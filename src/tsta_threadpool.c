#include "tsta_threadpool.h"

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define TSTA_TP_PAUSE() _mm_pause()
#else
#define TSTA_TP_PAUSE() ((void)0)
#endif

/* Spin iterations before a worker / waiter falls back to a condvar wait.
 * Bounded busy-wait avoids the futex syscall on the short, bursty task
 * pattern TSTA uses (submit a few blocks, barrier-wait, repeat). */
#define TSTA_TP_SPIN 1000

typedef struct tsta_threadpool_task {
  void (*function)(void* arg);
  unsigned char* payload;
} tsta_threadpool_task_t;

typedef struct tsta_threadpool_worker_context {
  tsta_threadpool_t* pool;
  void* task_buffer;
} tsta_threadpool_worker_context_t;

struct tsta_threadpool {
  /* Shared state, protected by mutex_pool. */
  _Alignas(64) pthread_mutex_t mutex_pool;
  pthread_cond_t not_empty;
  pthread_cond_t finished;
  tsta_threadpool_task_t* task_queue;
  unsigned char* task_storage;
  size_t task_size;
  size_t task_stride;
  int queue_capacity;
  int queue_size;
  int queue_front;
  int queue_rear;
  int shutdown;
  /* Hot barrier counters, atomically updated; cache-line isolated from the
   * mutex/queue above and the read-only data below to cut false sharing.
   *   pending: submitted, not yet finished (waited on by tsta_threadpool_wait)
   *   queued:  submitted, not yet dequeued (spun on by workers)             */
  _Alignas(64) int pending;
  int queued;
  /* Read-only after create. */
  _Alignas(64) unsigned char* worker_storage;
  tsta_threadpool_worker_context_t* worker_contexts;
  pthread_t* threads;
  int thread_count;
};

static size_t
tsta_threadpool_align_size(size_t size)
{
  size_t alignment = (size_t)_Alignof(max_align_t);
  if (alignment == 0)
    return size;
  return ((size + alignment - 1) / alignment) * alignment;
}

static void*
tsta_threadpool_worker(void* arg)
{
  tsta_threadpool_worker_context_t* context =
      (tsta_threadpool_worker_context_t*)arg;
  tsta_threadpool_t* pool = context->pool;
  void* task_buffer = context->task_buffer;

  for (;;) {
    /* Phase 1: bounded busy-wait for a queued task (no syscall). */
    int spins = 0;
    while (spins < TSTA_TP_SPIN
           && __atomic_load_n(&pool->queued, __ATOMIC_ACQUIRE) == 0
           && !__atomic_load_n(&pool->shutdown, __ATOMIC_ACQUIRE)) {
      TSTA_TP_PAUSE();
      spins++;
    }

    /* Phase 2: blocking acquire. */
    pthread_mutex_lock(&pool->mutex_pool);
    while (pool->queue_size == 0
           && !__atomic_load_n(&pool->shutdown, __ATOMIC_ACQUIRE))
      pthread_cond_wait(&pool->not_empty, &pool->mutex_pool);

    if (__atomic_load_n(&pool->shutdown, __ATOMIC_ACQUIRE)
        && pool->queue_size == 0) {
      pthread_mutex_unlock(&pool->mutex_pool);
      break;
    }

    tsta_threadpool_task_t task = pool->task_queue[pool->queue_front];
    pool->queue_front = (pool->queue_front + 1) % pool->queue_capacity;
    pool->queue_size--;
    /* Copy into the worker's private buffer while holding the lock: the ring
     * slot becomes free immediately and the producer may overwrite it. */
    memcpy(task_buffer, task.payload, pool->task_size);
    pthread_mutex_unlock(&pool->mutex_pool);
    __atomic_sub_fetch(&pool->queued, 1, __ATOMIC_ACQ_REL);

    /* Execute outside the lock; each worker has its own task buffer. */
    task.function(task_buffer);

    /* Completion: the last in-flight task signals the waiter. */
    if (__atomic_sub_fetch(&pool->pending, 1, __ATOMIC_ACQ_REL) == 0) {
      pthread_mutex_lock(&pool->mutex_pool);
      pthread_cond_signal(&pool->finished);
      pthread_mutex_unlock(&pool->mutex_pool);
    }
  }
  return NULL;
}

tsta_threadpool_t*
tsta_threadpool_create(int thread_count, int queue_capacity, size_t task_size)
{
  tsta_threadpool_t* pool;
  size_t task_stride;
  int mutex_initialized = 0;
  int not_empty_initialized = 0;
  int finished_initialized = 0;

  if (thread_count <= 0 || queue_capacity <= 0 || task_size == 0)
    return NULL;

  task_stride = tsta_threadpool_align_size(task_size);
  if (task_stride == 0 || (size_t)queue_capacity > SIZE_MAX / task_stride
      || (size_t)thread_count > SIZE_MAX / task_stride)
    return NULL;

  pool = (tsta_threadpool_t*)calloc(1, sizeof(*pool));
  if (!pool)
    return NULL;

  pool->threads =
      (pthread_t*)calloc((size_t)thread_count, sizeof(*pool->threads));
  pool->worker_contexts = (tsta_threadpool_worker_context_t*)calloc(
      (size_t)thread_count, sizeof(*pool->worker_contexts));
  pool->task_queue = (tsta_threadpool_task_t*)calloc(
      (size_t)queue_capacity, sizeof(*pool->task_queue));
  pool->task_storage =
      (unsigned char*)calloc((size_t)queue_capacity, task_stride);
  pool->worker_storage =
      (unsigned char*)calloc((size_t)thread_count, task_stride);
  if (!pool->threads || !pool->worker_contexts || !pool->task_queue
      || !pool->task_storage || !pool->worker_storage)
    goto fail;

  pool->task_size = task_size;
  pool->task_stride = task_stride;
  pool->queue_capacity = queue_capacity;
  pool->thread_count = thread_count;

  for (int i = 0; i < queue_capacity; ++i)
    pool->task_queue[i].payload =
        pool->task_storage + (size_t)i * pool->task_stride;

  if (pthread_mutex_init(&pool->mutex_pool, NULL) != 0)
    goto fail;
  mutex_initialized = 1;
  if (pthread_cond_init(&pool->not_empty, NULL) != 0)
    goto fail;
  not_empty_initialized = 1;
  if (pthread_cond_init(&pool->finished, NULL) != 0)
    goto fail;
  finished_initialized = 1;

  for (int i = 0; i < thread_count; ++i) {
    pool->worker_contexts[i].pool = pool;
    pool->worker_contexts[i].task_buffer =
        pool->worker_storage + (size_t)i * pool->task_stride;
    if (pthread_create(&pool->threads[i], NULL, tsta_threadpool_worker,
                       &pool->worker_contexts[i])
        != 0) {
      pthread_mutex_lock(&pool->mutex_pool);
      __atomic_store_n(&pool->shutdown, 1, __ATOMIC_RELEASE);
      pthread_cond_broadcast(&pool->not_empty);
      pthread_mutex_unlock(&pool->mutex_pool);
      for (int j = 0; j < i; ++j)
        pthread_join(pool->threads[j], NULL);
      goto fail;
    }
  }
  return pool;

fail:
  if (finished_initialized)
    pthread_cond_destroy(&pool->finished);
  if (not_empty_initialized)
    pthread_cond_destroy(&pool->not_empty);
  if (mutex_initialized)
    pthread_mutex_destroy(&pool->mutex_pool);
  free(pool->worker_storage);
  free(pool->task_storage);
  free(pool->task_queue);
  free(pool->worker_contexts);
  free(pool->threads);
  free(pool);
  return NULL;
}

int
tsta_threadpool_destroy(tsta_threadpool_t* pool)
{
  if (!pool)
    return -1;

  pthread_mutex_lock(&pool->mutex_pool);
  while (__atomic_load_n(&pool->pending, __ATOMIC_ACQUIRE) > 0)
    pthread_cond_wait(&pool->finished, &pool->mutex_pool);
  __atomic_store_n(&pool->shutdown, 1, __ATOMIC_RELEASE);
  pthread_cond_broadcast(&pool->not_empty);
  pthread_mutex_unlock(&pool->mutex_pool);

  for (int i = 0; i < pool->thread_count; ++i)
    pthread_join(pool->threads[i], NULL);

  pthread_mutex_destroy(&pool->mutex_pool);
  pthread_cond_destroy(&pool->not_empty);
  pthread_cond_destroy(&pool->finished);
  free(pool->worker_storage);
  free(pool->task_storage);
  free(pool->task_queue);
  free(pool->worker_contexts);
  free(pool->threads);
  free(pool);
  return 0;
}

int
tsta_threadpool_submit(tsta_threadpool_t* pool,
                       void (*function)(void*),
                       const void* arg,
                       size_t arg_size)
{
  if (!pool || !function || !arg || arg_size != pool->task_size)
    return -1;

  pthread_mutex_lock(&pool->mutex_pool);
  if (__atomic_load_n(&pool->shutdown, __ATOMIC_ACQUIRE)) {
    pthread_mutex_unlock(&pool->mutex_pool);
    return -1;
  }
  if (pool->queue_size == pool->queue_capacity) {
    pthread_mutex_unlock(&pool->mutex_pool);
    return 1;
  }

  tsta_threadpool_task_t* slot = &pool->task_queue[pool->queue_rear];
  memcpy(slot->payload, arg, pool->task_size);
  slot->function = function;
  pool->queue_rear = (pool->queue_rear + 1) % pool->queue_capacity;
  pool->queue_size++;
  __atomic_add_fetch(&pool->queued, 1, __ATOMIC_ACQ_REL);
  __atomic_add_fetch(&pool->pending, 1, __ATOMIC_ACQ_REL);
  /* Signal only on the empty->non-empty transition; busy workers will see
   * queued via their spin loop. */
  if (pool->queue_size == 1)
    pthread_cond_signal(&pool->not_empty);
  pthread_mutex_unlock(&pool->mutex_pool);
  return 0;
}

void
tsta_threadpool_wait(tsta_threadpool_t* pool)
{
  if (!pool)
    return;

  /* Bounded busy-wait first: the common case (short tasks) drains before
   * the futex round-trip cost of a condvar wait. */
  int spins = 0;
  while (spins < TSTA_TP_SPIN
         && __atomic_load_n(&pool->pending, __ATOMIC_ACQUIRE) > 0) {
    TSTA_TP_PAUSE();
    spins++;
  }
  if (__atomic_load_n(&pool->pending, __ATOMIC_ACQUIRE) == 0)
    return;

  pthread_mutex_lock(&pool->mutex_pool);
  while (__atomic_load_n(&pool->pending, __ATOMIC_ACQUIRE) > 0)
    pthread_cond_wait(&pool->finished, &pool->mutex_pool);
  pthread_mutex_unlock(&pool->mutex_pool);
}
