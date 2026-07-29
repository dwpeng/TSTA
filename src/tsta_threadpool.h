#ifndef TSTA_THREADPOOL_H
#define TSTA_THREADPOOL_H

#include <stddef.h>

typedef struct tsta_threadpool tsta_threadpool_t;

tsta_threadpool_t*
tsta_threadpool_create(int thread_count, int queue_capacity, size_t task_size);
int tsta_threadpool_destroy(tsta_threadpool_t* pool);
int tsta_threadpool_submit(tsta_threadpool_t* pool,
                           void (*function)(void*),
                           const void* arg,
                           size_t arg_size);
void tsta_threadpool_wait(tsta_threadpool_t* pool);

#endif
