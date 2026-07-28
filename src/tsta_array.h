#ifndef TSTA_ARRAY_H
#define TSTA_ARRAY_H

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;

#ifndef roundup64
#define roundup64(x)                                                           \
  ({                                                                           \
    u64 __x = (x);                                                             \
    __x--;                                                                     \
    __x |= __x >> 1;                                                           \
    __x |= __x >> 2;                                                           \
    __x |= __x >> 4;                                                           \
    __x |= __x >> 8;                                                           \
    __x |= __x >> 16;                                                          \
    __x |= __x >> 32;                                                          \
    __x++;                                                                     \
    __x;                                                                       \
  })
#endif

#define __define_array_struct(type, name)                                      \
  typedef struct {                                                             \
    type *data;                                                                \
    u64 length;                                                                \
    u64 capacity;                                                              \
  } name

#define _min_array_capacity 1
#define array_auto_resize(arr, type)                                           \
  if ((arr)->length == (arr)->capacity) {                                      \
    if ((arr)->capacity == 0) {                                                \
      (arr)->data = (type *)malloc(sizeof(type) * _min_array_capacity);        \
      (arr)->capacity = _min_array_capacity;                                   \
      memset((arr)->data, 0, sizeof(type) * _min_array_capacity);              \
    } else {                                                                   \
      (arr)->capacity = roundup64((u64)((arr)->capacity + 1ULL));              \
      (arr)->data =                                                            \
          (type *)realloc((arr)->data, sizeof(type) * (arr)->capacity);        \
      memset((arr)->data + (arr)->length, 0,                                   \
             sizeof(type) * ((arr)->capacity - (arr)->length));                \
    }                                                                          \
  }

#define __array_method(type, name, method_name)                                \
  static inline void method_name##_init__(name *a) {                           \
    a->data = (type *)malloc(_min_array_capacity * sizeof(type));              \
    a->length = 0;                                                             \
    a->capacity = _min_array_capacity;                                         \
  }                                                                            \
  static inline name *method_name##_new() {                                    \
    name *arr = (name *)malloc(sizeof(name));                                  \
    arr->data = NULL;                                                          \
    arr->length = 0;                                                           \
    arr->capacity = 0;                                                         \
    return arr;                                                                \
  }                                                                            \
  static inline name *method_name##_with_capacity(u64 capacity) {              \
    name *arr = (name *)malloc(sizeof(name));                                  \
    arr->data = (type *)malloc(sizeof(type) * capacity);                       \
    memset(arr->data, 0, sizeof(type) * capacity);                             \
    arr->length = 0;                                                           \
    arr->capacity = capacity;                                                  \
    return arr;                                                                \
  }                                                                            \
  static inline void method_name##_shrink(name *a) {                           \
    if (a->capacity > a->length) {                                             \
      a->data = (type *)realloc(a->data, sizeof(type) * a->length);            \
      a->capacity = a->length;                                                 \
    }                                                                          \
  }                                                                            \
  static inline void method_name##_free(name *a) {                             \
    if (!a)                                                                    \
      return;                                                                  \
    if (a->capacity) {                                                         \
      free(a->data);                                                           \
      a->data = NULL;                                                          \
    }                                                                          \
    a->length = 0;                                                             \
    a->capacity = 0;                                                           \
    free(a);                                                                   \
  }                                                                            \
  static inline void method_name##_push(name *a, type value) {                 \
    array_auto_resize(a, type);                                                \
    a->data[a->length++] = value;                                              \
  }                                                                            \
  static inline int method_name##_reserve(name *a, u64 min_capacity) {         \
    type *new_data;                                                            \
    u64 new_capacity;                                                          \
    u64 old_capacity;                                                          \
    if (!a)                                                                    \
      return -1;                                                               \
    if (min_capacity <= a->capacity)                                           \
      return 0;                                                                \
    old_capacity = a->capacity;                                                \
    new_capacity = old_capacity > 0 ? old_capacity : _min_array_capacity;      \
    while (new_capacity < min_capacity) {                                      \
      if (new_capacity > UINT64_MAX / 2) {                                     \
        new_capacity = min_capacity;                                           \
        break;                                                                 \
      }                                                                        \
      new_capacity *= 2;                                                       \
    }                                                                          \
    new_data = (type *)realloc(a->data, sizeof(type) * new_capacity);          \
    if (!new_data)                                                             \
      return -1;                                                               \
    if (new_capacity > old_capacity) {                                         \
      memset(new_data + old_capacity, 0,                                       \
             sizeof(type) * (new_capacity - old_capacity));                    \
    }                                                                          \
    a->data = new_data;                                                        \
    a->capacity = new_capacity;                                                \
    return 0;                                                                  \
  }                                                                            \
  static inline int method_name##_append(name *a, type value) {                \
    if (method_name##_reserve(a, a->length + 1) != 0)                          \
      return -1;                                                               \
    a->data[a->length++] = value;                                              \
    return 0;                                                                  \
  }                                                                            \
  static inline type method_name##_get(name *a, u64 index) {                   \
    return a->data[index];                                                     \
  }                                                                            \
  static inline void method_name##_set(name *a, u64 index, type value) {       \
    a->data[index] = value;                                                    \
  }                                                                            \
  static inline void method_name##_clear(name *a) { a->length = 0; }           \
  static inline void method_name##_zero(name *a) {                             \
    memset(a->data, 0, sizeof(type) * a->length);                              \
  }                                                                            \
  static inline void method_name##_resize(name *a, u64 new_size) {             \
    assert(new_size > 0);                                                      \
    new_size = roundup64(new_size);                                            \
    if (new_size > a->capacity) {                                              \
      a->capacity = new_size;                                                  \
      a->data = (type *)realloc(a->data, sizeof(type) * a->capacity);          \
      memset(a->data + a->length, 0,                                           \
             sizeof(type) * (a->capacity - a->length));                        \
    }                                                                          \
  }                                                                            \
  static inline type method_name##_pop(name *a) {                              \
    if (a->length == 0)                                                        \
      return (type){0};                                                        \
    return a->data[--a->length];                                               \
  }                                                                            \
  static inline type method_name##_ipop(name *a, u64 i) {                      \
    type el = a->data[i];                                                      \
    a->data[i] = a->data[--a->length];                                         \
    return el;                                                                 \
  }                                                                            \
  static inline type *method_name##_nextref(name *a) {                         \
    array_auto_resize(a, type);                                                \
    return &a->data[a->length++];                                              \
  }                                                                            \
  static inline type method_name##_last(name *a) {                             \
    return a->data[a->length - 1];                                             \
  }                                                                            \
  static inline type method_name##_first(name *a) { return a->data[0]; }       \
  static inline name *method_name##_concat(name *a, name *b) {                 \
    if (a->length + b->length > a->capacity) {                                 \
      u64 new_cap = a->length + b->length;                                     \
      method_name##_resize(a, new_cap);                                        \
    }                                                                          \
    memcpy(a->data + a->length, b->data, sizeof(type) * b->length);            \
    a->length += b->length;                                                    \
    return a;                                                                  \
  }                                                                            \
  static inline u64 method_name##_memory_usage(name *a) {                      \
    if (!a)                                                                    \
      return 0;                                                                \
    return sizeof(name) + sizeof(type) * a->capacity;                          \
  }

#define define_array(type, name, method_name)                                  \
  __define_array_struct(type, name);                                           \
  __array_method(type, name, method_name)

#define array_malloc(arr, type, size)                                          \
  do {                                                                         \
    (arr)->data = (type *)malloc(sizeof(type) * size);                         \
    (arr)->length = 0;                                                         \
    (arr)->capacity = size;                                                    \
  } while (0)

#define array_realloc(arr, type, size)                                         \
  do {                                                                         \
    (arr)->data = (type *)realloc((arr)->data, sizeof(type) * size);           \
    (arr)->capacity = size;                                                    \
  } while (0)

#define array_free(a, block)                                                   \
  do {                                                                         \
    block;                                                                     \
    {                                                                          \
      if ((a)->capacity) {                                                     \
        free((a)->data);                                                       \
        (a)->data = NULL;                                                      \
      }                                                                        \
      (a)->length = 0;                                                         \
      (a)->capacity = 0;                                                       \
      free(a);                                                                 \
    }                                                                          \
  } while (0)
#endif
