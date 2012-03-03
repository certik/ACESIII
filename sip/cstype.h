#ifndef _UNIXTHREADS_H_
#define _UNIXTHREADS_H_
typedef pthread_mutex_t CRITICAL_SECTION_t;

#define InitializeCriticalSection(_mut_ptr_) \
        pthread_mutex_init(_mut_ptr_, NULL)
#define DeleteCriticalSection(_mut_ptr_) \
        pthread_mutex_destroy(_mut_ptr_)
#define EnterCriticalSection(_mut_ptr_) \
        pthread_mutex_lock(_mut_ptr_)
#define LeaveCriticalSection(_mut_ptr_) \
        pthread_mutex_unlock(_mut_ptr_)

#endif /* _UNIXTHREADS_H_ */
