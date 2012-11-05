#ifndef __SHARK_SCHED_IMPL_HPP
#define __SHARK_SCHED_IMPL_HPP

#include <cstddef>   // std

#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)
#include <pthread.h>
#include <functional>     // std::function
#include <shark/globals.hpp>
#elif defined(SHARK_OMP_SCHED)
#include <omp.h>
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif

namespace shark {
#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)
	extern volatile int rdy_count;
	extern pthread_mutex_t mutex;
	extern pthread_cond_t cond_go;
	extern pthread_cond_t cond_rdy;
	extern std::function<void(int)> work;
	extern pthread_t* children;

	void* ThreadMain(void* arg);
#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif
}

#endif
