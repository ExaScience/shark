/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

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
	extern volatile int epoch;
	extern pthread_spinlock_t rdy_count_lock;
	extern std::function<void(int)> work;
	extern pthread_t* children;

	void* ThreadMain(void* arg);
#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif
}

#endif
