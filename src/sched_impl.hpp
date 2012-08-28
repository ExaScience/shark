#ifndef __SHARK_SCHED_IMPL_HPP
#define __SHARK_SCHED_IMPL_HPP

#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)
#include <pthread.h>
#include <functional>     // std::function
#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif

namespace shark {
#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)
	extern pthread_barrier_t barrier;
	extern std::function<void(int)> work;
	extern pthread_t* children;

	void* ThreadMain(void* arg);
#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif
}

#endif
