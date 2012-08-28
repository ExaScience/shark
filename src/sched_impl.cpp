
#include "sched_impl.hpp"

using namespace std;
using namespace shark;

#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)

pthread_barrier_t shark::barrier;
function<void(int)> shark::work;
pthread_t* shark::children;

void* shark::ThreadMain(void* arg) {
	int threadid = reinterpret_cast<size_t>(arg);
	while(true) {
		pthread_barrier_wait(&barrier);
		if(!work) break;
		work(threadid);
		pthread_barrier_wait(&barrier);
	}
	return NULL;
}

#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif

