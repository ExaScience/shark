
#include "sched_impl.hpp"

using namespace std;
using namespace shark;

#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)

volatile int shark::rdy_count;
pthread_mutex_t shark::mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t shark::cond_go = PTHREAD_COND_INITIALIZER;
pthread_cond_t shark::cond_rdy = PTHREAD_COND_INITIALIZER;
function<void(int)> shark::work;
pthread_t* shark::children;

void* shark::ThreadMain(void* arg) {
	int threadid = reinterpret_cast<size_t>(arg);
	while(true) {
		pthread_mutex_lock(&mutex);
		rdy_count--;
		pthread_cond_signal(&cond_rdy);
		pthread_cond_wait(&cond_go, &mutex);
		if(!work) {
			rdy_count--;
			pthread_cond_signal(&cond_rdy);
			pthread_mutex_unlock(&mutex);
			break;
		}
		pthread_mutex_unlock(&mutex);
		work(threadid);
	}
	return NULL;
}

#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif

