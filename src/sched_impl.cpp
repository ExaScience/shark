/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include "sched_impl.hpp"

using namespace std;
using namespace shark;

#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)

volatile int shark::rdy_count;
volatile int shark::epoch;
pthread_spinlock_t shark::rdy_count_lock;
function<void(int)> shark::work;
pthread_t* shark::children;

void* shark::ThreadMain(void* arg) {
	int threadid = reinterpret_cast<size_t>(arg);
	int next_epoch = 0;
	while(true) {
		pthread_spin_lock(&rdy_count_lock);
		rdy_count--;
		pthread_spin_unlock(&rdy_count_lock);
		// busy wait
		next_epoch++;
		while(epoch != next_epoch)
			;
		if(!work)
			break;
		work(threadid);
	}
	pthread_spin_lock(&rdy_count_lock);
	rdy_count--;
	pthread_spin_unlock(&rdy_count_lock);
	return NULL;
}

#elif defined(SHARK_OMP_SCHED)
#elif defined(SHARK_TBB_SCHED)
#elif defined(SHARK_COBRA_SCHED)
#endif

