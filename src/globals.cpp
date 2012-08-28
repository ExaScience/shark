
#include <cassert>                     // assert
#include <utility>                     // std::move

#if defined(SHARK_NO_COMM)
#include <cstdlib>                     // std::exit
#include <chrono>                      // std::chrono
#endif

#include <shark/globals.hpp>
#include "comm_impl.hpp"
#include "sched_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

#if defined(SHARK_SER_SCHED)
const string shark::sched("ser");
int shark::nthrds(1);

#elif defined(SHARK_PTHREAD_SCHED)
const string shark::sched("pthread");
int shark::nthrds(8);

#elif defined(SHARK_OMP_SCHED)
const string shark::sched("omp");
int shark::nthrds(omp_get_num_procs());

#elif defined(SHARK_TBB_SCHED)
const string shark::sched("tbb");
int shark::nthrds(tbb::task_scheduler_init::default_num_threads());

#elif defined(SHARK_COBRA_SCHED)
const string shark::sched("cobra");
int shark::nthrds(8);

#else
#error "No scheduler nthrds"

#endif

#include "types"

void shark::Init(int* argc, char*** argv) {
#if defined(SHARK_MPI_COMM)
	MPI_Init(argc, argv);
	{
		unique_ptr<GroupImpl> impl(new GroupImpl());
		impl->comm = MPI_COMM_WORLD;
		Group::w.reset(new Group(move(impl)));
	}
#elif defined(SHARK_NO_COMM)
	unused(argc, argv);
	Group::w.reset(new Group(unique_ptr<GroupImpl>(new GroupImpl())));
#endif
}

void shark::Finalize() {
#if defined(SHARK_PTHREAD_SCHED)
	pthread_barrier_wait(&barrier);
	for(int k = 0; k < nthrds-1; k++)
		pthread_join(children[k], nullptr);
	delete[] children;
	pthread_barrier_destroy(&barrier);
#elif defined(SHARK_TBB_SCHED)
	delete tbb_init;
#elif defined(SHARK_COBRA_SCHED)
	delete sch;
#endif

#if defined(SHARK_MPI_COMM)
	Group::w.reset();
	MPI_Finalize();
#elif defined(SHARK_NO_COMM)
	Group::w.reset();
#endif
}

void shark::SetupThreads() {
#if defined(SHARK_SER_SCHED)
	assert(nthrds == 1);
#elif defined(SHARK_PTHREAD_SCHED)
	work = 0;
	pthread_barrier_init(&barrier, nullptr, nthrds);
	children = new pthread_t[nthrds-1];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	for(int k = 0; k < nthrds-1; k++)
		pthread_create(&children[k], &attr, ThreadMain, reinterpret_cast<void*>(static_cast<size_t>(k+1)));
	pthread_attr_destroy(&attr);
#elif defined(SHARK_OMP_SCHED)
	omp_set_num_threads(nthrds);
#elif defined(SHARK_TBB_SCHED)
	tbb_init = new tbb::task_scheduler_init(nthrds);
#elif defined(SHARK_COBRA_SCHED)
	sch = new cobra::scheduler(nthrds);
#endif
}

void shark::Abort(int errorcode) {
#if defined(SHARK_MPI_COMM)
	MPI_Abort(MPI_COMM_WORLD, errorcode);
#elif defined(SHARK_NO_COMM)
	exit(errorcode);
#else
#error "No comm abort"
#endif
}

double shark::Wtime() {
#if defined(SHARK_MPI_COMM)
	return MPI_Wtime();
#elif defined(SHARK_NO_COMM)
	return chrono::duration_cast<chrono::duration<double>>(chrono::high_resolution_clock::now().time_since_epoch()).count();
#else
#error "No comm wtime"
#endif
}

#if defined(SHARK_PTHREAD_SCHED)
void shark::ThreadWork(function<void(int)> w) {
	work = w;
	pthread_barrier_wait(&barrier);
	work(0);
	pthread_barrier_wait(&barrier);
	work = 0;
}
#endif

