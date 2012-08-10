
#include <cassert>                     // assert
#include <utility>                     // std::move

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
	MPI_Init(argc, argv);
	{
		unique_ptr<GroupImpl> impl(new GroupImpl());
		impl->comm = MPI_COMM_WORLD;
		Group::w.reset(new Group(move(impl)));
	}
#define SYMB(d,T) mpi_type<vec<d,T>>::init(); 
#include "inst_dimtype"
#undef SYMB
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
#define SYMB(d,T) mpi_type<vec<d,T>>::destroy();
#include "inst_dimtype"
#undef SYMB
	Group::w.reset();
	MPI_Finalize();
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
	MPI_Abort(MPI_COMM_WORLD, errorcode);
}

double shark::Wtime() {
	return MPI_Wtime();
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

