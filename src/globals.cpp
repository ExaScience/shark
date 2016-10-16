/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <cassert>                     // assert
#include <cstdlib>                     // std::getenv
#include <utility>                     // std::move
#include <iostream>
#include <iomanip>

#if defined(SHARK_NO_COMM)
#include <cstdlib>                     // std::exit
#include <chrono>                      // std::chrono
#endif

#if defined(SHARK_SER_SCHED)
#elif defined(SHARK_PTHREAD_SCHED)
#include <pthread.h>
#include <functional>                  // std::function

namespace shark {
	extern volatile int rdy_count;
	extern volatile int epoch;
	extern pthread_spinlock_t rdy_count_lock;
	extern std::function<void(int)> work;
	extern pthread_t* children;

	void* ThreadMain(void* arg);

}
#elif defined(SHARK_OMP_SCHED)
#include <omp.h>
#elif defined(SHARK_TBB_SCHED)
#include <tbb/task_scheduler_init.h>

namespace shark {
	extern tbb::task_scheduler_init* tbb_init;
}
#elif defined(SHARK_COBRA_SCHED)
#endif

#include <shark/globals.hpp>
#include <shark/version.hpp>
#include "comm_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

#define xstr(a) str(a)
#define str(a) #a

const string shark::version(xstr(SHARK_VERSION));

#undef xstr
#undef str

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

#if defined(SHARK_MPI_COMM)
#ifdef SHARK_STPACK
const string shark::comm("mpi-stpack");
#else
const string shark::comm("mpi");
#endif
#elif defined(SHARK_GPI_COMM)
#ifdef SHARK_STPACK
const string shark::comm("gpi-stpack");
#else
const string shark::comm("gpi");
#endif
#elif defined(SHARK_NO_COMM)
const string shark::comm("nocomm");
#else
#error "No comm"

#endif

verbosity_mask shark::log_mask = verbosity_mask();
std::ostream* shark::log_out(nullptr);

#include "types"

#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
namespace {
	inline MPI_Aint offset(void* addr) {
		return static_cast<MPI_Aint>(static_cast<char*>(addr) - static_cast<char*>(nullptr));
	}

	template<int ndim>
	void init_part_type() {
		part<ndim>* p = nullptr;
		int bl[] = { mpi_type<vec<ndim,part_position>>::count(), mpi_type<vec<ndim,double>>::count(), 1, 1 };
		MPI_Aint disp[] = { offset(&p->x), offset(&p->v), offset(&p->w), offset(p+1) };
		MPI_Datatype types[] = { mpi_type<vec<ndim,part_position>>::t, mpi_type<vec<ndim,double>>::t, MPI_DOUBLE, MPI_UB };
		MPI_Type_create_struct(4, bl, disp, types, &mpi_type<part<ndim>>::t);
		MPI_Type_commit(&mpi_type<part<ndim>>::t);
	}

	template<int ndim>
	void free_part_type() {
		MPI_Type_free(&mpi_type<part<ndim>>::t);
	}
}

void shark::Init(int* argc, char*** argv, bool thread_multiple) {
        if (thread_multiple) {
                int provided;
                MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
                assert(provided == MPI_THREAD_MULTIPLE);
        } else {
                MPI_Init(argc, argv);
        }

#if defined(SHARK_GPI_COMM)
	gaspi_proc_init(GASPI_BLOCK);
#endif

	{
		unique_ptr<GroupImpl> impl(new GroupImpl());
		impl->comm = MPI_COMM_WORLD;
#if defined(SHARK_GPI_COMM)
		impl->grp = GASPI_GROUP_ALL;
#endif
		Group::w.reset(new Group(move(impl)));
	}
	{
		unique_ptr<GroupImpl> impl(new GroupImpl());
		impl->comm = MPI_COMM_SELF;
#if defined(SHARK_GPI_COMM)
                gaspi_rank_t my_rank;
                gaspi_proc_rank(&my_rank);
                gaspi_group_create(&impl->grp);
                gaspi_group_add(impl->grp, my_rank);
#endif
		Group::s.reset(new Group(move(impl)));
	}
        {
		part_position* pp = nullptr;
		int bl[] = { 1, 1, 1 };
		MPI_Aint disp[] = { offset(&pp->pos), offset(&pp->off), offset(pp+1) };
		MPI_Datatype types[] = { MPI_INT, MPI_FLOAT, MPI_UB };
		MPI_Type_create_struct(3, bl, disp, types, &mpi_type<part_position>::t);
		MPI_Type_commit(&mpi_type<part_position>::t);
	}
#define SYMBD(d) init_part_type<d>();
#include "inst_dim"
#undef SYMBD
#elif defined(SHARK_NO_COMM)
void shark::Init(int*, char***, bool) {
	Group::w.reset(new Group(unique_ptr<GroupImpl>(new GroupImpl())));
	Group::s.reset(new Group(unique_ptr<GroupImpl>(new GroupImpl())));
#endif
	char* nthrds_str = getenv("SHARK_NTHRDS");
	if(nthrds_str != nullptr)
		nthrds = stoi(nthrds_str);

#ifdef NDEBUG
        log_out = new std::ostream(0);
#else
        log_out = &std::cout;
#endif


}

void shark::Finalize() {
#if defined(SHARK_PTHREAD_SCHED)
	// Set expected
	rdy_count = nthrds-1;
	// Release threads
	epoch++;
	// Busy-wait for expected
	while(rdy_count > 0)
		;
	// Cleanup
	for(int k = 0; k < nthrds-1; k++)
		pthread_join(children[k], nullptr);
	delete[] children;
	pthread_spin_destroy(&rdy_count_lock);
#elif defined(SHARK_TBB_SCHED)
	delete tbb_init;
#elif defined(SHARK_COBRA_SCHED)
	delete sch;
#endif

#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
#define SYMBD(d) free_part_type<d>();
#include "inst_dim"
#undef SYMBD
	MPI_Type_free(&mpi_type<part_position>::t);
	Group::w.reset();
	Group::s.reset();
	MPI_Finalize();
#if defined(SHARK_GPI_COMM)
	gaspi_proc_term(GASPI_BLOCK);
#endif
#elif defined(SHARK_NO_COMM)
#endif
}

void shark::SetupThreads() {
#if defined(SHARK_SER_SCHED)
	assert(nthrds == 1);
#elif defined(SHARK_PTHREAD_SCHED)
	work = 0;
	pthread_spin_init(&rdy_count_lock, PTHREAD_PROCESS_PRIVATE);
	children = new pthread_t[nthrds-1];
	// Set expected
	rdy_count = nthrds-1;
	epoch = 0;
	// Release threads
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	for(int k = 0; k < nthrds-1; k++)
		pthread_create(&children[k], &attr, ThreadMain, reinterpret_cast<void*>(static_cast<size_t>(k+1)));
	pthread_attr_destroy(&attr);
	// Busy wait for expected
	while(rdy_count > 0)
		;
#elif defined(SHARK_OMP_SCHED)
	omp_set_num_threads(nthrds);
#elif defined(SHARK_TBB_SCHED)
	tbb_init = new tbb::task_scheduler_init(nthrds);
#elif defined(SHARK_COBRA_SCHED)
	sch = new cobra::scheduler(nthrds);
#endif
}

void shark::Abort(int errorcode) {
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
	MPI_Abort(MPI_COMM_WORLD, errorcode);
#elif defined(SHARK_NO_COMM)
	exit(errorcode);
#else
#error "No comm abort"
#endif
}

double shark::Wtime() {
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
	return MPI_Wtime();
#elif defined(SHARK_NO_COMM)
	return chrono::duration_cast<chrono::duration<double>>(chrono::high_resolution_clock::now().time_since_epoch()).count();
#else
#error "No comm wtime"
#endif
}

#if defined(SHARK_PTHREAD_SCHED)

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

// from shark/common.hpp
void shark::ThreadWork(function<void(int)> w) {
	work = w;
	// Set expected
	rdy_count = nthrds-1;
	// Release threads
	epoch++;
	work(0);
	// Busy wait for expected
	while(rdy_count > 0)
		;
	work = 0;
}

#elif defined(SHARK_TBB_SCHED)

tbb::task_scheduler_init* shark::tbb_init;

#endif

#ifdef SHARK_PROFILING

shark::Counter::Counter(std::string name)
    : name(name), count(1), diff_wall(0), total_counter(false)
{
    shark::perf_data.read(start_wall);
}

shark::Counter::Counter() 
    : name(std::string()), count(0), diff_wall(0), total_counter(true)
{} 

shark::Counter::~Counter() {
    if(total_counter) return;

    shark::perf_data.read(stop_wall);
    diff_wall = stop_wall - start_wall;
    shark::perf_data[name] += *this;
}

void shark::Counter::operator+=(const shark::Counter &other) {
    if (name.empty()) name = other.name;
    diff_wall += other.diff_wall;
    count+= other.count;
}

std::string shark::Counter::as_string(const shark::Counter &total) {
    std::ostringstream os;
    int percent = 100.0 * diff_wall / (total.diff_wall + 0.000001);
    os << ">> " << name << ":\t" <<
        count << "\t" << std::fixed << std::setw(11)
        << std::setprecision(4) 
        << diff_wall << "\t(" << percent << "%)\n";
    return os.str();
}

shark::TotalsCounter shark::perf_data;

shark::TotalsCounter::TotalsCounter() {}

void shark::TotalsCounter::print() {
    if (world().procid != 0) return;
    std::cout << "\nTotals:\n";
    for(auto &t : data)
        std::cout << t.second.as_string(data["main"]);
}

void shark::TotalsCounter::read(double &w) {
    w = Wtime();
}

#endif // SHARK_PROFILING

