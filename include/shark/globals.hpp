/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_GLOBALS_HPP
#define __SHARK_GLOBALS_HPP

#include <string>                      // std::string
#include <map>
#include <iostream>
#include <sstream>

#include "common.hpp"
#include "group.hpp"

namespace shark {

	/**
	 * Shark version.
	 */
	extern const std::string version;

	/**
	 * The number of threads per process.
	 */
	extern int nthrds;

	/**
	 * The name of the scheduler and the comm library
	 */
	extern const std::string sched;
	extern const std::string comm;

	/**
	 * Initialize shark.
	 * This has to be called at the beginning of the program.
	 * @param argc pointer to the argument count
	 * @param argv pointer to the argument values
	 */
	void Init(int* argc, char*** argv, bool thread_multiple = false);

	/**
	 * Setup threads.
	 */
	void SetupThreads();

	/**
	 * Finalize shark.
	 * This has to be called at the end of the program.
	 */
	void Finalize();

	/**
	 * Abort shark execution.
	 */
	void Abort(int errorcode);

	/**
	 * Returns an elapsed time on the calling process (local).
	 */
	double Wtime();

	/**
	 * The group of all processes.
	 * This is available after initialization.
	 */
	INLINE const Group& world();

	// Inline function implementations

	inline const Group& world() {
		return *Group::w;
	}

	/**
	 * The group of only the current process
	 */
	INLINE const Group& self();

	// Inline function implementations

	inline const Group& self() {
		return *Group::s;
	}

#ifdef SHARK_PROFILING

#define SHARK_COUNTER(name) Counter c(name)

        struct Counter {
            std::string name;
          
            unsigned count;
            double start_wall, stop_wall, diff_wall; // wallclock time

            bool total_counter;

            Counter(std::string name);
            Counter();

            ~Counter();

            void operator+=(const Counter &other);

            std::string as_string(const Counter &total);
        };

        struct TotalsCounter {
        private:
            std::map<std::string, Counter> data;

        public:
            //c-tor starts PAPI
            TotalsCounter();

            //prints results
            void print();

            Counter &operator[](const std::string &name) {
                return data[name];
            }

            void read(double &w);
        };

        extern TotalsCounter perf_data;
#else 

#define SHARK_COUNTER(name) 

#endif //SHARK_PROFILING
}

#endif
