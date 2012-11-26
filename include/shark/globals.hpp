/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_GLOBALS_HPP
#define __SHARK_GLOBALS_HPP

#include <string>                      // std::string
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
	 * The name of the scheduler.
	 */
	extern const std::string sched;

	/**
	 * Initialize shark.
	 * This has to be called at the beginning of the program.
	 * @param argc pointer to the argument count
	 * @param argv pointer to the argument values
	 */
	void Init(int* argc, char*** argv);

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
}

#endif
