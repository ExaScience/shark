/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_FUTURE_HPP
#define __SHARK_FUTURE_HPP

#include <memory>                      // std::unique_ptr
#include <set>
#include <functional>

#include "common.hpp"

namespace shark {

	class Handle {
	public:
                static std::set<std::shared_ptr<Handle>> ongoing;
                static void make_progress();

                Handle();
		virtual ~Handle();
		virtual bool test() = 0;
		virtual void wait() = 0;
	};

	class DoneHandle: public Handle {
	public:
		DoneHandle();
		virtual ~DoneHandle();
		virtual bool test();
		virtual void wait();
	};

	/**
	 * A future value of type T
	 */
	template<typename T>
	struct Future {
		friend class shark::Group;
		bool done;
		std::shared_ptr<Handle> h;
		std::unique_ptr<T> val;

		Future(std::shared_ptr<Handle> h);
		Future(std::shared_ptr<Handle> h, std::unique_ptr<T>&& val);

	public:
		/**
		 * Create an inactive Future
		 */
		Future();
		~Future();

		/**
		 * No copy semantics
		 */
		Future(const Future<T>&) = delete;
		Future<T>& operator=(const Future<T>&) = delete;


		/**
		 * Move semantics
		 */
		Future(Future<T>&& f);
		Future& operator=(Future<T>&& f);

		/**
		 * Test for completion. Returns immediately.
		 */
		bool test();

		/**
		 * Wait for completion and return value.
		 */
		const T& wait();
		INLINE operator const T&();
	};

	// Inline members
	
	template<typename T>
	inline Future<T>::operator const T&() {
		return wait();
	}

	/**
	 * A future with no value
	 */
	template<>
	struct Future<void> {
		template<int,typename> friend class shark::ndim::GlobalArray;

		bool done;
		std::shared_ptr<Handle> h;

		Future(std::shared_ptr<Handle> h);
	
	public:
		/**
		 * Create an inactive Future
		 */
		Future();
		~Future();

		/**
		 * No copy semantics
		 */
		Future(const Future<void>&) = delete;
		Future<void>& operator=(const Future<void>&) = delete;


		/**
		 * Move semantics
		 */
		Future(Future<void>&& f);
		Future& operator=(Future<void>&& f);

		/**
		 * Test for completion. Returns immediately.
		 */
		bool test();

		/**
		 * Wait for completion.
		 */
		void wait();
	};

}

#endif
