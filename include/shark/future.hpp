/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_FUTURE_HPP
#define __SHARK_FUTURE_HPP

#include <memory>                      // std::unique_ptr

#include "common.hpp"

namespace shark {

	/**
	 * A future value of type T
	 */
	template<typename T>
	class Future {
		bool done;
		std::unique_ptr<T> val;

		Future(bool done);
		Future(bool done, const T& val);
		Future(bool done, T&& val);

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

		bool test();
		const T& wait();

		INLINE operator const T&();
	};

	// Inline members
	
	template<typename T>
	inline Future<T>::operator const T&() {
		return wait();
	}

	class Handle {
	public:
		virtual ~Handle() = 0;
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
	 * A future with no value
	 */
	template<>
	class Future<void> {
		template<int,typename> friend class shark::ndim::GlobalArray;

		bool done;
		std::unique_ptr<Handle> h;

		Future(bool done, Handle* h);
	
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

		bool test();
		void wait();

	};

}

#endif
