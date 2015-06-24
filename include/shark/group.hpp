/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_GROUP_HPP
#define __SHARK_GROUP_HPP

#include <memory>                      // std::unique_ptr
#include "common.hpp"
#include "future.hpp"

namespace shark {

	class GroupImpl;

	/**
	 * A group of processes
	 */
	class Group {
		template<int> friend class shark::ndim::Domain;
		template<int,typename> friend class shark::ndim::GlobalArray;
		template<int,typename> friend class shark::ndim::SparseArray;
		friend const Group& world();
		friend void Init(int*, char***);
		friend void Finalize();

	private:
		static std::unique_ptr<Group> w;

		const std::unique_ptr<GroupImpl> impl;
		
	public:
		/**
		 * The identifier of the calling process in the group (zero-based).
		 */
		const int procid;

		/**
		 * The number of processes in the group.
		 */
		const int nprocs;

	private:
		Group(std::unique_ptr<GroupImpl>&& impl);
	
	public:
		~Group();
		Group(const Group& other) = delete;
		Group& operator=(const Group& other) = delete;


		bool operator==(const Group& other) const;

		/**
		 * Synchronize all processes in the group (collective).
		 */
		void sync() const;

		template<typename T>
		T external_sum(T&& val) const;

		template<typename T>
		Future<T> external_isum(T&& val) const;

	};

}

#endif
