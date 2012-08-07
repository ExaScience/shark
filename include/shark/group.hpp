#ifndef __SHARK_GROUP_HPP
#define __SHARK_GROUP_HPP

#include <memory>                      // std::unique_ptr
#include "common.hpp"

namespace shark {

	class GroupImpl;

	/**
	 * A group of processes
	 */
	class Group {
		template<int> friend class shark::ndim::Domain;
		template<int,typename> friend class shark::ndim::GlobalArray;
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
	};

}

#endif
