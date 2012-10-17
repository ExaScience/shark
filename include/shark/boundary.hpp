#ifndef __SHARK_BOUNDARY_HPP
#define __SHARK_BOUNDARY_HPP

#include <array>                       // std::array
#include <memory>                      // std::unique_ptr

#include "common.hpp"
#include "access.hpp"
#include "coords.hpp"
#include "coords_range.hpp"

namespace shark {

	namespace ndim {

		/**
		 * A boundary configuration for a global array
		 */
		template<int ndim, typename T>
		class Boundary {
			friend class GlobalArray<ndim,T>;

			class type {
			public:
				virtual ~type();
				virtual type* clone() const;
			};

			class periodic_type: public type {
			public:
				periodic_type();
				virtual ~periodic_type();
				virtual periodic_type* clone() const;
			};

			class fixed_type: public type {
			public:
				virtual void set(Access<ndim,T>& a, coords_range<ndim> r) const = 0;
			};

			class general_type: public type {
			public:
				virtual void set(Access<ndim,T>& a, coords_range<ndim> r, long k) const = 0;
			};

			template<typename Func>
			class fun_fixed_type: public fixed_type {
				Func f;
			public:
				fun_fixed_type(const Func& f);
				virtual ~fun_fixed_type();
				virtual fun_fixed_type<Func>* clone() const;
				virtual void set(Access<ndim,T>& a, coords_range<ndim> r) const;
			};

			template<typename Func>
			class fun_general_type: public general_type {
				Func f;
			public:
				fun_general_type(const Func& f);
				virtual ~fun_general_type();
				virtual fun_general_type<Func>* clone() const;
				virtual void set(Access<ndim,T>& a, coords_range<ndim> r, long k) const;
			};

			std::unique_ptr<type> t;

			Boundary(type* t);

		public:
			/**
			 * Create an unmanaged boundary
			 */
			Boundary();
			static std::array<Boundary<ndim,T>,ndim> all_unmanaged();

			~Boundary();

			/**
			 * Copy semantics
			 */
			Boundary(const Boundary<ndim,T>& other);
			Boundary<ndim,T>& operator=(const Boundary<ndim,T>& other);

			/**
			 * Create a periodic boundary
			 */
			static Boundary<ndim,T> periodic();
			static std::array<Boundary<ndim,T>,ndim> all_periodic();

			/**
			 * Create a boundary with a constant value over space and time
			 */
			static Boundary<ndim,T> constant(const T& val);
			static std::array<Boundary<ndim,T>,ndim> all_constant(const T& val);

			/**
			 * Create a boundary with values fixed over time
			 */
			template<typename Func>
			static Boundary<ndim,T> fixed(const Func& f);
			template<typename Func>
			static std::array<Boundary<ndim,T>,ndim> all_fixed(const Func& f);

			/**
			 * Create a general boundary
			 */
			template<typename Func>
			static Boundary<ndim,T> general(const Func& f);
			template<typename Func>
			static std::array<Boundary<ndim,T>,ndim> all_general(const Func& f);
		};

		// Generic members

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T>::fun_fixed_type<Func>::fun_fixed_type(const Func& f): f(f) { }

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T>::fun_fixed_type<Func>::~fun_fixed_type() { }

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T>::fun_fixed_type<Func>* Boundary<ndim,T>::fun_fixed_type<Func>::clone() const {
			return new fun_fixed_type<Func>(f);
		}

		template<int ndim, typename T> template<typename Func>
		void Boundary<ndim,T>::fun_fixed_type<Func>::set(Access<ndim,T>& a, coords_range<ndim> r) const {
			r.for_each([this,&a](coords<ndim> ii) {
				a(ii) = f(ii);
			});
		}

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T>::fun_general_type<Func>::fun_general_type(const Func& f): f(f) { }

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T>::fun_general_type<Func>::~fun_general_type() { }

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T>::fun_general_type<Func>* Boundary<ndim,T>::fun_general_type<Func>::clone() const {
			return new fun_general_type<Func>(f);
		}

		template<int ndim, typename T> template<typename Func>
		void Boundary<ndim,T>::fun_general_type<Func>::set(Access<ndim,T>& a, coords_range<ndim> r, long k) const {
			r.for_each([this,&a,k](coords<ndim> ii) {
				a(ii) = f(k, ii);
			});
		}

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T> Boundary<ndim,T>::fixed(const Func& f) {
			return Boundary<ndim,T>(new fun_fixed_type<Func>(f));
		}

		template<int ndim, typename T> template<typename Func>
		std::array<Boundary<ndim,T>,ndim> Boundary<ndim,T>::all_fixed(const Func& f) {
			std::array<Boundary<ndim,T>,ndim> bds;
			for(int d = 0; d < ndim; d++)
				bds[d].t.reset(new fun_fixed_type<Func>(f));
			return bds;
		}

		template<int ndim, typename T> template<typename Func>
		Boundary<ndim,T> Boundary<ndim,T>::general(const Func& f) {
			return Boundary<ndim,T>(new fun_general_type<Func>(f));
		}

		template<int ndim, typename T> template<typename Func>
		std::array<Boundary<ndim,T>,ndim> Boundary<ndim,T>::all_general(const Func& f) {
			std::array<Boundary<ndim,T>,ndim> bds;
			for(int d = 0; d < ndim; d++)
				bds[d].t.reset(new fun_general_type<Func>(f));
			return bds;
		}

	}

}

#endif
