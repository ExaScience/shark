
#ifndef __SHARK_EXPR_HPP
#define __SHARK_EXPR_HPP

#include <cassert>
#include "common.hpp"
#include "coords.hpp"

namespace shark {

	namespace ndim {

		template<int,typename>
		class ConstAcc;

		template<int ndim,typename T>
		class ConstExp {
			friend class ConstAcc<ndim,T>;
			const Domain<ndim>& dom;
			const T val;
		public:
			static const int number_of_dimensions = ndim;
			typedef ConstAcc<ndim,T> accessor;
			ConstExp(const Domain<ndim>& dom, const T& val);
			INLINE const Domain<ndim>& domain() const;
		};

		template<int ndim,typename T>
		class ConstAcc {
			const ConstExp<ndim,T>& c;
		public:
			ConstAcc(const ConstExp<ndim,T>& c);
			INLINE const T& operator()(coords<ndim>) const;
		};

		// Inline Const member functions

		template<int ndim, typename T>
		inline const Domain<ndim>& ConstExp<ndim,T>::domain() const {
			return dom;
		}

		template<int ndim, typename T>
		inline const T& ConstAcc<ndim,T>::operator()(coords<ndim>) const {
			return c.val;
		}

		/**
		 * Constant value over domain
		 */
		template<int ndim,typename T>
		ConstExp<ndim,T> constant(const Domain<ndim>& dom, const T& val);

		template<typename>
		class NegAcc;

		template<typename S>
		class NegExp {
			friend class NegAcc<S>;
			const S& src;
		public:
			static const int number_of_dimensions = S::number_of_dimensions;
			typedef NegAcc<S> accessor;
			NegExp(const S& src): src(src) { }
			const Domain<number_of_dimensions>& domain() const {
				return src.domain();
			}
		};

		template<typename S>
		class NegAcc {
			const typename S::accessor a;
		public:
			NegAcc(const NegExp<S>& exp): a(exp.src) { }
			auto operator()(coords<S::number_of_dimensions> ii) const -> decltype(-a(ii)) {
				return -a(ii);
			}
		};

		/**
		 * Negation
		 */
		template<typename S>
		NegExp<S> operator-(const S& src) {
			return NegExp<S>(src);
		}

		template<typename,typename>
		class AddAcc;

		template<typename S1, typename S2>
		class AddExp {
			friend class AddAcc<S1,S2>;
			const S1& src1;
			const S2& src2;
		public:
			static_assert(S1::number_of_dimensions == S2::number_of_dimensions, "source dimensionality");
			static const int number_of_dimensions = S1::number_of_dimensions;
			typedef AddAcc<S1,S2> accessor;
			AddExp(const S1& src1, const S2& src2): src1(src1), src2(src2) { }
			const Domain<number_of_dimensions>& domain() const {
				assert(src1.domain() == src2.domain());
				return src1.domain();
			}
		};

		template<typename S1, typename S2>
		class AddAcc {
			const typename S1::accessor a1;
			const typename S2::accessor a2;
		public:
			AddAcc(const AddExp<S1,S2>& exp): a1(exp.src1), a2(exp.src2) { }
			auto operator()(coords<S1::number_of_dimensions> ii) const -> decltype(a1(ii) + a2(ii)) {
				return a1(ii) + a2(ii);
			}
		};

		/**
		 * Addition
		 */
		template<typename S1, typename S2>
		AddExp<S1,S2> operator+(const S1& src1, const S2& src2) {
			return AddExp<S1,S2>(src1, src2);
		}
	}

}

#endif
