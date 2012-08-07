
#ifndef __SHARK_EXPR_HPP
#define __SHARK_EXPR_HPP

#include <cassert>
#include <type_traits>
#include "common.hpp"
#include "coords.hpp"

namespace shark {

	namespace ndim {

		// Using SFINAE to determine whether S can be considered a source
		template<typename S>
		class is_source {
			template<typename T>
			static std::true_type test(decltype(T::number_of_dimensions)*);
			template<typename>
			static std::false_type test(...);
		public:
			// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=6709
			typedef decltype(test<S>(0)) result_type;
			static const bool value = result_type::value;
		};

		/**
		 * Nullary expressions
		 */

		template<int, typename>
		class NullaryAcc;

		template<int ndim, typename Func>
		class NullaryExp {
			friend class NullaryAcc<ndim,Func>;
			const Domain<ndim>& dom;
			const Func f;
		public:
			static const int number_of_dimensions = ndim;
			typedef NullaryAcc<ndim,Func> accessor;
			NullaryExp(const Domain<ndim>& dom, const Func& f);
			~NullaryExp();
			INLINE const Domain<ndim>& domain() const;
		};

		template<int ndim, typename Func>
		NullaryExp<ndim,Func>::NullaryExp(const Domain<ndim>& dom, const Func& f): dom(dom), f(f) { }

		template<int ndim, typename Func>
		NullaryExp<ndim,Func>::~NullaryExp() { }

		template<int ndim, typename Func>
		inline const Domain<ndim>& NullaryExp<ndim,Func>::domain() const {
			return dom;
		}

		template<int ndim, typename Func>
		class NullaryAcc {
			const Func& f;
		public:
			NullaryAcc(const NullaryExp<ndim,Func>& exp);
			~NullaryAcc();
			INLINE auto operator()(coords<ndim> ii) const -> decltype(f(ii));
		};

		template<int ndim, typename Func>
		NullaryAcc<ndim,Func>::NullaryAcc(const NullaryExp<ndim,Func>& exp): f(exp.f) { }

		template<int ndim, typename Func>
		NullaryAcc<ndim,Func>::~NullaryAcc() { }

		template<int ndim, typename Func>
		inline auto NullaryAcc<ndim,Func>::operator()(coords<ndim> ii) const -> decltype(f(ii)) {
			return f(ii);
		}

		/**
		 * Unary expressions
		 */

		template<typename, typename>
		class UnaryAcc;

		template<typename S, typename Func>
		class UnaryExp {
			friend class UnaryAcc<S,Func>;
			const S& src;
			const Func f;
		public:
			static const int number_of_dimensions = S::number_of_dimensions;
			typedef UnaryAcc<S,Func> accessor;
			UnaryExp(const S& src, const Func& f);
			~UnaryExp();
			INLINE const Domain<number_of_dimensions>& domain() const;
		};

		template<typename S, typename Func>
		UnaryExp<S,Func>::UnaryExp(const S& src, const Func& f): src(src), f(f) { }

		template<typename S, typename Func>
		UnaryExp<S,Func>::~UnaryExp() { }

		template<typename S, typename Func>
		inline const Domain<UnaryExp<S,Func>::number_of_dimensions>& UnaryExp<S,Func>::domain() const {
			return src.domain();
		}

		template<typename S, typename Func>
		class UnaryAcc {
			const typename S::accessor a;
			const Func& f;
		public:
			UnaryAcc(const UnaryExp<S,Func>& exp);
			~UnaryAcc();
			INLINE auto operator()(coords<S::number_of_dimensions> ii) const -> decltype(f(a, ii));
		};

		template<typename S, typename Func>
		UnaryAcc<S,Func>::UnaryAcc(const UnaryExp<S,Func>& exp): a(exp.src), f(exp.f) { }

		template<typename S, typename Func>
		UnaryAcc<S,Func>::~UnaryAcc() { }

		template<typename S, typename Func>
		inline auto UnaryAcc<S,Func>::operator()(coords<S::number_of_dimensions> ii) const -> decltype(f(a, ii)) {
			return f(a, ii);
		}

		/**
		 * Binary expressions
		 */

		template<typename, typename, typename>
		class BinaryAcc;

		template<typename S1, typename S2, typename Func>
		class BinaryExp {
			friend class BinaryAcc<S1,S2,Func>;
			const S1& src1;
			const S2& src2;
			const Func f;
		public:
			static_assert(S1::number_of_dimensions == S2::number_of_dimensions, "source dimensionality");
			static const int number_of_dimensions = S1::number_of_dimensions;
			typedef BinaryAcc<S1,S2,Func> accessor;
			BinaryExp(const S1& src1, const S2& src2, const Func& f);
			~BinaryExp();
			INLINE const Domain<number_of_dimensions>& domain() const;
		};

		template<typename S1, typename S2, typename Func>
		BinaryExp<S1,S2,Func>::BinaryExp(const S1& src1, const S2& src2, const Func& f): src1(src1), src2(src2), f(f) {
			assert(src1.domain() == src2.domain());
		}

		template<typename S1, typename S2, typename Func>
		BinaryExp<S1,S2,Func>::~BinaryExp() { }

		template<typename S1, typename S2, typename Func>
		inline const Domain<BinaryExp<S1,S2,Func>::number_of_dimensions>& BinaryExp<S1,S2,Func>::domain() const {
			return src1.domain();
		}

		template<typename S1, typename S2, typename Func>
		class BinaryAcc {
			const typename S1::accessor a1;
			const typename S2::accessor a2;
			const Func& f;
		public:
			BinaryAcc(const BinaryExp<S1,S2,Func>& exp);
			~BinaryAcc();
			INLINE auto operator()(coords<S1::number_of_dimensions> ii) const -> decltype(f(a1, a2, ii));
		};

		template<typename S1, typename S2, typename Func>
		BinaryAcc<S1,S2,Func>::BinaryAcc(const BinaryExp<S1,S2,Func>& exp): a1(exp.src1), a2(exp.src2), f(exp.f) { }

		template<typename S1, typename S2, typename Func>
		BinaryAcc<S1,S2,Func>::~BinaryAcc() { }

		template<typename S1, typename S2, typename Func>
		inline auto BinaryAcc<S1,S2,Func>::operator()(coords<S1::number_of_dimensions> ii) const -> decltype(f(a1, a2, ii)) {
			return f(a1, a2, ii);
		}

		/**
		 * Constant value over domain
		 */

		template<int ndim, typename T>
		class Const {
			const T val;
		public:
			Const(const T& val);
			~Const();
			INLINE const T& operator()(coords<ndim> ii) const;
		};

		template<int ndim, typename T>
		inline const T& Const<ndim,T>::operator()(coords<ndim>) const {
			return val;
		}

		template<int ndim, typename T>
		NullaryExp<ndim,Const<ndim,T>> constant(const Domain<ndim>& dom, const T& val);

		/**
		 * Negation
		 */

		template<typename S>
		class Neg {
		public:
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(-a(ii)) {
				return -a(ii);
			}
		};

		template<typename S>
		typename std::enable_if<is_source<S>::value, UnaryExp<S,Neg<S>>>::type operator-(const S& src) {
			return UnaryExp<S,Neg<S>>(src, Neg<S>());
		}

		/**
		 * Scalar multiplication
		 */

		template<typename S, typename T>
		class Scale {
		private:
			const T val;
		public:
			Scale(const T& val): val(val) { }
			~Scale() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii)*val) {
				return a(ii) * val;
			}
		};

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,Scale<S,T>>>::type operator*(const S& src, const T& val) {
			return UnaryExp<S,Scale<S,T>>(src, Scale<S,T>(val));
		}

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,Scale<S,T>>>::type operator*(const T& val, const S& src) {
			return UnaryExp<S,Scale<S,T>>(src, Scale<S,T>(val));
		}

		/**
		 * Addition
		 */
		template<typename S1, typename S2>
		class Add {
		public:
			auto operator()(const typename S1::accessor& a1, const typename S2::accessor& a2, coords<S1::number_of_dimensions> ii) const -> decltype(a1(ii) + a2(ii)) {
				return a1(ii) + a2(ii);
			}
		};

		template<typename S1, typename S2>
		typename std::enable_if<is_source<S1>::value && is_source<S2>::value, BinaryExp<S1,S2,Add<S1,S2>>>::type operator+(const S1& src1, const S2& src2) {
			return BinaryExp<S1,S2,Add<S1,S2>>(src1, src2, Add<S1,S2>());
		}

		/**
		 * Subtraction
		 */
		template<typename S1, typename S2>
		class Sub {
		public:
			auto operator()(const typename S1::accessor& a1, const typename S2::accessor& a2, coords<S1::number_of_dimensions> ii) const -> decltype(a1(ii) - a2(ii)) {
				return a1(ii) - a2(ii);
			}
		};

		template<typename S1, typename S2>
		typename std::enable_if<is_source<S1>::value && is_source<S2>::value, BinaryExp<S1,S2,Sub<S1,S2>>>::type operator-(const S1& src1, const S2& src2) {
			return BinaryExp<S1,S2,Sub<S1,S2>>(src1, src2, Sub<S1,S2>());
		}

		/**
		 * Multiplication
		 */
		template<typename S1, typename S2>
		class Mul {
		public:
			auto operator()(const typename S1::accessor& a1, const typename S2::accessor& a2, coords<S1::number_of_dimensions> ii) const -> decltype(a1(ii) * a2(ii)) {
				return a1(ii) * a2(ii);
			}
		};

		template<typename S1, typename S2>
		typename std::enable_if<is_source<S1>::value && is_source<S2>::value, BinaryExp<S1,S2,Mul<S1,S2>>>::type operator*(const S1& src1, const S2& src2) {
			return BinaryExp<S1,S2,Mul<S1,S2>>(src1, src2, Mul<S1,S2>());
		}
	}

}

#endif
