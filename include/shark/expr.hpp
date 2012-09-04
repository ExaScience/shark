
#ifndef __SHARK_EXPR_HPP
#define __SHARK_EXPR_HPP

#include <cassert>
#include <utility>         // std::declval
#include <type_traits>     // std::true_type, std::enable_if
// Temporary fix to enable declval for libstdc++ 4.4
#if defined(__GLIBCXX__) && __GLIBCXX__ < 20100414
#include <type_traits>
namespace std {
	template<typename T>
	typename add_rvalue_reference<T>::type declval();
}
#endif
#include "abs.hpp"
#include "common.hpp"
#include "coords.hpp"
#include "vec.hpp"

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

		template<typename S>
		struct source {
			typedef typename std::remove_reference<
				//typename std::result_of<typename S::accessor(coords<S::number_of_dimensions>)>::type
				decltype(std::declval<typename S::accessor>()(std::declval<coords<S::number_of_dimensions>>()))
			>::type element_type;
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
			const coords_range<ndim> r;
			const Func f;
		public:
			static const int number_of_dimensions = ndim;
			typedef NullaryExp<ndim,Func> storage;
			typedef NullaryAcc<ndim,Func> accessor;
			NullaryExp(const Domain<ndim>& dom, const Func& f);
			NullaryExp(const Domain<ndim>& dom, coords_range<ndim> r, const Func& f);
			~NullaryExp();
			INLINE const Domain<ndim>& domain() const;
			INLINE coords_range<ndim> region() const;
		};

		template<int ndim, typename Func>
		NullaryExp<ndim,Func>::NullaryExp(const Domain<ndim>& dom, const Func& f): dom(dom), r(dom.total()), f(f) { }

		template<int ndim, typename Func>
		NullaryExp<ndim,Func>::NullaryExp(const Domain<ndim>& dom, coords_range<ndim> r, const Func& f): dom(dom), r(r), f(f) { }

		template<int ndim, typename Func>
		NullaryExp<ndim,Func>::~NullaryExp() { }

		template<int ndim, typename Func>
		inline const Domain<ndim>& NullaryExp<ndim,Func>::domain() const {
			return dom;
		}

		template<int ndim, typename Func>
		inline coords_range<ndim> NullaryExp<ndim,Func>::region() const {
			return r;
		}

		template<int ndim, typename Func>
		class NullaryAcc {
			const Func& f;
		public:
			NullaryAcc(const NullaryExp<ndim,Func>& exp);
			~NullaryAcc();
			// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=54359
			//INLINE auto operator()(coords<ndim> ii) const -> decltype(f(ii));
			//INLINE typename std::result_of<Func(coords<ndim>)>::type
			INLINE decltype(std::declval<Func>()(std::declval<coords<ndim>>()))
			operator()(coords<ndim> ii) const;
		};

		template<int ndim, typename Func>
		NullaryAcc<ndim,Func>::NullaryAcc(const NullaryExp<ndim,Func>& exp): f(exp.f) { }

		template<int ndim, typename Func>
		NullaryAcc<ndim,Func>::~NullaryAcc() { }

		template<int ndim, typename Func>
		//inline typename std::result_of<Func(coords<ndim>)>::type 
		inline decltype(std::declval<Func>()(std::declval<coords<ndim>>()))
		NullaryAcc<ndim,Func>::operator()(coords<ndim> ii) const {
			return f(ii);
		}

		template<int ndim, typename Func>
		NullaryExp<ndim,Func> nullary(const Domain<ndim>& dom, const Func& f) {
			return NullaryExp<ndim,Func>(dom,f);
		}

		/**
		 * Unary expressions
		 */

		template<typename, typename>
		class UnaryAcc;

		template<typename S, typename Func>
		class UnaryExp {
			friend class UnaryAcc<S,Func>;
			const typename S::storage src;
			const Func f;
		public:
			static const int number_of_dimensions = S::number_of_dimensions;
			typedef UnaryExp<S,Func> storage;
			typedef UnaryAcc<S,Func> accessor;
			UnaryExp(const S& src, const Func& f);
			~UnaryExp();
			// http://llvm.org/bugs/show_bug.cgi?id=13729
			//INLINE const Domain<number_of_dimensions>& domain() const;
			//INLINE coords_range<number_of_dimensions> region() const;
			INLINE const Domain<S::number_of_dimensions>& domain() const;
			INLINE coords_range<S::number_of_dimensions> region() const;
		};

		template<typename S, typename Func>
		UnaryExp<S,Func>::UnaryExp(const S& src, const Func& f): src(src), f(f) { }

		template<typename S, typename Func>
		UnaryExp<S,Func>::~UnaryExp() { }

		template<typename S, typename Func>
		inline const Domain<S::number_of_dimensions>& UnaryExp<S,Func>::domain() const {
			return src.domain();
		}

		template<typename S, typename Func>
		inline coords_range<S::number_of_dimensions> UnaryExp<S,Func>::region() const {
			return src.region();
		}

		template<typename S, typename Func>
		class UnaryAcc {
			const typename S::accessor a;
			const Func& f;
		public:
			UnaryAcc(const UnaryExp<S,Func>& exp);
			~UnaryAcc();
			// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=54359
			//INLINE auto operator()(coords<S::number_of_dimensions> ii) const -> decltype(f(a, ii));
			//INLINE typename std::result_of<Func(typename S::accessor, coords<S::number_of_dimensions>)>::type
			INLINE decltype(std::declval<Func>()(std::declval<typename S::accessor>(), std::declval<coords<S::number_of_dimensions>>()))
			operator()(coords<S::number_of_dimensions> ii) const;
		};

		template<typename S, typename Func>
		UnaryAcc<S,Func>::UnaryAcc(const UnaryExp<S,Func>& exp): a(exp.src), f(exp.f) { }

		template<typename S, typename Func>
		UnaryAcc<S,Func>::~UnaryAcc() { }

		template<typename S, typename Func>
		//inline typename std::result_of<Func(typename S::accessor, coords<S::number_of_dimensions>)>::type
		inline decltype(std::declval<Func>()(std::declval<typename S::accessor>(), std::declval<coords<S::number_of_dimensions>>()))
		UnaryAcc<S,Func>::operator()(coords<S::number_of_dimensions> ii) const {
			return f(a, ii);
		}

		template<typename S, typename Func>
		UnaryExp<S,Func> unary(const S& src, const Func& f) {
			return UnaryExp<S,Func>(src, f);
		}

		/**
		 * Binary expressions
		 */

		template<typename, typename, typename>
		class BinaryAcc;

		template<typename S1, typename S2, typename Func>
		class BinaryExp {
			friend class BinaryAcc<S1,S2,Func>;
			const typename S1::storage src1;
			const typename S2::storage src2;
			const Func f;
		public:
			static_assert(S1::number_of_dimensions == S2::number_of_dimensions, "source dimensionality");
			static const int number_of_dimensions = S1::number_of_dimensions;
			typedef BinaryExp<S1,S2,Func> storage;
			typedef BinaryAcc<S1,S2,Func> accessor;
			BinaryExp(const S1& src1, const S2& src2, const Func& f);
			~BinaryExp();
			// http://llvm.org/bugs/show_bug.cgi?id=13729
			//INLINE const Domain<number_of_dimensions>& domain() const;
			//INLINE coords_range<number_of_dimensions> region() const;
			INLINE const Domain<S1::number_of_dimensions>& domain() const;
			INLINE coords_range<S1::number_of_dimensions> region() const;
		};

		template<typename S1, typename S2, typename Func>
		BinaryExp<S1,S2,Func>::BinaryExp(const S1& src1, const S2& src2, const Func& f): src1(src1), src2(src2), f(f) {
			assert(src1.domain() == src2.domain());
		}

		template<typename S1, typename S2, typename Func>
		BinaryExp<S1,S2,Func>::~BinaryExp() { }

		template<typename S1, typename S2, typename Func>
		inline const Domain<S1::number_of_dimensions>& BinaryExp<S1,S2,Func>::domain() const {
			return src1.domain();
		}

		template<typename S1, typename S2, typename Func>
		inline coords_range<S1::number_of_dimensions> BinaryExp<S1,S2,Func>::region() const {
			return src1.region().overlap(src2.region());
		}

		template<typename S1, typename S2, typename Func>
		class BinaryAcc {
			const typename S1::accessor a1;
			const typename S2::accessor a2;
			const Func& f;
		public:
			BinaryAcc(const BinaryExp<S1,S2,Func>& exp);
			~BinaryAcc();
			// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=54359
			// INLINE auto operator()(coords<S1::number_of_dimensions> ii) const -> decltype(f(a1, a2, ii));
			// INLINE typename std::result_of<Func(typename S1::accessor, typename S2::accessor, coords<S1::number_of_dimensions>)>::type
			INLINE decltype(std::declval<Func>()(std::declval<typename S1::accessor>(), std::declval<typename S2::accessor>(), std::declval<coords<S1::number_of_dimensions>>()))
			operator()(coords<S1::number_of_dimensions> ii) const;
		};

		template<typename S1, typename S2, typename Func>
		BinaryAcc<S1,S2,Func>::BinaryAcc(const BinaryExp<S1,S2,Func>& exp): a1(exp.src1), a2(exp.src2), f(exp.f) { }

		template<typename S1, typename S2, typename Func>
		BinaryAcc<S1,S2,Func>::~BinaryAcc() { }

		template<typename S1, typename S2, typename Func>
		//inline typename std::result_of<Func(typename S1::accessor, typename S2::accessor, coords<S1::number_of_dimensions>)>::type
		inline decltype(std::declval<Func>()(std::declval<typename S1::accessor>(), std::declval<typename S2::accessor>(), std::declval<coords<S1::number_of_dimensions>>()))
		BinaryAcc<S1,S2,Func>::operator()(coords<S1::number_of_dimensions> ii) const {
			return f(a1, a2, ii);
		}

		template<typename S1, typename S2, typename Func>
		BinaryExp<S1,S2,Func> binary(const S1& src1, const S2& src2, const Func& f) {
			return BinaryExp<S1,S2,Func>(src1, src2, f);
		}

		/**
		 * Group sources as a vector
		 */
		template<typename...>
		class AsVecAcc;

		template<typename...>
		class AsVecExp;

		template<typename S>
		class AsVecExp<S> {
			friend class AsVecAcc<S>;
			const typename S::storage src;
		public:
			static const int number_of_dimensions = S::number_of_dimensions;
			typedef AsVecExp<S> storage;
			typedef AsVecAcc<S> accessor;
			AsVecExp(const S& s);
			~AsVecExp();
			// http://llvm.org/bugs/show_bug.cgi?id=13729
			//INLINE const Domain<number_of_dimensions>& domain() const;
			//INLINE coords_range<number_of_dimensions> region() const;
			INLINE const Domain<S::number_of_dimensions>& domain() const;
			INLINE coords_range<S::number_of_dimensions> region() const;
		};

		template<typename S>
		AsVecExp<S>::AsVecExp(const S& src): src(src) { }

		template<typename S>
		AsVecExp<S>::~AsVecExp() { }

		template<typename S>
		inline const Domain<S::number_of_dimensions>& AsVecExp<S>::domain() const {
			return src.domain();
		}

		template<typename S>
		inline coords_range<S::number_of_dimensions> AsVecExp<S>::region() const {
			return src.region();
		}

		template<typename S, typename... Ss>
		class AsVecExp<S,Ss...> {
			friend class AsVecAcc<S,Ss...>;
			const typename S::storage src;
			const AsVecExp<Ss...> srcs;
		public:
			static_assert(S::number_of_dimensions == AsVecExp<Ss...>::number_of_dimensions, "source dimensionality");
			static const int number_of_dimensions = S::number_of_dimensions;
			typedef AsVecExp<S,Ss...> storage;
			typedef AsVecAcc<S,Ss...> accessor;
			AsVecExp(const S& src, const Ss&... srcs);
			~AsVecExp();
			// http://llvm.org/bugs/show_bug.cgi?id=13729
			//INLINE const Domain<number_of_dimensions>& domain() const;
			//INLINE coords_range<number_of_dimensions> region() const;
			INLINE const Domain<S::number_of_dimensions>& domain() const;
			INLINE coords_range<S::number_of_dimensions> region() const;
		};

		template<typename S, typename... Ss>
		AsVecExp<S,Ss...>::AsVecExp(const S& src, const Ss&... srcs): src(src), srcs(srcs...) {
			assert(src.domain() == this->srcs.domain());
		}

		template<typename S, typename... Ss>
		AsVecExp<S,Ss...>::~AsVecExp() { }

		template<typename S, typename... Ss>
		inline const Domain<S::number_of_dimensions>& AsVecExp<S,Ss...>::domain() const {
			return src.domain();
		}

		template<typename S, typename... Ss>
		inline coords_range<S::number_of_dimensions> AsVecExp<S,Ss...>::region() const {
			return src.region().overlap(srcs.region());
		}

		template<typename... Ss>
		AsVecExp<Ss...> as_vec(const Ss&... ss) {
			return AsVecExp<Ss...>(ss...);
		}

		template<typename S>
		class AsVecAcc<S> {
			const typename S::accessor a;
		public:
			typedef vec<1, typename source<S>::element_type> vec_type;
			AsVecAcc(const AsVecExp<S>& exp);
			~AsVecAcc();
			template<int d, typename T> INLINE void fill(T& t, coords<S::number_of_dimensions> ii) const;
			INLINE vec_type operator()(coords<S::number_of_dimensions> ii) const;
		};

		template<typename S>
		AsVecAcc<S>::AsVecAcc(const AsVecExp<S>& exp): a(exp.src) { }

		template<typename S>
		AsVecAcc<S>::~AsVecAcc() { }

		template<typename S> template<int d, typename T>
		inline void AsVecAcc<S>::fill(T& t, coords<S::number_of_dimensions> ii) const {
			t[d] = a(ii);
		}

		template<typename S>
		inline typename AsVecAcc<S>::vec_type AsVecAcc<S>::operator()(coords<S::number_of_dimensions> ii) const {
			vec_type v;
			fill<0>(v, ii);
			return v;
		}

		template<typename S, typename... Ss>
		class AsVecAcc<S,Ss...> {
			const typename S::accessor a;
			const AsVecAcc<Ss...> aa;
		public:
			typedef vec<1+sizeof...(Ss), typename std::common_type<typename source<S>::element_type, typename source<Ss>::element_type...>::type> vec_type;
			AsVecAcc(const AsVecExp<S,Ss...>& exp);
			~AsVecAcc();
			template<int d, typename T> INLINE void fill(T& t, coords<S::number_of_dimensions> ii) const;
			INLINE vec_type operator()(coords<S::number_of_dimensions> ii) const;
		};

		template<typename S, typename... Ss>
		AsVecAcc<S,Ss...>::AsVecAcc(const AsVecExp<S,Ss...>& exp): a(exp.src), aa(exp.srcs) { }

		template<typename S, typename... Ss>
		AsVecAcc<S,Ss...>::~AsVecAcc() { }

		template<typename S, typename... Ss> template<int d, typename T>
		inline void AsVecAcc<S,Ss...>::fill(T& t, coords<S::number_of_dimensions> ii) const {
			t[d] = a(ii);
			aa.template fill<d+1>(t, ii);
		}
		
		template<typename S, typename... Ss>
		inline typename AsVecAcc<S,Ss...>::vec_type AsVecAcc<S,Ss...>::operator()(coords<S::number_of_dimensions> ii) const {
			vec_type v;
			fill<0>(v, ii);
			return v;
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

		template<int ndim, typename T>
		NullaryExp<ndim,Const<ndim,T>> constant(const Domain<ndim>& dom, coords_range<ndim> r, const T& val);

		/**
		 * Coordinate value over domain
		 */
		template<int ndim>
		class CoordVec {
		public:
			INLINE vec<ndim,coord> operator()(coords<ndim> ii) const;
		};

		template<int ndim>
		inline vec<ndim,coord> CoordVec<ndim>::operator()(coords<ndim> ii) const {
			return ii.to_vec();
		}

		template<int ndim>
		NullaryExp<ndim,CoordVec<ndim>> coord_vec(const Domain<ndim>& dom);

		template<int ndim, typename T>
		class SCoordVec {
			const coords<ndim> lower;
			const vec<ndim,T> step;
		public:
			SCoordVec(coords<ndim> lower, vec<ndim,T> step);
			~SCoordVec();
			INLINE vec<ndim,T> operator()(coords<ndim> ii) const;
		};

		template<int ndim, typename T>
		inline vec<ndim,T> SCoordVec<ndim,T>::operator()(coords<ndim> ii) const {
			return (ii - lower).to_vec() * step;
		}

		template<int ndim, typename T>
		SCoordVec<ndim,T>::SCoordVec(coords<ndim> lower, vec<ndim,T> step): lower(lower), step(step) { }

		template<int ndim, typename T>
		SCoordVec<ndim,T>::~SCoordVec() { }

		template<int ndim, typename T>
		NullaryExp<ndim,SCoordVec<ndim,T>> coord_vec(const Domain<ndim>& dom, coords_range<ndim> r, vec<ndim,T> one) {
			return NullaryExp<ndim,SCoordVec<ndim,T>>(dom, SCoordVec<ndim,T>(r.lower, one / (r.upper - r.lower - coords<ndim>::one()).to_vec()));
		}

		template<int ndim, typename T>
		NullaryExp<ndim,SCoordVec<ndim,T>> coord_vec(const Domain<ndim>& dom, vec<ndim,T> one) {
			return coord_vec(dom, dom.total(), one);
		}

		template<int d, int ndim>
		class CoordVal {
		public:
			static_assert(d >= 0 && d < ndim, "0 <= d < ndim");
			INLINE coord operator()(coords<ndim> ii) const;
		};

		template<int d, int ndim>
		inline coord CoordVal<d,ndim>::operator()(coords<ndim> ii) const {
			return ii[d];
		}

		template<int d, int ndim>
		NullaryExp<ndim,CoordVal<d,ndim>> coord_val(const Domain<ndim>& dom) {
			return NullaryExp<ndim,CoordVal<d,ndim>>(dom,CoordVal<d,ndim>());
		}

		template<int d, int ndim, typename T>
		class SCoordVal {
			const coord lower;
			const T step;
		public:
			static_assert(d >= 0 && d < ndim, "0 <= d < ndim");
			SCoordVal(coord lower, T step);
			~SCoordVal();
			INLINE T operator()(coords<ndim> ii) const;
		};

		template<int d, int ndim, typename T>
		inline T SCoordVal<d,ndim,T>::operator()(coords<ndim> ii) const {
			return (ii[d] - lower) * step;
		}

		template<int d, int ndim, typename T>
		SCoordVal<d,ndim,T>::SCoordVal(coord lower, T step): lower(lower), step(step) { }

		template<int d, int ndim, typename T>
		SCoordVal<d,ndim,T>::~SCoordVal() { }

		template<int d, int ndim, typename T>
		NullaryExp<ndim,SCoordVal<d,ndim,T>> coord_val(const Domain<ndim>& dom, coord lower, coord upper, T one) {
			return NullaryExp<ndim,SCoordVal<d,ndim,T>>(dom, SCoordVal<d,ndim,T>(lower, one / (upper - lower - 1)));
		}

		template<int d, int ndim, typename T>
		NullaryExp<ndim,SCoordVal<d,ndim,T>> coord_val(const Domain<ndim>& dom, T one) {
			return coord_val<d>(dom, coord(0), dom.n[d], one);
		}

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
		 * Absolute value
		 */

		template<typename S>
		class Abs {
		public:
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(abs(a(ii))) {
				return abs(a(ii));
			}
		};

		template<typename S>
		typename std::enable_if<is_source<S>::value, UnaryExp<S,Abs<S>>>::type abs(const S& src) {
			return UnaryExp<S,Abs<S>>(src, Abs<S>());
		}

		/**
		 * Minimal/maximal element
		 */
		template<typename S>
		class MinElement {
		public:
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii).min()) {
				return a(ii).min();
			}
		};

		template<typename S>
		class MaxElement {
		public:
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii).max()) {
				return a(ii).max();
			}
		};

		template<typename S>
		typename std::enable_if<is_source<S>::value, UnaryExp<S,MinElement<S>>>::type min_element(const S& src) {
			return UnaryExp<S,MinElement<S>>(src, MinElement<S>());
		}

		template<typename S>
		typename std::enable_if<is_source<S>::value, UnaryExp<S,MaxElement<S>>>::type max_element(const S& src) {
			return UnaryExp<S,MaxElement<S>>(src, MaxElement<S>());
		}

		/**
		 * Addition
		 */

		template<typename S, typename T>
		class AddC {
		private:
			const T val;
		public:
			AddC(const T& val): val(val) { }
			~AddC() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii)+val) {
				return a(ii) + val;
			}
		};

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,AddC<S,T>>>::type operator+(const S& src, const T& val) {
			return UnaryExp<S,AddC<S,T>>(src, AddC<S,T>(val));
		}

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,AddC<S,T>>>::type operator+(const T& val, const S& src) {
			return UnaryExp<S,AddC<S,T>>(src, AddC<S,T>(val));
		}

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

		template<typename S, typename T>
		class SubCL {
		private:
			const T val;
		public:
			SubCL(const T& val): val(val) { }
			~SubCL() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(val-a(ii)) {
				return val - a(ii);
			}
		};

		template<typename S, typename T>
		class SubCR {
		private:
			const T val;
		public:
			SubCR(const T& val): val(val) { }
			~SubCR() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii)-val) {
				return a(ii) - val;
			}
		};

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,SubCL<S,T>>>::type operator-(const T& val, const S& src) {
			return UnaryExp<S,SubCL<S,T>>(src, SubCL<S,T>(val));
		}

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,SubCR<S,T>>>::type operator-(const S& src, const T& val) {
			return UnaryExp<S,SubCR<S,T>>(src, SubCR<S,T>(val));
		}

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
		template<typename S, typename T>
		class MulC {
		private:
			const T val;
		public:
			MulC(const T& val): val(val) { }
			~MulC() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii)*val) {
				return a(ii) * val;
			}
		};

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,MulC<S,T>>>::type operator*(const S& src, const T& val) {
			return UnaryExp<S,MulC<S,T>>(src, MulC<S,T>(val));
		}

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,MulC<S,T>>>::type operator*(const T& val, const S& src) {
			return UnaryExp<S,MulC<S,T>>(src, MulC<S,T>(val));
		}

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

		/**
		 * Division
		 */

		template<typename S, typename T>
		class DivCL {
		private:
			const T val;
		public:
			DivCL(const T& val): val(val) { }
			~DivCL() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(val/a(ii)) {
				return val / a(ii);
			}
		};

		template<typename S, typename T>
		class DivCR {
		private:
			const T val;
		public:
			DivCR(const T& val): val(val) { }
			~DivCR() { }
			auto operator()(const typename S::accessor& a, coords<S::number_of_dimensions> ii) const -> decltype(a(ii)/val) {
				return a(ii) / val;
			}
		};

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,DivCL<S,T>>>::type operator/(const T& val, const S& src) {
			return UnaryExp<S,DivCL<S,T>>(src, DivCL<S,T>(val));
		}

		template<typename S, typename T>
		typename std::enable_if<is_source<S>::value && !is_source<T>::value, UnaryExp<S,DivCR<S,T>>>::type operator/(const S& src, const T& val) {
			return UnaryExp<S,DivCR<S,T>>(src, DivCR<S,T>(val));
		}

		template<typename S1, typename S2>
		class Div {
		public:
			auto operator()(const typename S1::accessor& a1, const typename S2::accessor& a2, coords<S1::number_of_dimensions> ii) const -> decltype(a1(ii) / a2(ii)) {
				return a1(ii) / a2(ii);
			}
		};

		template<typename S1, typename S2>
		typename std::enable_if<is_source<S1>::value && is_source<S2>::value, BinaryExp<S1,S2,Div<S1,S2>>>::type operator/(const S1& src1, const S2& src2) {
			return BinaryExp<S1,S2,Div<S1,S2>>(src1, src2, Div<S1,S2>());
		}

		/**
		 * Comparisons
		 */
		template<typename S1, typename S2>
		class Eq {
		public:
			auto operator()(const typename S1::accessor& a1, const typename S2::accessor& a2, coords<S1::number_of_dimensions> ii) const -> decltype(a1(ii) == a2(ii)) {
				return a1(ii) == a2(ii);
			}
		};

		template<typename S1, typename S2>
		typename std::enable_if<is_source<S1>::value && is_source<S2>::value, BinaryExp<S1,S2,Eq<S1,S2>>>::type operator==(const S1& src1, const S2& src2) {
			return BinaryExp<S1,S2,Eq<S1,S2>>(src1, src2, Eq<S1,S2>());
		}

		/**
		 * Sum reduction
		 */
		template<typename R, typename S>
		typename std::enable_if<is_source<S>::value, R>::type sum(const R& zero, const S& src) {
			const typename S::accessor s(src);
			return src.domain().sum(src.region(), zero, [&s](R& acc, coords<S::number_of_dimensions> i) {
				acc += s(i);
			});
		}

		/**
		 * Matrix norms
		 */
		template<typename S>
		typename std::enable_if<is_source<S>::value, typename source<S>::element_type>::type norm1(const S& src) {
			return sum(typename source<S>::element_type(), abs(src));
		}

		template<typename S>
		typename std::enable_if<is_source<S>::value, typename source<S>::element_type>::type norm2(const S& src) {
			return sqrt(sum(typename source<S>::element_type(), src * src));
		}

		/**
		 * Dot product
		 */
		template<typename S1, typename S2>
		typename std::enable_if<is_source<S1>::value && is_source<S2>::value, typename source<S1>::element_type>::type dot(const S1& src1, const S2& src2) {
			return sum(typename source<S1>::element_type(), src1 * src2);
		}

	}

}

#endif
