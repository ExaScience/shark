/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <utility>    // std::move
#include <string>     // std::string
#include <iostream>   // std::ostream
#include <memory>     // std::unique_ptr
#include <shark.hpp>

using namespace std;
using namespace shark;

class tester {
	ostream& out;
	string current;
	long fails;
	long tests;
	long curfails;
	long curchecks;

public:
	tester(ostream& out): out(out), fails(0), tests(0) { }
	~tester() { }
	void begin_test(const std::string& name);
	void end_test();
	void add_result(test_result r);
};

void tester::begin_test(const std::string& name) {
	current = name;
	curfails = 0;
	curchecks = 0;
	if(world().procid == 0)
		out << "Begin test " << current << endl;
	world().sync();
}

void tester::end_test() {
	world().sync();
	bool fail = curfails > 0;
	tests++;
	if(fail) fails++;
	if(world().procid == 0) {
		out << "End test " << current << ": ";
		if(!fail)
			out << "Success (" << curchecks << " checks)" << endl;
		else
			out << "FAILED (" << curfails << " of " << curchecks << " checks failed)" << endl;
	}
}

void tester::add_result(test_result r) {
	curchecks += r.checks;
	curfails += r.fails;
}

using namespace shark::ndim;

template<int ndim, typename S>
class suite1 {
	static_assert(S::number_of_dimensions == ndim, "Source dimensionality");

	const Domain<ndim>& dom;
	const S src;
	typedef typename source<S>::element_type T;

	coords_range<ndim> subrange();
	unique_ptr<T[]> src_buf(coords_range<ndim> r);

	void test_reshape_domain(tester& t, const Domain<ndim>& alt_dom);

public:
	suite1(const Domain<ndim>& dom, const S& src): dom(dom), src(src) { }
	~suite1() { }
	void test_basic(tester& t);
	void test_move(tester& t);
	void test_ghost(tester& t);
	void test_ghost_corner(tester& t);
	void test_ghost_periodic(tester& t);
	void test_get(tester& t);
	void test_put(tester& t);
	void test_accumulate(tester& t);
	void test_gather(tester& t);
	void test_reshape(tester& t);
	void run(tester& t) {
		test_basic(t);
		test_move(t);
		test_ghost(t);
		test_ghost_corner(t);
		test_ghost_periodic(t);
		test_get(t);
		test_put(t);
		test_accumulate(t);
		test_gather(t);
		test_reshape(t);
	}
};

template<int ndim, typename S>
suite1<ndim,S> make_suite1(const Domain<ndim>& dom, const S& src) {
	return suite1<ndim,S>(dom, src);
}

template<int ndim, typename S>
coords_range<ndim> suite1<ndim,S>::subrange() {
	coords_range<ndim> total = dom.total();
	coords_range<ndim> r;
	for(int d = 0; d < ndim; d++) {
		coord q = (total.upper[d] - total.lower[d])/4;
		r.lower[d] = total.lower[d] + q;
		r.upper[d] = total.lower[d] + q * 3;
	}
	return r;
}

template<int ndim, typename S>
unique_ptr<typename suite1<ndim,S>::T[]> suite1<ndim,S>::src_buf(coords_range<ndim> r) {
	const typename S::accessor s(src);
	coords<ndim+1> ld = r.stride();
	unique_ptr<T[]> ptr(new T[ld[0]]);
	r.for_each([&s,r,ld,&ptr](coords<ndim> ii) {
		ptr[(ii - r.lower).offset(ld)] = s(ii);
	});
	return ptr;
}

template<int ndim, typename S>
void suite1<ndim,S>::test_basic(tester& t) {
	t.begin_test("test_basic");
	{
		GlobalArray<ndim,T> ga(dom);
		ga = src;
		t.add_result(check(ga == src));
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_move(tester& t) {
	t.begin_test("test_move");
	{
		GlobalArray<ndim,T> ga(dom);
		ga = constant(dom, T());
		{
			GlobalArray<ndim,T> tmp(dom);
			tmp = src;
			ga = std::move(tmp);
		}
		t.add_result(check(ga == src));
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_ghost(tester& t) {
	t.begin_test("test_ghost");
	coords<ndim> gw;
	for(int d = 0; d < ndim; d++)
		gw[d] = 2;
	coords_range<ndim> inner = dom.total();
	for(int d = 0; d < ndim; d++) {
		inner.lower[d] += gw[d];
		inner.upper[d] -= gw[d];
	}
	{
		GlobalArray<ndim,T> ga(dom, gw);
		ga = src;
		ga.update();
		t.add_result(check(inner, ga == src));
		for(int d = 0; d < ndim; d++) {
			coords<ndim> sd = {{}}, su = {{}};
			sd[d] = -1;
			su[d] = 1;
			t.add_result(check(inner,
				shift(ga,sd) == shift(src,sd) &&
				shift(ga,su) == shift(src,su) &&
				shift(ga,sd+sd) == shift(src,sd+sd) &&
				shift(ga,su+su) == shift(src,su+su)));
		}
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_ghost_corner(tester& t) {
	t.begin_test("test_ghost_corner");
	coords<ndim> gw;
	for(int d = 0; d < ndim; d++)
		gw[d] = 2;
	coords_range<ndim> inner = dom.total();
	for(int d = 0; d < ndim; d++) {
		inner.lower[d] += gw[d];
		inner.upper[d] -= gw[d];
	}
	{
		GlobalArray<ndim,T> ga(dom, gw, true);
		ga = src;
		ga.update();
		t.add_result(check(inner, ga == src));
		t.add_result(check(inner, shift(ga,gw) == shift(src,gw) && shift(ga,-gw) == shift(src,-gw)));
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_ghost_periodic(tester& t) {
	t.begin_test("test_ghost_periodic");
	coords<ndim> gw;
	for(int d = 0; d < ndim; d++)
		gw[d] = 2;
	typename GlobalArray<ndim,T>::bounds bd;
	for(int d = 0; d < ndim; d++)
		bd[d] = Boundary<ndim,T>::periodic();
	{
		GlobalArray<ndim,T> ga(dom, gw, false, bd);
		ga = src;
		ga.update();
		t.add_result(check(ga == src));
		for(int d = 0; d < ndim; d++) {
			coords<ndim> sd = {{}}, su = {{}};
			sd[d] = -1;
			su[d] = 1;
			t.add_result(check(
				shift(ga,sd) == shift(src,sd,dom.total()) &&
				shift(ga,su) == shift(src,su,dom.total()) && 
				shift(ga,sd+sd) == shift(src,sd+sd,dom.total()) &&
				shift(ga,su+su) == shift(src,su+su,dom.total())));
		}
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_get(tester& t) {
	t.begin_test("test_get");
	coords_range<ndim> r = subrange();
	coords<ndim+1> ld = r.stride();
	{
		GlobalArray<ndim,T> ga(dom);
		ga = src;
		// Everyone does a "get" and checks the result
		unique_ptr<T[]> ptr(new T[ld[0]]);
		dom.group.sync();
		ga.get(r, ptr.get());
		dom.group.sync();
		{
			test_result tr = test_result();
			const typename S::accessor s(src);
			r.for_each([&tr,&s,r,ld,&ptr](coords<ndim> ii) {
				if(s(ii) != ptr[(ii - r.lower).offset(ld)])
					tr.fails++;
				tr.checks++;
			});
			// Add everything together
			t.add_result(dom.group.external_sum(move(tr)));
		}
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_put(tester& t) {
	t.begin_test("test_put");
	coords_range<ndim> r = subrange();
	{
		GlobalArray<ndim,T> ga(dom);
		ga = constant(dom, T());
		// First one does a put
		dom.group.sync();
		if (dom.group.procid == 0)
			ga.put(r, src_buf(r).get());
		dom.group.sync();
		// Everyone helps check
		t.add_result(check(r, ga == src));
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_accumulate(tester& t) {
	t.begin_test("test_accumulate");
	coords_range<ndim> r = subrange();
	{
		GlobalArray<ndim,T> ga(dom);
		ga = constant(dom, T());
		// First and last one contribute
		dom.group.sync();
		if (dom.group.procid == 0)
			ga.accumulate(r, src_buf(r).get());
		if (dom.group.procid == dom.group.nprocs-1)
			ga.accumulate(r, src_buf(r).get());
		dom.group.sync();
		// Everyone helps check
		t.add_result(check(r, ga == src+src));
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_gather(tester& t) {
	t.begin_test("test_gather");
	coords_range<ndim> r = subrange();
	{
		// Data
		GlobalArray<ndim,T> ga(dom);
		ga = src;
		// Mark local region of interest
		SparseArray<ndim,T> sa(dom, 6);
		r.for_each([&sa](coords<ndim> i) {
			sa.set(i);
		});
		// Retrieve
		ga.gather(sa);
		// Check
		{
			test_result tr = test_result();
			const typename S::accessor s(src);
			r.for_each([&sa,&s,&tr](coords<ndim> ii) {
				if(s(ii) != sa.get(ii))
					tr.fails++;
				tr.checks++;
			});
			// Add everything together
			t.add_result(dom.group.external_sum(move(tr)));
		}
	}
	t.end_test();
}

template<int ndim, typename S>
void suite1<ndim,S>::test_reshape_domain(tester& t, const Domain<ndim>& alt_dom) {
	GlobalArray<ndim,T> ga(dom);
	ga = src;
	ga.reshape(alt_dom);
	const typename GlobalArray<ndim,T>::accessor a(ga);
	const typename S::accessor s(src);
	t.add_result(alt_dom.sum(test_result(), [&a,&s](test_result& tr, coords<ndim> ii) {
		if(s(ii) != a(ii))
			tr.fails++;
		tr.checks++;
	}));
}

template<int ndim, typename S>
void suite1<ndim,S>::test_reshape(tester& t) {
	t.begin_test("test_reshape");
	{
		typename Domain<ndim>::pcoords row_np;
		row_np[0] = dom.group.nprocs;
		for(int d = 1; d < ndim; d++)
			row_np[d] = 1;
		test_reshape_domain(t, Domain<ndim>(dom.group, dom.n, row_np));
		typename Domain<ndim>::pcoords col_np;
		for(int d = 0; d < ndim-1; d++)
			col_np[d] = 1;
		col_np[ndim-1] = dom.group.nprocs;
		test_reshape_domain(t, Domain<ndim>(dom.group, dom.n, col_np));
	}
	t.end_test();
}

string usagestr("Usage: ");

int usage(string msg=string()) {
        if(!msg.empty())
                cerr << msg << endl;
        cerr << usagestr << endl;
        return msg.empty() ? 0 : 1;
}


int main(int argc, char* argv[]) {
	Init(&argc, &argv);

	usagestr += argv[0];
	usagestr += " [-v]";

	for(int k = 1; k < argc; k++) {
		string arg(argv[k]);
		if(arg == "-v") {
			if(world().procid == 0) {
				log_out = &cerr;
				log_mask.set();
			}
		} else if(arg == "-h" || arg == "--help")
			return usage();
		else if(arg[0] == '-')
			return usage("Unknown option: " + arg);
		else
			return usage("Unexpected argument: " + arg);
	}

	SetupThreads();
	tester t(cerr);
	if(world().procid == 0)
		cerr << "Testing <1,int>" << endl;
	{
		coords<1> n = {{1000}};
		Domain<1> dom(world(), n);
		if(world().procid == 0)
			dom.outputDistribution(cerr);
		auto f = coord_val<0>(dom, 1000);
		make_suite1(dom, f).run(t);
	}
	if(world().procid == 0)
		cerr << endl << "Testing <3,double>" << endl;
	{
		coords<3> n = {{100,100,100}};
		Domain<3> dom(world(), n);
		if(world().procid == 0)
			dom.outputDistribution(cerr);
		auto f = sin(coord_val<0>(dom, 2*M_PI, false)) * cos(coord_val<1>(dom, 2*M_PI, false)) * coord_val<2>(dom, 1.0, false);
		make_suite1(dom, f).run(t);
	}
	if(world().procid == 0)
		cerr << endl << "Testing <3,vec<3,double>>" << endl;
	{
		coords<3> n = {{100,100,100}};
		Domain<3> dom(world(), n);
		if(world().procid == 0)
			dom.outputDistribution(cerr);
		const vec<3,double> one = {{1.0, 1.0, 1.0}};
		const vec<3,double> mid = {{0.5, 0.5, 0.5}};
		auto f = abs(coord_vec(dom, one) - mid);
		make_suite1(dom, f).run(t);
	}
	Finalize();
}
