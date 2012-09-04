
#include <utility>    // std::move
#include <string>     // std::string
#include <iostream>   // std::ostream
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
	void begin_test(std::string&& name);
	void end_test();
	void add_result(test_result r);
};

void tester::begin_test(std::string&& name) {
	current = move(name);
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

	const Domain<ndim>& dom;
	const S src;

public:
	suite1(const Domain<ndim>& dom, const S& src): dom(dom), src(src) { }
	~suite1() { }
	void test_basic(tester& t);
	void test_ghost(tester& t);
	void test_ghost_corner(tester& t);
	void run(tester& t) {
		test_basic(t);
		test_ghost(t);
		test_ghost_corner(t);
	}
};

template<int ndim, typename S>
suite1<ndim,S> make_suite1(const Domain<ndim>& dom, const S& src) {
	return suite1<ndim,S>(dom, src);
}

template<int ndim, typename S>
void suite1<ndim,S>::test_basic(tester& t) {
	t.begin_test("test_basic");
	{
		GlobalArray<ndim,double> ga(dom);
		ga = src;
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
	for(int d = 0; d < ndim; d++)
		if(!dom.pd[d]) {
			inner.lower[d] += gw[d];
			inner.upper[d] -= gw[d];
		}
	{
		GlobalArray<ndim,double> ga(dom, gw);
		ga = src;
		ga.update();
		t.add_result(check(inner, ga == src));
		for(int d = 0; d < ndim; d++) {
			coords<ndim> sd = {{}};
			coords<ndim> su = {{}};
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
	for(int d = 0; d < ndim; d++)
		if(!dom.pd[d]) {
			inner.lower[d] += gw[d];
			inner.upper[d] -= gw[d];
		}
	{
		GlobalArray<ndim,double> ga(dom, gw, true);
		ga = src;
		ga.update();
		t.add_result(check(inner, ga == src));
		t.add_result(check(inner, shift(ga,gw) == shift(src,gw) && shift(ga,-gw) == shift(src,-gw)));
	}
	t.end_test();
}


int main(int argc, char* argv[]) {
	Init(&argc, &argv);
	SetupThreads();
	tester t(cerr);
	{
		coords<3> n = {{100,100,100}};
		Domain<3> dom(world(), n);
		auto f = sin(coord_val<0>(dom)) * cos(coord_val<1>(dom));
		make_suite1(dom, f).run(t);
	}
	Finalize();
}
