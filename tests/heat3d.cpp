/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <thread>
#include <cmath>
#include <iostream>
#include <chrono> 
#include <sstream>
#include <functional>
#include <getopt.h>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types3d;

coords split(const std::string &s) {
    std::stringstream ss(s);
    std::string item;
    coords ret;


    for(int i=0; i<3; ++i) {
        std::getline(ss, item, ',');
        ret[i] = atoi(item.c_str());
   }

   return ret;
}

double roundOff(double n) {
    double d = n * 100.0;
    int i = d + 0.5;
    d = (float)i / 100.0;
    return d;
}

std::pair<double, std::string> convertSize(size_t size) {               
    std::vector<std::string> SIZES = { "B", "KB", "MB", "GB" };
    unsigned div = 0;
    size_t rem = 0;

    while (size >= 1024 && div < SIZES.size()) {
        rem = (size % 1024);
        div++;
        size /= 1024;
    }

    double size_d = (double)size + (double)rem / 1024.0;
    return make_pair(size_d, SIZES[div]);
}

void init_pyramid(GlobalArrayD& ga)
{
	const coords_range outer = { ga.domain().total().lower - ga.ghost_width(), ga.domain().total().upper + ga.ghost_width() };
	const vecD one = {{1.0, 1.0, 1.0}};
	const vecD mid = {{0.5, 0.5, 0.5}};

	ga <<= 1.0 - max_element(abs(coord_vec(ga.domain(), outer, one) - mid)) / 0.5;
}

static int degree = 1;

template<typename S>
struct Heat_3d_7pt {
    const double nu = 0.25 * 0.25 * 0.25;
    inline double operator()(const typename S::accessor& u, coords ii) const {
#if 1
        const coords  left = {{ -1,  0,  0}};
        const coords right = {{  1,  0,  0}};
        const coords above = {{  0, -1,  0}};
        const coords below = {{  0,  1,  0}};
        const coords  back = {{  0,  0, -1}};
        const coords front = {{  0,  0,  1}};

        return u(ii) + nu * (u(ii + left) + u(ii + right) + u(ii + above) + u(ii + below) + u(ii + back) + u(ii + front) - 6.0 * u(ii));
#else 
        coords  left;  left[0] = ii[0]-1; left[1] = ii[1]-0; left[2] = ii[2]+0;
        coords right;  right[0] = ii[0]+1; right[1] = ii[1]+0; right[2] = ii[2]+0;
        coords above;  above[0] = ii[0]+0; above[1] = ii[1]-1; above[2] = ii[2]+0;
        coords below;  below[0] = ii[0]+0; below[1] = ii[1]+1; below[2] = ii[2]+0;
        coords  back;  back[0] = ii[0]+0; back[1] = ii[1]+0; back[2] = ii[2]-1;
        coords front;  front[0] = ii[0]+0; front[1] = ii[1]+0; front[2] = ii[2]+1;

        return u(ii) + nu * (u(left) + u(right) + u(above) + u(below) + u(back) + u(front) - 6.0 * u(ii));
#endif

    }
};

template<typename S>
struct Heat_3d_19pt {
    const double nu = 0.25 * 0.25 * 0.25;
    inline double operator()(const typename S::accessor& u, coords ii) const {
        const coords  left = {{ -1,  0,  0}};
        const coords right = {{  1,  0,  0}};
        const coords above = {{  0, -1,  0}};
        const coords below = {{  0,  1,  0}};
        const coords  back = {{  0,  0, -1}};
        const coords front = {{  0,  0,  1}};

        double ret = u(ii);
        
        for(int d=0; d<degree; d++) {
            ret += nu * (
                u(ii + left*1) + u(ii + right*1) + u(ii + above*1) + u(ii + below*1) + u(ii + back*1) + u(ii + front*1)
              + u(ii + left*2) + u(ii + right*2) + u(ii + above*2) + u(ii + below*2) + u(ii + back*2) + u(ii + front*2)
              + u(ii + left*3) + u(ii + right*3) + u(ii + above*3) + u(ii + below*3) + u(ii + back*3) + u(ii + front*3)
              -  18.0 * u(ii));
        }

        return ret;

    }
};

template<typename S>
void heat_overlap(GlobalArrayD& ga, GlobalArrayD& gb)
{
    auto f = ga.iupdate();
    {
        SHARK_COUNTER("inner");
        gb.region(gb.inner()) = unary(ga, S());
    }
    {
        SHARK_COUNTER("wait");
        f.wait();
    }
    {
        SHARK_COUNTER("outer");
        for (auto r : gb.outer()) gb.region(r, true) = unary(ga, S());
    }
}

template<typename S>
void heat(GlobalArrayD& ga, GlobalArrayD& gb)
{
    ga.update();
    {
    SHARK_COUNTER("total");
    gb <<= unary(ga, S());
    }
}

template<typename S>
void heat_only_update(GlobalArrayD& ga, GlobalArrayD&)
{
    ga.update();
}

template<typename S>
void heat_only_S(GlobalArrayD& ga, GlobalArrayD& gb)
{
    SHARK_COUNTER("total");
    gb <<= unary(ga, S());
}

template<typename S>
void heat_only_S_inner_outer(GlobalArrayD& ga, GlobalArrayD& gb)
{
    {
        SHARK_COUNTER("inner");
        gb.region(gb.inner()) = unary(ga, S());
    }

    {
        SHARK_COUNTER("outer");
        for (auto r : gb.outer()) gb.region(r, true) = unary(ga, S());
    }
 }

template<typename S>
void heat_nowait(GlobalArrayD& ga, GlobalArrayD& gb)
{
    ga.iupdate();
    {
    SHARK_COUNTER("total");
    gb <<= unary(ga, S());
    }
}

template<typename S>
void heat_only_ghost_nowait(GlobalArrayD& ga, GlobalArrayD&)
{
    ga.iupdate();
}

template<typename S>
void heat_ghost_sleep_coarse(GlobalArrayD& ga, GlobalArrayD&)
{
    auto f = ga.iupdate();
    std::this_thread::sleep_for(std::chrono::seconds(1));
    f.wait();
}

template<typename S>
void heat_ghost_sleep_fine(GlobalArrayD& ga, GlobalArrayD&)
{
    auto f = ga.iupdate();
    for(int i=0; i<10; ++i) {
	    f.test(); // progress...
	    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    f.wait();
}


template<typename S>
void heat_loop(GlobalArrayD& ga, GlobalArrayD& gb, int iter, unsigned mode) {
        SHARK_COUNTER("main");
	for(int k = 0; k < iter; k+=2)
	{
            switch (mode) {
                case 0 : heat<S>(ga, gb); heat<S>(gb, ga); break;
                case 1 : heat_overlap<S>(ga, gb); heat_overlap<S>(gb, ga); break;
                case 2 : heat_only_update<S>(ga, gb); heat_only_update<S>(gb, ga); break;
                case 3 : heat_only_S<S>(ga, gb); heat_only_S<S>(gb, ga); break;
                case 4 : heat_only_S_inner_outer<S>(ga, gb); heat_only_S_inner_outer<S>(gb, ga); break;
                case 5 : heat_nowait<S>(ga, gb); heat_nowait<S>(gb, ga); break;
                case 6 : heat_only_ghost_nowait<S>(ga, gb); heat_only_ghost_nowait<S>(gb, ga); break;
                case 7 : heat_ghost_sleep_coarse<S>(ga, gb); heat_ghost_sleep_coarse<S>(gb, ga); break;
                case 9 : heat_ghost_sleep_fine<S>(ga, gb); heat_ghost_sleep_fine<S>(gb, ga); break;
                default: Abort(1);
            }
	}
}

int main(int argc, char **argv)
{
        static const vector<const char *> mode_names = 
                                    { "iupdate+wait+total",
                                      "iupdate+inner+wait+outer",
                                      "iupdate+wait",
                                      "total",
                                      "inner+outer",
                                      "iupdate+total",
                                      "iupdate",
                                      "iupdate+sleep+wait",
                                      "iupdate+inner+(wait+outer)xn",
                                      "iupdate+sleep_progress+wait",
                                    };
	Init(&argc, &argv);

	//Default values that might be overridden by options
	int iter = -1; 
	int n = 511;	//incremented by one internally
        int npt = 19;
	bool block = false;
        unsigned mode = 0;
        coords bs = coords();

	//Process arguments
	int ch;

	while((ch = getopt(argc, argv, "n:i:t:bm:d:hp:c:")) != -1)
	{
                int v = -1;
                if (optarg) v = atoi(optarg);
		switch(ch)
		{
			case 'n': n = v;        break;
			case 'i': iter = v;     break;
			case 't': nthrds = v;   break;
			case 'p': npt = v;      break;
			case 'b': block = true; break;
			case 'm': mode = v;     break;
			case 'd': degree = v;     break;
			case 'c': bs = split(optarg); break;
			case 'h':
			case '?':
			default:
				cout << "Usage: " << argv[0] << " [-n <grid-cells>] [-i <iterations>] [-t <threads>]\n" 
                                     << "         [-m <mode>] [-b (blocked domain)] [-p <7|19 point stencil>]\n"
                                     << "         [-c <cache blocking size as v1,v2,v3>]\n"
                                     << "\n\tmodes:\n";
                                for(unsigned i = 0; i<mode_names.size(); ++i) {
                                     cout << "\t\t" << i << ": " << mode_names[i] << endl;
                                }
				Abort(1);
		}
	}

        assert((npt == 19) || (npt == 7));


	// if the number of iterations is not provided we calculate a
	// reasonable number for benchmarking at given size
	if (iter < 0)  iter = 4e6 * world().nprocs / n / n;

	if(world().procid == 0)
	{
		cout << "comm: " << comm << endl;
		cout << "sched: " << sched << endl;
                auto s = convertSize(n*n*n*sizeof(double));
		cout << "n: " << n  << " (" << s.first << s.second << ")" << endl;
		cout << "niter: " << iter << endl;
                cout << "degree: " << degree << endl;
		cout << "block: " << boolalpha << block << noboolalpha << endl;
		cout << "mode: " << mode_names[mode] << " (" << mode << ")" << endl;
		cout << "npt: " << npt << " point stencil" << endl;
		cout << "nprocs: " << world().nprocs << endl;
		cout << "nthrds: " << nthrds << endl;
                cout << "cache block sizes: " <<  bs << endl;
	}

	//Declare auxiliary variables and run	
	double nu = 0.25 * 0.25 * 0.25;
	SetupThreads();

	{
		const coords size  = {{n,n,n}};
		const coords ghost = {{3,3,3}};

		const array<int,3> pcoords = {{ 0, 0, block ? 0 : 1 }};
		Domain d(world(), size, bs, pcoords);

		if(world().procid == 0)
			d.outputDistribution(cout);

		typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};

		GlobalArrayD ga(d, ghost, false, bd);
		GlobalArrayD gb(d, ghost, false, bd);
		world().sync();

		init_pyramid(ga);

		double e0 = norm1(ga) / n / n / n;
		double f0 = norm2(ga) / n / n / n;

		if(world().procid == 0)
		{
			cout << "e0: " << e0 << endl;
			cout << "f0: " << f0 << endl;
		}

		double starttime = Wtime();

                if (npt == 19) {
                    heat_loop<Heat_3d_19pt<GlobalArrayD>>(ga, gb, iter, mode);
                } else if (npt == 7) {
                    heat_loop<Heat_3d_7pt<GlobalArrayD>>(ga, gb, iter, mode);
                } else abort();

		double endtime = Wtime();
		double e1 = norm1(ga) / n / n / n;
		double f1 = norm2(ga) / n / n / n;

		if(world().procid == 0)
		{
			cout << "e1: " << e1 << endl;
			cout << "f1: " << f1 << endl;
			cout << "t: " << (nu / n / n / n)   * iter << endl;
                        cout << "iter: " << iter << endl;
			cout << "runtime: " << endtime - starttime << endl;
                        cout << "MLUP/s: " << (double)iter*n*n*n/(endtime - starttime)/1e6 << endl;
                        cout << "ns/LUP: " << 1e9*(endtime - starttime)/((double)iter*n*n*n) << endl;
                        perf_data.print();
		}
	}

	Finalize();

	return 0;
}

