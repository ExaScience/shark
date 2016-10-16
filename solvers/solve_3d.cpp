#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <cassert>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types3d;

struct Solver;

std::ostream &operator<<(std::ostream &os, const Solver &opts);

struct Solver {
    //Default values that might be overridden by options
    int nx = 65;	//incremented by one internally
    int ny = 65;	//incremented by one internally
    int nz = 65;	//incremented by one internally
    double tol = 1e-06;
    int maxit = 1000;
    int it = 0;
    int reps = 1;
    int degree = 1;
    ostream* out = NULL;
    string solvername = "cg";

    // gmres
    int restart = 30;
    enum GS {CLASSICAL, MODIFIED} GStype;
    // pipe gmres
    enum BASIS { MONOMIAL, NEWTON, CHEBYSHEV } basis;
    int l = 2;
    double lmin=0.0, lmax=0.0;

    void onlydot();
    void cg();
    void pipe_cg();
    void gmres();
    void pipe_gmres();
    void solve();

    double starttime;
    double mintime = 1.e+10;
    double maxtime = 1.e-10;

#ifdef SHARK_PROFILING
    TotalsCounter min_perf_data;
    Counter *main_counter;
#endif
    int r = 0;

    void start();
    void stop(int k);
    void print(double errorNorm);
    void setn(int n) { nx = ny = nz = n; }
    void parse_args(int argc, char **argv);
    
    // the GAs
    GlobalArrayD x_exact;
    GlobalArrayD b;
    GlobalArrayD x;
    Domain *dom;

    void init();
    void finalize();
};

#include "laplaceOperator3d.hpp"

void Solver::solve() {
    if      (solvername == "dot") return onlydot();
    else if (solvername == "cg") return cg();
    else if (solvername == "pipecg") return pipe_cg();
    else if (solvername == "gmres") return gmres();
    else if (solvername == "pipegmres") return pipe_gmres();
    Abort(1);
}

void Solver::start() {
    applyOperator(b, x_exact, degree);
    x <<= constant(x.domain(), 0.0); 
    starttime = Wtime();
#ifdef SHARK_PROFILING
    main_counter = new Counter("main");
#endif
}

void Solver::stop(int k) {
    double totaltime = Wtime() - starttime;

#ifdef SHARK_PROFILING
    delete main_counter;
#endif

    if (k>it) it = k;
    if (totaltime > maxtime)	maxtime = totaltime;

    if (totaltime < mintime) {
        mintime = totaltime;
#ifdef SHARK_PROFILING
        min_perf_data = perf_data;
#endif
    }

#ifdef SHARK_PROFILING
    perf_data = TotalsCounter();
#endif

    if (world().procid == 0) cout << "rep " << r++ << "/" << reps << ": " << totaltime << endl;
}

void Solver::print(double errorNorm) {
    if (world().procid != 0) return;
    cout << "runtime: " << mintime << endl;
    cout << "maxtime: " << maxtime << "\t(" << (int)(100.*maxtime/mintime) << "%)" << endl;
    assert(maxtime / mintime < 2.);
    cout << "err: " << errorNorm << endl;
    cout << "it: " << it << endl;
#ifdef SHARK_PROFILING
    min_perf_data.print();
#endif
}

void Solver::parse_args(int argc, char **argv)
{

    //Process arguments
    int ch;
    while((ch = getopt(argc, argv, "x:y:z:n:e:m:t:r:d:l:c:vh?")) != -1) {
        switch(ch) {
            case 'x': { istringstream iss(optarg); iss >> nx; } break;
            case 'y': { istringstream iss(optarg); iss >> ny; } break;
            case 'z': { istringstream iss(optarg); iss >> nz; } break;
            case 'n': { setn(atoi(optarg)); } break;
            case 'e': { istringstream iss(optarg); iss >> tol; } break;
            case 'm': { istringstream iss(optarg); iss >> maxit; } break;
            case 't': { istringstream iss(optarg); iss >> shark::nthrds; } break;
            case 'r': { istringstream iss(optarg); iss >> reps; } break;
            case 'd': { istringstream iss(optarg); iss >> degree; } break;
            case 'l': { istringstream iss(optarg); iss >> l; } break;
            case 'c': { solvername = optarg; } break;
            case 'v': { out = world().procid == 0 ? &cout : NULL; } break; 
            case 'h':
            case '?':
            default:
                      cout << "Usage: " << argv[0] << " [-x] [-y] [-z] [-e <tolerance>] [-m <max-iterations>] [-t <threads>]  [-c cg|cg3|cgd|cg3d|p1cg|cgsgmres|mgsgmres|p1gmres] \n" << endl;
                      Abort(1);
        }
    }

    if(world().procid == 0) cout << *this;
}

void Solver::init() {
    coords size = {{nx,ny,nz}};
    coords gw = {{1,1,1}};
    typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};
    dom = new Domain(world(), size);
    x_exact = GlobalArrayD(*dom, gw);
    x_exact <<= constant(*dom, 1.0);
    b = GlobalArrayD(*dom, gw);
    x = GlobalArrayD(*dom, gw);
}

void Solver::finalize() {
    GlobalArrayD error(x, true);
    error = error - x_exact;
    double errorNorm = norm2(error);
    print(errorNorm);
    delete dom;
}

std::ostream &operator<<(std::ostream &os, const Solver &opts)
{
    os << "Opts:" << endl;
    os << "  nx: " << opts.nx << endl;
    os << "  ny: " << opts.ny << endl;
    os << "  nz: " << opts.nz << endl;
    os << "  solver: " << opts.solvername << endl;
    os << "  reps: " << opts.reps << endl;
    os << "  degree: " << opts.degree << endl;
    os << "  maxit: " << opts.maxit << endl;
    os << "  tol: " << opts.tol << endl;
    os << "SHARK:" << endl;
    os << "  sched: " << shark::sched << endl;
    os << "  nprocs: " << shark::world().nprocs << endl;
    os << "  nthrds: " << shark::nthrds << endl;
    return os;
}

#include "cg.hpp"
#include "pipe_cg.hpp"
#include "gmres.hpp"
#include "pipe_gmres.hpp"

int main(int argc, char **argv)
{
    Init(&argc, &argv);
    SetupThreads();

    {
        Solver s;
        s.parse_args(argc, argv);
        s.init();
        s.solve();
        s.finalize();
    }

    Finalize();

    return 0;
}
