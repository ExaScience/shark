
#include <stdlib.h>     /* srand, rand */

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <shark.hpp>
#include <mpi.h>
#include <omp.h>

#include "bpmf.h"

using namespace std;
using namespace Eigen;
using namespace shark;
using namespace shark::types2d;

typedef SparseMatrix<double> SparseMatrixD;

const int num_feat = 100;
unsigned num_p = 0;
unsigned num_m = 0;

const int alpha = 2;
int nsims = 20;
const int burnin = 5;

double mean_rating = .0;

SparseMatrixD M;
typedef Eigen::Triplet<double> T;
vector<T> probe_vec;

VectorXd mu_u(num_feat);
VectorXd mu_m(num_feat);

MatrixXd Lambda_u(num_feat, num_feat);
MatrixXd Lambda_m(num_feat, num_feat);

MatrixXd sample_u;
MatrixXd sample_m;

// parameters of Inv-Whishart distribution (see paper for details)
MatrixXd WI_u(num_feat, num_feat);
const int b0_u = 2;
const int df_u = num_feat;
VectorXd mu0_u(num_feat);

MatrixXd WI_m(num_feat, num_feat);
const int b0_m = 2;
const int df_m = num_feat;
VectorXd mu0_m(num_feat);

void loadNetflix_trainMatrix(const char* fname)
{
    std::vector<T> lst;
    lst.reserve(100000000);

    FILE *f = fopen(fname, "r");
    assert(f && "Could not open file");

    // skip header
    char buf[2048];
    fscanf(f, "%s\n", buf);

    // data
    unsigned i, j;
    double v;
    double tmp;
    while (!feof(f))
    {
        if (!fscanf(f, "%d,%d,%lg,%d-%d-%d\n", &i, &j, &v, &tmp, &tmp, &tmp)) continue;

        i--;
        j--;

        num_p = std::max(num_p, j);
        num_m = std::max(num_m, i);
        mean_rating += v;
        lst.push_back(T(j,i,v));
    }

    num_p++;
    num_m++;
    mean_rating /= lst.size();
    fclose(f);

    M = SparseMatrix<double>(num_p, num_m);

    M.setFromTriplets(lst.begin(), lst.end());
}

void loadNetflix_probeVector(const char* fname)
{
    FILE *f = fopen(fname, "r");
    assert(f && "Could not open file");

    // skip header
    char buf[2048];
    fscanf(f, "%s\n", buf);

    // data
    unsigned i, j;
    double v;
    double tmp;
    while (!feof(f))
    {
        if (!fscanf(f, "%d,%d,%lg,%d-%d-%d\n", &i, &j, &v, &tmp, &tmp, &tmp)) continue;
        i--;
        j--;
        probe_vec.push_back(T(j,i,v));
    }

    fclose(f);
}

void loadChemo(const char* fname)
{
    std::vector<T> lst;
    lst.reserve(100000);
    
    FILE *f = fopen(fname, "r");
    assert(f && "Could not open file");

    // skip header
    char buf[2048];
    fscanf(f, "%s\n", buf);

    // data
    unsigned i, j;
    double v;

    while (!feof(f)) {
        if (!fscanf(f, "%d,%d,%lg\n", &i, &j, &v)) continue;
        i--;
        j--;

        if ((rand() % 5) == 0) {
            probe_vec.push_back(T(i,j,log10(v)));
        } 
#ifndef TEST_SAMPLE
        else // if not in test case -> remove probe_vec from lst
#endif
        {
            num_p = std::max(num_p, i);
            num_m = std::max(num_m, j);
            mean_rating += v;
            lst.push_back(T(i,j,log10(v)));
        }
    }
    num_p++;
    num_m++;
    mean_rating /= lst.size();
    fclose(f);

    M = SparseMatrix<double>(num_p, num_m);

    M.setFromTriplets(lst.begin(), lst.end());
}

void init()
{
    mean_rating = M.sum() / M.nonZeros();

    Lambda_u.setIdentity();
    Lambda_m.setIdentity();

    sample_m = MatrixXd(num_feat, num_m);
    sample_m.setZero();

    // parameters of Inv-Whishart distribution (see paper for details)
    WI_u.setIdentity();
    mu0_u.setZero();

    WI_m.setIdentity();
    mu0_m.setZero();
}

pair<double,double> eval_probe_vec_netflix(const vector<T> &probe_vec, const MatrixXd &sample_m, const MatrixXd &sample_u, double mean_rating , GlobalArrayD& users)
{
    unsigned n = probe_vec.size();
    unsigned correct = 0;
    unsigned done = 0;
    double diff = .0;
    for(auto t : probe_vec) {
    	if ( (users.domain().local().lower[1] <= t.row()) && (t.row() < users.domain().local().upper[1]))
    	{
    		double prediction = sample_m.col(t.col()).dot(sample_u.col(t.row() - users.domain().local().lower[1])) + mean_rating;
    		correct += (t.value() - prediction) * (t.value() - prediction) ;
    		diff += abs(t.value() - prediction);
    		done++;
    	}
    }
   
    return std::make_pair((double)correct / done, diff / done);
}

pair<double,double> eval_probe_vec(const vector<T> &probe_vec, const MatrixXd &sample_m, const MatrixXd &sample_u, double mean_rating , GlobalArrayD& users)
{
    unsigned n = probe_vec.size();
    unsigned correct = 0;
    unsigned done = 0;
    double diff = .0;
    for(auto t : probe_vec) {
    	if ( (users.domain().local().lower[1] <= t.row()) && (t.row() < users.domain().local().upper[1]))
    	{
    		double prediction = sample_m.col(t.col()).dot(sample_u.col(t.row() - users.domain().local().lower[1])) + mean_rating;
    		correct += (t.value() < log10(200)) == (prediction < log10(200));
    		diff += abs(t.value() - prediction);
    		done++;
    	}
    }

    return std::make_pair((double)correct / done, diff / done);
}

void sample_movies(MatrixXd &s, int mm, const SparseMatrixD &mat, double mean_rating,
    const MatrixXd &samples, int alpha, const MatrixXd &mu_u, const MatrixXd &Lambda_u, GlobalArrayD& users)
{
    int i = 0;
    MatrixXd E(num_feat,mat.col(mm).nonZeros());
    VectorXd rr(mat.col(mm).nonZeros());
    //cout << "movie " << endl;
    for (SparseMatrixD::InnerIterator it(mat,mm); it; ++it, ++i) {

      //   cout <<world().procid  <<  "M[" << it.row() << "," << it.col() << "] = " << it.value() << endl;

    	if ( (users.domain().local().lower[1] > it.row()) || (it.row() > users.domain().local().upper[1]))
    	{
    		vector<double> communication_buffer(num_feat);

    		coords_range missing_user;

    		missing_user.lower[1] = it.row();
    		missing_user.upper[1] = it.row() + 1;

    	    missing_user.lower[0] = 0;
    	    missing_user.upper[0] = num_feat;

			//#pragma omp critical
    	    	users.get(missing_user,communication_buffer.data());

    	    E.col(i) = Map<VectorXd> (communication_buffer.data(), num_feat);
    	}
    	else
    	{
       		E.col(i) = samples.col(it.row() - users.domain().local().lower[1]);
    	}

        rr(i) = it.value() - mean_rating;
    }

    auto MM = E * E.transpose();
    MatrixXd MMs = alpha * MM.array();
    assert(MMs.cols() == num_feat && MMs.rows() == num_feat);
    MatrixXd covar = (Lambda_u + MMs).inverse();
    MatrixXd MMrr = (E * rr) * alpha;  
    auto U = Lambda_u * mu_u;
    auto mu = covar * (MMrr + U);

    MatrixXd chol = covar.llt().matrixL();
#ifdef TEST_SAMPLE
    auto r(num_feat); r.setConstant(0.25);
#else
    auto r = nrandn(num_feat);
#endif
    s.col(mm) = chol * r + mu;

#ifdef TEST_SAMPLE
      cout << "movie " << mm << ":" << result.cols() << " x" << result.rows() << endl;
      cout << "mean rating " << mean_rating << endl;
      cout << "E = [" << E << "]" << endl;
      cout << "rr = [" << rr << "]" << endl;
      cout << "MM = [" << MM << "]" << endl;
      cout << "Lambda_u = [" << Lambda_u << "]" << endl;
      cout << "covar = [" << covar << "]" << endl;
      cout << "mu = [" << mu << "]" << endl;
      cout << "chol = [" << chol << "]" << endl;
      cout << "rand = [" << r <<"]" <<  endl;
      cout << "result = [" << result << "]" << endl;
#endif

}

void sample_users(MatrixXd &s, int mm, const SparseMatrixD &mat, double mean_rating,
    const MatrixXd &samples, int alpha, const MatrixXd &mu_u, const MatrixXd &Lambda_, int idx)
{
    int i = 0;
    MatrixXd E(num_feat,mat.col(mm).nonZeros());
    VectorXd rr(mat.col(mm).nonZeros());
    //cout << "movie " << endl;

    for (SparseMatrixD::InnerIterator it(mat,mm); it; ++it, ++i)
    {
        // cout << "M[" << it.row() << "," << it.col() << "] = " << it.value() << endl;

    	E.col(i) = samples.col(it.row());
        rr(i) = it.value() - mean_rating;
    }

    auto MM = E * E.transpose();
    MatrixXd MMs = alpha * MM.array();
    assert(MMs.cols() == num_feat && MMs.rows() == num_feat);
    MatrixXd covar = (Lambda_u + MMs).inverse();
    MatrixXd MMrr = (E * rr) * alpha;
    auto U = Lambda_u * mu_u;
    auto mu = covar * (MMrr + U);

    MatrixXd chol = covar.llt().matrixL();
#ifdef TEST_SAMPLE
    auto r(num_feat); r.setConstant(0.25);
#else
    auto r = nrandn(num_feat);
#endif
    s.col(idx) = chol * r + mu;

#ifdef TEST_SAMPLE
      cout << "movie " << mm << ":" << result.cols() << " x" << result.rows() << endl;
      cout << "mean rating " << mean_rating << endl;
      cout << "E = [" << E << "]" << endl;
      cout << "rr = [" << rr << "]" << endl;
      cout << "MM = [" << MM << "]" << endl;
      cout << "Lambda_u = [" << Lambda_u << "]" << endl;
      cout << "covar = [" << covar << "]" << endl;
      cout << "mu = [" << mu << "]" << endl;
      cout << "chol = [" << chol << "]" << endl;
      cout << "rand = [" << r <<"]" <<  endl;
      cout << "result = [" << result << "]" << endl;
#endif

}

int main(int argc, char *argv[])
{
	Init(&argc, &argv);

	bool random_assignment = false;
	bool block = false;
    int ch;

    string fname;
    string probevectorname;
    string problem;

    while((ch = getopt(argc, argv, "n:t:b:p:r:v:i:")) != -1)
    {
    	switch(ch)
    	{
    		case 'i':
				{
					istringstream iss(optarg);
					iss >> nsims;
				}
				break;
    	    case 't':
    			{
    				istringstream iss(optarg);
    				iss >> nthrds;
    			}
    			break;
    		case 'n':
    			{
    				istringstream iss(optarg);
    				iss >> fname;
    			}
    			break;
    		case 'v':
    		    {
    		    	istringstream iss(optarg);
    		    	iss >> probevectorname;
    		    	}
    		    break;
    		case 'p':
    			{
    				istringstream iss(optarg);
    				iss >> problem;
    			}
    			break;
    		case 'r':
    			{
    				istringstream iss(optarg);
    				iss >> random_assignment;
    			}
    			break;
    		case 'b':
    				block = true;
    				break;
    		case '?':
    				cout << "Usage: " << argv[0] << " [-n <path_to_the_data>] [-t <threads>] [-p <problem>] [-r <random_assignement>] \n" << endl;
    				Abort(1);
    		}
    	}

    assert(fname.c_str() && "filename missing");

    Eigen::initParallel();
    Eigen::setNbThreads(16);

	if ( problem.compare("chemo") == 0 ) loadChemo(fname.c_str());
	if ( problem.compare("Netflix") == 0 )
	{
		loadNetflix_trainMatrix(fname.c_str());
		loadNetflix_probeVector(probevectorname.c_str());
	}

	init();

	SetupThreads();

	typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};

	coords gw = { { 1,1 } };
	const array<int,2> pcoords = {{ block ? 0 : 1, 0}};
	coords size = { {num_feat + 1, num_p + 1} };
	Domain d(world(), size, pcoords);

	GlobalArrayD users(d,gw, false, bd);

	if(world().procid == 0)
	{
		cout << "num_feat: " << num_feat<<endl;
		cout << "sched: " << sched << endl;
		cout << "block: " << boolalpha << block << noboolalpha << endl;
		cout << "nprocs: " << world().nprocs << endl;
		cout << "nthrds: " << nthrds << endl;

		d.outputDistribution(cout);
	}

	int my_users_min = users.domain().local().lower[1];
	int my_users_max = users.domain().local().upper[1];

	std::vector<int> movie_bids(world().nprocs,0); //contains the indices of the movies that this process will sample
	std::vector<int> assigned_movies(num_m, 0); //contains the indice of the process that will sample the movie assigned_movies[0] = 1 the process number 1 will sample the movies number 0 ...

	vector<double> recv(num_feat); //first element is the movie index to update

	if (world().nprocs != 0)
		if (random_assignment == 0 ) //assign a movie to the process that owns the biggest number of users that rated that movie
		{
			for (int mm = 0; mm < num_m; mm++)
			{
				int max = 0;

				for(int k = 0; k < world().nprocs; k++) movie_bids[k] = 0;

				for (SparseMatrixD::InnerIterator it(M,mm); it; ++it)
				{
					for(int k = 0; k < world().nprocs; k++)
					{
						if ( ( d.local(k).lower[1] <= it.row() ) && ( it.row() <= d.local(k).upper[1] ) ) movie_bids[k]++;
						if ( movie_bids[k] > movie_bids[max] ) max = k;
					}

					assigned_movies[mm] = max;
				}
			}
		}
		else //random assignment of the movies to the processes
		{
			int ratio = num_m / world().nprocs;

			for (int u = 0; u < num_m; u++)
			{
				for(int k = 0; k < world().nprocs; k++)
					if ( (int) (u/ratio) == k )
						assigned_movies[u] = k;

				if ( (int) (u/ratio) == world().nprocs ) assigned_movies[u] = 0;
			}
		}

        auto start = chrono::duration_cast<chrono::duration<double>>(chrono::high_resolution_clock::now().time_since_epoch()).count();

        SparseMatrixD Mt = M.transpose();

        sample_u = MatrixXd(num_feat, (my_users_max - my_users_min));
        sample_u = Map<MatrixXd>(users.ptr, num_feat, (my_users_max - my_users_min));

        for(int i=0; i<nsims; ++i)
        {
    	  // Sample from user hyperparams
          tie(mu_u, Lambda_u) = CondNormalWishart(sample_u, mu0_u, b0_u, WI_u, df_u);

          // Sample from movie hyperparams
          tie(mu_m, Lambda_m) = CondNormalWishart(sample_m, mu0_m, b0_m, WI_m, df_m);

		  #pragma omp parallel for
          	  for(int uu = my_users_min; uu < my_users_max - 1; uu++)
          	  {
          		  sample_users(sample_u, uu, Mt, mean_rating, sample_m, alpha, mu_u, Lambda_u, uu - my_users_min);
          	  }

          MPI_Status status;

          //receive info about a movie only once, this is used in the case this process have more than one user that have rated that movie
          vector<int> flagsm(num_m, 0);

          for(int mm = 0; mm < num_m; mm++)
          {
           //this movie is assigned to me so I sample it and then I send the new movie data to all the users that are concerned
       	   if (assigned_movies[mm] == world().procid)
       	   {
       		   sample_movies(sample_m, mm, M, mean_rating, sample_u, alpha, mu_m, Lambda_m, users);

       		   vector<int> flagsp(world().nprocs, 0);

    		   for (SparseMatrixD::InnerIterator it(M,mm); it; ++it)
        	   {
    			   for(int k = 0; k < world().nprocs; k++)
    				   if ( ( users.domain().local(k).lower[1] <= it.row() ) && ( it.row() < users.domain().local(k).upper[1] ) && flagsp[k] == 0 && k != world().procid )
    			  	   {
    		               flagsp[k] = 1;

    		               MPI_Send(sample_m.col(mm).data(), num_feat, MPI_DOUBLE, k , 0, MPI_COMM_WORLD);
    			       }
        	   }
       	   }
           else //the movie is not assigned to me but I may have users that rated it
           {
        	   for (SparseMatrixD::InnerIterator it(M,mm); it; ++it)
        	   {
         		   if ( (users.domain().local().lower[1] <= it.row()) && (it.row() < users.domain().local().upper[1]) && flagsm[mm] == 0 )
           		   {
           				flagsm[mm] = 1;

           			    MPI_Recv(recv.data(), num_feat, MPI_DOUBLE, assigned_movies[mm] , 0, MPI_COMM_WORLD, &status);

           				sample_m.col(mm) = Map<VectorXd> (recv.data(), num_feat);
           		   }
        	   }
           }
        }

       if ( problem.compare("chemo") == 0 )
       {
    	   auto eval = eval_probe_vec(probe_vec, sample_m, sample_u, mean_rating, users);
    	   double norm_u = sample_u.norm();
    	   double norm_m = sample_m.norm();
    	   auto end = chrono::duration_cast<chrono::duration<double>>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    	   auto elapsed = end - start;
    	   double samples_per_sec = (i + 1) * (num_p + num_m) / elapsed;

    	   if (world().procid == 0)
    	      printf("Iteration %d:\t num_correct: %3.2f%%\tavg_diff: %3.2f\tFU(%6.2f)\tFM(%6.2f)\tSamples/sec: %6.2f\n",
    	             i, 100*eval.first, eval.second, norm_u, norm_m, samples_per_sec);
       }
       else if ( problem.compare("Netflix") == 0 )
       {
    	   auto eval = eval_probe_vec_netflix(probe_vec, sample_m, sample_u, mean_rating, users);
    	   double norm_u = sample_u.norm();
       	   double norm_m = sample_m.norm();
       	   auto end = chrono::duration_cast<chrono::duration<double>>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    	   auto elapsed = end - start;
    	   double samples_per_sec = (i + 1) * (num_p + num_m) / elapsed;

    	   if (world().procid == 0)
    	      printf("Iteration %d:\t num_correct: %3.2f%%\t mrse: %3.2f%%\t avg_diff: %3.2f\tFU(%6.2f)\tFM(%6.2f)\tSamples/sec: %6.2f\n",
    	        i, 100*eval.first, sqrt(eval.first), eval.second, norm_u, norm_m, samples_per_sec);
       }

      }

      return 0;
}
