
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <unistd.h>
#include <sys/time.h>

#include <shark.hpp>
#include "../src/comm_impl.hpp"
#include <unsupported/Eigen/SparseExtra>


#ifdef _OPENMP
#include <omp.h>
#else
#include <tbb/tbb.h>
#endif

#include "bpmf.h"

using namespace std;
using namespace Eigen;

using namespace shark;
using namespace shark::types2d;

const int num_feat = 30;

const double alpha = 2;
int nsims = 20;
const int burnin = 5;

double mean_rating = .0;

int num_p = 0;
int num_m = 0;

typedef SparseMatrix<double> SparseMatrixD;
SparseMatrixD M, Mt, P, Pt;

long int total_exchanged_movies = 0;
long int total_exchanged_movies_elements = 0;
long int total_exchanged_users = 0;
long int total_exchanged_users_elements = 0;

VectorXd mu_u(num_feat);
VectorXd mu_m(num_feat);
MatrixXd Lambda_u(num_feat, num_feat);
MatrixXd Lambda_m(num_feat, num_feat);
MatrixXd sample_u;
MatrixXd sample_m;

// parameters of Inv-Whishart distribution (see paper for details)

const int b0_u = 2;
const int df_u = num_feat;
VectorXd mu0_u(num_feat);

MatrixXd WI_u(num_feat, num_feat);
MatrixXd WI_m(num_feat, num_feat);

const int b0_m = 2;
const int df_m = num_feat;
VectorXd mu0_m(num_feat);

void init() {
    mean_rating = M.sum() / M.nonZeros();
    Lambda_u.setIdentity();
    Lambda_m.setIdentity();

    sample_u = MatrixXd(num_feat,M.rows());
    sample_m = MatrixXd(num_feat,M.cols());
    sample_u.setZero();
    sample_m.setZero();

    // parameters of Inv-Whishart distribution (see paper for details)
    WI_u.setIdentity();
    mu0_u.setZero();

    WI_m.setIdentity();
    mu0_m.setZero();
}

inline double sqr(double x) { return x*x; }
std::pair<double,double> eval_probe_vec(int n, VectorXd & predictions, const MatrixXd &sample_m, const MatrixXd &sample_u, double mean_rating,
		GlobalArrayD& users)
{
    double se = 0.0, se_avg = 0.0;
    unsigned idx = 0;

    auto range = users.domain().local();
    for (int k=range.lower[1]; k<std::min(range.upper[1], (long)Pt.outerSize()); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(Pt,k); it; ++it)
        {
            const double pred = sample_m.col(it.row()).dot(sample_u.col(it.col() - users.domain().local().lower[1])) + mean_rating;
            se += sqr(it.value() - pred);

            const double pred_avg = (n == 0) ? pred : (predictions[idx] + (pred - predictions[idx]) / n);
            se_avg += sqr(it.value() - pred_avg);
            predictions[idx++] = pred_avg;
        }
    }


    const double rmse = sqrt( se / idx );
    const double rmse_avg = sqrt( se_avg / idx );
    return std::make_pair(rmse, rmse_avg);
}

void sample_movies(MatrixXd &s, int mm, const SparseMatrixD &mat, double mean_rating,
const MatrixXd &samples, int alpha, const VectorXd &mu_u, const MatrixXd &Lambda_u, GlobalArrayD& users,double* communication_time)
{
	int i = 0;
        MatrixXd MM(num_feat, num_feat); MM.setZero();
	VectorXd rr(num_feat); rr.setZero();

	for (SparseMatrixD::InnerIterator it(mat,mm); it; ++it, ++i)
	{
	       // cout << "M[" << it.row() << "," << it.col() << "] = " << it.value() << endl;
		if ( (users.domain().local().lower[1] > it.row()) || (it.row() >= users.domain().local().upper[1]))
		{
			vector<double> communication_buffer(num_feat);

			double s = tick();

			coords_range missing_user;

			missing_user.lower[1] = it.row();
			missing_user.upper[1] = it.row() + 1;

		    missing_user.lower[0] = 0;
		    missing_user.upper[0] = num_feat;

		    //#pragma omp critical
			   	users.get(missing_user,communication_buffer.data());

			total_exchanged_users_elements ++;


			*communication_time += tick() - s;

			auto col = Map<VectorXd> (communication_buffer.data(), num_feat);

			MM.noalias() += col * col.transpose();
			rr.noalias() += col * ((it.value() - mean_rating) * alpha);
		 }
		 else
		 {
			 auto col = samples.col(it.row() - users.domain().local().lower[1]);
			 MM.noalias() += col * col.transpose();
			 rr.noalias() += col * ((it.value() - mean_rating) * alpha);
		 }


	}

	Eigen::LLT<MatrixXd> chol = (Lambda_u + alpha * MM).llt();

	VectorXd tmp = rr + Lambda_u * mu_u;
	chol.matrixL().solveInPlace(tmp);
	tmp += nrandn(num_feat);
	chol.matrixU().solveInPlace(tmp);
	s.col(mm) = tmp;
}

void sample_users(MatrixXd &s, int mm, const SparseMatrixD &mat, double mean_rating,
    const MatrixXd &samples, int alpha, const VectorXd &mu_u, const MatrixXd &Lambda_u, int idx)
{
    int i = 0;
    MatrixXd MM(num_feat, num_feat); MM.setZero();
    VectorXd rr(num_feat); rr.setZero();
    for (SparseMatrixD::InnerIterator it(mat,mm); it; ++it, ++i) {
        // cout << "M[" << it.row() << "," << it.col() << "] = " << it.value() << endl;
        auto col = samples.col(it.row());
        MM.noalias() += col * col.transpose();
        rr.noalias() += col * ((it.value() - mean_rating) * alpha);
    }

    Eigen::LLT<MatrixXd> chol = (Lambda_u + alpha * MM).llt();
    if(chol.info() != Eigen::Success) {
      throw std::runtime_error("Cholesky Decomposition failed!");
    }

    VectorXd tmp = rr + Lambda_u * mu_u;
    chol.matrixL().solveInPlace(tmp);
    tmp += nrandn(num_feat);
    chol.matrixU().solveInPlace(tmp);
    s.col(idx) = tmp;
}

int main(int argc, char *argv[])
{
	Init(&argc, &argv);
        {

	bool random_assignment = false;
	bool block = false;
    int ch;
    string fname;
    string probename;

    while((ch = getopt(argc, argv, "n:t:b:r:p:i:")) != -1)
    {
        switch(ch)
        {
            case 'i': nsims = atoi(optarg); break;
            case 't': nthrds = atoi(optarg); break;
            case 'n': fname = optarg; break;
            case 'p': probename = optarg; break;
            case 'r': random_assignment = true; break;
            case 'b': block = true; break;
            case '?':
            default:
                cout << "Usage: " << argv[0] << " [-t <threads>]  [-r (random_assignement)] [ -b (block assignment) ]" 
                     << "[ -i <niters> ] -n <samples.mtx> -p <probe.mtx>"
                     << endl;
                Abort(1);
        }
    }

    if (fname.empty() || probename.empty()) { 
        cout << "Usage: " << argv[0] << " [-t <threads>] [-p <problem>] [-r <random_assignement>]" 
             << " -n <samples.mtx> -p <probe.mtx>"
             << endl;
        Abort(1);
    }

    assert(fname.c_str() && "filename missing");

    loadMarket(M, fname);

    Mt = M.transpose();

    loadMarket(P, probename);

    Pt = P.transpose();

    assert(M.nonZeros() > P.nonZeros());

    num_m = M.cols();
    num_p = M.rows();

	init();

	SetupThreads();

	typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};

	double sampling_users_time = 0, sampling_movies_time = 0, wishart_users_time = 0,  wishart_movies_time = 0;
    double get_user_time = 0, receiving_time = 0, sending_time = 0;

	coords gw = { { 1,1 } };
	const array<int,2> pcoords = {{ block ? 0 : 1, 0}};
	coords size = { {num_feat + 1, num_p + 1} };
	Domain d(world(), size, pcoords);

	GlobalArrayD users(d,gw, false, bd);

	VectorXd predictions(VectorXd::Zero( P.nonZeros() ));

	long double  average_sampling_sec =0;
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

	if (world().nprocs != 0) {
		if (random_assignment == 0) //assign a movie to the process that owns the biggest number of users that rated that movie
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
        }

		auto start = tick();

    	std::vector<int> movies_assigned_to_me(world().nprocs,0);

    	for (int u = 0; u < num_m; u++)
    		for(int k = 0; k < world().nprocs; k++)
    			if (assigned_movies[u] == k) movies_assigned_to_me[k]++;

    	cout << "I am process " << world().procid << " I will sample " << movies_assigned_to_me[world().procid] << " movies " <<endl;

        SparseMatrixD Mt = M.transpose();

        sample_u = Map<MatrixXd>(&users.da(users.domain().local().lower), num_feat, (my_users_max - my_users_min));

        for(int i=0; i<nsims; ++i)
        {
          double s = tick();

    	  // Sample from user hyperparams
          tie(mu_u, Lambda_u) = CondNormalWishart(sample_u, mu0_u, b0_u, WI_u, df_u);


          wishart_users_time  += tick() - s;

          s = tick();

          // Sample from movie hyperparams
          tie(mu_m, Lambda_m) = CondNormalWishart(sample_m, mu0_m, b0_m, WI_m, df_m);

          wishart_movies_time += tick() - s;


          s = tick();

		  #pragma omp parallel for
          	  for(int uu = my_users_min; uu < my_users_max - 1; uu++)
          	  {
          		  sample_users(sample_u, uu, Mt, mean_rating, sample_m, alpha, mu_u, Lambda_u, uu - my_users_min);
          	  }


          sampling_users_time += tick() - s;

          MPI_Status status;

          vector<int> flagsm(num_m, 0); //receive info about a movie only once, this is used in the case this process have more than one user that have rated that movie

          for(int mm = 0; mm < num_m; mm++)
          {
       	   if (assigned_movies[mm] == world().procid) //this movie is assigned to me so I sample it and then I send the new movie data to all the users that are concerned
       	   {
       		  double s = tick();

       		  sample_movies(sample_m, mm, M, mean_rating, sample_u, alpha, mu_m, Lambda_m, users, &get_user_time);


        	   sampling_movies_time += tick() - s;

       		   s = tick();

       		   vector<int> flagsp(world().nprocs, 0);

    		   for (SparseMatrixD::InnerIterator it(M,mm); it; ++it)
        	   {
    			   for(int k = 0; k < world().nprocs; k++)
    				   if ( ( users.domain().local(k).lower[1] <= it.row() ) && ( it.row() < users.domain().local(k).upper[1] ) && flagsp[k] == 0 && k != world().procid )
    			  	   {
    		               flagsp[k] = 1;

    		               MPI_Send(sample_m.col(mm).data(), num_feat, MPI_DOUBLE, k , 0, MPI_COMM_WORLD);

    		               total_exchanged_movies += num_feat;

    		               total_exchanged_movies_elements ++;
    			  	   }
        	   }


    		   sending_time += tick() - s;
       	   }
           else //the movie is not assigned to me but I may have users that rated it
           {
        	   double s = tick();

        	   for (SparseMatrixD::InnerIterator it(M,mm); it; ++it)
        	   {
         		   if ( (users.domain().local().lower[1] <= it.row()) && (it.row() < users.domain().local().upper[1]) && flagsm[mm] == 0 )
           		   {
           				flagsm[mm] = 1;

           			    MPI_Recv(recv.data(), num_feat, MPI_DOUBLE, assigned_movies[mm] , 0, MPI_COMM_WORLD, &status);

           				sample_m.col(mm) = Map<VectorXd> (recv.data(), num_feat);
           		   }
        	   }

    		   receiving_time += tick() - s;
           }
        }

          auto eval = eval_probe_vec( (i < burnin) ? 0 : (i - burnin), predictions, sample_m, sample_u, mean_rating, users);
           double norm_u = sample_u.norm();
           double norm_m = sample_m.norm();
           auto end = tick();
           auto elapsed = end - start;

           double samples_per_sec = (i + 1) * (M.rows() + M.cols()) / elapsed;

           if (world().procid == 0)  printf("Iteration %d:\t RMSE: %3.2f\tavg RMSE: %3.2f\tFU(%6.2f)\tFM(%6.2f)\tSamples/sec: %6.2f\n",
                   i, eval.first, eval.second, norm_u, norm_m, samples_per_sec);

           average_sampling_sec += samples_per_sec;
      }

      auto end = tick();
      auto elapsed = end - start;

      if (world().procid == 0)
      {
    	  cout << "Total time: " << elapsed <<endl <<flush;
		  cout << "Average Samples/sec: " << average_sampling_sec / nsims << endl <<flush;
		  cout << "Sampling users time: " << sampling_users_time << endl <<flush;
		  cout << "Sampling movies time: " << sampling_movies_time << endl <<flush;
		  cout << "Wishart users times " << wishart_users_time << endl <<flush;
		  cout << "Wishart movies times " << wishart_movies_time << endl <<flush;
		  cout << "Receiving time: " << receiving_time << endl <<flush;
		  cout << "Sending time: " << sending_time << endl <<flush;
		  cout << "Get user time: " << get_user_time << endl <<flush;
		  cout << "Sampling users time: " << sampling_users_time << endl <<flush;
	  }

	  cout << "Total exchanged movies elements: " << total_exchanged_movies_elements / nsims << endl;
	  cout << "Total exchanged users elements: " << total_exchanged_users_elements / nsims << endl;

        }
      Finalize();

      return 0;
}
