/*
 * BPMF.cpp
 *
 *  Created on: 5 févr. 2015
 *  Author: chakro23
 *
 * Copyright (c) 2015, imec
 *  
 */

#include "BPMF.h"
#include "Data.h"

#include <cstdlib>
#include <iostream>
#include <random>
#include <iostream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <cstddef>
#include <mpi.h>

#include <Dense>
#include <Cholesky>

using namespace std;
using namespace shark;
using namespace Eigen;

const int num_features = 20;

typedef Matrix<double, num_features, num_features> MatrixFeat;

/* Compute mean of a column of the current matrix */

double columnsMean(MatrixXd& gi, int column)
{
	double sum = 0.0;

	for(int k = 0; k < gi.rows(); k++)
		sum+= gi(k,column);

	return (double) sum/gi.rows();
}

/* Compute mean of a column of the current matrix */

double columnsMean(GlobalArrayD& gi, int column)
{
	const Domain& dom = gi.domain();
	const coords size = dom.local().upper - dom.local().lower;

	double sum = 0.0;

	for(int k = 0; k < size[0]; k++)
		sum+= gi.ptr[k * size[1] + column];

	return (double) sum/size[0];
}

void covariance(MatrixXd& ga, MatrixFeat&out_matrix)
{
	double sum, mean;
	int nb_rows = ga.rows();
	int nb_columns = ga.cols();

	std::vector<double> xi(nb_rows,0), yi(nb_rows,0);

	for (int i = 0 ; i < nb_columns; i++)
	{
		sum = 0.0;

		for (int j = 0 ; j < nb_rows; j++)
		{
			sum += ga(j,i);
			mean = sum/nb_rows;
		}

		double inner_product = 0;

		for (int j = 0 ; j < nb_rows; j++)
		{
			xi[j] = ga(j,i) - mean; //vector = column(i), vector[i]
			inner_product += xi[j] * xi[j];
		}

		out_matrix(i, i) = inner_product / (nb_rows - 1);

		for (int j = i + 1; j < nb_columns; j++)
		{
			sum = 0.0;

			inner_product = 0;

			for (int k = 0 ; k < nb_rows; k++)
			{
				sum += ga(k,j);
				mean = sum/nb_rows;
			}

			for (int k = 0 ; k < nb_rows; k++)
			{
				yi[k] = ga(k,j) - mean;
				inner_product += xi[k] * yi[k];
			}

			double val =  inner_product / (nb_rows - 1);

			out_matrix(i,j) = val;
			out_matrix(j,i) = val;
		}
	}
}

void covariance(GlobalArrayD& ga, MatrixFeat& out_matrix)
{
	const Domain& dom = ga.domain();
	const coords size = dom.local().upper - dom.local().lower;

	double sum, mean;
	int nb_rows = size[0];
	int nb_columns = size[1];

	std::vector<double> xi(nb_rows,0), yi(nb_rows,0);

	for (int i = 0 ; i < nb_columns; i++)
	{
		sum = 0.0;

		for (int j = 0 ; j < nb_rows; j++)
		{
			sum += ga.ptr[i + nb_columns * j];
			mean = sum/nb_rows;
		}

		double inner_product = 0;

		for (int j = 0 ; j < nb_rows; j++)
		{
			xi[j] = ga.ptr[i + nb_columns * j] - mean; //vector = column(i), vector[i]
			inner_product += xi[j] * xi[j];
		}

		out_matrix(i, i) =  inner_product / (nb_rows - 1);

		for (int j = i + 1; j < nb_columns; j++)
		{
			sum = 0.0;

			inner_product = 0;

			for (int k = 0 ; k < nb_rows; k++)
			{
				sum += ga.ptr[j + nb_columns* k];
				mean = sum/nb_rows;
			}

			for (int k = 0 ; k < nb_rows; k++)
			{
				yi[k] = ga.ptr[j + nb_columns * k] - mean;
				inner_product += xi[k] * yi[k];
			}

			double val =  inner_product / (nb_rows - 1);

			out_matrix(i,j) = val;
			out_matrix(j,i) = val;
		}
	}
}

void mult(MatrixFeat& A ,std::vector<double>& B, std::vector<double>& C)
{
	if (C.size() != A.rows()) C.resize(A.rows());

	for (int i = 0; i < A.rows(); i++)
	{
		double product = 0;

		for (int j = 0; j < A.cols(); j++)
			product += A(i,j) * B[j];

		C[i] = product;
	}
}

void mult(MatrixFeat& A ,MatrixFeat& B, MatrixFeat& C)
{
	if (C.rows() != A.rows() &&  C.cols() != B.cols())  C.resize(A.rows(), B.cols());

	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < B.cols(); j++)
		{
			double product = 0;

			for (int k = 0; k < A.cols(); k++)
				product += A(i,k) * B(k,j);

			C(i,j) = product;
		}
	}
}

MatrixFeat  wishart(MatrixFeat& scale, double df)
{
	LLT<MatrixXd> lltOfscale(scale); // compute the Cholesky decomposition of A
	MatrixFeat A = lltOfscale.matrixL();

	int p = scale.rows();

	std::vector<double> z(p*p,0);

	std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.5,1);

	for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
			z[i * p + j] = distribution(generator);

	std::vector<double> y(p);

	for (int i = 0; i < p; i++)
	{
		std::default_random_engine generator2;
		std::gamma_distribution<double> distribution2((df - (i + 1)) / 2 , 2.0);
		y[i] = distribution2(generator2);
	}

	std::vector<double> B(p * p);

	B[0] = y[0];

	if (p > 1)
	{
		// rest of diagonal:
		for (int j = 1; j < p; j++)
		{
			double inner_product = 0;

			for (int k = 0; k < j; k++)
				inner_product += z[k* p + j] * z[k* p + j];

			B[j * p + j] = y[j] + inner_product;
		}

		// first row and column:
		for (int j = 1; j < p; j++)
		{
			B[j] = z[j] * sqrt(y[0]);
			B[j * p] = B[j]; // mirror
		}
	}

	if (p > 2)
	{
		for (int j = 2; j < p; j++)
		{
			for (int i = 1; i <= j - 1; i++)
			{
				double inner_product = 0;

				for (int k = 0; k < i ; k++)
					inner_product += z[k * p + i] * z[k * p + j];

				B[i * p + j] = z[i * p + j] * sqrt(y[i] + inner_product);
				B[j * p + i] = B[i * p + j]; // mirror
			}
		}
	}

	//optimized A^t*A*B

	A.transposeInPlace();

	std::vector<double> row(A.cols(),0);

	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < A.cols(); j++)
		{
			double product = 0;

			for (int k = 0; k < A.cols(); k++)
				product += A(i,k) * A(k,j);

			row[j] = product;
		}

		for (int j = 0; j < A.cols(); j++)
		{
			double product = 0;

			for (int k = 0; k < A.cols(); k++)
				product += row[k] * B[k * p + j];

			A(i,j) = product;
		}
	}

	return A;
}

void sample_user_hyper_parameters(GlobalArrayD& users, int num_users, MatrixFeat& alpha_users, vector<double>& mu_u)
{
	vector<double> x_bar(num_features,0), mu0_u_x_bar(num_features,0), mu0_u(num_features,0), mu_temp(num_features,0), normalRdn(num_features,0);
	MatrixFeat S_bar, e1e2, WI_post, WU_I, lam;
	int b0_u = 2;
	double df_upost;
	int df_u = num_features;

	WU_I.setIdentity();

	for (int f = 0; f < num_features; f++)
		x_bar[f] = columnsMean(users,f);

	covariance(users,S_bar);

	for (int i = 0; i < num_features; i++) mu0_u_x_bar[i] = mu0_u[i] - x_bar[i];

	for (int i = 0; i < num_features; i++)
		for (int j = 0; j < num_features; j++)
			e1e2(i, j) = mu0_u_x_bar[i] * mu0_u_x_bar[j] * (num_users * b0_u / (b0_u + num_users + 0.0));

	WI_post = WU_I.inverse() + (S_bar * num_users) + e1e2;
	WI_post = WI_post.inverse();
	WI_post = WI_post + (WI_post.transpose() * 0.5);

	df_upost = df_u + num_users;

	alpha_users = wishart(WI_post, df_upost);

	for (int i = 0; i < num_features; i++)
		mu_temp[i] = ( (mu0_u[i] * b0_u) + (x_bar[i] * num_users) ) * (1 / (b0_u + num_users + 0.0));

	alpha_users = alpha_users * (b0_u + num_users);

	LLT<MatrixXd> lltOfusers(alpha_users.inverse());
	lam = lltOfusers.matrixL();
	lam.transposeInPlace();

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.5,1);

	for (int f = 0; f < num_features; f++)
		normalRdn[f] = distribution(generator);

	for (int i = 0; i < num_features; i++)
	{
		double inner_product = 0;

		for (int j = 0; j < num_features; j++)
			inner_product += lam(i,j) * normalRdn[j];

		mu_u[i] = inner_product + mu_temp[i];
	}
}

void sample_movie_hyper_parameters(MatrixXd movies, int num_movies, MatrixFeat& alpha_movies, vector<double>& mu_m)
{
	vector<double> x_bar(num_features,0), mu0_m_x_bar(num_features,0), mu0_m(num_features,0), mu_temp(num_features,0), normalRdn(num_features,0);
	MatrixFeat S_bar, e1e2, WI_post, WM_I, lam;
	int b0_m = 2;
	int df_m = num_features;
	double df_mpost;

	WM_I.setIdentity();

	for (int f = 0; f < num_features; f++)
		x_bar[f] = columnsMean(movies,f);

	covariance(movies, S_bar);

	for (int i = 0; i < num_features; i++) mu0_m_x_bar[i] = mu0_m[i] - x_bar[i];

	for (int i = 0; i < num_features; i++)
		for (int j = 0; j < num_features; j++)
			e1e2(i, j) = mu0_m_x_bar[i] * mu0_m_x_bar[j] * (num_movies * b0_m / (b0_m + num_movies + 0.0));

	WI_post = WM_I.inverse() + (S_bar * num_movies) + e1e2;
	WI_post = WI_post.inverse();
	WI_post = WI_post + (WI_post.transpose() * 0.5);

	df_mpost = df_m + num_movies;

	alpha_movies = wishart(WI_post, df_mpost);

	for (int i = 0; i < num_features; i++)
		mu_temp[i] = ( (mu0_m[i] * b0_m) + (x_bar[i] * num_movies) ) * (1 / (b0_m + num_movies + 0.0));

	alpha_movies = alpha_movies * (b0_m + num_movies);

	LLT<MatrixXd> lltOfmovies(alpha_movies.inverse());
	lam = lltOfmovies.matrixL();
	lam.transposeInPlace();

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.5,1);

	for (int f = 0; f < num_features; f++)
		normalRdn[f] = distribution(generator);

	for (int i = 0; i < num_features; i++)
	{
		double inner_product = 0;

		for (int j = 0; j < num_features; j++)
			inner_product += lam(i,j) * normalRdn[j];

		mu_m[i] = inner_product + mu_temp[i];
	}
}

int main(int argc, char** argv)
{
	Init(&argc, &argv);

	int iterations = 20;
	bool block = false;

	int ch;
	string problem;

	while((ch = getopt(argc, argv, "n:i:t:b:p:")) != -1)
	{
		switch(ch)
		{
			case 'i':
				{
					istringstream iss(optarg);
					iss >> iterations;
				}
				break;
			case 't':
				{
					istringstream iss(optarg);
					iss >> nthrds;
				}
				break;
			case 'p':
				{
					istringstream iss(optarg);
					iss >> problem;
				}
				break;
			case 'b':
				block = true;
				break;
			case '?':
			default:
				cerr << "Usage: " << argv[0] << " [-i <iterations>] [-t <threads>] [-p <problem>] \n" << endl;
				Abort(1);
		}
	}

	int num_users, num_movies, num_features = 20, numRates = 0;
	double globalMean;

	std::vector<double> rate_values;
	std::vector<int> user_indices;
	std::vector<int> movie_indices;

	Data data;

	printf("Start loading Data...\n");

	if (problem.compare("chemo") == 0 ) data.readData("chemo", rate_values, user_indices, movie_indices, num_users, num_movies, numRates, globalMean);
	else if (problem.compare("Netflix") == 0 ) data.readData("Netflix", rate_values, user_indices, movie_indices, num_users, num_movies, numRates, globalMean);
	else if (problem.compare("MovieLens") == 0 ) data.readData("MovieLens", rate_values, user_indices, movie_indices, num_users, num_movies, numRates, globalMean);

	printf("Data are: num_users %d num_movies %d numRates %d globalMean %f \n", num_users, num_movies, numRates, globalMean);

	printf("Finish loading Data...\n");

	coords gw = { { 1,1 } };
	const array<int,2> pcoords = {{ 0, block ? 0 : 1 }};
	coords size = { {num_users , num_features} }; //the output matrix would be a global array that is distributed horizontally
	Domain d(world(), size);

	if (world().procid == 0)
	{
		cerr << "sched: " << sched << endl;
		cerr << "iteration: " << iterations << endl;
		cerr << "block: " << boolalpha << block << noboolalpha << endl;
		cerr << "nprocs: " << world().nprocs << endl;
		cerr << "nthrds: " << nthrds << endl;

		d.outputDistribution(cerr);
	}

	GlobalArrayD users(d,gw); //numUsers*num_features
	MatrixXd movies(num_movies,num_features);

	std::vector<double> mu_u(num_features,0), mu_m(num_features,0), mean(num_features,0), tmp3(num_features,0);
	std::vector<double> tmp4(num_features,0), tmp5(num_features,0),  normalRdn(num_features,0);

	MatrixFeat alpha_users, alpha_movies, tmp, lam, tmp2, covar;

	const coords inner = users.domain().local().upper - users.domain().local().lower;

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.5,1);

	for (unsigned i = 0; i < inner[0]; i++)
		for (unsigned j = 0; j < inner[1]; j++)
			users.ptr[i * inner[1] + j] = distribution(generator);

	for (unsigned i = 0; i < movies.rows(); i++)
	   for (unsigned j = 0; j < movies.cols(); j++)
	   	   movies(i, j) = distribution(generator);

	for (int f = 0; f < num_features; f++)
	{
		mu_u[f] = columnsMean(users,f);
		mu_m[f] = columnsMean(movies,f);
	}

	alpha_users.setIdentity();
	covariance(users,tmp);
	alpha_users = tmp.inverse();

	alpha_movies.setIdentity();
	covariance(movies,tmp);
	alpha_movies = tmp.inverse();

	//constantes

	int beta=2;
	int decay = -1; //constant decay: Niu et al, Hogwild!: A lock-free approach to parallelizing stochastic gradient descent, NIPS 2011.
	double initLRate = 0.001;
	double maxLRate = -1;
	double regU = 0.001;
	double regI = 0.001;
	double regB = 0.001;

	int my_users_min = users.domain().local().lower[0];
	int my_users_max = users.domain().local().upper[0];

	//TODO attribuer les films aux processus qui possèdent le plus d'utilisteurs qui ont noté ce film
	// la matrice de départ est inchangée donc on peut la diagonaliser une fois et attribuer les films à échantilloner aux processus qui ont le plus grands nombre d'utilisateurs qui ont marqué ce film
	// diagonaliser la matrice et attribuer à chaque processus va garder les indices des films à faire

	std::vector<int> movie_bids(world().nprocs,0); //contains the indices of the movies that this process will sample
	std::vector<int> assigned_movies(num_movies); //contains the indice of the process that will sample the movie assigned_movies[0] = 1 the process number 1 will sample the movies number 0 ...

	double* communication_buffer = (double*) malloc (num_features * sizeof(double));
	double send[num_features], recv[num_features]; //first element is the movie index to update

	MatrixXd MM;
	std::vector<double> rr;

	for (int u = 0; u < num_movies; u++)
	{
		int max = 0;

		for(int i = 0; i < movie_indices.size(); i++)
			if (movie_indices[i] == u)
			{
			    for(int k = 0; k < world().nprocs; k++)
				{
					if ( ( d.local(k).lower[0] <= user_indices[i] ) && ( user_indices[i] <= d.local(k).upper[0] ) ) movie_bids[k]++;
					if ( movie_bids[k] > movie_bids[max] ) max = k;
				}
			}

		assigned_movies[u] = max;
	}

	for (int iter = 1; iter <= iterations; iter++)
	{
		sample_user_hyper_parameters(users, num_users, alpha_users, mu_u);

		sample_movie_hyper_parameters(movies, num_movies, alpha_movies, mu_m);

		for (int gibbs = 0; gibbs < 1; gibbs++)
		{
			for (int u = my_users_min; u < my_users_max; u++) //do it for my users range
			{
				// Get the list of items rated by user uu:
				// features of items rated by user uu:

				int count = 0, idx = 0;

				for (int i = 0; i < user_indices.size(); i++)
					if (user_indices[i] == u)	count++;

				if (MM.rows() != count) MM.resize(count,num_features);
				if (rr.size() != count) rr.resize(count);

				for (int i = 0; i < user_indices.size(); i++)
				{
					if (user_indices[i] == u)
					{
						rr[idx] = rate_values[i] - globalMean;

						for (int f = 0; f < num_features; f++)
							MM(idx, f) = movies(movie_indices[i],f) ;

						idx++;
					}
				}

				if( MM.rows() != 0 )
				{
					for (int i = 0; i < MM.cols(); i++) //tM.rows()
					{
						for (int j = 0; j < MM.cols(); j++) //M.cols()
						{
							double product1 = 0; // tmp2 = tMM * MM

							for (int k = 0; k < MM.rows(); k++) //tM.cols()
								product1 += MM(k,i) * MM(k,j);

							tmp2(i,j) = product1;
						}

						double product2 = 0; // tmp3 = tMM * rr

						for (int j = 0; j < MM.rows(); j++) //tM.cols()
							product2 += MM(j,i) * rr[j];

						tmp3[i] = product2;
					}

					for (int i = 0 ; i < num_features; i++)
						for (int j = 0 ; j < num_features; j++)
							tmp(i,j) = alpha_users(i,j) + (tmp2(i,j) * beta);

					covar = tmp.inverse();

					mult(alpha_users, mu_u, tmp4);

					for (int j = 0 ; j < num_features; j++)
						tmp5[j] = tmp3[j] * beta + tmp4[j];

					mult(covar, tmp5, mean);

					LLT<MatrixXd> lltOfcovar(covar);
					lam = lltOfcovar.matrixL();
					lam.transposeInPlace();

					for (int i = 0; i < num_features; i++)
					{
						double product = 0;

						for (int j = 0; j < num_features; j++)
						{
							normalRdn[j] = distribution(generator);
							product += lam(i,j) * normalRdn[j];
						}

						users.ptr[ ( u - my_users_min ) * num_features + i] = mean[i] + product ;
					}
				}
			}

		    printf("Process %d Finished sampling users\n" , world().procid);

			//sampling movies

			for (int j = 1; j < num_movies; j++) //movies 0 does'nt exist
			{
				int receiver_idx = -1 , sender_idx;

				int count = 0, idx = 0;

				for (int i = 0; i < movie_indices.size(); i++)
					if (movie_indices[i] == j)	count++;

				if ( MM.rows() != count ) MM.resize(count,num_features);
				if ( rr.size() != count ) rr.resize(count);

				if (assigned_movies[j] == world().procid) // process this only if this movie is assigned to me
				{
					if (count != 0) //this movie has users that rated it
					{
						for (int i = 0; i < movie_indices.size(); i++)
							if (movie_indices[i] == j)
							{
								rr[idx] = rate_values[i] - globalMean;

								if( (my_users_min <= user_indices[i]) && (user_indices[i] <= my_users_max) ) //if I have this user in my range of users
								{
									for (int f = 0; f < num_features; f++)
										MM(idx, f) = users.ptr[(user_indices[i] - my_users_min) * num_features + f];
								}
								else //get this user's information from the process that owns it using get()
								{
									printf("%d is here???? %d \n", world().procid , i);

									coords_range missing_user;

									missing_user.lower[0] = user_indices[i];
									missing_user.upper[0] = user_indices[i] + 1;

									missing_user.lower[1] = 0;
									missing_user.upper[1] = num_features;

								//	world().sync();

									users.get(missing_user,communication_buffer);

								//	world().sync();

									for (int f = 0; f < num_features; f++)
										MM(idx, f) = communication_buffer[f];
								}

								idx++;
							}

						printf("%d is here!!!!!!!!!!\n", world().procid );

						for (int i = 0; i < MM.cols(); i++) //tM.rows()
						{
							for (int j = 0; j < MM.cols(); j++) //M.cols()
							{
								double product1 = 0; // tmp2 = tMM * MM

								for (int k = 0; k < MM.rows(); k++) //tM.cols()
									product1 += MM(k,i) * MM(k,j);

								tmp2(i,j) = product1;
							}

							double product2 = 0; // tmp3 = tMM * rr

							for (int j = 0; j < MM.rows(); j++) //tM.cols()
								product2 += MM(j,i) * rr[j];

							tmp3[i] = product2;
						}

						printf("%d is here\n", world().procid );

						for (int i = 0 ; i < num_features; i++)
							for (int l = 0 ; l < num_features; l++)
								tmp(i,l) = alpha_movies(i,l) + (tmp2(i,l) * beta);

						covar = tmp.inverse();

						mult(alpha_movies, mu_m, tmp4);

						for (int l = 0 ; l < num_features; l++)
							tmp5[l] = tmp3[l] * beta + tmp4[l];

						mult(covar, tmp5, mean);

						LLT<MatrixXd> lltOfcovar(covar);
						lam = lltOfcovar.matrixL();
						lam.transposeInPlace();

						for (int i = 0; i < num_features; i++)
						{
							double product = 0;

							for (int l = 0; l < num_features; l++)
							{
								normalRdn[l] = distribution(generator);
								product += lam(i,l) * normalRdn[l];
							}

							communication_buffer[i] = mean[i] + product;

							movies(j,i) = mean[i] + product;
						}

						vector<int> flags(world().nprocs, 0); //notify a process about a movie only once, this is used in the case a process have more than one user that have rated that movie

						printf("%d is here\n", world().procid );

						for (int i = 0; i < movie_indices.size(); i++)
							if ( movie_indices[i] == j && ( (user_indices[i] < my_users_min) || (user_indices[i] > my_users_max) ) )
							{
								for(int k = 0; k < world().nprocs; k++)
									if ( ( d.local(k).lower[0] <= user_indices[i] ) && ( user_indices[i] <= d.local(k).upper[0] ) && flags[k] == 0 )
									{
									   receiver_idx = k;  // k it is the indice of the process that I have to notify

									   for (int i = 0; i < num_features; i++) send[i] = communication_buffer[i];

									   flags[k] = 1;

									   printf("%d sending info for movie %d which assigned to %d\n", world().procid, j, assigned_movies[j]);

									   MPI_Send(send, num_features, MPI_DOUBLE, receiver_idx , 0, MPI_COMM_WORLD);

									   printf("%d waiting info for movie %d which assigned to %d\n", world().procid, j, assigned_movies[j]);
									}
							}
					}
				}
				else //the movie is not assigned to me but I may have users that rated it
				{
					if (count != 0) //this movie has users that rated it
					{
						for (int i = 0; i < movie_indices.size(); i++)
						   if ( movie_indices[i] == j && ( (my_users_min <= user_indices[i]) && (user_indices[i] < my_users_max) ) )
						   {
								sender_idx = assigned_movies[j];

								MPI_Status status;

								printf("%d waiting info for movie %d which assigned to %d\n", world().procid, j, assigned_movies[j]);

								MPI_Recv(recv, num_features, MPI_DOUBLE, sender_idx , 0, MPI_COMM_WORLD, &status);

								printf("%d received info for movie %d which assigned to %d\n", world().procid, j, assigned_movies[j]);

								for (int i = 0; i < num_features; i++) movies(recv[0], i) = recv[i];
							}
					}
				}

				printf("Process %d Finished sampling the movies %d\n", world().procid ,j);
			}

			printf("Process %d Finished sampling the movies \n", world().procid );

		} // end of gibbs

		int my_rates = 0;

		long double errs = 0, loss = 0;
		long double last_errs , last_loss, lRate;

		for (int i = 0 ; i < user_indices.size(); i++)
		{
			if ( (my_users_min <= user_indices[i]) && (user_indices[i] <= my_users_max))
			{
					int user = user_indices[i];
					int movie = movie_indices[i];
					double rate = rate_values[i];

					long double pred = 0.0;

					for (int j = 0; j < num_features ; j++)
						pred += users.ptr[ (user-my_users_min) * num_features + j] * movies(movie, j);

					pred += globalMean;

					long double euj = rate - pred;

					errs += euj * euj;

					//cerr << "pred " << pred << " rate " << rate << " euj " << euj << " errs " << errs << endl;

					my_rates++;
			}
		}

		errs *= 0.5;

		// check if converged

		bool cond1 = (errs < 1e-5);
		bool cond2 = (last_errs >= errs && last_errs - errs < 1e-5);
		bool converged = cond1 || cond2;

		// if not converged, update learning rate

		if (!converged && lRate >= 0)
		{
			if (last_loss != 0.0)
			{
				if (abs(last_loss) > abs(loss))	lRate *= 1.05;
				else	lRate *= 0.5;
			}
			else if (decay > 0 && decay < 1) lRate *= decay;
			else if (decay == 0) lRate = initLRate / (1 + initLRate * ((regU + regI) / 2.0) * iter);

			// limit to max-learn-rate after update
			if (maxLRate > 0 && lRate > maxLRate) lRate = maxLRate;
		}

		cerr << "iter " << iter << " errs = " << errs << " last_errs = " << last_errs << " delta_errs = " << last_errs - errs
			 << " m_rmse " << sqrt(errs/my_rates) << " learn_rate = " << lRate << endl;

		last_errs = errs;

		if (converged == true)	break;
	}
}

