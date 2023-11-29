#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <mpi.h>
#include <time.h>
#include "harmonics.h"
#include "model.h"

#define BLOCK_SIZE 64

int lmax = -1;
long npoint;
int spacing = 1;
char * data_filename;
char * model_filename;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--data FILENAME              input file containing experimental data points\n");
        printf("--model FILENAME             output file containing the model\n");
        printf("--npoint N                   number of points to read\n");
        printf("--lmax N                     order of the model\n");
		printf("--spacing N					 spacing between two selected points\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[6] = {
                {"data", required_argument, NULL, 'd'},
                {"npoint", required_argument, NULL, 'n'},
                {"lmax", required_argument, NULL, 'l'},
                {"model", required_argument, NULL, 'm'},
				{"spacing", required_argument, NULL, 's'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'd':
                        data_filename = optarg;
                        break;
                case 'm':
                        model_filename = optarg;
                        break;
                case 'n':
                        npoint = atoi(optarg);
                        break;
				case 's' :
						spacing = atoi(optarg);
						break;
                case 'l':
                        lmax = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        if (data_filename == NULL || model_filename == NULL || lmax < 0 || npoint <= 0)
                usage(argv);
}

/**************************** LINEAR ALGEBRA *********************************/


/*
 * Return the euclidean norm of x[0:n] using tricks for a greater precision
 */
double norm(int n, double const *x)
{
   double rdwarf = 3.834e-20, rgiant = 1.304e19;
   double s1 = 0, s2 = 0, s3 = 0;
   double x1max = 0, x3max = 0;
   double agiant = rgiant / ((double) n);
   for (int i = 0; i < n; i++) {
       double xabs = fabs(x[i]);
       if (xabs > rdwarf && xabs < agiant) {  // sum for intermediate components
           s2 += xabs * xabs;
           continue;
       }
       if (xabs <= rdwarf) {                         // sum for small components 
           if (xabs > x3max) {
               double d3 = x3max / xabs;
               s3 = 1 + s3 * (d3 * d3);
               x3max = xabs;
               continue;
           }
           if (xabs != 0) {
               double d4 = xabs / x3max;
               s3 += d4 * d4;
           }
           continue;
       }
       if (xabs <= x1max) {                          // sum for large components
           double d2 = xabs / x1max;
           s1 += d2 * d2;
           continue;
       }
       double d1 = x1max / xabs;
       s1 = 1 + s1 * (d1 * d1);
       x1max = xabs;
   }
   if (s1 == 0) {                                         // calculation of norm
       if (s2 == 0)
           return x3max * sqrt(s3);
       if (s2 >= x3max)
           return sqrt(s2 * (1 + x3max / s2 * (x3max * s3)));
       if (s2 < x3max)
           return sqrt(x3max * (s2 / x3max + x3max * s3));
   }
   return x1max * sqrt(s1 + s2 / x1max / x1max);
}

/*
 * Apply a real elementary reflector H to a real m-by-n matrix C. H is
 * represented in the form
 *
 *       H = I - tau * v * v**T
 *
 * where tau is a real scalar and v is a real vector.
 *
 * C is a 2D array of dimension (ldc, n).  On exit, C is overwritten with H*C.
 * It is required that ldc >= m.
 */
void multiply_householder(int m, int n, double *v, double tau, double *c, int ldc)
{
	for (int j = 0; j < n; j++) {
		double sum = 0;
		for (int i = 0; i < m; i++)
			sum += c[j * ldc + i] * v[i];
		for (int i = 0; i < m; i++)
			c[j * ldc + i] -= tau * v[i] * sum;
	}
}


/*
 * Overwrite vector c with transpose(Q) * c where Q is a
 * real m-by-m orthogonal matrix defined as the product of k elementary
 * reflectors
 *
 *       Q = H(1) * H(2) ... H(k)
 * 
 * A is a 2D array of dimension (m, k), which contains a QR factorisation
 * computed by QR_factorize().  A is not modified.
 *
 * tau is an array of dimension k. tau[i] must contain the scalar factor of the
 * elementary reflector H(i), as returned by QR_factorize().  tau is read-only.
 *
 * c is a vector of dimension m.  On exit, c is overwritten by transpose(Q)*c.
 */
void multiply_Qt(int m, int k, double * A, double * tau, double * c)
{
	for (int i = 0; i < k; i++) {
		/* Apply H(i) to A[i:m] */
		double aii = A[i + i * m];
		A[i + i * m] = 1;
		multiply_householder(m-i, 1, &A[i*m + i], tau[i], &c[i], m);
		A[i + i * m] = aii;
	}
}

/*
 * Solve the triangular linear system U*x == b
 *
 * U is a 2D array of dimension (ldu, n) with non-zero diagonal entries. Only
 * the upper-triangle is read by this function. b and x are n element vectors.
 * On exit, b is overwritten with x.
 */
void triangular_solve(int n, const double *U, int ldu, double *b) 
{
	for (int k = n - 1; k >= 0; k--) {
           b[k] /= U[k * ldu + k];
           for (int i = 0; i < k; i++)
               b[i] -= b[k] * U[i + k*ldu];
    }
}

void print_matrix(int m, int n, double* A) {
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			printf("%.1lf\t", A[j * m + i]);
		} printf("\n");
	} printf("\n");
}


/*****************************************************************************/

long max(long a, long b) {
	return (a > b) ? a : b; 
}

long isqrt(long a) {
	return (long) sqrt((double)a);
}

/* Computes the integers p and q such thtat a = p * q and the interval between p and q is minimal 
*  If a is a sqare, p = q = sqrt(a)
*/

void approx_sqrt(long a, long* p, long* q) {
	int sa = isqrt(a);
	if(sa * sa == a) {
		*p = sa;
		*q = sa;
		return;
	} else {
		int tsa = sa;
		while(a % tsa != 0) tsa--;
		*p = tsa;
		tsa = sa+1;
		while(a % tsa != 0) tsa++;
		*q = tsa;
		return;
	}
}

int main(int argc, char ** argv)
{
	long i_one = 1, i_zero = 0;
	//double zero=0.0E+0, one=1.0E+0;

	long rank = 0, nprocs = 0;

    /* BLACS setup --- obtain system context */
	long myrow, mycol;
    long sysctx;

    Cblacs_get(0, 0, &sysctx);
	Cblacs_setup(&rank, &nprocs);

	long nprow, npcol;
	approx_sqrt(nprocs, &nprow, &npcol);

	long root = 0;
	MPI_Status *status = NULL;
	assert(nprocs % nprow == 0); assert(nprocs % npcol == 0);
	

    /* obtain BLACS grid context */
    long ctx = sysctx;

	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);
	Cblacs_gridinfo(ctx, &nprow, &npcol, &myrow, &mycol);

	double *A = NULL, *A_dist = NULL;
	
	process_command_line_options(argc, argv);

	long nvar = (lmax + 1) * (lmax + 1);
	double start = 0, FLOP = 0, FLOPS = 0, t = 0;

	struct data_points data;
	struct spherical_harmonics model;

	npoint /= spacing;

	if(rank == 0) {
		
		printf("Linear Least Squares with dimension %ld x %ld\n", npoint, nvar);
		if (nvar > npoint)
			errx(1, "not enough data points");

		long matrix_size = sizeof(double) * nvar * npoint;
		char hsize[16];
		human_format(hsize, matrix_size);
		printf("Matrix size: %sB\n", hsize);

		A = malloc(matrix_size);
		if (A == NULL)
			err(1, "cannot allocate matrix");

		double * P = malloc((lmax + 1) * (lmax + 1) * sizeof(*P));
		double * v = malloc(npoint * sizeof(*v));
		if (P == NULL || v == NULL)
			err(1, "cannot allocate data points\n");

		printf("Reading data points from %s\n", data_filename);
		load_data_points(data_filename, npoint, spacing, &data);
		printf("Successfully read %ld data points\n", npoint);

		printf("Building matrix\n");
		setup_spherical_harmonics(lmax, &model);


		for (long i = 0; i < npoint; i++) {
			computeP(&model, P, sin(data.phi[i]));
			//printf("%ld\n", i);
			for (long l = 0; l <= lmax; l++) {
				/* zonal term */
				A[i + npoint * CT(l, 0)] = P[PT(l, 0)];

				/* tesseral terms */
				for (long m = 1; m <= l; m++) {
					//printf("%d\n", m);
					A[i + npoint * CT(l, m)] = P[PT(l, m)] * cos(m * data.lambda[i]);
					A[i + npoint * ST(l, m)] = P[PT(l, m)] * sin(m * data.lambda[i]);
				}
			}
		}

		FLOP = 2. * nvar * nvar * npoint;
		char hflop[16];
		human_format(hflop, FLOP);
		printf("Least Squares (%sFLOP)\n", hflop);
		start = wtime();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("coucou\n");
	long m = npoint;
	long n = nvar;
	assert(m >= n);
	//assert(m % nprow == 0); assert(n % npcol == 0);

	long mb = BLOCK_SIZE;
	long nb = BLOCK_SIZE;
	long info = 0;
	long ia = 1;//(myrow*mb) + 1;
	long ja = 1;//(mycol*nb) + 1;
	long ib = 1;
	long jb = 1;
	//long ia = myrow, ja = mycol;
	long rsrc = i_zero, csrc = i_zero;
	long descA[9], descA_dist[9]; 
	
	double *tau = malloc(n * sizeof(double));
	
	long mp = numroc_( &m, &mb, &myrow, &rsrc, &nprow );
	long nq = numroc_( &n, &nb, &mycol, &csrc, &npcol );
	
	//printf ("#%d : myrow = %d, mycol = %d, m = %d, n = %d, mb = %d, nb = %d, mp = %d, nq = %d, ia = %d, ja = %d\n", rank, myrow, mycol, m, n, mb, nb, mp, nq, ia, ja);

	/* Computing lwork */

	long iroff = (ia-1) % mb,  icoff = (ja-1) % nb;
    long iarow = indxg2p_( &ia, &mb, &myrow, &rsrc, &nprow );
    long iacol = indxg2p_( &ja, &nb, &mycol, &csrc, &npcol );
	long t1 = m+iroff, t2 =  n+icoff;
    long mp0 = numroc_( &t1, &mb, &myrow, &iarow, &nprow );
    long nq0 = numroc_( &t2, &nb, &mycol, &iacol, &npcol );

	long lwork = nb * ( mp0 + nq0 + nb );
	double *work = malloc(lwork * sizeof(double));

	/* lwork computed */

	//printf("coucou1\n");
	A_dist = malloc(mp * nq * sizeof(double));

	long lld = max( numroc_( &m, &m, &myrow, &i_zero, &nprow ), 1 );
	long lld_dist = max( mp, 1 );

	
	descinit_( descA,		&m, &n, &m	, &n	, &rsrc, &csrc, &ctx,	 	&lld		, &info );
	descinit_( descA_dist, 	&m, &n, &mb	, &nb	, &rsrc, &csrc, &ctx, 		&lld_dist	, &info );

	/* 2D bloc cyclic matrix scattering */
	//printf("rank = %ld. m = %ld, n = %ld\n", rank, m, n);
	pdgemr2d_(&m, &n, A, &ia, &ja, descA, A_dist, &ib, &jb, descA_dist, &ctx);
	/* Distributed QR factorization */
	//printf("Data distribution done\n");
	pdgeqrf_(&m, &n, A_dist, &i_one, &i_one, descA_dist, tau, work, &lwork, &info);
	/* 2D bloc cyclic matrix gathering */
	pdgemr2d_(&m, &n, A_dist, &ib, &jb, descA_dist, A, &ia, &ja, descA, &ctx);


	free(A_dist);
	free(work);
	
	long count = n/npcol;

	/* Sending tau parts to the main process */

	if((rank < npcol) && (rank != 0)) 
		MPI_Send(tau, count, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
	
	/* Receiving tau parts */
	if(rank == root) 
		for(long i = 1; i < npcol; i++) 
			MPI_Recv(&tau[(n * i) / npcol], count, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, status);
		

	if (rank != 0) free(tau);

	Cblacs_gridexit( ctx );
	Cblacs_exit( i_zero );

	if (rank == 0) {
		double *b = data.V;

		multiply_Qt(m, n, A, tau, b);
		triangular_solve(n, A, m, b);   

		t = wtime()  - start;
		FLOPS = FLOP / t;
		char hflops[16];
		human_format(hflops, FLOPS);
		printf("Completed in %.1f s (%s FLOPS)\n", t, hflops);
		fprintf(stderr, "%ld, %ld, %ld, %lf, %lf, %lf\n", nprocs, npoint, nvar, t, FLOP, FLOPS);

		double res = 0;
		for (int j = nvar; j < npoint; j++)
			res += data.V[j] * data.V[j];
		printf("residual sum of squares %g\n", res);


		printf("Saving model in %s\n", model_filename);
		FILE *g = fopen(model_filename, "w");
		if (g == NULL)
			err(1, "cannot open %s for writing\n", model_filename);
		for (int l = 0; l <= lmax; l++) {
			fprintf(g, "%d\t0\t%.18g\t0\n", l, data.V[CT(l, 0)]);
			for (int m = 1; m <= l; m++)
				fprintf(g, "%d\t%d\t%.18g\t%.18g\n", l, m, data.V[CT(l, m)], data.V[ST(l, m)]);
		}
	}
}
