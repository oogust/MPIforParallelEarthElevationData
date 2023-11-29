#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <mpi.h>

#include "harmonics.h"

int lmax = -1;
char * data_filename;
char * model_filename;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--data FILENAME              CSV file containing experimental data points\n");
	printf("--model FILENAME             file containing the model\n");        
        printf("--npoint N                   number of points to read\n");
        printf("--lmax N                     order of the model\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[5] = {
                {"data", required_argument, NULL, 'd'},
                {"lmax", required_argument, NULL, 'l'},
                {"model", required_argument, NULL, 'm'},
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
                case 'l':
                        lmax = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        if (data_filename == NULL || model_filename == NULL || lmax < 0)
                usage(argv);
}


int main(int argc, char ** argv)
{
        int rank, nprocs;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Status *status = NULL;
        int root = 0;

	process_command_line_options(argc, argv);

	/* preparations and memory allocation */
	double * P = malloc((lmax + 1) * (lmax + 1) * sizeof(*P));
	if (P == NULL)
		err(1, "cannot allocate data points\n");

	if (rank == 0) printf("Reading order-%d model from %s\n", lmax, model_filename);
	struct spherical_harmonics model;
	load_spherical_harmonics(model_filename, lmax, &model);
	if (rank == 0) printf("Successfully read model\n");

	/* evaluate the model on corresponding input points */
        FILE *f = fopen(data_filename, "r");
        if (f == NULL)
                err(1, "cannot open %s", data_filename);
        
        if (rank == 0) printf("Reading data points from %s\n", data_filename);
        double max_error = 0;
        double error_l1 = 0;
        double error_l2 = 0;
        int npoint = 0;
        double start = 0.;
        if (rank == root) start = wtime();

        int j = rank;
        int bufferLength = 256;
	char buffer[bufferLength];
        while(j--) fgets(buffer, bufferLength, f);
        npoint += rank;

        for (;;) {
                /* read data point */
                double lambda, phi, v;
                int k = fscanf(f, "%lg %lg %lg", &lambda, &phi, &v);

		j = nprocs;
		while(j--) fgets(buffer, bufferLength, f);
        
                if (k == EOF) {
                        if (ferror(f))
                                err(1, "read error");
                        if (feof(f))
                                break;
                }
                npoint += nprocs;
                if (k != 3)
                        errx(1, "parse error on line %d", npoint);

                /* compute error */
		computeP(&model, P, sin(phi));
                double vmodel = evaluate(&model, P, lambda);
                double error = fabs(v - vmodel);
                if (error > max_error)
                        max_error = error;
                error_l1 += error;
                error_l2 += error * error;
        }
        if(rank == root) {
                for(int i = 1; i < nprocs; i++) {
                        double error_l1_tmp = 0;
                        double error_l2_tmp = 0;
                        double max_error_tmp = 0;

                        MPI_Recv(&error_l1_tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, status);
                        MPI_Recv(&error_l2_tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, status);
                        MPI_Recv(&max_error_tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, status);
                        error_l1 += error_l1_tmp;
                        error_l2 += error_l2_tmp;
                        max_error = (max_error_tmp > max_error) ? max_error_tmp : max_error;
                }
        } else {
                MPI_Send(&error_l1, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
                MPI_Send(&error_l2, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
                MPI_Send(&max_error, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        }
        if(rank == root) {
                printf("Completed in %.1f s\n", wtime() - start);
                fprintf(stderr, "%d, %lf\n", nprocs, wtime() - start);
	        printf("Average error %g meters\n", error_l1 / npoint);
                printf("Max error %g meters\n", max_error);
	        printf("Standard Deviation %g\n", sqrt(error_l2 / npoint));
        }
        MPI_Finalize();
}