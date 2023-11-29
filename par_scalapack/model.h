#ifndef _MODEL_H
#define _MODEL_H

#define BLOCK_CYCLIC_2D 1
#define DLEN_ 9 
#define DTYPE_ 1
#define CTXT_ 2
#define M_ 3 
#define N_ 4 
#define MB_ 5
#define NB_ 6
#define RSRC_ 7 
#define CSRC_ 8
#define LLD_ 9

extern void Cblacs_get(long icontxt, long what, long *val);
extern void Cblacs_setup(long * rank, long * nprocs);
extern void Cblacs_gridinit(long * ictx, const char * order, long nprow, long npcol);
extern void Cblacs_gridinfo(long ictx, long * nprow, long * npcol, long * myrow, long * mycol);
extern void Cblacs_gridexit( long ictx );
extern void Cblacs_exit( long doneflag );

extern void* mr2d_malloc(long size);

/*
* This function computes the local number of rows or columns of a 
* block-cyclically distributed matrix contained in a process row or 
* process column, respectively, indicated by the calling sequence argument iproc.
*/
long numroc_(long* N, long* NB, long* IPROC, long *ISRCPROC, long* NPROCS );

/*
* This function computes the process row or column index of a 
* global element of a block-cyclically distributed matrix.
*/
long indxg2p_( long* INDXGLOB, long* NB, long* IPROC, long* ISRCPROC, long* NPROCS  );

/*
* The pdgemr2d routine copies the indicated matrix or submatrix 
* of A to the indicated matrix or submatrix of B. It provides a 
* truly general copy from any block cyclicly-distributed matrix 
* or submatrix to any other block cyclicly-distributed matrix or submatrix.
*/
void pdgemr2d_( long* m, long* n, double* a, long* ia, long* ja, long* desca, double* b, long* ib, long* jb, long* descb, long* ictxt);

/*
* These subroutines compute the QR factorization of a general matrix A, 
* where, in this description:
* A represents the global general submatrix Aia:ia+m-1, ja:ja+n-1 to be factored.
* For PDGEQRF, Q is an orthogonal matrix.
* For m >= n, R is an upper triangular matrix.
* If m = 0 or n = 0, no computation is performed and the subroutine returns 
* after doing some parameter checking.
*/
void pdgeqrf_(long* M, long *N, double* A, long* IA, long *JA, long* DESCA, double *TAU, double *WORK, long* LWORK, long *INFO);

/* 
* This subroutine initializes a type-1 array descriptor with error checking 
*/
void descinit_ (long *desc, const long *m, const long *n, const long *mb, const long *nb, const long *irsrc, const long *icsrc, const long *ictxt, const long *lld, long *info);


#endif