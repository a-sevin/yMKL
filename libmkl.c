#include"yapi.h"

// Macro that simplifies calls to y_error() when reporting errors
// to Yorick.
// The do{... } while(0) is there to avoid forgetting { ... } around it.
// __FILE__ and __LINE__ are replaced at compilation by actual filename and file line.
// (note: __FUNCTION__ also exists ..)
#define ERROR2YORICK(fmt, args...) do { char errmsg[128]; sprintf(errmsg, "[%s@%d]: " fmt, __FILE__, __LINE__, ##args); y_error(errmsg); } while(0)

#include<mkl_lapack.h>
#include<mkl_blas.h>
#include<stdlib.h>
#include<stdio.h>


int Y_potri_cpu(int argc) {
/** @brief wrapper routine for mkl_lapack potri method
 *  @param[in] argc : command line arguments
 *  can work as a subroutine (return discarded)
 *    - first   : the matrix to inverse
 */
  if (yarg_subroutine()) {
    int yType;
    long ntot;
    long dims[Y_DIMSIZE];
    MKL_INT info=0;
    char uplo='L';
    MKL_INT one = 1;
    MKL_INT row;

    void *h_mat=ygeta_any(argc - 1, &ntot, dims, &yType);
    MKL_INT n=dims[1];
    MKL_INT m=dims[2];

    if (yType == Y_FLOAT) {
      spotrf( &uplo, &n, (float*) h_mat, &n, &info );
      if(info) ERROR2YORICK("Cholesky failed");
      spotri( &uplo, &n, (float*) h_mat, &n, &info );
      if(info) ERROR2YORICK("Cholesky inverse failed");

      // Copy lower part to the upper
      MKL_INT size=n-1;
      float *matL = (float*)h_mat+1;
      float *matU = (float*)h_mat+n;
      do{
	scopy(&size, matL, &one, matU, &n);
	size--;
	matL+=n+1;
	matU+=n+1;
      } while (size>0);
    } else if (yType == Y_DOUBLE) {
      dpotrf( &uplo, &n, (double*) h_mat, &n, &info);
      if(info) ERROR2YORICK("Cholesky failed");
      dpotri( &uplo, &n, (double*) h_mat, &n, &info );
      if(info) ERROR2YORICK("Cholesky inverse failed");

      // Copy lower part to the upper
      MKL_INT size=n-1;
      double *matL = (double*)h_mat+1;
      double *matU = (double*)h_mat+n;
      do{
	dcopy(&size, matL, &one, matU, &n);
	size--;
	matL+=n+1;
	matU+=n+1;
      } while (size>0);
    } else {
      y_error("carma_potri not implemented for this type");
    }
  } else {
    y_error("carma_potri must be call as a subroutine");
  }
  return 0;
}

void fComputeMatrixEigenVectors(const MKL_INT matrix_rank, float *matrix,float *EV)  
{  
  char JOBZ='V';  
  char UPLO='L';  
  int info=0;  
  MKL_INT n=matrix_rank;  
  MKL_INT lda=matrix_rank;  
  MKL_INT lwork=-1;  
  MKL_INT liwork=-1;  
  MKL_INT iwkopt;  
  MKL_INT *iwork;  
  float wkopt;  
  float *work;  
  fsyevd(&JOBZ,&UPLO,&n,matrix,&lda,EV,&wkopt,&lwork,&iwkopt,&liwork,&info);  
  if(info>0)  ERROR2YORICK("The algorithm failed to query workspace.");
  lwork=(MKL_INT)wkopt;  
  liwork=iwkopt;  
  work=malloc(lwork*sizeof(float));  
  iwork=malloc(liwork*sizeof(MKL_INT));
  fsyevd(&JOBZ,&UPLO,&n,matrix,&lda,EV,work,&lwork,iwork,&liwork,&info);  
  if(info>0)  ERROR2YORICK("The algorithm failed to compute eigenvalues.");
  free(work);  
  free(iwork);  
}  

void dComputeMatrixEigenVectors(const MKL_INT matrix_rank, double *matrix,double *EV)  
{  
  char JOBZ='V';  
  char UPLO='L';  
  MKL_INT info=0;  
  MKL_INT n=matrix_rank;  
  MKL_INT lda=matrix_rank;  
  MKL_INT lwork=-1;  
  MKL_INT liwork=-1;  
  MKL_INT iwkopt;  
  MKL_INT *iwork;  
  double wkopt;  
  double *work;  
  dsyevd(&JOBZ,&UPLO,&n,matrix,&lda,EV,&wkopt,&lwork,&iwkopt,&liwork,&info);  
  if(info>0)  ERROR2YORICK("The algorithm failed to query workspace.");
  lwork=(MKL_INT)wkopt;  
  liwork=iwkopt;  
  work=malloc(lwork*sizeof(double));  
  iwork=malloc(liwork*sizeof(MKL_INT));
  dsyevd(&JOBZ,&UPLO,&n,matrix,&lda,EV,work,&lwork,iwork,&liwork,&info);  
  if(info>0)  ERROR2YORICK("The algorithm failed to compute eigenvalues.");
  free(work);  
  free(iwork);  
}  

int Y_syevd_cpu(int argc) {
/** @brief wrapper routine for mkl_lapack potri method
 *  @param[in] argc : command line arguments
 *  can work as a subroutine (return discarded)
 *    - first   : the matrix to inverse
 */
  if (yarg_subroutine()) {
    int yType;
    long ntot;
    long dims[Y_DIMSIZE];
    char uplo='L';
    MKL_INT one = 1;
    MKL_INT row;

    void *h_mat=ygeta_any(argc - 1, &ntot, dims, &yType);
    void *h_ev=ygeta_any(argc - 2, &ntot, dims, &yType);
    MKL_INT n=dims[1];

    if (yType == Y_FLOAT) {
      fComputeMatrixEigenVectors(n, (float*)h_mat, (float*)h_ev);
    } else if (yType == Y_DOUBLE) {
      dComputeMatrixEigenVectors(n, (double*)h_mat, (double*)h_ev);
    } else {
      y_error("carma_potri not implemented for this type");
    }
  } else {
    int yType;
    long ntot;
    long dims[Y_DIMSIZE];
    char uplo='L';
    MKL_INT one = 1;
    MKL_INT row;

    void *h_mat=ygeta_any(argc - 1, &ntot, dims, &yType);
    dims[0]=1;
    MKL_INT n=dims[1];

    if (yType == Y_FLOAT) {
      float *h_ev=ypush_f(dims);
      fComputeMatrixEigenVectors(n, (float*)h_mat, h_ev);
    } else if (yType == Y_DOUBLE) {
      double *h_ev=ypush_d(dims);
      dComputeMatrixEigenVectors(n, (double*)h_mat, h_ev);
    } else {
      y_error("carma_potri not implemented for this type");
    }
  }
  return 0;
}


void Y_gemm_cpu(int argc)
/** @brief wrapper routine for mkl_blas gemm method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first   : (1) the C matrix  / (2) the A matrix
 *    - second  : (1) the A matrix  / (2) the B matrix
 *    - third   : (1) the B matrix
 *    - fourth  : (1) the alpha coeff
 *    - fifth   : (1) the beta coeff
 *    - sixth   : (1) the opA
 *    - seventh : (1) the opB
 *  in case (2) the destination is pushed on the stack as a yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
  int yType=yarg_typeid(argc - 1);
  if (yarg_subroutine()) {
    char opA = 'n';
    if (argc > 3)
      opA = ygets_c(argc - 4);
    char opB = 'n';
    if (argc > 4)
      opB = ygets_c(argc - 5);

    if (yType == Y_FLOAT) {
      float alpha = 1.0f;
      if (argc > 5)
        alpha = ygets_f(argc - 6);
      float beta = 0.0f;
      if (argc > 6)
        beta = ygets_f(argc - 7);

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE],dimC[Y_DIMSIZE];
      float *matA = ygeta_f(argc - 1, &ntot, dimA);
      float *matB = ygeta_f(argc - 2, &ntot, dimB);
      float *matC = ygeta_f(argc - 3, &ntot, dimC);

      MKL_INT m = dimC[1];
      MKL_INT n =  dimC[2];
      MKL_INT k = (opA == 'n') ? dimA[2] : dimA[1];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimC[1];
      sgemm(&opA, &opB, &m, &n, &k, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);
    }
    if (yType == Y_DOUBLE) {
      double alpha = 1.0;
      if (argc > 5)
        alpha = ygets_d(argc - 6);
      double beta = 0.0;
      if (argc > 6)
        beta = ygets_d(argc - 7);

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE],dimC[Y_DIMSIZE];
      double *matA = ygeta_d(argc - 1, &ntot, dimA);
      double *matB = ygeta_d(argc - 2, &ntot, dimB);
      double *matC = ygeta_d(argc - 3, &ntot, dimC);

      MKL_INT m = dimC[1];
      MKL_INT n =  dimC[2];
      MKL_INT k = (opA == 'n') ? dimA[2] : dimA[1];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimC[1];
      dgemm(&opA, &opB, &m, &n, &k, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);
    }
  } else {
    // called as a function : need to allocate space
    char opA = 'n';
    if (argc > 2)
      opA = ygets_c(argc - 3);
    char opB = 'n';
    if (argc > 3)
      opB = ygets_c(argc - 4);

    if (yType == Y_FLOAT) {
      float alpha = 1.0f;
      if (argc > 4)
        alpha = ygets_f(argc - 5);
      float beta = 0.0f;
      if (argc > 5)
        beta = ygets_f(argc - 6);

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE],dimC[Y_DIMSIZE];
      float *matA = ygeta_f(argc - 1, &ntot, dimA);
      float *matB = ygeta_f(argc - 2, &ntot, dimB);

      dimC[0] = 2;
      dimC[1] = (opA == 'n') ? dimA[1] : dimA[2];
      dimC[2] = (opB == 'n') ? dimB[2] : dimB[1];
      float *matC = ypush_f(dimC);

      MKL_INT m = dimC[1];
      MKL_INT n =  dimC[2];
      MKL_INT k = (opA == 'n') ? dimA[2] : dimA[1];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimC[1];
      sgemm(&opA, &opB, &m, &n, &k, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);

    } else if (yType == Y_DOUBLE) {
      double alpha = 1.0;
      if (argc > 4)
        alpha = ygets_d(argc - 5);
      double beta = 0.0;
      if (argc > 5)
        beta = ygets_d(argc - 6);

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE],dimC[Y_DIMSIZE];
      double *matA = ygeta_d(argc - 1, &ntot, dimA);
      double *matB = ygeta_d(argc - 2, &ntot, dimB);

      dimC[0] = 2;
      dimC[1] = (opA == 'n') ? dimA[1] : dimA[2];
      dimC[2] = (opB == 'n') ? dimB[2] : dimB[1];
      double *matC = ypush_d(dimC);

      MKL_INT m = dimC[1];
      MKL_INT n =  dimC[2];
      MKL_INT k = (opA == 'n') ? dimA[2] : dimA[1];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimC[1];
      dgemm(&opA, &opB, &m, &n, &k, &alpha, matA,
	    &lda, matB,&ldb, &beta, matC, &ldc);
    }
  }
}

void Y_symm_cpu(int argc)
/** @brief wrapper routine for mkl_blas symm method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first   : (1) the C matrix  / (2) the A matrix
 *    - second  : (1) the A matrix  / (2) the B matrix
 *    - third   : (1) the B matrix
 *    - fourth  : (1) the alpha coeff
 *    - fifth   : (1) the beta coeff
 *    - sixth   : (1) the opA
 *    - seventh : (1) the opB
 *  in case (2) the destination is pushed on the stack as a yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
  int yType=yarg_typeid(argc - 1);
  if (yarg_subroutine()) {
    char side = 'L'; //or 'R'
    if (argc > 3)
      side = ygets_c(argc - 4);

    if (yType == Y_FLOAT) {
      float alpha = 1.0f;
      if (argc > 4)
        alpha = ygets_f(argc - 5);
      float beta = 0.0f;
      if (argc > 5)
        beta = ygets_f(argc - 6);

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE],dimC[Y_DIMSIZE];
      float *matA = ygeta_f(argc - 1, &ntot, dimA);
      float *matB = ygeta_f(argc - 2, &ntot, dimB);
      float *matC = ygeta_f(argc - 3, &ntot, dimC);

      MKL_INT m = dimC[1];
      MKL_INT n =  dimC[2];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimC[1];
      ssymm(&side, "L", &m, &n, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);
    }
    if (yType == Y_DOUBLE) {
      double alpha = 1.0f;
      if (argc > 4)
        alpha = ygets_d(argc - 5);
      double beta = 0.0f;
      if (argc > 5)
        beta = ygets_d(argc - 6);

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE],dimC[Y_DIMSIZE];
      double *matA = ygeta_d(argc - 1, &ntot, dimA);
      double *matB = ygeta_d(argc - 2, &ntot, dimB);
      double *matC = ygeta_d(argc - 3, &ntot, dimC);

      MKL_INT m = dimC[1];
      MKL_INT n =  dimC[2];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimC[1];
      dsymm(&side, "L", &m, &n, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);
    }
  } else {
    // called as a function : need to allocate space
    char side = 'L'; // or'R'
    if (argc > 2)
      side = ygets_c(argc - 3);

    if (yType == Y_FLOAT) {
      float alpha = 1.0f;
      if (argc > 3)
        alpha = ygets_f(argc - 4);
      float beta = 0.0f;

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE];
      float *matA = ygeta_f(argc - 1, &ntot, dimA);
      float *matB = ygeta_f(argc - 2, &ntot, dimB);

      float *matC = ypush_f(dimB);

      MKL_INT m = dimB[1];
      MKL_INT n =  dimB[2];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimB[1];
      ssymm(&side, "L", &m, &n, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);

    } else if (yType == Y_DOUBLE) {
      double alpha = 1.0;
      if (argc > 3)
        alpha = ygets_d(argc - 4);
      double beta = 0.0;

      long ntot;
      long dimA[Y_DIMSIZE],dimB[Y_DIMSIZE];
      double *matA = ygeta_d(argc - 1, &ntot, dimA);
      double *matB = ygeta_d(argc - 2, &ntot, dimB);

      double *matC = ypush_d(dimB);

      MKL_INT m = dimB[1];
      MKL_INT n =  dimB[2];

      MKL_INT lda = dimA[1];
      MKL_INT ldb = dimB[1];
      MKL_INT ldc = dimB[1];
      dsymm(&side, "L", &m, &n, &alpha, matA,
    		  &lda, matB,&ldb, &beta, matC, &ldc);
    }
  }
}

