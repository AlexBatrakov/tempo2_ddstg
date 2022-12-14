#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "tempo2.h"
#include "T2accel.h"
#include "TKmatrix.h"
#include <assert.h>

/*
 * Solves L.x=b by forward substitusion for triangular L.
 */
void TKforwardSubVec( double **L, double *b,int ndata,double *x) {

    for (int i=0; i < ndata; ++i){
        double s=0;
        for (int j = 0; j < i; ++j) {
            s+=L[i][j]*x[j];
        }
        x[i] = (b[i]-s)/L[i][i];
    }
}

/*
 * Solves L.X=M by forward substitusion for triangular L.
 */
void TKforwardSub( double **L, double **M,int ndata,int np,double **X) {
    for (int ip=0; ip < np; ++ip){
        for (int i=0; i < ndata; ++i){
            double s=0;
            for (int j = 0; j < i; ++j) {
                s+=L[i][j]*X[j][ip];
            }
            X[i][ip] = (M[i][ip]-s)/L[i][i];
        }
    }
}



void TKmultMatrix( double **idcm, double **u,int ndata,int ndata2,int npol,double **uout)
{
#ifdef ACCEL_MULTMATRIX
    if(useT2accel){
        logdbg("Using ACCELERATED multmatrix (M.Keith 2012)");
        accel_multMatrix(idcm[0],u[0],ndata,ndata2,npol,uout[0]);
    } else {
#endif
        logdbg("Using SLOW multmatrix");
        int i,j,k;
        for (j=0;j<npol;j++)
        {
            for (i=0;i<ndata;i++)
            {
                uout[i][j]=0.0;
            }
        }
        for (k=0;k<ndata2;k++) {
            for (i=0;i<ndata;i++){
                for (j=0;j<npol;j++){
                    uout[i][j]+=idcm[i][k]*u[k][j];
                }
            }
        }

#ifdef ACCEL_MULTMATRIX
    }
#endif
}


void TKmultMatrix_sq( double **idcm, double **u,int ndata,int npol,double **uout)
{
    TKmultMatrix(idcm,u,ndata,ndata,npol,uout);
}
void TKmultMatrixVec( double **idcm, double *b,int ndata,int ndata2,double *bout)
{
#ifdef ACCEL_MULTMATRIX
    if (useT2accel){
        logdbg("using accelerated multmatrix (m.keith 2012)");
        accel_multMatrixVec(idcm[0],b,ndata,ndata2,bout);
        //accel_multMatrix(idcm[0],b,ndata,ndata2,1,bout);
    }else{
#endif
        int i,j;
        for (i=0;i<ndata;i++)
        {
            bout[i] = 0.0;
        }

        for (j=0;j<ndata2;j++)
        {
            for (i=0;i<ndata;i++)
                bout[i]+=idcm[i][j]*b[j];
        }

#ifdef ACCEL_MULTMATRIX
    }
#endif
}


void TKmultMatrixVec_sq( double **idcm, double *b,int ndata,double *bout)
{
    TKmultMatrixVec(idcm,b,ndata,ndata,bout);
}



longdouble** malloc_2dLL(size_t rows,size_t cols){
    size_t i;
    longdouble* memory;
    longdouble** m;
    logdbg("Allocate %d x %d longdouble array (%.3f kb)",rows,cols, (double)(rows*cols*sizeof(longdouble)/1024.0));

    m=(longdouble**) calloc(rows,sizeof(longdouble*));
    memory=(longdouble*)malloc(sizeof(longdouble)*rows*cols);

    if(memory==NULL || m==NULL){
        logdbg("Could not allocate %d x %d longdouble array (%.3f kb)",rows,cols, (double)(rows*cols*sizeof(double)/1024.0));
        logerr("Cannot allocate enough memory for array");
        exit(1);
    }
    logdbg("Allocated mem=0x%016llx m=0x%016llx",reinterpret_cast<uint64_t>(memory),reinterpret_cast<uint64_t>(m));
    for(i=0;i<rows;i++){
        m[i]=memory+cols*i;
    }

    logdbg("Accessible memory 0x%016llx -> 0x%016llx",reinterpret_cast<uint64_t>(memory),reinterpret_cast<uint64_t>(memory+rows*cols));

    return m;

}

void free_2dLL(longdouble** m){
    free(m[0]);
    free(m);
}


/**
 * Allocate uinv in a "BLAS/LAPACK" compatile way
 */
double **malloc_uinv(size_t n){
    return malloc_blas(n,n);
}

size_t get_blas_rows(double** uinv){
    size_t* dim=(size_t*)(uinv[0]-2);
    return *dim;
}

size_t get_blas_cols(double** uinv){
    size_t* dim=(size_t*)(uinv[0]-1);
    return *dim;
}
/**
 * Allocate uinv in a "BLAS/LAPACK" compatile way
 * store the dimensions of the array in two secret elements before
 * the main memory allocation. Useful for checks.
 * WARNING: assumes that sizeof(int) <= sizeof(double)
 */
double **malloc_blas(size_t rows,size_t cols){
    double *memory;
    size_t *dimN;
    size_t *dimM;
    double** uinv;
    size_t i;

    logdbg("Allocate %d x %d double array (%.3f kb)",rows,cols, (double)(rows*cols*sizeof(double)/1024.0));

    assert(rows*cols!=0);

    if (sizeof(int) > sizeof(double)){
        logerr("Error, somehow you have a system with sizeof(int) > sizeof(double)");
        exit(1);
    }
    memory=(double*)calloc((rows*cols+2),sizeof(double)); // we add 2 to store dimensions
    uinv=(double**)malloc(sizeof(double*)*rows);
    if(memory==NULL || uinv==NULL){
        logdbg("Could not allocate %d x %d double array (%.3f kb)",rows,cols, (double)(rows*cols*sizeof(double)/1024.0));
        logerr("Cannot allocate enough memory for array");
        exit(1);
    }
    logdbg("Allocated mem=0x%016llx uinv=0x%016llx",reinterpret_cast<uint64_t>(memory),reinterpret_cast<uint64_t>(uinv));
    memory+=2; // the first two bytes are for the dimensions.
    dimN=(size_t*)(memory-2);
    dimM=(size_t*)(memory-1);

    *dimN = rows;
    *dimM = cols;

    logdbg("Secret memory rows=0x%016llx cols=0x%016llx",reinterpret_cast<uint64_t>(dimN),reinterpret_cast<uint64_t>(dimM));

    logdbg("Accessible memory 0x%016llx -> 0x%016llx",reinterpret_cast<uint64_t>(memory),reinterpret_cast<uint64_t>((memory+rows*cols)));
    for(i=0;i<rows;i++){
        uinv[i]=memory+cols*i;
    }
    return uinv;
}

void free_blas(double** m){
    if(debugFlag){
        logdbg("free 0x%016llx",reinterpret_cast<uint64_t>(m[0])-2);
        logdbg("m was %d x %d",get_blas_rows(m),get_blas_cols(m));
    }
    fflush(stdout);
    free(m[0]-2);
    logdbg("free 0x%016llx",reinterpret_cast<uint64_t>(m));
    fflush(stdout);
    free(m);
    logdbg("leaving free_blas");
    fflush(stdout);
}


void free_uinv(double** uinv){
    free_blas(uinv);
}



float** malloc_2df(size_t rows,size_t cols){
    size_t i;
    float* memory;
    float** m;
    logdbg("Allocate %d x %d float array (%.3f kb)",rows,cols, (double)(rows*cols*sizeof(float)/1024.0));

    m=(float**) calloc(rows,sizeof(float*));
    //memory=(float*)malloc(sizeof(float)*rows*cols);
    posix_memalign((void**)&memory,16,sizeof(float)*rows*cols);

    for(i=0;i<rows*cols;i++){
        memory[i]=0;
    }

    if(memory==NULL || m==NULL){
        logdbg("Could not allocate %d x %d float array (%.3f kb)",rows,cols, (double)(rows*cols*sizeof(float)/1024.0));
        logerr("Cannot allocate enough memory for array");
        exit(1);
    }
    logdbg("Allocated mem=0x%016llx m=0x%016llx",reinterpret_cast<uint64_t>(memory),reinterpret_cast<uint64_t>(m));
    for(i=0;i<rows;i++){
        m[i]=memory+cols*i;
    }

    logdbg("Accessible memory 0x%016llx -> 0x%016llx",reinterpret_cast<uint64_t>(memory),reinterpret_cast<uint64_t>(memory+rows*cols));

    return m;

}

void free_2df(float** m){
    free(m[0]);
    free(m);
}


