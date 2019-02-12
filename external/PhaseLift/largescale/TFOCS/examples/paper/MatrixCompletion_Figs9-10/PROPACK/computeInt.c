/* 
 * Stephen Becker, July 2010
 *
 * Not tested!  If this is slowing you down,
 * try setting opt.eta = 0, which forces reorthogonalization
 * against all previous vectors.
 *
 * */

#include "mex.h"

//#include "math.h"

#define min(X,Y) (((X) < (Y)) ? (X) : (Y))
#define max(X,Y) (((X) > (Y)) ? (X) : (Y))


/* to use "mwSize", need matrix.h, but doesn't work for R2006a */
/* #include "matrix.h" */
/* So, use the following definitions instead: */
#ifndef mwSize
    #define mwSize size_t
#endif
#ifndef mwIndex
    #define mwIndex size_t  /* should make it compatible w/ 64-bit systems */
#endif

void printUsage() {
    //mexPrintf("usage:\tcomputeInt(IO,mu,extra,eta,j)\n");
    mexPrintf("usage:\tcomputeInt(IO,mu,eta)\n");
}

void mexFunction(
         int nlhs,       mxArray *plhs[],
         int nrhs, const mxArray *prhs[]
         )
{
    /* Declare variable */
    mwSize lengthIO, a,b;
    mwIndex i, j, r, s, k;
    double *IO, *mu, *INT;
    //int extra;
    double eta;
    /* Check for proper number of input and output arguments */    
    if ( nrhs != 3 ) {
        printUsage();
        //mexErrMsgTxt("Needs 5 input arguments");
        mexErrMsgTxt("Needs 3 input arguments");
    } 
    a = mxGetN( prhs[0] );
    b = mxGetM( prhs[0] );
    lengthIO = (a>b) ? a : b;

    IO = mxGetPr( prhs[0] );
    mu = mxGetPr( prhs[1] );
    //extra = (int) *mxGetPr( prhs[2] );
    eta = (double) *mxGetPr( prhs[2] );
    //j = (mwIndex) *mxGetPr( prhs[4] );

    a = mxGetN( prhs[1] );
    b = mxGetM( prhs[1] );
    j = (a>b) ? a : b;

    plhs[0] = mxCreateDoubleMatrix(j,1,mxREAL);
    INT = mxGetPr( plhs[0] );

    //mexPrintf("a is %d, b is %d, lengthIO is %d, j is %d\n",
            //a, b, lengthIO, j );

    for ( i = 0 ; i < lengthIO ; i++ ) {
        // IO[i] - 1 to convert to 0-based
        for ( r=IO[i]-1 ; r >= 0 ; r-- ) {  // -1 to get 0-based
            if ( ( mu[r]  < eta ) | (INT[r] == 1) ) 
                break;
            else
                INT[r] = 1;
        }

        // IO[i] - 1 to convert to 0-based
        for ( s = IO[i]-1 ; s < j ; s++ ) { // -1 to get 0-based
            if ( ( mu[s] < eta ) | (INT[s] == 1) ) 
                break;
            else
                INT[s] = 1;
        }
    }
     
}
