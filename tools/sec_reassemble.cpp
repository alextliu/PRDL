/* sec_reassemble
 * 
 * reassemble short-time sections
 * 
 * outMatrix = sec_reassemble(inMatirx, hopsize, nSection, nChannel);
 * 
 * compile with:
 *  mex('-v','-R2018a','sec_reassemble.cpp');
*/

#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    mwSize window_len, hopsize, overlap_len, nSection, nChannel, elChannel, nOutMatrix;
    mxComplexDouble *inMatrix, *outMatrix;
    mwIndex i,j,k, iChannel_in, iChannel_out, jSection, jHop;

    // Assign input arguments to C
    inMatrix = mxGetComplexDoubles(prhs[0]);
    window_len = mxGetM(prhs[0]);       // window length
    hopsize = mxGetScalar(prhs[1]);     // hop size
    nSection = mxGetScalar(prhs[2]);    // no. of short-time sections for each channel
    nChannel = mxGetScalar(prhs[3]);    // no. of channels

    overlap_len = window_len - hopsize; // overlap length
    nOutMatrix = nSection * hopsize + overlap_len;       // length of first dimension of outMatrix
    elChannel = nSection * window_len;    // no. of elements of each channel in inMatrix

    // Assign output arguments back to MATLAB
    plhs[0] = mxCreateNumericMatrix(nOutMatrix, nChannel, mxDOUBLE_CLASS, mxCOMPLEX);
    outMatrix = mxGetComplexDoubles(plhs[0]);

    for (i = 0; i < nChannel; i++) {
        iChannel_in = i * elChannel;    // initial index of channel i in inMatrix
        iChannel_out = i * nOutMatrix;              // initial index of channel i in outMatrix
        for (j = 0; j < nSection; j++) {
            jSection = j * window_len + iChannel_in;    // initial index of section j of channel i in inMatrix
            jHop = j * hopsize + iChannel_out;          // initial index of section j of channel i in outMatrix
            for (k = 0; k < window_len; k++) {
                outMatrix[k + jHop].real += inMatrix[k + jSection].real;
                outMatrix[k + jHop].imag += inMatrix[k + jSection].imag;
            }
        }
    }



}