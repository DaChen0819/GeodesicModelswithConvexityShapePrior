/* function to compute the Bhattacharyya coefficients between two histograms*/
/* mex  -compatibleArrayDims BhattacharyyaCoefficient.cpp */

#include "mex.h"
#include <time.h>
#include <cmath>
#include <iostream>

#define Sub2Ind(ix,iy) (iy)*iNx+ (ix)
#define rint(A) floor((A)+(((A) < 0)? -0.5 : 0.5))


/* This function computes a Gaussian function */
void ComputeGaussianKernel(double* pGaussianKernel, double fStdParzenWindow, int iSizeGaussKernel){
    
    const int iGaussianCenter = (int) ((iSizeGaussKernel-1)/2);
    double f1 = 2.0* (fStdParzenWindow*fStdParzenWindow);
    double f2 = fStdParzenWindow*std::sqrt(2.0*3.1415926);
    
    for (int i=0; i< iSizeGaussKernel; i++)
        pGaussianKernel[i] = (double)(exp(-(i-iGaussianCenter)*(i-iGaussianCenter)/f1)/f2);
    
}
/****************************************/


/****************************************/
/* This function estimates the intensity distribution from an histogram using a Gaussian kernel*/
void ComputeParzenWindPDF(double*        pProb,
                          const double*  pHisto,
                          double         fArea,
                          int            iNumBins,
                          const double*  pGaussianKernel,
                          int            iSizeGaussKernel){

    double   fSumGauss;
    
    const int iNg = iSizeGaussKernel;
    const int iNg2 = (int) ((iSizeGaussKernel-1)/2);
    const int iGaussianCenter = (int) ((iSizeGaussKernel-1)/2);
    
    /* Center of p */
    for (int i=iNg2; i<iNumBins - iNg2; i++){
        for (int iGaussianPos=0; iGaussianPos< iNg; iGaussianPos++){
            pProb[i] += pHisto[i+(iGaussianPos-iGaussianCenter)]* pGaussianKernel[iGaussianPos];
        }
        pProb[i] /= fArea;
    }
    
    
    /* Borders of p */
    for (int i=0; i< iNg2; i++){
        fSumGauss = 0.0;
        for (int iGaussianPos=0;iGaussianPos<iNg; iGaussianPos++){
            int iShift = i+(iGaussianPos-iGaussianCenter);
            if (iShift>= 0  &&  iShift< iNumBins)
                pProb[i] += pHisto[iShift]* pGaussianKernel[iGaussianPos];
            fSumGauss += pGaussianKernel[iGaussianPos];
        }
        pProb[i] /= fArea;
        pProb[i] /= fSumGauss;
    }
    
    for (int i=iNumBins-iNg2; i<iNumBins; i++){
        fSumGauss = 0.0;
        for (int iGaussianPos=0; iGaussianPos< iNg; iGaussianPos++){
            int iShift = i+(iGaussianPos-iGaussianCenter);
            if (iShift>= 0  &&  iShift< iNumBins)
                pProb[i] += pHisto[iShift]*pGaussianKernel[iGaussianPos];
            fSumGauss += pGaussianKernel[iGaussianPos];
        }
        pProb[i] /= fArea;
        pProb[i] /= fSumGauss;
    }
}
/****************************************/

/****************************************/
void ComputeParzenWindFunction_ConvTerm(const double  *pProbIn,
                                        const double  *pProbOut,
                                        double  *pConvTerm,
                                        double  fAreaIn,
                                        double  fAreaOut,
                                        int     iNumBins,
                                        const double  *pGaussianKernel,
                                        int     iSizeGaussKernel){
    
    const int iNg  = iSizeGaussKernel;  //Gaussian window size.
    const int iNg2 = (int) ((iSizeGaussKernel-1)/2); // half of the gaussian window size.
    const int iGaussianCenter = (int) ((iSizeGaussKernel-1)/2);
    
    /* Initialization */
    for (int i=0; i< iNumBins; i++) {
        pConvTerm[i] = 0.0;
    }
    
    /* Center of ConvTerm */
    for(int i=iNg2; i<iNumBins-iNg2; i++){
        for(int iGaussianPos=0; iGaussianPos< iNg; iGaussianPos++){
            int iShift= i+(iGaussianPos-iGaussianCenter);
             double termIn=std::sqrt(pProbOut[iShift]/pProbIn[iShift])/fAreaIn;
             double termOut=std::sqrt(pProbIn[iShift]/pProbOut[iShift])/fAreaOut;
             pConvTerm[i]+=(termIn-termOut)*pGaussianKernel[iGaussianPos];
        }
    }
    
    /* Boders of ConvTerm*/
    for(int i=0; i< iNg2; i++){
        double fSumGauss = 0.0;
        for(int iGaussianPos=0; iGaussianPos< iNg; iGaussianPos++){
            int iShift = i+(iGaussianPos-iGaussianCenter);
            if(iShift>= 0  &&  iShift<iNumBins){
                double termIn=std::sqrt(pProbOut[iShift]/pProbIn[iShift])/fAreaIn;
                double termOut=std::sqrt(pProbIn[iShift]/pProbOut[iShift])/fAreaOut;
                pConvTerm[i] += (termIn-termOut)*pGaussianKernel[iGaussianPos];
                fSumGauss += pGaussianKernel[iGaussianPos];
            }
        }
        pConvTerm[i] /= fSumGauss;
    }
    
    for(int i=iNumBins-iNg2; i<iNumBins; i++){
        double fSumGauss = 0.0;
        for(int iGaussianPos=0; iGaussianPos< iNg; iGaussianPos++){
            int iShift = i+(iGaussianPos-iGaussianCenter);
            if(iShift>= 0  &&  iShift< iNumBins){
                double termIn=std::sqrt(pProbOut[iShift]/pProbIn[iShift])/fAreaIn;
                double termOut=std::sqrt(pProbIn[iShift]/pProbOut[iShift])/fAreaOut;
                pConvTerm[i] += (termIn-termOut)*pGaussianKernel[iGaussianPos];
                fSumGauss += pGaussianKernel[iGaussianPos];
            }
        }
        pConvTerm[i] /= fSumGauss;
    }
}
/****************************************/

double ComputeBhattacharyyaCoefficient(const double* pProbIn, const double *pProbOut,int iNumBins){
    
    double bhaCoefficient=0.0;
    
    for (int i=0; i< iNumBins; i++)
        bhaCoefficient+=std::sqrt(pProbIn[i]*pProbOut[i]);
    
    return bhaCoefficient;
}


/****************************************/
void ComputeBhattacharyyaGradientFlow(const double* pInputImage,
                                      const double* pInputShape,
                                      int     iNx,
                                      int     iNy,
                                      double* pGradientFlow,
                                      double* pProbIn,
                                      double* pProbOut,
                                      double* pHistoIn,
                                      double* pHistoOut,
                                      int     iNumBins,
                                      const double* pGaussianKernel,
                                      int     iSizeGaussKernel,
                                      double  *pConvTerm){
    
    double fNormalizationIn = 0.0;
    double fNormalizationOut = 0.0;
    
    // initialize the histograms
    for (int i=0; i< iNumBins; i++) {
        pHistoIn[i] = 0.0;
        pHistoOut[i] = 0.0;
    }
    
    /* Compute histograms inside and outside the boundary of {u>0.5} */
    for (int ix=0; ix< iNx; ix++){
        for (int iy=0; iy< iNy; iy++){
            int iIndexIntensity = (int)(rint(pInputImage[Sub2Ind(ix,iy)]));
            if (iIndexIntensity<=0)
                iIndexIntensity=0;
            if (pInputShape[Sub2Ind(ix,iy)] >= 0.5)
                pHistoIn[iIndexIntensity] += 1.0;
            else
                pHistoOut[iIndexIntensity]+= 1.0;
            
            fNormalizationIn  += pInputShape[Sub2Ind(ix,iy)];
            fNormalizationOut += 1.0-pInputShape[Sub2Ind(ix,iy)];
        }
    }
    
    for (int i=0; i< iNumBins; i++){
        if( pHistoIn[i] < 1 )
            pHistoIn[i] = 1;
                
        if( pHistoOut[i] < 1 )
            pHistoOut[i] = 1;
    }
    
    
    /* Estimate the intensity distributions inside and outside from histograms, using the Parzen estimation method */
    ComputeParzenWindPDF(pProbIn,pHistoIn,fNormalizationIn,iNumBins,pGaussianKernel,iSizeGaussKernel);
    ComputeParzenWindPDF(pProbOut,pHistoOut,fNormalizationOut,iNumBins,pGaussianKernel,iSizeGaussKernel);
    ComputeParzenWindFunction_ConvTerm(pProbIn,pProbOut,pConvTerm,fNormalizationIn,fNormalizationOut,iNumBins,pGaussianKernel,iSizeGaussKernel);
    double bhaCoefficient=ComputeBhattacharyyaCoefficient(pProbIn, pProbOut,iNumBins);
    
    
    for (int ix=0; ix<iNx; ix++){
        for (int iy=0; iy< iNy; iy++){
            int iIndexIntensity =(int)rint(pInputImage[Sub2Ind(ix,iy)]);
            
            if(iIndexIntensity<=0)
                iIndexIntensity=0;
            
            pGradientFlow[Sub2Ind(ix,iy)]=0.5*bhaCoefficient*(1.0/fNormalizationOut-1.0/fNormalizationIn)+0.5*pConvTerm[iIndexIntensity];
            
        }
    }
}


void mexFunction(int nlhs, mxArray* plhs[],
				 int nrhs, const mxArray* prhs[] ) {
    
    
    if(nrhs!=3)
        mexErrMsgTxt("3 input arguments are required.");
    
    if(nlhs!=3)
        mexErrMsgTxt("3 output arguments are required.");
    
    
    // given image.
    double* pInputImage = mxGetPr(prhs[0]);
    
    //current shape.
    double* pInputShape = mxGetPr(prhs[1]);
    
    //image size
    const mwSize iNx = mxGetM(prhs[0]);
    const mwSize iNy = mxGetN(prhs[0]);
    
    
    
    //parameters:
    //1. number of parameters
    //2. number of bins for histogram
    //3. sigma of Gaussian kernel
    double* pVecParameters=mxGetPr(prhs[2]);
    const int iNumParameters=(int)pVecParameters[0];
    const double iNumBins=pVecParameters[1];
    
    //std of the Gaussian kernel
    const double fStdParzenWindow = pVecParameters[2];
    
    
    // the size of the truncated gaussian window.
    int iSizeGaussKernel = (int)rint(6.0*fStdParzenWindow);
    iSizeGaussKernel = (iSizeGaussKernel%2 == 1)?iSizeGaussKernel:iSizeGaussKernel+1;
    
    
    //gradient flow
    plhs[0] = mxCreateDoubleMatrix(iNx,iNy,mxREAL);
    double* pGradientFlow = mxGetPr(plhs[0]);
    
    
    //prob  inside the contour
    plhs[1] = mxCreateDoubleMatrix(iNumBins+1, 1, mxREAL);
    double* pProbIn = mxGetPr(plhs[1]);
    
    //prob  outside the contour
    plhs[2] = mxCreateDoubleMatrix(iNumBins+1, 1, mxREAL);
    double* pProbOut = mxGetPr(plhs[2]);
    
    //histogram inside the contour
    double* pHistoIn = (double*)calloc((unsigned)(iNumBins+1),sizeof(double));
    
    //histogram outside the contour
    double* pHistoOut = (double*)calloc((unsigned)(iNumBins+1),sizeof(double));
    
    //histogram outside the contour
    double* pConvTerm = (double*)calloc((unsigned)(iNumBins+1),sizeof(double));
    
    //gaussian kernel.
    double* pGaussianKernel=(double*)calloc((unsigned)(iSizeGaussKernel),sizeof(double));
    if(!pGaussianKernel)
        mexPrintf("Memory allocation failure\n");
    
    // construct the gaussian kernel.
    
    clock_t timing=-clock();
    
    ComputeGaussianKernel(pGaussianKernel,fStdParzenWindow,iSizeGaussKernel);
    ComputeBhattacharyyaGradientFlow(pInputImage,pInputShape,iNx,iNy,pGradientFlow,pProbIn,pProbOut,pHistoIn,pHistoOut,iNumBins,pGaussianKernel,iSizeGaussKernel,pConvTerm);
    
    timing+=clock();
    bool showData=false;
    if(iNumParameters==3)
        showData=pVecParameters[2]>0.5;
    else if(iNumParameters==2)
        showData=false;
    else
        mexErrMsgTxt("incorrect input parameters.");
       
    if(showData){  
        std::cout<<"number of parameters = "<<iNumParameters<<"."<<std::endl;
        std::cout<<"number of bins = "<<iNumBins<<"."<<std::endl;
        std::cout<<"Gaussian sigma = "<<fStdParzenWindow<<"."<<std::endl;
        std::cout<<"Image Size=["<<iNx<<", "<<iNy<<"]."<<std::endl;
        std::cout<<"Computing Time for the computation of Bhattacharyya Flow: "<<timing/double(CLOCKS_PER_SEC)<<"."<<std::endl;
    }

    
    delete [] pConvTerm;
    delete [] pHistoIn;
    delete [] pHistoOut;
    delete [] pGaussianKernel;
}





