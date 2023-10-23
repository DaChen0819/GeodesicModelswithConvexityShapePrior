#include <fstream>
#include <iostream>
#include <stdio.h>
#include "mex.h"
#include "JMM_CPPLibs/Output/MexIO.h"
#include "JMM_CPPLibs/Macros/ExportArrow.h"
#include "JMM_CPPLibs/Output/FileIO.h"
#include "JMM_CPPLibs/LinearAlgebra/ArrayType.h"
#include "JMM_CPPLibs/LinearAlgebra/VectorType.h"
#include "pathTransformation.h"
#include "ConvertRealPathsToDigitalPaths.h"

typedef IO_<MexIO> IO;
typedef typename IO::Msg Msg;
typedef typename IO::WarnMsg WarnMsg;


void mexFunction(int nlhs, mxArray *plhs[],
                                     int nrhs, const mxArray *prhs[] ){
    if(nrhs!=1 || nlhs!=1){
        IO::WarnMsg() << "Exactly one input and one output are expected (in structures).\n";
        return;
    }
    try{
        IO io(prhs[0],plhs);
        io.arrayOrdering = IO::ArrayOrdering::YXZ_ColumnMajor;
        Run<IO>(io);
    }
    catch(const std::exception & e){
        IO::WarnMsg() << "Convert Real Paths to Digital Paths exception.\n " << e.what() << "\n";
    }
}

