#ifndef DispatchAndRun_h
#define DispatchAndRun_h
 
#include "JMM_CPPLibs/Macros/String.h"
#include "JMM_CPPLibs/Macros/PPCat.h"

#include "Specializations/Curvature2.h"
#include "Experimental/PrescribedCurvature2.h"
#include "Experimental/ConvexityCurvature2.h"

// ------- Custom invocation, with multiple models.  ---------
#define HFMSpecializationMacro(modelName) \
{ \
using StencilDataType = Stencil ## modelName ;\
using HFMI = StencilDataType::HFMI; \
if(model== #modelName){ \
    io.currentSetter=IO::SetterTag::Compute;\
    StencilDataType stencil; \
    HFMI(io, stencil).Run();\
    io.currentSetter=IO::SetterTag::User; return;} \
}

void Run(IO & io){
    const std::string rawModel = io.GetString("model");
    if(rawModel == "ElasticaConvex2"){
        const std::string model="ElasticaConvex2<5>";
        HFMSpecializationMacro(ElasticaConvex2<5>);
    }
    else if(rawModel == "ReedsSheppForwardConvex2"){
        const std::string model="ReedsSheppForwardConvex2";
        HFMSpecializationMacro(ReedsSheppForwardConvex2);
    }
    else if(rawModel == "DubinsConvex2"){
        const std::string model="DubinsConvex2";
        HFMSpecializationMacro(DubinsConvex2);
    }
    else{
        ExceptionMacro("Unrecognized model : " << rawModel);
    }
    
}


#endif 
