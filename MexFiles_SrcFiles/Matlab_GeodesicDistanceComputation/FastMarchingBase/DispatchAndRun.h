#ifndef DispatchAndRun_h
#define DispatchAndRun_h

#include "JMM_CPPLibs/Macros/String.h"
#include "JMM_CPPLibs/Macros/PPCat.h"

#include "Specializations/Isotropic.h"
#include "Specializations/Riemannian.h"
#include "Experimental/AsymmetricQuadratic.h"
#include "Specializations/Curvature2.h"
#include "Specializations/Curvature3.h"
#include "Specializations/QuadLinLag2.h"
#include "Experimental/AsymRander.h"

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
    if(rawModel=="Isotropic2"){
        const std::string model="Isotropic<2>";
        HFMSpecializationMacro(Isotropic<2>);
    }
    else if(rawModel=="Isotropic3"){
        const std::string model="Isotropic<3>";
        HFMSpecializationMacro(Isotropic<3>);
    }
    else if(rawModel=="Riemann2"){
        const std::string model="Riemann<2>";
        HFMSpecializationMacro(Riemann<2>);
    }
    else if(rawModel=="Riemann3"){
        const std::string model="Riemann<3>";
        HFMSpecializationMacro(Riemann<3>);
    }
    else{
        ExceptionMacro("Unrecognized model : " << rawModel);
    }
}


#endif 
