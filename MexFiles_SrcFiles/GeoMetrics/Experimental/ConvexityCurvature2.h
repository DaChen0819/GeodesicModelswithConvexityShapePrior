// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef ConvexityCurvature2_h
#define ConvexityCurvature2_h

#include "Specializations/CommonTraits.h"
#include "Specializations/Curvature2.h"


// --------------- 2D  Convexity ReedsSheppForward model -------------
struct TraitsReedsSheppForwardConvex2 : TraitsR2S1 {
    typedef EulerianStencil<DifferenceType,0,4> StencilType;
};

struct StencilReedsSheppForwardConvex2 final
: HamiltonFastMarching<TraitsReedsSheppForwardConvex2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppForwardConvex2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType xi=1;
    std::string parameterizationOrder;
    
    typedef Traits::BasisReduction<2> ReductionType;
    Voronoi1Vec<ReductionType> reduc;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = index[2]*param.dependScale;
        const typename ReductionType::VectorType v{cos(theta),sin(theta)};
        
        std::array<DifferenceType,3> physicalForward;
        reduc(&physicalForward[0],v/param.gridScale);
        
        for(int k=0;k<3;k++)
            stencil.forward[0][k]=physicalForward[k];
        
        // adding convexity constraint to the augular dimension.
        if(parameterizationOrder == "Clockwise"){
            stencil.forward[0][3].offset=OffsetType{0,0,-1};
            stencil.forward[0][3].baseWeight = 1./square(xi*param.dependScale);
        }
        else if(parameterizationOrder == "AntiClockwise"){
            stencil.forward[0][3].offset=OffsetType{0,0,1};
            stencil.forward[0][3].baseWeight=1./square(xi*param.dependScale);
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc.eps=that->io.template Get<ScalarType>("eps",reduc.eps);
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/dims.back());
        
        if(that->io.HasField("parameterizationOrder")){
            parameterizationOrder=that->io.GetString("parameterizationOrder");
            if(parameterizationOrder!= "Clockwise" && parameterizationOrder!="AntiClockwise"){
                std::cout<<"Please check the input for parameterizationOrder: Clockwise or AntiClockwise."<<std::endl;
                std::cout<<"parameterizationOrder is set to Clockwise."<<std::endl;
            }
        }
        else this->parameterizationOrder="Clockwise";
    }
};


// --------------- 2D  Convexity Dubins model ---------------

struct TraitsDubinsConvex2 : TraitsR2S1 {
//    typedef EulerianStencil<DifferenceType,0,6,3> StencilType;// with relaxization.
    typedef EulerianStencil<DifferenceType,0,6,2> StencilType;
};

struct StencilDubinsConvex2 final
: HamiltonFastMarching<TraitsDubinsConvex2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDubinsConvex2> HFM;
    typedef HFM::StencilDataType Superclass;
    Redeclare1Type(Superclass,OffsetType)
    
    HFM::ParamDefault param;
    ScalarType xi=1.0;
    std::string parameterizationOrder;
    
    typedef Traits::BasisReduction<2> Reduction2DType;
    Voronoi1Vec<Reduction2DType> reduc2;
    
    typedef Traits::BasisReduction<3> Reduction3DType;
    Voronoi1Vec<Reduction3DType> reduc3;
    
    virtual void SetStencil(const IndexType& index, StencilType& stencil) override {
        const ScalarType theta = index[2]*param.dependScale;
        
        std::array<DifferenceType,3> physicalForward;
        const typename Reduction2DType::VectorType v2{cos(theta),sin(theta)};
        reduc2(&stencil.forward[1][0],v2/param.gridScale);
        
        // adding convexity to the angular dimension.
        if(parameterizationOrder == "Clockwise"){
            const VectorType v3{cos(theta)/param.gridScale,sin(theta)/param.gridScale,-1./(xi*param.dependScale)};
            reduc3(&stencil.forward[0][0],v3);
        }
        else if(parameterizationOrder == "AntiClockwise"){
            const VectorType v3{cos(theta)/param.gridScale,sin(theta)/param.gridScale,1./(xi*param.dependScale)};
            reduc3(&stencil.forward[0][0],v3);
        }
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc2.eps=that->io.template Get<ScalarType>("eps",reduc2.eps);
        reduc3.eps=that->io.template Get<ScalarType>("eps",reduc3.eps);
        
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/dims.back());
        
        if(that->io.HasField("parameterizationOrder")){
            parameterizationOrder=that->io.GetString("parameterizationOrder");
            if(parameterizationOrder!= "Clockwise" && parameterizationOrder!="AntiClockwise"){
                std::cout<<"Please check the input for parameterizationOrder: Clockwise or AntiClockwise."<<std::endl;
                std::cout<<"parameterizationOrder is set to Clockwise."<<std::endl;
            }
        }
        else this->parameterizationOrder="Clockwise";
    }
};
// ------------- 2D  Convexity Elastica model ----------------

// Number of Fejer numerical integration points.
template<int nFejer>
struct TraitsElasticaConvex2 : TraitsR2S1 {
    static const DiscreteType nForward=6*nFejer - 3*((nFejer%2)==1);
    typedef EulerianStencil<DifferenceType,0,nForward> StencilType;
};

template<int nFejer>
struct StencilElasticaConvex2 final
: HamiltonFastMarching<TraitsElasticaConvex2<nFejer> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsElasticaConvex2<nFejer> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare8Types(Superclass,Traits,ScalarType,IndexType,VectorType,StencilType,ParamInterface,HFMI,OffsetType)
    Redeclare1Constant(Traits,mathPi)
    
    typename HFM::ParamDefault param;
    ScalarType eps=0.1, xi=1.0;
//    ScalarType relaxation=0.0;
    static const std::array<ScalarType,nFejer> fejerWeights;
    
    //In the convexity-shaped elastica model, we apply two ways to add convexity to the angular dimension.
    bool adaptiveFejer=false;
    
    template<size_t n> using BasisReduction = typename Traits::template BasisReduction<n>;
    Voronoi1Vec<BasisReduction<2> > reduc2;
    Voronoi1Vec<BasisReduction<3> > reduc3;
    std::string parameterizationOrder;
    
    virtual void SetStencil(const IndexType & index, StencilType & stencil) override {
        const ScalarType theta = (2.*mathPi*index[2])/this->dims[2];
        auto & forward = stencil.forward[0];
        for(int l=0; l<nFejer/2; ++l){
            const ScalarType phi = mathPi*(l+0.5)/nFejer;
            const VectorType v{
                sin(phi)*cos(theta)/this->param.gridScale,
                sin(phi)*sin(theta)/this->param.gridScale,
                cos(phi)/(xi*this->param.dependScale)
            };

            const int s0=2*l,s1=2*l+1;
            reduc3(&forward[6*s0], v);
            for(int i=0; i<6; ++i) {
                forward[6*s0+i].baseWeight*=fejerWeights[l];
                forward[6*s1+i] = forward[6*s0+i];
                forward[6*s1+i].offset[2]*=-1;
            }
        }
        if(nFejer%2==1){
            typedef typename Traits::template BasisReduction<2> ReductionType;
            typename ReductionType::VectorType
            v{cos(theta)/this->param.gridScale,sin(theta)/this->param.gridScale};
            const int s = nFejer-1;
            reduc2(&forward[6*s], v);
            for(int i=0; i<3; ++i){
                forward[6*s+i].baseWeight*=fejerWeights[nFejer/2];}
        } // if odd
        
        // adding convexity to the angular dimension.
        // set the weights as 0 if the last coordinate of the corresponding offsets do not obey the criteria.
        if(parameterizationOrder == "Clockwise"){
            for(int i=0;i<6*nFejer - 3*((nFejer%2)==1);i++){
                OffsetType offset=stencil.forward[0][i].offset;
                if(offset[2]>0.0)
                    stencil.forward[0][i].baseWeight=0;
            }
        }
        else if(parameterizationOrder == "AntiClockwise"){
            for(int i=0;i<6*nFejer - 3*((nFejer%2)==1);i++){
                OffsetType offset=stencil.forward[0][i].offset;
                if(offset[2]<0.0)
                    stencil.forward[0][i].baseWeight=0.0;
            }
        }
    }
    
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
        reduc3.eps=that->io.template Get<ScalarType>("eps",reduc3.eps);
        reduc2.eps=reduc3.eps;
        xi=that->io.template Get<ScalarType>("xi",xi);
        param.Setup(that,2*mathPi/this->dims.back());
        
        if(that->io.HasField("parameterizationOrder")){
            parameterizationOrder=that->io.GetString("parameterizationOrder");
            if(parameterizationOrder!= "Clockwise" && parameterizationOrder!="AntiClockwise"){
                std::cout<<"Please check the input for parameterizationOrder: Clockwise or AntiClockwise."<<std::endl;
                std::cout<<"parameterizationOrder is set to Clockwise."<<std::endl;
            }
        }
        else this->parameterizationOrder="Clockwise";
        
        if(that->io.HasField("adaptiveFejer")){
            adaptiveFejer=that->io.template Get<ScalarType>("adaptiveFejer")>0.5;
        }
        else{
            adaptiveFejer=false;
        }
//        relaxation=that->io.template Get<ScalarType>("relaxation",0.0);
    }
};

template<> const std::array<double, 1> StencilElasticaConvex2<1>::fejerWeights = {{2.}};
template<> const std::array<double, 2> StencilElasticaConvex2<2>::fejerWeights = {{1.,1.}};
template<> const std::array<double, 3> StencilElasticaConvex2<3>::fejerWeights =
{{0.444444, 1.11111, 0.444444}};
template<> const std::array<double, 4> StencilElasticaConvex2<4>::fejerWeights =
{{0.264298, 0.735702, 0.735702, 0.264298}};
template<> const std::array<double, 5> StencilElasticaConvex2<5>::fejerWeights =
{{0.167781, 0.525552, 0.613333, 0.525552, 0.167781}};
template<> const std::array<double, 6> StencilElasticaConvex2<6>::fejerWeights =
{{0.118661, 0.377778, 0.503561, 0.503561, 0.377778, 0.118661}};
template<> const std::array<double, 7> StencilElasticaConvex2<7>::fejerWeights =
{{0.0867162, 0.287831, 0.398242, 0.454422, 0.398242, 0.287831, 0.0867162}};
template<> const std::array<double, 8> StencilElasticaConvex2<8>::fejerWeights =
{{0.0669829, 0.222988, 0.324153, 0.385877, 0.385877, 0.324153, 0.222988, 0.0669829}};
template<> const std::array<double, 9> StencilElasticaConvex2<9>::fejerWeights =
{{0.0527366, 0.179189, 0.264037, 0.330845, 0.346384, 0.330845, 0.264037, 0.179189, 0.0527366}};

#endif /* Curvature2_h */
