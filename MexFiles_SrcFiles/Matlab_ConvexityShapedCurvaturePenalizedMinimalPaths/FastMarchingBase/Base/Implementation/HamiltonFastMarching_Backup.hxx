// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HamiltonFastMarching_hxx
#define HamiltonFastMarching_hxx

// ----- Printing some structures ------ //
template<typename Traits> void
HamiltonFastMarching<Traits>::FlowDataType::PrintSelf(std::ostream & os) const {
    os << "{" << flow << "," << value << "," << width << "}";
}

template<typename Traits> void
HamiltonFastMarching<Traits>::FullIndexType::PrintSelf(std::ostream & os) const {
	os << "{" << index << "," << linear << "}";
}

// --------- Construction -------

template<typename T> HamiltonFastMarching<T>::
HamiltonFastMarching(StencilDataType & _stencilData):
dom(_stencilData.dims), stencilData(_stencilData){
    values.dims = stencilData.dims;
    values.resize(values.dims.Product(),Traits::Infinity());
    
    acceptedFlags.dims=values.dims;
    acceptedFlags.resize(values.size(),false);
    
    activeNeighs.dims = values.dims;
    activeNeighs.resize(values.size());
    stencilData.Initialize(this);
};

template<typename T> auto HamiltonFastMarching<T>::
MaxStencilWidth() const -> DiscreteType {
    return std::accumulate(
    stencilData.reversedOffsets.values.begin(),
	stencilData.reversedOffsets.values.end(),0,
    [](DiscreteType a, const OffsetType & offset)->DiscreteType {
    return std::max(a,std::accumulate(offset.begin(),offset.end(),0,
    [](DiscreteType b, DiscreteType c)->DiscreteType {return b+std::abs(c);}
                                      ));});
}


// ------------ Running Fast Marching -------------
template<typename T>
struct HamiltonFastMarching<T>::QueueElement {
    DiscreteType linearIndex;
    ScalarType value;
    bool operator < (const QueueElement & other) const {return value>other.value;}
};


template<typename T>
void HamiltonFastMarching<T>::Run(){
    RunInit();
    do {} while(!RunOnce());
}

template<typename T>
void HamiltonFastMarching<T>::RunInit(){
    assert(values.CheckDims());
    assert(values.size()>0);
    
    thetaScale = (2.0*Traits::mathPi)/((ScalarType)this->values.dims[2]);
    
    if(acceptedFlags.empty()){
        acceptedFlags.dims = values.dims;
        acceptedFlags.resize(values.size(),false);
    } else {
        assert(acceptedFlags.dims==values.dims);
        assert(acceptedFlags.CheckDims());
    }
    
    if(totalCurvature.empty()){
        totalCurvature.dims = values.dims;
        totalCurvature.resize(values.size(),Traits::Infinity());
    }
    
    for(const auto & [index,value] : seeds){
        const DiscreteType linearIndex = values.Convert(index);
        values[linearIndex] = value;
//        acceptedFlags[linearIndex] = true;
        totalCurvature[linearIndex] = 0.0;
        queue.push({linearIndex,value});
    }
}

template<typename T>
bool HamiltonFastMarching<T>::ConvexSimpleGeodesic(IndexType acceptedIndex, VectorType& flowSum){
    flowSum=VectorType::Constant(0);
    ScalarType weightSum=0;
    DiscreteFlowType flow;
    ScalarType ProcessorCurvature=0.0;
    Recompute(acceptedIndex, flow);
    
    for(const DiscreteFlowElement & fl : flow){
        if(fl.weight==0)
            continue;
        weightSum+=fl.weight;
        flowSum+=fl.weight * VectorType::CastCoordinates(fl.offset);
        IndexType neigh = acceptedIndex;
        for(int i=0; i<Dimension; ++i)
            neigh[i]+=fl.offset[i];
        const auto transform = this->dom.Periodize(neigh,acceptedIndex);
        assert(transform.IsValid());
        (void)transform;
        ProcessorCurvature+=fl.weight*totalCurvature(neigh);
    }
    if(weightSum>0)
        flowSum/=weightSum;
    else
        return false; // seeds.
    
    ScalarType curvature=0.0;
    curvature=(flowSum[Dimension-1]*thetaScale);
    totalCurvature(acceptedIndex)=ProcessorCurvature/weightSum+curvature;
    if(std::abs(totalCurvature(acceptedIndex))>2.01*Traits::mathPi)
        return true;
    
    return false;
}


template<typename T>
bool HamiltonFastMarching<T>::RunOnce(){
    QueueElement top = queue.top();
    queue.pop();
    
    if(acceptedFlags[top.linearIndex])
        return queue.empty();
	
   
    
    const FullIndexType accepted = {values.Convert(top.linearIndex),top.linearIndex};
    acceptedFlags[accepted.linear]=true;
    stencilData.EraseCache(accepted.linear);
    
    // ------ opposite gradient flow --------
    VectorType realFlow;
    if(ConvexSimpleGeodesic(accepted.index,realFlow) && useSimplicity)
        return queue.empty();
        
        
//    if(parameterizationOrder == "Clockwise" && realFlow[Dimension-1]>0.0)
//        return true;
//        
//    if(parameterizationOrder == "AntiClockwise" && realFlow[Dimension-1]<0.0)
//        return true;
    
    int dec = PostProcess(accepted.index);
    if(dec & Decision::kRecompute){
        DiscreteFlowType flow;
        const RecomputeType rec = Recompute(accepted.index, flow);
		top.value = rec.value;
        dec|=PostProcessWithRecompute(accepted.index, rec, flow);
    }
    if(dec & Decision::kTerminate) return true;
    if(dec & Decision::kContinue)  return queue.empty();
    
    
    const auto offsets = stencilData.ReversedOffsets(accepted);
    for(OffsetCRef offset : offsets){
        if(offsetRemoval && parameterizationOrder == "Clockwise" && offset[Dimension-1]>1e-5)
            continue;
        if(offsetRemoval && parameterizationOrder == "AntiClockwise" && offset[Dimension-1]<-1e-5)
            continue;
        ConditionalUpdate(accepted.index, offset, top.value);
    }

    return queue.empty();
}

template<typename T> int HamiltonFastMarching<T>::
PostProcess(IndexCRef acceptedIndex) {
    int result =
	(order>1
	 || factoring.NeedsRecompute(acceptedIndex)
	 || !extras.postProcessWithRecompute.empty())
	? Decision::kRecompute : Decision::kAccept;
	for(ExtraAlgorithmInterface * p : extras.postProcess) {
		result|=p->PostProcess(acceptedIndex);}
    return result;
}

template<typename T> int HamiltonFastMarching<T>::
PostProcessWithRecompute(IndexCRef acceptedIndex, const RecomputeType & rec,
						 const DiscreteFlowType & flow){
    values(acceptedIndex)=rec.value;
    int result = Decision::kAccept;
    for(ExtraAlgorithmInterface * p : extras.postProcessWithRecompute)
        result|=p->PostProcessWithRecompute(acceptedIndex, rec, flow);
    return result;
}

template<typename T> void
HamiltonFastMarching<T>::ConditionalUpdate(IndexCRef acceptedIndex,
                                           OffsetType offset,
                                           ScalarType acceptedValue){
    FullIndexType updated;
    const auto transform = VisibleOffset(acceptedIndex,-offset, updated.index);
	
    if(!transform.IsValid()) return;
    transform.PullVector(offset);    
    updated.linear = values.Convert(updated.index);
    if(acceptedFlags[updated.linear]) return;
    // Next line forbids updating of seeds or given data.
    // Can be specialized to allow e.g. for sequential computation of Voronoi diagrams.
    if(activeNeighs[updated.linear].none() &&
	   values[updated.linear]!=Traits::Infinity()) return;
    if(values[updated.linear]<=acceptedValue) return;
    
    PhysicalIndexType updatedPhysicalIndex;
    for(int k=0;k<Dimension-1;k++)
        updatedPhysicalIndex[k]=updated.index[k];
    
    
    if(!physicalActiveRegion.empty() && physicalActiveRegion(updatedPhysicalIndex)<0.5)
        return;
    
    if(ScribblesCrossingDetection(acceptedIndex,updated.index))
         return;
    

    if(AdaptiveCrossingDetection(acceptedIndex,updated.index))
        return;
	
    Update(updated, offset, acceptedValue); // also pushes in queue
}


template<typename T> void HamiltonFastMarching<T>::
Update(FullIndexCRef updated, OffsetCRef offset, ScalarType acceptedValue){
    auto & active = activeNeighs[updated.linear];

/*    // Alternatively, only insert in queue if value is strictly decreased.
    ScalarType & val = values[updatedLinearIndex];
    if(result.first>=val) return;
 */
    const ScalarType updatedValue =
    stencilData.HopfLaxUpdate(updated,offset,acceptedValue,active);
    values[updated.linear] = updatedValue;
    queue.push({updated.linear,updatedValue});    
}

// ----------------- Boundary conditions -------------------

template<typename T> auto HamiltonFastMarching<T>::
VisibleOffset(const IndexType & acceptedIndex, const OffsetType & offset,
			  IndexType & updatedIndex) const -> DomainTransformType {
    updatedIndex=acceptedIndex+IndexDiff::CastCoordinates(offset);
    DomainTransformType result = dom.Periodize(updatedIndex,acceptedIndex);
    if(!result.IsValid()) return result;
    for(ExtraAlgorithmInterface * p : extras.visible){
        if(!p->Visible(acceptedIndex,offset,updatedIndex)) result.Invalidate();
	}
    return result;
}

// ---------- Recompute ----------

template<typename T> auto HamiltonFastMarching<T>::
Recompute(IndexCRef updatedIndex, DiscreteFlowType & discreteFlow)
const -> RecomputeType {
    assert(discreteFlow.empty());
    for(ExtraAlgorithmInterface * p : extras.beforeRecompute) {
		p->BeforeRecompute(updatedIndex);}
    const DiscreteType updatedLinearIndex = values.Convert(updatedIndex);
    const ActiveNeighFlagType active = activeNeighs[updatedLinearIndex];
    if(active.none()) return {values[updatedLinearIndex],0.};
	
	factoring.SetIndex(updatedIndex);
	
	// Used in criteria for ditching the high order scheme
	const ScalarType oldValue = values[updatedLinearIndex];
    
    auto GetValueCorr = [this,&updatedIndex,&oldValue]
	(OffsetType offset, int & ord) -> ScalarType {
        //order code : 0 -> invalid, else requested/used order
		
        IndexType acceptedIndex = updatedIndex+IndexDiff::CastCoordinates(offset);
        const auto transform = dom.Periodize(acceptedIndex,updatedIndex);
        if(!transform.IsValid()) {ord=0; return -Traits::Infinity();}
        
        // ----------------- used for winding geodesic curves only -----------------//
        if(useExtraCrossingCheck){
            if(ScribblesCrossingDetection(updatedIndex,acceptedIndex)){
                ord=0;
                return -Traits::Infinity();
            }

            if(AdaptiveCrossingDetection(updatedIndex,acceptedIndex)){
                ord=0;
                return -Traits::Infinity();
            }
        }
        // --------------------------------------------------------------------//
        
        
		const DiscreteType acceptedLinearIndex = values.Convert(acceptedIndex);
		if(!acceptedFlags[acceptedLinearIndex]) {ord=0; return -Traits::Infinity();}
		
        const ScalarType acceptedValue = values[acceptedLinearIndex];

		ord=std::min(order,ord);
		while(ord>=2){ // Single iteration
			OffsetType offset2 = offset;
			transform.PullVector(offset2);
			IndexType acceptedIndex2;
			const auto transform2 = VisibleOffset(acceptedIndex, offset2, acceptedIndex2);
			if(!transform2.IsValid()) break;
			const DiscreteType acceptedLinearIndex2 = values.Convert(acceptedIndex2);
			if(!acceptedFlags[acceptedLinearIndex2]) break;
			const ScalarType acceptedValue2 = values(acceptedIndex2);
			// Ditch if non-causal. Implied by next test, if maxRatioOrder2 <=1
			if(acceptedValue2>acceptedValue) break;
			
			// Estimate only reasonable if the scheme is strictly causal
			const ScalarType offsetNormApprox = oldValue-acceptedValue;
			if(strictlyCausal){ // Ditch if not a sufficiently small correction
				const ScalarType diff2 = oldValue-2*acceptedValue+acceptedValue2;
				if(std::abs(diff2) > maxRatioOrder2*offsetNormApprox) break;
			}
			
			while(ord>=3){ // Single iteration
				OffsetType offset3 = offset2;
				transform2.PullVector(offset3);
				IndexType acceptedIndex3;
				const auto transform3 =
				VisibleOffset(acceptedIndex2, offset3, acceptedIndex3);
				if(!transform3.IsValid()) break;
				const DiscreteType acceptedLinearIndex3 = values.Convert(acceptedIndex3);
				if(!acceptedFlags[acceptedLinearIndex3]) break;
				const ScalarType acceptedValue3 = values(acceptedIndex3);
				//Ditch if non-causal
				if(acceptedValue3>acceptedValue2) break;
				
				if(strictlyCausal){// Ditch if not a sufficiently small correction
					const ScalarType diff3 = oldValue-3*acceptedValue+3*acceptedValue2-acceptedValue3;
					if(std::abs(diff3) > maxRatioOrder3*offsetNormApprox) break;
				}
				
				ord=3;
				return (6./11.)*
				(3.*acceptedValue-1.5*acceptedValue2+(1./3.)*acceptedValue3
				 +factoring.Correction(offset,3)
				 );
			}
			ord=2;
			return (2./3.)*
			(2.*acceptedValue-0.5*acceptedValue2
			 +factoring.Correction(offset,2)
			 );
		}
        ord=1;
		return acceptedValue
		+factoring.Correction(offset,1);
    };
	
	assert(discreteFlow.empty());
    return stencilData.HopfLaxRecompute(GetValueCorr,updatedIndex,active,discreteFlow);
    
}

// ------- Getting values around a point ------

template<typename Traits> template<bool useFactoring, bool smallCorrection>
void HamiltonFastMarching<Traits>::SetIndex(IndexCRef index) const {

	auto & tmp = getNeighborValue_tmp;
	tmp.index = index;
    
	if(smallCorrection){
		tmp.value = values(index);
	}
    
	if(useFactoring){factoring.SetIndex(index);}
}

template<typename Traits> template<bool useFactoring, bool smallCorrection, int maxOrder>
auto HamiltonFastMarching<Traits>::
GetNeighborValue(OffsetType offset,int& ord)
const -> ScalarType {
	//order code : 0 -> invalid, else requested/used order
	assert(ord<=maxOrder);
	
	const auto & index = getNeighborValue_tmp.index;
	const auto & oldValue = getNeighborValue_tmp.value;
	
	IndexType acceptedIndex = index+IndexDiff::CastCoordinates(offset);
	const auto transform = dom.Periodize(acceptedIndex,index);
	if(!transform.IsValid()) {ord=0; return -Traits::Infinity();}
	const DiscreteType acceptedLinearIndex = values.Convert(acceptedIndex);

	if(!acceptedFlags[acceptedLinearIndex]) {ord=0; return -Traits::Infinity();}
	const ScalarType acceptedValue = values[acceptedLinearIndex];
	
	ord=std::min(order,ord);
	while(maxOrder>=2 && ord>=2){ // Single iteration
		OffsetType offset2 = offset;
		transform.PullVector(offset2);
		IndexType acceptedIndex2;
		const auto transform2 = VisibleOffset(acceptedIndex, offset2, acceptedIndex2);
		if(!transform2.IsValid()) break;
		const DiscreteType acceptedLinearIndex2 = values.Convert(acceptedIndex2);
		if(!acceptedFlags[acceptedLinearIndex2]) break;
		const ScalarType acceptedValue2 = values(acceptedIndex2);
		// Ditch if non-causal. (Implied by next test, smallCorrection, if maxRatioOrder2 <=1)
		if(acceptedValue2>acceptedValue) break;
		
		// The estimate below is only reasonable if the scheme is strictly causal
		const ScalarType offsetNormApprox = oldValue-acceptedValue;
		if(smallCorrection){ // Ditch if not a sufficiently small correction
			const ScalarType diff2 = oldValue-2*acceptedValue+acceptedValue2;
			if(std::abs(diff2) > maxRatioOrder2*offsetNormApprox) break;
		}
		
		while(maxOrder>=3 && ord>=3){ // Single iteration
			OffsetType offset3 = offset2;
			transform2.PullVector(offset3);
			IndexType acceptedIndex3;
			const auto transform3 =
			VisibleOffset(acceptedIndex2, offset3, acceptedIndex3);
			if(!transform3.IsValid()) break;
			const DiscreteType acceptedLinearIndex3 = values.Convert(acceptedIndex3);
			if(!acceptedFlags[acceptedLinearIndex3]) break;
			const ScalarType acceptedValue3 = values(acceptedIndex3);
			//Ditch if non-causal
			if(acceptedValue3>acceptedValue2) break;
			
			if(smallCorrection){// Ditch if not a sufficiently small correction
				const ScalarType diff3 =
				oldValue-3*acceptedValue+3*acceptedValue2-acceptedValue3;
				if(std::abs(diff3) > maxRatioOrder3*offsetNormApprox) break;
			}
            
			ord=3;
			ScalarType result =
			3.*acceptedValue -1.5*acceptedValue2 +(1./3.)*acceptedValue3;
			if(useFactoring) {result+= factoring.Correction(offset,3);}
			return (6./11.)* result;
		}
		ord=2;
		ScalarType result = 2.*acceptedValue -0.5*acceptedValue2;
		if(useFactoring) {result+= factoring.Correction(offset,2);}
		return (2./3.)*result;
	}
	ord=1;
	ScalarType result = acceptedValue;
	if(useFactoring) result+=factoring.Correction(offset,1);
	return result;
}

template<typename Traits> auto HamiltonFastMarching<Traits>::
GeodesicFlow(const IndexType & index) const -> FlowDataType {
    DiscreteFlowType discreteFlow;
    const RecomputeType & rec = Recompute(index,discreteFlow);
    FlowDataType result;
    result.value =  rec.value;
    result.width = rec.width;
    
    result.flow.fill(0.);
    for(const auto & [offset,weight] : discreteFlow){
        assert(weight>=0);
        for(int i=0; i<Dimension; ++i){
            result.flow[i]+=weight*ScalarType(offset[i]);}
    }
    return result;
}


template<typename T> bool HamiltonFastMarching<T>::
AdaptiveCrossingDetection(IndexType acceptedIndex,IndexType updatedIndex) const{
    
    if(taggedAdaptiveCutRegion.empty() || adaptiveCut.empty())
        return false;
    
    PhysicalIndexType acceptedPhysicalIndex, updatedPhysicalIndex;
    for(int k=0;k<Dimension-1;k++){
        acceptedPhysicalIndex[k]=acceptedIndex[k];
        updatedPhysicalIndex[k]=updatedIndex[k];
    }
    
    int acceptedLabel = (int)taggedAdaptiveCutRegion(acceptedPhysicalIndex);
    int updatedLabel  = (int)taggedAdaptiveCutRegion(updatedPhysicalIndex);
    
    if(acceptedLabel==0 || updatedLabel==0)
        return false;
    
    // 2 means anti-clockwise part; 3 means clockwise part; 4 means overlap part.
    if(adaptiveCut.find(acceptedPhysicalIndex)!=adaptiveCut.end() && adaptiveCut.find(updatedPhysicalIndex)!=adaptiveCut.end()
       && acceptedLabel==4 && updatedLabel==4){
        return false;
    }
    
    if(acceptedLabel==updatedLabel && acceptedLabel!=4){
        return false;
    }
    if(updatedLabel==2 && acceptedLabel==3){
        return true;
    }
    else if(updatedLabel==3 && acceptedLabel==2){
        return true;
    }
    else if(updatedLabel==4 && acceptedLabel>0.5){
        if(adaptiveCut.find(acceptedPhysicalIndex)!=adaptiveCut.end())
            return true;
        
        if(adaptiveCut.find(updatedPhysicalIndex)!=adaptiveCut.end())
            return true;
        
        PointType acceptedPoint=dom.PointFromIndex(acceptedIndex);
        PointType updatedPoint=dom.PointFromIndex(updatedIndex);
        VectorType tangent=VectorType::Constant(0.0);
        const ScalarType scale=0.5;
        for(int k=0;k<Dimension-1;k++)
            tangent[k]=(ScalarType)(updatedPoint[k]-acceptedPoint[k]);
        
        tangent=tangent.Normalized();
        
        PhysicalIndexType currentPhysicalIndex;
        for(int k=0; k<Dimension-1;k++)
            currentPhysicalIndex[k]=acceptedIndex[k];
        
        double count=1.0;
        while(currentPhysicalIndex!=updatedPhysicalIndex && count<=15000){
            PointType currentPoint;
            for(int k=0;k<Dimension-1;k++){
                currentPoint[k]=acceptedPoint[k]+tangent[k]*scale*count;
                currentPhysicalIndex[k]=(DiscreteType)std::floor(currentPoint[k]);
            }
            count++;
            if(adaptiveCut.find(currentPhysicalIndex)!=adaptiveCut.end())
                return true;
        }
        
    }
    else if(acceptedLabel==4 && updatedLabel>0.5){
        
        if(adaptiveCut.find(acceptedPhysicalIndex)!=adaptiveCut.end())
            return true;
        
        if(adaptiveCut.find(updatedPhysicalIndex)!=adaptiveCut.end())
            return true;
        
        PointType acceptedPoint=dom.PointFromIndex(acceptedIndex);
        PointType updatedPoint=dom.PointFromIndex(updatedIndex);
        VectorType tangent=VectorType::Constant(0.0);
        const ScalarType scale=0.5;
        for(int k=0;k<Dimension-1;k++)
            tangent[k]=(ScalarType)(updatedPoint[k]-acceptedPoint[k]);
        tangent=tangent.Normalized();
        
        PhysicalIndexType currentPhysicalIndex;
        for(int k=0; k<Dimension-1;k++)
            currentPhysicalIndex[k]=acceptedIndex[k];
        
        double count=1.0;
        while(currentPhysicalIndex!=updatedPhysicalIndex && count<=15000){
            PointType currentPoint;
            for(int k=0;k<Dimension-1;k++){
                currentPoint[k]=acceptedPoint[k]+tangent[k]*scale*count;
                currentPhysicalIndex[k]=(DiscreteType)std::floor(currentPoint[k]);
            }
            count++;
            if(adaptiveCut.find(currentPhysicalIndex)!=adaptiveCut.end())
                return true;
        }
    }
    
    return false;
}

template<typename T> bool HamiltonFastMarching<T>::
ScribblesCrossingDetection(IndexType acceptedIndex,IndexType updatedIndex) const{
    
    
    if(scribbleRegion.empty() || scribbles.empty())
        return false;
    
    
    PhysicalIndexType acceptedPhysicalIndex, updatedPhysicalIndex;
    for(int k=0;k<Dimension-1;k++){
        acceptedPhysicalIndex[k]= acceptedIndex[k];
        updatedPhysicalIndex[k] = updatedIndex[k];
    }
    
     // scribbles to prevent paths to cross over.
    
    int acceptedLabel=(int)scribbleRegion(acceptedPhysicalIndex);
    int updatedLabel=(int)scribbleRegion(updatedPhysicalIndex);
    
    if(acceptedLabel==0 || updatedLabel==0)
        return false;
    
    if(scribbles.find(acceptedPhysicalIndex)!=scribbles.end())
        return true;
    
    if(scribbles.find(updatedPhysicalIndex)!=scribbles.end())
        return true;
    
    if(updatedLabel>0.5 && acceptedLabel>0.5){
        VectorType tangent=VectorType::Constant(0.0);
        const ScalarType scale=0.5;
        
        PointType acceptedPoint=dom.PointFromIndex(acceptedIndex);
        PointType updatedPoint=dom.PointFromIndex(updatedIndex);
        for(int k=0;k<Dimension-1;k++)
            tangent[k]=(ScalarType)(updatedPoint[k]-acceptedPoint[k]);
        tangent=tangent.Normalized();
        
        PhysicalIndexType currentPhysicalIndex;
        for(int k=0; k<Dimension-1;k++)
            currentPhysicalIndex[k]=acceptedIndex[k];
        
        double count=1.0;
        while(currentPhysicalIndex!=updatedPhysicalIndex && count<=15000){
            PointType currentPoint;
            for(int k=0;k<Dimension-1;k++){
                currentPoint[k]=acceptedPoint[k]+tangent[k]*scale*count;
                currentPhysicalIndex[k]=(DiscreteType)std::floor(currentPoint[k]);
            }
            count++;
            if(scribbles.find(currentPhysicalIndex)!=scribbles.end())
                return true;
        }
    }
    return false;
}


#endif /* HamiltonFastMarching_hxx */
