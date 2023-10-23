// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HamiltonFastMarching_hxx
#define HamiltonFastMarching_hxx

// ----- Printing some structures ------

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
    
    tubularNeighbourhood.dims=values.dims;
    tubularNeighbourhood.resize(values.size(),0.0);
    
    normalizedGeoFlows.dims=values.dims;
    normalizedGeoFlows.resize(values.size(),VectorType::Constant(0.0));
    
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
    
    if(acceptedFlags.empty()){
        acceptedFlags.dims = values.dims;
        acceptedFlags.resize(values.size(),false);
    }
    else {
        assert(acceptedFlags.dims==values.dims);
        assert(acceptedFlags.CheckDims());
    }
    
    
    for(const auto & [index,value] : seeds){
        const DiscreteType linearIndex = values.Convert(index);
        values[linearIndex] = value;
        listSeedIndices.insert(index);
        queue.push({linearIndex,value});
    }
    if(Dimension==2){
        OffsetType offset;
        for(int i=-1;i<=1;i++){
            offset[0]=i;
            for(int j=-1;j<=1;j++){
                offset[1]=j;
                if(std::abs(i)+std::abs(j)==0)
                    continue;
                neighbourOffsetList.push_back(offset);
            }
        }
    }
    
    if(Dimension==3){
        OffsetType offset;
        for(int i=-1;i<=1;i++){
            offset[0]=i;
            for(int j=-1;j<=1;j++){
                offset[1]=j;
                for(int k=-1;k<=1;k++){
                    if(std::abs(i)+std::abs(j)+std::abs(k)==0)
                        continue;
                    neighbourOffsetList.push_back(offset);
                }
            }
        }
    }

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
    
    DiscreteFlowType discreteFlow;
    const RecomputeType rec = Recompute(accepted.index,discreteFlow);
    if(listSeedIndices.find(accepted.index)==listSeedIndices.end()){
        VectorType geodesicFlow;
        for(const auto & [offset,weight]:discreteFlow){
            assert(weight>=0);
            for(int i=0; i<Dimension; ++i)
                geodesicFlow[i]+=weight*ScalarType(offset[i]);
        }
        normalizedGeoFlows(accepted.index)=geodesicFlow.Normalized();
    }
    const bool isInsideTubeNeigh=IsInsideTubularNeighbourhood(accepted.index);
    int dec = PostProcess(accepted.index);
    if(dec & Decision::kRecompute){
		top.value = rec.value;
        dec|=PostProcessWithRecompute(accepted.index, rec,discreteFlow);
    }
    
    if(dec & Decision::kTerminate){
        if(isInsideTubeNeigh)
            tubularNeighbourhood[accepted.linear]=1.0;
        
        return true;
    }
    
    if(dec & Decision::kContinue){
        if(isInsideTubeNeigh)
            tubularNeighbourhood[accepted.linear]=1.0;
        
        return queue.empty();
    }
    if(isInsideTubeNeigh)
        tubularNeighbourhood[accepted.linear]=1.0;
    
    const auto offsets = stencilData.ReversedOffsets(accepted);
    for(OffsetCRef offset : offsets){
        ConditionalUpdate(accepted.index, offset, top.value);}

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

template<typename T> bool  HamiltonFastMarching<T>::
IsInsideTubularNeighbourhood(IndexCRef acceptedIndex){
    // if the adaptive neighbourhood is not used
    if(!acuteAnglePreservedShape)
        return true;
    
    if(listSeedIndices.find(acceptedIndex)!=listSeedIndices.end())
        return true;
    VectorType acceptedGeoFlow=normalizedGeoFlows(acceptedIndex);
    VectorType neighGeoFlow=VectorType::Constant(0.0);
    VectorType discreteDirection=VectorType::Constant(0.0);
    const int numNeighs=(int)std::pow((double)Dimension,(double)Dimension)-1;
    for(int k=0; k<numNeighs;k++){
        IndexType neighIndex;
        OffsetType offset=neighbourOffsetList[k];
        for(int idx=0;idx<Dimension;idx++)
            neighIndex[idx]=acceptedIndex[idx]+offset[idx];
        
        if(!normalizedGeoFlows.InRange(neighIndex))
            continue;
        if(!acceptedFlags(neighIndex))
            continue;
        
        if(listSeedIndices.find(neighIndex)!=listSeedIndices.end())
            continue;
        
        neighGeoFlow=normalizedGeoFlows(neighIndex);
        const double scalarProduct=acceptedGeoFlow.ScalarProduct(neighGeoFlow);
        if(scalarProduct<angleScalarProductTolerance){
            if(values(acceptedIndex)>minRadiusTolerance && values(neighIndex)>minRadiusTolerance)
                return false;
        }
        // experimentally setting
        for(int h=0;h<Dimension;h++)
            discreteDirection[h]=(double)(neighIndex[h]-acceptedIndex[h]);
        
        discreteDirection=discreteDirection.Normalized();
        if(scalarProduct<-0.1 && discreteDirection.ScalarProduct(acceptedGeoFlow)<0 && values(acceptedIndex)>1.9){
            return false;
        }
    }
    return true;
}

template<typename T> void
HamiltonFastMarching<T>::ConditionalUpdate(IndexCRef acceptedIndex,
                                           OffsetType offset,
                                           ScalarType acceptedValue){
    FullIndexType updated;
    const auto transform = VisibleOffset(acceptedIndex,-offset, updated.index);
	
    if(!transform.IsValid())
        return;
    transform.PullVector(offset);    
    updated.linear = values.Convert(updated.index);
    
    if(!mask.empty() && mask(updated.index)<0.5)
        return;
    
    if(acceptedFlags[updated.linear])
        return;
    
    // Next line forbids updating of seeds or given data.
    // Can be specialized to allow e.g. for sequential computation of Voronoi diagrams.
    if(activeNeighs[updated.linear].none() &&
	   values[updated.linear]!=Traits::Infinity())
        return;
    if(values[updated.linear]<=acceptedValue)
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


#endif /* HamiltonFastMarching_hxx */
