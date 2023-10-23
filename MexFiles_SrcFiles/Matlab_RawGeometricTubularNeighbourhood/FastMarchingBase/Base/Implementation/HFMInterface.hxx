// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef HFMInterface_hxx
#define HFMInterface_hxx

// ---- Default stencil setup ------

template<typename T> template<typename Dummy> struct
HFMInterface<T>::SpecializationsDefault_<true,Dummy> {
    typedef HFMInterface<T> HFMI;
    typedef typename HFMI::HFM HFM;
    Redeclare6Types(HFM,ScalarType,Traits,MultiplierType,PointType,IndexType,DiscreteType)
    Redeclare2Constants(Traits,mathPi,Dimension)
    
    static const int DimDep =Traits::nStencilDependencies, DimIndep = Dimension-DimDep;
    template<typename E> using DataSource = typename HFM::template DataSource<E>;
    template<typename E> using DataSource_Value = typename HFMI::template DataSource_Value<E>;
    template<typename E> using DataSource_Array = typename HFMI::template DataSource_Array<E>;
    template<typename E> using DataSource_Indep = typename HFMI::template DataSource_Indep<E>;
    template<typename E> using DataSource_Dep = typename HFMI::template DataSource_Dep<E>;
    template<typename E> static std::unique_ptr<DataSource<E> > GetField(KeyCRef name,HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<E>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>(io.template Get<E>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<E, Dimension>(name)));}
        else if(nDims==DimIndep) {  return ResultType(new DataSource_Indep<E>(io.template GetArray<E, DimIndep>(name)));}
        else if(nDims==DimDep) {    return ResultType(new DataSource_Dep<E>(io.template GetArray<E, DimDep>(name)));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.\n");}
    }
    template<typename E> static std::unique_ptr<DataSource<E> > GetIntegralField(KeyCRef name,HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<ScalarType>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>((E)io.template Get<ScalarType>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<ScalarType, Dimension>(name).template Cast<E>()));}
        else if(nDims==DimIndep) {  return ResultType(new DataSource_Indep<E>(io.template GetArray<ScalarType, DimIndep>(name).template Cast<E>()));}
        else if(nDims==DimDep) {    return ResultType(new DataSource_Dep<E>(io.template GetArray<ScalarType, DimDep>(name).template Cast<E>()));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.");}
    }
    static const DiscreteType UnorientedDim = DimIndep==0 ? 1 : DimIndep; // Avoid zero-dimensional
    typedef std::array<ScalarType, UnorientedDim> UnorientedPointType;
    typedef std::array<DiscreteType, UnorientedDim> UnorientedIndexType;
    static UnorientedIndexType StripUnoriented(const IndexType & p){
        UnorientedIndexType q;
        auto depIt = Traits::stencilDependencies.begin();
        auto qIt = q.begin();
        for(int i=0; i<Dimension; ++i){
            if(depIt!=Traits::stencilDependencies.end() && i==*depIt) ++depIt;
            else { *qIt = p[i]; ++qIt;}
        }
        return q;
    }

    static PointType PadUnoriented(const UnorientedPointType & p){
        PointType q=PointType::Constant(0);
        auto depIt = Traits::stencilDependencies.begin();
        auto pIt = p.begin();
        for(int i=0; i<Dimension; ++i) {
            if(depIt!=Traits::stencilDependencies.end() && i==*depIt) ++depIt;
            else {q[i]=*pIt; ++pIt;}
        }
        assert(DimIndep==0 || pIt==p.end());
        return q;
    }
    static void EquivalentPoints(HFMI*that,PointType p, std::vector<PointType> & result){
        const IndexType dims = that->stencil.dims;
        const auto & dom = that->pFM->dom;
        IndexType index=IndexType::Constant(0);
        const auto dep = Traits::stencilDependencies;
        while(true){
            const PointType r = dom.PointFromIndex(index);
            for(DiscreteType i : dep){p[i]=r[i];}
            result.push_back(p);

            // Next index
            int k;
            for(k=0; k<Traits::nStencilDependencies; ++k){
                if(index[dep[k]]<dims[dep[k]]-1) {++index[dep[k]]; break;}
                else {index[dep[k]]=0;}
            }
            if(k==Traits::nStencilDependencies) break;
        }
    }
    static void PadAdimEquiv(HFMI*that, UnorientedPointType p, std::vector<PointType> & equiv){
        equiv.clear();
        EquivalentPoints(that, that->stencil.Param().ADim( PadUnoriented(p) ), equiv);
    }
    
};

template<typename T> template<typename Dummy> struct
HFMInterface<T>::SpecializationsDefault_<false,Dummy> {
    typedef HFMInterface<T> HFMI;
    typedef typename HFMI::HFM HFM;
    Redeclare3Types(HFM,ScalarType,PointType,IndexType)
    Redeclare1Constant(HFM,Dimension)
    
    template<typename E> using DataSource = typename HFM::template DataSource<E>;
    template<typename E> using DataSource_Value = typename HFMI::template DataSource_Value<E>;
    template<typename E> using DataSource_Array = typename HFMI::template DataSource_Array<E>;
    template<typename E> static std::unique_ptr<DataSource<E> > GetField(KeyCRef name, HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<E>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>(io.template Get<E>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<E, Dimension>(name)));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.\n");}
    }
    template<typename E> static std::unique_ptr<DataSource<E> > GetIntegralField(KeyCRef name,HFMI*that) {
        typedef std::unique_ptr<DataSource<E> > ResultType;
        auto & io = that->io;
        const size_t nDims = io.template GetDimensions<ScalarType>(name).size();
        if(nDims==0) {              return ResultType(new DataSource_Value<E>((E)io.template Get<ScalarType>(name)));}
        else if(nDims==Dimension) { return ResultType(new DataSource_Array<E>(io.template GetArray<ScalarType, Dimension>(name).template Cast<E>()));}
        else {ExceptionMacro("Field " << name << " has incorrect depth.");}
    }
    // Dummy functions 
    typedef PointType UnorientedPointType;
    typedef IndexType UnorientedIndexType;
    static UnorientedIndexType StripUnoriented(const IndexType & p){return p;}
    static UnorientedPointType PadUnoriented(const PointType & p){return p;}
    static void EquivalentPoints(HFMI*, PointType, std::vector<PointType> &){assert(false);};
    static void PadAdimEquiv(HFMI*, UnorientedPointType, std::vector<PointType> &){assert(false);}
};


// ------- Data sources ------
template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Inverse : DataSource<E> {
    typedef DataSource<E> Superclass;
    typedef typename Superclass::ReturnType ReturnType;
    std::unique_ptr<Superclass> pSource;
    DataSource_Inverse(std::unique_ptr<Superclass> p):pSource(std::move(p)){}
    
    template<bool,typename> struct Inv_; // Better with C++17 'if constexpr'
	typedef Inv_<std::numeric_limits<ReturnType>::is_specialized, void> Inv;
    template<typename Dummy> struct Inv_<true,Dummy>{
        ReturnType operator()(const ReturnType & t){return 1./t;}};
    template<typename Dummy> struct Inv_<false,Dummy> {
        ReturnType operator()(const ReturnType & t){
            ReturnType s; for(int i=0; i<t.size(); ++i) s[i]=1./t[i]; return s;}};

    virtual ReturnType operator()(const IndexType & index) const {
        return Inv()((*pSource)(index));}
    virtual bool CheckDims(const IndexType & dims) const {
        if(!pSource) return false; return pSource->CheckDims(dims);}
};

template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Value : DataSource<E> {
    typedef typename DataSource<E>::ReturnType ReturnType;
    E value;
    DataSource_Value(E _value):value(std::move(_value)){};
    virtual bool CheckDims(const IndexType &) const {return true;}
    virtual ReturnType operator()(const IndexType &) const {return value;}
};

template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Array : DataSource<E> {
    typedef typename DataSource<E>::ReturnType ReturnType;
    Array<E,Dimension> values;
    DataSource_Array(Array<E, Dimension> && _values):values(std::move(_values)){};
    virtual bool CheckDims(const IndexType & dims) const {
        return dims==values.dims;}
    virtual ReturnType operator()(const IndexType & index) const {return values(index);}
};


template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Indep<E,true> : DataSource<E> {
    typedef typename HFMInterface<T>::HFM HFM;
    typedef typename HFM::template DataSource<E>::ReturnType ReturnType;
    Redeclare2Types(HFM,Traits,IndexType)
    Redeclare1Constant(HFM,Dimension)
    typedef typename HFMInterface<T>::template Array<E, Dimension-Traits::nStencilDependencies> ArrayType;
    ArrayType values;
    DataSource_Indep(ArrayType && _values):values(std::move(_values)){};
    typedef typename ArrayType::IndexType ShortIndexType;
    ShortIndexType Convert(const IndexType & index) const {
        ShortIndexType result;
        int i=0;
        for(int j=0; j<Dimension; ++j){
            if(i<Traits::nStencilDependencies && Traits::stencilDependencies[i]==j){++i;}
            else result[j-i]=index[j];
        }
        return result;
    }
    virtual bool CheckDims(const IndexType & dims) const {return Convert(dims)==values.dims;}
    virtual ReturnType operator()(const IndexType & index) const {return values(Convert(index));}
};

template<typename T> template<typename E> struct // Dummy specialization
HFMInterface<T>::DataSource_Indep<E,false> {};

template<typename T> template<typename E> struct
HFMInterface<T>::DataSource_Dep<E,true> : DataSource<E> {
    typedef typename HFMInterface<T>::HFM HFM;
    typedef typename HFM::template DataSource<E>::ReturnType ReturnType;
    Redeclare2Types(HFM,Traits,IndexType)
    Redeclare1Constant(HFM,Dimension)
    typedef typename HFMInterface<T>::template Array<E, Traits::nStencilDependencies> ArrayType;
    ArrayType values;
    DataSource_Dep(ArrayType && _values):values(std::move(_values)){};
    typedef typename ArrayType::IndexType ShortIndexType;
    ShortIndexType Convert(const IndexType & index) const {
        ShortIndexType result;
        for(int i=0; i<Traits::nStencilDependencies; ++i){
            result[i]=index[Traits::stencilDependencies[i]];}
        return result;
    }
    virtual bool CheckDims(const IndexType & dims) const {return Convert(dims)==values.dims;}
    virtual ReturnType operator()(const IndexType & index) const {return values(Convert(index));}
};

template<typename T> template<typename E> struct // Dummy specialization
HFMInterface<T>::DataSource_Dep<E,false> {};

template<typename T> template<typename E>
struct HFMInterface<T>::TimeDependentSource : DataSource<E> {
    typedef T Traits;
    const ScalarType * pCurrentTime=nullptr;
    typedef std::pair<ScalarType,std::unique_ptr<DataSource<E> > > InterpolationPoint;
    std::vector<InterpolationPoint> interpolationData;
    typedef typename DataSource<E>::ReturnType ReturnType;
    virtual bool CheckDims(const IndexType & index) const {
        // Also checks that times are sorted and pointers well assigned
        if(interpolationData.empty()) return false;
        ScalarType lastTime = -T::Infinity();
        for(const auto & [time,source] : interpolationData){
            if(time<=lastTime) return false;
            lastTime=time;
            if(source==nullptr) return false;
            if(!source->CheckDims(index)) return false;
        }
        return true;
    };
    virtual ReturnType operator()(const IndexType & index) const {
        assert(pCurrentTime!=nullptr);
        const ScalarType t=*pCurrentTime; 
        if(t<=interpolationData.front().first) return (*interpolationData.front().second)(index);
        if(t>=interpolationData.back().first) return (*interpolationData.back().second)(index);
        const auto it1 = std::lower_bound(interpolationData.begin(), interpolationData.end(), *pCurrentTime,
                                          [](const InterpolationPoint & a, ScalarType b)->bool{return a.first<b;});
        const auto it0 = it1-1;
        const ScalarType t0 = it0->first, t1=it1->first;
        const E r0 = (*it0->second)(index), r1=(*it1->second)(index);
        return r0 + ((t-t0)/(t1-t0))*(r1-r0);
    };
};

// ---- Getting fields -------

template<typename T> template<typename E> auto
HFMInterface<T>::GetField(KeyCRef s, bool mayDependOnTime) -> std::unique_ptr<DataSource<E> > {
    if(io.HasField(s)){
    std::unique_ptr<DataSource<E> > result = SpecializationsDefault::template GetField<E>(s,this);
    if(!result->CheckDims(stencil.dims)){
        ExceptionMacro("Error : Field " << s << " has inconsistent dimensions.");}
        return std::move(result);
    } else if(mayDependOnTime && io.HasField(s+"_times")) {
        typedef TimeDependentSource<E> ResultType;
        std::unique_ptr<ResultType> pResult(new ResultType);
        pResult->pCurrentTime = &(pTime->currentTime);
        pTime->currentTime = -Traits::Infinity();
        std::vector<ScalarType> times = io.template GetVector<ScalarType>(s+"_times");
        for(int i=1; i<times.size(); ++i){
            if(times[i-1]>=times[i]) ExceptionMacro("Time dependent field error : times are not strictly increasing.");}
        pResult->interpolationData.reserve(times.size());
        for(int i=0; i<times.size(); ++i){
            pResult->interpolationData.push_back({times[i], GetField<E>(s+"_"+std::to_string(i))});}
        return std::move(pResult);
    } else ExceptionMacro("Error : Field " << s << " not found.");
    
}

template<typename T>
int HFMInterface<T>::FieldElementSize(KeyCRef key) const {
	int size = io.GetElementSize(key, Dimension);
	if(size<0) size = io.GetElementSize(key, 0);
	if(size<0) ExceptionMacro("Error : field " << key << " has invalid depth.");
	return size;
}

template<typename T> template<typename E> auto
HFMInterface<T>::GetIntegralField(KeyCRef s) -> std::unique_ptr<DataSource<E> > {
    std::unique_ptr<DataSource<E> > result = SpecializationsDefault::template GetIntegralField<E>(s,this);
    if(!result->CheckDims(stencil.dims)){
        ExceptionMacro("Error : Field " << s << " has inconsistent dimensions.\n");}
    return std::move(result);}


// ---------- Running  ---------
template<typename T> void HFMInterface<T>::
Run() {
    Run_SetupIO();
	Run_SetupStencil();
	Run_SetupSolver(); 
    Run_SetupExtraAlgorithms();
    if(Run_RunSolver())
        return;
	
	const clock_t top = clock();
//    Run_ExtractGeodesics();
	io.Set<ScalarType>("GeodesicCPUTime",(clock()-top)/double(CLOCKS_PER_SEC));

    Run_ExportData();
}

template<typename T> void HFMInterface<T>::
Run_SetupIO() {
    // Setup input array ordering policy
    io.verbosity = (int)io.Get<ScalarType>("verbosity",io.verbosity);
    
	// "Special case: Internally, the input arrays are converted to RowMajor"
    io.arrayOrdering = enumFromString<IO::ArrayOrdering>(io.GetString("arrayOrdering",enumToRealString(io.arrayOrdering)));
    if(static_cast<int>(io.arrayOrdering)==-1)
        ExceptionMacro("Error : unrecognized array ordering");
    pTime = std::unique_ptr<TimeDependentFields<T> >(new TimeDependentFields<T>);
    pTime->Setup(this);
}

template<typename T> void HFMInterface<T>::
Run_SetupStencil() {
    // Setup stencils (hence metric)
    stencil.Setup(this);
    
    const clock_t top = clock();
    pFM=std::unique_ptr<HFM>(new HFM(stencil));
	
    io.Set<ScalarType>("StencilCPUTime",ScalarType(clock()-top)/CLOCKS_PER_SEC);
}

template<typename T> void HFMInterface<T>::
Run_SetupSolver() {
    // ------- Some exports that are independent of the fast marching results -------
    io.Set<ScalarType>("MaxStencilWidth",pFM->MaxStencilWidth());
    pFM->order = io.template Get<ScalarType>("order",1.);
      
    // set mask if required
    if(io.HasField("physicalMask"))
        pFM->mask = io.GetArray<ScalarType, Dimension>("physicalMask");
    else if(io.HasField("mask"))
        pFM->mask = io.GetArray<ScalarType, Dimension>("mask");
    else if(io.HasField("activeRegion"))
        pFM->mask = io.GetArray<ScalarType, Dimension>("activeRegion");
    
    if(io.HasField("acuteAnglePreservedShape"))
        pFM->acuteAnglePreservedShape=io.template Get<ScalarType>("acuteAnglePreservedShape")>0.5;
    
    if(io.HasField("minRadiusTolerance"))
        pFM->minRadiusTolerance=io.template Get<ScalarType>("minRadiusTolerance")>0.5;
    

    if(io.HasField("getStencils")) {
        const auto & pts = io.GetVector<PointType>("getStencils");
        std::ostringstream oss;
        oss << "{";
        for(const PointType & p : pts){
            IndexType index = pFM->dom.IndexFromPoint(stencil.Param().ADim(p));
            if(!pFM->dom.PeriodizeNoBase(index).IsValid()){
                WarnMsg() << "GetStencil error : point " << p << "is out of range.\n";
                continue;}
            typename HFM::StencilType indStencil;
            stencil.SetStencil(index,indStencil);
            oss << "{" << index << "," << indStencil << "},";
        }
        oss << "}";
        io.SetString("stencils", oss.str());
    }
    
    if(io.HasField("pointToIndex")){
        auto pts = io.GetVector<PointType>("pointToIndex");
        for(PointType & p : pts)
            p=PointType::CastCoordinates(pFM->dom.IndexFromPoint(stencil.Param().ADim(p)));
        io.SetVector("indexFromPoint", pts);
    }
    if(io.HasField("indexToPoint")){
        auto pts = io.GetVector<PointType>("indexToPoint");
        for(PointType & p : pts) p=stencil.Param().ReDim(pFM->dom.PointFromIndex(IndexType::CastCoordinates(p)));
        io.SetVector("pointFromIndex", pts);
    }

	pFM->factoring.Setup(this); // Sets the seeds and factor
    
    // Reinserting data from previous run
    if(io.HasField("values")){
        pFM->values = io.GetArray<ScalarType, Dimension>("values");}
    
    if(io.HasField("activeNeighs")){
        pFM->activeNeighs = io.GetArray<ScalarType, Dimension>("activeNeighs")
		.template Cast<ActiveNeighFlagType>();}
}

template<typename T> template<typename Alg> Alg* HFMInterface<T>::
SetupSingleAlgorithm(){
    auto pAlg = std::unique_ptr<Alg>(new Alg);
    pAlg->Setup(this);
    if(! pAlg->ImplementIn(pFM.get())) return nullptr;
    Alg * result = pAlg.get();
    extras.push_back(std::move(pAlg));
    return result;
}

template<typename T> bool HFMInterface<T>::
Run_RunSolver() {
	if(io.HasField("values") && io.HasField("activeNeighs")){
		if(io.verbosity>=1) Msg() << "Bypassing fast marching solver based on cached data.\n";
		std::fill(pFM->acceptedFlags.begin(),pFM->acceptedFlags.end(),true);
	} else if(pFM->seeds.empty()) {
		WarnMsg() << "No seeds to run fast marching solver,"
		<< " and missing values or activeNeighs to compute geodesics.\n";
		return true;
    } else {
        const clock_t top = clock();
        pFM->Run();
        const clock_t elapsed = clock()-top;
        const ScalarType FMCPUTime = ScalarType(elapsed)/CLOCKS_PER_SEC;
        io.Set<ScalarType>("FMCPUTime",FMCPUTime);
		if(io.verbosity>=1){
			Msg() << "Fast marching solver completed in " << FMCPUTime << " s.\n";}
    }
    return false;
}

template<typename T> void HFMInterface<T>::
ExportGeodesics(KeyCRef suffix, const std::vector<PointType> & tips){
    if(pGeodesicSolver==nullptr){
        typedef std::unique_ptr<GeodesicSolverInterface> GeoSolverPtr;
        std::string geodesicSolverType =
        HFM::DomainType::periodizeUsesBase ? "ODE" : 
        io.GetString("geodesicSolver", "Discrete");
        if(geodesicSolverType=="Discrete"){
            pGeodesicSolver = GeoSolverPtr(new GeodesicDiscreteSolver<Traits>(*pFM));
        }
        else if(geodesicSolverType=="ODE"){
            pGeodesicSolver = GeoSolverPtr(new GeodesicODESolver<Traits>(*pFM));
        }
        else{
            ExceptionMacro("Error : Unrecognized geodesic solver " << geodesicSolverType << ".\n");
        }
        pGeodesicSolver->Setup(this);
    }
    
    const auto & geodesics = pGeodesicSolver->Run(this,tips);

    std::vector<ScalarType> geodesicLengths;
    for(const auto & geo : geodesics) {geodesicLengths.push_back((ScalarType)geo.size());}
    std::vector<PointType> geodesicPoints;
    geodesicPoints.reserve((size_t)std::accumulate(geodesicLengths.begin(), geodesicLengths.end(), 0.));
   
    for(const auto & geo : geodesics)
        for(const PointType & p : geo){
            PointType newPoint;
            for(int k=0;k<Dimension;k++) //HFM::Traits::nStencilDependencies
                newPoint[k]=p[k]-0.5;
            geodesicPoints.push_back(stencil.Param().ReDim(newPoint));
        }
	const std::string sep = suffix.empty() ? "" : "_";
    io.SetVector("geodesicPoints"+sep+suffix, geodesicPoints);
    io.SetVector("geodesicLengths"+sep+suffix, geodesicLengths);
}

template<typename T> void HFMInterface<T>::
Run_ExtractGeodesics() {
    
    if(!io.HasField("stopWhenAnyAccepted")){
        if(io.HasField("tips")){
            std::vector<PointType> tips = io.GetVector<PointType>("tips");
            for(PointType & tip : tips){
                PointType tipTem = stencil.Param().ADim(tip);
                IndexType tipIdx;
                for (int k=0;k<Dimension;k++)
                    tipIdx[k]=tipTem[k];
                
                tip=this->pFM->dom.PointFromIndex(tipIdx);
            }
            ExportGeodesics("",tips);
        }
    }
    
    if(HFM::hasBundle && io.HasField("tips_Unoriented")){
        typedef typename SpecializationsDefault::UnorientedPointType UnorientedPointType;
        const auto indepTips = io.GetVector<UnorientedPointType>("tips_Unoriented");
        
        std::vector<PointType> tips;
        std::vector<PointType> equiv;
        for(const UnorientedPointType & indepTip : indepTips){
            SpecializationsDefault::PadAdimEquiv(this,indepTip,equiv);
            ScalarType valMin=Traits::Infinity();
            PointType pMin=equiv.front();
            for(const PointType & q : equiv){
                PointType r=q;
                if(!pFM->dom.PeriodizeNoBase(r).IsValid()) continue; // Out of range
                const ScalarType val =pFM->values(pFM->dom.IndexFromPoint(r));
                if(valMin<val) continue; // Too large value
                pMin=q; valMin=val;
            }
            tips.push_back(pMin);
        }
        ExportGeodesics("Unoriented", tips);
    }
}


template<typename T> void HFMInterface<T>::
Run_ExportData() {
    if(io.Get<ScalarType>("exportValues",0.,2)) {
        io.SetArray("values", pFM->values);}
    
    io.SetArray("tubularNeighbourhood", pFM->tubularNeighbourhood);
    
    if(io.Get<ScalarType>("exportActiveNeighs",0.,3)) {
		io.SetArray("activeNeighs",pFM->activeNeighs.template Cast<ScalarType>());}
    
    if(io.Get<ScalarType>("exportGeodesicFlow",0.,2)) {
        Array<VectorType, Dimension> flow;
        flow.dims=pFM->values.dims;
        flow.resize(pFM->values.size());
        for(DiscreteType i=0; i<flow.size(); ++i){
            flow[i] = pFM->GeodesicFlow(flow.Convert(i)).flow;
            flow[i] = - stencil.Param().ReDim(flow[i]); // Flows from seeds to tips.
        }
        io.SetArray<VectorType>("geodesicFlow", flow);
    }
	if(io.Get<ScalarType>("exportActiveOffsets",0.,3)) {
		// Possible improvement : do simultaneously with exportGeodesicFlow
		// to avoid recomputing twice.
		Redeclare3Types(HFM,DiscreteFlowType,OffsetType,ShortType)
		using ElemType = std::array<ScalarType,DiscreteFlowType::max_size()>;
		Array<ElemType,Dimension> flow;
        flow.dims=pFM->values.dims;
        flow.resize(pFM->values.size());
        for(DiscreteType i=0; i<flow.size(); ++i){
			DiscreteFlowType discreteFlow;
			pFM->Recompute(flow.Convert(i),discreteFlow);
			ElemType & elem = flow[i];
			
			for(int j=0; j<discreteFlow.size(); ++j){
				const OffsetType & offset = discreteFlow[j].offset;
				unsigned long offset_tolong = 0;
				for(int k=0; k<Dimension; ++k){
					offset_tolong += std::make_unsigned_t<ShortType>(offset[k]) << (8*sizeof(ShortType)*(Dimension-k-1));}
				elem[j] = ScalarType(offset_tolong);
			}
			for(int j=(int)discreteFlow.size(); j<discreteFlow.max_size(); ++j) {elem[j]=0;}
        }
        io.SetArray<ElemType>("activeOffsets", flow);
	}
    for(const auto & pAlg : extras) pAlg->Finally(this);
}

#endif /* HFMInterface_hxx */
