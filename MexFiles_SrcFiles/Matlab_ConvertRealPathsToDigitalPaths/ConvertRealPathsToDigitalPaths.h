#ifndef _ConvertRealPathsToDigitalPaths_H_
#define _ConvertRealPathsToDigitalPaths_H_

#include <limits>
#include <iostream>
#include "couple.h"
#include "triplet.h"
#include "path_Simplified.h"


template<typename IO>
int Run(IO& io){
    
    typedef typename LinearAlgebra::Array<double, 2>  ImageType;
    typedef typename LinearAlgebra::Point<double, 2>  PointType;
    typedef typename LinearAlgebra::Point<int, 2>     IndexType;
    typedef PathTransformationBase                    PathTransType;
    
    std::vector<PointType> realPath = io.template GetVector<PointType>("realPath");
 
    CRealPath2D realPath_;
    typename std::vector<PointType>::const_iterator itRealPoint;
    for(itRealPoint=realPath.begin();itRealPoint!=realPath.end();itRealPoint++){
        CCouple<float> pt_;
        PointType pt=*itRealPoint;
        pt_.x=(float)pt[0];
        pt_.y=(float)pt[1];
        realPath_.push_back(pt_);
    }
    
    
    int numRows;
    int numCols;
    std::vector<double> imageSize = io.template GetVector<double>("imageSize");
    if(imageSize.size()==1){
        numRows=(int)imageSize[0];
        numCols=(int)imageSize[0];
    }
    else if(imageSize.size()>1){
        numRows=(int)imageSize[0];
        numCols=(int)imageSize[1];
    }
    
    
    const bool  bLoopRemoval=io.template Get<double>("loopRemoval",0.)>0.5;
    std::unique_ptr<PathTransType> pPathTrans = std::unique_ptr<PathTransType>(new PathTransType());
   
    pPathTrans->SetInput(realPath_);
    pPathTrans->bLoopRemoval=bLoopRemoval;
    pPathTrans->SetDomainSize(numRows, numCols);
    pPathTrans->Update();
    CIntegerPath2D digitalPath_=pPathTrans->GetConvertedDigitalPath();
    std::vector<PointType> digitalPath;
    CIntegerPath2D::const_iterator itIntPoint;
    for(itIntPoint=digitalPath_.begin();itIntPoint!=digitalPath_.end();itIntPoint++){
        CCouple<int> intPoint=*itIntPoint;
        PointType pt;
        pt[0]=(double)intPoint.x;
        pt[1]=(double)intPoint.y;
        digitalPath.push_back(pt);
    }
    io.SetVector("digitalPath",digitalPath);

    
    IndexType dims;
    dims[0]=numCols;
    dims[1]=numRows;
    ImageType binaryImage;
    binaryImage.dims=dims;
    binaryImage.resize(numRows*numCols,0.0);
    std::vector<PointType>::const_iterator itPoint;
    for(itPoint=digitalPath.begin();itPoint!=digitalPath.end();itPoint++){
        IndexType index= IndexType::CastCoordinates(*itPoint);

        if(binaryImage.InRange(index))
            binaryImage(index)=1.0;
    }
    io.SetArray("binaryImage", binaryImage);
    
    return EXIT_SUCCESS;
}
#endif
