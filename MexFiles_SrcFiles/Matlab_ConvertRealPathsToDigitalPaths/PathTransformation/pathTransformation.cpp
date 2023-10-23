/*
This file is an adaption form Miller's work: <piecewisegeodesiccombination.cpp>.
 <piecewisegeodesiccombination.cpp> is a part of
[MBC15] J. Mille, S. Bougleux and L. Cohen. Combination of piecewise-geodesic
        paths for interactive segmentation. International Journal of Computer
        Vision, 112(1):1-22, 2015.
*/

#include <iostream>
#include "path_Simplified.h"
#include "path_Simplified.cpp"
#include "pathTransformation.h"

PathTransformationBase::PathTransformationBase()
{
    numRows=1;
    numCols=1;
    bLoopRemoval=false;
}



void PathTransformationBase::ConvertRealPathToDigitalPathSimple(const CRealPath2D &realPath, CIntegerPath2D &intPath) const{
    CRealPath2D::const_iterator itRealPoint;
    const CCouple<int> *pPoint;
    CArray2D<CCouple<int> *> arrayPtrPoints;
    CCouple<float> pfHalf(0.5f, 0.5f);
    CCouple<int> piLastAdded, piCurrent;
    CCouple<int> viAbsDiff;

    arrayPtrPoints.Init(numCols, numRows);
    arrayPtrPoints.Fill(NULL);

    intPath.clear();

    itRealPoint = realPath.begin();
    piLastAdded = (CCouple<int>)(*itRealPoint + pfHalf);
    if (arrayPtrPoints.Element(piLastAdded)==NULL){
        intPath.push_back(piLastAdded);
        arrayPtrPoints.Element(piLastAdded) = &intPath.back();
    }

    for (itRealPoint++; itRealPoint!=realPath.end(); itRealPoint++){
        piCurrent = (CCouple<int>)(*itRealPoint + pfHalf);
        if (piCurrent!=piLastAdded)
        {
            viAbsDiff.x = abs(piCurrent.x-piLastAdded.x);
            viAbsDiff.y = abs(piCurrent.y-piLastAdded.y);

            if (viAbsDiff.x+viAbsDiff.y==1)
            {
                if (arrayPtrPoints.Element(piCurrent)==NULL)
                {
                    intPath.push_back(piCurrent);
                    arrayPtrPoints.Element(piCurrent) = &intPath.back();
                    piLastAdded = piCurrent;
                }
                else {
                    pPoint = arrayPtrPoints.Element(piCurrent);
                    while (pPoint!=&intPath.back())
                    {
                        arrayPtrPoints.Element(intPath.back()) = NULL;
                        intPath.pop_back();
                    }
                    piLastAdded = piCurrent;
                }
            }
            else {
                // Draw digital line from piLastAdded to piLastCurrent
                CIntegerPath2D pathDigitalLineSeg;
                CIntegerPath2D::const_iterator itPointSeg;

                if (viAbsDiff.x==1 && viAbsDiff.y==1)
                {
                    pathDigitalLineSeg.push_back(piLastAdded);
                    pathDigitalLineSeg.push_back(CCouple<int>(piLastAdded.x, piCurrent.y));
                    pathDigitalLineSeg.push_back(piCurrent);
                }
                else {
                    pathDigitalLineSeg.SetDigitalLineSegment4connected(piLastAdded, piCurrent);
                }

                for (itPointSeg=++pathDigitalLineSeg.begin(); itPointSeg!=pathDigitalLineSeg.end(); itPointSeg++)
                {
                    piCurrent = *itPointSeg;
                    if (arrayPtrPoints.Element(piCurrent)==NULL)
                    {
                        intPath.push_back(piCurrent);
                        arrayPtrPoints.Element(piCurrent) = &intPath.back();
                        piLastAdded = piCurrent;
                    }
                    else {
                        pPoint = arrayPtrPoints.Element(piCurrent);
                        while (pPoint!=&intPath.back())
                        {
                            arrayPtrPoints.Element(intPath.back()) = NULL;
                            intPath.pop_back();
                        }
                        piLastAdded = piCurrent;
                    }
                }
            }
        }
    }
}

void PathTransformationBase::ConvertRealPathToDigitalPath(const CRealPath2D &realPath, CIntegerPath2D &intPath) const{
    CRealPath2D::const_iterator itRealPoint;
    CCouple<float> pfHalf(0.5f, 0.5f);
    CCouple<int> piLastAdded, piCurrent;
    CCouple<int> viAbsDiff;

    intPath.clear();

    itRealPoint = realPath.begin();
    piLastAdded = (CCouple<int>)(*itRealPoint + pfHalf);
    intPath.push_back(piLastAdded);

    for (itRealPoint++; itRealPoint!=realPath.end(); itRealPoint++)
    {
        piCurrent = (CCouple<int>)(*itRealPoint + pfHalf);
        if (piCurrent!=piLastAdded)
        {
            viAbsDiff.x = abs(piCurrent.x-piLastAdded.x);
            viAbsDiff.y = abs(piCurrent.y-piLastAdded.y);

            if (viAbsDiff.x+viAbsDiff.y==1)
            {
                intPath.push_back(piCurrent);
                piLastAdded = piCurrent;
            }
            else {
                // Draw digital line from piLastAdded to piLastCurrent
                CIntegerPath2D pathDigitalLineSeg;
                CIntegerPath2D::const_iterator itPointSeg;

                if (viAbsDiff.x==1 && viAbsDiff.y==1)
                {
                    pathDigitalLineSeg.push_back(piLastAdded);
                    pathDigitalLineSeg.push_back(CCouple<int>(piLastAdded.x, piCurrent.y));
                    pathDigitalLineSeg.push_back(piCurrent);
                }
                else {
                    pathDigitalLineSeg.SetDigitalLineSegment4connected(piLastAdded, piCurrent);
                }

                for (itPointSeg=++pathDigitalLineSeg.begin(); itPointSeg!=pathDigitalLineSeg.end(); itPointSeg++)
                {
                    piCurrent = *itPointSeg;
                    intPath.push_back(piCurrent);
                    piLastAdded = piCurrent;
                }
            }
        }
    }

}

void PathTransformationBase::SetDomainSize(const int nRows, const int nCols){
    numRows=nRows;
    numCols=nCols;
}

void PathTransformationBase::SetInput(const CRealPath2D& inputRealPath){
    realPath=CRealPath2D(inputRealPath);
}

void PathTransformationBase::Update(){
    
    if (realPath.empty())
        std::cerr<<"ConvertRealPathsToDigitalPaths: No real path is given. "<<std::endl;
    
    if(bLoopRemoval){
        ConvertRealPathToDigitalPathSimple(realPath, integerPath);
    }
    else{
        ConvertRealPathToDigitalPath(realPath, integerPath);
    }
}

const CIntegerPath2D& PathTransformationBase::GetConvertedDigitalPath() const{
    return integerPath;
}

