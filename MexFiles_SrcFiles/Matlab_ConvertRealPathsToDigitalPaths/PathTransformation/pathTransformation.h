/*

This file is an adaption form Miller's work: <piecewisegeodesiccombination.h>.
 <piecewisegeodesiccombination.h> is a part of
[MBC15] J. Mille, S. Bougleux and L. Cohen. Combination of piecewise-geodesic
        paths for interactive segmentation. International Journal of Computer
        Vision, 112(1):1-22, 2015.
*/

#ifndef _PATHTTRANSFORMATION_H_
#define _PATHTTRANSFORMATION_H_

#include <list>
#include "path_Simplified.h"


// class PathTransformationBase
// Contains pure virtual member functions and thus cannot be instantiated
class PathTransformationBase
{
  // Member variables
public:
    
    bool bLoopRemoval;
    
  protected:
    // Pointer to input image
	// Initialized as null in default constructor. Needs to be set by calling AttachImage()
    int numRows;
    int numCols;
    
	CRealPath2D realPath;
	CIntegerPath2D integerPath;
    
  public:
    
    bool bShowInfo;
    

  public:
	// Default constructor
    PathTransformationBase();

    // Destructor
    virtual ~PathTransformationBase() {}

    virtual void SetDomainSize(const int, const int);
    virtual void SetInput(const CRealPath2D& );
    virtual void Update();
    // Create 4-connected simple digital path from a real path
    // Possible loops are eliminated
    // Params: real-valued curve [in], digital 4-connected simple curve [out]
    void ConvertRealPathToDigitalPathSimple(const CRealPath2D &, CIntegerPath2D &) const;

    // Create 4-connected digital path from a real path
    // The path is not necessarily simple: it can have self-tangencies or self-intersections
    // Params: real-valued curve [in], digital 4-connected curve [out]
    void ConvertRealPathToDigitalPath(const CRealPath2D &, CIntegerPath2D &) const;

    virtual const CIntegerPath2D& GetConvertedDigitalPath() const;
    
};
#include "pathTransformation.cpp"
#endif
