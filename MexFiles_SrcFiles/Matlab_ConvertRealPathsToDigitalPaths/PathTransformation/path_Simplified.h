/*
Copyright 2015 Julien Mille

This file is simplification of the file <path.h>.
 <path.h> is a part of the following work:

[MBC15] J. Mille, S. Bougleux and L. Cohen. Combination of piecewise-geodesic
        paths for interactive segmentation. International Journal of Computer
        Vision, 112(1):1-22, 2015.

*/

#ifndef _PATH_SIMPLIFIED_H_
#define _PATH_SIMPLIFIED_H_

#include <list>

#include "couple.h"
#include "arraynd.h"

using namespace std;

class CIntegerPath2D;
class CRealPath2D;

class CIntegerPath2D : public list<CCouple<int> >
{
  friend class CRealPath2D;

  // Member variables
  public:
    float fSortingValue;

  // Member functions
  public:
	CIntegerPath2D(){}
	CIntegerPath2D(const CRealPath2D &);

    // CIntegerPath2D is used to implement the discretized 4-connected admissible paths
    // Before searching combination of admissible paths, we sort them in admissible
    // sets with respect to some value (exteriority). We need to overload comparison
    // operators to allow the STL sort() funtion
    bool operator <(const CIntegerPath2D &path) const {
        return fSortingValue<path.fSortingValue;
    }
    bool operator >(const CIntegerPath2D &path) const {
        return fSortingValue>path.fSortingValue;
    }

	float GetLength() const;

	// Compute area with discrete Green's theorem (should be a digital 4-connected closed curve)
	float GetArea4connected() const;

    void SetDigitalLineSegment4connected(const CCouple<int> &, const CCouple<int> &);
    void AddDigitalPath4connected(const CIntegerPath2D &);

};

class CRealPath2D : public list<CCouple<float> >
{
  friend class CIntegerPath2D;

  // Member variables
  public:
    float fSortingValue;

  // Member functions
  public:
	CRealPath2D()
	{
//		rgbVertices.Set(0,255,0);
//		rgbEdges.Set(0,0,0);
//		bDisplayVertices = true;
//		bDisplayEdges = true;
	}
	CRealPath2D(const CIntegerPath2D &);

    // Admissible paths are implemented in curve with real point coordinates
    // Before searching combination of admissible paths, we sort them in admissible
    // sets with respect to some value (exteriority). We need to overload comparison
    // operators to allow the STL sort() funtion
    bool operator <(const CRealPath2D &path) const {
        return fSortingValue<path.fSortingValue;
    }
    bool operator >(const CRealPath2D &path) const {
        return fSortingValue>path.fSortingValue;
    }

	float GetLength() const;
	void InitLine(const CCouple<float> &, const CCouple<float> &);
	void LaplacianSmooth();

	void Translate(const CCouple<float> &);

    float GetSignedAreaWithLineSegment() const;

    CRealPath2D Resampled(float) const;

};

#endif
