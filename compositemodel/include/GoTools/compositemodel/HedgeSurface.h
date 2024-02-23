/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _HEDGESURFACE_H
#define _HEDGESURFACE_H

#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/ClassType.h"
//#include "GoTools/compositemodel/RevEngRegion.h"

namespace Go {

class RevEngRegion;
  class CurveOnSurface;

  enum
    {
     SURF_TYPE_UNDEF, LINEARSWEPT_SURF, ROTATIONALSWEPT_SURF
    };
  
class HedgeSurface : public ftSurface
{
public:
  
  // Constructor
  // HedgeSurface();

  HedgeSurface();

  HedgeSurface(shared_ptr<ParamSurface> sf, RevEngRegion *region);

  HedgeSurface(shared_ptr<ParamSurface> sf, std::vector<RevEngRegion*>& region);

  // Destructor
  ~HedgeSurface();

  void setLinearSweepInfo(shared_ptr<SplineCurve> profile,
			  Point startpt, Point endpt)
  {
    surf_code_ = LINEARSWEPT_SURF;
    profile_ = profile;
    sweep1_ = startpt;
    sweep2_ = endpt;
  }
  
  void setRotationalSweepInfo(shared_ptr<SplineCurve> profile,
			      Point location, Point axis)
  {
    surf_code_ = ROTATIONALSWEPT_SURF;
    profile = profile_;
    sweep1_ = location;
    sweep2_ = axis;
  }
  
  int dimension()
  {
    return surface()->dimension();
  }
  
  int numPoints();

  ClassType instanceType(int& code);

  bool isPlane()
  {
    int code;
    return (instanceType(code) == Class_Plane);
  }

  bool isCylinder()
  {
    int code;
    return (instanceType(code) == Class_Cylinder);
  }

  bool isSphere()
  {
    int code;
    return (instanceType(code) == Class_Sphere);
  }

  bool isTorus()
  {
    int code;
    return (instanceType(code) == Class_Torus);
  }

  bool isCone()
  {
    int code;
    return (instanceType(code) == Class_Cone);
  }

  bool isSpline()
  {
    int code;
    return (instanceType(code) == Class_SplineSurface);
  }

  std::vector<RevEngRegion*> getRegions()
  {
    return regions_;
  }

  int numRegions()
  {
    return (int)regions_.size();
  }

  RevEngRegion* getRegion(int ix)
  {
    if (ix < 0 || ix >= (int)regions_.size())
      return 0;
    else
      return regions_[ix];
  }

  void addRegion(RevEngRegion* reg);
  
  bool removeRegion(RevEngRegion* reg);
  
  BoundingBox regionsBox()
  {
    return bbox_;
  }

  bool updateSurfaceWithAxis(Point axis[3], int ix, double tol, double angtol);
    
  bool updatePlaneWithAxis(Point axis[3], int ix, double tol, double angtol);
    
  bool updateCylinderWithAxis(Point axis[3], int ix, double tol, double angtol);
    
  bool checkAccuracyAndUpdate(shared_ptr<ParamSurface> surf, double tol,
			      double angtol);

  bool isCompatible(HedgeSurface* other, double angtol, double approx_tol,
		    ClassType& type, double& score);

  bool hasBaseSf();

  void ensureSurfaceBounded();

  bool isTangential(HedgeSurface* surf);

  void doTrim(std::vector<shared_ptr<CurveOnSurface> >& int_cvs,
	      shared_ptr<BoundedSurface>& bdsf,
	      double tol,
	      std::vector<shared_ptr<HedgeSurface> >& added_sfs);
  void limitSurf(double diag = -1.0);
  void trimWithPoints(double aeps);

  void store(std::ostream& os) const;
  void read(std::istream& is);
    
private:
  std::vector<RevEngRegion*> regions_;
  BoundingBox bbox_;
  int surf_code_;
  shared_ptr<SplineCurve> profile_;
  Point sweep1_;
  Point sweep2_;
};
}

#endif // _HEDGESURFACE_H
