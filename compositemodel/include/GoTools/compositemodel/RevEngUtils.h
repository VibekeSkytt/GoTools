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

#ifndef _REVENGUTILS_H
#define _REVENGUTILS_H

#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/BoundingBox.h"
#include <vector>

namespace Go {

  class RevEngPoint;
  class Line;
  class Plane;
  class Cylinder;
  class Torus;
  class Sphere;
  class Cone;
  
  namespace RevEngUtils
  {
    void principalAnalysis(Point& curr, std::vector<Point>& points, 
			   double lambda[3], double eigenvec[3][3]);

    void principalAnalysis(std::vector<RevEngPoint*>& points, 
			   double lambda[3], double eigenvec[3][3]);

    void TaubinCurvature(Point curr, std::vector<Point>& points,
			 Point& tvec, Point& normal, Point& mincvec,
			 double& minc, Point& maxcvec, double& maxc);

    void computeMonge(Point& curr, std::vector<Point>& points,
		      Point& vec1, Point& vec2, Point& normal, Point& mincvec,
		      double& minc, Point& maxcvec, double& maxc,
		      double& currdist, double& avdist);

    void smoothSurf(shared_ptr<SplineSurface>& surf, int fixed);

    shared_ptr<SplineSurface> surfApprox(std::vector<double>& data, int dim,
					 std::vector<double>& param, int order1,
					 int order2, int nmb_coef1, int nmb_coef2,
					 bool close1, bool close2,
					 int max_iter, double tol, double& maxd, 
					 double& avd, int& num_out,
					 std::vector<double>& parvals,
					 double del=0.0);
    
    shared_ptr<SplineSurface> surfApprox(std::vector<double>& data, int dim,
					 std::vector<double>& param, int order1,
					 int order2, int nmb_coef1, int nmb_coef2,
					 double del=0.0);

    shared_ptr<SplineSurface> surfApprox(std::vector<double>& data, int dim,
					 std::vector<double>& param, int order1,
					 int order2, int nmb_coef1, int nmb_coef2,
					 double umin, double umax, double vmin,
					 double vmax);

    void parameterizeWithPlane(std::vector<RevEngPoint*>& pnts, const BoundingBox& bbox,
			       const Point& vec1, const Point& vec2,
			       std::vector<double>& data, std::vector<double>& param);
    
    void parameterizeWithPlane(std::vector<Point>& pnts, const BoundingBox& bbox,
			       const Point& vec1, const Point& vec2,
			       std::vector<double>& data, std::vector<double>& param);
    bool parameterizeOnPrimary(std::vector<RevEngPoint*>& points,
			       shared_ptr<ParamSurface> surf,
			       std::vector<double>& data, 
			       std::vector<double>& param,
			       int& inner1, int& inner2, bool& close1, bool& close2);

  void computeAxis(std::vector<Point>& points,
		     Point& axis, Point& Cx, Point& Cy);

    void computeAxis(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> >& points,
		     Point& axis, Point& Cx, Point& Cy);

    void
    computeSphereProp(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		      std::vector<RevEngPoint*>::iterator> >& points,
		      Point& centre, double& radius);

    void coneAxis(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> >& points,
		     Point& axis, Point& Cx, Point& Cy);

    void coneApex(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		  std::vector<RevEngPoint*>::iterator> >& points,
		  Point axis, Point& apex, double& phi);

    void computeCylPosRadius(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			     std::vector<RevEngPoint*>::iterator> >& points,
			     Point& low, Point& high, Point& axis, Point& Cx, 
			     Point& Cy, Point& pos, double& radius);

    void computeCircPosRadius(std::vector<Point>& points,
			      const Point& axis, const Point& Cx, const Point& Cy,
			      Point& pos, double& radius);
    
    void computeRadius(std::vector<Point>& points,
		       Point& axis, Point& Cx, Point& Cy, double& radius);
    
    void computePlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		      std::vector<RevEngPoint*>::iterator> >& points,
		      Point& pos, Point& norm);

    void computeLine(vector<Point>& points, Point& pos, Point& dir);
    
    void computePlane(std::vector<Point>& points, Point normal,Point mainaxis[3],
		      Point& pos, Point& norm, Point& Cx, Point& Cy);

    void rotateToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		       std::vector<RevEngPoint*>::iterator> >& points,
		       Point& xvec, Point& axis, Point& mid, std::vector<Point>& rotated);

    void projectToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			std::vector<RevEngPoint*>::iterator> >& points,
			Point& axis, Point& mid, std::vector<Point>& projected,
			double& maxdist, double& avdist, double dlen=-1.0);

    void projectToPlane(std::vector<RevEngPoint*>& points,
			Point& axis, Point& mid, std::vector<Point>& projected,
			double& maxdist, double& avdist, double dlen=-1.0);

     void rotateToPlane(std::vector<Point>& points,
			Point& xvec, Point& axis, Point& mid, std::vector<Point>& rotated);

   void distToSurf(std::vector<RevEngPoint*>::iterator start,
		    std::vector<RevEngPoint*>::iterator end,
		    shared_ptr<ParamSurface> surf, double tol,
		   double& maxdist, double& avdist,
		   int& inside, int& inside2,
		   std::vector<RevEngPoint*>& in,
		   std::vector<RevEngPoint*>& out,
		   std::vector<double>& parvals,
		   std::vector<std::pair<double,double> >& distang,
		   double angtol=-1.0);
    
    void distToSurf(std::vector<Point>& points,
		    shared_ptr<ParamSurface> surf, double tol,
		    double& maxdist, double& avdist, int& inside,
		    std::vector<double>& distance);
    
     void distToSurf(std::vector<RevEngPoint*>& points,
		    shared_ptr<ParamSurface> surf, double tol,
		     double angtol, double& maxdist, double& avdist, 
		     int& inside, int& inside2,
		     std::vector<std::pair<double,double> >& dist_ang);
    
   void distToCurve(std::vector<Point>& points,
		    shared_ptr<ParamCurve> curve, double tol,
		    double& maxdist, double& avdist, int& inside);

     void distToCurve(std::vector<Point>& points,
		      shared_ptr<ParamCurve> curve, double tol,
		      double& maxdist, double& avdist, int& inside,
		      std::vector<double>& dist);

     void distToCurve(std::vector<Point>& points,
		      shared_ptr<ParamCurve> curve, double tol,
		      double& maxdist, double& avdist, int& inside,
		      std::vector<double>& parvals, std::vector<double>& dist);

    shared_ptr<ElementarySurface> elemsurfWithAxis(shared_ptr<ElementarySurface> sf_in,
						 std::vector<RevEngPoint*>& points,
						   Point mainaxis[3], double diag);
    
   shared_ptr<Plane> planeWithAxis(std::vector<RevEngPoint*>& points,
				    Point axis, Point init_loc,
				    Point mainaxis[3]);
    
    shared_ptr<Cylinder> cylinderWithAxis(std::vector<RevEngPoint*>& points,
					  Point axis, Point low, 
					  Point high, Point mainaxis[3]);
    
    shared_ptr<Cylinder> cylinderWithAxis(std::vector<RevEngPoint*>& points,
					  Point axis, Point Cx, Point pos);
    
    shared_ptr<Torus> torusWithAxis(std::vector<RevEngPoint*>& points,
				    Point axis, Point loc, 
				    Point mainaxis[3]);
    
    shared_ptr<Torus> torusWithAxis(std::vector<RevEngPoint*>& points,
				    Point axis, Point Cx, Point pos);
    
    shared_ptr<Sphere> sphereWithAxis(std::vector<RevEngPoint*>& points,
				      Point axis, 
				      Point mainaxis[3]);
    
    shared_ptr<Sphere> sphereWithAxis(std::vector<RevEngPoint*>& points,
				      Point axis, Point Cx, Point pos);
    
    shared_ptr<Cone> coneWithAxis(vector<RevEngPoint*>& points,
				  Point axis, Point low, 
				  Point high, Point mainaxis[3]);
    
    shared_ptr<Cone> coneWithAxis(std::vector<RevEngPoint*>& points,
				      Point axis, Point Cx, Point pos,
				      double len);
    
    shared_ptr<ParamSurface> doMergePlanes(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					   std::vector<RevEngPoint*>::iterator> > points,
					   const BoundingBox& bbox,
					   std::vector<int>& nmbpts,
					   bool set_bound = true);
    shared_ptr<ParamSurface> doMergeCylinders(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					      std::vector<RevEngPoint*>::iterator> > points,
					      const BoundingBox& bbox,
					      std::vector<int>& nmbpts,
					      bool set_bound = true);
    shared_ptr<ParamSurface> doMergeSpheres(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					    std::vector<RevEngPoint*>::iterator> > points,
					    const BoundingBox& bbox,
					    std::vector<int>& nmbpts, Point& normal);
    
    shared_ptr<ParamSurface> doMergeTorus(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					  std::vector<RevEngPoint*>::iterator> > points,
					  const BoundingBox& bbox,
					  std::vector<int>& nmbpts);

    shared_ptr<SplineCurve> midCurve(shared_ptr<SplineCurve>& cv1,
				     shared_ptr<SplineCurve>& cv2);
    
    void  curveApprox(std::vector<Point>& points,
		      shared_ptr<ParamCurve> cvin,
		      int ik, int in, 
		      shared_ptr<SplineCurve>& curve);
    
    void  curveApprox(std::vector<Point>& points,
		      std::vector<double>& param,
		      int ik, int in, 
		      shared_ptr<SplineCurve>& curve);
    
    shared_ptr<SplineCurve> createCurve(std::vector<RevEngPoint*>& points, int degree,
					double tol, int maxiter);

    void extractLinearPoints(std::vector<RevEngPoint*>& points, 
			     std::vector<Point>& rotated, double len,
			     Point& pos, Point& axis, double rad,
			     Point& axis2, bool plane,
			     double tol, double angtol,
			     std::vector<std::pair<double,double> >& dist_ang,
			     std::vector<RevEngPoint*>& linear, bool start,
			     std::vector<RevEngPoint*>& remaining);

    bool extractLinearPoints(std::vector<Point>& points,
			     shared_ptr<Line>& line, Point norm, double tol, 
			     bool start, double& splitpar,
			     std::vector<Point>& linear, 
			     std::vector<Point>& remaining);
    
    void distributePointsToRegions(std::vector<RevEngPoint*>& points,
				   std::vector<shared_ptr<ElementarySurface> >& sfs,
				   shared_ptr<ElementarySurface> curr_sf,
				   double tol, double angtol,
				   std::vector<std::vector<RevEngPoint*> >& sfs_pts,
				   std::vector<RevEngPoint*>& curr_pts,
				   std::vector<RevEngPoint*>& remainin);

    void identifyEndPoints(std::vector<RevEngPoint*> edge_pts, shared_ptr<CurveOnSurface>& sfcv,
			   RevEngPoint*& first_pt, double& t1, RevEngPoint*& last_pt, double& t2);
    
    void setLoopSeq(std::vector<shared_ptr<CurveOnSurface> >& cvs);

    void identifyConGroups(std::vector<RevEngPoint*>& init,
			   std::vector<std::vector<RevEngPoint*> >& groups);
  }
  
} // namespace Go


#endif // _REVENGUTILS_H

