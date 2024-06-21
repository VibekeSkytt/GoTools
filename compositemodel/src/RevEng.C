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

#include "GoTools/compositemodel/RevEng.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngEdge.h"
//#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "sislP.h"
#include <vector>
#include <set>
#include <fstream>
#include <iostream> // @@ debug

using namespace Go;
using std::vector;
using std::pair;
using std::istream;
using std::ostream;

typedef MatrixXD<double, 3> Matrix3D;

#define MAX_COLORS 12
int colors[MAX_COLORS][3] = {
  {255, 0, 0},
  {0, 255, 0},
  {0, 0, 255},
  {255, 255, 0},
  {255, 0, 255},
  {0, 255, 255},
  {128, 255, 0},
  {255, 128, 0},
  {128, 0, 255},
  {255, 0, 128},
  {0, 128, 255},
  {0, 255, 128},
};

#define DEBUG_DIV
#define DEBUG_EDGE
#define DEBUG_BLEND
//#define DEBUG_MONGE
//#define DEBUG_ENHANCE
#define DEBUG_SEG
#define DEBUG
#define DEBUGONE
#define DEBUG_CHECK
#define DEBUG_PLANAR
#define DEBUG_AXIS
#define DEBUG_GROW
#define DEBUG_VALIDATE
#define DEBUG_EDGE
//#define DEBUG_TRIANG
#define DEBUG_TRIM
#define DEBUG_MODEL

//===========================================================================
RevEng::RevEng(shared_ptr<ftPointSet> tri_sf, double mean_edge_len)
  : tri_sf_(tri_sf), mean_edge_len_(mean_edge_len)
//===========================================================================
{
  // Set default parameters
  model_character_ = MEDIUM_ROUGH;
  initParameters();
  max_next_ = std::min(80, tri_sf_->size()/200);
  max_next_ = std::max(2*min_next_, max_next_);
}


//===========================================================================
RevEng::RevEng()
//===========================================================================
{
  // Empty infrastructure for reading stage
  model_character_ = MEDIUM_ROUGH;
  initParameters();
}


//===========================================================================
RevEng::~RevEng()
//===========================================================================
{
}


int compare_x_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[0] < p2->getPoint()[0]);
  // if (p1->getPoint()[0] < p2->getPoint()[0])
  //   return -1;
  // else if (p1->getPoint()[0] > p2->getPoint()[0])
  //   return 1;
  // else
  //   return 0;
}

int compare_y_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[1] < p2->getPoint()[1]);
  // if (p1->getPoint()[1] < p2->getPoint()[1])
  //   return -1;
  // else if (p1->getPoint()[1] > p2->getPoint()[1])
  //   return 1;
  // else
  //   return 0;
}

int compare_z_par(const RevEngPoint* p1, const RevEngPoint* p2)
{
  return (p1->getPoint()[2] < p2->getPoint()[2]);
  // if (p1->getPoint()[2] < p2->getPoint()[2])
  //   return -1;
  // else if (p1->getPoint()[2] > p2->getPoint()[2])
  //   return 1;
  // else
  //   return 0;
}

//===========================================================================
void RevEng::setBoundingBox()
//===========================================================================
{
  int nmbpt = tri_sf_->size();
  vector<RevEngPoint*> all_pts(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      Point xyz2 = Point(xyz[0], xyz[1], xyz[2]);
      all_pts[ki] = pt;
      if (ki == 0)
	bbox_ = BoundingBox(xyz2, xyz2);
      else
	bbox_.addUnionWith(xyz2);
    }

}

//===========================================================================
void RevEng::enhancePoints()
//===========================================================================
{
  setBoundingBox();
#ifdef DEBUG_ENHANCE  
  std::cout << "Bounding box, min: " << bbox_.low() << ", max: " << bbox_.high() << std::endl;
#endif

  // Update parameters based on surface roughness
  updateParameters();
  int nmbpt = tri_sf_->size();

#ifdef DEBUG_TRIANG
  vector<vector<RevEngPoint*> > conn_groups;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->visited())
	continue;
      vector<RevEngPoint*> curr_group;
      pt->fetchConnected(0, nmbpt, curr_group);
      conn_groups.push_back(curr_group);
    }
  
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->unsetVisited();
    }

  std::ofstream oftri("init_ptgroups.g2");
  for (size_t kj=0; kj<conn_groups.size(); ++kj)
    {
      oftri << "400 1 0 0" << std::endl;
      oftri << conn_groups[kj].size() << std::endl;
      for (size_t kr=0; kr<conn_groups[kj].size(); ++kr)
	oftri << conn_groups[kj][kr]->getPoint() << std::endl;
    }
#endif
  
  int writepoints = 0;
  vector<double> tri_ang(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);

      // Compute surface normal from triangulation
      pt->computeTriangNormal(100.0*mean_edge_len_);
      //double avlen = pt->getMeanEdgLen();
      tri_ang[ki] = pt->getTriangAngle();
    }
  std::sort(tri_ang.begin(), tri_ang.end());
#ifdef DEBUG_ENHANCE  
  std::cout << "Triangle angles: " << tri_ang[0] << " " << tri_ang[nmbpt/4];
  std::cout << " " << tri_ang[nmbpt/2] << " " << tri_ang[3*nmbpt/4];
  std::cout << " " << tri_ang[nmbpt-1] << std::endl;
  std::cout << "norm_ang_lim_ : " << norm_ang_lim_ << std::endl;
#endif
  
  double wgt_nmb = 1.0/(double)nmbpt;
  double av_close = 0.0;
  int max_close = 0;
  int min_close = nmbpt;
  vector<double> lambda_3(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->nmbMonge() > 0)
	continue;  // Already enhanced

      // Compute surface normal from triangulation
      pt->computeTriangNormal(100.0*mean_edge_len_);
      if (pt->isOutlier())
	continue;

      if (pt->getNmbNeighbour() == 0)
	{
	  pt->setOutlier();
	  continue;
	}

      //double avlen = pt->getMeanEdgLen();

      //Fetch nearby points
      vector<RevEngPoint*> nearpts;
      double local_len = pt->getMeanEdgLen(10.0*mean_edge_len_);
      double radius = rfac_*(local_len + mean_edge_len_);
      double radius2 = 0.5*radius;
      radius = std::min(radius, 20.0*mean_edge_len_);
      //radius *= 1.5; // TEST 
      //double radius = 0.5*rfac_*(local_len + mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_, max_next_, nearpts);

      Point mincvec(0.0, 0.0, 0.0), maxcvec(0.0, 0.0, 0.0);
      //      Point mincvec2(0.0, 0.0, 0.0), maxcvec2(0.0, 0.0, 0.0);

      av_close += wgt_nmb*(double)nearpts.size();
      max_close = std::max(max_close, (int)nearpts.size());
      min_close = std::min(min_close, (int)nearpts.size());
      
      if (nearpts.size() >= 3)
	{
#ifdef DEBUG_DIV
	  std::ofstream of("nearpts.g2");
	  if (writepoints)
	    {
	      of << "400 1 0 4 255 0 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << pt->getPoint() << std::endl << std::endl;
	      of << "400 1 0 4 0 255 0 255" << std::endl;
	      of << nearpts.size() << std::endl;
	      for (size_t kr=0; kr<nearpts.size(); ++kr)
		of << nearpts[kr]->getPoint() << std::endl;
	    }
#endif
	  
	  // Compute eigenvectors and values of covariance matrix
	  nearpts.push_back(pt);
	  double lambda[3];
	  double eigenvec[3][3];
	  RevEngUtils::principalAnalysis(nearpts, lambda, eigenvec);
	  Point eigen1(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
	  Point eigen2(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
	  Point eigen3(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
	  lambda_3[ki] = lambda[2];
	  Point tnorm = pt->getTriangNormal();
	  if (tnorm.length() < 1.0e-10)
	    {
	      int stop_norm = 1;
	    }
	  else if (eigen3*tnorm <  0.0)
	    {
	      eigen2 *= -1;
	      eigen3 *= -1;
	    }

	  for (size_t kr=0; kr<nearpts.size(); ++kr)
	    {
	      if (pt->pntDist(nearpts[kr]) <= radius2)
		nearpts[kr]->addCovarianceEigen(eigen1, lambda[0], eigen2, lambda[1],
						eigen3, lambda[2]);
	    }
#ifdef DEBUG_DIV      
	  if (writepoints)
	    {
	      for (int ki=0; ki<3; ++ki)
		{
		  Vector3D vec(eigenvec[ki][0], eigenvec[ki][1], eigenvec[ki][2]);
		  of << "410 1 0 4 0 0 0 255" << std::endl;
		  of << "1" << std::endl;
		  Vector3D curr = pt->getPoint();
		  of << curr << " " << curr+0.1*vec << std::endl;
		}
	    }
#endif
	  // Compute normal and curvature using Monge patch
	  // Point normal;//, mincvec, maxcvec;
	  // double minc, maxc;
	  // double currdist, avdist;
	  // RevEngUtils::computeMonge(curr, nearpts, eigen1, eigen3, normal, mincvec, minc,
	  // 				maxcvec, maxc, currdist, avdist);
	  computeMonge(pt, nearpts, eigen1, eigen3, radius2);
	  // Orient vectors with respect to triangulation normal
	  // The normal vectors should be OK. Curvature vectors are not necessarily
	  // consistent with regard to orientation

	  // double minc2, maxc2;
	  // RevEngUtils::TaubinCurvature(curr, nearpts, eigen1, eigen3, mincvec2, minc2,
	  // 				   maxcvec2, maxc2);
      
	  // if (normal*tnorm < 0.0)
	  // 	normal *= -1;
	  // pt->addMongeInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist,
	  // 		       zero_si_);
	  // 	}
	  // of01 << curr << " " << curr+mincvec << std::endl;
	  // //of02 << curr << " " << curr+mincvec2 << std::endl;
	  // of03 << curr << " " << curr+maxcvec << std::endl;
	  //of04 << curr << " " << curr+maxcvec2 << std::endl;

	  int stop_break = 1;
	}
    }

#ifdef DEBUG
  vector<Vector3D> pts1;
  vector<Vector3D> lin1;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (!pt->isOutlier())
	{
	  pts1.push_back(pt->getPoint());
	  vector<ftSamplePoint*> next = pt->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *pt2 = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (!pt2->isOutlier())
		{
		  lin1.push_back(pt->getPoint());
		  lin1.push_back(pt2->getPoint());
		}
	    }
	}
    }
  std::ofstream ofout1("not_outliers.g2");
  ofout1 << "400 1 0 4 0 0 0 255" << std::endl;
  ofout1 << pts1.size() << std::endl;
  for (size_t kj=0; kj<pts1.size(); ++kj)
    ofout1 << pts1[kj] << std::endl;
  ofout1 << "410 1 0 4 55 100 100 255" << std::endl;
  ofout1 << lin1.size()/2 << std::endl;
  for (size_t kj=0; kj<lin1.size(); kj+=2)
    ofout1 << lin1[kj] << " " << lin1[kj+1] << std::endl;
#endif
  std::sort(lambda_3.begin(), lambda_3.end());
#ifdef DEBUG_ENHANCE  
  std::cout << "No close, min: " << min_close << ", max: " << max_close << ", average: " << av_close << std::endl;
  std::cout << "lambda3, min: " << lambda_3[0] << ", max: " << lambda_3[nmbpt-1] << ", medium: " << lambda_3[nmbpt/2] << std::endl;
#endif

  if (false)
    {
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (pt->isOutlier())
	continue;
      double rp[2];
      setRp(pt, rp);
      pt->setRp(rp);
    }
    }

  for (int ki=0; ki<nmbpt; ++ki)
    {
      // Check orientation of curvature
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      if (pt->isOutlier())
	continue;
      Point tnorm = pt->getTriangNormal();
      if (tnorm.length() < 1.0e-10)
	{
	  // Fetch triangulation normal in neighbour
	  vector<ftSamplePoint*> next = pt->getNeighbours();
	  double mindist = std::numeric_limits<double>::max();
	  int ix = -1;
	  for (size_t kr=0; kr<next.size(); ++kr)
	    {
	      double dist = pt->pntDist(next[kr]);
	      RevEngPoint *nextpt = dynamic_cast<RevEngPoint*>(next[kr]);
	      Point nextnorm = nextpt->getTriangNormal();
	      if (dist < mindist && nextnorm.length() > 1.0e-10)
		{
		  mindist = dist;
		  ix = (int)kr;
		}
	    }
	  if (ix >= 0)
	    {
	      RevEngPoint *nextpt = dynamic_cast<RevEngPoint*>(next[ix]);
	      tnorm = nextpt->getTriangNormal();
	      Point PCAnorm = pt->getPCANormal();
	      if (tnorm*PCAnorm < 0.0)
		pt->turnPCA();
	      Point Mongenorm = pt->getMongeNormal();
	      if (tnorm*Mongenorm < 0.0)
		pt->turnMongeNorm();
	    }
	  
	}
    }
  
  //setClassificationParams();
  // of01 << "410 1 0 4 0 0 255 255" << std::endl;
  // of01 << nmbpt << std::endl;
  // of03 << "410 1 0 4 0 255 0 255" << std::endl;
  // of03 << nmbpt << std::endl;


#ifdef DEBUG_ENHANCE  
  std::cout << "Start curvature filter" << std::endl;
#endif
  curvatureFilter();
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->resetPointAssociation();
    }
  
#ifdef DEBUG_ENHANCE  
  std::cout << "Finish curvature filter" << std::endl;
#endif 
  int stop_break = 1;

}

//===========================================================================
void RevEng::computeMonge(RevEngPoint* pt, std::vector<RevEngPoint*>& points,
			  Point& vec1, Point& vec2, double radius)
//===========================================================================
{
  // Transform points to coordinate system given by vec1 (x-axis) and vec2 (y-axis)
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D zaxis(0, 0, 1);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, zaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;
  //rotmat.identity();

  // Perform rotation and sort parameter values and z-value
  int nmbpts = (int)points.size();
  vector<double> par(2*nmbpts);
  vector<double> zval(nmbpts);
  Vector3D curr = pt->getPoint();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      Vector3D dv = points[ki]->getPoint() - curr;
      Vector3D dvrot = rotmat*dv;
      //Vector3D dvrot = mat2*dvrot0;
      par[2*ki] = curr[0] + dvrot[0];
      par[2*ki+1] = curr[1] + dvrot[1];
      zval[ki] = curr[2] + dvrot[2];
    }

  // Approximate z-component by biquadratic Bezier function in x and y
  int order = 3;
  shared_ptr<SplineSurface> mongesf = RevEngUtils::surfApprox(zval, 1, par, order,
							      order, order, order);

  vector<double> coefs2(3*order*order);
  std::vector<double>::iterator cf = mongesf->coefs_begin();
  for (int ka=0; ka<order; ++ka)
    {
      double vpar = mongesf->basis_v().grevilleParameter(ka);
      for (int kb=0; kb<order; ++kb, ++cf)
	{
	  double upar = mongesf->basis_u().grevilleParameter(kb);
	  coefs2[(ka*order+kb)*3] = upar;
	  coefs2[(ka*order+kb)*3+1] = vpar;
	  coefs2[(ka*order+kb)*3+2] = *cf;
	}
    }
  shared_ptr<SplineSurface> tmp(new SplineSurface(order, order, order, order, 
						  mongesf->basis_u().begin(),
						  mongesf->basis_v().begin(), &coefs2[0], 3));
#ifdef DEBUG_DIV
  int writesurface = 0;
  if (writesurface)
    {
      std::ofstream of("approx_sf.g2");
      tmp->writeStandardHeader(of);
      tmp->write(of);
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << 1 << std::endl;
      of << curr << std::endl;
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << nmbpts << std::endl;
      for (int ka=0; ka<nmbpts; ++ka)
	{
	  Point tmppt(par[2*ka], par[2*ka+1], zval[ka]);
	  of << tmppt << std::endl;
	}
    }
#endif
  
  // Compute surface normal 
  double avdist = 0.0;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      Point pos;
      mongesf->point(pos, par[2*ki], par[2*ki+1]);
      avdist += fabs(zval[ki] - pos[0]);
    }
  avdist /= (double)nmbpts;

#ifdef DEBUG_MONGE
  std::ofstream of2("Monge_curvature.g2");
  std::ofstream of3("Monge_curvature2.g2");
#endif
  vector<Point> monge1, monge2, monge3, monge4;
   for (size_t kr=0; kr<points.size(); ++kr)
    {
      if (pt->pntDist(points[kr]) > radius)
	continue;

      Point triang_norm = points[kr]->getTriangNormal();
      vector<Point> der(3);
      mongesf->point(der, par[2*kr], par[2*kr+1], 1);
      Vector3D norm(-der[1][0], -der[2][0], 1.0);
      norm.normalize();

      // Accuracy of approximation
      double currdist = fabs(zval[kr] - der[0][0]);
  
      // Compute principal curvatures in curr
      SISLSurf *sislsf = GoSurf2SISL(*mongesf, false);
      int left1 = 0, left2 = 0;
      int stat = 0;
      double minc, maxc;
      double d1[2], d2[2];
      s2542(sislsf, 0, 0, 0, &par[0], &left1, &left2, &minc, &maxc, d1, d2, &stat);
      Vector3D du(1.0, 0.0, der[1][0]);
      Vector3D dv(0.0, 1.0, der[2][0]);
      Vector3D cvec1 = d1[0]*du + d1[1]*dv;
      Vector3D cvec2 = d2[0]*du + d2[1]*dv;
      if (sislsf) freeSurf(sislsf);

      // Vector3D origin(par[0], par[1], zval[0]);
      // of << "410 1 0 4 0 0 0 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+norm << std::endl;

      // of << "410 1 0 4 0 55 155 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+cvec1 << std::endl;

  
      // of << "410 1 0 4 155 55 0 255" << std::endl;
      // of << "1" << std::endl;
      // of << origin << " " << origin+cvec2 << std::endl;

  
  
      // Transform results to original coordinate system
      Matrix3D mat3, mat4, rotmat2;
      mat4.setToRotation(zaxis, vec2_3);
      mat3.setToRotation(xaxis, vec1_2);
      rotmat2 = mat3*mat4;
      //rotmat2.identity();
      //Vector3D norm0 = mat4*norm;
      Vector3D norm2 = rotmat2*norm;
      Point normal = Point(norm2[0], norm2[1], norm2[2]);
      if (triang_norm.length() > 1.0e-10 && normal*triang_norm < 0.0)
	normal *= -1;
  
      Vector3D cvec3 = rotmat2*cvec1;
      Point mincvec = Point(cvec3[0], cvec3[1], cvec3[2]); 
      Vector3D cvec4 = rotmat2*cvec2;
      Point maxcvec = Point(cvec4[0], cvec4[1], cvec4[2]);
      points[kr]->addMongeInfo(normal, mincvec, minc, maxcvec, maxc, currdist, avdist,
      		       zero_si_);

      // Vector3D xyz = points[kr]->getPoint();
      // Point xyz2 = Point(xyz[0], xyz[1], xyz[2]);
      // Vector3D der2(der[0][0], der[0][1], der[0][2]);
      // Vector3D der3 = rotmat2*der2;
      // Point der4(der3[0], der3[1], der3[2]);  // Not a 3D point!!!
      // monge1.push_back(xyz2);
      // monge1.push_back(xyz2+mincvec);
      // monge2.push_back(xyz2);
      // monge2.push_back(xyz2+maxcvec);
      // monge3.push_back(der4);
      // monge3.push_back(der4+mincvec);
      // monge4.push_back(der4);
      // monge4.push_back(der4+maxcvec);
    }

   // int writeMonge = 0;

   // if (writeMonge)
   //   {
   //     of2 << "410 1 0 4 255 0 0 255" << std::endl;
   //     of2 << monge1.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge1.size(); kr+=2)
   // 	 of2 << monge1[kr] << " " << monge1[kr+1] << std::endl;
   //     of2 << "410 1 0 4 0 0 255 255" << std::endl;
   //     of2 << monge2.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge2.size(); kr+=2)
   // 	 of2 << monge2[kr] << " " << monge2[kr+1] << std::endl;
   
   //     of3 << "400 1 0 4 0 255 0 255" << std::endl;
   //     of3 << monge3.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge3.size(); kr+=2)
   // 	 of3 << monge3[kr] << std::endl;
   //     of3 << "410 1 0 4 255 0 0 255" << std::endl;
   //     of3 << monge3.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge3.size(); kr+=2)
   // 	 of3 << monge3[kr] << " " << monge3[kr+1] << std::endl;
   //     of3 << "410 1 0 4 0 0 255 255" << std::endl;
   //     of3 << monge4.size()/2 << std::endl;
   //     for (size_t kr=0; kr<monge4.size(); kr+=2)
   // 	 of3 << monge4[kr] << " " << monge4[kr+1] << std::endl;
   //   }
  int stop_break = 1;
}
 
//===========================================================================
void RevEng::setRp(RevEngPoint* first, double rp[2])
//===========================================================================
{
  if (first->nmbMonge() == 0)
    return;

  // Fetch associated triangles
  vector<RevEngPoint*> next;
  first->setVisited();
  next.push_back(first);

  vector<vector<int> > tri;
  first->getAttachedTriangles(tri);
  rp[0] = computeRp(first, tri);
  size_t ki=0;
  //for (int kb=1; kb<2; ++kb)
  for (int kb=1; kb<1; ++kb)  // To save time
    {
      size_t nn = next.size();
      for (; ki<nn; ++ki)
	{
	  vector<ftSamplePoint*> curr_next = next[ki]->getNeighbours();
	  for (size_t kj=0; kj<curr_next.size(); ++kj)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(curr_next[kj]);
	      if (pt->visited())
		continue;
	      next.push_back(pt);
	      vector<vector<int> > tri2;
	      pt->getAttachedTriangles(tri2);
	      for (size_t kr=0; kr<tri2.size(); ++kr)
		{
		  auto it = std::find(tri.begin(), tri.end(), tri2[kr]);
		  if (it == tri.end())
		    tri.push_back(tri2[kr]);
		}
	    }
	}
      rp[kb] = computeRp(first, tri);
    }

  for (size_t kr=0; kr<next.size(); ++kr)
    next[kr]->unsetVisited();
}

//===========================================================================
double RevEng::computeRp(RevEngPoint* first, vector<vector<int> >& tri)
//===========================================================================
{
  double rp = 0.0;
  Point norm = first->getMongeNormal();
  Vector3D pnt0 = first->getPoint();
  Point pnt(pnt0[0], pnt0[1], pnt0[2]);
  size_t ntri = tri.size();
  Point cvecmax = first->maxCurvatureVec();
  double kmax = first->maxPrincipalCurvature();
  Point cvecmin = first->minCurvatureVec();
  double kmin = first->minPrincipalCurvature();
  for (size_t ki=0; ki<ntri; ++ki)
    {
      RevEngPoint* pt1 = dynamic_cast<RevEngPoint*>((*tri_sf_)[tri[ki][0]]);
      RevEngPoint* pt2 = dynamic_cast<RevEngPoint*>((*tri_sf_)[tri[ki][1]]);
      RevEngPoint* pt3 = dynamic_cast<RevEngPoint*>((*tri_sf_)[tri[ki][2]]);
      Point tnorm0 = pt1->getTriangNormal() + pt2->getTriangNormal() + pt3->getTriangNormal();
      Vector3D xyz1 = pt1->getPoint();
      Point pos1(xyz1[0], xyz1[1], xyz1[2]);
      Vector3D xyz2 = pt2->getPoint();
      Point pos2(xyz2[0], xyz2[1], xyz2[2]);
      Vector3D xyz3 = pt3->getPoint();
      Point pos3(xyz3[0], xyz3[1], xyz3[2]);
      Point btri = (pos1 + pos2 + pos3)/3.0;
      Point vec1 = pos2 - pos1;
      Point vec2 = pos3 - pos1;
      Point tnorm = vec1.cross(vec2);
      if (tnorm0*tnorm < 0.0)
	tnorm *= -1;
      tnorm.normalize_checked();
      Point vec3 = btri - pnt;
      vec3.normalize_checked();
      double fac1 = vec3*cvecmax;
      double fac2 = vec3*cvecmin;
      double curv = fac1*kmax + fac2*kmin;
      double div = 1.0 - pnt.dist2(btri)*curv*curv;
      double tmp = (div > 1.0e-12) ? 1.0 - ((tnorm*norm)/sqrt(div)) : 0.0;
      rp += fabs(tmp);
    }
  rp /= (double)ntri;
  return rp;
}

//===========================================================================
void RevEng::curvatureFilter()
//===========================================================================
{
  int nmbpt = tri_sf_->size();
  double radius_fac = 0.2; //0.7;
  vector<vector<RevEngPoint*> > nearpts(nmbpt);
  bool smoothcurv = true; //false;
  if (smoothcurv)
    {
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      
      // Fetch nearby points
      double local_len = pt->getMeanEdgLen();
      double radius = 0.5*radius_fac*rfac_*(local_len + mean_edge_len_);
      radius = std::min(radius, 20.0*mean_edge_len_);
      pt->fetchClosePoints2(radius, min_next_/2, max_next_/2, nearpts[ki]);
      if (nearpts[ki].size() == 0)
	continue;
      vector<double> H0(nearpts[ki].size()+1);
      vector<double> K0(nearpts[ki].size()+1);
      H0[0] = pt->meanCurvature0();
      K0[0] = pt->GaussCurvature0();
      for (size_t kj=0; kj<nearpts[ki].size(); ++kj)
	{
	  H0[kj+1] = nearpts[ki][kj]->meanCurvature0();
	  K0[kj+1] = nearpts[ki][kj]->GaussCurvature0();
	}
      std::sort(H0.begin(), H0.end());
      std::sort(K0.begin(), K0.end());
      pt->setMeanCurvature(0.5*H0[H0.size()/2] + H0[(H0.size()+1)/2]);
      pt->setGaussCurvature(0.5*K0[K0.size()/2] + K0[(K0.size()+1)/2]);
    }
  
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->updateCurvature();
    }

  int nmbsmooth = (model_character_ <= MEDIUM_ROUGH) ? 2 : 5;
  for (int ka=0; ka<nmbsmooth; ++ka)
    {
      for (int ki=0; ki<nmbpt; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      
	  // Fetch nearby points
	  double Hmean = pt->meanCurvature0();
	  double Kmean = pt->GaussCurvature0();
	  for (size_t kj=0; kj<nearpts[ki].size(); ++kj)
	    {
	      Hmean += nearpts[ki][kj]->meanCurvature0();
	      Kmean += nearpts[ki][kj]->GaussCurvature0();
	    }
	  Hmean /= (double)(nearpts[ki].size()+1);
	  Kmean /= (double)(nearpts[ki].size()+1);
	  pt->setMeanCurvature(Hmean);
	  pt->setGaussCurvature(Kmean);
	}
      for (int ki=0; ki<nmbpt; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
	  pt->updateCurvature();
	}
   }
    }
}

//===========================================================================
void RevEng::setClassificationParams()
//===========================================================================
{
#ifdef DEBUG_ENHANCE  
  int class_type = getClassificationType();
  std::cout << "Classification type: " << class_type << std::endl;
#endif
  // std::cout << "New classification type: " << std::endl;
  // std::cin >> class_type;
  // if (class_type == 1 || class_type == 3)
  //   setClassificationType(class_type);

  // if (class_type == CLASSIFICATION_CURVATURE)
  //   {
  //     double zero_H = getMeanCurvatureZero();
  //     std::cout << "Mean curvature zero limit: " << zero_H << ", give limit: " << std::endl;
  //     std::cin >> zero_H;
  //     setMeanCurvatureZero(zero_H);
      
  //     double zero_K = getGaussCurvatureZero();
  //     std::cout << "Gauss curvature zero limit: " << zero_K << ", give limit: " << std::endl;
  //     std::cin >> zero_K;
  //     setGaussCurvatureZero(zero_K);
  //   }
  // else if (class_type == CLASSIFICATION_POINTASSOCIATION)
  //   {
  //     std::cout << "ffac = " << ffac_ << ", give ffac: " << std::endl;
  //     std::cin >> ffac_;
  //     std::cout << "sfac = " << sfac_ << ", give sfac: " << std::endl;
  //     std::cin >> sfac_;
  //   }

  int nmbpt = tri_sf_->size();
#ifdef DEBUG_ENHANCE
  std::ofstream of01("minc1.g2");
  std::ofstream of03("maxc1.g2");
  of01 << "410 1 0 4 200 50 0 255" << std::endl;
  of01 << nmbpt << std::endl;
  of03 << "410 1 0 4 0 50 200 255" << std::endl;
  of03 << nmbpt << std::endl;

  std::ofstream of("triangnorm.g2");
  std::ofstream ofM("Mongenorm.g2");
  std::ofstream ofP("PCAnorm.g2");
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << nmbpt << std::endl;
  ofM << "410 1 0 4 255 0 0 255" << std::endl;
  ofM << nmbpt << std::endl;
  ofP << "410 1 0 4 0 200 55 255" << std::endl;
  ofP << nmbpt << std::endl;
#endif
  vector<Vector3D> triangplane;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Point minc = pt->minCurvatureVec();
      Point maxc = pt->maxCurvatureVec();
      Vector3D xyz = pt->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = pt->getTriangNormal();
      double ang = pt->getTriangAngle();
      Point Mnorm = pt->getMongeNormal();
      Point Pnorm = pt->getPCANormal();

#ifdef DEBUG_ENHANCE
      double avlen = pt->getMeanEdgLen();
      double fac = (pt->nmbMonge() == 0) ? 0.0 : 5.0;
      of01 << xyz2 << " " << xyz2+fac*avlen*minc << std::endl;
      of03 << xyz2 << " " << xyz2+fac*avlen*maxc << std::endl;

      of << xyz2 << " " << xyz2+5.0*avlen*norm << std::endl;
      ofM << xyz2 << " " << xyz2+5.0*avlen*Mnorm << std::endl;

      ofP << xyz2 << " " << xyz2+5.0*avlen*Pnorm << std::endl;
#endif
      if (pt->isOutlier())
	continue;
      
      if (ang <= norm_plane_lim_)
	triangplane.push_back(xyz);
    }

#ifdef DEBUG_ENHANCE
  std::ofstream of4("triangplane.g2");
  of4 << "400 1 0 4 200 0 200 255" << std::endl;
  of4 << triangplane.size() << std::endl;
  for (size_t kj=0; kj<triangplane.size(); ++kj)
    of4 << triangplane[kj] << std::endl;
#endif
  
  // // Surface variation
  // std::sort(sfvar.begin(), sfvar.end());
  // double varmed = 0.5*(sfvar[nmbpts/2]+sfvar[(nmbpts+1)/2]);
  // double varQ1 = sfvar[nmbpts/4];
  // double varQ2 = sfvar[3*nmbpts/4];

  // double stepmean = (sfvar[sfvar.size() - 1] - sfvar[0])/(double)(nmbpts-1);
  // double stepfac = 20.0; //50.0; //100.0;
  // double steplim1 = stepfac*stepmean;
  // vector<int> largestep;
  // for (int ki=0; ki<nmbpts-1; ++ki)
  //   {
  //     double diff = sfvar[ki+1] - sfvar[ki];
  //     if (diff > steplim1)
  // 	largestep.push_back(ki);
  //   }

  // // Curvedness
  // std::sort(curvedness.begin(), curvedness.end());
  // double stepmean2 = (curvedness[curvedness.size() - 1] - curvedness[0])/(double)(nmbpts-1);
  // double steplim2 = stepfac*stepmean2;
  // double curvedmed = 0.5*(curvedness[nmbpts/2] + curvedness[nmbpts/2+1]);
  // double curvedQ1 = curvedness[nmbpts/4];
  // double curvedQ2 = curvedness[3*nmbpts/4];
  // vector<int> largecurved;
  // for (int ki=0; ki<nmbpts-1; ++ki)
  //   {
  //     double diff = curvedness[ki+1] - curvedness[ki];
  //     if (diff > steplim2)
  // 	largecurved.push_back(ki);
  //   }

  // Testing
  // std::ofstream of1("sfvar_1.g2");
  // of1 << "400 1 0 4 155 0 100 255" << std::endl;
  // of1 << nmbpts-largestep[0]-1 << std::endl;
  // std::ofstream of2("curvedness_1.g2");
  // of2 << "400 1 0 4 0 155 100 255" << std::endl;
  // of2 << nmbpts-largecurved[0]-1 << std::endl;
//   vector<Vector3D> pts_v, pts_c;
//   for (int ka=0; ka<nmbpts; ++ka)
//     {
//       RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
//       double var = pt->getSurfaceVariation();
//       double curved = pt->getCurvedness();
//       Vector3D xyz = pt->getPoint();
//       // if (var > sfvar[largestep[0]])
//       // 	of1 << xyz << std::endl;
//       // if (curved > curvedness[largecurved[0]])
//       // 	of2 << xyz << std::endl;
//       if (var > varh)
// 	pts_v.push_back(xyz);
//       if (curved > curvh)
// 	pts_c.push_back(xyz);
//     }
//   std::ofstream of3("sfvar.g2");
//   of3 << "400 1 0 4 155 0 100 255" << std::endl;
//   of3 << pts_v.size() << std::endl;
//   for (size_t kr=0; kr<pts_v.size(); ++kr)
//     of3 << pts_v[kr] << std::endl;
//   std::ofstream of5("curvedness.g2");
//   of5 << "400 1 0 4 0 155 100 255" << std::endl;
//   of5 << pts_c.size() << std::endl;
//   for (size_t kr=0; kr<pts_c.size(); ++kr)
//     of5 << pts_c[kr] << std::endl;

// int stop_break = 1;
}

//===========================================================================
void RevEng::setEdgeClassificationParams()
//===========================================================================
{
  int edge_class_type = getEdgeClassificationType();
#ifdef DEBUG_ENHANCE  
  std::cout << "Edge classification type: " << edge_class_type << std::endl;
#endif
  // std::cout << "New classification type: " << std::endl;
  // std::cin >> edge_class_type;
  if (edge_class_type == CURVATURE_EDGE || edge_class_type == PCATYPE_EDGE ||
      edge_class_type == CNESS_EDGE || edge_class_type == RPFAC_EDGE)
    {
      setEdgeClassificationType(edge_class_type);
      if (edge_class_type == CURVATURE_EDGE)
	{
	  //double cfac = getCfac();
	  // std::cout << "Edge classification factor with curvature: " << cfac << std::endl;
	  // std::cout << "New classification factor: " << std::endl;
	  // std::cin >> cfac;
	  // if (cfac >= 5.0 && cfac <= 10.0)
	  //   setCfac(cfac);
	  
	  vector<Vector3D> curvaturecorners;
	  int nmbpts = tri_sf_->size();
	  for (int ka=0; ka<nmbpts; ++ka)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
	      double avlen = pt->getMeanEdgLen();
	      double maxpc = std::max(fabs(pt->maxPrincipalCurvature()),
				      fabs(pt->minPrincipalCurvature()));
	      Vector3D xyz = pt->getPoint();
	      double crvrad = 1.0/maxpc; //fabs(maxpc);
	      if (crvrad < cfac_*avlen)
		curvaturecorners.push_back(xyz);
	    }
#ifdef DEBUG_ENHANCE
	  std::ofstream of3("curvaturecorners0.g2");
	  of3 << "400 1 0 4 10 10 10 255" << std::endl;
	  of3 << curvaturecorners.size() << std::endl;
	  for (size_t kj=0; kj<curvaturecorners.size(); ++kj)
	    of3 << curvaturecorners[kj] << std::endl;
#endif
	}
      else if (edge_class_type == PCATYPE_EDGE)
	{
	  double pca_lim = getPCAlim();
	  // std::cout << "Edge classification factor with PCA: " << pca_lim << std::endl;
	  // std::cout << "New classification factor: " << std::endl;
	  // std::cin >> pca_lim;
	  // setPCAlim(pca_lim);
	  
	  vector<Vector3D> pts_v;
	  int nmbpts = tri_sf_->size();
	  for (int ka=0; ka<nmbpts; ++ka)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
	      double var = pt->getSurfaceVariation();
	      Vector3D xyz = pt->getPoint();
	      if (var > pca_lim)
		pts_v.push_back(xyz);
	    }
#ifdef DEBUG_ENHANCE
	  std::ofstream of3("sfvar.g2");
	  of3 << "400 1 0 4 155 0 100 255" << std::endl;
	  of3 << pts_v.size() << std::endl;
	  for (size_t kr=0; kr<pts_v.size(); ++kr)
	    of3 << pts_v[kr] << std::endl;
#endif
	}
      else if (edge_class_type == CNESS_EDGE)
	{
	  double cness_lim = getCnesslim();
	  // std::cout << "Edge classification factor with cness: " << cness_lim << std::endl;
	  // std::cout << "New classification factor: " << std::endl;
	  // std::cin >> cness_lim;
	  // setCnesslim(cness_lim);
	  
	  vector<Vector3D> pts_c;
	  int nmbpts = tri_sf_->size();
	  for (int ka=0; ka<nmbpts; ++ka)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
	      double curved = pt->getCurvedness();
	      Vector3D xyz = pt->getPoint();
	      if (curved > cness_lim)
		pts_c.push_back(xyz);
	    }
#ifdef DEBUG_ENHANCE
	  std::ofstream of5("curvedness.g2");
	  of5 << "400 1 0 4 0 155 100 255" << std::endl;
	  of5 << pts_c.size() << std::endl;
	  for (size_t kr=0; kr<pts_c.size(); ++kr)
	    of5 << pts_c[kr] << std::endl;
#endif
	}
      // else if (edge_class_type == RPFAC_EDGE)
      // 	{
	  //double RP_fac = getRPfac();
	  // std::cout << "Edge classification factor with RP: " << RP_fac << std::endl;
	  // std::cout << "New classification factor: " << std::endl;
	  // std::cin >> RP_fac;
	  // setRPfac(RP_fac);
	// }
    }
}

//===========================================================================
double RevEng::getPCAlim()
//===========================================================================
{
  if (pca_lim_ < 0.0)
    {
      int nmbpts = tri_sf_->size();
      double varh = 0.0;
      double hfac = 14.0/(15.0*nmbpts); //6.0/(nmbpts*7.0);
      double sfvar;
      for (int ki=0; ki<nmbpts; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
	  sfvar = pt->getSurfaceVariation();
	  varh += hfac*sfvar;
	}
      pca_lim_ = varh;
     }

  return pca_lim_;
}

//===========================================================================
double RevEng::getCnesslim()
//===========================================================================
{
  if (cness_lim_ < 0.0)
    {
      int nmbpts = tri_sf_->size();
      //double varh = 0.0;
      double hfac = 14.0/(15.0*nmbpts); //6.0/(nmbpts*7.0);
      double curvedness;
      double curvh;
      for (int ki=0; ki<nmbpts; ++ki)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
	  curvedness = pt->getCurvedness();
	  curvh += hfac*curvedness;
	}
      cness_lim_ = curvh;
     }

  return cness_lim_;
}

//===========================================================================
void RevEng::edgeClassification()
//===========================================================================
{
  //double zero_si = getShapeIndexZero();
  // std::cout << "Shape index zero limit: " << zero_si << ", give limit: " << std::endl;
  // std::cin >> zero_si;
  // setShapeIndexZero(zero_si);
  
  int nmbpts = tri_sf_->size();
  vector<Vector3D> triangcorners;
  vector<Vector3D> curvaturecorners;
  vector<Vector3D> smoothnesscorners;
  vector<Vector3D> PCAcorners;
  vector<Vector3D> Rpcorners;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      
      if (pt->getTriangAngle() > norm_ang_lim_)
	triangcorners.push_back(xyz);
      
       // Curvature edge classification
      double avlen = pt->getMeanEdgLen();
      double maxpc = std::max(fabs(pt->maxPrincipalCurvature()),
			      fabs(pt->minPrincipalCurvature()));
      double crvrad = 1.0/maxpc; 
      int c1_edge = (crvrad < cfac_*avlen) ? C1_EDGE : C1_NOT_EDGE;
      if (c1_edge == C1_EDGE)
	curvaturecorners.push_back(xyz);

      // Edge classification with curvedness
      double curved = pt->getCurvedness();
      int c2_edge = (curved > cness_lim_) ? C2_EDGE : C2_NOT_EDGE;
      if (c2_edge == C2_EDGE)
	smoothnesscorners.push_back(xyz);

      // Edge classification with surface variation
      double var = pt->getSurfaceVariation();
      int pca_edge = (var > pca_lim_) ? PCA_EDGE : PCA_NOT_EDGE;
      if (pca_edge == PCA_EDGE)
	PCAcorners.push_back(xyz);

      double rp = pt->getRp(rpix_);
      int rp_edge =  (rp >= rpfac_) ? RP_EDGE : RP_NOT_EDGE;
      if (rp_edge == RP_EDGE)
	Rpcorners.push_back(xyz);

      // Store classification in point
      pt->setEdgeClassification(c1_edge, c2_edge, pca_edge, rp_edge);
    }
  
  // Specify/clean edge classification
  bool closeedge = false;
  int nmbedge = 2;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      if (pt->isEdge(edge_class_type_))
	{
	  if (pt->isolatedEdge(edge_class_type_, nmbedge, closeedge))
	    pt->setEdgeUndef();

	  else 
	    // If the angular difference between triangle normals is less
	    // then the limit, classify the point as CLOSE_EDGE.
	    pt->adjustWithTriangNorm(norm_ang_lim_);
	}
    }
  
  vector<Vector3D> edgepts;
  for (int ka=0; ka<nmbpts; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      if (pt->isEdge(edge_class_type_))
	edgepts.push_back(pt->getPoint());
    }


#ifdef DEBUG_ENHANCE
  std::ofstream of2("triangcorners.g2");
  of2 << "400 1 0 4 0 0 0 255" << std::endl;
  of2 << triangcorners.size() << std::endl;
  for (size_t kj=0; kj<triangcorners.size(); ++kj)
    of2 << triangcorners[kj] << std::endl;

  std::ofstream of3("curvaturecorners.g2");
  of3 << "400 1 0 4 10 10 10 255" << std::endl;
  of3 << curvaturecorners.size() << std::endl;
  for (size_t kj=0; kj<curvaturecorners.size(); ++kj)
    of3 << curvaturecorners[kj] << std::endl;

  std::ofstream of5("smoothnesscorners.g2");
  of5 << "400 1 0 4 10 10 10 255" << std::endl;
  of5 << smoothnesscorners.size() << std::endl;
  for (size_t kj=0; kj<smoothnesscorners.size(); ++kj)
    of5 << smoothnesscorners[kj] << std::endl;

  std::ofstream of6("PCAcorners.g2");
  of6 << "400 1 0 4 10 10 10 255" << std::endl;
  of6 << PCAcorners.size() << std::endl;
  for (size_t kj=0; kj<PCAcorners.size(); ++kj)
    of6 << PCAcorners[kj] << std::endl;

  std::ofstream of7("Rpcorners.g2");
  of7 << "400 1 0 4 10 10 10 255" << std::endl;
  of7 << Rpcorners.size() << std::endl;
  for (size_t kj=0; kj<Rpcorners.size(); ++kj)
    of7 << Rpcorners[kj] << std::endl;
#endif

#ifdef DEBUG_EDGE
   if (edgepts.size() > 0)
    {
      std::ofstream ofedg("edgepts.g2");
      ofedg << "400 1 0 4 255 0 0 255" << std::endl;
      ofedg << edgepts.size() << std::endl;
      for (size_t kr=0; kr<edgepts.size();++kr)
	ofedg << edgepts[kr] << std::endl;
    }
#endif
   int stop_break = 1;
 }

//===========================================================================
void RevEng::classifyPoints()
//===========================================================================
{
  // Fetch relevant values for all points

  vector<vector<Vector3D> > class_pts(9);
  vector<vector<Vector3D> > class_shape(10);
  vector<vector<Vector3D> > class_pfs(10);
  int nmbpts = tri_sf_->size();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();

      // Curvature surface classification
      int ctype = C1_UNDEF;
      double gausscurv = pt->GaussCurvature();
      double meancurv = pt->meanCurvature();
      if (meancurv < -zero_H_)
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_PEAK;
	      class_pts[0].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_SRIDGE;
	      class_pts[2].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_RIDGE;
	      class_pts[1].push_back(xyz);
	    }
	}
      else if (meancurv > zero_H_)
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_PIT;
	      class_pts[6].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_SVALLEY;
	      class_pts[8].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_VALLEY;
	      class_pts[7].push_back(xyz);
	    }
	}
      else
	{
	  if (gausscurv > zero_K_)
	    {
	      ctype = C1_NONE;
	      class_pts[3].push_back(xyz);
	    }
	  else if (gausscurv < -zero_K_)
	    {
	      ctype = C1_MINSURF;
	      class_pts[5].push_back(xyz);
	    }
	  else
	    {
	      ctype = C1_FLAT;
	      class_pts[4].push_back(xyz);
	    }
	}

     // Surface classification with shape index
      int si_type = SI_UNDEF;
      double shapeindex = pt->getShapeIndex();
      if (shapeindex < -0.875)
	{
	  si_type = SI_SCUP;
	  class_shape[0].push_back(xyz);
	}
      else  if (shapeindex < -0.625)
	{
	  si_type = SI_TRO;
	  class_shape[1].push_back(xyz);
	}
      else  if (shapeindex < -0.375)
	{
	  si_type = SI_RUT;
	  class_shape[2].push_back(xyz);
	}
       else  if (shapeindex < -0.125)
	 {
	   si_type = SI_SRUT;
	   class_shape[3].push_back(xyz);
	 }
       else  if (shapeindex < 0.125)
	 {
	   si_type = SI_SAD;
	   class_shape[4].push_back(xyz);
	 }
       else  if (shapeindex < 0.375)
	 {
	   si_type = SI_SRID;
	   class_shape[5].push_back(xyz);
	 }
      else  if (shapeindex < 0.625)
	{
	  si_type = SI_RID;
	  class_shape[6].push_back(xyz);
	}
      else  if (shapeindex < 0.875)
	{
	  si_type = SI_DOM;
	  class_shape[7].push_back(xyz);
	}
      else  if (shapeindex <= 1.0)
	{
	  si_type = SI_SCAP;
	  class_shape[8].push_back(xyz);
	}
      else
	{
	  si_type = SI_PLANE;
	  class_shape[9].push_back(xyz);
	}

      int ps_type = PS_UNDEF;
      double fpa = pt->getfpa();
      double spa = pt->getspa();
      if (fpa <= ffac_)
	{
	  ps_type = PS_PLANE;
	  class_pfs[0].push_back(xyz);
	}
      else if (spa <= -1+sfac_)
	{
	  ps_type = PS_UC;
	  class_pfs[1].push_back(xyz);
	}
      else if (spa < -0.5-sfac_)
	{
	  ps_type = PS_EC;
	  class_pfs[2].push_back(xyz);
	}
      else if (spa <= -0.5+sfac_)
	{
	  ps_type = PS_PC;
	  class_pfs[3].push_back(xyz);
	}
      else if (spa < -sfac_)
	{
	  ps_type = PS_HC;
	  class_pfs[4].push_back(xyz);
	}
      else if (spa <= sfac_)
	{
	  ps_type = PS_HS;
	  class_pfs[5].push_back(xyz);
	}
      else if (spa < 0.5-sfac_)
	{
	  ps_type = PS_HX; 
	  class_pfs[6].push_back(xyz);
	}
      else if (spa <= 0.5+sfac_)
	{
	  ps_type = PS_PX;
	  class_pfs[7].push_back(xyz);
	}
      else if (spa < 1-sfac_)
	{
	  ps_type = PS_EX;
	  class_pfs[8].push_back(xyz);
	}
      else
	{
	  ps_type = PS_UX;
	  class_pfs[9].push_back(xyz);
	}

      // Store classification in point
      pt->setClassification(ctype, si_type, ps_type);
   }

#ifdef DEBUG_SEG
  std::ofstream of("curvature_segments.g2");
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      {
	of << "400 1 0 4 ";
	for (int kc=0; kc<3; ++kc)
	  of << colors[3*ka+kb][kc] << " ";
	of << "255" << std::endl;
	of << class_pts[3*ka+kb].size() << std::endl;
	for (size_t kr=0; kr<class_pts[3*ka+kb].size(); ++kr)
	  of << class_pts[3*ka+kb][kr] << std::endl;
      }


  // Shape index
  std::ofstream of2("shapeindex_segments.g2");
  for (int ka=0; ka<10; ++ka)
    {
      of2 << "400 1 0 4 ";
      for (int kc=0; kc<3; ++kc)
	of2 << colors[ka][kc] << " ";
      of2 << "255" << std::endl;
      of2 << class_shape[ka].size() << std::endl;
      for (size_t kr=0; kr<class_shape[ka].size(); ++kr)
	of2 << class_shape[ka][kr] << std::endl;
      }
  
   // Point classification association
  std::ofstream of4("ptsclass_segments.g2");
  for (int ka=0; ka<10; ++ka)
    {
      of4 << "400 1 0 4 ";
      for (int kc=0; kc<3; ++kc)
	of4 << colors[ka][kc] << " ";
      of4 << "255" << std::endl;
      of4 << class_pfs[ka].size() << std::endl;
      for (size_t kr=0; kr<class_pfs[ka].size(); ++kr)
	of4 << class_pfs[ka][kr] << std::endl;
      }
#endif  
  int stop_break = 1;
}

struct
{
  bool operator()(shared_ptr<RevEngRegion> a, shared_ptr<RevEngRegion> b)
  {
    return (a->numPoints() > b->numPoints());
    //   return 1;
    // else if (a->numPoints() == b->numPoints())
    //   return 2;
    // else
    //   return 3;
  }
} sort_region;

//===========================================================================
void RevEng::setApproxTolerance()
//===========================================================================
{
  double eps = getInitApproxTol();
#ifdef DEBUG
  std::cout << "Approx tol: " << eps << std::endl;
#endif
  // std::cout << "New tolerance: " << std::endl;
  // std::cin >> eps;
  setApproxTol(eps);
}

//===========================================================================
double RevEng::getInitApproxTol()
//===========================================================================
{
  int nmbpts = tri_sf_->size();
  vector<double> pointdist;
  double maxdist = 0.0;
  double avptdist = 0.0;
  double minptdist = std::numeric_limits<double>::max();
  double maxptdist = 0.0;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (!pt->isEdge(edge_class_type_))
	{
	  double ptdist = pt->getPointDistance();
	  pointdist.push_back(ptdist);
	  minptdist = std::min(minptdist, ptdist);
	  maxptdist = std::max(maxptdist, ptdist);
	  avptdist += ptdist;
	}
    }
  if (pointdist.size() > 0)
    avptdist /= (double)pointdist.size();
  
  std::sort(pointdist.begin(), pointdist.end());
  
  std::sort(pointdist.begin(), pointdist.end());
  double dlim = 0.93; //0.75;
  int dix = (int)(dlim*(double)pointdist.size());
  double eps = pointdist[dix];
  if (model_character_ == MEDIUM_ROUGH)
    eps *= 1.5;
  else if (model_character_ == ROUGH)
    eps *= 2.0;

  // Just to test
  //eps = 0.1;

#ifdef DEBUG_DIV
  std::cout << "Maxptdist: " << maxdist << ", avdist: " << avptdist;
  std::cout << ", medptdist: " << pointdist[pointdist.size()/2];
  std::cout << ", eps: " << eps << std::endl;
  
  std::ofstream ofd("ptdist.g2");
#endif
  double ptd1 = minptdist;
  double ptd2 = maxptdist;
  int nmbd = 12;
  double pdel = (ptd2 - ptd1)/(double)nmbd;
  vector<vector<Vector3D> > ptrs(nmbd);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      Vector3D xyz = pt->getPoint();
      double dd = pt->getPointDistance();
      int ix = (int)((dd-ptd1)/pdel);
      ix = std::min(ix, nmbd-1);
      ptrs[ix].push_back(xyz);
    }

#ifdef DEBUG_DIV
  for (int ki=0; ki<nmbd; ++ki)
    {
      if (ptrs[ki].size() > 0)
	{
	  ofd << "400 1 0 4 " << colors[ki][0] << " " << colors[ki][1] << " ";
	  ofd << colors[ki][2] << " 255"  << std::endl;
	  ofd << ptrs[ki].size();
	  for (size_t kr=0; kr<ptrs[ki].size(); ++kr)
	    ofd << ptrs[ki][kr] << std::endl;
	}
    }
#endif

  return eps;
 }

//===========================================================================
void RevEng::segmentIntoRegions()
//===========================================================================
{
  // Collect continous regions
  //int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  int nmbpts = tri_sf_->size();
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      if (pt->hasRegion())
	continue;
      if (pt->closeEdge(edge_class_type_))
	continue;

      if (pt->nmbSameClassification(classification_type_) == 0)
	single_points_.push_back(pt);
      else
	{
	  shared_ptr<RevEngRegion> region(new RevEngRegion(classification_type_,
							   edge_class_type_));
	  pt->setRegion(region.get());
	  region->collect(pt);
	  if (region->numPoints() == 1)
	    {
	      RevEngPoint *pt_single = region->getPoint(0);
	      pt_single->unsetRegion();
	      single_points_.push_back(pt_single);
	    }
	  else
	    regions_.push_back(region);
	}
    }

  std::sort(regions_.begin(), regions_.end(), sort_region);
#ifdef DEBUG
  if (regions_.size() > 0)
    {
      std::cout << "Regions 1" << std::endl;
      std::ofstream of("regions1.g2");
      std::ofstream ofm("mid_regions1.g2");
      std::ofstream ofs("small_regions1.g2");
      writeRegionStage(of, ofm, ofs);
    }
#endif

  // Sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

  min_point_region_ = setSmallRegionNumber();
#ifdef DEBUG
  std::cout << "Min point region: " << min_point_region_ << std::endl;
#endif
#ifdef DEBUG_PLANAR
  std::ofstream ofpc("cand_planar.g2");
#endif  
  double lim_cone = 0.1*M_PI;
  size_t numreg = regions_.size();
  for (size_t ki=0; ki<numreg; ++ki)
    {
      double avH, avK, MAH, MAK;
      regions_[ki]->getAvCurvatureInfo(avH, avK, MAH, MAK);
      if (regions_[ki]->numPoints() > min_point_region_ &&
	  (regions_[ki]->planartype() || MAH <= zero_H_) && 
	  (!regions_[ki]->feasiblePlane(zero_H_, zero_K_)))
	{
#ifdef DEBUG_PLANAR
	  regions_[ki]->writeRegionPoints(ofpc);
#endif
	  vector<vector<RevEngPoint*> > other_groups;
	  vector<RevEngPoint*> single;
	  regions_[ki]->splitPlanar(lim_cone, min_point_region_/2, other_groups, single);
	  if (other_groups.size() > 0)
	    {
	      vector<HedgeSurface*> dummy;
	      surfaceExtractOutput((int)ki, other_groups, dummy);
	    }
	  if (single.size())
	    single_points_.insert(single_points_.end(), single.begin(), single.end());

	  int stop_cand_plane = 1;
	}
    }
  
  // Set adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG_PLANAR
  std::ofstream ofp("planar_reg.g2");
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->feasiblePlane(zero_H_, zero_K_))
	regions_[ki]->writeRegionPoints(ofp);
    }
#endif
  // Integrate single points when appropriate
  vector<RevEngPoint*> remaining_single;
  for (int ka=0; ka<(int)single_points_.size(); ++ka)
    {
      bool merged = single_points_[ka]->mergeWithAdjacent(mean_edge_len_);
      if (!merged)
	{
	  single_points_[ka]->setOutlier();
	  remaining_single.push_back(single_points_[ka]);
	}
    }
  std::swap(single_points_, remaining_single);

  // Merge adjacent planar regions
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->feasiblePlane(zero_H_, zero_K_))
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  bool merged = regions_[ki]->mergePlanarReg(zero_H_, zero_K_, 
						     approx_tol_, mainaxis_,
						     grown_regions);
	  if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	    updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	  regions_[ki]->updateRegionAdjacency();
	}
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
#ifdef DEBUG
  checkConsistence("Regions1_2");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 1_2" << std::endl;
      std::ofstream of("regions1_2.g2");
      std::ofstream ofm("mid_regions1_2.g2");
      std::ofstream ofs("small_regions1_2.g2");
      writeRegionStage(of, ofm, ofs);
     }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // Simplify regions structure
#ifdef DEBUG
  std::cout << "Number of regions pre integrate: " << regions_.size() << std::endl;
#endif
  int num_reg = (int)regions_.size();
  for (int ka=num_reg-1; ka>=0; --ka)
    {
      bool to_be_removed = regions_[ka]->integrateInAdjacent(mean_edge_len_,
							     min_next_, max_next_,
							     approx_tol_, 0.5,
							     max_nmb_outlier_);
      if (to_be_removed)
	{
	  if (ka < num_reg-1)
	    std::swap(regions_[ka], regions_[num_reg-1]);
	  num_reg--;
	}
    }
  if (num_reg < (int)regions_.size())
    regions_.erase(regions_.begin()+num_reg, regions_.end());
#ifdef DEBUG
  std::cout << "Number of regions post integrate: " << regions_.size() << std::endl;
#endif  
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
#ifdef DEBUG
  checkConsistence("Regions2");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 2" << std::endl;
      std::ofstream of("regions2.g2");
      std::ofstream ofm("mid_regions2.g2");
      std::ofstream ofs("small_regions2.g2");
      writeRegionStage(of, ofm, ofs);
     }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // Update minimum number of points in region for surface generation.
  // First sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

  min_point_region_ = setSmallRegionNumber();
#ifdef DEBUG
  std::cout << "Min point region (2): " << min_point_region_ << std::endl;
#endif
}


//===========================================================================
void RevEng::initialSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // Create surfaces in simple regions (planes and cylinders) and extract
  // deviant points
  int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  double angtol = 5.0*anglim_;
  //int regsize = (int)regions_.size();
  std::sort(regions_.begin(), regions_.end(), sort_region);
#ifdef DEBUG
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      RevEngRegion *first = regions_[ki]->getPoint(0)->region();
      int num = regions_[ki]->numPoints();
      for (int ka=1; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != first)
	  std::cout << "Inconsistent region pointers, pre initPlaneCyl: " << ki << " " << ka << std::endl;
    }
#endif
  for (int kr=0; kr<(int)regions_.size(); ++kr)
    {
#ifdef DEBUG
      std::ofstream of0("init_reg.g2");
      regions_[kr]->writeRegionInfo(of0);
#endif
      
      if (regions_[kr]->numPoints() < min_point_region_)
	continue;

      vector<vector<RevEngPoint*> > out_groups;
      vector<RevEngPoint*> single;
      vector<shared_ptr<HedgeSurface> > sfs;
      vector<HedgeSurface*> prev_sfs;
      bool repeat = false;
      regions_[kr]->initPlaneCyl(min_point_in, min_point_region_,
				 approx_tol_, angtol, mainaxis_,
				 zero_H_, zero_K_, sfs, out_groups, single, repeat);
      if (single.size() > 0)
	single_points_.insert(single_points_.end(), single.begin(), single.end());
      if (out_groups.size() > 0 || prev_sfs.size() > 0)
	surfaceExtractOutput((int)kr, out_groups, prev_sfs);

#ifdef DEBUG_CHECK
      bool connect = regions_[kr]->isConnected();
      if (!connect)
	std::cout << "initPlaneCyl, disconnected region " << kr << std::endl;
#endif
      if (sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), sfs.begin(), sfs.end());
      if (repeat)
	--kr;
    }
  
  std::sort(regions_.begin(), regions_.end(), sort_region);

#ifdef DEBUG
  checkConsistence("Regions3");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 3" << std::endl;
      std::ofstream of("regions3.g2");
      std::ofstream ofm("mid_regions3.g2");
      std::ofstream ofs("small_regions3.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions3_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
    }
  
  if (single_points_.size() > 0)
    {
      std::ofstream of_single("single_pts3.g2");
      of_single << "400 1 0 4 0 0 0 255" << std::endl;
      of_single << single_points_.size() << std::endl;
      for (size_t kr=0; kr<single_points_.size(); ++kr)
	of_single << single_points_[kr]->getPoint() << std::endl;
     }

  if (surfaces_.size() > 0)
    {
      std::ofstream of("regsurf3.g2");
      writeRegionWithSurf(of);
    }
#endif
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 3. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post initPlaneCyl: " << ki << " " << ka << std::endl;
    }
#endif

  // Integrate single points when appropriate
  vector<RevEngPoint*> remaining_single;
  for (int ka=0; ka<(int)single_points_.size(); ++ka)
    {
      if (single_points_[ka]->isOutlier())
	continue;
      bool merged = single_points_[ka]->mergeWithAdjacent(mean_edge_len_);
      if (!merged)
	{
	  single_points_[ka]->setOutlier();
	  remaining_single.push_back(single_points_[ka]);
	}
    }
  std::swap(single_points_, remaining_single);

  // Sort regions according to number of points
  std::sort(regions_.begin(), regions_.end(), sort_region);

#ifdef DEBUG
  if (surfaces_.size() > 0)
    {
      std::ofstream of("regsurf4_0.g2");
      writeRegionWithSurf(of);
    }
#endif
	  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
}

//===========================================================================
void RevEng::growSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  double angtol = 5.0*anglim_;

  bool joinreg = true;
  if (joinreg)
    {
#ifdef DEBUG
      std::cout << "Pre join. Number of regions: " << regions_.size() << std::endl;
#endif
      // int num_reg = (int)regions_.size();
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
#ifdef DEBUG
      std::ofstream of0("init_join.g2");
      regions_[ki]->writeRegionInfo(of0);
#endif
	  if (regions_[ki]->numPoints() < min_point_region_)
	    continue;   // Grow into larger
	  
	  if (regions_[ki]->hasSurface())
	    {
	      growSurface(ki);
#ifdef DEBUG_CHECK
	      int num = regions_[ki]->numPoints();
	      for (int ka=0; ka<num; ++ka)
		if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
		  std::cout << "Inconsistent region pointers, post grow: " << ki << " " << ka << std::endl;
#endif
	    }
	  // else
	  //   {
	  //     vector<RevEngRegion*> adapted_regions;
	  //     regions_[ki]->joinRegions(mainaxis_, approx_tol_,
	  // 				angtol, adapted_regions);
	  //     for (size_t kj=0; kj<adapted_regions.size(); ++kj)
	  // 	{
	  // 	  size_t kr=0;
	  // 	  for (kr=0; kr<regions_.size(); ++kr)
	  // 	    if (adapted_regions[kj] == regions_[kr].get())
	  // 	      break;

	  // 	  if (kr < regions_.size())
	  // 	    {
	  // 	      // std::swap(regions_[kr], regions_[num_reg-1]);
	  // 	      // num_reg--;
	  // 	      regions_.erase(regions_.begin()+kr);
	  // 	    }
	  // 	}
	  //   }
#ifdef DEBUG
      std::ofstream of02("post_join.g2");
      regions_[ki]->writeRegionInfo(of02);
#endif
      int stop_grow = 1;
	}
      // if (num_reg < (int)regions_.size())
      //   regions_.erase(regions_.begin()+num_reg, regions_.end());
#ifdef DEBUG
      std::cout << "Post join. Number of regions: " << regions_.size() << std::endl;

      checkConsistence("Regions4");

      if (regions_.size() > 0)
	{
	  std::cout << "Regions 4" << std::endl;
	  std::ofstream of("regions4.g2");
	  std::ofstream ofm("mid_regions4.g2");
	  std::ofstream ofs("small_regions4.g2");
	  writeRegionStage(of, ofm, ofs);
	  if (surfaces_.size() > 0)
	    {
	      std::ofstream of("regsurf4.g2");
	      writeRegionWithSurf(of);
	    }
	  std::ofstream of0("regions4_helix.g2");
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
		{
		  regions_[ki]->writeRegionInfo(of0);
		  if (regions_[ki]->hasSurface())
		    regions_[ki]->writeSurface(of0);
		}
	    }
	}
#endif
    }
#ifdef DEBUG_CHECK
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post grow: " << ki << " " << ka << std::endl;
    }

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 4. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }


  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);

	}
    }
  
#ifdef DEBUG_CHECK
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post mergeAdjacentSimilar: " << ki << " " << ka << std::endl;
    }

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 4_2. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  
#ifdef DEBUG
      std::cout << "Post merge similar. Number of regions: " << regions_.size() << std::endl;

      checkConsistence("Regions4_2");

      if (regions_.size() > 0)
	{
	  std::cout << "Regions 4_2" << std::endl;
	  std::ofstream of("regions4_2.g2");
	  std::ofstream ofm("mid_regions4_2.g2");
	  std::ofstream ofs("small_regions4_2.g2");
	  writeRegionStage(of, ofm, ofs);
	  std::ofstream of0("regions4_2_helix.g2");
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
		{
		  regions_[ki]->writeRegionInfo(of0);
		  if (regions_[ki]->hasSurface())
		    regions_[ki]->writeSurface(of0);
		}
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf4_2.g2");
	  writeRegionWithSurf(of);
	}
#endif
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG
      std::ofstream ofp("init_growplane.g2");
      regions_[ki]->writeRegionInfo(ofp);
#endif
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;   // Not a stable source
	  
      if (!regions_[ki]->hasSurface())
	continue;

      vector<RevEngRegion*> grown_regions;
      vector<HedgeSurface*> adj_surfs;
      vector<vector<RevEngPoint*> > added_groups;
      regions_[ki]->growPlaneOrCyl(mainaxis_, min_point_region_, approx_tol_,
				   angtol, grown_regions, adj_surfs,
				   added_groups);
      updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
      vector<HedgeSurface*> dummy_surfs;
      if (added_groups.size() > 0)
	surfaceExtractOutput((int)ki, added_groups, dummy_surfs);
      int stop_growp = 1;
}

#ifdef DEBUG
      std::cout << "Grow plane. Number of regions: " << regions_.size() << std::endl;

      checkConsistence("Regions4_3");

      if (regions_.size() > 0)
	{
	  std::cout << "Regions 4_3" << std::endl;
	  std::ofstream of("regions4_3.g2");
	  std::ofstream ofm("mid_regions4_3.g2");
	  std::ofstream ofs("small_regions4_3.g2");
	  writeRegionStage(of, ofm, ofs);
	  std::ofstream of0("regions4_3_helix.g2");
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
		{
		  regions_[ki]->writeRegionInfo(of0);
		  if (regions_[ki]->hasSurface())
		    regions_[ki]->writeSurface(of0);
		}
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf4_3.g2");
	  writeRegionWithSurf(of);
	}
#endif
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  // Simplify regions structure
#ifdef DEBUG
  std::cout << "Number of regions pre integrate: " << regions_.size() << std::endl;
#endif
  int num_reg = (int)regions_.size();
  for (int ka=num_reg-1; ka>=0; --ka)
    {
      HedgeSurface* hedge = (regions_[ka]->hasSurface()) ?
	regions_[ka]->getSurface(0) : 0;
      bool to_be_removed = regions_[ka]->integrateInAdjacent(mean_edge_len_,
							     min_next_, max_next_,
							     approx_tol_, 0.5,
							     max_nmb_outlier_);
      if (to_be_removed)
	{
	  if (hedge)
	    {
	      size_t kr;
	      for (kr=0; kr<surfaces_.size(); ++kr)
		if (surfaces_[kr].get() == hedge)
		  break;
	      if (kr < surfaces_.size())
		surfaces_.erase(surfaces_.begin()+kr);
	    }

	  if (ka < num_reg-1)
	    std::swap(regions_[ka], regions_[num_reg-1]);
	  num_reg--;
	}
    }
  if (num_reg < (int)regions_.size())
    regions_.erase(regions_.begin()+num_reg, regions_.end());
#ifdef DEBUG
  std::cout << "Number of regions post integrate: " << regions_.size() << std::endl;
  
  checkConsistence("Regions5");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 5" << std::endl;
      std::ofstream of("regions5.g2");
      std::ofstream ofm("mid_regions5.g2");
      std::ofstream ofs("small_regions5.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf5.g2");
	  writeRegionWithSurf(of);
	}
     }
#endif
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, regions 5. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  
  int stop_break2 = 1;
}


//===========================================================================
void RevEng::updateAxesAndSurfaces()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // // First extract low accuracy points from groups with surface
  double angtol = 5.0*anglim_;
  // for (size_t ki=0; ki<regions_.size(); ++ki)
  //   {
  //     if (!regions_[ki]->hasSurface())
  // 	continue;

  //     vector<vector<RevEngPoint*> > added_groups;
  //     vector<HedgeSurface*> dummy_surfs;
  //     regions_[ki]->removeLowAccuracyPoints(min_point_region_, 
  // 					    approx_tol_, angtol, added_groups);
  //     if (added_groups.size() > 0)
  // 	surfaceExtractOutput((int)ki, added_groups, dummy_surfs);
  //   }

  std::sort(regions_.begin(), regions_.end(), sort_region);
  
  vector<int> reg_size(surfaces_.size());
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    reg_size[ki] = surfaces_[ki]->numPoints();

  std::sort(reg_size.begin(), reg_size.end());

  Point axis[3];
  int min_num = reg_size[(int)reg_size.size()/4];
  min_num = std::min(min_num, reg_size[reg_size.size()-1]/10);
  min_num = std::max(min_num, reg_size[reg_size.size()-1]/100);
  double max_ang = 0.1*M_PI;

  Point plane_axis[3];
  int num_pts1[3];
  computeAxisFromPlane(mainaxis_, min_num, max_ang, plane_axis, num_pts1);

  Point cyl_axis[3];
  int num_pts2[3];
  computeAxisFromCylinder(plane_axis, min_num, max_ang, cyl_axis, num_pts2);

  // Update main axes. Prioritize information from planes
  for (int ka=0; ka<3; ++ka)
    {
      num_pts1[ka] *= 2;
      int all_pts = num_pts1[ka] + num_pts2[ka];
      if (all_pts == 0)
	continue;
      double fac1 = (double)num_pts1[ka]/(double)(all_pts);
      double fac2 = (double)num_pts2[ka]/(double)(all_pts);
      if (cyl_axis[ka]*plane_axis[ka] < 0.0)
	cyl_axis[ka] *= -1.0;
      mainaxis_[ka] = fac1*plane_axis[ka] + fac2*cyl_axis[ka];
      mainaxis_[ka].normalize();
    }

  // Ensure orthogonality
  for (int ka=0; ka<3; ++ka)
    for (int kb=ka+1; kb<3; ++kb)
      if (num_pts1[ka] + num_pts2[ka] < num_pts1[kb] + num_pts2[kb])
	{
	  std::swap(num_pts1[ka], num_pts1[kb]);
	  std::swap(num_pts2[ka], num_pts2[kb]);
	  std::swap(mainaxis_[ka], mainaxis_[kb]);
	}

  Point tmp_axis[3];
  tmp_axis[2] = mainaxis_[0].cross(mainaxis_[1]);
  tmp_axis[1] = (num_pts1[2]+num_pts2[2] > 0) ?
	    mainaxis_[2].cross(mainaxis_[0]) : mainaxis_[1];
  tmp_axis[0] = (num_pts1[1]+num_pts2[1] > 0 && num_pts1[2]+num_pts2[2] > 0) ?
		 mainaxis_[1].cross(mainaxis_[2]) : mainaxis_[0];
  for (int ka=0; ka<3; ++ka)
    {
      if (tmp_axis[ka]*mainaxis_[ka] < 0.0)
	tmp_axis[ka] *= -1.0;
      if (num_pts1[ka] + num_pts2[ka] > 0)
	tmp_axis[ka] = 0.5*(tmp_axis[ka] + mainaxis_[ka]);
      tmp_axis[ka].normalize();
    }

  if (num_pts1[0]+num_pts2[0] > 0 && num_pts1[1]+num_pts2[1] > 0 &&
      num_pts1[2]+num_pts2[2] > 0)
    {
      mainaxis_[0] = tmp_axis[1].cross(tmp_axis[2]);
      mainaxis_[1] = tmp_axis[2].cross(tmp_axis[0]);
      mainaxis_[2] = tmp_axis[0].cross(tmp_axis[1]);
      for (int ka=0; ka<3; ++ka)
	{
	  if (mainaxis_[ka]*tmp_axis[ka] < 0.0)
	    mainaxis_[ka] *= -1.0;
	  mainaxis_[ka] = 0.5*(tmp_axis[ka] + mainaxis_[ka]);
	  mainaxis_[ka].normalize();
	}
    }
  else
    {
      for (int ka=0; ka<3; ++ka)
	{
	  mainaxis_[ka] = tmp_axis[ka];
	}
    }
  // Point tmp_axis = mainaxis_[0].cross(mainaxis_[1]);
  // tmp_axis.normalize();
  // if (mainaxis_[2]*tmp_axis < 0.0)
  //   tmp_axis *= -1;
  // mainaxis_[2] = 0.5*(mainaxis_[2] + tmp_axis);
  // tmp_axis = mainaxis_[2].cross(mainaxis_[0]);
  // tmp_axis.normalize();
  // if (mainaxis_[1]*tmp_axis < 0.0)
  //   tmp_axis *= -1;
  // mainaxis_[1] = 0.5*(mainaxis_[1] + tmp_axis);
  // tmp_axis = mainaxis_[1].cross(mainaxis_[2]);
  // tmp_axis.normalize();
  // if (mainaxis_[0]*tmp_axis < 0.0)
  //   tmp_axis *= -1;
  // mainaxis_[0] = 0.5*(mainaxis_[0] + tmp_axis);

  mainaxis_[2] = mainaxis_[0].cross(mainaxis_[1]);
  mainaxis_[1] = mainaxis_[2].cross(mainaxis_[0]);
  for (int ka=0; ka<3; ++ka)
    mainaxis_[ka].normalize();

  // Update surfaces
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())// && (!regions_[ki]->hasRevEdges()))
	{
	  bool updated = axisUpdate(ki, max_ang, angtol);
	  // if (updated)
	  //   {
	  //     growSurface(ki);
	  //     int stop_break0 = 1;
	  //   }
	}
    }
  
#ifdef DEBUG
  std::ofstream of0("regions5_2_0_helix.g2");
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	{
	  regions_[ki]->writeRegionInfo(of0);
	  if (regions_[ki]->hasSurface())
	    regions_[ki]->writeSurface(of0);
	}
    }
#endif
	  
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  {
	    updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	    // if (!regions_[ki]->hasRevEdges())
	    //   {
		bool updated = axisUpdate(ki, max_ang, angtol);
		if (!updated)
		  int stop_break1 = 1;
	      // }
	  }
	}
    }
#ifdef DEBUG
  std::cout << "Post merge similar. Number of regions: " << regions_.size() << std::endl;

  checkConsistence("Regions5_2");

  if (regions_.size() > 0)
    {
      std::cout << "Regions 5_2" << std::endl;
      std::ofstream of("regions5_2.g2");
      std::ofstream ofm("mid_regions5_2.g2");
      std::ofstream ofs("small_regions5_2.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions5_2_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf5_2.g2");
	  writeRegionWithSurf(of);
	}
    }
#endif
#ifdef DEBUG
  std::ofstream ofu1("unresolved.g2");
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      if (regions_[kr]->hasSurface())
	continue;
      if (regions_[kr]->hasAssociatedBlend())
	continue;
      regions_[kr]->writeRegionPoints(ofu1);
    }
#endif
   int stop_break = 1;
}


//===========================================================================
bool RevEng::axisUpdate(int ix, double max_ang, double angtol)
//===========================================================================
{
  if (!regions_[ix]->hasSurface())
    return false;

  HedgeSurface *hsurf = regions_[ix]->getSurface(0);
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(hsurf->surface());
  if (!elem.get())
    return false;

#ifdef DEBUG_AXIS
  std::ofstream of("axis_adapt.g2");
  elem->writeStandardHeader(of);
  elem->write(of);
  regions_[ix]->writeRegionPoints(of);
#endif
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base;
  regions_[ix]->getAdjacentElemInfo(adj_elem, adj_elem_base);
  Point adj_axis, adj_pos;
  double adj_ang = M_PI;
  Point vec = elem->direction();
  double pihalf = 0.5*M_PI;
  double min_perpen = M_PI;
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      if (adj_elem[ki].first->instanceType() != Class_Plane &&
	  adj_elem[ki].first->instanceType() != Class_Cylinder)
	continue;
      if (adj_elem[ki].second->hasBlendEdge())
	continue; // Derived information
      int sfflag = adj_elem[ki].second->getSurfaceFlag();
      double anglim = (sfflag == ACCURACY_OK) ? 2.0*max_ang : max_ang;
      Point dir = adj_elem[ki].first->direction();
      double ang = vec.angle(dir);
      if (fabs(pihalf - ang) < min_perpen)
	min_perpen = fabs(pihalf - ang);
      ang = std::min(ang, M_PI-ang);
      if (ang < anglim && ang < adj_ang)
	{
	  adj_ang = ang;
	  adj_axis = dir;
	  if (adj_elem[ki].first->instanceType() == Class_Cylinder)
	    adj_pos = adj_elem[ki].first->location();
	}
    }

  int ka_min = -1;
  double ang_main = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = vec.angle(mainaxis_[ka]);
      ang = std::min(ang, M_PI-ang);
      if (ang < max_ang && ang < ang_main)
	{
	  ang_main = ang;
	  ka_min = ka;
	}
    }
  
  bool updated = false;
  if (adj_axis.dimension() != 0 || ka_min >= 0) 
    updated = regions_[ix]->updateSurfaceWithAxis(min_point_region_, adj_axis, mainaxis_,
						  ka_min, approx_tol_, angtol, adj_pos);
  // if (updated == false && elem->instanceType() == Class_Plane && min_perpen < anglim_)
  //   {
  //     regions_[ix]->setPlaneParam(min_point_region_, mainaxis_, approx_tol_, angtol);
  //   }
  
  return updated;
}

//===========================================================================
void RevEng::recognizeEdges()
//===========================================================================
{
 // Ensure some limitation of surface size
  Point low = bbox_.low();
  Point high = bbox_.high();
  double diag = low.dist(high);
  double blendfac = 2.0;

  // Ensure bounded surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    surfaces_[ki]->limitSurf(diag);

  double angtol = 5.0*anglim_;
  double pihalf = 0.5*M_PI;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;
      if (regions_[ki]->hasBlendEdge())
	continue;

      vector<RevEngEdge*> rev_edgs1 = regions_[ki]->getAllRevEdges();
      
      int code;
      int classtype = regions_[ki]->getSurface(0)->instanceType(code);
      if (classtype != Class_Plane && classtype != Class_Cylinder && classtype != Class_Cone)
	continue;  // Preliminary

      shared_ptr<ParamSurface> surf1 = regions_[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem1 =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
      Point dir1 = elem1->direction();
      if (elem1->instanceType() == Class_Plane &&
	  dir1*regions_[ki]->getPoint(0)->getTriangNormal() < 0.0)
	dir1 *= -1;
      for (size_t kj=ki+1; kj<regions_.size(); ++kj)
	{
	  if (!regions_[kj]->hasSurface())
	    continue;
	  if (regions_[ki]->hasAssociatedBlend())
	    continue;
	  if (regions_[ki]->hasBlendEdge())
	    continue;
	  
	  vector<RevEngEdge*> rev_edgs2 = regions_[kj]->getAllRevEdges();
	  if (rev_edgs1.size() > 0 && rev_edgs2.size() > 0)
	    {
	      size_t kr, kh;
	      for (kr=0; kr<rev_edgs1.size(); ++kr)
		{
		  for (kh=0; kh<rev_edgs2.size(); ++kh)
		    if (rev_edgs1[kr] == rev_edgs2[kh])
		      break;
		  if (kh < rev_edgs2.size())
		    break;
		}
	      if (kr < rev_edgs1.size())
		continue;
	    }
	  
 	  int code2;
	  int classtype2 = regions_[kj]->getSurface(0)->instanceType(code2);
	  if (classtype2 != Class_Plane && classtype2 != Class_Cylinder &&
	      classtype2 != Class_Cone)
	    continue;  // Preliminary
	  if (classtype == Class_Cylinder && classtype2 == Class_Cylinder)
	    continue;
	  if (classtype == Class_Cone && classtype2 == Class_Cone)
	    continue;

	  shared_ptr<ParamSurface> surf2 = regions_[kj]->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> elem2 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
	  if (regions_[ki]->isAdjacent(regions_[kj].get())
	      || regions_[ki]->isNextToAdjacent(regions_[kj].get()))
	    {
#ifdef DEBUG_EDGE
	      std::ofstream of1("adj_regs.g2");
	      regions_[ki]->writeRegionPoints(of1);
	      regions_[ki]->writeSurface(of1);
	      regions_[kj]->writeRegionPoints(of1);
	      regions_[kj]->writeSurface(of1);
#endif
	      Point dir2 = elem2->direction();
	      if (elem2->instanceType() == Class_Plane &&
		  dir2*regions_[kj]->getPoint(0)->getTriangNormal() < 0.0)
		dir2 *= -1;
	      double ang = dir1.angle(dir2);
	      ang = std::min(ang, M_PI-ang);
	      bool compute_edge = false;
	      if (classtype == Class_Plane && classtype2 == Class_Plane)
		{
		  if (ang > blendfac*angtol)
		    compute_edge = true;
		}
	      else if (classtype == Class_Plane || classtype2 == Class_Plane)
		{
		  if (ang < blendfac*angtol)
		    compute_edge = true;
		  else if (fabs(pihalf-ang) < blendfac*angtol)
		    {
		      // Check for near tangential cases
		      Point norm = (classtype == Class_Plane) ? elem1->direction() :
			elem2->direction();
		      double dlen = fabs((elem1->location() - elem2->location())*norm);
		      double rad = (classtype == Class_Plane) ? elem2->radius(0.0, 0.0) :
			elem1->radius(0.0, 0.0);   // Not appropriate for a cone
		      if (dlen > rad + blendfac*approx_tol_ ||
			  dlen < rad - blendfac*approx_tol_)
			compute_edge = true;
		    }
		}
	      else
		{
		  // One cone and one cylinder. Check for same axis and
		  // a significant cone angle
		  Point loc1 = elem1->location();
		  Point loc2 = elem2->location();
		  Point tmp = loc1 + ((loc2-loc1)*dir1)*dir1;
		  if (ang >= blendfac*angtol || loc2.dist(tmp) > approx_tol_)
		    compute_edge = false;
		  else
		    {
		      double phi = 0.0;
		      shared_ptr<Cone> cone1 =
			dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
		      shared_ptr<Cone> cone2 =
			dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
		      if (cone1.get())
			phi = cone1->getConeAngle();
		      else if (cone2.get())
			phi = cone2->getConeAngle();
		      double conefac = 4.0;
		      if (fabs(phi) > conefac*angtol)
			compute_edge = true;
		      else
			compute_edge = false;
		    }
		}
	      

	      if (compute_edge)
		{
		  // Make sure that cone domains do not cover the apex
		  int ka;
		  shared_ptr<ElementarySurface> elem;
		  shared_ptr<RevEngRegion> reg;
		  for (ka=0, elem=elem1, reg=regions_[ki]; ka<2;
		       ++ka, elem=elem2, reg=regions_[kj])
		    {
		      if (elem->instanceType() == Class_Cone)
			{
			  shared_ptr<Cone> cone =
			    dynamic_pointer_cast<Cone,ElementarySurface>(elem);
			  double apar;
			  int adir;
			  cone->getDegenerateParam(apar, adir);
			  if (adir > 0)
			    {
			      RectDomain dom = cone->getParameterBounds();
			      double dom2[4];
			      dom2[0] = dom.umin();
			      dom2[1] = dom.umax();
			      dom2[2] = dom.vmin();
			      dom2[3] = dom.vmax();
			      double dom3[4];
			      reg->getDomain(dom3);
			      adir--;
			      // Assumes that the relevant part of the cone
			      // does not cover the apex
			      double midp = 0.5*(dom3[2*adir]+dom3[2*adir+1]);
			      double del = 0.01*(dom3[2*adir+1]-dom3[2*adir]);
			      if (midp < apar)
				dom2[2*adir+1] = apar - del;
			      else
				dom2[2*adir] = apar + del;
			      cone->setParameterBounds(dom2[0], dom2[2],
							dom2[1], dom2[3]);
			    }
			}
		    }

		  vector<shared_ptr<RevEngEdge> > edges =
		    defineEdgesBetween(ki, elem1, dir1,kj, elem2, dir2);
		  if (edges.size() > 0)
		    edges_.insert(edges_.end(), edges.begin(),edges.end());
		}
	    }
	}
    }
}

//===========================================================================
void RevEng::firstEdges()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  recognizeEdges();

#ifdef DEBUG
  checkConsistence("Regions6");

   if (regions_.size() > 0)
    {
      std::cout << "Regions6" << std::endl;
      std::ofstream of("regions6.g2");
      std::ofstream ofm("mid_regions6.g2");
      std::ofstream ofs("small_regions6.g2");
      writeRegionStage(of, ofm, ofs);
     }

   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges6.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	   int num_blend = edges_[kr]->numBlendRegs();
	   for (int ka=0; ka<num_blend; ++ka)
	     edges_[kr]->getBlendReg(ka)->writeRegionPoints(ofe);
	 }
     }
  
#endif

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

   
#ifdef DEBUG
   std::cout << "Extend blend region collection" << std::endl;
#endif
   for (size_t ki=0; ki<edges_.size(); ++ki)
     {
       extendBlendAssociation(ki);
     }
   
  // Just in case composed regions have been segmented
  int min_num_point = min_point_region_/10;
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;
      if (regions_[ki]->numPoints() < min_num_point)
	continue;
#ifdef DEBUG_EDGE
      std::ofstream of("planar_merge_cand.g2");
      regions_[ki]->writeRegionInfo(of);
#endif
      if (regions_[ki]->feasiblePlane(zero_H_, zero_K_))
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  bool merged = regions_[ki]->mergePlanarReg(zero_H_, zero_K_, 
						     approx_tol_, mainaxis_,
						     grown_regions);
	  if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	    updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);
	  regions_[ki]->updateRegionAdjacency();
	}
    }

  
 #ifdef DEBUG
  checkConsistence("Regions6_2");

   if (regions_.size() > 0)
    {
      std::cout << "Regions6_2" << std::endl;
      std::ofstream of("regions6_2.g2");
      std::ofstream ofm("mid_regions6_2.g2");
      std::ofstream ofs("small_regions6_2.g2");
      writeRegionStage(of, ofm, ofs);
     }

   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges6_2.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	   int num_blend = edges_[kr]->numBlendRegs();
	   for (int ka=0; ka<num_blend; ++ka)
	     edges_[kr]->getBlendReg(ka)->writeRegionPoints(ofe);
	 }
     }
  
#endif
   int stop_break = 1;
}

//===========================================================================
void RevEng::extendBlendAssociation(size_t ix)
//===========================================================================
{
  // Regions in blend area
  vector<RevEngRegion*> blend_regs;
  edges_[ix]->getAllBlendRegs(blend_regs);
  size_t num_blend_regs = blend_regs.size();

#ifdef DEBUG_BLEND
  std::ofstream of1("blend_regs.g2");
  for (size_t kj=0; kj<blend_regs.size(); ++kj)
    blend_regs[kj]->writeRegionPoints(of1);
#endif
  // Intersection curve
  vector<shared_ptr<CurveOnSurface> > cvs;
  edges_[ix]->getCurve(cvs);
  
  // Width
  double width = edges_[ix]->getDistance();

  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      vector<RevEngRegion*> new_blends;
      blend_regs[ki]->neighbourBlends(cvs, 2.0*width, approx_tol_, new_blends);
#ifdef DEBUG_BLEND
      std::ofstream of2("blend_regs2.g2");
      for (size_t kj=0; kj<new_blends.size(); ++kj)
	new_blends[kj]->writeRegionPoints(of2);
#endif
      for (size_t kj=0; kj<new_blends.size(); ++kj)
	{
	  size_t kr;
	  for (kr=num_blend_regs; kr<blend_regs.size(); ++kr)
	    if (blend_regs[kr] == new_blends[kj])
	      break;
	  if (kr == blend_regs.size())
	    blend_regs.push_back(new_blends[kj]);
	}
      int stop_break = 1;
    }

  for (size_t kj=num_blend_regs; kj<blend_regs.size(); ++kj)
    {
      edges_[ix]->addBlendRegion(blend_regs[kj]);
      blend_regs[kj]->setAssociatedBlend(edges_[ix].get());
    }
  int stop_break2 = 1;
}

//===========================================================================
bool RevEng::setBlendEdge(size_t ix)
//===========================================================================
{
 Point low = bbox_.low();
  Point high = bbox_.high();
  double diag = low.dist(high);
  double blendfac = 2.0;

  //double pihalf = 0.5*M_PI;
  double angtol = 5.0*anglim_;
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  regions_[ix]->getAdjacentElemInfo(adj_elem, adj_elem_base);
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      ClassType type1 = adj_elem[ki].first->instanceType();
      if (type1 != Class_Plane && type1 != Class_Cylinder)
	continue;
      Point dir1 = adj_elem[ki].first->direction();
      if (type1 == Class_Plane &&
	  dir1*adj_elem[ki].second->getMeanNormalTriang() < 0.0)
	dir1 *= -1;
      for (size_t kj=ki+1; kj<adj_elem.size(); ++kj)
	{
	  ClassType type2 = adj_elem[kj].first->instanceType();
	  if (type2 != Class_Plane && type2 != Class_Cylinder)
	    continue;
	  if (type1 == Class_Cylinder && type2 == Class_Cylinder)
	    continue;
	  Point dir2 = adj_elem[kj].first->direction();
	  if (type2 == Class_Plane &&
	      dir2*adj_elem[kj].second->getMeanNormalTriang() < 0.0)
	    dir2 *= -1;
	  double ang = dir1.angle(dir2);
	  ang = std::min(ang, M_PI-ang);
	  bool compute_edge = false;
	  if (type1 == Class_Plane && type2 == Class_Plane)
	    {
	      if (ang > blendfac*angtol)
		compute_edge = true;
	    }
	  else
	    {
	      if (ang < blendfac*angtol)
		compute_edge = true;
	    }

	  if (compute_edge)
	    {
	      // Check if the two regions already has a common edge
	      vector<RevEngEdge*> edges1 = adj_elem[ki].second->getAllRevEdges();
	      vector<RevEngEdge*> edges2 = adj_elem[kj].second->getAllRevEdges();
	      size_t kr, kh;
	      for (kr=0; kr<edges1.size(); ++kr)
		{
		  for (kh=0; kh<edges2.size(); ++kh)
		    if (edges1[kr] == edges2[kh])
		      break;
		  if (kh < edges2.size())
		    break;
		}

	      if (kr < edges1.size())
		{
		  // An edge exist already. Extend
		  continue;  // For the time being
		}
	      else
		{
		  // Define new edge
		  size_t ix1, ix2;
		  for (ix1=0; ix1<regions_.size(); ++ix1)
		    if (regions_[ix1].get() == adj_elem[ki].second)
		      break;
		  for (ix2=0; ix2<regions_.size(); ++ix2)
		    if (regions_[ix2].get() == adj_elem[kj].second)
		      break;
		  if (ix1 == regions_.size() || ix2 == regions_.size())
		    continue;

		  // Make sure that the adjacent surfaces is bounded
		  if (!adj_elem[ki].first->isBounded())
		    regions_[ix1]->getSurface(0)->limitSurf(diag);
		  if (!adj_elem[kj].first->isBounded())
		    regions_[ix2]->getSurface(0)->limitSurf(diag);
		  vector<shared_ptr<RevEngEdge> > edges =
		    defineEdgesBetween(ix1, adj_elem[ki].first,
				       dir1, ix2, adj_elem[kj].first, dir2);
		  if (edges.size() > 0)
		    {
		      edges_.insert(edges_.end(), edges.begin(), edges.end());
		      return true;
		    }
		}
	    }
	}
    }
  return false;
}

//===========================================================================
vector<shared_ptr<RevEngEdge> >
RevEng::defineEdgesBetween(size_t ix1, shared_ptr<ElementarySurface>& surf1,
			   Point& dir1, size_t ix2,
			   shared_ptr<ElementarySurface>& surf2, Point& dir2,
			   double lenlim0, bool check_common)
//===========================================================================
{
  vector<shared_ptr<RevEngEdge> > edges;
  if (regions_[ix1]->hasAssociatedBlend() || regions_[ix2]->hasAssociatedBlend())
    return edges;
  if ((regions_[ix1]->getSurfaceFlag() >= ACCURACY_POOR  &&
       (!regions_[ix1]->hasBlendEdge())) ||
      (regions_[ix2]->getSurfaceFlag() >= ACCURACY_POOR  &&
       (!regions_[ix2]->hasBlendEdge())))
    return edges;

  //double angtol = 5.0*anglim_;
  double blendlim = 50.0*mean_edge_len_;
  double lenlim = 10.0*mean_edge_len_; //blendlim;
  if (lenlim0 > 0.0)
    lenlim= std::min(lenlim, lenlim0);
  double int_tol = 1.0e-6;
  //bool plane1 = (surf1->instanceType() == Class_Plane);
  //bool plane2 = (surf2->instanceType() == Class_Plane);
  bool adjacent = regions_[ix1]->isAdjacent(regions_[ix2].get());

  shared_ptr<BoundedSurface> bd1, bd2;
  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
  BoundedUtils::getSurfaceIntersections(surf1, surf2, int_tol_,
					int_cvs1, bd1,
					int_cvs2, bd2);
  if (int_cvs1.size() == 0)
    return edges;

  bool keep_length = (int_cvs1.size() == 1 && (regions_[ix1]->hasBlendEdge() ||
						regions_[ix2]->hasBlendEdge()));
		  
  // Limit intersection curves to relevant intervals
  vector<pair<double,double> > t1_t2, t3_t4;
  bool OK1 =
    regions_[ix1]->getCurveRestriction(int_cvs1, approx_tol_,
				       anglim_, t1_t2);
  bool OK2 =
    regions_[ix2]->getCurveRestriction(int_cvs2, approx_tol_,
				       anglim_, t3_t4);

  if (!keep_length)
    {
      for (int ka=(int)int_cvs1.size()-1; ka>=0; --ka)
	{
	  double t1 = std::max(t1_t2[ka].first, t3_t4[ka].first);
	  double t2 = std::min(t1_t2[ka].second, t3_t4[ka].second);
	  if (t2 > t1 && (t1 > int_cvs1[ka]->startparam() ||
			  t2 < int_cvs1[ka]->endparam()))
	    {
	      double pmin = std::max(t1, int_cvs1[ka]->startparam());
	      double pmax = std::min(t2, int_cvs1[ka]->endparam());
	      shared_ptr<CurveOnSurface> sub1(int_cvs1[ka]->subCurve(pmin,pmax));
	      int_cvs1[ka] = sub1;
	      shared_ptr<CurveOnSurface> sub2(int_cvs2[ka]->subCurve(pmin,pmax));
	      int_cvs2[ka] = sub2;
	    }

	  if (t2 <= t1 || int_cvs1[ka]->estimatedCurveLength() < lenlim)
	    {
	      int_cvs1.erase(int_cvs1.begin()+ka);
	      int_cvs2.erase(int_cvs2.begin()+ka);
	    }
	}
      if (int_cvs1.size() == 0)
	return edges;
    }
  
  vector<RevEngRegion*> common_reg =
    regions_[ix1]->commonAdjacent(regions_[ix2].get());
  for (size_t kj=0; kj<common_reg.size(); )
    {
      if (common_reg[kj]->hasAssociatedBlend() || common_reg[kj]->hasBlendEdge())
  	common_reg.erase(common_reg.begin()+kj);
      else
  	++kj;
    }

  if (!adjacent && check_common)
    {
      // Check if any of the regions between is more significant
      int num1 = regions_[ix1]->numPoints();
      int num2 = regions_[ix2]->numPoints();
      int fac = 2;
      for (size_t kj=0; kj<common_reg.size(); ++kj)
	{
	  if (!common_reg[kj]->hasSurface())
	    continue;
	  int num3 = common_reg[kj]->numPoints();
	  if (num3 > fac*std::min(num1, num2))
	    return edges;  // Could be a small passage. 
	}
    }
		  
#ifdef DEBUG_EDGE
  std::ofstream of1("adj_regs_cv.g2");
  for (size_t kr=0; kr<int_cvs1.size(); ++kr)
    {
      shared_ptr<ParamCurve> cv = int_cvs1[kr]->spaceCurve();
      cv->writeStandardHeader(of1);
      cv->write(of1);
    }

  if (common_reg.size() > 0)
    {
      std::ofstream of2("between_regs.g2");
      for (size_t kr=0; kr<common_reg.size(); ++kr)
	{
	  common_reg[kr]->writeRegionPoints(of2);
	  if (common_reg[kr]->hasSurface())
	    common_reg[kr]->writeSurface(of2);
	}
    }
#endif

  vector<RevEngRegion*> regs1;
  regs1.push_back(regions_[ix2].get());
  regs1.insert(regs1.end(), common_reg.begin(), common_reg.end());
  vector<RevEngPoint*> bd_pts1 =
    regions_[ix1]->extractBdPoints(); //regs1);

  vector<RevEngRegion*> regs2;
  regs2.push_back(regions_[ix1].get());
  regs2.insert(regs2.end(), common_reg.begin(), common_reg.end());
  vector<RevEngPoint*> bd_pts2 =
    regions_[ix2]->extractBdPoints();//regs2);
  if (bd_pts1.size() == 0 || bd_pts2.size() == 0)
    return edges;
  
#ifdef DEBUG_EDGE
  std::ofstream of3("bd_pts.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << bd_pts1.size() << std::endl;
  for (size_t kr=0; kr<bd_pts1.size(); ++kr)
    of3 << bd_pts1[kr]->getPoint() << std::endl;
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << bd_pts2.size() << std::endl;
  for (size_t kr=0; kr<bd_pts2.size(); ++kr)
    of3 << bd_pts2[kr]->getPoint() << std::endl;
#endif

  // if (int_cvs1.size() > 1)
  //   {
  //     // Ensure consistent curve parameterization. Sort
  //     vector<Point> startpt(int_cvs1.size()), endpt(int_cvs1.size());
  //     for (size_t kr=0; kr<int_cvs1.size(); ++kr)
  // 	{
  // 	  startpt[kr] = int_cvs1[kr]->ParamCurve::point(int_cvs1[kr]->startparam());
  // 	  endpt[kr] = int_cvs1[kr]->ParamCurve::point(int_cvs1[kr]->endparam());
  // 	}
  //     int stop_breakp = 1;
  //   }
  int num_in_lim1=0, num_in_lim2=0;
  vector<pair<double, double> > t5_t6, t7_t8;
  vector<double> wwd1, wwd2;
  regions_[ix1]->estimateBlendDimensions(int_cvs1, bd_pts1,
					 approx_tol_, blendlim,
					 t5_t6, wwd1, num_in_lim1);

  regions_[ix2]->estimateBlendDimensions(int_cvs2, bd_pts2,
					 approx_tol_, blendlim,
					 t7_t8, wwd2, num_in_lim2);
  if (num_in_lim1 == 0 || num_in_lim2 == 0)
    return edges;
  if (t5_t6.size() == 0 || t7_t8.size() == 0)
    return edges;

  // Unify intersection curve limitations
  vector<pair<double,double> > tmin_tmax;
  vector<double> width;
  if (keep_length)
    {
      tmin_tmax.push_back(std::make_pair(int_cvs1[0]->startparam(),
					 int_cvs1[0]->endparam()));
      double wwd = 0.0;
      double w1min = std::numeric_limits<double>::max();
      double w2min = std::numeric_limits<double>::max();
      for (size_t kr=0; kr<wwd1.size(); ++kr)
	{
	  wwd += wwd1[kr];
	  w1min = std::min(w1min, wwd1[kr]);
	}
      for (size_t kr=0; kr<wwd2.size(); ++kr)
	{
	  wwd += wwd2[kr];
	  w2min = std::min(w2min, wwd2[kr]);
	}
      wwd /= (double)(wwd1.size()+wwd2.size());
      wwd = std::max(wwd, std::max(w1min, w2min));
      width.push_back(wwd);
    }
  else
    {
      size_t kk1, kk2;
      for (kk1=0, kk2=0; kk1<wwd1.size() && kk2<wwd2.size(); )
	{
	  for (; kk1<wwd1.size() && t7_t8[kk2].first >= t5_t6[kk1].second; ++kk1);
	  for (; kk2<wwd2.size() && t5_t6[kk1].first >= t7_t8[kk2].second; ++kk2);
	  tmin_tmax.push_back(std::make_pair(std::max(t5_t6[kk1].first,t7_t8[kk2].first),
					     std::min(t5_t6[kk1].second,t7_t8[kk2].second)));
	  width.push_back(0.5*(wwd1[kk1]+wwd2[kk2]));
	  if (t5_t6[kk1].second < t7_t8[kk2].second)
	    ++kk1;
	  else if (t7_t8[kk2].second < t5_t6[kk1].second)
	    ++kk2;
	  else
	    {
	      ++kk1;
	      ++kk2;
	    }
	}
      if (width.size() == 0)
	return edges;
    }

  // Don't. A rotational surface will still have a seam
  // bool connect_first = false;
  // if (width.size() > 1)
  //   {
  //     // Check for a connected piece across the seam of a closed intersection curve
  //     vector<Point> der1(2), der2(2);
  //     int ncvs = (int)int_cvs1.size();
  //     int_cvs1[0]->point(der1, int_cvs1[0]->startparam(), 1);
  //     int_cvs1[ncvs-1]->point(der2, int_cvs1[ncvs-1]->endparam(), 1);
  //     double dd = der1[0].dist(der2[0]);
  //     double angd = der1[1].angle(der2[1]);
  //     if (dd <= approx_tol_ && angd <= angtol)
  // 	connect_first = true;
  //   }

  size_t ki;
  //for (ki=connect_first ? 1 : 0; ki<width.size(); ++ki)
  for (ki=0; ki<width.size(); ++ki)
    {
      double tmin = tmin_tmax[ki].first;
      double tmax = tmin_tmax[ki].second;
      double width2 = width[ki];
      if (tmax <= tmin+int_tol)
	continue;
      
      vector<shared_ptr<CurveOnSurface> > cvs1, cvs2;
      for (size_t kr=0; kr<int_cvs1.size(); ++kr)
	{
	  double tp1 = std::max(int_cvs1[kr]->startparam(), tmin);
	  double tp2 = std::min(int_cvs1[kr]->endparam(), tmax);
	  if (tp2 <= tp1+int_tol)
	    continue;
	  if (tp1 > int_cvs1[kr]->startparam()+int_tol ||
	      tp2 < int_cvs1[kr]->endparam()-int_tol)
	    {
	      shared_ptr<CurveOnSurface> sub1(int_cvs1[kr]->subCurve(tp1, tp2));
	      shared_ptr<CurveOnSurface> sub2(int_cvs2[kr]->subCurve(tp1, tp2));
	      cvs1.push_back(sub1);
	      cvs2.push_back(sub2);
	    }
	  else if (fabs(tp1-int_cvs1[kr]->startparam()) < int_tol &&
		   fabs(tp2-int_cvs1[kr]->endparam()) < int_tol)
	    {
	      cvs1.push_back(int_cvs1[kr]);
	      cvs2.push_back(int_cvs2[kr]);

	    }
	}
    

      // if (connect_first && ki == width.size()-1 &&
      // 	  tmin_tmax[0].second-tmin_tmax[0].first > int_tol)
      // 	{
      // 	  for (size_t kr=0; kr<int_cvs1.size(); ++kr)
      // 	    {
      // 	      double tp1 = std::max(int_cvs1[kr]->startparam(), tmin_tmax[0].first);
      // 	      double tp2 = std::min(int_cvs1[kr]->endparam(), tmin_tmax[0].second);
      // 	      if (tp2 <= tp1+int_tol)
      // 		continue;
      // 	      shared_ptr<CurveOnSurface> sub1, sub2;
      // 	      if (tp1 > int_cvs1[kr]->startparam()+int_tol ||
      // 		  tp2 < int_cvs1[kr]->endparam()-int_tol)
      // 		{
      // 		  sub1 = shared_ptr<CurveOnSurface>(int_cvs1[kr]->subCurve(tp1, tp2));
      // 		  sub2 = shared_ptr<CurveOnSurface>(int_cvs2[kr]->subCurve(tp1, tp2));
      // 		}
      // 	      else if (fabs(tp1-int_cvs1[kr]->startparam()) < int_tol &&
      // 		       fabs(tp2-int_cvs1[kr]->endparam()) < int_tol)
      // 		{
      // 		  sub1 = int_cvs1[kr];
      // 		  sub2 = int_cvs2[kr];
      // 		}
      // 	      if (sub1.get())
      // 		{
      // 		  double dd1, dd2;
      // 		  cvs1[cvs1.size()-1]->appendCurve(sub1.get(), 0, dd1, false);
      // 		  cvs2[cvs2.size()-1]->appendCurve(sub2.get(), 0, dd1, false);
      // 		}
      // 	    }
      // 	  tmax = cvs1[cvs1.size()-1]->endparam();
      // 	  width2 = 0.5*(width2 + width[0]);
      // 	}
	  

      if (cvs1.size() == 0)
	continue;

#ifdef DEBUG_EDGE
      std::ofstream of1e("one_edgcv.g2");
      for (size_t kr=0; kr<cvs1.size(); ++kr)
	{
	  cvs1[kr]->spaceCurve()->writeStandardHeader(of1e);
	  cvs1[kr]->spaceCurve()->write(of1e);
	}
#endif
      for (size_t kj=0; kj<cvs1.size(); ++kj)
	{
	  shared_ptr<RevEngEdge> edg = defineOneEdge(ix1, surf1, dir1, ix2,
						     surf2, dir2, cvs1[kj],
						     cvs2[kj], width2,
						     common_reg);
	  if (edg.get())
	    edges.push_back(edg);
	}
    }

  return edges;
}
	

//===========================================================================
shared_ptr<RevEngEdge> 
RevEng::defineOneEdge(size_t ix1, shared_ptr<ElementarySurface>& surf1,
		      Point& dir1, size_t ix2,
		      shared_ptr<ElementarySurface>& surf2, Point& dir2,
		      shared_ptr<CurveOnSurface>& int_cv1,
		      shared_ptr<CurveOnSurface>& int_cv2,
		      double width, vector<RevEngRegion*>& common_reg)
//===========================================================================
{
  shared_ptr<RevEngEdge> dummy_edg;
  bool plane1 = (surf1->instanceType() == Class_Plane);
  bool plane2 = (surf2->instanceType() == Class_Plane);
  double angtol = 5.0*anglim_;
  double tol10 = 10.0*approx_tol_;
  int min_pt_blend = 20;

  if (plane1)
    {
      // Make sure that the plane normal poins out
      Point avnorm = regions_[ix1]->getMeanNormal();
      if (dir1*avnorm < 0.0)
	dir1 *= -1;
    }
  if (plane2)
    {
      // Make sure that the plane normal poins out
      Point avnorm = regions_[ix2]->getMeanNormal();
      if (dir2*avnorm < 0.0)
	dir2 *= -1;
    }
  int state = (plane1 || plane2) ? 1 : 2;


  vector<vector<RevEngPoint*> > near_pts(common_reg.size()+2);
  vector<RevEngPoint*> curr_near1, curr_near2;
  double tmin1 = int_cv1->startparam();
  double tmax1 = int_cv1->endparam();
  regions_[ix1]->getNearPoints(int_cv1, tmin1, tmax1, width,
			       angtol, curr_near1);
  double tmin2 = int_cv2->startparam();
  double tmax2 = int_cv2->endparam();
  regions_[ix2]->getNearPoints(int_cv2, tmin2, tmax2, width,
			       angtol, curr_near2);
  if (curr_near1.size() == 0 && curr_near2.size() == 0)
    return dummy_edg;

  RevEngPoint *distant1 = 0, *distant2 = 0;

  if (curr_near1.size() > 0)
    distant1 = getDistantPoint(int_cv1, std::max(tmin1, tmin2),
			       std::min(tmax1, tmax2), curr_near1);
  if (curr_near2.size() > 0)
    distant2 = getDistantPoint(int_cv2, std::max(tmin1, tmin2),
			       std::min(tmax1, tmax2), curr_near2);
  if (!distant1)
    {
      Point mid;
      int_cv1->point(mid, 0.5*(int_cv1->startparam()+int_cv1->endparam()));
      double distmid;
      distant1 = regions_[ix1]->closestPoint(mid, distmid);
    }
  if (!distant2)
    {
      Point mid;
      int_cv2->point(mid, 0.5*(int_cv2->startparam()+int_cv2->endparam()));
      double distmid;
      distant2 = regions_[ix2]->closestPoint(mid, distmid);
    }
  if ((!distant1) || (!distant2))
    return dummy_edg;

  vector<Point> der(2);
  Vector3D xyz1 = distant1->getPoint();
  Vector3D xyz2 = distant2->getPoint();
  Point loc1(xyz1[0], xyz1[1], xyz1[2]);
  Point loc2(xyz2[0], xyz2[1], xyz2[2]);
  Point distpt = (plane1) ? loc1 : loc2;
  double dist;
  double tpar;
  Point close;
  int_cv1->closestPoint(distpt, int_cv1->startparam(),
			int_cv1->endparam(), tpar, close, dist);
  int_cv1->point(der, tpar, 1);
		  
  Point loc1_2 = loc1 - ((loc1-der[0])*der[1])*der[1];
  Point loc2_2 = loc2 - ((loc2-der[0])*der[1])*der[1];
  Point vec = loc1_2 - loc2_2;
  bool outer1, outer2;
  Point pos2;
  if (plane1)
    outer1 = (dir1*vec < 0.0);
  else if (plane2)
    {
      // Cylinder or cone combined with plane
      Point pos = surf1->location();
       double vval = (der[0]-pos)*dir1;
     double rad = surf1->radius(0.0,vval);
      pos2 = pos + vval*dir1;
      Point loc2_3 = loc2 - ((loc2-der[0])*dir1)*dir1;
      outer1 = (pos2.dist(loc2_3) < rad);
    }
  else if (surf1->instanceType() == Class_Cylinder)
    {
      // Cylinder combined with cone
      Point vec2 = loc1 - der[0];
      Point axs = surf1->direction();
      outer1 = (vec*axs >= 0.0);
    }
  else
    {
      // Cone combined with cylinder
      Vector2D uv = distant1->getPar();
      double rad1 = surf1->radius(0.0, uv[1]);
      double rad2 = surf2->radius(0.0, 0.0);
      outer1 = (rad1 >= rad2);
    }
  
  if (plane2)
    outer2 = (dir2*vec > 0.0);
  else if (plane1)
    {
      Point pos = surf2->location();
      double vval = (der[0]-pos)*dir2;
      double rad = surf2->radius(0.0,vval);
      pos2 = pos + vval*dir2;
      Point loc1_3 = loc1 - ((loc1-der[0])*dir2)*dir2;
      outer2 = (pos2.dist(loc1_3) < rad);
    }
  else if (surf2->instanceType() == Class_Cylinder)
    {
      // Cylinder combined with cone
      Point vec2 = loc2 - der[0];
      Point axs = surf2->direction();
      outer2 = (vec*axs >= 0.0);
    }
  else
    {
      // Cone combined with cylinder
      Vector2D uv = distant2->getPar();
      double rad1 = surf2->radius(0.0, uv[1]);
      double rad2 = surf1->radius(0.0, 0.0);
      outer2 = (rad1 >= rad2);
    }

  double mindist1=0.0, mindist2=0.0;
  near_pts[0] = (curr_near1.size() == 0 || state == 2) ? curr_near1 :
    regions_[ix2]->removeOutOfSurf(curr_near1, tol10,
				   angtol, outer2, mindist1);
  near_pts[1] = (curr_near2.size() == 0 || state == 2) ? curr_near2 :
    regions_[ix1]->removeOutOfSurf(curr_near2, tol10,
				   angtol, outer1, mindist2);

  int nmb_near = (int)(near_pts[0].size() + near_pts[1].size());
  for (size_t kr=0; kr<common_reg.size(); ++kr)
    {
      vector<RevEngPoint*> adj_near1, adj_near2;
      double dummy_min = 0.0;
      double tmin3 = int_cv1->startparam();
      double tmax3 = int_cv1->endparam();
      common_reg[kr]->getNearPoints(int_cv1, tmin3, tmax3, width,
				    angtol, adj_near1);
      if (state == 1)
	{
	  adj_near2 =
	    regions_[ix1]->removeOutOfSurf(adj_near1, tol10,
					   angtol, outer1, dummy_min);
	  near_pts[2+kr] = 
	    regions_[ix2]->removeOutOfSurf(adj_near2, tol10,
					   angtol, outer2, dummy_min);
	}
      else
	near_pts[2+kr] = adj_near1;
      
      nmb_near += (int)near_pts[2+kr].size();
    }
  if (nmb_near < min_pt_blend)
    return dummy_edg;
  
#ifdef DEBUG_EDGE
  std::ofstream of4("blend_pts.g2");
  for (size_t kr=0; kr<near_pts.size(); ++kr)
    {
      of4 << "400 1 0 0" << std::endl;
      of4 << near_pts[kr].size() << std::endl;
      for (size_t kh=0; kh<near_pts[kr].size(); ++kh)
	of4 << near_pts[kr][kh]->getPoint() << std::endl;
    }
#endif

  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > blend_groups;
  for (size_t kr=0; kr<near_pts.size(); ++kr)
    blend_groups.push_back(std::make_pair(near_pts[kr].begin(),
					  near_pts[kr].end()));
  double radius = 0.0, ylen = 0.0;
  if (plane1 && plane2)
    {
      Point lin1 = der[1].cross(dir1);
      Point lin2 = der[1].cross(dir2);
      lin1.normalize();
      if (lin1*dir2 < 0.0)
	lin1 *= -1;
      lin2.normalize();
      if (lin2*dir1 < 0.0)
	lin2 *= -1;
      radius = computeCylinderRadius(near_pts, width, der[0],
				     der[1], lin1, lin2);
      Point vec = lin1 + lin2;
      vec.normalize();
      double alpha = lin1.angle(lin2);
      double fac = 1.0/sin(0.5*alpha);
      Point centre = der[0] - fac*radius*vec;
      //double ang = lin1.angle(lin2);
      double xlen = der[0].dist(centre);
      ylen = sqrt(xlen*xlen - radius*radius);
    }
  else
    {
      int sgn = 1;
      if ((plane1 && surf1->direction()*dir1 < 0.0) ||
	  (plane2 && surf2->direction()*dir2 < 0.0))
	sgn = -1;
      double d2 = 0.0;
      radius = computeTorusRadius(near_pts, int_cv1, surf1, surf2, 
				  width, outer1, outer2, sgn, d2);


      Point centre, normal, Cx;
      double Rrad;
      getTorusParameters(surf1, surf2, der[0], radius, d2, outer1, 
			 outer2, sgn, Rrad, centre, normal, Cx);
      double xlen = der[0].dist(centre);
      ylen = fabs(xlen - Rrad);
    }
  ylen *= 1.5;  // Allow some slack
  ylen= std::max(ylen, 1.1*std::max(mindist1, mindist2));

  // Reduce near points from adjacent groups according to the updated width
  if (ylen < width)
    {
      for (size_t kr=2; kr<near_pts.size(); ++kr)
	{
	  vector<RevEngPoint*> near_pts2;
	  common_reg[kr]->getNearPoints2(near_pts[kr], int_cv1, ylen, near_pts2);
	  std::swap(near_pts[kr], near_pts2);
	}
    }

#ifdef DEBUG_EDGE
  std::ofstream of5("blend_pts2.g2");
  for (size_t kr=0; kr<near_pts.size(); ++kr)
    {
      of5 << "400 1 0 0" << std::endl;
      of5 << near_pts[kr].size() << std::endl;
      for (size_t kh=0; kh<near_pts[kr].size(); ++kh)
	of5 << near_pts[kr][kh]->getPoint() << std::endl;
    }
#endif

  // An intersection/blend curve is found
  // Segment identified adjacent groups according to curve
  vector<RevEngRegion*> adj_regs;
  vector<size_t> remove_ix;
  for (size_t kr=0; kr<common_reg.size(); ++kr)
    {
      if (near_pts[kr+2].size() == 0)
	continue;

      if ((int)near_pts[kr+2].size() < common_reg[kr]->numPoints())
	{
	  int num_init = common_reg[kr]->numPoints();
	  vector<vector<RevEngPoint*> > out_groups;
	  vector<HedgeSurface*> out_sfs;
	  vector<vector<RevEngPoint*> > near_groups;
	  common_reg[kr]->extractSpesPoints(near_pts[kr+2], near_groups);
	  common_reg[kr]->updateInfo();
	  common_reg[kr]->splitRegion(out_groups);
	  if (common_reg[kr]->hasSurface() && common_reg[kr]->numPoints() < num_init/2)
	    {
	      int num_sf = common_reg[kr]->numSurface();
	      for (int ka=0; ka<num_sf; ++ka)
		out_sfs.push_back(common_reg[kr]->getSurface(ka));
	      common_reg[kr]->clearSurface();
	    }
	  vector<RevEngPoint*> all_near;
	  for (size_t kh=0; kh<near_groups.size(); ++kh)
	    all_near.insert(all_near.end(), near_groups[kh].begin(),
			    near_groups[kh].end());  // Can be disconnected

	  // Make new region
	  shared_ptr<RevEngRegion> curr_adj(new RevEngRegion(common_reg[kr]->getClassificationType(),
							     common_reg[kr]->getEdgeClassificationType(),
							     near_pts[kr+2]));
	  regions_.push_back(curr_adj);
	  adj_regs.push_back(curr_adj.get());

	  size_t kh=0;
	  for (kh=0; kh<regions_.size(); ++kh)
	    if (common_reg[kr] == regions_[kh].get())
	      break;
	  surfaceExtractOutput((int)kh, out_groups, out_sfs);
	}
      else
	{
	  adj_regs.push_back(common_reg[kr]);
	  remove_ix.push_back(kr);
	}
    }

  for (int ka=(int)remove_ix.size()-1; ka>=0; --ka)
    common_reg.erase(common_reg.begin() + remove_ix[ka]);

  for (size_t kr=0; kr<adj_regs.size(); ++kr)
    {
      vector<vector<RevEngPoint*> > curr_out;
      vector<HedgeSurface*> dummy_surfs;
      adj_regs[kr]->splitRegion(curr_out);
      if (curr_out.size() > 0)
	{
	  size_t kh=0;
	  for (kh=0; kh<regions_.size(); ++kh)
	    if (adj_regs[kr] == regions_[kh].get())
	      break;
	  surfaceExtractOutput((int)kh, curr_out, dummy_surfs);
	}
    }

  int edge_type = BLEND_NOT_SET;
  vector<shared_ptr<CurveOnSurface> > int_cvs1(1, int_cv1);
  vector<shared_ptr<CurveOnSurface> > int_cvs2(1, int_cv2);
  shared_ptr<RevEngEdge> edge(new RevEngEdge(edge_type, regions_[ix1].get(),
					     int_cvs1, outer1, regions_[ix2].get(),
					     int_cvs2, outer2, radius,
					     std::min(width, ylen)));
  regions_[ix1]->addRevEdge(edge.get());
  regions_[ix2]->addRevEdge(edge.get());
  if (adj_regs.size() > 0)
    edge->addBlendRegions(adj_regs);
  for (size_t kr=0; kr<adj_regs.size(); ++kr)
    adj_regs[kr]->setAssociatedBlend(edge.get());

  return edge;
}

					  
//===========================================================================
RevEngPoint* RevEng::getDistantPoint(shared_ptr<CurveOnSurface>& cv,
				     double tmin, double tmax,
				     vector<RevEngPoint*>& points)
//===========================================================================
{
  RevEngPoint *distant = 0;
  double ptdist = 0.0;
  double eps = 1.0e-9;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      Vector3D xyz = points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      cv->closestPoint(pt, cv->startparam(), cv->endparam(),
		       tpar, close, dist);
      if (tpar > tmin+eps && tpar < tmax-eps && dist > ptdist)
	{
	  ptdist = dist;
	  distant = points[ki];
	}
    }
  return distant;
}


//===========================================================================
void RevEng::updateRegionStructure()
//===========================================================================
{
  int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use

  // Segmentation of composed regions
#ifdef DEBUG
  std::cout << "Segment composed regions" << std::endl;
#endif
  double frac_norm_lim = 0.2; //0.75;
  double angtol = 5.0*anglim_;

  size_t reg_size = regions_.size();
  for (size_t ki=0; ki<reg_size; ++ki)
    {
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;

      if (regions_[ki]->hasSurface())
	continue;

      if (regions_[ki]->hasSweepInfo())
	continue;
      
#ifdef DEBUG_DIV
      std::ofstream ofs("region_to_segm.g2");
      regions_[ki]->writeRegionInfo(ofs);
      std::ofstream of2("unitsphere_segm.g2");
      regions_[ki]->writeUnitSphereInfo(of2);
      Point origo(0.0, 0.0, 0.0);
      of2 << "410 1 0 4 0 0 0 255" << std::endl;
      of2 << "3" << std::endl;
      for (int ka=0; ka<3; ++ka)
	of2 << origo << " " << mainaxis_[ka] << std::endl;

      std::ofstream ofa("adjacent_to_segm.g2");
      regions_[ki]->writeAdjacentPoints(ofa);
      //regions_[ki]->sortByAxis(mainaxis_, min_point_in, approx_tol_);
#endif
	  
      bool segmented = false;
      if (regions_[ki]->getFracNorm() > frac_norm_lim)
	{
	  // Grow sub regions according to surface type
#ifdef DEBUG_DIV
	  //std::cout << "Grow planes" << std::endl;
#endif
	  segmented = segmentByPlaneGrow((int)ki, min_point_in, angtol);
	}
      if (!segmented)
	{
	  // Check if a segmentation into several cylinder like
	  // regions is feasible
	  double avH, avK, MAH, MAK;
	  regions_[ki]->getAvCurvatureInfo(avH, avK, MAH, MAK);
	  double fac = 5.0;
	  if (MAH > fac*MAK)
	    {
	      segmented = segmentByAxis((int)ki, min_point_in);
	    }
	  if (!segmented)
	    segmented = segmentByContext((int)ki, min_point_in, angtol, true);
	}
    }

  // Update adjcency pointers
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
#ifdef DEBUG
  checkConsistence("Regions6_2");

   if (regions_.size() > 0)
    {
      std::cout << "Regions6_2" << std::endl;
      std::ofstream of("regions6_2.g2");
      std::ofstream ofm("mid_regions6_2.g2");
      std::ofstream ofs("small_regions6_2.g2");
      writeRegionStage(of, ofm, ofs);
    }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges6_2.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif
   
#ifdef DEBUG_DIV
  std::cout << "Merge adjacent regions" << std::endl;
#endif
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);

	}
    }
  
#ifdef DEBUG_DIV
  std::cout << "Merge adjacent regions" << std::endl;
#endif
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);

	}
    }
#ifdef DEBUG
  checkConsistence("Regions7");

   if (regions_.size() > 0)
    {
      std::cout << "Regions7" << std::endl;
      std::ofstream of("regions7.g2");
      std::ofstream ofm("mid_regions7.g2");
      std::ofstream ofs("small_regions7.g2");
      writeRegionStage(of, ofm, ofs);
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges7.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif   
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
  std::cout << "Pre join. Number of regions: " << regions_.size() << std::endl;
#endif
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;   // Grow into larger
      
      if (regions_[ki]->hasSurface())
	growSurface(ki);
    }
  
#ifdef DEBUG
  checkConsistence("Regions8");

   if (regions_.size() > 0)
    {
      std::cout << "Regions8" << std::endl;
      std::ofstream of("regions8.g2");
      std::ofstream ofm("mid_regions8.g2");
      std::ofstream ofs("small_regions8.g2");
      writeRegionStage(of, ofm, ofs);
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges8.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif
   
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }


}

//===========================================================================
void RevEng::computeAxisFromCylinder(Point initaxis[3], int min_num, double max_ang,
				     Point axis[3], int num_points[3])
//===========================================================================
{
  vector<vector<pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> > > points(3);
  num_points[0] =  num_points[1] = num_points[2] = 0;
  for (int ka=0; ka<3; ++ka)
    axis[ka] = initaxis[ka];

  Point axis0[3];
  for (int ka=0; ka<3; ++ka)
    axis0[ka] = initaxis[ka];
  double max_ang2 = 0.5*M_PI - 1.0e-9;
  for (int kc=0; kc<2; ++kc)
    {
      Point dummy(0.0, 0.0, 0.0);
      vector<double> min_angle(3, max_ang2);
      vector<Point> min_axis(3, dummy);
   
      for (size_t ki=0; ki<surfaces_.size(); ++ki)
	{
	  int sfcode;
	  ClassType type = surfaces_[ki]->instanceType(sfcode);
	  if ((type != Class_Cylinder) && (type != Class_Cone))
	    continue;
      
	  shared_ptr<ElementarySurface> curr =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surfaces_[ki]->surface());
	  if (!curr.get())
	    continue;

	  int num_reg = surfaces_[ki]->numRegions();
	  if (num_reg == 0)
	    continue; // Should not happen

	  RevEngRegion *reg0 = surfaces_[ki]->getRegion(0);
	  if (reg0->hasBlendEdge())
	    continue;  // Derived information
      
	  if (reg0->hasAssociatedBlend())
	    continue;  // Outdated information
      
	  int num_pts = surfaces_[ki]->numPoints();
	  if (num_pts < min_num)
	    continue;

	  Point vec = curr->direction();
	  int ix = -1;
	  double min_angle0 = std::numeric_limits<double>::max();
	  for (int ka=0; ka<3; ++ka)
	    {
	      double ang = vec.angle(axis0[ka]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang < min_angle0)
		{
		  min_angle0 = ang;
		  ix = ka;
		}
	      if (ang < max_ang)
		{
		  for (int kb=0; kb<num_reg; ++kb)
		    {
		      RevEngRegion *reg = surfaces_[ki]->getRegion(kb);
		      points[ka].push_back(std::make_pair(reg->pointsBegin(),
							  reg->pointsEnd()));
		    }
		  num_points[ka] += num_pts;
		}
	    }
	  if (ix >= 0)
	    {
	      min_angle[ix] = min_angle0;
	      min_axis[ix] = vec;
	    }
	}
      if (num_points[0]+num_points[1]+num_points[2] == 0 && kc == 0)
	{
	  for (int kb=0; kb<3; ++kb)
	    if (min_axis[kb].length() > 0.0)
	      axis0[kb] = min_axis[kb];
	}
      else if (kc == 0)
	break;
    }

  Point Cx, Cy;
  for (int ka=0; ka<3; ++ka)
    {
      if (points[ka].size() > 0)
	{
	  RevEngUtils::computeAxis(points[ka], axis[ka], Cx, Cy);
	}
      else
	axis[ka] = initaxis[ka];
    }

}

//===========================================================================
void RevEng::computeAxisFromPlane(Point initaxis[3], int min_num, double max_ang,
				  Point axis[3], int num_points[3])
//===========================================================================
{
  vector<vector<shared_ptr<Plane> > > planes(3);
  vector<vector<int> > num(3);
  vector<vector<double> > avd(3);
  num_points[0] =  num_points[1] = num_points[2] = 0;
  Point axis0[3];
  for (int ka=0; ka<3; ++ka)
    axis0[ka] = initaxis[ka];
  double max_ang2 = 0.5*M_PI - 1.0e-9;
  for (int kc=0; kc<2; ++kc)
    {
      Point dummy(0.0, 0.0, 0.0);
      vector<double> min_angle(3, max_ang2);
      vector<Point> min_axis(3, dummy);
      for (size_t ki=0; ki<surfaces_.size(); ++ki)
	{
	  int sfcode;
	  ClassType type = surfaces_[ki]->instanceType(sfcode);
	  if (type != Class_Plane)
	    continue;

	  shared_ptr<Plane> curr =
	    dynamic_pointer_cast<Plane,ParamSurface>(surfaces_[ki]->surface());
	  if (!curr.get())
	    continue;

	  int num_pts = surfaces_[ki]->numPoints();
	  if (num_pts < min_num)
	    continue;
      
	  // Check for a more accurate base surface
	  double avdist = 0.0;
	  int num_reg = surfaces_[ki]->numRegions();
	  double fac = 1.0/(double)num_pts;
	  if (num_reg > 0)
	    {
	      // Should always be the case
	      RevEngRegion* reg = surfaces_[ki]->getRegion(0);
	      shared_ptr<ParamSurface> base = reg->getBase();
	      if (base.get() && base->instanceType() == Class_Plane)
		{
		  double av1 = 0.0, av2 = 0.0;
		  //int nmb = 0;
		  for (int ka=0; ka<num_reg; ++ka)
		    {
		      reg = surfaces_[ki]->getRegion(ka);
		      double maxds, avds, maxdb, avdb;
		      int num_in, num2_in, num_inb, num2_inb;
		      reg->getAccuracy(maxds, avds, num_in, num2_in);
		      reg->getBaseDist(maxdb, avdb, num_inb, num2_inb);
		      int nn = reg->numPoints();
		      av1 += (double)nn*fac*avds;
		      av2 += (double)nn*fac*avdb;
		    }
		  if (av2 < av1)
		    {
		      curr = dynamic_pointer_cast<Plane,ParamSurface>(base);
		      avdist = av2;
		    }
		  else
		    avdist = av1;
		}
	    }

	  Point normal = curr->direction();
	  int ix = -1;
	  double min_angle0 = std::numeric_limits<double>::max();
	  for (int kb=0; kb<3; ++kb)
	    {
	      double ang = normal.angle(axis0[kb]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang < min_angle0)
		{
		  min_angle0 = ang;
		  ix = kb;
		}
	      if (ang <= max_ang)
		{
		  planes[kb].push_back(curr);
		  num[kb].push_back(num_pts);
		  avd[kb].push_back(avdist);
		  num_points[kb] += num_pts;
		}
	    }
	  if (ix >= 0)
	    {
	      min_angle[ix] = min_angle0;
	      min_axis[ix] = normal;
	    }
	}
      if (num_points[0]+num_points[1]+num_points[2] == 0 && kc == 0)
	{
	  for (int kb=0; kb<3; ++kb)
	    if (min_axis[kb].length() > 0.0)
	      axis0[kb] = min_axis[kb];
	}
      else if (kc == 0)
	break;
    }
  
  for (int kb=0; kb<3; ++kb)
    {
      if (planes[kb].size() == 0)
	{
	  axis[kb] = initaxis[kb];
	  continue;
	}

      axis[kb] = Point(0.0, 0.0, 0.0);
      double fac = 1.0/(double)num_points[kb];
      for (size_t ki=0; ki<planes[kb].size(); ++ki)
	{
	  Point normal = planes[kb][ki]->direction();
	  if (normal*axis0[kb] < 0.0)
	    normal *= -1.0;

	  double wgt = fac*(1.0 - avd[kb][ki])*num[kb][ki];
	  wgt = std::max(wgt, 0.0);
	  axis[kb] += wgt*normal;
	}
      axis[kb].normalize_checked();
    }
}

//===========================================================================
bool RevEng::recognizeOneSurface(int& ix, int min_point_in, double angtol,
				 int pass)
//===========================================================================
{
  bool firstpass = (pass <= 2);
  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) /*&&
								  regions_[ix]->feasiblePlane(zero_H_, zero_K_)*/)
    {
      vector<shared_ptr<HedgeSurface> > plane_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = regions_[ix]->extractPlane(mainaxis_,
					      approx_tol_, min_point_in,
					      min_point_region_, angtol, 
					      prefer_elementary_,
					      plane_sfs, prev_surfs, out_groups);

      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (plane_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;  // Result accepted
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass))
    {
      vector<shared_ptr<HedgeSurface> > cyl_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool repeat = false;
      
      // Cylinder from context information
      bool found =
	regions_[ix]->contextCylinder(mainaxis_, approx_tol_,
				      min_point_in, min_point_region_, 
				      angtol, prefer_elementary_,
				      cyl_sfs, prev_surfs,
				      out_groups);

      if (found == false && regions_[ix]->feasibleCylinder(zero_H_, zero_K_))
	(void)regions_[ix]->extractCylinder(approx_tol_, min_point_in,
					    min_point_region_, angtol,
					    prefer_elementary_,
					    cyl_sfs, prev_surfs,
					    out_groups, repeat);

      
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (repeat)
	{
	  --ix;
	  return false;  // New try
	}
	

      if (cyl_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }
  
  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && (!firstpass) &&
      regions_[ix]->hasSweepInfo() && regions_[ix]->sweepType() == 1)
    {
      vector<shared_ptr<HedgeSurface> > linsweep_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = regions_[ix]->extractLinearSweep(approx_tol_, min_point_in,
						 min_point_region_, angtol, 
						 prefer_elementary_,
						 linsweep_sfs, 
						 prev_surfs);
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (linsweep_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), linsweep_sfs.begin(), linsweep_sfs.end());
	   
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }
      
  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && pass > 1)
    {
      vector<shared_ptr<HedgeSurface> > tor_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;

      // Torus from context information
      bool found =  regions_[ix]->contextTorus(mainaxis_, approx_tol_,
					       min_point_in, min_point_region_, 
					       angtol, prefer_elementary_,
					       tor_sfs, prev_surfs,
					       out_groups);
      if (regions_[ix]->numPoints() == 0)
	{
	  regions_.erase(regions_.begin()+ix);
	  --ix;
	  return false;
	}

      if (!found)
	found = regions_[ix]->extractTorus(mainaxis_, approx_tol_, min_point_in,
					   min_point_region_,
					   angtol, prefer_elementary_,
					   tor_sfs, prev_surfs,
					   out_groups);

      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (tor_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), tor_sfs.begin(), tor_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && pass > 1)
    {
      vector<shared_ptr<HedgeSurface> > sph_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = regions_[ix]->extractSphere(mainaxis_, approx_tol_,
					       min_point_in, min_point_region_,
					       angtol, prefer_elementary_,
					       sph_sfs, prev_surfs,
					       out_groups);

      
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (sph_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), sph_sfs.begin(), sph_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && pass > 1)
    {
      vector<shared_ptr<HedgeSurface> > cone_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found =
	regions_[ix]->extractCone(approx_tol_, min_point_in, min_point_region_, 
				  angtol, prefer_elementary_,
				  cone_sfs, prev_surfs,
				  out_groups);
	   
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);
	
      if (cone_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), cone_sfs.begin(), cone_sfs.end());
      if (found && regions_[ix]->getMaxSfDist() < 0.5*approx_tol_)
	return true;
    }

  if (regions_[ix]->tryOtherSurf(prefer_elementary_, firstpass) && (!firstpass)) 
    {
      vector<shared_ptr<HedgeSurface> > spl_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
      bool found = 
	regions_[ix]->extractFreeform(approx_tol_, min_point_in,
				      min_point_region_, angtol,
				      prefer_elementary_,
				      spl_sfs, prev_surfs,
				      out_groups);

      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (spl_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), spl_sfs.begin(), spl_sfs.end());
    }

  if (!regions_[ix]->hasSurface())
    {
      vector<shared_ptr<HedgeSurface> > adj_sfs;
      vector<HedgeSurface*> prev_surfs;
      vector<vector<RevEngPoint*> > out_groups;
       bool found =
	regions_[ix]->adjacentToCylinder(mainaxis_, approx_tol_, min_point_in,
					 min_point_region_, angtol,
					 prefer_elementary_,
					 adj_sfs, prev_surfs,
					 out_groups);
      
      if (out_groups.size() > 0 || prev_surfs.size() > 0)
	surfaceExtractOutput(ix, out_groups, prev_surfs);

      if (adj_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), adj_sfs.begin(), adj_sfs.end());
    }
#ifdef DEBUG
  if (regions_[ix]->hasSurface())
    {
      std::ofstream of("one_surface.g2");
      regions_[ix]->writeSurface(of);
    }
#endif
  return (regions_[ix]->hasSurface());
}

//===========================================================================
void RevEng::recognizeSurfaces(int min_point_in, int pass)
//===========================================================================
{
  double angfac = 5.0;
  double angtol = angfac*anglim_;
  int pass2 = pass + 1;
  std::sort(regions_.begin(), regions_.end(), sort_region);
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, pre recognizeSurfaces. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, pre recognizeSurfaces: " << ki << " " << ka << std::endl;
    }
#endif

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG_CHECK
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, pre. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      
 #endif
    }
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG_CHECK
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      for (size_t kj=0; kj<regions_.size(); ++kj)
	{
	  int nump = regions_[kj]->numPoints();
	  for (int ka=0; ka<nump; ++ka)
	    if (regions_[kj]->getPoint(ka)->region() != regions_[kj].get())
	      std::cout << "Inconsistent region pointers, recognizeSurfaces: " << ki << " " << kj << " " << ka << std::endl;
	}
#endif
      if (regions_[ki]->numPoints() < min_point_region_)
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;  // Treated separately
#ifdef DEBUGONE
      //int classtype = regions_[ki]->getClassification();
      //std::cout << "No " << ki <<  ", classtype: " << classtype << std::endl;
      
      std::ofstream of1("region.g2");
      regions_[ki]->writeRegionInfo(of1);
      std::ofstream of2("unitsphere.g2");
      regions_[ki]->writeUnitSphereInfo(of2);
#endif
      if (regions_[ki]->hasBlendEdge())
	continue;
      if (regions_[ki]->hasRevEdges())
	continue;
      
      bool pot_blend = false;
      if (!regions_[ki]->hasSurface())
	pot_blend = regions_[ki]->potentialBlend(angtol);
      if (pot_blend)
	{
#ifdef DEBUG
	  std::cout << "Potential blend" << std::endl;
#endif
	  bool done = setBlendEdge(ki);
	  if (done)
	    continue;
	}
      bool found = recognizeOneSurface(ki, min_point_in, angtol, pass2);
#ifdef DEBUG_CHECK
      bool connect = regions_[ki]->isConnected();
      if (!connect)
	std::cout << "recognizeOneSurface, disconnected region " << ki << std::endl;
#endif
        
#ifdef DEBUG_CHECK
  for (int kj=0; kj<(int)regions_.size(); ++kj)
    {
      std::set<RevEngPoint*> tmpset2(regions_[kj]->pointsBegin(), regions_[kj]->pointsEnd());
      if ((int)tmpset2.size() != regions_[kj]->numPoints())
	std::cout << "Point number mismatch 2. " << kj << " " << ki << " " << tmpset2.size() << " " << regions_[kj]->numPoints() << std::endl;
    }
#endif
  if ((!regions_[ki]->hasSurface()) ||
      (regions_[ki]->hasSurface() &&
       (regions_[ki]->getSurfaceFlag() == ACCURACY_POOR ||
	regions_[ki]->getSurfaceFlag() == NOT_SET)))
	{
	  // Still no surface. Try to divide composite regions into smaller
	  // pieces
	  bool split = segmentComposite(ki, min_point_in, angtol);
#ifdef DEBUG_CHECK
	  if (ki >= 0)
	    {
	      bool connect = regions_[ki]->isConnected();
	      if (!connect)
		std::cout << "segmentComposite, disconnected region " << ki << std::endl;
	    }
#endif
	  int stop_split = 1;
	}
   }

#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, post recognizeSurfaces. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, post recognizeSurfaces: " << ki << " " << ka << std::endl;
    }
#endif

  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

}

//===========================================================================
void RevEng::validateSurfaces()
//===========================================================================
{
#ifdef DEBUG_VALIDATE
  std::ofstream of0("init_helix.g2");
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	{
	  regions_[ki]->writeRegionInfo(of0);
	  if (regions_[ki]->hasSurface())
	    regions_[ki]->writeSurface(of0);
	}
    }
#endif
  int min_point_in = 50;   // Should be a class parameter
  double angfac = 5.0;
  double angtol = angfac*anglim_;
  double fracfac = 1.5;
  std::sort(regions_.begin(), regions_.end(), sort_region);
  bool do_adjust = true;
  int no_adjust = -1;
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface() &&
	  regions_[ki]->getSurfaceFlag() > ACCURACY_OK)
	{
	  if (regions_[ki]->hasBlendEdge())
	    continue;
	  shared_ptr<ParamSurface> surf = regions_[ki]->getSurface(0)->surface();
	  int classtype = surf->instanceType();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
	  if (!elem.get())
	    continue;
	  
	  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;;
	  regions_[ki]->getAdjacentElemInfo(adj_elem, adj_elem_base);
	  
#ifdef DEBUG_VALIDATE
	  std::ofstream of1("validatesurf.g2");
	  regions_[ki]->writeRegionInfo(of1);
	  regions_[ki]->writeSurface(of1);
	  std::ofstream of2("adj_validate.g2");
	  for (size_t kj=0; kj<adj_elem.size(); ++kj)
	    {
	      adj_elem[kj].second->writeRegionPoints(of2);
	      adj_elem[kj].second->writeSurface(of2);
	    }
#endif
	  // Check if the region is connected
	  if (!regions_[ki]->isConnected())
	    {
	      vector<vector<RevEngPoint*> > sep_groups;
	      regions_[ki]->splitRegion(sep_groups);
	      if (sep_groups.size() > 0)
		{
		  vector<HedgeSurface*> dummy_sfs;
		  surfaceExtractOutput(ki, sep_groups, dummy_sfs);
		  regions_[ki]->updateInfo(approx_tol_, angtol);
		  regions_[ki]->checkReplaceSurf(mainaxis_, min_point_region_,
						 approx_tol_, angtol);
		  --ki;
		  continue;
		}
	    }
	  
	  double in_frac1 = regions_[ki]->getFracNorm();
	  double in_frac2 = regions_[ki]->getFracNormTriang();
	  if (classtype == Class_Plane && in_frac1 > fracfac*in_frac2)
	    {
	      bool update = integratePlanarInHelix(ki, elem, adj_elem);
	      if (update)
		continue;
	    }

	  // Check for significant axis
	  if (adj_elem.size() > 0)
	    {
	      if (no_adjust != ki)
		{
		  bool update = adjustWithAdjacent(ki, min_point_in,
						   angtol, adj_elem);
		  if (update)
		    {
		      no_adjust = ki;
		      --ki;
		      continue;
		    }
		}

	      Point axis, axis2, pos;
	      bool found = regions_[ki]->identifySignificantAxis(adj_elem, pos,
								 axis, axis2);
	      if (found)
		{
		  regions_[ki]->analyseRotated(pos, axis, axis2);
		  int stop_break0 = 1;
		}
	    }
	  int stop_break = 1;
	}
    }
#ifdef DEBUG
  checkConsistence("Regions12");

   std::cout << "Finished validate surf, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions12" << std::endl;
      std::ofstream of("regions12.g2");
      std::ofstream ofm("mid_regions12.g2");
      std::ofstream ofs("small_regions12.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions12_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf12.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges12.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif
}

//===========================================================================
bool RevEng::adjustWithAdjacent(int& ix, int min_point_in, double angtol,
				vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj)
//===========================================================================
{
  // Identify reliable adjacent surfaces
  vector<RevEngRegion*> adj_reg;
  for (size_t ki=0; ki<adj.size(); ++ki)
    if (adj[ki].second->getSurfaceFlag() < ACCURACY_POOR)
      adj_reg.push_back(adj[ki].second);

  if (adj_reg.size() == 0)
    return false;

  vector<vector<RevEngPoint*> > added_groups;
  int num_pt = regions_[ix]->numPoints();
  bool changed = regions_[ix]->segmentByAdjSfContext(mainaxis_, min_point_in,
						     min_point_region_,
						     approx_tol_, angtol,
						     adj_reg, added_groups);

  if (added_groups.size() > 0)
    {
      vector<HedgeSurface*> dummy_sfs;
      surfaceExtractOutput(ix, added_groups, dummy_sfs);
    }

  // Check if remaining region is connected
  if (!regions_[ix]->isConnected())
    {
      vector<vector<RevEngPoint*> > sep_groups;
      regions_[ix]->splitRegion(sep_groups);
      if (sep_groups.size() > 0)
	{
	  vector<HedgeSurface*> dummy_sfs;
	  surfaceExtractOutput(ix, sep_groups, dummy_sfs);
	}
    }

  if (regions_[ix]->numPoints() < num_pt)
    {
      regions_[ix]->updateInfo(approx_tol_, angtol);
      regions_[ix]->checkReplaceSurf(mainaxis_, min_point_region_,
				     approx_tol_, angtol);
      return true;
    }
  else
    return false;
}


//===========================================================================
bool RevEng::integratePlanarInHelix(int& ix, shared_ptr<ElementarySurface> surf,
				    vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj)
//===========================================================================
{
  // Check for a feasible adjacent regions
  double pihalf = 0.5*M_PI;
  double anglim = 0.1*M_PI;
  vector<RevEngRegion*> adj_reg;
  Point norm = surf->direction();
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if ((adj[ki].first->instanceType() == Class_Cylinder ||
	   adj[ki].first->instanceType() == Class_Cone) &&
	  adj[ki].second->getSurfaceFlag() == PROBABLE_HELIX)
	{
	  // Check for roughly orthogonal axis/normal
	  Point axis = adj[ki].first->direction();
	  double ang = norm.angle(axis);
	  if (fabs(pihalf-ang) <= anglim)
	    adj_reg.push_back(adj[ki].second);
	}
    }

  if (adj_reg.size() == 0)
    return false;  // No feasible candidate

  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_reg.size(); ++kj)
      if (adj_reg[kj]->numPoints() > adj_reg[ki]->numPoints())
	std::swap(adj_reg[ki], adj_reg[kj]);

  double angtol = 5.0*anglim_;
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    {
      vector<RevEngRegion*> grown_regions;
      vector<HedgeSurface*> adj_surfs;
      bool done = adj_reg[ki]->includeAdjacent(regions_[ix].get(),
					       mainaxis_, approx_tol_,
					       angtol, grown_regions,
					       adj_surfs);
      if (done)
	{
#ifdef DEBUG_VALIDATE
	  std::ofstream of("update_val.g2");
	  adj_reg[ki]->writeRegionPoints(of);
	  adj_reg[ki]->writeSurface(of);
#endif
	  int ix2 = ix+1;  // To force the function to remove
	  updateRegionsAndSurfaces(ix2, grown_regions, adj_surfs);
	  --ix;
	  return true;
	}
    }
  return false;
}

//===========================================================================
bool RevEng::segmentComposite(int& ix, int min_point_in, double angtol)
//===========================================================================
{
  // To be updated
  double frac_norm_lim = 0.2; //0.75;
  bool segmented = false;
  
  if (regions_[ix]->numPoints() < min_point_region_)
    return false;

  if (regions_[ix]->hasSurface() &&
      (!(regions_[ix]->getSurfaceFlag() == ACCURACY_POOR ||
	 regions_[ix]->getSurfaceFlag() == NOT_SET)))
    return false;

  if (regions_[ix]->hasSweepInfo())
    return false;

#ifdef DEBUG_DIV
  std::cout << "Segment composite, ix=" << ix << std::endl;
  std::ofstream ofs("region_to_segm.g2");
  regions_[ix]->writeRegionInfo(ofs);
  std::ofstream of2("unitsphere_segm.g2");
  regions_[ix]->writeUnitSphereInfo(of2);
  Point origo(0.0, 0.0, 0.0);
  of2 << "410 1 0 4 0 0 0 255" << std::endl;
  of2 << "3" << std::endl;
  for (int ka=0; ka<3; ++ka)
    of2 << origo << " " << mainaxis_[ka] << std::endl;

  std::ofstream ofa("adjacent_to_segm.g2");
  regions_[ix]->writeAdjacentPoints(ofa);
  //regions_[ix]->sortByAxis(mainaxis_, min_point_in, approx_tol_);
#endif

  size_t num_points = regions_[ix]->numPoints();
  
  bool plane_grow = false;
  if (plane_grow && regions_[ix]->getFracNorm() > frac_norm_lim)
    {
      // Grow sub regions according to surface type
#ifdef DEBUG_DIV
      //std::cout << "Grow planes" << std::endl;
#endif
      segmented = segmentByPlaneGrow(ix, min_point_in, angtol);
#ifdef DEBUG_DIV
      if (segmented)
	std::cout << "Segmented by grow planes" << std::endl;
#endif
    }

  bool repeat = false;
  vector<shared_ptr<HedgeSurface> > hedgesfs;
  vector<shared_ptr<RevEngRegion> > added_reg;
  vector<vector<RevEngPoint*> > separate_groups;
  vector<RevEngPoint*> single_points;
  vector<HedgeSurface*> prevsfs;

  if (!segmented && regions_[ix]->hasDivideInfo())
    {
      for (int ki=0; ki<regions_[ix]->numDivideInfo(); ++ki)
	{
	  segmented = regions_[ix]->divideWithSegInfo(ki,
						      min_point_region_,
						      separate_groups,
						      single_points);
	  if (segmented)
	    break;
	}
    }
  
  if (!segmented)
    {
      // Check if a segmentation into several cylinder like
      // regions is feasible
      double avH, avK, MAH, MAK;
      regions_[ix]->getAvCurvatureInfo(avH, avK, MAH, MAK);
      double fac = 5.0;
      if (MAH > fac*MAK)
	{
	  segmented = regions_[ix]->extractCylByAxis(mainaxis_, min_point_in,
						     min_point_region_,
						     approx_tol_, angtol,
						     prefer_elementary_,
						     hedgesfs, added_reg,
						     separate_groups,
						     single_points);
	}
    }

#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset1(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset1.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 1. " << ki << " " << tmpset1.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
    
  if (!segmented)
    {
      vector<RevEngRegion*> adj_planar = regions_[ix]->fetchAdjacentPlanar();
      if (adj_planar.size() > 1)
	{
	  segmented = regions_[ix]->segmentByPlaneAxis(mainaxis_, min_point_in,
						       min_point_region_,
						       approx_tol_, angtol,
						       prefer_elementary_,
						       adj_planar, hedgesfs,
						       added_reg,
						       prevsfs, separate_groups);
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset2(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset2.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 2. " << ki << " " << tmpset2.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
    
  for (size_t kr=0; kr<adj_planar.size(); ++kr)
    {
      std::set<RevEngPoint*> tmpset(adj_planar[kr]->pointsBegin(),
				    adj_planar[kr]->pointsEnd());
      if ((int)tmpset.size() != adj_planar[kr]->numPoints())
	std::cout << "Point number mismatch (ByPlaneAxis). " << kr << " " << tmpset.size() << " " << adj_planar[kr]->numPoints() << std::endl;
    }
	  
#endif
	}
      
      if (!segmented)
	{
	  // Extend with cylindrical
	  vector<RevEngRegion*> adj_cyl =
	    regions_[ix]->fetchAdjacentCylindrical();
	  if (adj_cyl.size() > 0)
	    adj_planar.insert(adj_planar.end(), adj_cyl.begin(),
			      adj_cyl.end());
	  if (adj_planar.size() > 0)
	    segmented =
	      regions_[ix]->segmentByAdjSfContext(mainaxis_, min_point_in,
						  min_point_region_,
						  approx_tol_, angtol,
						  adj_planar, separate_groups);
#ifdef DEBUG_CHECK
	  for (int ki=0; ki<(int)regions_.size(); ++ki)
	    {
	      std::set<RevEngPoint*> tmpset4(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
	      if ((int)tmpset4.size() != regions_[ki]->numPoints())
		std::cout << "Point number mismatch, composite 4. " << ki << " " << tmpset4.size() << " " << regions_[ki]->numPoints() << std::endl;
	    }
	  for (size_t kr=0; kr<adj_planar.size(); ++kr)
	    {
	      std::set<RevEngPoint*> tmpset(adj_planar[kr]->pointsBegin(),
					    adj_planar[kr]->pointsEnd());
	      if ((int)tmpset.size() != adj_planar[kr]->numPoints())
		std::cout << "Point number mismatch (ByAdjSfContext). " << kr << " " << tmpset.size() << " " << adj_planar[kr]->numPoints() << std::endl;
	    }
#endif
	}
    }
  
  if (added_reg.size() > 0)
    repeat = true;
  else if (separate_groups.size() > 0)
    {
      int num = 0;
      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	num += (int)separate_groups[ki].size();
      if (num > (int)num_points/10)
	repeat = true;
    }
  
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset5(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset5.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 5. " << ki << " " << tmpset5.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
#ifdef DEBUG
  if (regions_[ix]->numPoints() < (int)num_points)
    {
      std::ofstream of("seg_by_context.g2");
      int num = regions_[ix]->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of <<  num << std::endl;
      for (int ka=0; ka<num; ++ka)
	of << regions_[ix]->getPoint(ka)->getPoint() << std::endl;

      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of <<  separate_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)separate_groups[ki].size(); ++ka)
	    of << separate_groups[ki][ka]->getPoint() << std::endl;
	}
      
      for (size_t ki=0; ki<added_reg.size(); ++ki)
	{
	  int num = added_reg[ki]->numPoints();
	  of << "400 1 0 4 0 0 255 255" << std::endl;
	  of <<  num << std::endl;
	  for (int ka=0; ka<num; ++ka)
	    of << added_reg[ki]->getPoint(ka)->getPoint() << std::endl;
	  if (added_reg[ki]->hasSurface())
	    added_reg[ki]->writeSurface(of);
	}
    }
#endif

  // Update adjacency information for current region
  regions_[ix]->updateRegionAdjacency();
  
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset6(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset6.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 6. " << ki << " " << tmpset6.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  if (added_reg.size() > 0 && regions_[ix]->hasSurface() == false)
    {
      // Swap sequence of regions
      int max_pt = 0;
      int max_ix = -1;
      for (size_t ki=0; ki<added_reg.size(); ++ki)
	if (added_reg[ki]->numPoints() > max_pt)
	  {
	    max_pt = added_reg[ki]->numPoints();
	    max_ix = (int)ki;
	  }
      if (max_ix >= 0)
	std::swap(regions_[ix], added_reg[max_ix]);
    }
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset6_2(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset6_2.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 6_2. " << ki << " " << tmpset6_2.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
  if (separate_groups.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, separate_groups, prev_surfs);
    }

#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset7(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if ((int)tmpset7.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 7. " << ki << " " << tmpset7.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
  if (!regions_[ix]->isConnected())
    std::cout << "Disconnected region (split), ix= " << ix << " " << regions_[ix].get() << std::endl;
#endif
  if (added_reg.size() > 0)
    regions_.insert(regions_.end(), added_reg.begin(), added_reg.end());
  if (hedgesfs.size() > 0)
    surfaces_.insert(surfaces_.end(), hedgesfs.begin(), hedgesfs.end());
  if (repeat && (!regions_[ix]->hasSurface()))
    {
#ifdef DEBUG_DIV
      std::cout << "Repeat ix = " << ix << std::endl;
#endif
      --ix;
    }
  
  return segmented;
}

//===========================================================================
void RevEng::adjustPointRegions(int min_point_in)
//===========================================================================
{
  double angfac = 5.0;
  double angtol = angfac*anglim_;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	continue;
      if (regions_[ki]->hasAssociatedBlend())
	continue;
      if (regions_[ki]->numPoints() < min_point_in)
	continue;
      
      vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base;
      regions_[ki]->getAdjacentElemInfo(adj_elem, adj_elem_base);
      if (adj_elem.size() == 0)
	continue;

      vector<RevEngRegion*> adj_groups(adj_elem.size());
      for (size_t kj=0; kj<adj_elem.size(); ++kj)
	adj_groups[kj] = adj_elem[kj].second;

      vector<vector<RevEngPoint*> > out_groups;
      bool changed = regions_[ki]->segmentByAdjSfContext(mainaxis_, min_point_in,
							 min_point_region_,
							 approx_tol_,
							 angtol, adj_groups,
							 out_groups);
      if (regions_[ki]->numPoints() == 0)
	{
	  vector<RevEngRegion*> adjacent;
	  regions_[ki]->getAdjacentRegions(adjacent);
	  for (size_t kj=0; kj<adjacent.size(); ++kj)
	    adjacent[kj]->removeAdjacentRegion(regions_[ki].get());

	  regions_.erase(regions_.begin()+ki);
	  --ki;
	}
      
      if (out_groups.size() > 0)
	{
	  vector<HedgeSurface*> dummy;
	  surfaceExtractOutput((int)ki, out_groups, dummy);
	}

    }
}

//===========================================================================
void RevEng::surfaceCreation(int pass)
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  int min_point_in = 50; //10; //20;
  adjustPointRegions(min_point_in);
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
  
  // First pass. Recognize elementary surfaces
  recognizeSurfaces(min_point_in, pass);

#ifdef DEBUG
  checkConsistence("Regions9");

   if (regions_.size() > 0)
    {
      std::cout << "Regions9" << std::endl;
      std::ofstream of("regions9.g2");
      std::ofstream ofm("mid_regions9.g2");
      std::ofstream ofs("small_regions9.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions9_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf9.g2");
	  writeRegionWithSurf(of);
	}
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges9.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
     }

   std::cout << "Merge adjacent regions, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
#endif
   double angtol = 5.0*anglim_;
   for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
	{
	  vector<RevEngRegion*> grown_regions;
	  vector<HedgeSurface*> adj_surfs;
	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
					     grown_regions, adj_surfs);
	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);

	}
    }
#ifdef DEBUG
  checkConsistence("Regions10");

   std::cout << "Finished merge adjacent regions, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions10" << std::endl;
      std::ofstream of("regions10.g2");
      std::ofstream ofm("mid_regions10.g2");
      std::ofstream ofs("small_regions10.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions10_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
        if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf10.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges10.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG_CHECK
  for (int ka=0; ka<(int)regions_.size(); ++ka)
    {
      std::set<RevEngPoint*> tmpset(regions_[ka]->pointsBegin(), regions_[ka]->pointsEnd());
      if ((int)tmpset.size() != regions_[ka]->numPoints())
	std::cout << "Point number mismatch, pre grow. " << ka << " " << tmpset.size() << " " << regions_[ka]->numPoints() << std::endl;

      bool con = regions_[ka]->isConnected();
      if (!con)
	std::cout << "Disconnected region, ka= " << ka << " " << regions_[ka].get() << std::endl;
    }
#endif

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface() && regions_[ki]->getSurfaceFlag() < ACCURACY_POOR)
	{
#ifdef DEBUG_GROW
      std::cout << "ki=" << ki << ", nmb reg: " << regions_.size() << ", nmb surf: " << surfaces_.size() << std::endl;
#endif
      growSurface(ki, pass);
#ifdef DEBUG_CHECK
  if (!regions_[ki]->isConnected())
    std::cout << "Disconnected region (grow), ki= " << ki << " " << regions_[ki].get() << std::endl;
#endif
	}
      // for (int kh=0; kh<(int)surfaces_.size(); ++kh)
      // 	{
      // 	  int numreg = surfaces_[kh]->numRegions();
      // 	  for (int ka=0; ka<numreg; ++ka)
      // 	    {
      // 	      RevEngRegion *reg = surfaces_[kh]->getRegion(ka);
      // 	      size_t kr;
      // 	      for (kr=0; kr<regions_.size(); ++kr)
      // 		if (reg == regions_[kr].get())
      // 		  break;
      // 	      if (kr == regions_.size())
      // 		std::cout << "Region4, surface 1. Obsolete region pointer, ki=" << ki << ", kh=" << kh << ". Region: " << reg << ", surface: " << surfaces_[kh].get() << std::endl;
      // 	    }
      // 	}
    }
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
  checkConsistence("Regions11");

   std::cout << "Finished grow with surf, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11" << std::endl;
      std::ofstream of("regions11.g2");
      std::ofstream ofm("mid_regions11.g2");
      std::ofstream ofs("small_regions11.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions11_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
        if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif
}

//===========================================================================
void RevEng::manageBlends1()
//===========================================================================
{
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  // Add possible missing edges
  recognizeEdges();
  
#ifdef DEBUG
  if (edges_.size() > 0)
    {
      std::ofstream ofe("edges10.g2");
      for (size_t kr=0; kr<edges_.size(); ++kr)
	{
	  vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	  for (size_t kh=0; kh<cvs.size(); ++kh)
	    {
	      cvs[kh]->writeStandardHeader(ofe);
	      cvs[kh]->write(ofe);
	    }
	  int num_blend = edges_[kr]->numBlendRegs();
	  for (int ka=0; ka<num_blend; ++ka)
	    edges_[kr]->getBlendReg(ka)->writeRegionPoints(ofe);
	}
    }
#endif
  
#ifdef DEBUG
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      int num = regions_[ki]->numPoints();
      for (int ka=0; ka<num; ++ka)
	if (regions_[ki]->getPoint(ka)->region() != regions_[ki].get())
	  std::cout << "Inconsistent region pointers, pre createBlendSurface: " << ki << " " << ka << std::endl;
    }
#endif
  
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
   std::cout << "Create blends" << std::endl;
#endif
   for (size_t ki=0; ki<edges_.size(); ++ki)
     bool found = createBlendSurface((int)ki);

   
#ifdef DEBUG
    std::cout << "Finished create blends, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_2" << std::endl;
      std::ofstream of("regions11_2.g2");
      std::ofstream ofm("mid_regions11_2.g2");
      std::ofstream ofs("small_regions11_2.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_2.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_2.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif

   // Update blend radii
   equalizeBlendRadii();
   
#ifdef DEBUG
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3_1" << std::endl;
      std::ofstream of("regions11_3_1.g2");
      std::ofstream ofm("mid_regions11_3_1.g2");
      std::ofstream ofs("small_regions11_3_1.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_4.g2");
	  writeRegionWithSurf(of);
	}
    }
#endif
   int stop_break = 1;
}

//===========================================================================
void RevEng::manageBlends2()
//===========================================================================
{
#ifdef DEBUG
   std::cout << "Set blend boundaries" << std::endl;
#endif
   for (size_t ki=0; ki<surfaces_.size(); ++ki)
     {
       int nreg = surfaces_[ki]->numRegions();
       if (nreg != 1)
	 continue;  // Not an issue currently
       RevEngRegion *reg = surfaces_[ki]->getRegion(0);
       if (reg->hasBlendEdge())
	 setBlendBoundaries(reg);
     }
   
#ifdef DEBUG
    std::cout << "Finished set blend boundaries, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3" << std::endl;
      std::ofstream of("regions11_3.g2");
      std::ofstream ofm("mid_regions11_3.g2");
      std::ofstream ofs("small_regions11_3.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_3.g2");
	  writeRegionWithSurf(of);
	}

      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_3.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_3.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif

   for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if (regions_[ka]->hasSurface() && regions_[ka]->hasBlendEdge())
   	 growBlendSurface(ka);
     }
   
#ifdef DEBUG
    std::cout << "Finished grow blend surface, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3_2" << std::endl;
      std::ofstream of("regions11_3_2.g2");
      std::ofstream ofm("mid_regions11_3_2.g2");
      std::ofstream ofs("small_regions11_3_2.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_3_2.g2");
	  writeRegionWithSurf(of);
	}

      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_3_2.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_3_2.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif

    for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if (regions_[ka]->hasSurface() && regions_[ka]->numTrimEdges() > 0 &&
	   (!regions_[ka]->hasBlendEdge()))
    	 {
    	   vector<vector<RevEngPoint*> > added_groups;
    	   vector<HedgeSurface*> dummy_surfs;
    	   double tol = 1.5*approx_tol_;
    	   double angtol = 5.0*anglim_;
    	   regions_[ka]->removeLowAccuracyPoints(min_point_region_, 
    						 tol, angtol, added_groups);
    	   if (added_groups.size() > 0)
    	     surfaceExtractOutput(ka, added_groups, dummy_surfs);
    	 }
     }
    
    for (int ka=0; ka<(int)regions_.size(); ++ka)
     {
       if (regions_[ka]->hasSurface() && regions_[ka]->numTrimEdges() > 0)
	 growMasterSurface(ka);
     }
   
#ifdef DEBUG
    std::cout << "Finished grow blend surface, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_3_3" << std::endl;
      std::ofstream of("regions11_3_3.g2");
      std::ofstream ofm("mid_regions11_3_3.g2");
      std::ofstream ofs("small_regions11_3_3.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_3_3.g2");
	  writeRegionWithSurf(of);
	}

      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_3_3.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges11_3_3.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif

#ifdef DEBUG_BLEND
   std::ofstream ofbb("blend_branch.g2");
   vector<RevEngPoint*> bbpts;
   for (size_t kr=0; kr<regions_.size(); ++kr)
     {
       if (!regions_[kr]->hasBlendEdge())
	 continue;
       vector<RevEngPoint*> currbb = regions_[kr]->extractBranchPoints();
       if (currbb.size() > 0)
	 bbpts.insert(bbpts.end(), currbb.begin(), currbb.end());
     }
   if (bbpts.size() > 0)
     {
       ofbb << "400 1 0 4 0 0 0 255" << std::endl;
       ofbb << bbpts.size() << std::endl;
       for (size_t kr=0; kr<bbpts.size(); ++kr)
	 ofbb << bbpts[kr]->getPoint() << std::endl;
     }
#endif

  // TESTING. Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

#ifdef DEBUG
   std::cout << "Torus corners" << std::endl;
#endif

   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       if (regions_[ki]->toBeRemoved())
	 continue;
       
      if (!regions_[ki]->hasSurface())
	continue;

      if (!regions_[ki]->hasBlendEdge())
	continue;

      bool done = defineTorusCorner(ki);
#ifdef DEBUG_BLEND
      std::ofstream ofsfs("curr_sfs.g2");
      for (size_t kj=0; kj<regions_.size(); ++kj)
	{
	  if (regions_[kj]->hasSurface())
	    {
	      shared_ptr<ParamSurface> curr_sf = regions_[kj]->getSurface(0)->surface();
	      curr_sf->writeStandardHeader(ofsfs);
	      curr_sf->write(ofsfs);
	    }
	}
      int stop_tor = 1;
#endif
     }
   vector<RevEngRegion*> removereg;
   vector<HedgeSurface*> removehedge;
   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       if (regions_[ki]->toBeRemoved())
	 {
	   removereg.push_back(regions_[ki].get());
	   int num_sfs = regions_[ki]->numSurface();
	   for (int ka=0; ka<num_sfs; ++ka)
	     removehedge.push_back(regions_[ki]->getSurface(ka));
	 }
     }
   if (removereg.size() > 0)
     {
      int dummy_ix = 0;
      updateRegionsAndSurfaces(dummy_ix, removereg, removehedge);
     }
   
#ifdef DEBUG
    std::cout << "Finished torus corners, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_4" << std::endl;
      std::ofstream of("regions11_4.g2");
      std::ofstream ofm("mid_regions11_4.g2");
      std::ofstream ofs("small_regions11_4.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_4.g2");
	  writeRegionWithSurf(of);
	}
      
      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_4.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
    }
#endif
   // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       // Check for possible missing corners
       if (regions_[ki]->toBeRemoved())
	 continue;
      if (regions_[ki]->hasSurface())
	continue;

      // Collect adjacent blend regions
      vector<RevEngRegion*> adj_blends;
      regions_[ki]->getAdjacentBlends(adj_blends);
#ifdef DEBUG_BLEND
      if (adj_blends.size() > 1)
	{
	  std::ofstream ofmb("candidate_blend_corner.g2");
	  regions_[ki]->writeRegionPoints(ofmb);
	  for (size_t kj=0; kj<adj_blends.size(); ++kj)
	    {
	      adj_blends[kj]->writeRegionPoints(ofmb);
	      adj_blends[kj]->writeSurface(ofmb);
	    }
	}
#endif
      if (adj_blends.size() >= 3)
	{
	  bool done = defineMissingCorner(ki, adj_blends);
	}
      
      int stop_break_corner = 1;
     }
   
   removereg.clear();
   removehedge.clear();
   for (size_t ki=0; ki<regions_.size(); ++ki)
     {
       if (regions_[ki]->toBeRemoved())
	 {
	   removereg.push_back(regions_[ki].get());
	   int num_sfs = regions_[ki]->numSurface();
	   for (int ka=0; ka<num_sfs; ++ka)
	     removehedge.push_back(regions_[ki]->getSurface(ka));
	 }
     }
   if (removereg.size() > 0)
     {
      int dummy_ix = 0;
      updateRegionsAndSurfaces(dummy_ix, removereg, removehedge);
     }
   
#ifdef DEBUG
    std::cout << "Finished missing corners, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions11_5" << std::endl;
      std::ofstream of("regions11_5.g2");
      std::ofstream ofm("mid_regions11_5.g2");
      std::ofstream ofs("small_regions11_5.g2");
      writeRegionStage(of, ofm, ofs);
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf11_5.g2");
	  writeRegionWithSurf(of);
	}
      
      vector<shared_ptr<ftEdge> > trim_edgs;
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int num = regions_[kr]->numTrimEdges();
	  if (num > 0)
	    {
	      vector<shared_ptr<ftEdge> > curr_edgs = regions_[kr]->getTrimEdges();
	      trim_edgs.insert(trim_edgs.end(), curr_edgs.begin(), curr_edgs.end());
	    }
	}
      if (trim_edgs.size() > 0)
	{
	  std::ofstream oft("trim_edges11_5.g2");
	  for (size_t kr=0; kr<trim_edgs.size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp = trim_edgs[kr]->geomCurve();
	      shared_ptr<CurveOnSurface> tmp2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      shared_ptr<ParamCurve> tmp3 = tmp2->spaceCurve();
	      tmp3->writeStandardHeader(oft);
	      tmp3->write(oft);
	    }
	}
    }
#endif
      
  int stop_break_blend = 1;
 
}


//===========================================================================
void RevEng::equalizeBlendRadii()
//===========================================================================
{
  // Collect information about blend radii
  vector<double> rad;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      
      if (!regions_[ki]->hasBlendEdge())
	continue;

      shared_ptr<ParamSurface> surf = regions_[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface, ParamSurface>(surf);
      if (!elem.get())
	continue;

      double radius = elem->radius(0.0, 0.0);
      double radius2 = elem->radius2(0.0, 0.0);
      if (radius2 > 0)
	rad.push_back(radius2);
      else
	rad.push_back(radius);
    }

  if (rad.size() == 0)
    return;
  
  std::sort(rad.begin(), rad.end());
  vector<size_t> ixs;
  ixs.push_back(0);
  double tmean = rad[0];
  size_t prev = 0;
  int nn = 1;
  size_t ki=1;
  double fac = 0.5;
  double curr_mean, range;
  for (; ki<rad.size(); ++ki)
    {
      curr_mean = (tmean + rad[ki])/(double)(nn+1);
      range = rad[ki]-rad[prev];
      if (range < fac*curr_mean)
	{
	  tmean += rad[ki];
	  ++nn;
	}
      else
	{
	  ixs.push_back(ki);
	  prev = ki;
	  tmean = rad[ki];
	  nn = 1;
	}
      int stop_break0 = 1;
    }

  vector<pair<double,double> > rad_range;
  for (ki=1; ki<ixs.size(); ++ki)
    rad_range.push_back(std::make_pair(rad[ixs[ki-1]], rad[ixs[ki]-1]));
  rad_range.push_back(std::make_pair(rad[ixs[ixs.size()-1]], rad[rad.size()-1]));

  vector<double> mean_rad(rad_range.size());
  for (ki=0; ki<rad_range.size(); ++ki)
    {
      mean_rad[ki] = 0.5*(rad_range[ki].first + rad_range[ki].second);
      if (mean_rad[ki] <= 0.1)
	{
	  double low = 0.01*(double)((int)(100.0*mean_rad[ki]));
	  double high = 0.01*(double)((int)(100.0*mean_rad[ki]+1));
	  mean_rad[ki] = (mean_rad[ki]-low < high-mean_rad[ki]) ? low : high;
	}
      else if (mean_rad[ki] <= 1.0)
	{
	  double low = 0.1*(double)((int)(10.0*mean_rad[ki]));
	  double high = 0.1*(double)((int)(10.0*mean_rad[ki]+1));
	  mean_rad[ki] = (mean_rad[ki]-low < high-mean_rad[ki]) ? low : high;
	}
      else if (mean_rad[ki] <= 10.0)
	{
	  double low = (double)((int)(mean_rad[ki]));
	  double high = (double)((int)(mean_rad[ki]+1));
	  double mid = 0.5*(low+high);
	  mean_rad[ki] = (mean_rad[ki]-low < std::min(fabs(mean_rad[ki]-mid), high-mean_rad[ki])) ? 
	    low : ((high-mean_rad[ki]) < fabs(mean_rad[ki]-mid) ? high : mid);;
	}
      else
	{
	  double low = (double)((int)(mean_rad[ki]));
	  double high = (double)((int)(mean_rad[ki]+1));
	  mean_rad[ki] = (mean_rad[ki]-low < high-mean_rad[ki]) ? low : high;
	}
    }
  
  for (ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      
      if (!regions_[ki]->hasBlendEdge())
	continue;

      shared_ptr<ParamSurface> surf = regions_[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elem =
	dynamic_pointer_cast<ElementarySurface, ParamSurface>(surf);
      if (!elem.get())
	continue;

      double radius = elem->radius(0.0, 0.0);
      double radius2 = elem->radius2(0.0, 0.0);
      if (radius2 > 0.0)
	radius = radius2;

      size_t kj;
      for (kj=0; kj<rad_range.size(); ++kj)
	if (radius >= rad_range[kj].first && radius <= rad_range[kj].second)
	  {
	    updateBlendRadius(ki, mean_rad[kj]);
	    break;
	  }
    }

  // Equialize adjacent blends between the same surfaces
  for (ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
      
      if (!regions_[ki]->hasBlendEdge())
	continue;

      RevEngEdge *edge1 = regions_[ki]->getBlendEdge();
      RevEngRegion* adj1[2];
      edge1->getAdjacent(adj1[0], adj1[1]);

      
      for (size_t kj=ki+1; kj<regions_.size(); ++kj)
	{
	  if (!regions_[kj]->hasSurface())
	    continue;
      
	  if (!regions_[kj]->hasBlendEdge())
	    continue;
	  
	  RevEngEdge *edge2 = regions_[kj]->getBlendEdge();
	  RevEngRegion* adj2[2];
	  edge2->getAdjacent(adj2[0], adj2[1]);

	  if ((adj1[0] == adj2[0] || adj1[0] == adj2[1]) &&
	      (adj1[1] == adj2[0] || adj1[1] == adj2[1]))
	    {
#ifdef DEBUG_BLEND
	      std::ofstream of("adj_blend.g2");
	      regions_[ki]->writeRegionPoints(of);
	      regions_[kj]->writeRegionPoints(of);
#endif
	      equalizeAdjacent(ki, kj);
	    }
	}
  
    }
  
    int stop_break = 1;

}

//===========================================================================
void RevEng::equalizeAdjacent(size_t ix1, size_t ix2)
//===========================================================================
{
  shared_ptr<ParamSurface> surf1 = regions_[ix1]->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = regions_[ix2]->getSurface(0)->surface();
  if ((!surf1) || (!surf2))
    return;

  if (surf1->instanceType() != surf2->instanceType())
    return;

  double angtol = 5.0*anglim_;
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  Point axis1 = elem1->direction();
  Point axis2 = elem1->direction();
  double ang = axis1.angle(axis2);
  ang = std::min(ang, M_PI-ang);
  if (ang > angtol)
    return;
  Point pos1 = elem1->location();
  Point pos2 = elem2->location();
  double rad1 = elem1->radius(0.0, 0.0);
  double rad2 = elem2->radius(0.0, 0.0);
  double fac = 0.1;
  if (fabs(rad1-rad2) > fac*std::max(rad1, rad2))
    return;  // Too large difference
  double rad = 0.5*(rad1 + rad2);
  
  if (axis1*axis2 < 0.0)
    axis2 *= -1;
  Point axis = 0.5*(axis1 + axis2);
  axis.normalize();

  Point Cx;
  Point Cx1 = elem1->direction2();
  Point Cx2 = elem2->direction2();
  double ang2 = Cx1.angle(Cx2);
  if (ang2 <= angtol || M_PI-ang2 <= angtol)
    {
      if (Cx1*Cx2 < 0.0)
	Cx2 *= -1;
      Cx = 0.5*(Cx1 + Cx2);
      Point Cy = axis.cross(Cx);
      Cx = Cy.cross(axis);
      Cx.normalize();
    }
  else
    {
      int ka = -1;
      double minang = std::numeric_limits<double>::max();
      for (int kb=0; kb<3; ++kb)
	{
	  double ang3 = mainaxis_[kb].angle(axis);
	  ang3 = std::min(ang3, M_PI-ang3);
	  if (ang3 < minang)
	    {
	      minang = ang3;
	      ka = kb;
	    }
	}
      if (ka < 0)
	return;
      int kb = (ka > 0) ? ka - 1 : 2;
      Cx = mainaxis_[kb].cross(axis);
    }

  RectDomain dom1 = elem1->getParameterBounds();
  RectDomain dom2 = elem2->getParameterBounds();
  shared_ptr<ElementarySurface> upd1, upd2;
  bool cyllike = true;
  if (surf1->instanceType() == Class_Cylinder)
    {
      Point pos2_2 = ((pos2 - pos1)*axis)*axis;
      if (pos2.dist(pos2_2) > approx_tol_)
	return;

     Point pos1_2 = ((pos1 - pos2)*axis)*axis;
      if (pos1.dist(pos1_2) > approx_tol_)
	return;

      Point pos3 = 0.5*(pos1_2 + pos2_2);  // Point on updated axis
      Point pos3_1 = pos3 + ((pos1 - pos3)*axis)*axis;
      Point pos3_2 = pos3 + ((pos2 - pos3)*axis)*axis;

      double rad = 0.5*(rad1 + rad2);
      upd1 = shared_ptr<Cylinder>(new Cylinder(rad, pos3_1, axis, Cx));
      upd2 = shared_ptr<Cylinder>(new Cylinder(rad, pos3_2, axis, Cx));

      double upar1 = 0.5*(dom1.umin() + dom2.umin());  // Could be a problem if
      // the direction of Cx is changed significantly
      double upar2 = 0.5*(dom1.umax() + dom2.umax());
      upd1->setParameterBounds(upar1, dom1.vmin(), upar2, dom1.vmax());
      upd2->setParameterBounds(upar1, dom2.vmin(), upar2, dom2.vmax());
      cyllike = true;
    }
  else if (surf1->instanceType() == Class_Torus)
    {
      if (pos1.dist(pos2) > approx_tol_)
	return;

      double minrad1 = elem1->radius2(0.0, 0.0);
      double minrad2 = elem2->radius2(0.0, 0.0);
      if (fabs(minrad1-minrad2) > fac*std::max(minrad1, minrad2))
	return;

      Point centre = 0.5*(pos1 + pos2);
      double minrad = 0.5*(minrad1 + minrad2);
      upd1 = shared_ptr<Torus>(new Torus(rad, minrad, centre, axis, Cx));
      upd2 = shared_ptr<Torus>(new Torus(rad, minrad, centre, axis, Cx));
      double vpar1 = 0.5*(dom1.vmin() + dom2.vmin());
      double vpar2 = 0.5*(dom1.vmax() + dom2.vmax());
      upd1->setParameterBounds(dom1.umin(), vpar1, dom1.umax(), vpar2);
      upd2->setParameterBounds(dom2.umin(), vpar1, dom2.umax(), vpar2);
      cyllike = false;
    }
  else
    return;  // Not supported
#ifdef DEBUG_BLEND
  std::ofstream of("updated_adjacent.g2");
  upd1->writeStandardHeader(of);
  upd1->write(of);
  upd2->writeStandardHeader(of);
  upd2->write(of);
#endif

  // Parameterize points
  double maxd1, avd1;
  int nmb_in1, nmb2_in1;
  vector<RevEngPoint*> in1, out1;
  vector<pair<double,double> > dist_ang1;
  vector<double> parvals1;
  RevEngUtils::distToSurf(regions_[ix1]->pointsBegin(), regions_[ix1]->pointsEnd(),
			  upd1, approx_tol_, maxd1, avd1, nmb_in1, nmb2_in1, in1, out1,
			  parvals1, dist_ang1, angtol);
  int sf_flag1 = regions_[ix1]->defineSfFlag(0, approx_tol_, nmb_in1, nmb2_in1, avd1,
					     cyllike);
  int num_pts1 = regions_[ix1]->numPoints();
  for (int ka=0; ka<num_pts1; ++ka)
    {
      RevEngPoint *curr = regions_[ix1]->getPoint(ka);
      curr->setPar(Vector2D(parvals1[2*ka],parvals1[2*ka+1]));
      curr->setSurfaceDist(dist_ang1[ka].first, dist_ang1[ka].second);
    }
  regions_[ix1]->updateInfo(approx_tol_, angtol);
  regions_[ix1]->setSurfaceFlag(sf_flag1);

  // Replace Surface
  regions_[ix1]->getSurface(0)->replaceSurf(upd1);
  
  double maxd2, avd2;
  int nmb_in2, nmb2_in2;
  vector<RevEngPoint*> in2, out2;
  vector<pair<double,double> > dist_ang2;
  vector<double> parvals2;
  RevEngUtils::distToSurf(regions_[ix2]->pointsBegin(), regions_[ix2]->pointsEnd(),
			  upd1, approx_tol_, maxd2, avd2, nmb_in2, nmb2_in2, in2, out2,
			  parvals2, dist_ang2, angtol);
  int sf_flag2 = regions_[ix2]->defineSfFlag(0, approx_tol_, nmb_in2, nmb2_in2, avd2,
					     cyllike);
  int num_pts2 = regions_[ix2]->numPoints();
  for (int ka=0; ka<num_pts2; ++ka)
    {
      RevEngPoint *curr = regions_[ix2]->getPoint(ka);
      curr->setPar(Vector2D(parvals2[2*ka],parvals2[2*ka+1]));
      curr->setSurfaceDist(dist_ang2[ka].first, dist_ang2[ka].second);
    }
  regions_[ix2]->updateInfo(approx_tol_, angtol);
  regions_[ix2]->setSurfaceFlag(sf_flag2);

  // Replace Surface
  regions_[ix2]->getSurface(0)->replaceSurf(upd2);
  
}

//===========================================================================
void RevEng::updateBlendRadius(size_t ix, double radius)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  double eps = 1.0e-6;
  RevEngEdge *edge = regions_[ix]->getBlendEdge();
  RevEngRegion* adj[2];
  edge->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return;  // Something wrong

  bool out1 = false, out2 = false;
  edge->getOuterInfo(out1, out2);

  // Intersection curve
  vector<shared_ptr<CurveOnSurface> > cvs;
  edge->getCurve(cvs);
  vector<Point> der(2);
  cvs[0]->point(der, 0.5*(cvs[0]->startparam()+cvs[0]->endparam()), 1);

  shared_ptr<ParamSurface> surf = regions_[ix]->getSurface(0)->surface();
  shared_ptr<Cylinder> init_cyl =
    dynamic_pointer_cast<Cylinder, ParamSurface>(surf);
  shared_ptr<Torus> init_tor =
    dynamic_pointer_cast<Torus, ParamSurface>(surf);
  if ((!init_cyl.get()) && (!init_tor.get()))
    return;
  shared_ptr<ParamSurface> surf1 = adj[0]->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adj[1]->getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  shared_ptr<ElementarySurface> upd_blend;
  double ylen;
  double init_rad;
  Point dir1 = elem1->direction();
  Point norm1 = adj[0]->getMeanNormalTriang();
  if (elem1->instanceType() == Class_Plane && dir1*norm1 < 0.0)
    dir1 *= -1;
  Point dir2 = elem2->direction();
  Point norm2 = adj[1]->getMeanNormalTriang();
  if (elem2->instanceType() == Class_Plane && dir2*norm2 < 0.0)
    dir2 *= -1;
  if (init_cyl.get())
    {
      if (elem1->instanceType() != Class_Plane || elem2->instanceType() != Class_Plane)
	return; 

      Point lin1, lin2;
      Point dir1_2 = dir1, dir2_2 = dir2;
      if (elem1->instanceType() == Class_Plane)
	lin1 = der[1].cross(dir1);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf1->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf1->normal(dir1_2, clo_u, clo_v);
	  lin1 = der[1].cross(dir1_2);
	}
      
      if (elem2->instanceType() == Class_Plane)
	lin2 = der[1].cross(dir2);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf2->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf2->normal(dir2_2, clo_u, clo_v);
	  lin2 = der[1].cross(dir2_2);
	}
     
      // Create new cylinder
      Point axis = init_cyl->direction();
      Point Cx = init_cyl->direction2();
      init_rad = init_cyl->getRadius();
      lin1.normalize();
      if (lin1*dir2_2 < 0.0)
	lin1 *= -1;
      lin2.normalize();
      if (lin2*dir1_2 < 0.0)
	lin2 *= -1;

      Point vec = lin1 + lin2;
      vec.normalize();
      double alpha = lin1.angle(lin2);
      double fac = 1.0/sin(0.5*alpha);
      int sgn = (out1 && out2) ? 1 : -1;
      Point loc = der[0] + sgn*fac*radius*vec;
      double xlen = der[0].dist(loc);
      ylen = sqrt(xlen*xlen - radius*radius);

      upd_blend = shared_ptr<Cylinder>(new Cylinder(radius, loc, axis, Cx));
    }
  else if (init_tor.get())
    {
      // int plane_ix = 0;
      // if (elem2->instanceType() == Class_Plane)
      // 	{
      // 	  std::swap(elem1, elem2);
      // 	  plane_ix = 1;
      // 	}
      // if (elem1->instanceType() != Class_Plane)
      // 	return; 
      // if (elem2->instanceType() != Class_Cylinder &&
      // 	  elem2->instanceType() != Class_Cone)
      // 	return;
      
      // init_rad = init_tor->radius2(0.0, 0.0);

      // bool rot_out = (plane_ix == 0) ? out2 : out1;
      // bool plane_out = (plane_ix == 0) ? out1 : out2;
      // int sgn1 = plane_out ? -1 : 1;
      // int sgn2 = rot_out ? -1 : 1;
      // Point normal0 = elem1->direction();
      // Point norm1 = adj[plane_ix]->getMeanNormalTriang();
      // if (normal0*norm1 < 0.0)
      // 	sgn1 *= -1;

      bool plane1 = (elem1->instanceType() == Class_Plane);
      bool plane2 = (elem2->instanceType() == Class_Plane);
      shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
      shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
      shared_ptr<ElementarySurface> rotational;
      if (plane1)
	rotational = elem2;
      else if (plane2)
	rotational = elem1;
      else if (cone1.get())
	rotational = elem2;  // elem2 is expected to be a cylinder
      else if (cone2.get())
	rotational = elem1;
      double alpha1 = 0.0, alpha2 = 0.0;
      if (cone1.get())
	alpha1 = cone1->getConeAngle();
      if (cone2.get())
	alpha2 = cone2->getConeAngle();
      double alpha = fabs(alpha1) + fabs(alpha2);
      double beta = (plane1 || plane2) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
      double phi2 = 0.5*(M_PI - beta);
      Point axis = rotational->direction();
      Point loc = rotational->location();
      double rd = (der[0]-loc)*axis;
      Point centr = loc + rd*axis;
      double hh = init_tor->location().dist(centr);
      init_rad = init_tor->radius2(0.0, 0.0);
      double d2 = hh/sin(phi2) - init_rad;
      int sgn = 1;
      if ((elem1->instanceType() == Class_Plane && elem1->direction()*dir1 < 0.0) ||
	  (elem2->instanceType() == Class_Plane && elem2->direction()*dir2 < 0.0))
	sgn = -1;
      double Rrad;
      Point centre, normal, Cx;
      getTorusParameters(elem1, elem2, der[0], radius, d2, out1, out2, sgn,
			 Rrad, centre, normal, Cx);
     
      upd_blend = shared_ptr<Torus>(new Torus(Rrad, radius, centre, normal, Cx));

      // bool setbounds = false;
      // if ((plane1 && out2) || (plane2 && out1) ||
      // 	  (plane1 && false && plane2 == false &&
      // 	   ((cone1.get() && out1 == false) || (cone2.get() && out2 == false))))
      // 	setbounds = true;
	
      // if (setbounds)
      // 	{
	  RectDomain dom = init_tor->getParameterBounds();
	  upd_blend->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
					dom.umax(), dom.vmax()-M_PI);
	// }
      double xlen = der[0].dist(centre);
      ylen = fabs(xlen - Rrad);
     }
  
  // RectDomain dom = init_cyl->getParameterBounds();
  // cyl->setParameterBounds(dom.umin(), dom.vmin(), dom.umax(), dom.vmax());
#ifdef DEBUG_BLEND
  std::ofstream of("blend2.g2");
  regions_[ix]->writeRegionPoints(of);
  surf->writeStandardHeader(of);
  surf->write(of);

  adj[0]->writeRegionPoints(of);
  adj[1]->writeRegionPoints(of);
  upd_blend->writeStandardHeader(of);
  upd_blend->write(of);
  RectDomain dom2 = upd_blend->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = upd_blend->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of);
  tmp_cvs[0]->write(of);
#endif

  if (radius < init_rad)
    {
      // Move outside points from the blend surface to the adjacent
      // surfaces. If the new radius is larger than the initial one,
      // moving of points will be performed by the next function
      vector<RevEngPoint*> points(regions_[ix]->pointsBegin(), regions_[ix]->pointsEnd());
      vector<vector<RevEngPoint*> > out_pts(2);
      adj[0]->sortBlendPoints(points, cvs, ylen, adj[1], out_pts[0], out_pts[1]);
      vector<RevEngPoint*> all_out;
      if (out_pts[0].size() > 0)
	all_out.insert(all_out.end(), out_pts[0].begin(), out_pts[0].end());
      if (out_pts[1].size() > 0)
	all_out.insert(all_out.end(), out_pts[1].begin(), out_pts[1].end());
      if (all_out.size() > 0)
	regions_[ix]->removePoints(all_out);
      
      for (int ka=0; ka<2; ++ka)
	{
	  vector<RevEngPoint*> to_blend;  // Not here
	  adj[ka]->updateWithPointsInOut(to_blend, out_pts[ka], approx_tol_, angtol);
	}
    }

  // Parameterize remaining points
  double maxd, avd;
  int nmb_in, nmb2_in;
  vector<RevEngPoint*> in, out;
  vector<pair<double,double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(regions_[ix]->pointsBegin(), regions_[ix]->pointsEnd(),
			  upd_blend, approx_tol_, maxd, avd, nmb_in, nmb2_in, in, out,
			  parvals, dist_ang, angtol);
  int sf_flag = regions_[ix]->defineSfFlag(0, approx_tol_, nmb_in, nmb2_in, avd, true);
  int num_pts = regions_[ix]->numPoints();
  for (int ka=0; ka<num_pts; ++ka)
    {
      RevEngPoint *curr = regions_[ix]->getPoint(ka);
      curr->setPar(Vector2D(parvals[2*ka],parvals[2*ka+1]));
      curr->setSurfaceDist(dist_ang[ka].first, dist_ang[ka].second);
    }
  regions_[ix]->setSurfaceFlag(sf_flag);

  // Replace Surface
  regions_[ix]->getSurface(0)->replaceSurf(upd_blend);
  

  int stop_break = 1;
}


//===========================================================================
bool RevEng::defineTorusCorner(size_t ix)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  double angtol = 5.0*anglim_;
  double blendfac = 2.0;
  
  shared_ptr<ParamSurface> surf1 = regions_[ix]->getSurface(0)->surface();
  if (surf1->instanceType() != Class_Cylinder)
    return false;

  if (!regions_[ix]->hasBlendEdge())
    return false;

  // Bound cylinder in the length direction
  double diag = bbox_.low().dist(bbox_.high());
  if (!surf1->isBounded())
    {
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
      cyl->setParamBoundsV(-0.5*diag, 0.5*diag);
    }

  RevEngEdge* blend_edge = regions_[ix]->getBlendEdge();
  RevEngRegion *adj1=0, *adj2=0;
  blend_edge->getAdjacent(adj1, adj2);
  
  vector<RevEngEdge*> revedg1 = regions_[ix]->getAllRevEdges();
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  Point dir1 = elem1->direction();

#ifdef DEBUG_BLEND
  std::ofstream of0("pot_tor_adj.g2");
  surf1->writeStandardHeader(of0);
  surf1->write(of0);
#endif  
  // Looking for a plane
  bool found_torus = false;
  for (size_t kj=0; kj<regions_.size(); ++kj)
    {
      if (kj == ix)
	continue;
      if (regions_[kj]->toBeRemoved())
	continue;
      if (!regions_[kj]->hasSurface())
	continue;
      if (regions_[kj].get() == adj1 || regions_[kj].get() == adj2)
	continue;

      shared_ptr<ParamSurface> surf2 = regions_[kj]->getSurface(0)->surface();
      if (surf2->instanceType() != Class_Plane)
	continue;

      if (!(regions_[ix]->isAdjacent(regions_[kj].get())
	    || regions_[ix]->isNextToAdjacent(regions_[kj].get())))
	continue;

      // Check if the plane and the cylinder already has a common RevEngEdge
      vector<RevEngEdge*> revedg2 = regions_[kj]->getAllRevEdges();
      size_t kr, kh;
      for (kr=0; kr<revedg1.size(); ++kr)
	{
	  for (kh=0; kh<revedg2.size(); ++kh)
	    if (revedg1[kr] == revedg2[kh])
	      break;
	  if (kh < revedg2.size())
	    break;
	}
      if (kr < revedg1.size())
	continue;
      
#ifdef DEBUG_BLEND
      std::ofstream of1("torus_blend1.g2");
      regions_[ix]->writeRegionPoints(of1);
      regions_[ix]->writeSurface(of1);
      regions_[kj]->writeRegionPoints(of1);
      regions_[kj]->writeSurface(of1);
#endif
      
      shared_ptr<ElementarySurface> elem2 =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
      if (!elem2->isBounded())
	elem2->setParameterBounds(-0.5*diag, -0.5*diag, 0.5*diag, 0.5*diag);
      Point dir2 = elem2->direction();
      double ang = dir1.angle(dir2);
      ang = fabs(pihalf - ang);
      if (ang > blendfac*angtol)
	{
	  double usz, vsz;
	  elem1->estimateSfSize(usz, vsz);
	  double lenlim = 0.9*usz;

	  vector<shared_ptr<RevEngEdge> > edges =
	    defineEdgesBetween(ix, elem1, dir1, kj, elem2, dir2, lenlim, false);
	  if (edges.size() > 0)
	    {
	      size_t ix2 = edges_.size();
	      edges_.insert(edges_.end(), edges.begin(), edges.end());
	      bool found = createTorusBlend(ix2);
	      // if (!found)
	      // 	found = createBlendSurface((int)ix2);
	      if (found)
		found_torus = true;
	      else
		edges_.erase(edges_.begin()+ix2, edges_.end());
	    }
	}
    }

  return found_torus;
}

//===========================================================================
bool RevEng::defineMissingCorner(size_t ix, vector<RevEngRegion*>& adj_blends)
//===========================================================================
{
  bool found_torus = false;
  double pihalf = 0.5*M_PI;
  double angtol = 5.0*anglim_;
  double diag = bbox_.low().dist(bbox_.high());

  // Get candidate opposite regions
  vector<RevEngRegion*> opposite_reg;
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      RevEngEdge* edg1 = adj_blends[ki]->getBlendEdge();
      if (!edg1)
	continue;
      RevEngRegion *adj1, *adj2;
      edg1->getAdjacent(adj1, adj2);
      if ((!adj1) || (!adj2))
	continue;
      for (size_t kj=ki+1; kj<adj_blends.size(); ++kj)
	{
	  RevEngEdge* edg2 = adj_blends[kj]->getBlendEdge();
	  if (!edg2)
	    continue;
	  RevEngRegion *adj3, *adj4;
	  edg2->getAdjacent(adj3, adj4);
	  if ((!adj3) || (!adj4))
	    continue;
	  if (adj1 == adj3 || adj1 == adj4)
	    opposite_reg.push_back(adj1);
	  if (adj2 == adj3 || adj2 == adj4)
	    opposite_reg.push_back(adj2);
	}
    }
#ifdef DEBUG_BLEND
  if (opposite_reg.size() > 0)
    {
      std::ofstream of1("opposite_reg.g2");
      for (size_t ki=0; ki<opposite_reg.size(); ++ki)
	opposite_reg[ki]->writeRegionPoints(of1);
    }
#endif

  double blendfac = 2.0;
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      RevEngEdge* edg1 = adj_blends[ki]->getBlendEdge();
      RevEngRegion *adj1, *adj2;
      edg1->getAdjacent(adj1, adj2);
      if ((!adj1) || (!adj2))
	continue;
      shared_ptr<ParamSurface> surf1 = adj_blends[ki]->getSurface(0)->surface();
      if (surf1->instanceType() != Class_Cylinder)
	continue;
      if (!surf1->isBounded())
	{
	  shared_ptr<Cylinder> cyl =
	    dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
	  cyl->setParamBoundsV(-0.5*diag, 0.5*diag);
	}
      shared_ptr<ElementarySurface> elem1 =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
      Point dir1 = elem1->direction();
      
      for (size_t kj=0; kj<opposite_reg.size(); ++kj)
	{
	  if (!opposite_reg[kj]->hasSurface())
	    continue;
	  if (opposite_reg[kj] == adj1 || opposite_reg[kj] == adj2)
	    continue;
	  shared_ptr<ParamSurface> surf2 =
	    opposite_reg[kj]->getSurface(0)->surface();
	  if (surf2->instanceType() != Class_Plane)
	    continue;
	  shared_ptr<ElementarySurface> elem2 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
	  if (!elem2->isBounded())
	    elem2->setParameterBounds(-0.5*diag, -0.5*diag, 0.5*diag, 0.5*diag);
	  Point dir2 = elem2->direction();
	  double ang = dir1.angle(dir2);
	  ang = fabs(pihalf - ang);
	  if (ang > blendfac*angtol)
	    {
	      double usz, vsz;
	      elem1->estimateSfSize(usz, vsz);
	      double lenlim = 0.9*usz;

	      size_t kr, kh;
	      for (kr=0; kr<regions_.size(); ++kr)
		if (regions_[kr].get() == adj_blends[ki])
		  break;
	      if (kr == regions_.size())
		continue;
	      for (kh=0; kh<regions_.size(); ++kh)
		if (regions_[kh].get() == opposite_reg[kj])
		  break;
	      if (kh == regions_.size())
		continue;
	      vector<shared_ptr<RevEngEdge> > edges =
		defineEdgesBetween(kr, elem1, dir1, kh, elem2, dir2, lenlim, false);
	      if (edges.size() > 0)
		{
		  size_t ix2 = edges_.size();
		  edges_.insert(edges_.end(), edges.begin(), edges.end());
		  bool found = createTorusBlend(ix2);
		  if (found)
		    found_torus = true;
		  else
		    edges_.erase(edges_.begin()+ix2, edges_.end());
		}
	    }
	  
	}
    }
  
  return found_torus;
}

//===========================================================================
bool getAdjacentToTorus(RevEngEdge* edge, vector<RevEngEdge*>& rev_edgs,
			double tol, RevEngEdge*& adj_edg1,
			RevEngEdge*& adj_edg2, double& rad1, double& rad2)
//===========================================================================
{
  adj_edg1 = adj_edg2 = 0;
  rad1 = rad2 = -1.0;
  
  Point pos1, pos2;
  edge->getCrvEndPoints(pos1, pos2);

  double mindist1 = std::numeric_limits<double>::max();
  double mindist2 = std::numeric_limits<double>::max();
  int min_ix1 = -1, min_ix2 = -1;
  for (size_t ki=0; ki<rev_edgs.size(); ++ki)
    {
      RevEngRegion *reg = rev_edgs[ki]->getBlendRegSurf();
      if (!reg)
	continue;

      double tpar1, tpar2, dist1, dist2;
      Point close1, close2;
      rev_edgs[ki]->closestPoint(pos1, tpar1, close1, dist1);
      rev_edgs[ki]->closestPoint(pos2, tpar2, close2, dist2);
      if (dist1 < mindist1)
	{
	  mindist1 = dist1;
	  min_ix1 = (int)ki;
	}
       if (dist2 < mindist2)
	{
	  mindist2 = dist2;
	  min_ix2 = (int)ki;
	}
      int stop_break = 1;
    }
  if (min_ix1 < 0 && min_ix2 < 0)
    return false;
  if (mindist1 > tol && mindist2 > tol)
    return false;   // Might need to tune tolerance

  if (min_ix1 >= 0 && mindist1 <= tol)
    {
      adj_edg1 = rev_edgs[min_ix1];
      RevEngRegion *reg = adj_edg1->getBlendRegSurf();
      if (reg->hasSurface())
	{
	  shared_ptr<ParamSurface> bsurf = reg->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> belem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(bsurf);
	  if (bsurf->instanceType() == Class_Cylinder /*||
							bsurf->instanceType() == Class_Cone*/)
	    rad1 = belem->radius(0.0, 0.0);
	  else if (bsurf->instanceType() == Class_Torus)
	    rad1 = belem->radius2(0.0, 0.0);
	}
    }
  if (min_ix2 >= 0 && mindist2 <= tol)
    {
      adj_edg2 = rev_edgs[min_ix2];
      RevEngRegion *reg = adj_edg2->getBlendRegSurf();
      if (reg->hasSurface())
	{
	  shared_ptr<ParamSurface> bsurf = reg->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> belem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(bsurf);
	  if (bsurf->instanceType() == Class_Cylinder /*||
							bsurf->instanceType() == Class_Cone*/)
	    rad2 = belem->radius(0.0, 0.0);
	  else if (bsurf->instanceType() == Class_Torus)
	    rad2 = belem->radius2(0.0, 0.0);
	}
    }

  if (rad1 < 0 && rad2 < 0)
    return false;

  return true;
}

//===========================================================================
bool RevEng::createTorusBlend(size_t ix)
//===========================================================================
{
  double eps = 1.0e-6;
  double angtol = 5.0*anglim_;
  double pihalf = 0.5*M_PI;
  double tol2 = 2*approx_tol_;  // Due to possible inaccuracies

  // Adjacent regions
  RevEngRegion* adj[2];
  edges_[ix]->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return false;  // Something wrong
  
  vector<shared_ptr<CurveOnSurface> > intcv1, intcv2;
  edges_[ix]->getCurve(intcv1, true);
  edges_[ix]->getCurve(intcv2, false);
  
  bool out1 = false, out2 = false;
  edges_[ix]->getOuterInfo(out1, out2);
  
  shared_ptr<ParamSurface> surf[2];
  surf[0] = adj[0]->getSurface(0)->surface();
  surf[1] = adj[1]->getSurface(0)->surface();
  ClassType type[2];
  type[0] = surf[0]->instanceType();
  type[1] = surf[1]->instanceType();
  int jx1 = (type[0] == Class_Plane) ? 0 :
    ((type[1] == Class_Plane) ? 1 : -1);
  if (jx1 == -1)
    return false;
  int jx2 = 1 - jx1;
  if (type[jx2] != Class_Cylinder)
    return false;
  bool outp = (jx1 == 0) ? out1 : out2;
  bool outr = (jx1 == 0) ? out2 : out1;
  shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf[jx2]);
  if (!cyl.get())
    return false;  // Should not happen
  double radius1 = cyl->getRadius();

  // Get adjacent blend edges and corresponding radii
  vector<RevEngEdge*> rev_edgs = adj[jx1]->getAllRevEdges();
  if (rev_edgs.size() == 0)
    return false;

  double rad1 = -1.0, rad2 = -1.0;
  RevEngEdge *revedg1, *revedg2;
  bool OK = getAdjacentToTorus(edges_[ix].get(), rev_edgs, tol2, 
			       revedg1, revedg2, rad1, rad2);
  if (!OK)
    return false;
  if ((!revedg1) || (!revedg2))
    return false;
  bool possible_suitcase =
    (fabs(rad1-rad2) > std::min(fabs(radius1-rad1), fabs(radius1-rad2)));

  double radius2 = (rad1 > 0.0 && rad2 > 0.0) ? 0.5*(rad1 + rad2) :
    ((rad1 > 0.0) ? rad1 : rad2);  // Should do extra checking if rad1 very
  // different from rad2 and both are larger than zero

  Point loc1 = cyl->getLocation();
  Point axis1 = cyl->getAxis();
  Point Cx = cyl->direction2();

  // The torus location coincides with the intersection point between the
  // cylinder axis and the plane
  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf[jx1]);
  Point dirp = plane->direction();
  Point avnorm = adj[jx1]->getMeanNormal();
  if (dirp*avnorm < 0.0)
    dirp *= -1;
  if (axis1*dirp < 0.0)
    axis1 *= -1;
  double alpha = axis1.angle(dirp);
  Point locp = plane->location();
  Point loct0 = loc1 + ((locp - loc1)*axis1)*axis1;
  double d2 = locp.dist(loct0);
  double xlen = d2*atan(alpha);
  int sgn1 = ((locp - loct0)*dirp < 0.0) ? -1 : 1;
  loct0 += sgn1*xlen*dirp;
  int sgn2 = outp ? 1 : -1;
  Point loct = loct0 + sgn2*radius2*dirp;
  int sgn3 = outr ? -1 : 1;
  double radiust = radius1 + sgn3*radius2;
  if (radius2 >= radiust)
    return false;  // Not expecting degenerate torus

  shared_ptr<Torus> torus(new Torus(radiust, radius2, loct, dirp, Cx));

  // Regions in blend area
  vector<RevEngRegion*> blend_regs;
  edges_[ix]->getAllBlendRegs(blend_regs);
#ifdef DEBUG_BLEND
  std::ofstream of("torus_blend2.g2");
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_regs[ki]->writeRegionPoints(of);
#endif
  
  // Parameterize on torus
  vector<RevEngPoint*> blend_pts;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_pts.insert(blend_pts.end(), blend_regs[ki]->pointsBegin(),
		     blend_regs[ki]->pointsEnd());
  double maxd, avd;
  int num_in, num2_in;
  vector<double> parvals;
  vector<pair<double,double> > dist_ang;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(blend_pts.begin(), blend_pts.end(), torus,
			  approx_tol_, maxd, avd, num_in, num2_in, inpt, outpt,
			  parvals, dist_ang, angtol);
  
  RectDomain dom = cyl->getParameterBounds();
  if (sgn3 < 0)
    torus->setParameterBounds(dom.umin(), 0.0, dom.umax(), pihalf);
  else
    torus->setParameterBounds(dom.umin(), M_PI, dom.umax(), 1.5*M_PI);
  // What if the cylinder and the plane are not perpendicular? It is not a
  // torus, but it is anyway not handled
  
#ifdef DEBUG_BLEND
  torus->writeStandardHeader(of);
  torus->write(of);
  vector<shared_ptr<ParamCurve> > c_cvs = torus->constParamCurves(0.0, true);
  c_cvs[0]->writeStandardHeader(of);
  c_cvs[0]->write(of);
#endif

  // Define longitudinal boundary edges of torus
  RectDomain tordom = torus->getParameterBounds();
  double torlim[4];
  torlim[0] = tordom.umin();
  torlim[1] = tordom.umax();
  torlim[2] = tordom.vmin();
  torlim[3] = tordom.vmax();
  vector<shared_ptr<CurveOnSurface> > torbound(4);
  for (int ka=0; ka<2; ++ka)
    {
      shared_ptr<Circle> circle = torus->getMajorCircle(torlim[2+ka]);
      Point mid(torlim[0], torlim[2+ka]);
      Point vec = Point(torlim[1],torlim[2+ka]) - Point(torlim[0],torlim[2+ka]);
      vec.normalize();
      shared_ptr<Line> line(new Line(mid, vec));
      line->setParamBounds(0.0, torlim[1]-torlim[0]);
      torbound[2+ka] =
	shared_ptr<CurveOnSurface>(new CurveOnSurface(torus, line, circle, false,
						      3, 2, torlim[2+ka], 2+ka, true));
    }

  // Restrict adjacent cylinders
  // The cylinder and the edges does not necessarily have the same
  // parameterization
  Point pt1, pt2;
  double midpar = 0.5*(tordom.vmin() + tordom.vmax());
  torus->point(pt1, tordom.umin(), midpar);
  torus->point(pt2, tordom.umax(), midpar);
  vector<RevEngRegion*> adj_reg(2, 0);
  bool upper2[2];
  double lim[2];
  vector<shared_ptr<ElementarySurface> > adj_surf(2);
  vector<double> parbound(8);
  int kx[2];
  for (int ka=0; ka<2; ++ka)
    {
      if ((ka == 0 && rad1 <= 0.0) || (ka == 1 && rad2 <= 0.0))
	{
	  kx[ka] = -1;
	  continue;
	}

      RevEngRegion *reg = (ka == 0) ? revedg1->getBlendRegSurf() :
	revedg2->getBlendRegSurf();
      shared_ptr<ParamSurface> bsurf = reg->getSurface(0)->surface();
      shared_ptr<ElementarySurface> belem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(bsurf);
      adj_reg[ka] = reg;
      adj_surf[ka] = belem;
      
     double reg_dom[4];
      reg->getDomain(reg_dom);
      Point pt3[2];
      int kx1 = (bsurf->instanceType() == Class_Cylinder) ? 0 : 1;
      for (int kb=0; kb<2; ++kb)
	{
	  double par[2];
	  par[kx1] = 0.5*(reg_dom[2*kx1] + reg_dom[2*kx1+1]);
	  par[1-kx1] = reg_dom[2*(1-kx1)+kb];
	  pt3[kb] = bsurf->point(par[0], par[1]);
	}
      
      double upar1, vpar1, bdist1, upar2, vpar2, bdist2;
      Point bclose1, bclose2;
      bsurf->closestPoint(pt1, upar1, vpar1, bclose1, bdist1, eps);
      bsurf->closestPoint(pt2, upar2, vpar2, bclose2, bdist2, eps);
      double uvpar;
      if (belem->instanceType() == Class_Cylinder)
	uvpar = (bdist1 <= bdist2) ? vpar1 : vpar2;
      else
	uvpar = (bdist1 <= bdist2) ? upar1 : upar2;
      Point tor_pt = (bdist1 <= bdist2) ? pt1 : pt2;
      lim[ka] = uvpar;
      kx[ka] = (bdist1 <= bdist2) ? 0 : 1;

      // Assume the part of the cylinder to remove is the smaller one
      RectDomain bdom = bsurf->containingDomain();
      parbound[4*ka] = bdom.umin();
      parbound[4*ka+1] = bdom.vmin();
      parbound[4*ka+2] = bdom.umax();
      parbound[4*ka+3] = bdom.vmax();
      if (tor_pt.dist(pt3[1]) < tor_pt.dist(pt3[0]))
	{
	  if (belem->instanceType() == Class_Cylinder)
	    parbound[4*ka+3] = uvpar;
	  else
	    parbound[4*ka+2] = uvpar;
	  upper2[ka] = true;
	}
      else 
	{
	  if (belem->instanceType() == Class_Cylinder)
	    parbound[4*ka+1] = uvpar;
	  else
	    parbound[4*ka] = uvpar;
	  upper2[ka] = false;
	}
#ifdef DEBUG_BLEND
      bsurf->writeStandardHeader(of);
      bsurf->write(of);
#endif
    }

  Point vec1 = (adj_surf[0].get()) ? adj_surf[0]->location() - plane->location() : Point(0.0, 0.0, 0.0);
  Point vec2 = (adj_surf[1].get()) ? adj_surf[1]->location() - plane->location() : Point(0.0, 0.0, 0.0);
  double sc1 = vec1*dirp;
  double sc2 = vec2*dirp;
  if (sc1*sc2 < 0.0) //kx[0] == kx[1] && kx[0]>=0)
    {
      // Need to wait with "setParameterBounds"
      vector<RevEngRegion*> blends(3);
      blends[0] = adj[jx2];
      blends[1] = revedg1->getBlendRegSurf();
      blends[2] = revedg2->getBlendRegSurf();
      return suitcaseCorner(blends, edges_[ix].get());
    }
  else if (possible_suitcase)
    return false;  // Not an expected configuration

  for (int ka=0; ka<2; ++ka)
    {
      if (adj_surf[ka].get())
	{
	  adj_surf[ka]->setParameterBounds(parbound[4*ka], parbound[4*ka+1],
					   parbound[4*ka+2], parbound[4*ka+3]);
	  adj_reg[ka]->adaptEdges();
	}
   }

  // Bound blend cylinder
  Point pt3;
  torus->point(pt3, 0.5*(dom.umin()+dom.vmin()), 0);
  double upar3, vpar3, dist3;
  Point close3;
  bool upper = false;
  cyl->closestPoint(pt3, upar3, vpar3, close3, dist3, eps);
  if (fabs(vpar3-dom.vmin()) < fabs(dom.vmax()-vpar3))
    cyl->setParamBoundsV(vpar3, dom.vmax());
  else
    {
      cyl->setParamBoundsV(dom.vmin(), vpar3);
      upper = true;
    }
  adj[jx2]->adaptEdges();
  
 #ifdef DEBUG_BLEND
  cyl->writeStandardHeader(of);
  cyl->write(of);
#endif

  shared_ptr<Circle> bcircle = cyl->getCircle(vpar3);
  shared_ptr<ElementaryCurve> bpar =
    cyl->getElementaryParamCurve(bcircle.get(), approx_tol_);
  shared_ptr<CurveOnSurface> blendbound(new CurveOnSurface(cyl, bpar, bcircle,
							   false, 3, 2, vpar3,
							   (upper) ? 3 : 2, true));
 
  // Define traverse boundary edges of torus and adjacent cylinders
  vector<shared_ptr<CurveOnSurface> > cylbound(2);
  for (int ka=0; ka<2; ++ka)
    {
      if (!adj_reg[ka])
	continue;

      bool udir = (adj_surf[ka]->instanceType() == Class_Cylinder) ? true : false;
      shared_ptr<Circle> circle = torus->getMinorCircle(torlim[ka]);
      Point mid(torlim[ka], torlim[2]);
      Point vec = Point(torlim[ka],torlim[3]) - Point(torlim[ka],torlim[2]);
      vec.normalize();
      shared_ptr<Line> line(new Line(mid, vec));
      line->setParamBounds(0.0, torlim[3]-torlim[2]);
      torbound[ka] =
	shared_ptr<CurveOnSurface>(new CurveOnSurface(torus, line, circle,
						      false, 3, 1, torlim[ka],
						      ka, true));
      vector<shared_ptr<ParamCurve> > tmpcv =
	adj_surf[ka]->constParamCurves(lim[ka], udir);
      if (tmpcv.size() != 1)
	continue;  // Should not happen
      shared_ptr<ElementaryCurve> space =
	dynamic_pointer_cast<ElementaryCurve,ParamCurve>(tmpcv[0]);
      shared_ptr<ElementaryCurve> par =
	adj_surf[ka]->getElementaryParamCurve(space.get(), approx_tol_);
      int bd;
      if (udir)
	bd = (upper2[ka]) ? 3 : 2;
      else
	bd = (upper2[ka]) ? 1 : 0;	
      cylbound[ka] =
	shared_ptr<CurveOnSurface>(new CurveOnSurface(adj_surf[ka], par, space,
						      false, 3, udir ? 1 : 2,
						      lim[ka],
						      bd, true));
    }

  // Define planar boundary curve
  shared_ptr<CurveOnSurface> planebound;
  shared_ptr<ParamCurve> pspace1 = torbound[3]->spaceCurve();
  shared_ptr<Circle> pspace = dynamic_pointer_cast<Circle,ParamCurve>(pspace1);
  shared_ptr<ElementaryCurve> ppar =
    plane->getElementaryParamCurve(pspace.get(), 10.0*tol2);  // Just to test
  planebound = shared_ptr<CurveOnSurface>(new CurveOnSurface(plane, ppar,
							     pspace, false, 1));

#ifdef DEBUG_BLEND
  std::ofstream ofc("torus_trim.g2");
  for (int ka=0; ka<4; ++ka)
    {
      if (!torbound[ka].get())
	continue;
      shared_ptr<ParamCurve> tmp = torbound[ka]->spaceCurve();
      tmp->writeStandardHeader(ofc);
      tmp->write(ofc);
    }
  for (int ka=0; ka<2; ++ka)
    {
      if (!cylbound[ka].get())
	continue;
      shared_ptr<ParamCurve> tmp = cylbound[ka]->spaceCurve();
      tmp->writeStandardHeader(ofc);
      tmp->write(ofc);
    }
  if (planebound.get())
    {
      shared_ptr<ParamCurve> tmp = planebound->spaceCurve();
      tmp->writeStandardHeader(ofc);
      tmp->write(ofc);
    }

  std::ofstream ofp("torus_par.g2");
  for (int ka=0; ka<4; ++ka)
    {
      if (!torbound[ka].get())
	continue;
      shared_ptr<ParamCurve> tmp = torbound[ka]->parameterCurve();
      if (!tmp.get())
	continue;
      SplineDebugUtils::writeSpaceParamCurve(tmp, ofp, 0.0);
    }
#endif
  
  // Move points as appropriate
  // From torus
  vector<vector<RevEngPoint*> > move1(4);
  vector<RevEngPoint*> torus_pts;
  for (size_t ki=0; ki<blend_pts.size(); ++ki)
    {
      if (parvals[2*ki+1] > M_PI)
	move1[0].push_back(blend_pts[ki]);  // To blend cylinder
      else if (parvals[2*ki+1] > pihalf)
	move1[1].push_back(blend_pts[ki]);  // To plane
      else if (parvals[2*ki] < dom.umin())
	move1[2].push_back(blend_pts[ki]);  // First adjacent cylinder
      else if (parvals[2*ki] > dom.umax())
	move1[3].push_back(blend_pts[ki]);  // Second adjacent cylinder
      else
	{
	  blend_pts[ki]->setPar(Vector2D(parvals[2*ki],parvals[2*ki+1]));
	  blend_pts[ki]->setSurfaceDist(dist_ang[ki].first, dist_ang[ki].second);
	  torus_pts.push_back(blend_pts[ki]);
	}
    }

  // From blend cylinder
  vector<RevEngPoint*> move_cyl;
  int num_pts_cyl = adj[jx2]->numPoints();
  for (int ka=0; ka<num_pts_cyl; ++ka)
    {
      RevEngPoint* curr = adj[jx2]->getPoint(ka);
      Vector2D par = curr->getPar();
      if ((upper && par[1] > vpar3) || (upper == false && par[1] < vpar3))
	move_cyl.push_back(curr);
    }

  // Remove identified points from cylinder
  if (move_cyl.size() > 0)
    adj[jx2]->removeAndUpdatePoints(move_cyl);

  // From adjacent cylinders
  vector<vector<RevEngPoint*> > move_adj(2);
  for (int kb=0; kb<2; ++kb)
    {
      if (!adj_reg[kb])
	continue;
      int num_pts_cyl = adj_reg[kb]->numPoints();
      for (int ka=0; ka<num_pts_cyl; ++ka)
	{
	  RevEngPoint* curr = adj_reg[kb]->getPoint(ka);
	  Vector2D par = curr->getPar();
	  if ((upper2[kb] && par[1] > lim[kb]) ||
	      (upper2[kb] == false && par[1] < lim[kb]))
	    move_adj[kb].push_back(curr);
	}
      if (move_adj[kb].size() > 0)
	adj_reg[kb]->removeAndUpdatePoints(move_adj[kb]);
    }

  // From adjacent plane
  vector<RevEngPoint*> move_plane;
  adj[jx1]->extractOutOfEdge(planebound, (jx1 == 0) ? intcv1 : intcv2,
			     radius2, approx_tol_, angtol, move_plane);

  // To plane
  bool OK1, OK2, OK3;
  if (move1[1].size() > 0)
    OK1 = adj[jx1]->addPointsToGroup(move1[1], approx_tol_, angtol);

  // To blend cylinder
  if (move1[0].size() > 0)
    OK2 = adj[jx2]->addPointsToGroup(move1[0], approx_tol_, angtol);

  // To adjacent cylinders
  vector<vector<RevEngPoint*> > added_groups;
  for (int ka=0; ka<2; ++ka)
    if (move1[2+ka].size() > 0)
      {
	if (adj_reg[ka])
	  OK3 = adj_reg[ka]->addPointsToGroup(move1[2+ka], approx_tol_, angtol);
	else
	  added_groups.push_back(move1[2+ka]);
      }


  // To torus
  if (move_plane.size() > 0)
    torus_pts.insert(torus_pts.end(), move_plane.begin(), move_plane.end());
  if (move_adj[0].size() > 0)
    torus_pts.insert(torus_pts.end(), move_adj[0].begin(), move_adj[0].end());
  if (move_adj[1].size() > 0)
    torus_pts.insert(torus_pts.end(), move_adj[1].begin(), move_adj[1].end());
  if (move_cyl.size() > 0)
    torus_pts.insert(torus_pts.end(), move_cyl.begin(), move_cyl.end());
  if (torus_pts.size() == 0)
    return false;

  // Define torus blend
  shared_ptr<RevEngRegion> blendreg(new RevEngRegion(classification_type_,
						     edge_class_type_,
						     torus_pts));
  blendreg->setRegionAdjacency();
  regions_.push_back(blendreg);
  shared_ptr<HedgeSurface> hedge;
  shared_ptr<ParamSurface> torus_tmp = torus;
  blendreg->setAssociatedSurface(torus_tmp, approx_tol_, angtol,
				 min_point_region_, hedge);
  if (hedge.get())
    surfaces_.push_back(hedge);
  
  if (added_groups.size() > 0)
    {
      vector<HedgeSurface*> dummy_hedge;
      surfaceExtractOutput((int)regions_.size()-1, added_groups, dummy_hedge);
    }

#ifdef DEBUG_BLEND
  std::ofstream oft("torus_corner.g2");
  blendreg->writeRegionPoints(oft);
  adj[jx1]->writeRegionPoints(oft);
  adj[jx2]->writeRegionPoints(oft);
  if (adj_reg[0])
    adj_reg[0]->writeRegionPoints(oft);
  if (adj_reg[1])
    adj_reg[1]->writeRegionPoints(oft);
#endif
  
  // Add edges to regions
  int status = 0;
  vector<shared_ptr<ftEdge> > tor_edg(4);
  for (int ka=0; ka<4; ++ka)
    {
      if (torbound[ka].get())
	{
	  tor_edg[ka] = shared_ptr<ftEdge>(new ftEdge(hedge.get(), torbound[ka],
						      torbound[ka]->startparam(),
						      torbound[ka]->endparam()));
	  blendreg->addTrimEdge(tor_edg[ka]);
	}
    }

  shared_ptr<ftEdge> plane_edg(new ftEdge(adj[jx1]->getSurface(0), planebound,
					  planebound->startparam(),
					  planebound->endparam()));
  adj[jx1]->addTrimEdge(plane_edg);
  plane_edg->setReversed(true);
  tor_edg[3]->connectTwin(plane_edg.get(), status);


   shared_ptr<ftEdge> cyl_edg(new ftEdge(adj[jx2]->getSurface(0), blendbound,
					 blendbound->startparam(),
					 blendbound->endparam()));
  adj[jx2]->addTrimEdge(cyl_edg);
  cyl_edg->setReversed(true);
  tor_edg[2]->connectTwin(cyl_edg.get(), status);

  for (int ka=0; ka<2; ++ka)
    {
      if (cylbound[ka].get())
	{
	  shared_ptr<ftEdge> adj_edg(new ftEdge(adj_reg[ka]->getSurface(0),
						cylbound[ka],
						cylbound[ka]->startparam(),
						cylbound[ka]->endparam()));
	  adj_reg[ka]->addTrimEdge(adj_edg);
	  adj_edg->setReversed(true);
	  if (kx[ka] >= 0 && tor_edg[kx[ka]].get())
	    tor_edg[kx[ka]]->connectTwin(adj_edg.get(), status);
	}
    }

  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      //blend_regs[ki]->removeFromAdjacent();
      blend_regs[ki]->setRemove();
    }
   edges_[ix]->clearBlendRegions();
 
  return true;
}

//===========================================================================
void RevEng::setBlendBoundaries(RevEngRegion *reg)
//===========================================================================
{
  double eps = 1.0e-6;
  if (!reg->hasSurface())
    return;

  shared_ptr<ParamSurface> surf = reg->getSurface(0)->surface();
  if (surf->instanceType() != Class_Cylinder &&
      surf->instanceType() != Class_Torus)
    return;   

  double angtol = 5.0*anglim_;
  
  // Adjacent regions
  RevEngEdge *edge = reg->getBlendEdge();
  RevEngRegion* adj[2];
  edge->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return;  // Something wrong

  bool out1 = false, out2 = false;
  edge->getOuterInfo(out1, out2);

  vector<shared_ptr<CurveOnSurface> > intcv1, intcv2;
  edge->getCurve(intcv1, true);
  edge->getCurve(intcv2, false);
  //double eps1 = approx_tol_; //std::max(1.0e-4, 0.1*approx_tol_);
  // for (size_t ki=0; ki<intcv1.size(); ++ki)
  //   intcv1[ki]->fixMismatchCurves(eps1);
  // for (size_t ki=0; ki<intcv2.size(); ++ki)
  //   intcv2[ki]->fixMismatchCurves(eps1);
  
#ifdef DEBUG_BLEND
  std::ofstream of("blend.g2");
  reg->writeRegionPoints(of);
  surf->writeStandardHeader(of);
  surf->write(of);

  adj[0]->writeRegionPoints(of);
  adj[1]->writeRegionPoints(of);
#endif

  // Extract longitudial boundary curves
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  //double radius = elem->radius(0.0, 0.0);

  vector<shared_ptr<ElementarySurface> > adj_elem(2);
  for (int ka=0; ka<2; ++ka)
    {
      shared_ptr<ParamSurface> tmp = adj[ka]->getSurface(0)->surface();
      adj_elem[ka] = dynamic_pointer_cast<ElementarySurface,ParamSurface>(tmp);
    }

#ifdef DEBUG_BLEND
  std::ofstream ofsf("adj_sfs.g2");
  adj_elem[0]->writeStandardHeader(ofsf);
  adj_elem[0]->write(ofsf);
  shared_ptr<ParamCurve> tmp_cv1 = intcv1[0]->spaceCurve();
  tmp_cv1->writeStandardHeader(ofsf);
  tmp_cv1->write(ofsf);
  adj_elem[1]->writeStandardHeader(ofsf);
  adj_elem[1]->write(ofsf);
  shared_ptr<ParamCurve> tmp_cv2 = intcv2[0]->spaceCurve();
  tmp_cv2->writeStandardHeader(ofsf);
  tmp_cv2->write(ofsf);
#endif

  Point posi1, posi2, pos;
  intcv1[0]->point(posi1, intcv1[0]->startparam());
  double pari = (edge->isClosed(approx_tol_)) ?
    0.5*(intcv1[0]->startparam() + intcv1[intcv1.size()-1]->endparam()) :
    intcv1[intcv1.size()-1]->endparam();
  intcv1[intcv1.size()-1]->point(posi2, pari);
  pos = 0.5*(posi1 + posi2);
  RectDomain surfdom = surf->containingDomain();
  double regdom[4];
  reg->getDomain(regdom);
  double regfac = 2.0;
  double rad;
  double tpar1, tpar2;
  bool udir;
  int constdir;
  double delfac = 0.6;
  double seamfac = 0.1;
  bool plane1 =  (adj_elem[0]->instanceType() == Class_Plane);
  bool plane2 =  (adj_elem[1]->instanceType() == Class_Plane);
  if (surf->instanceType() == Class_Cylinder)
    {
      Point norm1 = adj_elem[0]->direction();
      Point norm2 = adj_elem[1]->direction();
      double ang = norm1.angle(norm2);
      ang = std::min(ang, M_PI-ang);
      
      tpar1 = M_PI - 0.5*ang; 
      tpar2 = tpar1 + ang;
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
      cyl->setParamBoundsU(tpar1, tpar2);
      double tdelreg = regdom[3] - regdom[2];
      if (surfdom.vmax() - surfdom.vmin() > regfac*tdelreg)
	cyl->setParamBoundsV(std::max(surfdom.vmin(), regdom[2]-delfac*tdelreg),
			     std::min(surfdom.vmax(), regdom[3]+delfac*tdelreg));
      rad = cyl->getRadius();
      udir = false;
      constdir = 0;
    }
  else
    {
      shared_ptr<Torus> tor = dynamic_pointer_cast<Torus,ParamSurface>(surf);
      RectDomain tor_dom = tor->getParameterBounds();
      if (plane1 || plane2)
	{
	  int plane_ix = (plane1) ? 0 : 1;
	  if (adj_elem[plane_ix]->instanceType() != Class_Plane)
	    return;  // Not as expected
	  rad = tor->getMinorRadius();
	  bool plane_out = (plane_ix == 0) ? out1 : out2;
	  bool rot_out = (plane_ix == 0) ? out2 : out1;
	  double phi = 0.5*M_PI;
	  Point norm = adj_elem[plane_ix]->direction();
	  Point norm2 = adj[plane_ix]->getMeanNormalTriang();
	  double beta = 0.5; //(norm*norm2 > 0.0) ? 0.5 : 0.0;
	  if (adj_elem[1-plane_ix]->instanceType() == Class_Cone)
	    {
	      shared_ptr<Cone> cone =
		dynamic_pointer_cast<Cone,ElementarySurface>(adj_elem[1-plane_ix]);
	      double alpha = cone->getConeAngle();
	      Point axis = adj_elem[1-plane_ix]->direction();
	      int sgn = (norm*axis < 0.0) ? -1 : 1;
	      phi += sgn*alpha;
	    }
	  if (rot_out)
	    {
	      tpar2 = tor_dom.vmin()+(1.0+beta)*M_PI;
	      tpar1 = tpar2 - phi;
	    }
	  else
	    {
	      tpar1 = tor_dom.vmin()+(1.0-beta)*M_PI;
	      tpar2 = tpar1 + phi;
	    }
	  if (plane_out)
	    {
	      int sgn = 1; //(norm*norm2 < 0.0 && rot_out == false) ? -1 : 1;
	      tpar1 += 0.5*sgn*M_PI;
	      tpar2 += 0.5*sgn*M_PI;
	    }
	}
      else
	{
	  int cyl_ix = (adj_elem[0]->instanceType() == Class_Cylinder) ? 0 : 1;
	  int cone_ix = 1 - cyl_ix;
	  if (adj_elem[cyl_ix]->instanceType() != Class_Cylinder ||
	      adj_elem[cone_ix]->instanceType() != Class_Cone)
	    return;

	  shared_ptr<Cone> cone =
	    dynamic_pointer_cast<Cone,ElementarySurface>(adj_elem[cone_ix]);
	  double alpha = cone->getConeAngle();
	  double beta = M_PI - fabs(alpha);
	  double phi = 0.5*(M_PI - beta);
	  rad = tor->getMinorRadius();
	  Point axis = adj_elem[cyl_ix]->direction();
	  Point loc = adj_elem[cyl_ix]->location();
	  double rd = (pos-loc)*axis;
	  Point centr = loc + rd*axis;
	  Point loc2 = tor->location();
	  Point loc2_2 = centr + ((loc2 - centr)*axis)*axis;
	  double hh = loc2_2.dist(centr);
	  double d2 = hh/sin(phi) - rad;
	  double xlen = sqrt((rad+d2)*(rad+d2) - rad*rad);
	  tpar1 = tor_dom.vmin() + M_PI - 2.0*xlen;
	  tpar2 = tpar1 + 2*xlen;
	}
      
      double tdelreg = regdom[1] - regdom[0];
      double upar1 = tor_dom.umin();
      double upar2 = tor_dom.umax();

      double cp1, cp2;
      Point cpos1, cpos2;
      double seed[2];
      seed[0] = 0.5*(regdom[0]+regdom[1]);
      seed[1] = 0.5*(regdom[2]+regdom[3]);
      int seam1 = edge->closedSfAtEnd(approx_tol_, cp1, cpos1, true);
      int seam2 = edge->closedSfAtEnd(approx_tol_, cp2, cpos2, false);
      double dd = cpos1.dist(cpos2);
      double u1, v1, u2, v2, d1, d2;
      Point cl1, cl2;
      elem->closestPoint(cpos1, u1, v1, cl1, d1, eps, 0, seed);
      elem->closestPoint(cpos2, u2, v2, cl2, d2, eps, 0, seed);
      if (seam1 != seam2 && dd > approx_tol_)
	{
	  bool close1, close2;
	  tor->isClosed(close1, close2);  // Expects close1 = true
	  if (close1)
	    {
	      // Check for seam
	      double umid = 0.5*(regdom[0]+regdom[1]);
	      double eps2 = std::max(eps, 0.001*(upar2-upar1));
	      if (fabs(u1-upar1) < eps2 && umid > u2)
		u1 = upar2;
	      else if (fabs(upar2-u1) < eps2 && umid < u2)
		u1 = upar1;
	      if (fabs(u2-upar1) < eps2 && umid > u1)
		u2 = upar2;
	      else if (fabs(upar2-u2) < eps2 && umid < u1)
		u2 = upar1;
	    }
	      
	  upar1 = std::max(upar1, std::min(u1,u2));
	  upar2 = std::min(upar2, std::max(u1,u2));
	}
      
      if (upar2 - upar1 > regfac*tdelreg)
	{
	  if (u2 < u1)
	    std::swap(u1, u2);
	  
	  upar1 = std::max(upar1, std::min(regdom[0]-delfac*tdelreg, u1));
	  upar2 = std::min(upar2, std::max(regdom[1]+delfac*tdelreg, u2));
	  if (upar1 < seamfac)
	    upar1 = 0.0;
	  if (2*M_PI-upar2 < seamfac)
	    upar2 = 2*M_PI;
	}
      tor->setParameterBounds(upar1, tpar1, upar2, tpar2);
      udir = true;
      constdir = 1;
    }


  vector<shared_ptr<ParamCurve> > bdcv1 = elem->constParamCurves(tpar1, udir); 
  vector<shared_ptr<ParamCurve> > bdcv2 = elem->constParamCurves(tpar2, udir);
  if (bdcv1.size() != 1 || bdcv2.size() != 1)
    return;  // Something strange

  shared_ptr<ElementaryCurve> space[2];
  space[0] = dynamic_pointer_cast<ElementaryCurve,ParamCurve>(bdcv1[0]);
  space[1] = dynamic_pointer_cast<ElementaryCurve,ParamCurve>(bdcv2[0]);
  if ((!space[0].get()) || (!space[1].get()))
    return;
  BoundingBox bb1 = adj[0]->getBbox();
  BoundingBox bb2 = adj[1]->getBbox();
  double diag = std::min(bb1.low().dist(bb1.high()), bb2.low().dist(bb2.high()));
  if (space[0]->startparam() < -0.5*diag || space[0]->endparam() > 0.5*diag)
    {
      for (int ka=0; ka<2; ++ka)
	space[ka]->setParamBounds(std::max(space[0]->startparam(), -0.5*diag),
				  std::min(space[0]->endparam(), 0.5*diag));
    }
  
#ifdef DEBUG_BLEND
  space[0]->writeStandardHeader(of);
  space[0]->write(of);
  space[1]->writeStandardHeader(of);
  space[1]->write(of);
#endif
  RectDomain dom = elem->getParameterBounds();
  Point parpt1(dom.umin(), dom.vmin()), parpt2(dom.umin(), dom.vmax());
  Point parpt3(dom.umax(), dom.vmin()), parpt4(dom.umax(), dom.vmax());
  Point lpos1 = (udir) ? Point(0.0, dom.vmin()) : Point(dom.umin(), 0.0);
  Point lpos2 = (udir) ? Point(0.0, dom.vmax()) : Point(dom.umax(), 0.0);
  Point pvec = (udir) ? Point(1.0, 0.0) : Point(0.0, 1.0);

  double t1 = space[0]->startparam();
  double t2 = space[0]->endparam();
  shared_ptr<ElementaryCurve> par1(new Line(lpos1, pvec));
  par1->setParamBounds(t1, t2);
  shared_ptr<ElementaryCurve> par2(new Line(lpos2, pvec));
  par2->setParamBounds(t1, t2);
#ifdef DEBUG_BLEND
  // Check
  double tm1 = 0.75*t1 + 0.25*t2;
  double tm2 = 0.25*t1 + 0.75*t2;
  Point pp1, pp2, pp3, pp4;
  Point sp1, sp2, sp3, sp4;
  Point ssp1, ssp2, ssp3, ssp4;
  space[0]->point(sp1, tm1);
  space[0]->point(sp2, tm2);
  space[1]->point(sp3, tm1);
  space[1]->point(sp4, tm2);
  par1->point(pp1, tm1);
  par1->point(pp2, tm2);
  par2->point(pp3, tm1);
  par2->point(pp4, tm2);
  surf->point(ssp1, pp1[0], pp1[1]);
  surf->point(ssp2, pp2[0], pp2[1]);
  surf->point(ssp3, pp3[0], pp3[1]);
  surf->point(ssp4, pp4[0], pp4[1]);
#endif

  shared_ptr<ElementaryCurve> adj_par[2];
  Point close1, close2;
  double upar1, upar2, vpar1, vpar2, dist1, dist2;
  Point pos1, pos2;
  int ix[2];
  ix[0] = 0;
  ix[1] = 1;
  
  // Check for curve matching
  space[0]->point(pos1, t1);
  adj_elem[0]->closestPoint(pos1, upar1, vpar1, close1, dist1, eps);
  adj_elem[1]->closestPoint(pos1, upar2, vpar2, close2, dist2, eps);
  if (dist2 < dist1)
    std::swap(ix[0], ix[1]);

  double tol2 = 2.0*approx_tol_;
  for (int ka=0; ka<2; ++ka)
    {
      adj_par[ix[ka]] = adj_elem[ix[ka]]->getElementaryParamCurve(space[ka].get(),
								  tol2);
      if (!adj_par[ix[ka]].get())
	{
	  std::cout << "No parameter curve found" << std::endl;
	}

      // space[ka]->point(pos1, t1);
      // space[ka]->point(pos2, t2);
      // adj_elem[ix[ka]]->closestPoint(pos1, upar1, vpar1, close1, dist1, eps);
      // adj_elem[ix[ka]]->closestPoint(pos2, upar2, vpar2, close2, dist2, eps);
      // Point pp1(upar1,vpar1), pp2(upar2,vpar2);
      // Point mpar = (t2*pp1 -t1*pp2)/(t2-t1);
      // Point mvec = pp2 - pp1;
      // mvec.normalize();

      // adj_par[ix[ka]] = shared_ptr<ElementaryCurve>(new Line(mpar, mvec));
      // adj_par[ix[ka]]->setParamBounds(t1, t2);
    }
  
  // Create edges. Orientation in adjacent regions is not set
  int status = 0;
  vector<shared_ptr<ftEdge> > bdedg(2);
  shared_ptr<CurveOnSurface> sfcv1(new CurveOnSurface(elem, par1, space[0], false, 3,
						      constdir+1, tpar1, 2*constdir, true));
  bdedg[0] = shared_ptr<ftEdge>(new ftEdge(reg->getSurface(0), sfcv1, t1, t2));
  shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(elem, par2, space[1], false, 3,
						      constdir+1, tpar2, 2*constdir+1, true));
  //sfcv2->reverseParameterDirection();
  bdedg[1] = shared_ptr<ftEdge>(new ftEdge(reg->getSurface(0), sfcv2, t1, t2));
  reg->addTrimEdge(bdedg[0]);
  reg->addTrimEdge(bdedg[1]);

  shared_ptr<CurveOnSurface> adj_sfcv[2];
  shared_ptr<ftEdge> adj_edg[2];
  vector<RevEngPoint*> out[2];
  for (int ka=0; ka<2; ++ka)
    {
      adj_sfcv[ix[ka]] = shared_ptr<CurveOnSurface>(new CurveOnSurface(adj_elem[ix[ka]],
								       adj_par[ix[ka]],
								       space[ka],
								       false, 1));
      // if (!adj_sfcv[ix[ka]]->hasParameterCurve())
      // 	{
      // 	  vector<shared_ptr<CurveOnSurface> > tmp_cvs;
      // 	  tmp_cvs.push_back(adj_sfcv[ix[ka]]);
      // 	  vector<pair<double,double> > t1_t2;
      // 	  adj[ix[ka]]->getCurveRestriction(tmp_cvs, approx_tol_, anglim_, t1_t2);
      // 	  if (t1_t2.size() == 1)
      // 	    {
      // 	      if (t1_t2[0].first > adj_sfcv[ix[ka]]->startparam() ||
      // 		  t1_t2[0].second < adj_sfcv[ix[ka]]->endparam())
      // 		{
      // 		}
      // 	    }
      // 	}
      adj_edg[ix[ka]] = shared_ptr<ftEdge>(new ftEdge(adj[ix[ka]]->getSurface(0),
						      adj_sfcv[ix[ka]], t1, t2));

      adj[ix[ka]]->addTrimEdge(adj_edg[ix[ka]]);
      bdedg[ka]->setReversed(true);
      adj_edg[ix[ka]]->connectTwin(bdedg[ka].get(), status);
  
      // Identify points from the adjacent regions lying outside the corresponding
      // trimming curve
      adj[ix[ka]]->extractOutOfEdge(adj_sfcv[ix[ka]],
				    (ix[ka] == 0) ? intcv1 : intcv2,
				    rad, approx_tol_, angtol, out[ix[ka]]);
    }

#ifdef DEBUG_BLEND
  std::ofstream of2("out_points.g2");
  for (int ka=0; ka<2; ++ka)
    {
      if (out[ka].size() > 0)
	{
	  of2 << "400 1 0 4 0 255 0 255" << std::endl;
	  of2 << out[ka].size() << std::endl;
	  for (size_t kr=0; kr<out[ka].size(); ++kr)
	    of2 << out[ka][kr]->getPoint() << std::endl;
	}
    }
#endif


  // Add out-points to the blend regions.
  // NB! This can be too simple in a more complex configuration.
  // Let's wait for the problem to turn up
  if (out[0].size() > 0 || out[1].size() > 0)
    {
      vector<RevEngPoint*> points;
      for (int ka=0; ka<2; ++ka)
	if (out[ka].size() > 0)
	  points.insert(points.end(), out[ka].begin(), out[ka].end());
      bool integrated = reg->addPointsToGroup(points, approx_tol_, angtol);
      if (!integrated)
	{
	  vector<HedgeSurface*> dummy_sfs;
	  vector<RevEngPoint*> dummy_pts;
	  for (int ka=0; ka<2; ++ka)
	    {
	      if (out[ix[ka]].size() > 0)
		{
		  vector<vector<RevEngPoint*> > out_1;
		  adj[ix[ka]]->connectedGroups(out[ix[ka]], out_1, false, dummy_pts);
		  for (size_t kr=0; kr<regions_.size(); ++kr)
		    if (regions_[kr].get() == adj[ix[ka]])
		      {
			surfaceExtractOutput((int)kr, out_1, dummy_sfs);
			break;
		      }
		}
	    }
	}
    }

  // Identify points associated to the blend region that should be
  // moved to the adjacent regions
  vector<int> adj_ix(4);
  adj_ix[0] = 3;
  adj_ix[1] = 1;
  adj_ix[2] = 0;
  adj_ix[3] = 2;
  vector<vector<RevEngPoint*> > move2adj(4);
  vector<RevEngPoint*> remain;
  vector<RevEngPoint*> regpoints = reg->getPoints();
  extractOutPoints(regpoints, elem, adj_ix, approx_tol_, angtol,
		   move2adj, remain);

#ifdef DEBUG_BLEND
  std::ofstream of2_3("in_points.g2");
  for (int ka=0; ka<4; ++ka)
    {
      if (move2adj[ka].size() > 0)
	{
	  of2_3 << "400 1 0 4 100 155 0 255" << std::endl;
	  of2_3 << move2adj[ka].size() << std::endl;
	  for (size_t kr=0; kr<move2adj[ka].size(); ++kr)
	    of2_3 << move2adj[ka][kr]->getPoint() << std::endl;
	}
    }
#endif

  int kx = 2*(1-constdir);
  for (int ka=kx; ka<=kx+1; ++ka)
    {
      if (move2adj[ka].size() > 0)
	{
	  remain.insert(remain.end(), move2adj[ka].begin(), move2adj[ka].end());
	  move2adj[ka].clear();
	}
    }

  for (int ka=0; ka<=1; ++ka)
    {
      int kb = 2*constdir+ka;
      if (move2adj[kb].size() > 0)
	{
	  reg->removePoints(move2adj[kb]);
	  bool integrated = adj[ix[ka]]->addPointsToGroup(move2adj[kb],
							  approx_tol_, angtol);
	  if (!integrated)
	    MESSAGE("RevEng::setBlendBoundaris. Missing adjacent surface");
	}
    }
  
#ifdef DEBUG_BLEND
  std::ofstream of3("updated_blend.g2");
  reg->writeRegionPoints(of3);
  adj[0]->writeRegionPoints(of3);
  adj[1]->writeRegionPoints(of3);
#endif
 int stop_break = 1;
}

//===========================================================================

// Service functionality for suitcaseCorner

bool getBlendRegMatches(vector<RevEngRegion*>& adj_blends,
			vector<vector<shared_ptr<ftEdge> > >& trim_edgs,
			vector<vector<pair<double, double> > >& par_lim,
			vector<vector<size_t> >& match)
{
#ifdef DEBUG_BLEND
  std::ofstream of("int_pt.g2");
#endif
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_blends.size(); ++kj)
      {
	size_t kr, kh;
	for (kr=0; kr<trim_edgs[ki].size(); ++kr)
	  {
	    ftEdgeBase *edg1 = trim_edgs[ki][kr]->twin();
	    ftFaceBase *face1 = (edg1->geomEdge()) ?
	      edg1->geomEdge()->face() : 0;
	    if (!face1)
	      return false;
	    for (kh=0; kh<trim_edgs[kj].size(); ++kh)
	      {
		ftEdgeBase *edg2 = trim_edgs[kj][kh]->twin();
		ftFaceBase *face2 = (edg2->geomEdge()) ?
		  edg2->geomEdge()->face() : 0;
		if (!face2)
		  return false;
		if (face1 == face2)
		  break;
	      }
	    if (kh < trim_edgs[kj].size())
	      break;
	  }
	if (kr == trim_edgs[ki].size() || kh == trim_edgs[kj].size())
	  continue;

	shared_ptr<ParamCurve> cv1 = trim_edgs[ki][kr]->geomCurve();
	shared_ptr<ParamCurve> cv2 = trim_edgs[kj][kh]->geomCurve();
	double par1, par2, dist;
	Point ptc1, ptc2;
	ClosestPoint::closestPtCurves(cv1.get(), cv2.get(), par1, par2,
				      dist, ptc1, ptc2);
#ifdef DEBUG_BLEND
	of << "400 1 0 4 155 100 0 255" << std::endl;
	of << "1" << std::endl;
	of << ptc1 << std::endl;
	of << "400 1 0 4 155 100 0 255" << std::endl;
	of << "1" << std::endl;
	of << ptc2 << std::endl;
#endif
	par_lim[ki][kr] = std::make_pair(par1, dist);
	par_lim[kj][kh] = std::make_pair(par2, dist);
	vector<size_t> match0{ki, kr, kj, kh};
	match.push_back(match0);
      }
  return true;
}


bool getTrimEdgeMidpoint(vector<vector<shared_ptr<ftEdge> > >& trim_edgs,
			 vector<vector<pair<double, double> > >& par_lim,
			 vector<vector<pair<double, double> > >& midp)
{
  for (size_t ki=0; ki<par_lim.size(); ++ki)
    {
      int num = 0;
      double par = 0.0;
      for (size_t kj=0; kj<par_lim[ki].size(); ++kj)
	if (par_lim[ki][kj].second >= 0.0)
	  {
	    par += par_lim[ki][kj].first;
	    ++num;
	  }
      if (num != 2)
	return false;
      par /= (double)num;
      midp[ki].resize(par_lim[ki].size());
      for (size_t kj=0; kj<par_lim[ki].size(); ++kj)
	{
	  if (par_lim[ki][kj].second < 0.0)
	    {
	      midp[ki][kj] = std::make_pair(0.0, -1.0);
	      continue;
	    }
	  shared_ptr<ParamCurve> cv = trim_edgs[ki][kj]->geomCurve();
	  Point pt1 = cv->point(par_lim[ki][kj].first);
	  Point pt2 = cv->point(par);
	  double dist = pt1.dist(pt2);
	  midp[ki][kj] = std::make_pair(par,dist);
	  //int stop_break = 1;
	}
    }
  return true;
}

bool computeCornerBoundaryCurves(vector<RevEngRegion*>& adj_blends,
				 vector<vector<shared_ptr<ftEdge> > >& trim_edgs,
				 vector<vector<pair<double, double> > >& par_lim,
				 vector<vector<pair<double, double> > >& midp,
				 vector<vector<size_t> >& match, double tol,
				 vector<shared_ptr<CurveOnSurface> >& blend_bd,
				 vector<RevEngRegion*>& adjreg)
{
  double tol1 = 0.5*tol;

  // Want four boundary curves if possible
  vector<double> midd(midp.size());
  for (size_t ki=0; ki<midp.size(); ++ki)
    {
      double middist = 0.0;
      for (size_t kj=0; kj<midp[ki].size(); ++kj)
	middist = std::max(middist, midp[ki][kj].second);
      midd[ki] = middist;
    }
  vector<double> midd2(midd.begin(), midd.end());
  
  std::sort(midd2.begin(), midd2.end());
  double tol2;
  if (midd2.size() > 4 && midd2.size() < 2)
    return false;  // Currently not handled
  else if (midd2.size() == 4)
    tol2 = 2.0*midd2[3];
  else if (midd2.size() == 2)
    tol2 = tol1;
  else
    tol2 = 0.5*(midd2[0] + midd2[1]);
  
  adjreg.insert(adjreg.end(), adj_blends.begin(), adj_blends.end());

  vector<pair<double,double> > cvpar(adj_blends.size());
  vector<bool> parset(adj_blends.size(), false);

  // Start defining iso-parametric curves due to close intersection points
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      if (midd[ki] <= tol1)
	{
	  double par;
	  size_t kj;
	  for (kj=0; kj<midp[ki].size() && midp[ki][kj].second < 0.0; ++kj);
	  par = midp[ki][kj].first;
	  cvpar[ki] = std::make_pair(par, par);
	  parset[ki] = true;
	}
    }

  // Continue with straight curve between opposite intersection points
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      if (midd[ki] <= tol2)
	{
	  double par1, par2;
	  size_t kj, kr;
	  for (kj=0; kj<par_lim[ki].size() && par_lim[ki][kj].second < 0.0; ++kj);
	  par1 = par_lim[ki][kj].first;
	  for (kr=kj+1; kr<par_lim[ki].size() && par_lim[ki][kr].second < 0.0; ++kr);
	  par2 = par_lim[ki][kr].first;
	  cvpar[ki] = std::make_pair(par1, par2);
	  parset[ki] = true;
	}
    }

  if (adj_blends.size() != 3)
    return false;  // Waiting for a test case

  // Set iso-curve parameters for adjacent blends
  size_t kj;
  vector<pair<size_t, size_t> > missing;
  vector<double> div_par;
  for (kj=0; kj<parset.size(); ++kj)
    {
      if (parset[kj])
	{
	  size_t kh, kr;
	  for (kh=0; kh<par_lim[kj].size() && par_lim[kj][kh].second < 0.0; ++kh);
	  for (kr=kh+1; kr<par_lim[kj].size() && par_lim[kj][kr].second < 0.0; ++kr);

	  int ix1=-1, ix2=-1;
	  double par1, par2;
	  for (size_t kv=0; kv<match.size(); ++kv)
	    {
	      if (match[kv][0] == kj)
		{
		  if (match[kv][1] == kh)
		    {
		      ix1 = match[kv][2];
		      par1 = par_lim[ix1][match[kv][3]].first;
		      missing.push_back(std::make_pair(ix1, 1-match[kv][3]));
		      div_par.push_back(par1);
		    }		      
		  else if (match[kv][1] == kr)
		    {
		      ix2 = match[kv][2];
		      par2 = par_lim[ix2][match[kv][3]].first;
		      missing.push_back(std::make_pair(ix2, 1-match[kv][3]));
		      div_par.push_back(par2);
		    }
		}
	      else if (match[kv][2] == kj)
		{
		  if (match[kv][3] == kh)
		    {
		      ix1 = match[kv][0];
		      par1 = par_lim[ix1][match[kv][1]].first;
		      missing.push_back(std::make_pair(ix1, 1-match[kv][1]));
		      div_par.push_back(par1);
		    }
		  else if (match[kv][3] == kr)
		    {
		      ix2 = match[kv][0];
		      par2 = par_lim[ix2][match[kv][1]].first;
		      missing.push_back(std::make_pair(ix2, 1-match[kv][1]));
		      div_par.push_back(par2);
		    }
		}
	    }
	  if (ix1 < 0 || ix2 < 0)
	    return false;
	  cvpar[ix1] = std::make_pair(par1, par1);
	  parset[ix1] = true;
	  cvpar[ix2] = std::make_pair(par2, par2);
	  parset[ix2] = true;
	  break;
	}
    }
  if (kj == parset.size())
    return false;
  
  double eps = 1.0e-9;
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      if (!parset[ki])
	return false;  // Missing information
      if (fabs(cvpar[ki].first-cvpar[ki].second) < eps)
	{
	  // Iso-parametric curve
	  size_t kj;
	  for (kj=0; kj<midp[ki].size() && midp[ki][kj].second < 0.0; ++kj);
	  double par = 0.5*(cvpar[ki].first+cvpar[ki].second);
	  shared_ptr<ParamCurve> cv = trim_edgs[ki][kj]->geomCurve();
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  shared_ptr<ParamCurve> pcv = sfcv->parameterCurve();
	  shared_ptr<ParamSurface> surf = sfcv->underlyingSurface();
	  if ((!pcv.get()) || (!surf.get()))
	    return false;
	  double val;
	  int dir;
	  bool isconst = sfcv->isConstantCurve(tol1, dir, val);
	  if (!isconst)
	    return false;
	  RectDomain dom = surf->containingDomain();
	  double pmin = (dir == 1) ? dom.umin() : dom.vmin();
	  double pmax = (dir == 1) ? dom.umax() : dom.vmax();
	  shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(surf, (3-dir),
							      par, pmin, pmax,
							      -1));
	  blend_bd.push_back(sfcv2);
	}
      else
	{
	  // Straight curve in the parameter domain
	  size_t kj, kr;
	  for (kj=0; kj<par_lim[ki].size() && par_lim[ki][kj].second < 0.0; ++kj);
	  for (kr=kj+1; kr<par_lim[ki].size() && par_lim[ki][kr].second < 0.0; +kr);
	  double par1 = cvpar[ki].first;
	  double par2 = cvpar[ki].second;
	  shared_ptr<ParamCurve> cv1 = trim_edgs[ki][kj]->geomCurve();
	  shared_ptr<ParamCurve> cv2 = trim_edgs[ki][kr]->geomCurve();
	  shared_ptr<CurveOnSurface> sfcv1 =
			dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
	  shared_ptr<CurveOnSurface> sfcv2 =
			dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
	  shared_ptr<ParamCurve> pcv1 = sfcv1->parameterCurve();
	  shared_ptr<ParamCurve> pcv2 = sfcv2->parameterCurve();
	  if ((!pcv1.get()) || (!pcv2.get()))
	    return false;
	  Point ppt1 = pcv1->point(par1);
	  Point ppt2 = pcv2->point(par2);
	  shared_ptr<SplineCurve> pcv3(new SplineCurve(ppt1, ppt2));
	  shared_ptr<CurveOnSurface> sfcv3(new CurveOnSurface(sfcv1->underlyingSurface(),
							      pcv3, true));
	  bool space = sfcv3->ensureSpaceCrvExistence(tol1);
	  if (!space)
	    return false;
	  blend_bd.push_back(sfcv3);
	}
    }

#ifdef DEBUG_BLEND
  std::ofstream of1("space_bd.g2");
  for (size_t kr=0; kr<blend_bd.size(); ++kr)
    {
      shared_ptr<ParamCurve> space = blend_bd[kr]->spaceCurve();
      shared_ptr<ParamCurve> parcv = blend_bd[kr]->parameterCurve();
      space->writeStandardHeader(of1);
      space->write(of1);
    }
#endif

  if (missing.size() > 0)
    {
      // Add remaining curves
      size_t ki, kj;
      for (ki=0; ki<missing.size(); ++ki)
	for (kj=ki+1; kj<missing.size(); ++kj)
	  {
	    size_t kr;
	    for (kr=0; kr<match.size(); ++kr)
	      if (missing[ki].first == match[kr][0] &&
		  missing[ki].second == match[kr][1] &&
		  missing[kj].first == match[kr][2] &&
		  missing[kj].second == match[kr][3])
		break;
	    if (kr == match.size())
	      continue; // Should not happen
	    size_t ki1 = missing[ki].first, kr1 = missing[ki].second;
	    size_t kj1 = missing[kj].first, kh1 = missing[kj].second;
	    ftEdgeBase *edg = trim_edgs[ki1][kr1]->twin();
	    ftSurface *face = edg->geomEdge()->face()->asFtSurface();
	    if (!face)
	      continue;
	    shared_ptr<ParamSurface> surf = face->surface();
	    shared_ptr<ParamCurve> cv1 = trim_edgs[ki1][kr1]->geomCurve();
	    shared_ptr<ParamCurve> cv2 = trim_edgs[kj1][kh1]->geomCurve();
	    vector<Point> der1(2), der2(2);
	    cv1->point(der1, div_par[ki], 1);
	    cv2->point(der2, div_par[kj], 1);
	    double dd = der1[0].dist(der2[0]);
	    der1[1].normalize();
	    der2[1].normalize();
	    if ((der2[0]-der1[0])*der1[1] < 0.0)
	      der1[1] *= -1;
	    if ((der2[0]-der1[0])*der2[1] < 0.0)
	      der2[1] *= -1;
#ifdef DEBUG_BLEND
	    std::ofstream of3("missing_cvs.g2");
	    of3 << "400 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der1[0] << std::endl;
	    of3 << "410 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der1[0] << " " << der1[0]+0.3*dd*der1[1] << std::endl;
	    of3 << "400 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der2[0] << std::endl;
	    of3 << "410 1 0 4 255 0 0 255" << std::endl;
	    of3 << "1" << std::endl;
	    of3 << der2[0] << " " << der2[0]+0.3*dd*der2[1] << std::endl;
#endif
	    vector<double> knots(8, 0.0);
	    for (int ka=0; ka<4; ++ka)
	      knots[4+ka] = dd;
	    dd /= 3.0;
	    vector<double> coefs;
	    coefs.insert(coefs.end(), der1[0].begin(), der1[0].end());
	    Point cf1 = der1[0] + dd*der1[1];
	    coefs.insert(coefs.end(), cf1.begin(), cf1.end());
	    Point cf2 = der2[0] - dd*der2[1];
	    coefs.insert(coefs.end(), cf2.begin(), cf2.end());
	    coefs.insert(coefs.end(), der2[0].begin(), der2[0].end());
	    shared_ptr<ParamCurve> mcv(new SplineCurve(4, 4, &knots[0],
							&coefs[0], 3));
#ifdef DEBUG_BLEND
	    mcv->writeStandardHeader(of3);
	    mcv->write(of3);
	    surf->writeStandardHeader(of3);
	    surf->write(of3);
#endif

	    shared_ptr<SplineCurve> space_proj, par_proj;
	    CurveCreators::projectCurve(mcv, surf, tol1, space_proj,
					par_proj);

	    shared_ptr<CurveOnSurface> sfcv_proj(new CurveOnSurface(surf,
								    par_proj,
								    space_proj,
								    true, 1));
#ifdef DEBUG_BLEND
	    space_proj->writeStandardHeader(of3);
	    space_proj->write(of3);
#endif
	    blend_bd.push_back(sfcv_proj);
	    HedgeSurface *regface =  dynamic_cast<HedgeSurface*>(face);
	    if (!regface)
	      return false;
	    adjreg.push_back(regface->getRegion(0));
	  }
    }
  return true;
}

bool getCoonsBoundaryInfo(vector<shared_ptr<CurveOnSurface> >& blend_bd,
			  double tol,
			  vector<shared_ptr<ParamCurve> >& bdcvs,
			  vector<shared_ptr<ParamCurve> >& crosscvs)
{
  bdcvs.resize(blend_bd.size());
  crosscvs.resize(blend_bd.size());
  
  for (size_t ki=0; ki<blend_bd.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv(blend_bd[ki]->spaceCurve()->clone());
      if (tmp_cv->instanceType() == Class_Circle)
	{
	  shared_ptr<Circle> circ = dynamic_pointer_cast<Circle,ParamCurve>(tmp_cv);
	  bdcvs[ki] = shared_ptr<SplineCurve>(circ->createNonRationalSpline(tol));
	}
      else if (tmp_cv->instanceType() == Class_Line)
	{
	  shared_ptr<Line> line = dynamic_pointer_cast<Line,ParamCurve>(tmp_cv);
	  bdcvs[ki] = shared_ptr<SplineCurve>(line->createSplineCurve());
	}
      else if (tmp_cv->instanceType() == Class_SplineCurve)
	{
	  shared_ptr<SplineCurve> cv = dynamic_pointer_cast<SplineCurve,ParamCurve>(tmp_cv);
	  if (cv->rational())
	    {
	      shared_ptr<ParamCurve> tmp_par = blend_bd[ki]->parameterCurve();
	      shared_ptr<ParamSurface> tmp_sf = blend_bd[ki]->underlyingSurface();
	      bdcvs[ki] =
		shared_ptr<SplineCurve>(CurveCreators::liftParameterCurve(tmp_par,
									 tmp_sf,
									 tol));
	    }
	  else
	    bdcvs[ki] = tmp_cv;
	}
      else
	return false;
      shared_ptr<CurveOnSurface> tmp_sfcv(blend_bd[ki]->clone());
      tmp_sfcv->setSpaceCurve(bdcvs[ki]);
      crosscvs[ki] = CreatorsUtils::createCrossTangent(*tmp_sfcv);
    }

  // Ensure correct direction of cross tangent curves
  for (size_t ki=0; ki<crosscvs.size(); ++ki)
    {
      size_t kj = (ki == crosscvs.size()-1) ? 0 : ki+1;
      Point ctan = crosscvs[ki]->point(crosscvs[ki]->endparam());
      double tdel = bdcvs[kj]->endparam() - bdcvs[kj]->startparam();
      Point pos1 = bdcvs[kj]->point(bdcvs[kj]->startparam());
      Point pos2 = bdcvs[kj]->point(bdcvs[kj]->startparam() + 0.1*tdel);
      Point vec = pos2 - pos1;
      if (ctan*vec <  0.0)
	{
	  shared_ptr<SplineCurve> cross =
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(crosscvs[ki]);
	  for (auto it=cross->coefs_begin(); it!=cross->coefs_end(); ++it)
	    (*it) *= -1.0;
	}
    }
  return true;
}

void pairOfRegionEdges(RevEngRegion* blendreg, HedgeSurface *hedge,
		       shared_ptr<ParamCurve>& bdcv1,
		       shared_ptr<CurveOnSurface>& bdcv2,
		       RevEngRegion* adjreg, Point& mid, double tol)
{
  int stat = 0;
  double eps = 1.0e-9;
  shared_ptr<ftEdge> blend_edge(new ftEdge(hedge, bdcv1,
					   bdcv1->startparam(),
					   bdcv1->endparam()));
  blendreg->addTrimEdge(blend_edge);

  int pdir;
  double pval;
  if (adjreg->hasBlendEdge() && bdcv2->isConstantCurve(tol, pdir, pval))
    {
      shared_ptr<ElementarySurface> adj_elem =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(adjreg->getSurface(0)->surface());
      if (adj_elem.get())
	{
	  RectDomain adjdom = adj_elem->getParameterBounds();
	  double regdom[4];
	  adjreg->getDomain(regdom);
	  double tmin = regdom[2*(pdir-1)];
	  double tmax = regdom[2*(pdir-1)+1];
	  double upar, vpar, dist;
	  Point close;
	  adj_elem->closestPoint(mid, upar, vpar, close, dist, eps);
	  if (pdir == 1)
	    {
	      if (upar > pval)
		adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(),
					     pval, adjdom.vmax());
	      else
		adj_elem->setParameterBounds(pval, adjdom.vmin(),
					     adjdom.umax(), adjdom.vmax());
	    }
	  else
	    {
	      if (vpar > pval)
		adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(), 
					     adjdom.umax(), pval);
	      else
		adj_elem->setParameterBounds(adjdom.umin(), pval, 
					     adjdom.umax(), adjdom.vmax());
	    }
	  // if (fabs(tmax-pval) > fabs(pval-tmin))
	  //   {
	  //     if (pdir == 1)
	  // 	adj_elem->setParameterBounds(pval, adjdom.vmin(),
	  // 				     adjdom.umax(), adjdom.vmax());
	  //     else
	  // 	adj_elem->setParameterBounds(adjdom.umin(), pval, 
	  // 				     adjdom.umax(), adjdom.vmax());
	  //   }
	  // else
	  //   {
	  //     if (pdir == 1)
	  // 	adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(),
	  // 				     pval, adjdom.vmax());
	  //     else
	  // 	adj_elem->setParameterBounds(adjdom.umin(), adjdom.vmin(), 
	  // 				     adjdom.umax(), pval);
	  //   }
	}
      adjreg->adaptEdges();
    }

  HedgeSurface *other_hedge = adjreg->getSurface(0);
  shared_ptr<ftEdge> other_edge(new ftEdge(other_hedge, bdcv2,
					   bdcv2->startparam(),
					   bdcv2->endparam()));
  adjreg->addTrimEdge(other_edge);
  other_edge->setReversed(true);
  blend_edge->connectTwin(other_edge.get(), stat);
}

//===========================================================================
void RevEng::extractOutPoints(vector<RevEngPoint*>& points, shared_ptr<ParamSurface> surf,
			      vector<int>& cv_ix,
			      double tol, double angtol,
			      vector<vector<RevEngPoint*> >& move2adj,
			      vector<RevEngPoint*>& remain)
//===========================================================================
{
  double eps = 1.0e-9;
  double maxd, avd;
  int num_in, num2_in;
  vector<double> parvals;
  vector<pair<double,double> > dist_ang;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(points.begin(), points.end(), surf,
			  tol, maxd, avd, num_in, num2_in, inpt, outpt,
			  parvals, dist_ang, angtol);

  RectDomain dom = surf->containingDomain();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      int bd=-1, bd2=-1;
      Vector2D par = Vector2D(parvals[2*ki], parvals[2*ki+1]);
      if (dom.isOnBoundary(par, eps, bd, bd2))
	{
	  if (bd2 >= 0)
	    {
	      // Check best fit
	    }
	  int move_ix = (bd == 0) ? cv_ix[3] :
	    ((bd == 1) ? cv_ix[1] : ((bd == 2) ? cv_ix[0] : cv_ix[2]));
	  move2adj[move_ix].push_back(points[ki]);
	}
      else
	remain.push_back(points[ki]);
    }
}

//===========================================================================
bool RevEng::suitcaseCorner(vector<RevEngRegion*>& adj_blends,
			    RevEngEdge* rev_edge)
//===========================================================================
{
  vector<vector<shared_ptr<ftEdge> > > trim_edgs(adj_blends.size());
  vector<vector<pair<double, double> > > par_lim(adj_blends.size());
  for (size_t ki=0; ki<adj_blends.size(); ++ki)
    {
      trim_edgs[ki] = adj_blends[ki]->getTrimEdges();
      par_lim[ki].resize(trim_edgs[ki].size());
      for (size_t kj=0; kj<par_lim[ki].size(); ++kj)
	par_lim[ki][kj] = std::make_pair(-1.0, -1.0);   // Dummy
    }

  vector<vector<size_t> > match;
  bool found1 = getBlendRegMatches(adj_blends, trim_edgs, par_lim, match);
  if (!found1)
    return false;

  vector<vector<pair<double, double> > > midp(par_lim.size());
  bool found2 = getTrimEdgeMidpoint(trim_edgs, par_lim, midp);
  if (!found2)
    return false;

  // Compute boundary curves
  vector<shared_ptr<CurveOnSurface> > blend_bd;
  vector<RevEngRegion*> adjreg;
  bool found3 = computeCornerBoundaryCurves(adj_blends, trim_edgs, par_lim,
					    midp, match, approx_tol_,
					    blend_bd, adjreg);
  if (!found3)
    return false;
					    

  // Ensure consistent curve sequence and direction
  vector<shared_ptr<CurveOnSurface> > blend_bd0(blend_bd.begin(), blend_bd.end());
  RevEngUtils::setLoopSeq(blend_bd);
  
#ifdef DEBUG_BLEND
  std::ofstream of4("space_bd2.g2");
  for (size_t kr=0; kr<blend_bd.size(); ++kr)
    {
      shared_ptr<ParamCurve> space = blend_bd[kr]->spaceCurve();
      space->writeStandardHeader(of4);
      space->write(of4);
      of4 << "400 1 0 4 255 0 0 255" << std::endl;
      of4 << "1" << std::endl;
      of4 << space->point(space->startparam()) << std::endl;
    }
#endif

  // Ensure non-rational spline boundary curves and extract cross parameter curves
  double tol1 = 0.5*approx_tol_;
  vector<shared_ptr<ParamCurve> > crosscvs;
  vector<shared_ptr<ParamCurve> > bdcvs;
  bool found4 = getCoonsBoundaryInfo(blend_bd, tol1, bdcvs, crosscvs);
  if (!found4)
    return false;

#ifdef DEBUG_BLEND
  std::ofstream of4_2("space_bd3.g2");
  for (size_t kr=0; kr<bdcvs.size(); ++kr)
    {
      bdcvs[kr]->writeStandardHeader(of4_2);
      bdcvs[kr]->write(of4_2);
    }
#endif
  
  // Create surface
  shared_ptr<SplineSurface> coons;
  if (bdcvs.size() == 4)
    {
      coons =
	shared_ptr<SplineSurface>(CoonsPatchGen::createCoonsPatch(bdcvs,
								  crosscvs,
								  tol1,
								  anglim_));
      RevEngUtils::smoothSurf(coons, 2);
    }
  else
    {
      MESSAGE("RevEng::suitcaseCorner. Only 4 boundary curves are supported.");
      return false;
    }

  if (!coons.get())
    return false;
  
#ifdef DEBUG_BLEND
  std::ofstream of5("coons_patch.g2");
  coons->writeStandardHeader(of5);
  coons->write(of5);
#endif

  // Set trimming curves
  double eps = 1.0e-9;
  CurveLoop cvloop = SurfaceTools::outerBoundarySfLoop(coons, eps);
  vector<shared_ptr<ParamCurve> > loopcvs = cvloop.getCurves();

  // Define match between boundary curves of corner surface and
  // adjacent surfaces
  vector<int> cv_ix(loopcvs.size());
  for (size_t ki=0; ki<loopcvs.size(); ++ki)
    {
      Point mid = loopcvs[ki]->point(0.5*(loopcvs[ki]->startparam()+
					  loopcvs[ki]->endparam()));
      int ix = -1;
      double mindist = std::numeric_limits<double>::max();
      for (size_t kj=0; kj<blend_bd0.size(); ++kj)
	{
	  double tpar, dist;
	  Point close;
	  blend_bd0[kj]->closestPoint(mid, blend_bd0[kj]->startparam(),
				     blend_bd0[kj]->endparam(), tpar,
				     close, dist);
	  if (dist < mindist)
	    {
	      ix = (int)kj;
	      mindist = dist;
	    }
	}
      cv_ix[ki] = ix;
    }

  // Move points from blend region as appropriate. First fetch blend points
    vector<RevEngRegion*> blend_regs;
  rev_edge->getAllBlendRegs(blend_regs);
  vector<RevEngPoint*> blend_pts;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_pts.insert(blend_pts.end(), blend_regs[ki]->pointsBegin(),
		     blend_regs[ki]->pointsEnd());

  // double angtol = 5.0*anglim_;
  // double maxd, avd;
  // int num_in, num2_in;
  // vector<double> parvals;
  // vector<pair<double,double> > dist_ang;
  // vector<RevEngPoint*> inpt, outpt;
  // RevEngUtils::distToSurf(blend_pts.begin(), blend_pts.end(), coons,
  // 			  approx_tol_, maxd, avd, num_in, num2_in, inpt, outpt,
  // 			  parvals, dist_ang, angtol);

  double angtol = 5.0*anglim_;
  vector<vector<RevEngPoint*> > move2adj(4);
  vector<RevEngPoint*> remain;
  extractOutPoints(blend_pts, coons, cv_ix, approx_tol_, angtol,
		   move2adj, remain);
  // RectDomain dom = coons->containingDomain();
  // for (size_t ki=0; ki<blend_pts.size(); ++ki)
  //   {
  //     int bd=-1, bd2=-1;
  //     Vector2D par = Vector2D(parvals[2*ki], parvals[2*ki+1]);
  //     if (dom.isOnBoundary(par, eps, bd, bd2))
  // 	{
  // 	  if (bd2 >= 0)
  // 	    {
  // 	      // Check best fit
  // 	    }
  // 	  int move_ix = (bd == 0) ? cv_ix[3] :
  // 	    ((bd == 1) ? cv_ix[1] : ((bd == 2) ? cv_ix[0] : cv_ix[2]));
  // 	  move2adj[move_ix].push_back(blend_pts[ki]);
  // 	}
  //     else
  // 	remain.push_back(blend_pts[ki]);
  //   }
  // if (remain.size() == 0 && remain.size() < blend_pts.size())
  //   return false;

  bool OK;
  for (size_t ki=0; ki<move2adj.size(); ++ki)
    if (move2adj[ki].size() > 0)
      OK = adjreg[ki]->addPointsToGroup(move2adj[ki], approx_tol_, angtol);

  // Delete current blend regions
  vector<HedgeSurface*> prev_sfs;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      //blend_regs[ki]->removeFromAdjacent();
      blend_regs[ki]->setRemove();
    }
  
  // Create blend region
  shared_ptr<RevEngRegion> blendreg(new RevEngRegion(classification_type_,
						     edge_class_type_,
						     remain));  
  blendreg->setRegionAdjacency();
  regions_.push_back(blendreg);
  shared_ptr<HedgeSurface> hedge;
  shared_ptr<ParamSurface> blend_tmp = coons;
  blendreg->setAssociatedSurface(blend_tmp, approx_tol_, angtol,
				 min_point_region_, hedge);
  rev_edge->setBlendRegSurf(blendreg.get());
  blendreg->setBlendEdge(rev_edge);
  if (hedge.get())
    surfaces_.push_back(hedge);
  
  // Add edges to regions
  RectDomain dom = coons->containingDomain();
  double umid = 0.5*(dom.umin()+dom.umax());
  double vmid = 0.5*(dom.vmin()+dom.vmax());
  Point mid = coons->ParamSurface::point(umid, vmid);
  for (size_t ki=0; ki<loopcvs.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> other = blend_bd0[cv_ix[ki]];
      RevEngRegion *adj = adjreg[cv_ix[ki]];
      pairOfRegionEdges(blendreg.get(), hedge.get(), loopcvs[ki], other, adj, mid,
			tol1);
    }

  // Move points from adjacent regions to corner blend region as appropriate
  for (size_t ki=0; ki<adjreg.size(); ++ki)
    {
      blendreg->growInDomain(adjreg[ki], approx_tol_, angtol);
      if (adjreg[ki]->numPoints() == 0)
	adjreg[ki]->setRemove();
    }
  
  // Remove obsolete edge information
  RevEngRegion* adj[2];
  rev_edge->getAdjacent(adj[0], adj[1]);
  size_t kr;
  for (int ka=0; ka<2; ++ka)
    {
      for (kr=0; kr<adj_blends.size(); ++kr)
  	if (adj_blends[kr] == adj[ka])
  	  break;
      if (kr == adj_blends.size())
  	{
  	  rev_edge->eraseAdjacent(adj[ka]);
  	  break;
  	}
    }
  rev_edge->eraseCurves();
  
  return true;
}

//===========================================================================
bool RevEng::createBlendSurface(int ix)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  double diag = bbox_.low().dist(bbox_.high());
  double eps = 1.0e-6;
  
  // Adjacent regions
  RevEngRegion* adj[2];
  edges_[ix]->getAdjacent(adj[0], adj[1]);

  if ((!adj[0]) || (!adj[1]))
    return false;
  if ((!adj[0]->hasSurface()) || (!adj[1]->hasSurface()))
    return false;  // Something wrong

#ifdef DEBUG_BLEND
  std::ofstream of1("adj_groups.g2");
  adj[0]->writeRegionPoints(of1);
  adj[1]->writeRegionPoints(of1);
#endif

  if (adj[0]->hasAssociatedBlend() || adj[1]->hasAssociatedBlend())
    return false; 
  // Regions in blend area
  vector<RevEngRegion*> blend_regs;
  edges_[ix]->getAllBlendRegs(blend_regs);
  
  // Intersection curve
  double eps1 = approx_tol_; //std::max(1.0e-4, 0.1*approx_tol_);
  edges_[ix]->fixMismatchCurves(eps1);
  vector<shared_ptr<CurveOnSurface> > cvs;
  edges_[ix]->getCurve(cvs, true);
  // for (size_t ki=0; ki<cvs.size(); ++ki)
  //   cvs[ki]->fixMismatchCurves(eps1);
  // for (size_t ki=0; ki<cvs2.size(); ++ki)
  //   cvs2[ki]->fixMismatchCurves(eps1);

  // TEST. Check if the intersection curve must be updated
  shared_ptr<ParamSurface> surf1 = adj[0]->getSurface(0)->surface();
  if (!surf1->isBounded())
    adj[0]->getSurface(0)->limitSurf(diag);
  shared_ptr<ParamSurface> surf2 = adj[1]->getSurface(0)->surface();
  if (!surf2->isBounded())
    adj[1]->getSurface(0)->limitSurf(diag);
  shared_ptr<BoundedSurface> bd1, bd2;
  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
  BoundedUtils::getSurfaceIntersections(surf1, surf2, int_tol_,
					int_cvs1, bd1,
					int_cvs2, bd2);
#ifdef DEBUG_BLEND
  std::ofstream of_int("intcurves_edge.g2");
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = cvs[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of_int);
      tmp_cv->write(of_int);
    }
  for (size_t ki=0; ki<int_cvs1.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = int_cvs1[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of_int);
      tmp_cv->write(of_int);
    }
#endif
  //  #if 0
  if (edges_[ix]->isClosed(approx_tol_))
    edges_[ix]->replaceCurves(int_cvs1, int_cvs2);  // The seam might have moved
      
  else //if (int_cvs1.size() == cvs.size())
    {
      vector<shared_ptr<CurveOnSurface> > cvs2;
      edges_[ix]->getCurve(cvs2, false);
      if (cvs[0]->isClosed())
	{
	  cvs[0] = int_cvs1[0];
	  cvs2[0] = int_cvs2[0];
	}
      else
	{
	  size_t st = cvs.size() - 1;
	  Point pos1, pos2, pos3, pos4;
	  double tdel1 = cvs[0]->endparam() - cvs[0]->startparam();
	  double tdel2 = cvs[st]->endparam() - cvs[st]->startparam();
	  cvs[0]->point(pos1, cvs[0]->startparam());
	  cvs[0]->point(pos3, cvs[0]->startparam() + 0.1*tdel1);
	  cvs[st]->point(pos2, cvs[st]->endparam());
	  cvs[st]->point(pos4, cvs[st]->endparam() - 0.1*tdel2);
	  double tp1, tp2;
	  double td1 = std::numeric_limits<double>::max();
	  double td2 = std::numeric_limits<double>::max();
	  double td3 = std::numeric_limits<double>::max();
	  double td4 = std::numeric_limits<double>::max();
	  int ix1 = -1, ix2 = -1;
	  for (size_t kr=0; kr<int_cvs1.size(); ++kr)
	    {
	      double tp1_2, tp2_2, tp1_3, tp2_3, td1_2, td2_2, td1_3, td2_3;
	      Point cl1_2, cl2_2;
	      double tmin = int_cvs1[kr]->startparam();
	      double tmax = int_cvs1[kr]->endparam();
	      int_cvs1[kr]->closestPoint(pos1, tmin, tmax, tp1_2, cl1_2, td1_2);
	      int_cvs1[kr]->closestPoint(pos2, tmin, tmax, tp2_2, cl2_2, td2_2);
	      int_cvs1[kr]->closestPoint(pos3, tmin, tmax, tp1_3, cl1_2, td1_3);
	      int_cvs1[kr]->closestPoint(pos4, tmin, tmax, tp2_3, cl2_2, td2_3);
	      if (td1_2 < td1-eps || (td1_2 <= td1+eps && td1_3 < td3))
		{
		  ix1 = (int)kr;
		  tp1 = tp1_2;
		  td1 = td1_2;
		  td3 = td1_3;
		}
	      if (td2_2 < td2-eps || (td2_2 <= td2+eps && td2_3 < td4))
		{
		  ix2 = (int)kr;
		  tp2 = tp2_2;
		  td2 = td2_2;
		  td4 = td2_3;
		}
	    }
	  
	  if (ix2 - ix1 + 1 == (int)cvs.size() && tp1 < tp2)
	    {
	      // if (tp2 < tp1)
	      // 	{
	      // This is not expected. Should check if one of the closest points
	      // is found on the wrong side of a seam

		  // for (size_t kr=0; kr<int_cvs1.size(); ++kr)
		  //   {
		  //     int_cvs1[kr]->reverseParameterDirection();
		  //     int_cvs2[kr]->reverseParameterDirection();
		  //   }
		  // for (size_t kr=0; kr<int_cvs1.size()/2; ++kr)
		  //   {
		  //     std::swap(int_cvs1[kr], int_cvs1[int_cvs1.size()-kr-1]);
		  //     std::swap(int_cvs2[kr], int_cvs2[int_cvs2.size()-kr-1]);
		  //   }
		  // std::swap(ix1, ix2);
		  // Point cl1, cl2;
		  // int_cvs1[ix1]->closestPoint(pos1, int_cvs1[ix1]->startparam(),
		  // 			      int_cvs1[ix1]->endparam(), tp1, cl1, td1);
		  // int_cvs1[ix2]->closestPoint(pos2, int_cvs1[ix2]->startparam(),
		  // 			      int_cvs1[ix2]->endparam(), tp2, cl2, td2);
		// }

	      for (size_t kr=0; kr<int_cvs1.size(); )
		{
		  double tp3 = std::max(tp1, int_cvs1[kr]->startparam());
		  double tp4 = std::min(tp2, int_cvs1[kr]->endparam());
		  if (ix1 <= ix2 && ((int)kr < ix1 || (int)kr > ix2))
		    {
		      int_cvs1.erase(int_cvs1.begin()+kr);
		      int_cvs2.erase(int_cvs2.begin()+kr);
		      if ((int)kr < ix1)
			ix1--;
		      if ((int)kr < ix2)
			ix2--;
		    }
		  else if (tp4 > tp3 && tp3 < int_cvs1[kr]->endparam() &&
			   tp4 > int_cvs1[kr]->startparam() && 
			   (tp3 > int_cvs1[kr]->startparam() || tp4 < int_cvs1[kr]->endparam()))
		    {
		      shared_ptr<CurveOnSurface> sub1(int_cvs1[kr]->subCurve(tp3, tp4));
		      cvs[kr] = sub1;
		      shared_ptr<CurveOnSurface> sub2(int_cvs2[kr]->subCurve(tp3, tp4));
		      cvs2[kr] = sub2;
#ifdef DEBUG_BLEND
		      shared_ptr<ParamCurve> tmp_cv = sub1->spaceCurve();
		      tmp_cv->writeStandardHeader(of_int);
		      tmp_cv->write(of_int);
#endif
		      ++kr;
		    }
		  else
		    {
		      cvs[kr] = int_cvs1[kr];
		      cvs2[kr] = int_cvs2[kr];
		      ++kr;
		    }
		}
	    }
	}
      edges_[ix]->replaceCurves(cvs, cvs2);
    }
  
  cvs.clear();
  edges_[ix]->getCurve(cvs, true);
  //  #endif
  vector<Point> der(2);
  cvs[0]->point(der, 0.5*(cvs[0]->startparam()+cvs[0]->endparam()), 1);

#ifdef DEBUG_BLEND
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = cvs[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of1);
      tmp_cv->write(of1);
    }
#endif

  // Width
  double width = edges_[ix]->getDistance();

  bool out1 = false, out2 = false;
  edges_[ix]->getOuterInfo(out1, out2);

  ClassType classtype1 = surf1->instanceType();
  ClassType classtype2 = surf2->instanceType();
  if (classtype1 != Class_Plane && classtype1 != Class_Cylinder &&
      classtype1 != Class_Cone)
    return false;
  if (classtype2 != Class_Plane && classtype2 != Class_Cylinder &&
      classtype2 != Class_Cone)
    return false;
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  vector<shared_ptr<ElementarySurface> > elem_sfs(2);
  elem_sfs[0] = elem1;
  elem_sfs[1] = elem2;
  if (!(elem1.get() && elem2.get()))
    return false;
  Point dir1 = elem1->direction();
  Point norm1 = adj[0]->getMeanNormalTriang();
  if (elem1->instanceType() == Class_Plane && dir1*norm1 < 0.0)
    dir1 *= -1;
  Point dir2 = elem2->direction();
  Point norm2 = adj[1]->getMeanNormalTriang();
  if (elem2->instanceType() == Class_Plane && dir2*norm2 < 0.0)
    dir2 *= -1;

  if (classtype1 == Class_Plane && classtype2 == Class_Plane)
    {
      // Create cylinder
    }
  else if (classtype1 == Class_Plane || classtype2 == Class_Plane)
    {
      // Create torus
    }
  else if ((classtype1 == Class_Cylinder && classtype2 == Class_Cone) ||
	   (classtype2 == Class_Cylinder && classtype1 == Class_Cone))
    {
      // Create torus
    }
  else
    return false;  // Not supported

  // Collect points from adjacent surfaces
  vector<vector<RevEngPoint*> > blend_pts(3);
  double tmin = cvs[0]->startparam();
  double tmax = cvs[0]->endparam();
  for (int ka=0; ka<2; ++ka)
    adj[ka]->getNearPoints(cvs[0], tmin, tmax, 2*width, angtol, blend_pts[ka]);
 
  // Collect points from associated blend regions
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    blend_pts[2].insert(blend_pts[2].end(), blend_regs[ki]->pointsBegin(),
			blend_regs[ki]->pointsEnd());

#ifdef DEBUG_BLEND
  std::ofstream of2("blend_pts.g2");
  for (size_t ki=0; ki<3; ++ki)
    {
      if (blend_pts[ki].size() > 0)
	{
	  of2 << "400 1 0 4 0 255 0 255" << std::endl;
	  of2 << blend_pts[ki].size() << std::endl;
	  for (size_t kr=0; kr<blend_pts[ki].size(); ++kr)
	    of2 << blend_pts[ki][kr]->getPoint() << std::endl;
	}
    }
#endif

  shared_ptr<ElementarySurface> blend_surf;
  double ylen = 0.0;
  double axis_ang = dir1.angle(dir2);
  int ldir = -1;
  if ((classtype1 == Class_Plane && classtype2 == Class_Plane) ||
      ((classtype1 == Class_Plane || classtype2 == Class_Plane) &&
       fabs(0.5*M_PI - axis_ang) <= angtol))
    {
      // Create cylinder
      Point lin1, lin2;
      Point dir1_2 = dir1, dir2_2 = dir2;
      if (classtype1 == Class_Plane)
	lin1 = der[1].cross(dir1);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf1->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf1->normal(dir1_2, clo_u, clo_v);
	  lin1 = der[1].cross(dir1_2);
	}
      
      if (classtype2 == Class_Plane)
	lin2 = der[1].cross(dir2);
      else
	{
	  double clo_u, clo_v, clo_dist;
	  Point clo;
	  surf2->closestPoint(der[0], clo_u, clo_v, clo, clo_dist, eps);
	  surf2->normal(dir2_2, clo_u, clo_v);
	  lin2 = der[1].cross(dir2_2);
	}
      
      lin1.normalize();
      if (lin1*dir2_2 < 0.0)
	lin1 *= -1;
      lin2.normalize();
      if (lin2*dir1_2 < 0.0)
	lin2 *= -1;
      int sgn = (out1 && out2) ? 1 : -1;
      blend_surf = createCylinderBlend(blend_pts, width, der[0],
				       der[1], lin1, lin2, sgn);
      double radius = blend_surf->radius(0.0, 0.0);
      //double ang = dir1.angle(dir2);
      Point centre = blend_surf->location();
      double xlen = der[0].dist(centre);
      ylen = sqrt(xlen*xlen - radius*radius);
      
      ldir = 1;
    }
  else
    {
      int sgn = 1;
      if ((classtype1 == Class_Plane && dir1*elem1->direction() < 0.0) ||
	  (classtype2 == Class_Plane && dir2*elem2->direction() < 0.0))
	sgn = -1;

      blend_surf = torusBlend(blend_pts, cvs[0], elem1, elem2, width,
			      out1, out2, sgn);

      double Rrad = blend_surf->radius(0.0, 0.0);
      Point centre = blend_surf->location();
      double xlen = der[0].dist(centre);
      ylen = fabs(xlen - Rrad);

      ldir = 1;
    }
  
  if (!blend_surf.get())
    {
#ifdef DEBUG_BLEND
      std::cout << "No blend_surf" << std::endl;
#endif
      return false;
    }

#ifdef DEBUG_BLEND
  for (int ka=0; ka<2; ++ka)
    {
      int num = adj[ka]->numPoints();
      for (int kb=0; kb<num; ++kb)
	if (adj[ka]->getPoint(kb)->region() != adj[ka])
	  std::cout << "Inconsistent region pointers, pre sortBlendPoints: " << ka << " " << kb << std::endl;
    }
#endif
  
  // Associate blend points with the appropriate surface (region)
  // First remove blend points from the adjacent surfaces
  vector<vector<RevEngPoint*> > out_pts(2);
  for (int ka=0; ka<2; ++ka)
    {
      adj[ka]->sortBlendPoints(blend_pts[ka], cvs, ylen, true, out_pts[ka]);
      blend_pts[2].insert(blend_pts[2].end(), out_pts[ka].begin(),
			  out_pts[ka].end());
    }
#ifdef DEBUG_BLEND
  for (int ka=0; ka<2; ++ka)
    {
      int num = adj[ka]->numPoints();
      for (int kb=0; kb<num; ++kb)
	if (adj[ka]->getPoint(kb)->region() != adj[ka])
	  std::cout << "Inconsistent region pointers, post sortBlendPoints1: " << ka << " " << kb << std::endl;
    }
#endif
      
  // Move blend points to the adjacent surfaces if appropriate
  vector<vector<RevEngPoint*> > in_pts(2);
  adj[0]->sortBlendPoints(blend_pts[2], cvs, ylen, adj[1],
			  in_pts[0], in_pts[1]);

#ifdef DEBUG_BLEND
  for (int ka=0; ka<2; ++ka)
    {
      int num = adj[ka]->numPoints();
      for (int kb=0; kb<num; ++kb)
	if (adj[ka]->getPoint(kb)->region() != adj[ka])
	  std::cout << "Inconsistent region pointers, post sortBlendPoints2: " << ka << " " << kb << std::endl;
    }
#endif
      
  // Update adjacent surfaces with modified collection of points
  for (int ka=0; ka<2; ++ka)
    {
      adj[ka]->updateWithPointsInOut(out_pts[ka], in_pts[ka], approx_tol_, angtol);
    }

  // Delete current blend regions
  vector<HedgeSurface*> prev_sfs;
  for (size_t ki=0; ki<blend_regs.size(); ++ki)
    {
      blend_regs[ki]->removeFromAdjacent();
      int num_sfs = blend_regs[ki]->numSurface();
      for (int ka=0; ka<num_sfs; ++ka)
	prev_sfs.push_back(blend_regs[ki]->getSurface(ka));
    }
  edges_[ix]->clearBlendRegions();
  if (blend_regs.size() > 0)
    {
      int dummy_ix = 0;
      updateRegionsAndSurfaces(dummy_ix, blend_regs, prev_sfs);
    }
  
  if (blend_pts[2].size() == 0)
    {
#ifdef DEBUG_BLEND
      std::cout << "No points for blend_surf" << std::endl;
#endif
      return false;
    }
  
  // Define blend region
  shared_ptr<RevEngRegion> blendreg(new RevEngRegion(classification_type_,
						     edge_class_type_,
						     blend_pts[2]));
#ifdef DEBUG_BLEND
  std::ofstream of3("updated_blend_pts.g2");
  adj[0]->writeRegionPoints(of3);
  adj[1]->writeRegionPoints(of3);
  blendreg->writeRegionPoints(of3);
#endif
  
  blendreg->setRegionAdjacency();
  regions_.push_back(blendreg);
  shared_ptr<HedgeSurface> hedge;
  shared_ptr<ParamSurface> blend_surf_tmp = blend_surf;
  blendreg->setAssociatedSurface(blend_surf_tmp, approx_tol_, angtol,
				 min_point_region_, hedge);
  if (hedge.get())
    surfaces_.push_back(hedge);

  // // Check blend points
  // vector<vector<RevEngPoint*> > out_groups;
  // blendreg->extractOutPoints(ldir, approx_tol_, angtol,
  // 			     1.1*angtol, out_groups);
  // if (out_groups.size() > 0)
  //   {
  //     vector<HedgeSurface*> dummy_sfs;
  //     surfaceExtractOutput(regions_.size()-1, out_groups, dummy_sfs);
  //   }
   
  // Update edge with blend region (surface)
  edges_[ix]->setBlendRegSurf(blendreg.get());
  blendreg->setBlendEdge(edges_[ix].get());
  for (int ka=0; ka<2; ++ka)
    adj[ka]->updateRegionAdjacency();
  blendreg->setRegionAdjacency();
  //edges_[ix]->setAltRadius(radius1);

  for (int ka=0; ka<2; ++ka)
    {
      vector<vector<RevEngPoint*> > sep_groups;
      adj[ka]->splitRegion(sep_groups);
      if (sep_groups.size() > 0)
	{
	  size_t kh;
	  for (kh=0; kh<regions_.size(); ++kh)
	    if (regions_[kh].get() == adj[ka])
	      break;
	  vector<HedgeSurface*> dummy_sfs;
	  surfaceExtractOutput((kh<regions_.size()) ? (int)kh : 0, sep_groups,
			       dummy_sfs);
	}
   }

  return true;
}

//===========================================================================
double
RevEng::computeTorusRadius(vector<vector<RevEngPoint*> >& blend_pts,
			   shared_ptr<CurveOnSurface>& cv,
			   const Point& locp, const Point& normal,
			   shared_ptr<ElementarySurface> rotational,
			   double width, bool plane_out, bool rot_out)
//===========================================================================
{
  double alpha = 0;
  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ElementarySurface>(rotational);
  if (cone.get())
    alpha = cone->getConeAngle();
  Point axis = rotational->direction();
  
  // Compute minor radius of torus
  double lrad = 0.0;
  int lnmb = 0;
  double beta = 0.5*M_PI + alpha;
  //double phi = 0.5*M_PI - alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist;
	  Point close;
	  cv->closestPoint(ptpos, cv->startparam(), cv->endparam(),
				tpar, close, dist);

	  // Define line in the point between the adjacent surfaces
	  vector<Point> der(2);
	  cv->point(der, tpar, 1);
	  der[1].normalize();
	  Point dir1 = der[1].cross(normal);
	  Point vec1 = der[0] - centr;
	  if ((dir1*vec1 < 0.0 && (!rot_out)) || (dir1*vec1 > 0.0 && rot_out))
	    dir1 *= -1;
	  Point dir2 = (plane_out) ? normal : -normal;
	  Point dir3;
	  if (alpha > 0)
	    {
	      Matrix3D mat;
	      mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
	      Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
	      Vector3D dir3_2 = mat*dir2_2;
	      dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
	    }
	  else
	    dir3 = dir2;
	  Point dir4 = 0.5*(dir1 + dir3);
#ifdef DEBUG_BLEND
	  std::ofstream of("midlin.g2");
	  of << "410 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir1 << std::endl;
	  of << "410 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir4 << std::endl;
	  of << "410 1 0 4 0 0 255 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir2 << std::endl;
	      
#endif
	  shared_ptr<Line> line(new Line(der[0], dir4));
	  line->setParameterInterval(-2*width, 2*width);
	  double lpar, ldist;
	  Point lclose;
	  line->closestPoint(ptpos, line->startparam(), line->endparam(),
			     lpar, lclose, ldist);
	  if (ldist <= approx_tol_)
	    {
	      double dd2 = close.dist(lclose);
	      double xlen = dd2/div;
	      lrad += xlen;
	      lnmb++;
	    }
	  int stop_break = 1;
	}
    }

  if (lnmb > 0)
    lrad /= (double)lnmb;

  return lrad;
}

//===========================================================================
void RevEng::getTorusParameters(shared_ptr<ElementarySurface> planar,
				shared_ptr<ElementarySurface> rotational,
				double radius, int sgn1, int sgn2, double& Rrad, 
				Point& centre, Point& normal, Point& Cx)
//===========================================================================
{
  Point locp = planar->location();
  normal = planar->direction();
  Point loc = rotational->location();
  Point axis = rotational->direction();
  double rd = (locp - loc)*axis;
  Point centre0 = loc + rd*axis;
  centre = centre0 - sgn1*radius*normal;
  Rrad = rotational->radius(0.0, rd);
  double alpha = 0.0;
  if (rotational->instanceType() == Class_Cone)
    alpha = ((Cone*)(rotational.get()))->getConeAngle();
  double phi = 0.5*M_PI - alpha;
  double sd = radius/sin(phi);
  Cx = rotational->direction2();
  Rrad += (sgn2*sd);
}

//===========================================================================
shared_ptr<Torus>
RevEng::torusBlend(vector<vector<RevEngPoint*> >& blend_pts,
		   vector<shared_ptr<CurveOnSurface> >& cvs,
		   const Point& locp, const Point& normal,
		   shared_ptr<ElementarySurface> rotational,
		   double width, bool plane_out, bool rot_out)
//===========================================================================
{
  shared_ptr<Torus> torus;

  double alpha = 0;
  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ElementarySurface>(rotational);
  if (cone.get())
    alpha = cone->getConeAngle();
  Point axis = rotational->direction();
  
  // Compute minor radius of torus
  double lrad = 0.0;
  double d2 = 0.0;
  int lnmb = 0;
  double beta = 0.5*M_PI + alpha;
  double phi = 0.5*M_PI - alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist=std::numeric_limits<double>::max();
	  Point close;
	  int ix = -1;
	  for (size_t kr=0; kr<cvs.size(); ++kr)
	    {
	      double tpar0, dist0;
	      Point close0;
	      cvs[kr]->closestPoint(ptpos, cvs[kr]->startparam(), cvs[kr]->endparam(),
				tpar0, close0, dist0);
	      if (dist0 < dist)
		{
		  tpar = tpar0;
		  dist = dist0;
		  close = close0;
		  ix = (int)kr;
		}
	      if (ix < 0)
		continue;

	      // Define line in the point between the adjacent surfaces
	      vector<Point> der(2);
	      cvs[ix]->point(der, tpar, 1);
	      der[1].normalize();
	      Point dir1 = der[1].cross(normal);
	      Point vec1 = der[0] - centr;
	      if ((dir1*vec1 < 0.0 && (!rot_out)) || (dir1*vec1 > 0.0 && rot_out))
		dir1 *= -1;
	      Point dir2 = (plane_out) ? normal : -normal;
	      Point dir3;
	      if (alpha > 0)
		{
		  Matrix3D mat;
		  mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
		  Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
		  Vector3D dir3_2 = mat*dir2_2;
		  dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
		}
	      else
		dir3 = dir2;
	      Point dir4 = 0.5*(dir1 + dir3);
#ifdef DEBUG_BLEND
	      std::ofstream of("midlin.g2");
	      of << "410 1 0 4 0 255 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << der[0] << " " << der[0]+dir1 << std::endl;
	      of << "410 1 0 4 255 0 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << der[0] << " " << der[0]+dir4 << std::endl;
	      of << "410 1 0 4 0 0 255 255" << std::endl;
	      of << "1" << std::endl;
	      of << der[0] << " " << der[0]+dir2 << std::endl;
	      
#endif
	      shared_ptr<Line> line(new Line(der[0], dir4));
	      line->setParameterInterval(-2*width, 2*width);
	      double lpar, ldist;
	      Point lclose;
	      line->closestPoint(ptpos, line->startparam(), line->endparam(),
				 lpar, lclose, ldist);
	      if (ldist <= approx_tol_)
		{
		  double dd2 = close.dist(lclose);
		  double xlen = dd2/div;
		  lrad += xlen;
		  d2 += dd2;
		  lnmb++;
		}
	      int stop_break = 1;
	    }
	}
    }

  if (lnmb > 0)
    {
      lrad /= (double)lnmb;
      d2 /= (double)lnmb;
    }

  int sgn1 = (plane_out) ? -1 : 1;
  int sgn2 = rot_out ? -1 : 1;
  Point pos = centr - sgn1*lrad*normal;
  double Rrad = rotational->radius(0.0, rd);
  double sd = lrad/sin(phi);
  Point Cx = rotational->direction2();
  torus = shared_ptr<Torus>(new Torus(Rrad+sgn2*sd, lrad, pos, normal, Cx));
  if (rot_out)
    {
      RectDomain dom = torus->getParameterBounds();
      try {
	torus->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
				  dom.umax(), dom.vmax()-M_PI);
      }
      catch (...)
	{
	}
    }
  
#ifdef DEBUG_BLEND
  std::ofstream of2("tor_blend.g2");
  torus->writeStandardHeader(of2);
  torus->write(of2);
  RectDomain dom2 = torus->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = torus->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of2);
  tmp_cvs[0]->write(of2);
  double tf = (sgn2 == 1) ? 0.5 : 1.5;
  vector<shared_ptr<ParamCurve> > tmp_cvs2 =
    torus->constParamCurves(dom2.vmin()+tf*M_PI, true);
  tmp_cvs2[0]->writeStandardHeader(of2);
  tmp_cvs2[0]->write(of2);
#endif
  return torus;
}

//===========================================================================
double
RevEng::computeTorusRadius(vector<vector<RevEngPoint*> >& blend_pts,
			   shared_ptr<CurveOnSurface>& cv,
			   shared_ptr<ElementarySurface> elem1,
			   shared_ptr<ElementarySurface> elem2,
			   double width, bool out1, bool out2, int sgn,
			   double& d2)
//===========================================================================
{
  d2 = 0.0;
  double alpha1 = 0.0, alpha2 = 0.0;
  shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
  if (cone1.get())
    alpha1 = cone1->getConeAngle();
  shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
  if (cone2.get())
    alpha2 = cone2->getConeAngle();
  if (cone1.get() && cone2.get())
    return 0.0;   // Two cones are not handled
  double alpha = fabs(alpha1) + fabs(alpha2);
  shared_ptr<Plane> plane;
  if (elem1->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem1);
  else if (elem2->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem2);
  int state = (plane.get()) ? 1 : 2;
  shared_ptr<ElementarySurface> rotational;
  if (elem1->instanceType() == Class_Plane)
    rotational = elem2;
  else if (elem2->instanceType() == Class_Plane)
    rotational = elem1;
  else if (cone1.get())
    rotational = elem2;  // elem2 is expected to be a cylinder
  else if (cone2.get())
    rotational = elem1;
  Point axis = rotational->direction();
  Point normal = (plane.get()) ? plane->direction() : axis;
  if (state == 1)
    normal *= sgn;
  if ((state == 1 && elem2->instanceType() == Class_Plane) ||
      (state == 2 && elem2->instanceType() == Class_Cylinder))
    std::swap(out1,out2);  // Call by value means this swap stays local
  
  // Compute minor radius of torus
  double lrad = 0.0;
  int lnmb = 0;
  double beta = (plane.get()) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
  double phi = beta - 2.0*alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  Point locp = (plane.get()) ? plane->location() :
    cv->ParamCurve::point(0.5*(cv->startparam()+cv->endparam()));;
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist;
	  Point close;
	  cv->closestPoint(ptpos, cv->startparam(), cv->endparam(),
			   tpar, close, dist);

	  // Define line in the point between the adjacent surfaces
	  vector<Point> der(2);
	  cv->point(der, tpar, 1);
	  der[1].normalize();
	  Point dir1 = der[1].cross(normal);
	  Point vec1 = der[0] - centr;
	  if ((dir1*vec1 < 0.0 && (!out2)) || (dir1*vec1 > 0.0 && out2))
	    dir1 *= -1;
	  Point dir2 = (out1) ? normal : -normal;
	  Point dir3;
	  if (alpha > 0)
	    {
	      Matrix3D mat;
	      if (state == 1)
		mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
	      else
		mat.setToRotation(0.5*beta, der[1][0], der[1][1], der[1][2]);
	      Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
	      Vector3D dir3_2 = mat*dir2_2;
	      dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
	    }
	  else
	    dir3 = dir2;
	  Point dir4 = (state == 1 )? 0.5*(dir1 + dir3) : dir3;
	  if (state == 2 && (!out2))
	    dir4 *= -1;
#ifdef DEBUG_BLEND
	  std::ofstream of("midlin.g2");
	  of << "410 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir1 << std::endl;
	  of << "410 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir4 << std::endl;
	  of << "410 1 0 4 0 0 255 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir2 << std::endl;
	  of << "410 1 0 4 100 100 55 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir3 << std::endl;
	      
#endif
	  shared_ptr<Line> line(new Line(der[0], dir4));
	  line->setParameterInterval(-2*width, 2*width);
	  double lpar, ldist;
	  Point lclose;
	  line->closestPoint(ptpos, line->startparam(), line->endparam(),
			     lpar, lclose, ldist);
	  if (ldist <= approx_tol_)
	    {
	      double dd2 = close.dist(lclose);
	      double xlen = dd2/div;
	      lrad += xlen;
	      d2 += dd2;
	      lnmb++;
	    }
	  int stop_break = 1;
	}
    }


  if (lnmb > 0)
    {
      lrad /= (double)lnmb;
     d2 /= (double)lnmb;
    }

  return lrad;
}

//===========================================================================
void
RevEng::getTorusParameters(shared_ptr<ElementarySurface> elem1,
			   shared_ptr<ElementarySurface> elem2, Point pos,
			   double radius, double d2, bool out1, bool out2, int sgn,
			   double& Rrad, Point& centre, Point& normal, Point& Cx)
//===========================================================================
{
  double alpha1 = 0.0, alpha2 = 0.0;
  shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
  if (cone1.get())
    alpha1 = cone1->getConeAngle();
  shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
  if (cone2.get())
    alpha2 = cone2->getConeAngle();
  if (cone1.get() && cone2.get())
    return;   // Two cones are not handled
  double alpha = fabs(alpha1) + fabs(alpha2);
  shared_ptr<Plane> plane;
  if (elem1->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem1);
  else if (elem2->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem2);
  int state = (plane.get()) ? 1 : 2;
  shared_ptr<ElementarySurface> rotational;
  if (elem1->instanceType() == Class_Plane)
    rotational = elem2;
  else if (elem2->instanceType() == Class_Plane)
    rotational = elem1;
  else if (cone1.get())
    rotational = elem2;  // elem2 is expected to be a cylinder
  else if (cone2.get())
    rotational = elem1;
  Point axis = rotational->direction();
  normal = (plane.get()) ? plane->direction() : axis;
  if (state == 1)
    normal *= sgn;
  Point loc = rotational->location();
  if ((state == 1 && elem2->instanceType() == Class_Plane) ||
      (state == 2 && elem2->instanceType() == Class_Cylinder))
    std::swap(out1,out2);  // Call by value means this swap stays local
  
  double rd = (pos-loc)*axis;
  Point centr = loc + rd*axis;
  double beta = (plane.get()) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
  double phi = beta - 2.0*alpha;
  double phi2 = 0.5*(M_PI - beta);
  int sgn1 = (out1) ? -1 : 1;
  int sgn2 = out2 ? -1 : 1;
  if (state == 2)
    {
      sgn2 *= -1;
    }
  double hh = (state == 1) ? radius : sgn2*(radius + d2)*sin(phi2);
  centre = centr - sgn1*hh*normal;
  Rrad = rotational->radius(0.0, rd);
  double sd = (state == 1) ? radius/sin(phi) : (radius + d2)*cos(phi2);
  Cx = rotational->direction2();
  Rrad += (sgn2*sd);
  
#ifdef DEBUG_BLEND
  shared_ptr<Torus> torus(new Torus(Rrad, radius, centre, normal, Cx));
  if (sgn2 < 0)
    {
      RectDomain dom = torus->getParameterBounds();
      torus->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
				dom.umax(), dom.vmax()-M_PI);
    }
  std::ofstream of2("tor_blend.g2");
  torus->writeStandardHeader(of2);
  torus->write(of2);
  RectDomain dom2 = torus->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = torus->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of2);
  tmp_cvs[0]->write(of2);
  double tf = (sgn2 == 1) ? 0.5 : 1.5;
  vector<shared_ptr<ParamCurve> > tmp_cvs2 =
    torus->constParamCurves(dom2.vmin()+tf*M_PI, true);
  tmp_cvs2[0]->writeStandardHeader(of2);
  tmp_cvs2[0]->write(of2);
#endif
}

//===========================================================================
shared_ptr<Torus>
RevEng::torusBlend(vector<vector<RevEngPoint*> >& blend_pts,
		   shared_ptr<CurveOnSurface>& cv,
		   shared_ptr<ElementarySurface> elem1,
		   shared_ptr<ElementarySurface> elem2,
		   double width, bool out1, bool out2, int sgn)
//===========================================================================
{
  shared_ptr<Torus> torus;

  double alpha1 = 0.0, alpha2 = 0.0;
  shared_ptr<Cone> cone1 = dynamic_pointer_cast<Cone,ElementarySurface>(elem1);
  if (cone1.get())
    alpha1 = cone1->getConeAngle();
  shared_ptr<Cone> cone2 = dynamic_pointer_cast<Cone,ElementarySurface>(elem2);
  if (cone2.get())
    alpha2 = cone2->getConeAngle();
  if (cone1.get() && cone2.get())
    return torus;   // Two cones are not handled
  double alpha = fabs(alpha1) + fabs(alpha2);
  shared_ptr<Plane> plane;
  if (elem1->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem1);
  else if (elem2->instanceType() == Class_Plane)
    plane = dynamic_pointer_cast<Plane,ElementarySurface>(elem2);
  int state = (plane.get()) ? 1 : 2;
  shared_ptr<ElementarySurface> rotational;
  if (elem1->instanceType() == Class_Plane)
    rotational = elem2;
  else if (elem2->instanceType() == Class_Plane)
    rotational = elem1;
  else if (cone1.get())
    rotational = elem2;  // elem2 is expected to be a cylinder
  else if (cone2.get())
    rotational = elem1;
  Point axis = rotational->direction();
  Point normal = (plane.get()) ? plane->direction() : axis;
  if (state == 1)
    normal *= sgn;
  if ((state == 1 && elem2->instanceType() == Class_Plane) ||
      (state == 2 && elem2->instanceType() == Class_Cylinder))
    std::swap(out1,out2);  // Call by value means this swap stays local
  
  // Compute minor radius of torus
  double lrad = 0.0;
  double d2 = 0.0;
  int lnmb = 0;
  double beta = (plane.get()) ? 0.5*M_PI + alpha : M_PI-fabs(alpha);
  double phi = beta - 2.0*alpha;
  double fac = 1.0/sin(0.5*beta);
  double div = fac - 1.0;
  Point loc = rotational->location();
  Point locp = (plane.get()) ? plane->location() :
    cv->ParamCurve::point(0.5*(cv->startparam()+cv->endparam()));;
  double rd = (locp-loc)*axis;
  Point centr = loc + rd*axis;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  // Localize with respect to the guiding curve(s)
	  double tpar, dist;
	  Point close;
	  cv->closestPoint(ptpos, cv->startparam(), cv->endparam(),
			   tpar, close, dist);

	  // Define line in the point between the adjacent surfaces
	  vector<Point> der(2);
	  cv->point(der, tpar, 1);
	  der[1].normalize();
	  Point dir1 = der[1].cross(normal);
	  Point vec1 = der[0] - centr;
	  if ((dir1*vec1 < 0.0 && (!out2)) || (dir1*vec1 > 0.0 && out2))
	    dir1 *= -1;
	  Point dir2 = (out1) ? normal : -normal;
	  Point dir3;
	  if (alpha > 0)
	    {
	      Matrix3D mat;
	      if (state == 1)
		mat.setToRotation(-alpha, der[1][0], der[1][1], der[1][2]);
	      else
		mat.setToRotation(0.5*beta, der[1][0], der[1][1], der[1][2]);
	      Vector3D dir2_2(dir2[0], dir2[1], dir2[2]);
	      Vector3D dir3_2 = mat*dir2_2;
	      dir3 = Point(dir3_2[0], dir3_2[1], dir3_2[2]);
	    }
	  else
	    dir3 = dir2;
	  Point dir4 = (state == 1 )? 0.5*(dir1 + dir3) : dir3;
	  if (state == 2 && (!out2))
	    dir4 *= -1;
#ifdef DEBUG_BLEND
	  std::ofstream of("midlin.g2");
	  of << "410 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir1 << std::endl;
	  of << "410 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir4 << std::endl;
	  of << "410 1 0 4 0 0 255 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir2 << std::endl;
	  of << "410 1 0 4 100 100 55 255" << std::endl;
	  of << "1" << std::endl;
	  of << der[0] << " " << der[0]+dir3 << std::endl;
	      
#endif
	  shared_ptr<Line> line(new Line(der[0], dir4));
	  line->setParameterInterval(-2*width, 2*width);
	  double lpar, ldist;
	  Point lclose;
	  line->closestPoint(ptpos, line->startparam(), line->endparam(),
			     lpar, lclose, ldist);
	  if (ldist <= approx_tol_)
	    {
	      double dd2 = close.dist(lclose);
	      double xlen = dd2/div;
	      lrad += xlen;
	      d2 += dd2;
	      lnmb++;
	    }
	  int stop_break = 1;
	}
    }


  if (lnmb > 0)
    {
      lrad /= (double)lnmb;
      d2 /= (double)lnmb;
    }

  int sgn1 = (out1) ? -1 : 1;
  int sgn2 = out2 ? -1 : 1;
  if (state == 2)
    {
      sgn2 *= -1;
    }
  double phi2 = 0.5*(M_PI - beta);
  double hh = (state == 1) ? lrad : sgn2*(lrad + d2)*sin(phi2);
  Point pos = centr - sgn1*hh*normal;
  double Rrad = rotational->radius(0.0, rd);
  double sd = (state == 1) ? lrad/sin(phi) : (lrad + d2)*cos(phi2);
  Point Cx = rotational->direction2();
  torus = shared_ptr<Torus>(new Torus(Rrad+sgn2*sd, lrad, pos, normal, Cx));
  if (sgn2 > 0)
    {
      RectDomain dom = torus->getParameterBounds();
      try {
	torus->setParameterBounds(dom.umin(), dom.vmin()-M_PI,
				  dom.umax(), dom.vmax()-M_PI);
      }
      catch (...)
	{
	}
    }
  
#ifdef DEBUG_BLEND
  std::ofstream of2("tor_blend.g2");
  torus->writeStandardHeader(of2);
  torus->write(of2);
  RectDomain dom2 = torus->getParameterBounds();
  vector<shared_ptr<ParamCurve> > tmp_cvs = torus->constParamCurves(dom2.vmin(), true);
  tmp_cvs[0]->writeStandardHeader(of2);
  tmp_cvs[0]->write(of2);
  double tf = (sgn2 == 1) ? 0.5 : 1.5;
  vector<shared_ptr<ParamCurve> > tmp_cvs2 =
    torus->constParamCurves(dom2.vmin()+tf*M_PI, true);
  tmp_cvs2[0]->writeStandardHeader(of2);
  tmp_cvs2[0]->write(of2);
#endif
  return torus;
}

//===========================================================================
double
RevEng::computeCylinderRadius(vector<vector<RevEngPoint*> > blend_pts,
			    double width, const Point& pos, const Point& axis,
			    const Point& dir1, const Point& dir2)
//===========================================================================
{
  double eps = 1.0e-6;
  double upar2, vpar2, dist2;
  Point close, surfnorm, close2;
  double alpha = dir1.angle(dir2);
  Point linedir = width*dir1 + width*dir2;
  Point planenorm = linedir.cross(axis);
  planenorm.normalize();
  shared_ptr<Plane> plane(new Plane(pos, planenorm, linedir));
  double lrad = 0.0;
  int lnmb = 0;
  double fac = 1.0/sin(0.5*alpha);
  double div = fac - 1.0;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  plane->closestPoint(ptpos, upar2, vpar2, close2, dist2, eps);
	  if (dist2 <= approx_tol_)
	    {
	      Point pos2 = plane->ParamSurface::point(0.0, vpar2);
	      double dd2 = pos2.dist(close2);
	      double xlen = dd2/div;
	      lrad += xlen;
	      lnmb++;
	    }
	}
    }
  if (lnmb > 0)
    lrad /= (double)lnmb;
  return lrad;
}


//===========================================================================
shared_ptr<Cylinder>
RevEng::createCylinderBlend(vector<vector<RevEngPoint*> > blend_pts,
			    double rad1, const Point& pos, const Point& axis,
			    const Point& dir1, const Point& dir2, int sgn)
//===========================================================================
{
  double eps = 1.0e-6;
  double upar2, vpar2, dist2;
  Point close, surfnorm, close2;
  double alpha = dir1.angle(dir2);
  Point linedir = rad1*dir1 + rad1*dir2;
  Point planenorm = linedir.cross(axis);
  planenorm.normalize();
  shared_ptr<Plane> plane(new Plane(pos, planenorm, linedir));
  double lrad = 0.0;
  int lnmb = 0;
  double fac = 1/sin(0.5*alpha);
  double div = fac - 1.0;
  for (size_t kj=0; kj<blend_pts.size(); ++kj)
    {
      for (size_t ki=0; ki<blend_pts[kj].size(); ++ki)
	{
	  Vector3D xyz = blend_pts[kj][ki]->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);

	  plane->closestPoint(ptpos, upar2, vpar2, close2, dist2, eps);
	  if (dist2 <= approx_tol_)
	    {
	      Point pos2 = plane->ParamSurface::point(0.0, vpar2);
	      double dd2 = pos2.dist(close2);
	      double xlen = dd2/div;
	      lrad += xlen;
	      lnmb++;
	    }
	}
    }
  if (lnmb > 0)
    lrad /= (double)lnmb;
  
  Point Cx = sgn*(dir1 + dir2);
  Point Cy = axis.cross(Cx);
  // Point pos2;
  // double rad2;
  Point low = bbox_.low();
  Point high = bbox_.high();

  Point vec = dir1 + dir2;
  vec.normalize();
  Point pos3 = pos + sgn*fac*lrad*vec;
  shared_ptr<Cylinder> cyl(new Cylinder(lrad, pos3, axis, Cx));
  double diag = low.dist(high);
  cyl->setParamBoundsV(-0.5*diag,0.5*diag);
  
#ifdef DEBUG_BLEND
  std::ofstream of("blend_cyl2.g2");
  cyl->writeStandardHeader(of);
  cyl->write(of);
#endif

  
  return cyl;
}

//===========================================================================
void RevEng::buildSurfaces()
//===========================================================================
{
//   double frac = 0.75;   Fraction of points with a certain property
//   double angfac = 10.0;
//   double angtol = 5.0*anglim_;

//   // Update regions size limitation. Sort regions according to number of points
//   std::sort(regions_.begin(), regions_.end(), sort_region);

//   min_point_region_ = setSmallRegionNumber();
//   std::cout << "Min point region: " << min_point_region_ << std::endl;
  
//   First pass. Recognize elementary surfaces
//   int min_point_in = 50; 10; //20;
//   recognizeSurfaces(min_point_in, true);

//   Segmentation of composed regions
//   std::cout << "Segment composed regions" << std::endl;
//   size_t reg_size = regions_.size();
//   for (size_t ki=0; ki<reg_size; ++ki)
//     {
//       if (!regions_[ki]->hasSurface())
// 	{
// #ifdef DEBUG_DIV
// 	  std::ofstream ofs("region_to_segm.g2");
// 	  regions_[ki]->writeRegionInfo(ofs);
// #endif
	  
// 	  bool segmented = false;
// 	  if (regions_[ki]->hasBaseSf())
// 	    {
// 	      Grow sub regions according to surface type
// 	      shared_ptr<ParamSurface> primary = regions_[ki]->getBase();
// 	      if (primary->instanceType() == Class_Plane)
// 		{
// 		  std::cout << "Grow planes" << std::endl;
// 		  segmented = segmentByPlaneGrow((int)ki, min_point_in, angtol);
// 		}
// 	    }
// 	  if (!segmented)
// 	    {
// 	      segmented = segmentByContext((int)ki, min_point_in, angtol, true);
// 	    }
// 	}

//     }

// #ifdef DEBUG_DIV
//   std::cout << "Merge adjacent regions" << std::endl;
// #endif
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       if (regions_[ki]->hasSurface())
// 	{
// 	  vector<RevEngRegion*> grown_regions;
// 	  vector<HedgeSurface*> adj_surfs;
// 	  regions_[ki]->mergeAdjacentSimilar(approx_tol_, angtol,
// 					     grown_regions, adj_surfs);
// 	if (grown_regions.size() > 0 || adj_surfs.size() > 0)
// 	  updateRegionsAndSurfaces(ki, grown_regions, adj_surfs);

// 	}
//     }
// #ifdef DEBUG
//   checkConsistence("Regions3_01");

//    if (regions_.size() > 0)
//     {
//       std::cout << "Regions3_01" << std::endl;
//       std::ofstream of("regions3_01.g2");
//       std::ofstream ofs("small_regions3_01.g2");
//       writeRegionStage(of, ofs);
//      }
// #endif

  
//   Update axis
//   Point mainaxis[3];
//   bool only_surf = true;
//   std::sort(regions_.begin(), regions_.end(), sort_region);
//   int min_num = regions_[regions_.size()/8]->getNumInside();Points();
//   defineAxis(mainaxis, only_surf, min_num);

//   Just to test
//   adaptToMainAxis(mainaxis);
  
//   Second pass. Recognize also free form surfaces and utilize context information
//   recognizeSurfaces(min_point_in, false);
  
 

//   std::cout << "Number of surfaces: " << surfaces_.size() << ", give number: " << std::endl;
//   int sfix;
//   std::cin >> sfix;
//   while (sfix >=0 && sfix < (int)surfaces_.size())
//     {
//       std::ofstream ofpar("sfparpoints.txt");
//       vector<RevEngRegion*> regs = surfaces_[sfix]->getRegions();
//       for (size_t ki=0; ki<regs.size(); ++ki)
//   	{
//   	  int numpt = regs[ki]->numPoints();
//   	  for (int ka=0; ka<numpt; ++ka)
//   	    {
//   	      RevEngPoint *pt = regs[ki]->getPoint(ka);
//   	      ofpar << pt->getPar() << " " << pt->getPoint() << std::endl;
//   	    }
//   	}
//       std::cout << "New index: " << std::endl;
//       std::cin >> sfix;
//     }
  
//   bool doGrow = true; false;
//   if (doGrow)
//     {
//       std::cout << "Number of regions, pre grow with surf: " << regions_.size() << std::endl;
//       std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

//   Update adjacency between regions
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->clearRegionAdjacency();
//     }
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->setRegionAdjacency();
//     }

//   TESTING
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       if (regions_[ki]->hasSurface())
// 	{
// 	  regions_[ki]->adjustBoundaries(mean_edge_len_, approx_tol_, 10.0*anglim_);
// 	}
//     }

//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       std::cout << "ki=" << ki << ", nmb surf: " << surfaces_.size() << std::endl;
//       if (regions_[ki]->hasSurface())
// 	growSurface(ki);
//       for (int kh=0; kh<(int)surfaces_.size(); ++kh)
//       	{
//       	  int numreg = surfaces_[kh]->numRegions();
//       	  for (int ka=0; ka<numreg; ++ka)
//       	    {
//       	      RevEngRegion *reg = surfaces_[kh]->getRegion(ka);
//       	      size_t kr;
//       	      for (kr=0; kr<regions_.size(); ++kr)
//       		if (reg == regions_[kr].get())
//       		  break;
//       	      if (kr == regions_.size())
//       		std::cout << "Region4, surface 1. Obsolete region pointer, ki=" << ki << ", kh=" << kh << ". Region: " << reg << ", surface: " << surfaces_[kh].get() << std::endl;
//       	    }
//       	}
//     }
// #ifdef DEBUG_DIV
//       std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;
//       std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

//   for (int ki=0; ki<(int)regions_.size(); ++ki)
//     {
//       vector<RevEngRegion*> adjacent;
//       regions_[ki]->getAdjacentRegions(adjacent);
//       for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	{
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (adjacent[kj] == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Regions4. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	}
//     }
//   for (int ki=0; ki<(int)surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (reg == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Region4, surface 1. Obsolete region pointer, ki=" << ki << ", ka=" << ka << ". Region: " << reg << std::endl;
// 	  vector<RevEngRegion*> adjacent;
// 	  reg->getAdjacentRegions(adjacent);
// 	  for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	    {
// 	      size_t kr;
// 	      for (kr=0; kr<regions_.size(); ++kr)
// 		if (adjacent[kj] == regions_[kr].get())
// 		  break;
// 	      if (kr == regions_.size())
// 		std::cout << "Region4, surface. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	    }
// 	}
//     }

//   std::cout << "Regions4" << std::endl;
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  if (reg->hasSurface() == false)
// 	    std::cout << "Missing link surface-regions. ki=" << ki << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 	  else if (reg->getSurface(0) != surfaces_[ki].get())
// 	    std::cout << "Inconsistent link surface-regions. ki=" << ki << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 	}
//     }

//    if (regions_.size() > 0)
//     {
//       std::cout << "Regions4" << std::endl;
//       std::ofstream of("regions4.g2");
//       std::ofstream ofs("small_regions4.g2");
//       writeRegionStage(of, ofs);
//      }
// #endif
  
//   std::cout << "Number of surfaces: " << surfaces_.size() << ", give number: " << std::endl;
//   std::cin >> sfix;
//   while (sfix >=0 && sfix < (int)surfaces_.size())
//     {
//       std::ofstream ofpar("sfparpoints.txt");
//       vector<RevEngRegion*> regs = surfaces_[sfix]->getRegions();
//       for (size_t ki=0; ki<regs.size(); ++ki)
//   	{
//   	  int numpt = regs[ki]->numPoints();
//   	  for (int ka=0; ka<numpt; ++ka)
//   	    {
//   	      RevEngPoint *pt = regs[ki]->getPoint(ka);
//   	      ofpar << pt->getPar() << " " << pt->getPoint() << std::endl;
//   	    }
//   	}
//       std::cout << "New index: " << std::endl;
//       std::cin >> sfix;
//     }
  
//   Update adjacency between regions
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->clearRegionAdjacency();
//     }
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->setRegionAdjacency();
//     }

//   bool adjust_sf = false;
//   if (adjust_sf)
//     {
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       if (regions_[ki]->hasSurface())
// 	regions_[ki]->adjustWithSurf(mainaxis, min_point_region_,
// 				     approx_tol_, 10.0*anglim_);
//     }
  
// #ifdef DEBUG_DIV
//   std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;
//       std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;

//   if (regions_.size() > 0)
//     {
//       std::cout << "Regions5" << std::endl;
//       std::ofstream of("regions5.g2");
//       std::ofstream ofs("small_regions5.g2");
//       writeRegionStage(of, ofs);
//      }
// #endif
//     }
  
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       if (regions_[ki]->hasSurface())
// 	{
// 	  vector<RevEngRegion*> grown_regions;
// 	  vector<HedgeSurface*> adj_surfs;
// 	  vector<vector<RevEngPoint*> > out_groups;
// 	  regions_[ki]->adjustWithCylinder(mainaxis_, approx_tol_, anglim_,
// 					   min_point_region_, out_groups,
// 					   grown_regions, adj_surfs);
// 	for (size_t kr=0; kr<out_groups.size(); ++kr)
// 	  {
// 	    shared_ptr<RevEngRegion> reg(new RevEngRegion(regions_[ki]->getClassificationType(),
// 							  edge_class_type_,
// 							  out_groups[kr]));

// 	    reg->setPreviousReg(regions_[ki].get());
// 	    reg->setRegionAdjacency();
// 	    bool integrate = reg->integrateInAdjacent(mean_edge_len_,
// 						      min_next_, max_next_,
// 						      approx_tol_, 0.5,
// 						      max_nmb_outlier_,
// 						      regions_[ki].get());
// 	    if (!integrate)
// 	      regions_.push_back(reg);
// 	  }	  
// 	  if (grown_regions.size() > 0)
// 	    {
// 	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<regions_.size(); )
// 		    {
// 		      if (kj == ki)
// 			{
// 			  ++kj;
// 			  continue;
// 			}

// 		      if (grown_regions[kr] == regions_[kj].get())
// 			{
// 			  regions_.erase(regions_.begin()+kj);
// 			  if (kj < ki)
// 			    --ki;
// 			}
// 		      else
// 			++kj;
// 		    }
// 		}
// 	    }
// 	  for (size_t kr=0; kr<adj_surfs.size(); ++kr)
// 	    {
// 	      size_t kj;
// 	      for (kj=0; kj<surfaces_.size(); ++kj)
// 		if (surfaces_[kj].get() == adj_surfs[kr])
// 		  break;
// 	      if (kj < surfaces_.size())
// 		surfaces_.erase(surfaces_.begin()+kj);
// 	    }
// #ifdef DEBUG_DIV
// 	  for (size_t kh=0; kh<surfaces_.size(); ++kh)
// 	    {
// 	      int numreg = surfaces_[kh]->numRegions();
// 	      for (int ka=0; ka<numreg; ++ka)
// 		{
// 		  RevEngRegion *reg = surfaces_[kh]->getRegion(ka);
// 		  if (reg->hasSurface() == false)
// 		    std::cout << "Missing link surface-regions. ki=" << ki << ", kh=" << kh << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 		}
// #endif
// 	    }
// 	}
//     }
// #ifdef DEBUG_DIV
//   for (int ki=0; ki<(int)regions_.size(); ++ki)
//     {
//       vector<RevEngRegion*> adjacent;
//       regions_[ki]->getAdjacentRegions(adjacent);
//       for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	{
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (adjacent[kj] == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Regions5_2. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	}
//     }
//   for (int ki=0; ki<(int)surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (reg == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Region5_2, surface 1. Obsolete region pointer, ki=" << ki << ", ka=" << ka << ". Region: " << reg << std::endl;
// 	  vector<RevEngRegion*> adjacent;
// 	  reg->getAdjacentRegions(adjacent);
// 	  for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	    {
// 	      size_t kr;
// 	      for (kr=0; kr<regions_.size(); ++kr)
// 		if (adjacent[kj] == regions_[kr].get())
// 		  break;
// 	      if (kr == regions_.size())
// 		std::cout << "Region5_2, surface. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	    }
// 	}
//     }

//   std::cout << "Regions5_2" << std::endl;
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  if (reg->hasSurface() == false)
// 	    std::cout << "Missing link surface-regions. ki=" << ki << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 	  else if (reg->getSurface(0) != surfaces_[ki].get())
// 	    std::cout << "Inconsistent link surface-regions. ki=" << ki << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 	}
//     }

//    if (regions_.size() > 0)
//     {
//       std::cout << "Regions5_2" << std::endl;
//       std::ofstream of("regions5_2.g2");
//       std::ofstream ofs("small_regions5_2.g2");
//       writeRegionStage(of, ofs);
//       }
// #endif
   
//   Update adjacency between regions
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->clearRegionAdjacency();
//     }
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->setRegionAdjacency();
//     }

// #ifdef DEBUG_DIV
//   std::cout << "Number of regions, pre grow with surf: " << regions_.size() << std::endl;
//       std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;
// #endif
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       if (regions_[ki]->hasSurface())
// 	{
// 	  growSurface(ki);
// #ifdef DEBUG_CHECK
//   if (!regions_[ki]->isConnected())
//     std::cout << "Disconnected region (grow), ki= " << ki << " " << regions_[ki].get() << std::endl;
// #endif
//     }
// #ifdef DEBUG_DIV
//       std::cout << "Number of regions, post grow with surf: " << regions_.size() << std::endl;
//       std::cout << "Number of surfaces: " << surfaces_.size() << std::endl;
      
//   for (int ki=0; ki<(int)regions_.size(); ++ki)
//     {
//       vector<RevEngRegion*> adjacent;
//       regions_[ki]->getAdjacentRegions(adjacent);
//       for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	{
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (adjacent[kj] == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Regions6. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	}
//     }
//   for (int ki=0; ki<(int)surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (reg == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Region6, surface 1. Obsolete region pointer, ki=" << ki << ", ka=" << ka << ". Region: " << reg << std::endl;
// 	  vector<RevEngRegion*> adjacent;
// 	  reg->getAdjacentRegions(adjacent);
// 	  for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	    {
// 	      size_t kr;
// 	      for (kr=0; kr<regions_.size(); ++kr)
// 		if (adjacent[kj] == regions_[kr].get())
// 		  break;
// 	      if (kr == regions_.size())
// 		std::cout << "Region6, surface. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	    }
// 	}
//     }
// #endif

// #ifdef DEBUG_DIV
//   std::cout << "Regions6" << std::endl;
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  if (reg->hasSurface() == false)
// 	    std::cout << "Missing link surface-regions. ki=" << ki << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 	  else if (reg->getSurface(0) != surfaces_[ki].get())
// 	    std::cout << "Inconsistent link surface-regions. ki=" << ki << ", surf: " << surfaces_[ki].get() << ", region: " << reg << std::endl;
// 	}
//     }

//    if (regions_.size() > 0)
//     {
//       std::cout << "Regions6" << std::endl;
//       std::ofstream of("regions6.g2");
//       std::ofstream ofs("small_regions6.g2");
//       writeRegionStage(of, ofs);
//      }
//     }
// #endif
  
// #ifdef DEBUG_DIV
//   for (int ki=0; ki<(int)regions_.size(); ++ki)
//     {
//       vector<RevEngRegion*> adjacent;
//       regions_[ki]->getAdjacentRegions(adjacent);
//       for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	{
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (adjacent[kj] == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Pre merge. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	}
//     }

//   for (int ki=0; ki<(int)surfaces_.size(); ++ki)
//     {
//       int numreg = surfaces_[ki]->numRegions();
//       for (int ka=0; ka<numreg; ++ka)
// 	{
// 	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
// 	  size_t kr;
// 	  for (kr=0; kr<regions_.size(); ++kr)
// 	    if (reg == regions_[kr].get())
// 	      break;
// 	  if (kr == regions_.size())
// 	    std::cout << "Pre merge, surface 1. Obsolete region pointer, ki=" << ki << ", ka=" << ka << ". Region: " << reg << std::endl;
// 	  vector<RevEngRegion*> adjacent;
// 	  reg->getAdjacentRegions(adjacent);
// 	  for (size_t kj=0; kj<adjacent.size(); ++kj)
// 	    {
// 	      size_t kr;
// 	      for (kr=0; kr<regions_.size(); ++kr)
// 		if (adjacent[kj] == regions_[kr].get())
// 		  break;
// 	      if (kr == regions_.size())
// 		std::cout << "Pre merge, surface. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
// 	    }
// 	}
//     }

//    Merge surface pieces representing the surface
//   std::cout << "Pre merge: " << surfaces_.size() << " surfaces" << std::endl;
//   mergeSurfaces();
//   std::cout << "Post merge: " << surfaces_.size() << " surfaces" << std::endl;
// #endif
  
//   Update adjacency between regions
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->clearRegionAdjacency();
//     }
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       regions_[ki]->setRegionAdjacency();
//     }

//   std::cout << "Pre merge2: " << surfaces_.size() << " surfaces" << std::endl;
//   mergeSplineSurfaces();
//   std::cout << "Post merge2: " << surfaces_.size() << " surfaces" << std::endl;

//   adaptToMainAxis();
//   int stop_break = 1;
    }

//===========================================================================
void RevEng::defineAxis(Point axis[3], bool onlysurf, int min_num)
//===========================================================================
{
  int nmb_in_surf[3];
  double pi4 = 0.25*M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      axis[ka] = Point(0.0, 0.0, 0.0);
      nmb_in_surf[ka] = 0;
    }
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      shared_ptr<ElementarySurface> elem;
      if (regions_[ki]->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = regions_[ki]->getSurface(0)->surface();
	  elem = dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	}
      else if (regions_[ki]->hasBaseSf() && (!onlysurf))
	{
	  shared_ptr<ParamSurface> sf = regions_[ki]->getBase();
	  elem = dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	}
      if (elem.get() && (elem->instanceType() == Class_Plane ||
			 elem->instanceType() == Class_Cylinder))
	{
	  int ka=0;
	  Point dir = elem->direction();
	  int num = regions_[ki]->getNumInside(); //regions_[ki]->numPoints();
	  if (num < min_num)
	    continue;
	  for (ka=0; ka<3; ++ka)
	    if (nmb_in_surf[ka] > 0)
	      {
		double ang = axis[ka].angle(dir);
		ang = std::min(ang, M_PI-ang);
		if (ang < pi4)
		  break;
	      }
	  if (ka == 3)
	    {
	      for (int kb=0; kb<3; ++kb)
		if (nmb_in_surf[kb] == 0)
		  {
		    axis[kb] = dir;
		    nmb_in_surf[kb] = num;
		    break;
		  }
	    }
	  else
	    {
	      double fac1 = (double)nmb_in_surf[ka]/(double)(nmb_in_surf[ka]+num);
	      double fac2 = (double)num/(double)(nmb_in_surf[ka]+num);
	      if (axis[ka]*dir < 0.0)
		dir *= -1.0;
	      axis[ka] = fac1*axis[ka] + fac2*dir;
	      nmb_in_surf[ka] += num;
	    }
	}
    }
  for (int ka=0; ka<3; ++ka)
    axis[ka].normalize_checked();

  for (int ka=0; ka<3; ++ka)
    for (int kb=ka+1; kb<3; ++kb)
      if (nmb_in_surf[kb] > nmb_in_surf[ka])
	{
	  std::swap(axis[ka], axis[kb]);
	  std::swap(nmb_in_surf[ka], nmb_in_surf[kb]);
	}

  if (nmb_in_surf[0] > 0 && nmb_in_surf[1] == 0)
    {
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = axis[0].angle(mainaxis_[ka]);
	  ang = std::min(ang, M_PI-ang);
	  if (ang > pi4)
	    {
	      axis[1] = mainaxis_[ka];
	      break;
	    }
	}
    }
  if (nmb_in_surf[0] > 0)
    {
      axis[2] = axis[1].cross(axis[0]);
      axis[1] = axis[0].cross(axis[2]);
    }
  else
    {
      for (int ka=0; ka<3; ++ka)
	axis[ka] = mainaxis_[ka];
    }

#ifdef DEBUG_DIV
  std::ofstream ofax2("base_axis.g2");
  Point mid = 0.5*(bbox_.low() + bbox_.high());
  double len = 0.5*bbox_.low().dist(bbox_.high());
  for (int ka=0; ka<3; ++ka)
    {
      ofax2 << "410 1 0 4 0 0 0 255" << std::endl;
      ofax2 << "1" << std::endl;
      ofax2 << mid << " " << mid+len*axis[ka] << std::endl;
    }
#endif
}

//===========================================================================
void RevEng::surfaceExtractOutput(int idx,
				  vector<vector<RevEngPoint*> > out_groups,
				  vector<HedgeSurface*> prev_surfs)
//===========================================================================
{
  for (size_t kr=0; kr<out_groups.size(); ++kr)
    {
      for (size_t kh=0; kh<out_groups[kr].size(); ++kh)
	out_groups[kr][kh]->unsetRegion();  
    }
  
  int classtype = regions_[idx]->getClassification();
  for (size_t kr=0; kr<out_groups.size(); ++kr)
    {
      shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
						    edge_class_type_,
						    out_groups[kr]));
      reg->setPreviousReg(regions_[idx].get());
      reg->setRegionAdjacency();
#ifdef DEBUG_CHECK
      bool connect = reg->isConnected();
      connect = reg->isConnected();
      if (!connect)
	{
	  std::cout << "surfaceExtractOutput, disconnected region " << idx << std::endl;
	}
#endif
      bool integrate = reg->integrateInAdjacent(mean_edge_len_,
						min_next_, max_next_,
						approx_tol_, 0.5,
						max_nmb_outlier_,
						regions_[idx].get());
      if (!integrate)
	regions_.push_back(reg);
    }	  
  for (size_t kr=0; kr<prev_surfs.size(); ++kr)
    {
      size_t kj;
      for (kj=0; kj<surfaces_.size(); ++kj)
	if (surfaces_[kj].get() == prev_surfs[kr])
	  break;
      if (kj < surfaces_.size())
	surfaces_.erase(surfaces_.begin()+kj);
    }
}

//===========================================================================
bool RevEng::segmentByPlaneGrow(int ix, int min_point_in, double angtol)
//===========================================================================
{
#ifdef DEBUG
  std::ofstream ofreg("segment_reg.g2");
  regions_[ix]->writeRegionInfo(ofreg);
#endif

  vector<shared_ptr<HedgeSurface> > plane_sfs;
  vector<HedgeSurface*> prev_surfs;
  vector<vector<RevEngPoint*> > out_groups;
  regions_[ix]->segmentByPlaneGrow(mainaxis_, approx_tol_, angtol, min_point_in, 
				   plane_sfs, prev_surfs, out_groups);
#ifdef DEBUG
  for (size_t ki=0; ki<plane_sfs.size(); ++ki)
    {
      plane_sfs[ki]->surface()->writeStandardHeader(ofreg);
      plane_sfs[ki]->surface()->write(ofreg);
    }
#endif
  
  if (out_groups.size() > 0 || prev_surfs.size() > 0)
    surfaceExtractOutput(ix, out_groups, prev_surfs);
  if (plane_sfs.size() > 0)
    surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());

  bool segmented = (out_groups.size() > 0);
  return segmented;
}

//===========================================================================
bool RevEng::segmentByAxis(int ix, int min_point_in)
//===========================================================================
{
  double angtol = 5.0*anglim_;
  vector<shared_ptr<HedgeSurface> > hedgesfs;
  vector<shared_ptr<RevEngRegion> > added_reg;
  vector<vector<RevEngPoint*> > separate_groups;
  vector<RevEngPoint*> single_points;
  bool segmented = regions_[ix]->extractCylByAxis(mainaxis_, min_point_in,
						  min_point_region_,
					       approx_tol_, angtol,
					       prefer_elementary_,
					       hedgesfs, added_reg,
					       separate_groups,
					       single_points);
  if (segmented && single_points.size() > 0)
    single_points_.insert(single_points_.end(), single_points.begin(),
			  single_points.end());
#ifdef DEBUG
  if (segmented)
    {
      std::ofstream of("seg_by_axis.g2");
      int num = regions_[ix]->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of <<  num << std::endl;
      for (int ka=0; ka<num; ++ka)
	of << regions_[ix]->getPoint(ka)->getPoint() << std::endl;

      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of <<  separate_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)separate_groups[ki].size(); ++ka)
	    of << separate_groups[ki][ka]->getPoint() << std::endl;
	}
    }
#endif
  
  if (separate_groups.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, separate_groups, prev_surfs);
    }
  
  if (added_reg.size() > 0)
    regions_.insert(regions_.end(), added_reg.begin(), added_reg.end());
  if (hedgesfs.size() > 0)
    surfaces_.insert(surfaces_.end(), hedgesfs.begin(), hedgesfs.end());
  
  return segmented;
}

//===========================================================================
bool RevEng::segmentByContext(int ix, int min_point_in, double angtol, bool first)
//===========================================================================
{
#ifdef DEBUG_DIV
  vector<RevEngPoint*> branchpt = regions_[ix]->extractBranchPoints();
  if (branchpt.size() > 0)
    {
      std::ofstream ofb("branch_pts_seg.g2");
      ofb << "400 1 0 4 0 0 0 255" << std::endl;
      ofb << branchpt.size() << std::endl;
      for (size_t ki=0; ki<branchpt.size(); ++ki)
	ofb << branchpt[ki]->getPoint() << std::endl;
    }
  
#endif
  
  vector<vector<RevEngPoint*> > separate_groups;
  vector<RevEngRegion*> adj_planar = regions_[ix]->fetchAdjacentPlanar();
  vector<shared_ptr<HedgeSurface> > hedgesfs;
  vector<shared_ptr<RevEngRegion> > added_reg;
  vector<HedgeSurface*> prevsfs;
  bool segmented = false;
  if (adj_planar.size() > 1)
    {
      segmented = regions_[ix]->segmentByPlaneAxis(mainaxis_, min_point_in,
						   min_point_region_,
						   approx_tol_, angtol,
						   prefer_elementary_,
						   adj_planar, hedgesfs,
						   added_reg,
						   prevsfs, separate_groups);
    }

  if (!segmented)
    {
      // Extend with cylindrical
      vector<RevEngRegion*> adj_cyl = regions_[ix]->fetchAdjacentCylindrical();
      if (adj_cyl.size() > 0)
	adj_planar.insert(adj_planar.end(), adj_cyl.begin(), adj_cyl.end());
      if (adj_planar.size() > 0)
	segmented =
	  regions_[ix]->segmentByAdjSfContext(mainaxis_, min_point_in,
					      min_point_region_,
					       approx_tol_, angtol,
					       adj_planar, separate_groups);
    }
  
  if (!segmented)
    {
      // Search for context direction
      double angtol2 = 2.0*angtol;

      Point direction = regions_[ix]->directionFromAdjacent(angtol);
      vector<vector<RevEngPoint*> > separate_groups2;
      if (direction.dimension() == 3)
	segmented = regions_[ix]->segmentByDirectionContext(min_point_in, approx_tol_,
							    direction, angtol2,
							    separate_groups2);
      if (separate_groups2.size() > 0)
	separate_groups.insert(separate_groups.end(), separate_groups2.begin(),
			       separate_groups2.end());
      if (segmented && (!regions_[ix]->hasSurface()))
	{
	  double angtol = -1.0;
	  int pass = 1;
	  bool found = recognizeOneSurface(ix, min_point_in, angtol, pass);
	  int stop_break = 1;
	}
    }
  
#ifdef DEBUG
  if (segmented)
    {
      std::ofstream of("seg_by_context.g2");
      int num = regions_[ix]->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of <<  num << std::endl;
      for (int ka=0; ka<num; ++ka)
	of << regions_[ix]->getPoint(ka)->getPoint() << std::endl;

      for (size_t ki=0; ki<separate_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of <<  separate_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)separate_groups[ki].size(); ++ka)
	    of << separate_groups[ki][ka]->getPoint() << std::endl;
	}
    }
#endif
  
  if (separate_groups.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, separate_groups, prev_surfs);
    }

  if (added_reg.size() > 0)
    regions_.insert(regions_.end(), added_reg.begin(), added_reg.end());
  if (hedgesfs.size() > 0)
    surfaces_.insert(surfaces_.end(), hedgesfs.begin(), hedgesfs.end());
  
  return segmented;
}

//===========================================================================
void RevEng::growSurface(int& ix, int pass)
//===========================================================================
{
  vector<RevEngRegion*> grown_regions;
  int min_nmb = 5*min_point_region_;  // Should be set from distribution of how many
	  // points the regions have
  double angtol = 5.0*anglim_;
  vector<HedgeSurface*> adj_surfs;
  vector<RevEngEdge*> adj_edgs;
  regions_[ix]->growWithSurf(mainaxis_, min_nmb, min_point_region_,
			     approx_tol_, angtol, grown_regions,
			     adj_surfs, adj_edgs, (pass>1));
  updateRegionsAndSurfaces(ix, grown_regions, adj_surfs);
  for (size_t ki=0; ki<adj_edgs.size(); ++ki)
    {
      size_t kj;
      for (kj=0; kj<edges_.size(); ++kj)
	if (edges_[kj].get() == adj_edgs[ki])
	  break;
      if (kj < edges_.size())
	edges_.erase(edges_.begin()+kj);
    }
}

//===========================================================================
void RevEng::growBlendSurface(int& ix)
//===========================================================================
{
  // Collect associated blend surfaces
  RevEngEdge *edge = regions_[ix]->getBlendEdge();
  if (!edge)
    return;

#ifdef DEBUG_BLEND
  std::ofstream of("blend_grow.g2");
  regions_[ix]->writeRegionPoints(of);
#endif

  vector<RevEngRegion*> next_blend;
  RevEngRegion *adj[2];
  edge->getAdjacent(adj[0], adj[1]);
  for (int ka=0; ka<2; ++ka)
    {
      if (!adj[ka])
	continue;
      vector<RevEngEdge*> rev_edgs = adj[ka]->getAllRevEdges();
      for (size_t ki=0; ki<rev_edgs.size(); ++ki)
	{
	  RevEngRegion *blendreg = rev_edgs[ki]->getBlendRegSurf();
	  if (blendreg)
	    {
	      if (blendreg == regions_[ix].get())
		continue;
	      size_t kj=0;
	      for (kj=0; kj<next_blend.size(); ++kj)
		if (next_blend[kj] == blendreg)
		  break;
	      if (kj < next_blend.size())
		continue;
	      next_blend.push_back(blendreg);
	    }
	}
    }

#ifdef DEBUG_BLEND
  for (size_t kr=0; kr<next_blend.size(); ++kr)
    next_blend[kr]->writeRegionPoints(of);
#endif
  
  vector<RevEngRegion*> grown_regions;
  double angtol = 5.0*anglim_;
  vector<HedgeSurface*> adj_surfs;
  vector<vector<RevEngPoint*> > added_regs;
  regions_[ix]->growBlendSurf(next_blend, approx_tol_, angtol,
			      grown_regions, added_regs);

  updateRegionsAndSurfaces(ix, grown_regions, adj_surfs);
  if (added_regs.size() > 0)
    {
      vector<HedgeSurface*> prev_surfs;
      surfaceExtractOutput(ix, added_regs, prev_surfs);
    }
}

//===========================================================================
void RevEng::growMasterSurface(int& ix)
//===========================================================================
{

#ifdef DEBUG_BLEND
  std::ofstream of("master_grow.g2");
  regions_[ix]->writeRegionPoints(of);
#endif

  vector<RevEngRegion*> grown_regions;
  double angtol = 5.0*anglim_;
  vector<HedgeSurface*> adj_surfs;
  int small_lim = min_point_region_/20;
  regions_[ix]->joinToCurrent(approx_tol_, angtol, small_lim,
			      grown_regions);

  updateRegionsAndSurfaces(ix, grown_regions, adj_surfs);
}


//===========================================================================
void RevEng::updateRegionsAndSurfaces(int& ix, vector<RevEngRegion*>& grown_regions,
				      vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  if (grown_regions.size() > 0)
    {
      for (size_t kr=0; kr<grown_regions.size(); ++kr)
	{
// #ifdef DEBUG_CHECK
// 	  bool connect = grown_regions[kr]->isConnected();
// 	  connect = grown_regions[kr]->isConnected();
// 	  if (!connect)
// 	    std::cout << "updateRegionsAndSurfaces, disconnected region " << ix << kr << std::endl;
// #endif
	  size_t kj;
	  for (kj=0; kj<regions_.size(); )
	    {
	      if ((int)kj == ix)
		{
		  ++kj;
		  continue;
		}

	      if (grown_regions[kr] == regions_[kj].get())
		{
		  regions_.erase(regions_.begin()+kj);
		  if ((int)kj < ix)
		    --ix;
		}
	      else
		++kj;
	    }
	}
      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
	{
	  size_t kj;
	  for (kj=0; kj<surfaces_.size(); ++kj)
	    if (surfaces_[kj].get() == adj_surfs[kr])
	      break;
	  if (kj < surfaces_.size())
	    surfaces_.erase(surfaces_.begin()+kj);
	}
    }

  for (size_t kj=0; kj<regions_.size(); ++kj)
    regions_[kj]->setVisited(false);

#ifdef DEBUG_DIV
  std::ofstream ofpts("sfpoints2.g2");
  std::ofstream ofpar("sfparpoints2.txt");
  int numpt = regions_[ix]->numPoints();
  ofpts << "400 1 0 0" << std::endl;
  ofpts << numpt << std::endl;
  for (int ka=0; ka<numpt; ++ka)
    {
      RevEngPoint *pt = regions_[ix]->getPoint(ka);
      ofpar << pt->getPar() << " " << pt->getPoint() << std::endl;
      ofpts << pt->getPoint() << std::endl;
    }
#endif
  int stop_break = 1;
  
}

//===========================================================================
void RevEng::mergeSurfaces()
//===========================================================================
{
  // Sort according to number of associated points
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int nmb1 = surfaces_[ki]->numPoints();
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  int nmb2 = surfaces_[kj]->numPoints();
	  if (nmb2 > nmb1)
	    std::swap(surfaces_[ki], surfaces_[kj]);
	}
    }
#ifdef DEBUG
  std::ofstream ofm("merged_sfs.g2");
  std::ofstream of("surfs0.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
      
      int nreg = surfaces_[ki]->numRegions();
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  int nmb = reg->numPoints();
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << nmb << std::endl;
	  for (int kb=0; kb<nmb; ++kb)
	    {
	      of << reg->getPoint(kb)->getPoint() << std::endl;
	    }
	}
    }
#endif
  
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
#ifdef DEBUG
      std::ofstream ofn("surfsn.g2");
      for (size_t kh=0; kh<surfaces_.size(); ++kh)
	{
	  shared_ptr<ParamSurface> surf = surfaces_[kh]->surface();
	  surf->writeStandardHeader(ofn);
	  surf->write(ofn);
	}
#endif
      
      // Identify possible merge candidates
      vector<size_t> cand_ix;
      vector<double> cand_score;
      cand_ix.push_back(ki);
      cand_score.push_back(0.0);
      ClassType type;
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  double score;
	  if (surfaces_[ki]->isCompatible(surfaces_[kj].get(), anglim_,
					  approx_tol_, type, score))
	    {
	      cand_ix.push_back(kj);
	      cand_score.push_back(score);
	    }
	}

     if (cand_ix.size() > 1)
	{
	  // Sort accorading to compability
	  for (size_t kr=1; kr<cand_ix.size(); ++kr)
	    for (size_t kh=kr+1; kh<cand_ix.size(); ++kh)
	      {
		if (cand_score[kh] < cand_score[kr])
		  {
		    std::swap(cand_ix[kr], cand_ix[kh]);
		    std::swap(cand_score[kr], cand_score[kh]);
		  }
	      }
	  
#ifdef DEBUG
	  std::ofstream pre("pre_merge.g2");
	  for (size_t kr=0; kr<cand_ix.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
	      surf->writeStandardHeader(pre);
	      surf->write(pre);
	    }
#endif
	  shared_ptr<HedgeSurface> merged_surf = doMerge(cand_ix, cand_score,
							 type);
	  if (merged_surf.get())
	    {
#ifdef DEBUG
	      std::ofstream post("post_merge.g2");
	      for (size_t kr=0; kr<cand_ix.size(); ++kr)
		{
		  shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
		  surf->writeStandardHeader(post);
		  surf->write(post);
		}
	      shared_ptr<ParamSurface> surfm = merged_surf->surface();
	      surfm->writeStandardHeader(post);
	      surfm->write(post);

	      surfm->writeStandardHeader(ofm);
	      surfm->write(ofm);
#endif
	      int nreg = merged_surf->numRegions();
	      for (int ka=0; ka<nreg; ++ka)
		{
		  RevEngRegion *reg =  merged_surf->getRegion(ka);
		  int nmb = reg->numPoints();
#ifdef DEBUG
		  ofm << "400 1 0 4 255 0 0 255" << std::endl;
		  ofm << nmb << std::endl;
		  for (int kb=0; kb<nmb; ++kb)
		    {
		      ofm << reg->getPoint(kb)->getPoint() << std::endl;
		    }
#endif
		}
	     

	      std::sort(cand_ix.begin(), cand_ix.end());
	      std::swap(surfaces_[cand_ix[0]], merged_surf);
	      if (cand_ix[0] > ki)
		std::swap(surfaces_[ki], surfaces_[cand_ix[0]]);
	      for (size_t kr=cand_ix.size()-1; kr>=1; --kr)
		{
		  surfaces_.erase(surfaces_.begin()+cand_ix[kr]);
		}
	    }
	}
      int stop_break = 1;
    }

  // Limit primary surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      surfaces_[ki]->limitSurf();
    }
  
#ifdef DEBUG
  std::ofstream of1("surfs1.g2");
  shared_ptr<std::ofstream> ofp(new std::ofstream("planes_2.g2"));
  shared_ptr<std::ofstream> ofc1(new std::ofstream("cylinders_2.g2"));
  shared_ptr<std::ofstream> of0(new std::ofstream("spheres_2.g2"));
  shared_ptr<std::ofstream> ofc2(new std::ofstream("cones_2.g2"));
  shared_ptr<std::ofstream> oft(new std::ofstream("tori_2.g2"));
  shared_ptr<std::ofstream> ofsp(new std::ofstream("spline_2.g2"));
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
     shared_ptr<std::ofstream> ofs;
      if (surfaces_[ki]->isPlane())
	ofs = ofp;
      else if (surfaces_[ki]->isCylinder())
	ofs = ofc1;
       else if (surfaces_[ki]->isSphere())
	ofs = of0;
     else if (surfaces_[ki]->isCone())
	ofs = ofc2;
      else if (surfaces_[ki]->isTorus())
	ofs = oft;
      else if (surfaces_[ki]->isSpline())
	ofs = ofsp;
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of1);
      surf->write(of1);
      if (ofs.get())
	{
	  surf->writeStandardHeader(*ofs);
	  surf->write(*ofs);
	}
      
      int nreg = surfaces_[ki]->numRegions();
      int nmb = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  nmb += reg->numPoints();
	}
      of1 << "400 1 0 4 255 0 0 255" << std::endl;
      of1 << nmb << std::endl;
      *ofs << "400 1 0 4 255 0 0 255" << std::endl;
      *ofs << nmb << std::endl;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  int nmb2 = reg->numPoints();;
	  for (int kb=0; kb<nmb2; ++kb)
	    {
	      of1 << reg->getPoint(kb)->getPoint() << std::endl;
	      *ofs << reg->getPoint(kb)->getPoint() << std::endl;
	    }
	}
    }
#endif

  int stop_break2 = 1;
}

//===========================================================================
void RevEng::mergeSplineSurfaces()
//===========================================================================
{
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int code1;
      int type1 = surfaces_[ki]->instanceType(code1);
      if (type1 != Class_SplineSurface)
	continue;
      
#ifdef DEBUG
      std::ofstream ofn("surfsn.g2");
      for (size_t kh=0; kh<surfaces_.size(); ++kh)
	{
	  shared_ptr<ParamSurface> surf = surfaces_[kh]->surface();
	  surf->writeStandardHeader(ofn);
	  surf->write(ofn);
	}
#endif

      vector<RevEngRegion*> regions1 = surfaces_[ki]->getRegions();
      DirectionCone cone1 = surfaces_[ki]->surface()->normalCone();
      
      // Identify possible merge candidates
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
	{
	  int code2;
	  int type2 = surfaces_[kj]->instanceType(code2);
	  if (type2 != Class_SplineSurface)
	    continue;

	  // Check adjacency
	  vector<RevEngRegion*> regions2 = surfaces_[kj]->getRegions();
	  DirectionCone cone2 = surfaces_[kj]->surface()->normalCone();
	  DirectionCone cone3 = cone2;
	  cone3.addUnionWith(cone1);
	  // if (cone3.greaterThanPi())
	  //   continue;
	  bool can_merge = false;
	  int num_edge_between = 0;
	  for (size_t kr=0; kr<regions1.size(); ++kr)
	    for (size_t kh=0; kh<regions2.size(); ++kh)
	      {
		bool adjacent = regions1[kr]->hasAdjacentRegion(regions2[kh]);
		if (adjacent)
		  {
		    can_merge = true;
		    bool edge_between =
		      regions1[kr]->hasEdgeBetween(regions2[kh]);
		    if (edge_between)
		      ++num_edge_between;
		    int stop_break = 1;
		  }
	      }
	  if (can_merge && num_edge_between == 0)
	    {
#ifdef DEBUG	      
	      std::ofstream pre("pre_merge_spline.g2");
	      shared_ptr<ParamSurface> surf1 = surfaces_[ki]->surface();
	      surf1->writeStandardHeader(pre);
	      surf1->write(pre);
	      shared_ptr<ParamSurface> surf2 = surfaces_[kj]->surface();
	      surf2->writeStandardHeader(pre);
	      surf2->write(pre);
#endif
	      
	      shared_ptr<HedgeSurface> merged =
		doMergeSpline(surfaces_[ki], cone1, surfaces_[kj], cone2);
	      if (merged)
		{
#ifdef DEBUG
		  std::ofstream post("post_merge_spline.g2");
		  shared_ptr<ParamSurface> surfm = merged->surface();
		  surfm->writeStandardHeader(post);
		  surfm->write(post);
		  int nreg = merged->numRegions();
		  for (int ka=0; ka<nreg; ++ka)
		    {
		      RevEngRegion *reg =  merged->getRegion(ka);
		      int nmb = reg->numPoints();
		      post << "400 1 0 4 255 0 0 255" << std::endl;
		      post << nmb << std::endl;
		      for (int kb=0; kb<nmb; ++kb)
			{
			  post << reg->getPoint(kb)->getPoint() << std::endl;
			}
		    }
#endif
		  std::swap(surfaces_[ki], merged);
		  surfaces_.erase(surfaces_.begin()+kj);
		  kj--;
		}
	      int stop_break2 = 1;
	    }
	}
    }

#ifdef DEBUG
  std::ofstream of1("surfs2.g2");
  std::ofstream of1n("surfs2n.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of1);
      surf->write(of1);
      surf->writeStandardHeader(of1n);
      surf->write(of1n);
      int nreg = surfaces_[ki]->numRegions();
      int nmb = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  nmb += reg->numPoints();
	}
      of1 << "400 1 0 4 255 0 0 255" << std::endl;
      of1 << nmb << std::endl;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  surfaces_[ki]->getRegion(ka);
	  int nmb2 = reg->numPoints();;
	  for (int kb=0; kb<nmb2; ++kb)
	    {
	      of1 << reg->getPoint(kb)->getPoint() << std::endl;
	    }
	}
    }
  #endif
}

// //===========================================================================
// shared_ptr<HedgeSurface> RevEng::doMerge(vector<size_t>& cand_ix,
// 					 vector<double>& cand_score,
// 					 ClassType type)
// //===========================================================================
// {
//   size_t candsize = cand_ix.size();
//   shared_ptr<HedgeSurface> merged_surf;
//   double delscore = (cand_score[candsize-1]-cand_score[0])/(double)(candsize-1);

//   std::ofstream ofp("all_merge_points.g2");

//   // Collect regions and point clouds
//   vector<RevEngRegion*> regions;
//   vector<pair<vector<RevEngPoint*>::iterator,
// 	      vector<RevEngPoint*>::iterator> > points;
//   BoundingBox bbox(3);
//   vector<int> nmbpts;
//   vector<int> nmbreg;
//   for (size_t ki=0; ki<cand_ix.size(); ++ki)
//     {
//       HedgeSurface* surf = surfaces_[cand_ix[ki]].get();
//       vector<RevEngRegion*> reg = surf->getRegions();
//       regions.insert(regions.end(), reg.begin(), reg.end());
//       nmbreg.push_back((int)reg.size());
//       for (size_t kj=0; kj<reg.size(); ++kj)
// 	{
// 	  nmbpts.push_back(reg[kj]->numPoints());
// 	  points.push_back(std::make_pair(reg[kj]->pointsBegin(),
// 					  reg[kj]->pointsEnd()));
// 	  bbox.addUnionWith(reg[kj]->boundingBox());
// 	  ofp << "400 1 0 0" << std::endl;
// 	  int npt = reg[kj]->numPoints();
// 	  ofp << npt << std::endl;
// 	  for (int ka=0; ka<npt; ++ka)
// 	    ofp << reg[kj]->getPoint(ka)->getPoint() << std::endl;
// 	}
//     }

//   for (size_t kh=0; kh<cand_ix.size(); ++kh)
//     {
//       // Extract appropriate point clouds
//       vector<pair<vector<RevEngPoint*>::iterator,
// 		  vector<RevEngPoint*>::iterator> > curr_points;
//       vector<int> curr_nmbpts;
//       size_t ki, kj, kr;
//       for (ki=0, kr=0; ki<cand_ix.size(); ++ki)
// 	{
// 	  if (ki > 0 && cand_score[ki]-cand_score[ki-1] > delscore)
// 	    break;
// 	  for (kj=0; kj<nmbreg[ki]; ++kj, ++kr)
// 	    {
// 	      curr_points.push_back(points[kr]);
// 	      curr_nmbpts.push_back(nmbpts[kr]);
// 	    }
// 	}
//       size_t appn = ki;
      
//       // Try merge current candidate set
//       // Distinguish between the different surface types
//       shared_ptr<ParamSurface> surf;
//       if (type == Class_Plane)
// 	{
// 	  surf = doMergePlanes(curr_points, bbox, curr_nmbpts);
// 	}
//       else if (type == Class_Cylinder)
// 	{
// 	  surf = doMergeCylinders(curr_points, bbox, curr_nmbpts);
// 	}
//       else if (type == Class_Torus)
// 	{
// 	  surf = doMergeTorus(curr_points, bbox, curr_nmbpts);
// 	}
//       if (!surf.get())
// 	return merged_surf;

//       std::ofstream of0("merge_surface.g2");
//       surf->writeStandardHeader(of0);
//       surf->write(of0);

//       // Check accuracy
//       double dfac = 5.0;
//       std::ofstream of1("regions_merge.g2");
//       std::ofstream of2("in_out_merge.g2");
//       vector<vector<RevEngPoint*> > all_in;
//       double all_maxd = 0.0, all_avd = 0.0;
//       int all_inside = 0;
//       for (ki=0, kr=0; ki<cand_ix.size(); ++ki)
// 	{
// 	  for (kj=0; kj<nmbreg[ki]; ++kj, ++kr)
// 	    {
// 	      regions[kr]->writeRegionInfo(of1);
	      
// 	      double maxd, avd;
// 	      int num2;
// 	      vector<RevEngPoint*> in, out;
// 	      vector<pair<double,double> > distang;
// 	      RevEngUtils::distToSurf(points[kr].first, points[kr].second,
// 				      surf, approx_tol_, maxd, avd, num2, 
// 				      in, out, distang);

// 	      of2 << "400 1 0 4 155 50 50 255" << std::endl;
// 	      of2 << in.size() << std::endl;
// 	      for (size_t kn=0; kn<in.size(); ++kn)
// 		of2 << in[kn]->getPoint() << std::endl;
// 	      of2 << "400 1 0 4 50 155 50 255" << std::endl;
// 	      of2 << out.size() << std::endl;
// 	      for (size_t kn=0; kn<out.size(); ++kn)
// 		of2 << out[kn]->getPoint() << std::endl;
	  
// 	      double maxd_init, avd_init;
// 	      int num2_init;
// 	      regions[kr]->getAccuracy(maxd_init, avd_init, num2_init);
// 	      int num = regions[kr]->numPoints();

// 	      if (num2 < num/2 || avd > approx_tol_)
// 		{
// 		  if (ki < appn)
// 		    {
// 		      cand_ix.erase(cand_ix.begin() + ki);
// 		      cand_score.erase(cand_score.begin() + ki);

// 		      for ()
// 			{
// 			  regions.erase(regions.begin()+ki);
// 			  points.erase(points.begin()+ki);
// 			  nmbpts.erase(nmbpts.begin()+ki);
// 			}
// 		    }
// 		  else
// 		    {
// 		    }
// 		}
// 	      else
// 		{
// 		  all_in.push_back(in);
// 		  all_maxd = std::max(all_maxd, maxd);
// 		  all_avd += num*avd;
// 		  all_inside += num2;
// 		  ++ki;
// 		}
// 	  int stop_break = 1;
// 	}
//       if (cand_ix.size() <= 1)
// 	break;
//       candsize = cand_ix.size();
//     }
//   if (cand_ix.size() == candsize)
//     {
//       merged_surf =
// 	shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
//       for (size_t kj=0; kj<regions.size(); ++kj)
// 	regions[kj]->setHedge(merged_surf.get());
//     }

//   return merged_surf;
// }

//===========================================================================
shared_ptr<HedgeSurface> RevEng::doMerge(vector<size_t>& cand_ix,
					 vector<double>& cand_score,
					 ClassType type)
//===========================================================================
{
  //size_t candsize = cand_ix.size();
  shared_ptr<HedgeSurface> merged_surf;
  if (cand_ix.size() <= 1)
    return merged_surf;
  double delscore = 2.0*(cand_score[1]-cand_score[0]);

  vector<size_t> select_ix;
  select_ix.push_back(0);
  for (size_t ki=1; ki<cand_ix.size(); ++ki)
    {
      if (cand_score[ki]-cand_score[ki-1] > delscore)
	break;
      select_ix.push_back(ki);
    }
  size_t init_select = select_ix.size();

  shared_ptr<ParamSurface> prev_surf, surf;
  vector<vector<RevEngPoint*> > all_in, all_in2;
  BoundingBox bbox(3);
  double avd_all = 0.0, maxd_all = 0.0;
  int num_in_all = 0, num_all = 0, num_all2 = 0;
  vector<RevEngRegion*> regions, regions2;
  double tolfac = 2.0;
  vector<vector<pair<double,double> > > distang;
  vector<vector<double> > parvals;
  for (size_t kh=0; kh<init_select; ++kh)
    {
      surf = approxMergeSet(cand_ix, select_ix, type);
      if (!surf.get())
	break;

#ifdef DEBUG
      std::ofstream of1("merge_sf.g2");
      std::ofstream of2("in_out_merge.g2");
      surf->writeStandardHeader(of1);
      surf->write(of1);
#endif
      
      // Test accuracy
      all_in.clear();
      regions.clear();
      all_in2.clear();
      regions2.clear();
      distang.clear();
      parvals.clear();
      num_in_all = num_all = num_all2 = 0;
       vector<size_t> select_ix2;
      BoundingBox bb(3);
      size_t ki, kj;
      for (ki=0; ki<cand_ix.size(); ++ki)
	{
	  HedgeSurface* hsurf = surfaces_[cand_ix[ki]].get();
	  vector<RevEngRegion*> reg = hsurf->getRegions();
	  vector<RevEngPoint*> in, out;
	  double avd_sf = 0.0, maxd_sf = 0.0;
	  int num_in_sf = 0, num_sf = hsurf->numPoints();
	  vector<vector<pair<double,double> > > curr_distang(reg.size());
	  vector<vector<double> > curr_parvals(reg.size());
	  for (kj=0; kj<reg.size(); ++kj)
	    {
 	      double maxd, avd;
	      int num_in, num2_in;
	      int num = reg[kj]->numPoints();
	      RevEngUtils::distToSurf(reg[kj]->pointsBegin(),
				      reg[kj]->pointsEnd(),
				      surf, approx_tol_, maxd, avd, num_in, num2_in,
				      in, out, curr_parvals[kj], curr_distang[kj]);
	      maxd_sf = std::max(maxd_sf, maxd);
	      avd_sf += num*avd/(double)num_sf;
	      num_in_sf += num_in;
	    }
#ifdef DEBUG
	    of2 << "400 1 0 4 155 50 50 255" << std::endl;
	    of2 << in.size() << std::endl;
	    for (size_t kn=0; kn<in.size(); ++kn)
	      of2 << in[kn]->getPoint() << std::endl;
	    of2 << "400 1 0 4 50 155 50 255" << std::endl;
	    of2 << out.size() << std::endl;
	    for (size_t kn=0; kn<out.size(); ++kn)
	      of2 << out[kn]->getPoint() << std::endl;
#endif
	    
	    if (num_in_sf > num_sf/2 && avd_sf < approx_tol_)
	      {
		all_in.push_back(in);
		select_ix2.push_back(ki);
		bb.addUnionWith(hsurf->regionsBox());
		avd_all += num_sf*avd_sf;
		maxd_all = std::max(maxd_all, maxd_sf);
		num_in_all += num_in_sf;
		num_all += num_sf;
		regions.insert(regions.end(), reg.begin(), reg.end());
		distang.insert(distang.end(), curr_distang.begin(), curr_distang.end());
		parvals.insert(parvals.end(), curr_parvals.begin(), curr_parvals.end());
	      }
	    else
	      {
		if (num_in_sf > num_sf/2 && avd_sf < tolfac*approx_tol_)
		  {
		    all_in2.push_back(in);
		    regions2.insert(regions2.end(), reg.begin(), reg.end());
		    num_all2 += num_sf;
		  }
		if (ki < init_select)
		  init_select--;
		else
		  break;
	      }
	}
      
      avd_all /= (double)num_all;
      bbox = bb;
      if (select_ix == select_ix2)
	break;
      select_ix = select_ix2;
      
      if (ki < cand_ix.size())
	break;
      prev_surf = surf;
    }
  
  for (size_t ki=0; ki<cand_ix.size(); )
    {
      vector<size_t>::iterator found =
	std::find(select_ix.begin(), select_ix.end(), ki);
      if (found == select_ix.end())
	cand_ix.erase(cand_ix.begin()+ki);
      else
	++ki;
    }

  if (cand_ix.size() <= 1 || all_in.size() + all_in2.size() <= 1)
    return merged_surf;

  
  vector<RevEngRegion*> regions3 = regions;
  if (all_in2.size() > 0)
    {
      all_in.insert(all_in.end(), all_in2.begin(), all_in2.end());
      regions3.insert(regions3.end(), regions2.begin(), regions2.end());
      num_all2 += num_all;
    }
  
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points(all_in.size());
  vector<int> nmbpts(all_in.size());
  for (size_t kj=0; kj<all_in.size(); ++kj)
    {
      points[kj] = std::make_pair(all_in[kj].begin(),
				  all_in[kj].end());
      nmbpts[kj] = (int)all_in[kj].size();
    }

  BoundingBox bbox2(3);
  for (size_t kj=0; kj<regions3.size(); ++kj)
    {
      bbox2.addUnionWith(regions3[kj]->boundingBox());
    }

  shared_ptr<ParamSurface> surf2;
  if (type == Class_Plane)
    {
      surf2 = RevEngUtils::doMergePlanes(points, bbox2, nmbpts);
    }
  else if (type == Class_Cylinder)
    {
      surf2 = RevEngUtils::doMergeCylinders(points, bbox2, nmbpts);
    }
  else if (type == Class_Torus)
    {
      surf2 = RevEngUtils::doMergeTorus(points, bbox2, nmbpts);
    }

#ifdef DEBUG
  std::ofstream of3("merge_sf2.g2");
  std::ofstream of4("in_out_merge2.g2");
  surf2->writeStandardHeader(of3);
  surf2->write(of3);
#endif
  vector<RevEngPoint*> in2, out2;
  vector<vector<pair<double,double> > > distang2(regions3.size());
  vector<vector<double> > parvals2(regions3.size());
  double avd_all2 = 0.0, maxd_all2 = 0.0;
  int num_in_all2 = 0;
  for (size_t kj=0; kj<regions3.size(); ++kj)
    {
      double maxd, avd;
      int num_in, num2_in;
      RevEngUtils::distToSurf(regions3[kj]->pointsBegin(),
			      regions3[kj]->pointsEnd(),
			      surf2, approx_tol_, maxd, avd, num_in, num2_in,
			      in2, out2, parvals2[kj], distang2[kj]);
      maxd_all2 = std::max(maxd_all2, maxd);
      avd_all2 += regions3[kj]->numPoints()*avd/(double)num_all2;
      num_in_all2 += num_in;
    }
#ifdef DEBUG
  of4 << "400 1 0 4 155 50 50 255" << std::endl;
  of4 << in2.size() << std::endl;
  for (size_t kn=0; kn<in2.size(); ++kn)
    of4 << in2[kn]->getPoint() << std::endl;
  of4 << "400 1 0 4 50 155 50 255" << std::endl;
  of4 << out2.size() << std::endl;
  for (size_t kn=0; kn<out2.size(); ++kn)
    of4 << out2[kn]->getPoint() << std::endl;
#endif
  
  if (surf2.get() && num_in_all2 >= num_in_all && avd_all2 < avd_all &&
      num_in_all2 > num_all2/2 && avd_all2 < approx_tol_)
    {
      merged_surf =
	shared_ptr<HedgeSurface>(new HedgeSurface(surf2, regions3));
      for (size_t kj=0; kj<regions3.size(); ++kj)
	{
	  regions3[kj]->setHedge(merged_surf.get());
	  int numpt = regions3[kj]->numPoints();
	  for (int ka=0; ka<numpt; ++ka)
	    {
	      RevEngPoint *pt = regions3[kj]->getPoint(ka);
	      pt->setPar(Vector2D(parvals2[kj][2*ka],parvals2[kj][2*ka+1]));
	      pt->setSurfaceDist(distang2[kj][ka].first, distang2[kj][ka].second);
	    }
	}
    }
  else if (surf.get() && num_in_all > num_all/2 && avd_all < approx_tol_)
    {
      merged_surf =
	shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
      for (size_t kj=0; kj<regions.size(); ++kj)
	{
	  regions[kj]->setHedge(merged_surf.get());
	  int numpt = regions[kj]->numPoints();
	  for (int ka=0; ka<numpt; ++ka)
	    {
	      RevEngPoint *pt = regions[kj]->getPoint(ka);
	      pt->setPar(Vector2D(parvals[kj][2*ka],parvals[kj][2*ka+1]));
	      pt->setSurfaceDist(distang[kj][ka].first, distang[kj][ka].second);
	    }
	}
    }
#ifdef DEBUG
  if (merged_surf.get())
    {
      std::ofstream ofp("parpoints.txt");
      for (size_t kj=0; kj<regions.size(); ++kj)
	{
	  int numpt = regions[kj]->numPoints();
	  for (int ka=0; ka<numpt; ++ka)
	    {
	      RevEngPoint *pt = regions[kj]->getPoint(ka);
	      ofp << pt->getPar() << " " << pt->getPoint() << std::endl;
	    }
	}
    }
#endif
  return merged_surf;
}
 


//===========================================================================
  shared_ptr<ParamSurface> RevEng::approxMergeSet(vector<size_t>& cand_ix,
						  vector<size_t>& select_ix,
						  ClassType type)
						  
//===========================================================================
  {
  // Collect point clouds
#ifdef DEBUG
  std::ofstream ofp("all_merge_points.g2");
#endif
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  BoundingBox bbox(3);
  vector<int> nmbpts;
  for (size_t ki=0; ki<select_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[cand_ix[select_ix[ki]]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  nmbpts.push_back(reg[kj]->numPoints());
	  points.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox.addUnionWith(reg[kj]->boundingBox());
#ifdef DEBUE
	  ofp << "400 1 0 0" << std::endl;
	  int npt = reg[kj]->numPoints();
	  ofp << npt << std::endl;
	  for (int ka=0; ka<npt; ++ka)
	    ofp << reg[kj]->getPoint(ka)->getPoint() << std::endl;
#endif
	}
    }
  
  shared_ptr<ParamSurface> surf;
  if (type == Class_Plane)
    {
      surf = RevEngUtils::doMergePlanes(points, bbox, nmbpts);
    }
  else if (type == Class_Cylinder)
    {
      surf = RevEngUtils::doMergeCylinders(points, bbox, nmbpts);
    }
  else if (type == Class_Torus)
    {
      surf = RevEngUtils::doMergeTorus(points, bbox, nmbpts);
    }
  return surf;
  }
  
//===========================================================================
shared_ptr<HedgeSurface> RevEng::doMergeSpline(shared_ptr<HedgeSurface>& surf1,
					       DirectionCone& cone1,
					       shared_ptr<HedgeSurface>& surf2,
					       DirectionCone& cone2)
//===========================================================================
{
  ClassType type1 = Class_Unknown, type2 = Class_Unknown;

  // Parameterize
  // Collect all points
  vector<RevEngPoint*> all_pts;
  vector<RevEngRegion*> regions;
  int dim = surf1->surface()->dimension();
  BoundingBox bbox(dim);
  int num1 = surf1->numRegions();
  for (int ka=0; ka<num1; ++ka)
    {
      RevEngRegion *reg = surf1->getRegion(ka);
      regions.push_back(reg);
      bbox.addUnionWith(reg->getBbox());
      vector<RevEngPoint*> pts = reg->getPoints();
      all_pts.insert(all_pts.end(), pts.begin(), pts.end());
      if (reg->hasBaseSf())
	{
	  ClassType ctype = reg->getBase()->instanceType();
	  if (type1 == Class_Unknown)
	    type1 = ctype;
	  else if (type1 != ctype)
	    type1 = Class_Unknown;  // Can be a problem with several regions
	}
    }
  int num2 = surf2->numRegions();
  for (int ka=0; ka<num2; ++ka)
    {
      RevEngRegion *reg = surf2->getRegion(ka);
      regions.push_back(reg);
      bbox.addUnionWith(reg->getBbox());
      vector<RevEngPoint*> pts = reg->getPoints();
      all_pts.insert(all_pts.end(), pts.begin(), pts.end());
      if (reg->hasBaseSf())
	{
	  ClassType ctype = reg->getBase()->instanceType();
	  if (type2 == Class_Unknown)
	    type2 = ctype;
	  else if (type1 != ctype)
	    type2 = Class_Unknown;
	}
    }
  
  double avd_all = 0.0, maxd_all = 0.0;
  int num_in_all = 0;
  for (size_t ki=0; ki<regions.size(); ++ki)
    {
      double avdr, maxdr;
      int num_inr, num2_inr;
      regions[ki]->getAccuracy(maxdr, avdr, num_inr, num2_inr);
      maxd_all = std::max(maxd_all, maxdr);
      double wgt = (double)regions[ki]->numPoints()/(double)all_pts.size();
      avd_all += wgt*avdr;
      num_in_all += num_inr;
    }
  
  // Define base surface
  ClassType type = Class_Unknown;
  if (type1 == type2 && type1 != Class_Unknown)
    {
      type = type1;
    }
  else if (type1 != Class_Unknown && type2 != Class_Unknown)
    {
    }
  else if (type1 != Class_Unknown)
    {
      type = type1;
    }
  else if (type2 != Class_Unknown)
    {
      type = type2;
    }

  shared_ptr<ParamSurface> projectsf;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > all_pts2(1);
  all_pts2[0] = std::make_pair(all_pts.begin(),all_pts.end());
  vector<int> nmbpts(1, (int)all_pts.size());
  if (type == Class_Plane)
    {
      projectsf = RevEngUtils::doMergePlanes(all_pts2, bbox, nmbpts, false);
    }
  else if (type == Class_Cylinder)
    {
      projectsf = RevEngUtils::doMergeCylinders(all_pts2, bbox, nmbpts, false);
    }
  else if (type == Class_Torus)
    {
      projectsf = RevEngUtils::doMergeTorus(all_pts2, bbox, nmbpts);
    }

  #ifdef DEBUG
  if (projectsf.get())
    {
      std::ofstream of1("projectsf.g2");
      projectsf->writeStandardHeader(of1);
      projectsf->write(of1);
    }
#endif
  
  // Parameterize on surface
  vector<double> data;
  vector<double> param;
  int inner_u = 0, inner_v = 0;
  bool close1 = false, close2 = false;
  if (projectsf.get())
    RevEngUtils::parameterizeOnPrimary(all_pts, projectsf, data, param,
				       inner_u, inner_v, close1, close2);
  else
    {
     double lambda[3];
     double eigenvec[3][3];
     RevEngUtils::principalAnalysis(all_pts, lambda, eigenvec);
     Point xAxis(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
     Point yAxis(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
     RevEngUtils::parameterizeWithPlane(all_pts, bbox, xAxis, yAxis,
					data, param);
    }
  
  // Approximate
  double maxd, avd;
  int num_out;
  int order = 4;
  int ncoef1 = order + inner_u;
  int ncoef2 = order + inner_v;
  int max_iter = 6;
  double belt = 0.01;
  shared_ptr<HedgeSurface> merged;
  shared_ptr<SplineSurface> surf;
  vector<double> parvals;
  try {
    surf =
    RevEngUtils::surfApprox(data, dim, param, order, order, ncoef1, ncoef2,
			    close1, close2, max_iter, approx_tol_, maxd, avd,
			    num_out, parvals, belt);
  }
  catch (...)
    {
      return merged;
    }
  int num_in = (int)all_pts.size() - num_out;

#ifdef DEBUG
  if (surf.get())
    {
      std::ofstream of2("mergedsf.g2");
      surf->writeStandardHeader(of2);
      surf->write(of2);
      for (size_t ki=0; ki<regions.size(); ++ki)
	regions[ki]->writeRegionInfo(of2);
    }
#endif

  BoundingBox bb1 = surf1->boundingBox();
  BoundingBox bb2 = surf2->boundingBox();
  bb1.addUnionWith(bb2);
  BoundingBox bb3 = surf->boundingBox();
  double diag1 = bb1.low().dist(bb1.high());
  double diag3 = bb3.low().dist(bb3.high());
  double frac1 = (double)num_in_all/(double)all_pts.size();
  double frac2 = (double)num_in/(double)all_pts.size();
  if (frac2 > 2.0*frac1/3.0 && frac2 > 0.5 && avd <= approx_tol_ &&
      diag3 < 1.5*diag1)
    {
      merged = shared_ptr<HedgeSurface>(new HedgeSurface(surf, regions));
      for (size_t kj=0, kr=0; kj<regions.size(); ++kj)
	{
	  regions[kj]->setHedge(merged.get());
	  int num = regions[kj]->numPoints();
	  for (int ka=0; ka<num; ++ka, ++kr)
	    {
	      RevEngPoint *pt = regions[kj]->getPoint(ka);
	      Vector2D par(parvals[2*kr], parvals[2*kr+1]);
	      pt->setPar(par);
	    }
	}
    }
  return merged;
}

// //===========================================================================
// void RevEng::recognizePlanes()
// //===========================================================================
// {
//   std::ofstream planeout("plane.g2");
//   double anglim = 0.1;
//   size_t nmbsfs = surfaces_.size();
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       // Check type
//       bool planar = regions_[ki]->planartype();
//       DirectionCone normalcone = regions_[ki]->getNormalCone();
//       if (!planar /*&& (normalcone.angle() > anglim ||
// 		    normalcone.greaterThanPi())*/)
// 	continue;  // Not a probable plane

//       // Try to fit the point cloud with a plane
//       vector<shared_ptr<HedgeSurface> > plane_sfs;
//       vector<HedgeSurface*> prev_surfs;
//       bool found = regions_[ki]->extractPlane(approx_tol_, min_point_region_,
// 					      10.0*anglim_, plane_sfs, prev_surfs, planeout);
//       for (size_t kr=0; kr<prev_surfs.size(); ++kr)
// 	{
// 	  size_t kj;
// 	  for (kj=0; kj<surfaces_.size(); ++kj)
// 	    if (surfaces_[kj].get() == prev_surfs[kr])
// 	      break;
// 	  if (kj < surfaces_.size())
// 	    surfaces_.erase(surfaces_.begin()+kj);
// 	}
//       if (plane_sfs.size() > 0)
// 	{
// 	  surfaces_.insert(surfaces_.end(), plane_sfs.begin(), plane_sfs.end());

// 	  vector<RevEngRegion*> grown_regions;
// 	  int min_nmb = 5*min_point_region_;  // Should be set from distribution of how many
// 	  // points the regions have
// 	  vector<HedgeSurface*> adj_surfs;
// 	  regions_[ki]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
// 	  if (grown_regions.size() > 0)
// 	    {
// 	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<regions_.size(); )
// 		    {
// 		      if (kj == ki)
// 			{
// 			  ++kj;
// 			  continue;
// 			}

// 		      if (grown_regions[kr] == regions_[kj].get())
// 			{
// 			  regions_.erase(regions_.begin()+kj);
// 			  if (kj < ki)
// 			    --ki;
// 			}
// 		      else
// 			++kj;
// 		    }
// 		}
// 	      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<surfaces_.size(); ++kj)
// 		    if (surfaces_[kj].get() == adj_surfs[kr])
// 		      break;
// 		  if (kj < surfaces_.size())
// 		    surfaces_.erase(surfaces_.begin()+kj);
// 		}
// 	    }
// 	}

//       for (size_t kj=0; kj<regions_.size(); ++kj)
// 	regions_[kj]->setVisited(false);


//       // if not planar
//       // continue;

//       // apply method for recognition of plane
//       // Should RevEngPoint be given as input or should the method take a
//       // vector of Points as input to enforce independence on a triangulation?
//       // The class HedgeSurface must lie in compositemodel due to the
//       // connection to RevEngRegion
//       // Should the computation take place in HedgeSurface? Or in RevEngRegion?
//       // Or in this class?
      
//       // some points may be disassembled

//       // check if the number of points in the plane is large enough
//       // if not, should the region be removed?
//       int stop_break = 1;
//     }
  
//   std::ofstream ofp("planes2.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(ofp);
//       surf->write(ofp);
//     }

//   if (regions_.size() > 0)
//     {
//       std::ofstream ofr("regions3.g2");
//       for (size_t kr=0; kr<regions_.size(); ++kr)
// 	{
// 	  // BoundingBox bbox = regions_[kr]->boundingBox();
// 	  // if (bbox.low().dist(bbox.high()) < 0.1)
// 	  //   std::cout << "Small bounding box" << std::endl;
// 	  if (regions_[kr]->numPoints() < 5)
// 	    continue;
// 	  ofr << "400 1 0 0" << std::endl;
// 	  int nmb = regions_[kr]->numPoints();
// 	  ofr << nmb << std::endl;
// 	  for (int ki=0; ki<nmb; ++ki)
// 	    {
// 	      ofr << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
// 	    }
// 	}
//     }

//   if (surfaces_.size() - nmbsfs > 1)
//     mergePlanes(nmbsfs, surfaces_.size());

//   std::ofstream of("surfaces0.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(of);
//       surf->write(of);
//     }
// }

// //===========================================================================
// void RevEng::recognizeCylinders()
// //===========================================================================
// {
//   std::ofstream cylout("cylinder.g2");
//   int nmbsfs = surfaces_.size();
//   for (size_t ki=0; ki<regions_.size(); ++ki)
//     {
//       // Check type
//       bool cyl = regions_[ki]->cylindertype();
//       // if (!cyl)
//       // 	continue;
      
//       // So far, try to regognize cylinders
//       vector<shared_ptr<HedgeSurface> > cyl_sfs;
//       vector<HedgeSurface*> prev_surfs;
//       bool found = regions_[ki]->extractCylinder(approx_tol_, min_point_region_,
// 						 10.0*anglim_,
// 						 mean_edge_len_, cyl_sfs, prev_surfs,
// 						 cylout);
//       for (size_t kr=0; kr<prev_surfs.size(); ++kr)
// 	{
// 	  size_t kj;
// 	  for (kj=0; kj<surfaces_.size(); ++kj)
// 	    if (surfaces_[kj].get() == prev_surfs[kr])
// 	      break;
// 	  if (kj < surfaces_.size())
// 	    {
// 	      surfaces_.erase(surfaces_.begin()+kj);
// 	      nmbsfs--;
// 	    }
// 	}
//       if (cyl_sfs.size() > 0)
// 	{
// 	  surfaces_.insert(surfaces_.end(), cyl_sfs.begin(), cyl_sfs.end());
      
// 	  vector<RevEngRegion*> grown_regions;
// 	  int min_nmb = 10*min_point_region_;  // Should be set from distribution of how many
// 	  // points the regions have
// 	  vector<HedgeSurface*> adj_surfs;
// 	  regions_[ki]->growWithSurf(min_nmb, approx_tol_, grown_regions, adj_surfs);
// 	  if (grown_regions.size() > 0)
// 	    {
// 	      for (size_t kr=0; kr<grown_regions.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<regions_.size(); )
// 		    {
// 		      if (kj == ki)
// 			{
// 			  ++kj;
// 			  continue;
// 			}

// 		      if (grown_regions[kr] == regions_[kj].get())
// 			{
// 			  regions_.erase(regions_.begin()+kj);
// 			  if (kj < ki)
// 			    ki--;
// 			}
// 		      else
// 			++kj;
// 		    }
// 		}
// 	      for (size_t kr=0; kr<adj_surfs.size(); ++kr)
// 		{
// 		  size_t kj;
// 		  for (kj=0; kj<surfaces_.size(); ++kj)
// 		    if (surfaces_[kj].get() == adj_surfs[kr])
// 		      break;
// 		  if (kj < surfaces_.size())
// 		    {
// 		      surfaces_.erase(surfaces_.begin()+kj);
// 		      if (kj < nmbsfs)
// 			nmbsfs--;
// 		    }
// 		}
// 	    }
// 	}

//       for (size_t kj=0; kj<regions_.size(); ++kj)
// 	regions_[kj]->setVisited(false);
      
//       // if not planar
//       // continue;

//       // apply method for recognition of plane
//       // Should RevEngPoint be given as input or should the method take a
//       // vector of Points as input to enforce independence on a triangulation?
//       // The class HedgeSurface must lie in compositemodel due to the
//       // connection to RevEngRegion
//       // Should the computation take place in HedgeSurface? Or in RevEngRegion?
//       // Or in this class?
      
//       // some points may be disassembled

//       // check if the number of points in the plane is large enough
//       // if not, should the region be removed?
//     }

  
//   std::ofstream ofp("cylinders2.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(ofp);
//       surf->write(ofp);
//     }

//   if (regions_.size() > 0)
//     {
//       std::ofstream ofr("regions4.g2");
//       for (size_t kr=0; kr<regions_.size(); ++kr)
// 	{
// 	  // BoundingBox bbox = regions_[kr]->boundingBox();
// 	  // if (bbox.low().dist(bbox.high()) < 0.1)
// 	  //   std::cout << "Small bounding box" << std::endl;
// 	  if (regions_[kr]->numPoints() < 5)
// 	    continue;
// 	  ofr << "400 1 0 0" << std::endl;
// 	  int nmb = regions_[kr]->numPoints();
// 	  ofr << nmb << std::endl;
// 	  for (int ki=0; ki<nmb; ++ki)
// 	    {
// 	      ofr << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
// 	    }
// 	}
//     }

//   // Now the surfaces between nmbsfs and surfaces_.size() are cylinder surfaces
//   // and can be merged if they represent the same cylinder
//   if (surfaces_.size() + nmbsfs > 1)
//     mergeCylinders(nmbsfs, surfaces_.size());
//   //    mergeCylinders(0, surfaces_.size());

//   std::ofstream of("surfaces.g2");
//   for (size_t ki=0; ki<surfaces_.size(); ++ki)
//     {
//       shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
//       surf->writeStandardHeader(of);
//       surf->write(of);
//     }
// }

//===========================================================================
void RevEng::mergePlanes(size_t first, size_t last)
//===========================================================================
{
  // Sort according to number of associated points
  for (size_t ki=first; ki<last; ++ki)
    {
      int nmb1 = surfaces_[ki]->numPoints();
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  int nmb2 = surfaces_[kj]->numPoints();
	  if (nmb2 > nmb1)
	    std::swap(surfaces_[ki], surfaces_[kj]);
	}
    }

#ifdef DEBUG_DIV
  std::ofstream of("planes0.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
    }
#endif
  
  // Find similar planes
  for (size_t ki=first; ki<last; ++ki)
    {
      vector<size_t> cand_ix;
      int code;
      ClassType type1 = surfaces_[ki]->instanceType(code);
      if (type1 != Class_Plane && type1 != Class_BoundedSurface)
	continue;

      shared_ptr<ParamSurface> surf1 = surfaces_[ki]->surface();
      shared_ptr<Plane> pla1 =
	dynamic_pointer_cast<Plane,ParamSurface>(surf1);
      if (!pla1.get())
	{
	  shared_ptr<BoundedSurface> bdsf1 =
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
	  if (bdsf1.get())
	    {
	      surf1 = bdsf1->underlyingSurface();
	      pla1 = dynamic_pointer_cast<Plane,ParamSurface>(surf1);
	    }
	}
      if (!pla1.get())
	continue; 

      cand_ix.push_back(ki); 
      Point norm1 = pla1->getNormal();
      Point pos1 = pla1->getPoint();
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  ClassType type2 = surfaces_[kj]->instanceType(code);
	  if (type2 != Class_Plane && type2 != Class_BoundedSurface)
	    continue;
	  
	  shared_ptr<ParamSurface> surf2 = surfaces_[kj]->surface();
	  shared_ptr<Plane> pla2 =
	    dynamic_pointer_cast<Plane,ParamSurface>(surf2);
	  if (!pla2.get())
	    {
	      shared_ptr<BoundedSurface> bdsf2 =
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	      if (bdsf2.get())
		{
		  surf2 = bdsf2->underlyingSurface();
		  pla2 = dynamic_pointer_cast<Plane,ParamSurface>(surf2);
		}
	    }
 	  if (!pla2.get())
	    continue;  // Should not happen
	  Point norm2 = pla2->getNormal();
	  Point pos2 = pla2->getPoint();

	  double ang = norm1.angle(norm2);
	  ang = std::min(ang, M_PI-ang);
	  Point pos2_0 = pos2 - ((pos2-pos1)*norm1)*norm1;
	  Point pos1_0 = pos1 - ((pos2-pos1)*norm2)*norm2;
	  double pdist1 = pos2.dist(pos2_0);
	  double pdist2 = pos1.dist(pos1_0);

	  if (ang > 2.0*anglim_ || pdist1 > 5.0*approx_tol_ || pdist2 > 5.0*approx_tol_)
	    continue;
	  cand_ix.push_back(kj);
	}

      if (cand_ix.size() > 1)
	{
#ifdef DEBUG_DIV
	  std::ofstream pre("pre_merge.g2");
	  for (size_t kr=0; kr<cand_ix.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
	      surf->writeStandardHeader(pre);
	      surf->write(pre);
	    }
#endif
	  shared_ptr<HedgeSurface> merged_surf;// = doMergePlanes(cand_ix);
	  if (merged_surf.get())
	    {
#ifdef DEBUG_DIV
	      std::ofstream post("post_merge.g2");
	      for (size_t kr=0; kr<cand_ix.size(); ++kr)
		{
		  shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
		  surf->writeStandardHeader(post);
		  surf->write(post);
		}
	      shared_ptr<ParamSurface> surfm = merged_surf->surface();
	      surfm->writeStandardHeader(post);
	      surfm->write(post);
#endif  
	      std::swap(surfaces_[cand_ix[0]], merged_surf);
	      if (cand_ix[0] > ki)
		std::swap(surfaces_[ki], surfaces_[cand_ix[0]]);
	      for (size_t kr=cand_ix.size()-1; kr>=1; --kr)
		{
		  surfaces_.erase(surfaces_.begin()+cand_ix[kr]);
		}
	      last -= (cand_ix.size()-1);
	    }
	}
      int stop_break = 1;
    }

}



//===========================================================================
void RevEng::mergeCylinders(size_t first, size_t last)
//===========================================================================
{
  // Sort according to number of associated points
  for (size_t ki=first; ki<last; ++ki)
    {
      int nmb1 = surfaces_[ki]->numPoints();
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  int nmb2 = surfaces_[kj]->numPoints();
	  if (nmb2 > nmb1)
	    std::swap(surfaces_[ki], surfaces_[kj]);
	}
    }

#ifdef DEBUG_DIV
  std::ofstream cyl("sorted_cylinders.g2");
  for (size_t ki=first; ki<last; ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      surf->writeStandardHeader(cyl);
      surf->write(cyl);
    }
#endif
  
  // Find similar cylinders
  for (size_t ki=first; ki<last; ++ki)
    {
      vector<size_t> cand_ix;
      int code;
      ClassType type1 = surfaces_[ki]->instanceType(code);
      if (type1 != Class_Cylinder && type1 != Class_BoundedSurface)
	continue;

      shared_ptr<ParamSurface> surf1 = surfaces_[ki]->surface();
      shared_ptr<Cylinder> cyl1 =
	dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
      if (!cyl1.get())
	{
	  shared_ptr<BoundedSurface> bdsf1 =
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
	  if (bdsf1.get())
	    {
	      surf1 = bdsf1->underlyingSurface();
	      cyl1 = dynamic_pointer_cast<Cylinder,ParamSurface>(surf1);
	    }
	}
      if (!cyl1.get())
	continue; 

      cand_ix.push_back(ki); 
      Point axis1 = cyl1->getAxis();
      Point pos1 = cyl1->getLocation();
      double rad1 = cyl1->getRadius();
      double dlim = std::max(0.01*rad1, approx_tol_);
      for (size_t kj=ki+1; kj<last; ++kj)
	{
	  ClassType type2 = surfaces_[kj]->instanceType(code);
	  if (type2 != Class_Cylinder && type2 != Class_BoundedSurface)
	    continue;
	  
	  shared_ptr<ParamSurface> surf2 = surfaces_[kj]->surface();
	  shared_ptr<Cylinder> cyl2 =
	    dynamic_pointer_cast<Cylinder,ParamSurface>(surf2);
	  if (!cyl2.get())
	    {
	      shared_ptr<BoundedSurface> bdsf2 =
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	      if (bdsf2.get())
		{
		  surf2 = bdsf2->underlyingSurface();
		  cyl2 = dynamic_pointer_cast<Cylinder,ParamSurface>(surf2);
		}
	    }
 	  if (!cyl2.get())
	    continue;  // Should not happen
	  Point axis2 = cyl2->getAxis();
	  Point pos2 = cyl2->getLocation();
	  double rad2 = cyl2->getRadius();

	  double ang = axis1.angle(axis2);
	  ang = std::min(ang, M_PI-ang);
	  Point pos2_0 = pos1 + ((pos2-pos1)*axis1)*axis1;
	  double pdist = pos2.dist(pos2_0);

	  if (ang > anglim_ || pdist > 5.0*dlim ||
	      fabs(rad2-rad1) > dlim)
	    continue;
	  cand_ix.push_back(kj);
	}

      if (cand_ix.size() > 1)
	{
#ifdef DEBUG_DIV
	  std::ofstream pre("pre_merge.g2");
	  for (size_t kr=0; kr<cand_ix.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
	      surf->writeStandardHeader(pre);
	      surf->write(pre);
	    }
#endif
	  shared_ptr<HedgeSurface> merged_surf;// = doMergeCylinders(cand_ix);
	  if (merged_surf.get())
	    {
#ifdef DEBUG_DIV
	      std::ofstream post("post_merge.g2");
	      for (size_t kr=0; kr<cand_ix.size(); ++kr)
		{
		  shared_ptr<ParamSurface> surf = surfaces_[cand_ix[kr]]->surface();
		  surf->writeStandardHeader(post);
		  surf->write(post);
		}
	      shared_ptr<ParamSurface> surfm = merged_surf->surface();
	      surfm->writeStandardHeader(post);
	      surfm->write(post);
#endif
	      
	      std::swap(surfaces_[cand_ix[0]], merged_surf);
	      if (cand_ix[0] > ki)
		std::swap(surfaces_[ki], surfaces_[cand_ix[0]]);
	      for (size_t kr=cand_ix.size()-1; kr>=1; --kr)
		{
		  surfaces_.erase(surfaces_.begin()+cand_ix[kr]);
		}
	      last -= (cand_ix.size()-1);
	    }
	}
      int stop_break = 1;
    }
}


//===========================================================================
void RevEng::adaptToMainAxis(Point mainaxis[3])
//===========================================================================
{
  // For each axis
  double angtol = 0.1;
  for (int ka=0; ka<3; ++ka)
    {
      // Collect surfaces with almost the same axis
      vector<size_t> sf_ix;
      for (size_t ki=0; ki<surfaces_.size(); ++ki)
	{
	  shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
	  if (elem.get() && elem->instanceType() != Class_Sphere)
	    {
	      Point axis = elem->direction();
	      double ang = mainaxis[ka].angle(axis);
	      ang = std::min(ang, M_PI-ang);
	      if (ang <= angtol)
		sf_ix.push_back(ki);
	    }
	}

      // Check axis with cylinders
      vector<size_t> cyl_ix;
      for (size_t ki=0; ki<sf_ix.size(); ++ki)
	{
	  shared_ptr<ParamSurface> surf = surfaces_[sf_ix[ki]]->surface();
	  if (surf->instanceType() == Class_Cylinder)
	    cyl_ix.push_back(sf_ix[ki]);
	}

      if (cyl_ix.size() > 0)
	cylinderFit(cyl_ix, mainaxis, ka);
      int stop_break = 1;
    }
}

//===========================================================================
void RevEng::adaptToMainAxis()
//===========================================================================
{
  vector<SurfaceProperties> sfprop;
  collectAxis(sfprop);

  // Sort surfaces into groups with roughly the same axis
  double epsang = 0.05;
  vector<DirectionCone> axis_cone;
  vector<vector<size_t> > group_ixs;
  vector<int> num_pts;
  for (size_t ki=0; ki<sfprop.size(); ++ki)
    {
      size_t kr;
      for (kr=0; kr<axis_cone.size(); ++kr)
	{
	  Point axis = sfprop[ki].dir_;
	  Point centre = axis_cone[kr].centre();
	  if (axis*centre < 0.0)
	    axis *= -1;
	  double angle = axis_cone[kr].unionAngle(axis);
	  if (angle <= epsang)
	    {
	      axis_cone[kr].addUnionWith(axis);
	      group_ixs[kr].push_back(ki);
	      num_pts[kr] += sfprop[ki].num_points_;
	      break;
	    }
	}
      if (kr == axis_cone.size())
	{
	  DirectionCone cone(sfprop[ki].dir_);
	  axis_cone.push_back(cone);
	  vector<size_t> ixs;
	  ixs.push_back(ki);
	  group_ixs.push_back(ixs);
	  num_pts.push_back(sfprop[ki].num_points_);
	}
    }

  // Sort according to number of points assiciated to the surface
  for (size_t ki=0; ki<num_pts.size(); ++ki)
    for (size_t kj=ki+1; kj<num_pts.size(); ++kj)
      {
	if (num_pts[kj] > num_pts[ki])
	  {
	    std::swap(num_pts[ki], num_pts[kj]);
	    std::swap(axis_cone[ki], axis_cone[kj]);
	    std::swap(group_ixs[ki], group_ixs[kj]);
	  }
      }

  for (size_t ki=0; ki<num_pts.size(); ++ki)
    {
      // Identify cylinders with the "same" axis including location
      if (group_ixs[ki].size() < 2)
	continue;
      for (size_t kj=0; kj<group_ixs[ki].size(); ++kj)
	{
	  if (sfprop[group_ixs[ki][kj]].type_ != Class_Cylinder &&
	      sfprop[group_ixs[ki][kj]].prev_type_  != Class_Cylinder)
	    continue;
	  Point loc1 = sfprop[group_ixs[ki][kj]].loc_;
	  Point vec1 = sfprop[group_ixs[ki][kj]].dir_;
	  double rad1 = sfprop[group_ixs[ki][kj]].rad1_;
	  vector<int> adapt_ix;
	  adapt_ix.push_back(sfprop[group_ixs[ki][kj]].sfix_);
	  for (size_t kr=kj+1; kr<group_ixs[ki].size(); ++kr)
	    {
	      if (sfprop[group_ixs[ki][kr]].type_ != Class_Cylinder &&
		  sfprop[group_ixs[ki][kr]].prev_type_  != Class_Cylinder)
		continue;
	      Point loc2 = sfprop[group_ixs[ki][kr]].loc_;
	      Point loc2_0 = loc1 + ((loc2-loc1)*vec1)*vec1;
	      double rad2 = sfprop[group_ixs[ki][kr]].rad1_;
	      double pdist = loc2.dist(loc2_0);
	      double dlim = std::max(0.05*std::max(rad1,rad2), approx_tol_);
	      if (pdist < dlim)
		{
		  adapt_ix.push_back(sfprop[group_ixs[ki][kr]].sfix_);
		}
	    }

	  if (adapt_ix.size() > 1)
	    {
	      // Try simultanous fitting of cylinders
	      cylinderFit(adapt_ix, axis_cone[ki].centre());
	      int stop_break = 1;
	    }
	}
    }
  int stop_break = 1;
}

//===========================================================================
void RevEng::cylinderFit(vector<size_t>& sf_ix, Point mainaxis[3], int ix)
//===========================================================================
{
  // Collect data points
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  BoundingBox bbox(3);
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox.addUnionWith(reg[kj]->boundingBox());
	}
    }

  Point axis, Cx, Cy;
  RevEngUtils::computeAxis(points, axis, Cx, Cy);

  Point low = bbox.low();
  Point high = bbox.high();
  Point pos;
  double radius;
  RevEngUtils::computeCylPosRadius(points, low, high, axis, Cx, Cy, pos,
				   radius);

#ifdef DEBUG
  std::ofstream of("cylinderfit.g2");
#endif
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points0;
      BoundingBox bbox0(3);
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points0.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox0.addUnionWith(reg[kj]->boundingBox());
	}
      Point low0 = bbox0.low();
      Point high0 = bbox0.high();
      Point pos0;
      double radius0;
      RevEngUtils::computeCylPosRadius(points0, low0, high0, axis, Cx, Cy, 
				      pos0, radius0);

      Point pos1;
      double radius1;
      RevEngUtils::computeCylPosRadius(points0, low0, high0, mainaxis[ix],
				       mainaxis[(ix+1)%3], mainaxis[(ix+2)%3],
				       pos1, radius1);
      shared_ptr<Cylinder> cyl0(new Cylinder(radius0, pos0, axis, Cy));
#ifdef DEBUG
      cyl0->writeStandardHeader(of);
      cyl0->write(of);
#endif
      shared_ptr<Cylinder> cyl1(new Cylinder(radius1, pos1, mainaxis[ix],
					     mainaxis[(ix+2)%3]));
#ifdef DEBUG
      cyl1->writeStandardHeader(of);
      cyl1->write(of);
#endif
      
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  double maxd0, avd0;
	  int num_in0, num2_in0;
	  vector<RevEngPoint*> in0, out0;
	  vector<pair<double,double> > distang0;
	  vector<double> parvals0;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl0, approx_tol_, maxd0, avd0, num_in0, num2_in0,
				  in0, out0, parvals0, distang0);
#ifdef DEBUG
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in0.size() << std::endl;
	  for (size_t kr=0; kr<in0.size(); ++kr)
	    of << in0[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out0.size() << std::endl;
	  for (size_t kr=0; kr<out0.size(); ++kr)
	    of << out0[kr]->getPoint() << std::endl;
#endif
	  
	  double maxd1, avd1;
	  int num_in1, num2_in1;
	  vector<RevEngPoint*> in1, out1;
	  vector<pair<double,double> > distang1;
	  vector<double> parvals1;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl1, approx_tol_, maxd1, avd1, num_in1, num2_in1,
				  in1, out1, parvals1, distang1);
#ifdef DEBUG
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in1.size() << std::endl;
	  for (size_t kr=0; kr<in1.size(); ++kr)
	    of << in1[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out1.size() << std::endl;
	  for (size_t kr=0; kr<out1.size(); ++kr)
	    of << out1[kr]->getPoint() << std::endl;
#endif
	  int stop_break = 1;
	}
    }
}

//===========================================================================
void RevEng::cylinderFit(vector<int>& sf_ix, Point normal)
//===========================================================================
{
  // Collect data points
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  BoundingBox bbox(3);
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox.addUnionWith(reg[kj]->boundingBox());
	}
    }

  Point axis, Cx, Cy;
  RevEngUtils::computeAxis(points, axis, Cx, Cy);

  Point low = bbox.low();
  Point high = bbox.high();
  Point pos;
  double radius;
  RevEngUtils::computeCylPosRadius(points, low, high, axis, Cx, Cy, pos,
				   radius);

#ifdef DEBUG
  std::ofstream of("cylinderfit.g2");
#endif
  for (size_t ki=0; ki<sf_ix.size(); ++ki)
    {
      HedgeSurface* surf = surfaces_[sf_ix[ki]].get();
      vector<RevEngRegion*> reg = surf->getRegions();
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points0;
      BoundingBox bbox0(3);
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  points0.push_back(std::make_pair(reg[kj]->pointsBegin(),
					  reg[kj]->pointsEnd()));
	  bbox0.addUnionWith(reg[kj]->boundingBox());
	}
      Point low0 = bbox0.low();
      Point high0 = bbox0.high();
      Point pos0;
      double radius0;
      RevEngUtils::computeCylPosRadius(points0, low0, high0, axis, Cx, Cy, 
				      pos0, radius0);
      double radius2 = computeCylRadius(points0, pos, Cx, Cy);
      
      shared_ptr<Cylinder> cyl(new Cylinder(radius0, pos0, axis, Cy));
      shared_ptr<Cylinder> cyl2(new Cylinder(radius2, pos, axis, Cy));
      shared_ptr<Cylinder> cyl3(new Cylinder(radius0, pos, axis, Cy));
#ifdef DEBUG
      cyl->writeStandardHeader(of);
      cyl->write(of);
      cyl2->writeStandardHeader(of);
      cyl2->write(of);

      cyl2->writeStandardHeader(of);
      cyl2->write(of);
#endif
      for (size_t kj=0; kj<reg.size(); ++kj)
	{
	  double maxd, avd;
	  int num_in, num2_in;
	  vector<RevEngPoint*> in, out;
	  vector<pair<double,double> > distang;
	  vector<double> parvals;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl, approx_tol_, maxd, avd, num_in, num2_in, in, out,
				  parvals, distang);
#ifdef DEBUG
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in.size() << std::endl;
	  for (size_t kr=0; kr<in.size(); ++kr)
	    of << in[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out.size() << std::endl;
	  for (size_t kr=0; kr<out.size(); ++kr)
	    of << out[kr]->getPoint() << std::endl;
#endif
	  double maxd2, avd2;
	  int num_in2, num2_in2;
	  vector<RevEngPoint*> in2, out2;
	  vector<pair<double,double> > distang2;
	  vector<double> parvals2;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl2, approx_tol_, maxd2, avd2, num_in2, num2_in2,
				  in2, out2, parvals2, distang2);
#ifdef DEBUG
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in2.size() << std::endl;
	  for (size_t kr=0; kr<in2.size(); ++kr)
	    of << in2[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out2.size() << std::endl;
	  for (size_t kr=0; kr<out2.size(); ++kr)
	    of << out2[kr]->getPoint() << std::endl;
#endif
	  
	  double maxd3, avd3;
	  int num_in3, num2_in3;
	  vector<RevEngPoint*> in3, out3;
	  vector<pair<double,double> > distang3;
	  vector<double> parvals3;
	  RevEngUtils::distToSurf(reg[kj]->pointsBegin(), reg[kj]->pointsEnd(),
				  cyl3, approx_tol_, maxd3, avd3, num_in3, num2_in3,
				  in3, out3, parvals3, distang3);
#ifdef DEBUG
	  of << "400 1 0 4 155 50 50 255" << std::endl;
	  of << in3.size() << std::endl;
	  for (size_t kr=0; kr<in3.size(); ++kr)
	    of << in3[kr]->getPoint() << std::endl;
	  of << "400 1 0 4 50 155 50 255" << std::endl;
	  of << out3.size() << std::endl;
	  for (size_t kr=0; kr<out3.size(); ++kr)
	    of << out3[kr]->getPoint() << std::endl;
#endif
	  int stop_break = 1;
	}
    }
}

//===========================================================================
double RevEng::computeCylRadius(vector<pair<vector<RevEngPoint*>::iterator,
				vector<RevEngPoint*>::iterator> >& points,
				Point mid, Point vec1, Point vec2)
//===========================================================================
{
  int nmb = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    nmb += (int)(points[ki].second - points[ki].first);

  double wgt = 1.0/(double)nmb;

  // Make tranformation matrix
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D yaxis(0, 1, 0);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, yaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;

  // Rotate points and add to circle decription
  double r2 = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D pnt = pt->getPoint();
	  Vector3D rpnt = rotmat*pnt;
	  double curr = rpnt[0]*rpnt[0] + rpnt[1]*rpnt[1];
	  r2 += wgt*curr;
	}
    }

  double radius = sqrt(r2);
  return radius;
}

//===========================================================================
void RevEng::collectAxis(vector<SurfaceProperties>& sfprop)
//===========================================================================
{
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int code;
      ClassType type = surfaces_[ki]->instanceType(code);
      ClassType type2 = Class_Unknown;
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      shared_ptr<ElementarySurface> elemsf =
	dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);

      int nreg = surfaces_[ki]->numRegions();
      int num_pts = 0;
      shared_ptr<ParamSurface> primary;
      int num_pt_primary = 0;
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
	  int num = reg->numPoints();
	  num_pts += num;
	  if (reg->hasBaseSf() && num > num_pt_primary)
	    {
	      double maxdp, avdp;
	      int num_inp, num2_inp;
	      reg->getBaseDist(maxdp, avdp, num_inp, num2_inp);
	      if (num_inp > num/2 && avdp < approx_tol_)
		{
		  primary = reg->getBase();
		  num_pt_primary = num;
		}
	    }
	}
      if (primary.get())
	type2 = primary->instanceType();
      
      if (!elemsf.get())
	{
	  if (primary.get())
	    elemsf =
	      dynamic_pointer_cast<ElementarySurface,ParamSurface>(primary);
	}
      if (!elemsf.get())
	continue;
      
      Point loc, dir;
      double rad1, rad2;
      loc = elemsf->location();
      dir = elemsf->direction();
      rad1 = elemsf->radius(0.0, 0.0);   // Not good enough for cones
      rad2 = elemsf->radius2(0.0, 0.0);   // Not good enough for cones
      SurfaceProperties currprop((int)ki, type, num_pts, dir, loc, type2,
				 rad1, rad2);
      sfprop.push_back(currprop);
    }
}


//===========================================================================
void RevEng::trimSurfaces()
//===========================================================================
{
#ifdef DEBUG_TRIM
  std::ofstream of0("trimmedsfs.g2");
#endif  
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      if (!regions_[ki]->hasSurface())
	continue;
#ifdef DEBUG_TRIM
      std::ofstream of("trimreg.g2");
      regions_[ki]->writeRegionPoints(of);
#endif
      
      bool trimmed = regions_[ki]->trimSurface(approx_tol_);
#ifdef DEBUG_TRIM
      if (trimmed)
	{
	  shared_ptr<ParamSurface> tsurf = regions_[ki]->getSurface(0)->surface();
	  tsurf->writeStandardHeader(of0);
	  tsurf->write(of0);
	}
#endif
	
      int stop_break = 1;
    }
#ifdef DEBUG
   std::cout << "Finished trim surf, regions: " << regions_.size() << ", surfaces: " << surfaces_.size() << std::endl;
   if (regions_.size() > 0)
    {
      std::cout << "Regions13" << std::endl;
      std::ofstream of("regions13.g2");
      std::ofstream ofm("mid_regions13.g2");
      std::ofstream ofs("small_regions13.g2");
      writeRegionStage(of, ofm, ofs);
      std::ofstream of0("regions13_helix.g2");
      for (int ki=0; ki<(int)regions_.size(); ++ki)
	{
	  if (regions_[ki]->getSurfaceFlag() == PROBABLE_HELIX)
	    {
	      regions_[ki]->writeRegionInfo(of0);
	      if (regions_[ki]->hasSurface())
		regions_[ki]->writeSurface(of0);
	    }
	}
      if (surfaces_.size() > 0)
	{
	  std::ofstream of("regsurf13.g2");
	  writeRegionWithSurf(of);
	}
     }
   
   if (edges_.size() > 0)
     {
       std::ofstream ofe("edges13.g2");
       for (size_t kr=0; kr<edges_.size(); ++kr)
	 {
	   vector<shared_ptr<ParamCurve> > cvs = edges_[kr]->getSpaceCurves();
	   for (size_t kh=0; kh<cvs.size(); ++kh)
	     {
	       cvs[kh]->writeStandardHeader(ofe);
	       cvs[kh]->write(ofe);
	     }
	 }
     }
#endif
  
#if 0
#ifdef DEBUG
  std::ofstream of1("surfbd.g2");
#endif
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      // Restrict unbounded surfaces
      surfaces_[ki]->ensureSurfaceBounded();
#ifdef DEBUG
      surfaces_[ki]->surface()->writeStandardHeader(of1);
      surfaces_[ki]->surface()->write(of1);
#endif
    }

#ifdef DEBUG
  std::ofstream of3("trimsurfs.g2");
#endif
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      surfaces_[ki]->trimWithPoints(approx_tol_);

#ifdef DEBUG
      surfaces_[ki]->surface()->writeStandardHeader(of3);
      surfaces_[ki]->surface()->write(of3);
#endif
      int stop1 = 1;
    }
#endif
#if 0
  vector<shared_ptr<ParamSurface> > sub_sfs;
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      vector<shared_ptr<ParamSurface> > splitsfs =
	SurfaceModelUtils::checkClosedFaces(surf, 10.0*int_tol_);
      sub_sfs.insert(sub_sfs.end(), splitsfs.begin(), splitsfs.end());
    }
  
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs(surfaces_.size());
  vector<shared_ptr<BoundedSurface> > bd_sfs(surfaces_.size());
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf1 = surfaces_[ki]->surface();
      for (size_t kj=ki+1; kj<surfaces_.size(); ++kj)
  	{
	  if (surfaces_[ki]->isTangential(surfaces_[kj].get()))
	    continue;
	  shared_ptr<ParamSurface> sf2 = surfaces_[kj]->surface();

  	  // Intersect sf1 and sf2
  	  // Remember intersection curves
	  shared_ptr<BoundedSurface> bd1, bd2;
	  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
	  BoundedUtils::getSurfaceIntersections(sf1, sf2, int_tol_,
						int_cvs1, bd1, int_cvs2, bd2);
	  bd_sfs[ki] = bd1;
	  bd_sfs[kj] = bd2;
	  if (int_cvs1.size() > 0)
	    all_int_cvs[ki].insert(all_int_cvs[ki].end(), int_cvs1.begin(), int_cvs1.end());
	  if (int_cvs2.size() > 0)
	    all_int_cvs[kj].insert(all_int_cvs[kj].end(), int_cvs2.begin(), int_cvs2.end());
  	}
    }
  
  std::ofstream of2("intcvs.g2");
  for (size_t ki=0; ki<all_int_cvs.size(); ++ki)
    for (size_t kj=0; kj<all_int_cvs[ki].size(); ++kj)
      {
	shared_ptr<ParamCurve> cv = all_int_cvs[ki][kj]->spaceCurve();
	cv->writeStandardHeader(of2);
	cv->write(of2);
      }

  size_t nmb_sfs = surfaces_.size();
  for (size_t ki=0; ki<nmb_sfs; ++ki)
    {
      vector<shared_ptr<HedgeSurface> > added_sfs;
      surfaces_[ki]->doTrim(all_int_cvs[ki], bd_sfs[ki], int_tol_, added_sfs);
      if (added_sfs.size() > 0)
	surfaces_.insert(surfaces_.end(), added_sfs.begin(), added_sfs.end());
    }
  
  std::ofstream of3("trimsurfs.g2");
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      surfaces_[ki]->surface()->writeStandardHeader(of3);
      surfaces_[ki]->surface()->write(of3);
    }
#endif
  int stop_break = 1;
}

 //===========================================================================
shared_ptr<SurfaceModel> RevEng::createModel()
//===========================================================================
{
  // Ensure bounded surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = surfaces_[ki]->surface();
      if (!surf->isBounded())
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(0);
	  double dom[4];
	  reg->getDomain(dom);
	  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
	  shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
	  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ParamSurface>(surf);
	  if (plane.get())
	    plane->setParameterBounds(dom[0], dom[2], dom[1], dom[3]);
	  else if (cyl.get())
	    cyl->setParamBoundsV(dom[2], dom[3]);
	  else if (cone.get())
	    cone->setParamBoundsV(dom[2], dom[3]);
 	}
    }
  
  vector<shared_ptr<ftSurface> > tmpsfs(surfaces_.begin(), surfaces_.end());
  sfmodel_ = shared_ptr<SurfaceModel>(new SurfaceModel(approx_tol_, 10.0*int_tol_,
						       100*int_tol_, anglim_, 10*anglim_,
						       tmpsfs));
#ifdef DEBUG_MODEL
  int num_bd = sfmodel_->nmbBoundaries();
  if (num_bd > 0)
    {
      std::ofstream of1("model_boundaries.g2");
      for (int ka=0; ka<num_bd; ++ka)
	{
	  ftCurve bd = sfmodel_->getBoundary(ka);
	  int num_seg = bd.numSegments();
	  for (int kb=0; kb<num_seg; ++kb)
	    {
	      shared_ptr<ParamCurve> cv = bd.segment(kb).spaceCurve();
	      cv->writeStandardHeader(of1);
	      cv->write(of1);
	    }
	}
    }
  
#endif
  return sfmodel_;
}

 //===========================================================================
void RevEng::initParameters()
//===========================================================================
{
  // Set default parameters
  min_next_ = 10;  // Minimum number of neighbouring points
  max_next_ = 500; //std::min(80, tri_sf_->size()/200); //500;
  rfac_ = 6.0; //3.0;  // Factor for radius in which to search for neighbouring points
  cfac_ = 8.0;  // Edge points from curvature is given by
  // cfac_ times the average length of triangulation edges in a vertex
  norm_plane_lim_= 0.005; // Limit for when the cone angle corresponding
  // to triangle normals indicate an edge
  zero_H_ = 0.005; //0.001; //0.0001;  // When mean curvature is considered zero
  zero_K_ = 0.005; //0.001; //0.0001;  // When Gauss curvature is considered zero
  zero_si_ = 0.0075; //0.001; // When shape index is considered zero
  norm_ang_lim_ = 0.1*M_PI; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
  pca_lim_ = cness_lim_ = -1.0;
  min_point_region_ = 200; //50; //10;  // Should be updated with regard to the total
  // number of points
  approx_tol_ = 0.001;  // Very preliminary
  int_tol_ = 1.0e-6;
  anglim_ = 0.01;
  max_nmb_outlier_ = 3;
  rpix_ = 1;
  rpfac_ = 0.01;
  ffac_ = 0.01;
  sfac_ = 0.05;

  prefer_elementary_ = 1;
}

 //===========================================================================
void RevEng::updateParameters()
//===========================================================================
{
  if (model_character_ == SMOOTH)
    {
      rfac_ = 4.0;
    }
  else if (model_character_ == MEDIUM_ROUGH)
    {
      zero_H_ = 0.007;
      zero_K_ = 0.007;
    }
  else
    {
      rfac_ = 6.0;
      zero_H_ = 0.01;
      zero_K_ = 0.01;
      anglim_ = 0.02;
    }
}

 //===========================================================================
int RevEng::setSmallRegionNumber()
//===========================================================================
{
  vector<int> nmb_pt_reg(regions_.size());
  for (size_t ki=0; ki<regions_.size(); ++ki)
    nmb_pt_reg[ki] = regions_[ki]->numPoints();

  std::sort(nmb_pt_reg.begin(), nmb_pt_reg.end());
  int tot_num = tri_sf_->size();
  int num_reg = (int)regions_.size();
  int idel = tot_num/num_reg;
  int min_num = std::min(10, tot_num/20);
  int ixmax = (int)(0.99*num_reg);
  int Q4 = 3*num_reg/4;
  int max_num = std::max(min_num, nmb_pt_reg[ixmax]);
  int Q4_num = nmb_pt_reg[Q4];
  max_num = std::min(max_num, 10*idel);
  int ixdel = std::max(num_reg/100, 2);
  int prev = nmb_pt_reg[0], prev0 = 0;
  int ix;
  int fac = 2;
  int min_jump = std::min(idel, Q4_num); //2;
  for (ix=ixdel; ix<num_reg; ix+=ixdel)
    {
      int diff = nmb_pt_reg[ix] - prev;
      if (diff > fac*(std::max(min_jump, prev-prev0)))
	break;
      if (diff > 0)
	prev0 = prev;
      prev = nmb_pt_reg[ix];
    }
  ix = std::min(ix, ixmax);
      
  int num = std::max(min_num, std::min(nmb_pt_reg[ix], max_num/2));
  return num;
}


 //===========================================================================
void RevEng::checkConsistence(std::string text) const
//===========================================================================
{
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      vector<RevEngRegion*> adjacent;
      regions_[ki]->getAdjacentRegions(adjacent);
      for (size_t kj=0; kj<adjacent.size(); ++kj)
	{
	  size_t kr;
	  for (kr=0; kr<regions_.size(); ++kr)
	    if (adjacent[kj] == regions_[kr].get())
	      break;
	  if (kr == regions_.size())
	    std::cout << text << ", Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
	}
    }
  for (int ki=0; ki<(int)surfaces_.size(); ++ki)
    {
      int numreg = surfaces_[ki]->numRegions();
      for (int ka=0; ka<numreg; ++ka)
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
	  size_t kr;
	  for (kr=0; kr<regions_.size(); ++kr)
	    if (reg == regions_[kr].get())
	      break;
	  if (kr == regions_.size())
	    std::cout << text << ", surface 1. Obsolete region pointer, ki=" << ki << ", ka=" << ka << ". Region: " << reg << std::endl;
	  vector<RevEngRegion*> adjacent;
	  reg->getAdjacentRegions(adjacent);
	  for (size_t kj=0; kj<adjacent.size(); ++kj)
	    {
	      size_t kr;
	      for (kr=0; kr<regions_.size(); ++kr)
		if (adjacent[kj] == regions_[kr].get())
		  break;
	      if (kr == regions_.size())
		std::cout << text << ", surface. Obsolete region pointer, ki=" << ki << ", kj=" << kj << ". Region: " << adjacent[kj] << std::endl;
	    }
	}
    }
}

 //===========================================================================
void RevEng::storeClassified(ostream& os) const
//===========================================================================
{
  storeParams(os);
  int nmbpts = tri_sf_->size();
  os << nmbpts << std::endl;
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);
      pt->store(os);
    }
}

 //===========================================================================
void RevEng::readClassified(istream& is)
//===========================================================================
{
  readParams(is);
  int nmbpts;
  is >> nmbpts;
  tri_sf_ = shared_ptr<ftPointSet>(new ftPointSet());
  vector<vector<int> > next_ix(nmbpts);
  for (int ki=0; ki<nmbpts; ++ki)
    {
      shared_ptr<RevEngPoint> vertex(new RevEngPoint());
      vertex->read(is, zero_si_, next_ix[ki]);
      tri_sf_->addEntry(vertex);
    }

  // Add next information
  for (int ki=0; ki<nmbpts; ++ki)
    {
      ftSamplePoint* pt1 = (*tri_sf_)[ki];
      for (size_t kr=0; kr<next_ix[ki].size(); ++kr)
	{
	  int ix = next_ix[ki][kr];
	  ftSamplePoint* pt2 = (*tri_sf_)[ix];
	  pt1->addNeighbour(pt2);
	}
    }

  if (false)
    {
  for (int ki=0; ki<nmbpts; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]); 
      double rp[3];
      setRp(pt, rp);
      pt->setRp(rp);
    }
    }
  
  max_next_ = std::min(80, tri_sf_->size()/200);
  max_next_ = std::max(2*min_next_, max_next_);
  setBoundingBox();
}

 //===========================================================================
void RevEng::storeGrownRegions(ostream& os)
//===========================================================================
{
  storeClassified(os);
  os << surfaces_.size() << std::endl;
  for (int ka=0; ka<(int)surfaces_.size(); ++ka)
    {
      surfaces_[ka]->setId(ka);
      surfaces_[ka]->store(os);
    }
  os << regions_.size() << std::endl;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setId((int)ki);
      regions_[ki]->store(os);
    }

  os << single_points_.size() << std::endl;
  for (size_t ki=0; ki<single_points_.size(); ++ki)
    os << single_points_[ki]->getIndex() << " ";
  os << std::endl;
  
  os << edges_.size() << std::endl;
  for (int ka=0; ka<(int)edges_.size(); ++ka)
    {
      edges_[ka]->setId(ka);
      edges_[ka]->store(os);
    }
}

 //===========================================================================
void RevEng::readGrownRegions(istream& is)
//===========================================================================
{
  readClassified(is);
  int nmb = tri_sf_->size();
  vector<ftSamplePoint*> tmp_pts(nmb);
  for (int ka=0; ka<nmb; ++ka)
    tmp_pts[ka] = (*tri_sf_)[ka];
  std::set<ftSamplePoint*> tmp_set(tmp_pts.begin(), tmp_pts.end());
  std::cout << "Read grown, size1 = " << tmp_pts.size() << ", size2 = " << tmp_set.size() << std::endl;
  curvatureFilter();
  int num_sfs;
  is >> num_sfs;
  if (num_sfs > 0)
    surfaces_.resize(num_sfs);
  for (int ki=0; ki<num_sfs; ++ki)
    {
      surfaces_[ki] = shared_ptr<HedgeSurface>(new HedgeSurface());
      surfaces_[ki]->read(is);
    }
  
  int num_regions;
  is >> num_regions;
  regions_.resize(num_regions);
  for (int ki=0; ki<num_regions; ++ki)
    {
      vector<int> sf_id;
      regions_[ki] = shared_ptr<RevEngRegion>(new RevEngRegion(edge_class_type_));
      regions_[ki]->read(is, tri_sf_, sf_id);
      for (size_t kj=0; kj<sf_id.size(); ++kj)
	{
	  for (size_t kr=0; kr<surfaces_.size(); ++kr)
	    {
	      if (sf_id[kj] == surfaces_[kr]->getId())
		{
		  regions_[ki]->addHedge(surfaces_[kr].get());
		  surfaces_[kr]->addRegion(regions_[ki].get());
		  break;
		}
	    }
	}
    }

  for (int ki=0; ki<num_regions; ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }

  int num_single;
  is >> num_single;
  single_points_.resize(num_single);
  for (int ki=0; ki<num_single; ++ki)
    {
      int ix;
      is >> ix;
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ix]);
      single_points_[ki] = pt;
    }
  
  int num_edgs;
  is >> num_edgs;
  if (num_edgs > 0)
    edges_.resize(num_edgs);
  for (int ki=0; ki<num_edgs; ++ki)
    {
      edges_[ki] = shared_ptr<RevEngEdge>(new RevEngEdge());
      int id1, id2, id3;
      vector<int> blend_id;
      edges_[ki]->read(is, id1, id2, id3, blend_id);

      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int id = regions_[kr]->getId();
	  if (id == id1)
	    {
	      edges_[ki]->setReg1(regions_[kr].get());
	      regions_[kr]->addRevEdge(edges_[ki].get());
	      break;
	    }
	}
      
      for (size_t kr=0; kr<regions_.size(); ++kr)
	{
	  int id = regions_[kr]->getId();
	  if (id == id2)
	    {
	      edges_[ki]->setReg2(regions_[kr].get());
	      regions_[kr]->addRevEdge(edges_[ki].get());
	      break;
	    }
	}

      if (id3 >= 0)
	{
	  for (size_t kr=0; kr<regions_.size(); ++kr)
	    {
	      int id = regions_[kr]->getId();
	      if (id == id3)
		{
		  edges_[ki]->setBlendRegSurf(regions_[kr].get());
		  regions_[kr]->setBlendEdge(edges_[ki].get());
		  break;
		}
	    }
	}

      for (size_t kh=0; kh<blend_id.size(); ++kh)
	{
	  for (size_t kr=0; kr<regions_.size(); ++kr)
	    {
	      int id = regions_[kr]->getId();
	      if (blend_id[kh] == id)
		{
		  edges_[ki]->addBlendRegion(regions_[kr].get());
		  regions_[kr]->setAssociatedBlend(edges_[ki].get());
		  break;
		}
	    }
	}
    }
  
}
  

 //===========================================================================
void RevEng::storeParams(ostream& os) const
//===========================================================================
{
  os <<  model_character_ << " " << mean_edge_len_ << " " << min_next_ << " " << rfac_ << " " << cfac_;
  os << " " << pca_lim_ << " " << cness_lim_ << " " << norm_ang_lim_;
  os << " " << norm_plane_lim_ << " " << zero_H_ << " " << zero_K_;
  os << " " << zero_si_ << " " << min_point_region_ << " " << approx_tol_ ;
  os << " " << anglim_ << " " << max_nmb_outlier_ << " ";
  os << edge_class_type_ << " " << classification_type_ << std::endl;
  os << mainaxis_[0] << " " << mainaxis_[1] << " " << mainaxis_[2] << std::endl;
}

 //===========================================================================
void RevEng::readParams(istream& is)
//===========================================================================
{
  is >>  model_character_ >> mean_edge_len_ >> min_next_ >> rfac_ >> cfac_ >> pca_lim_;
  is >> cness_lim_ >> norm_ang_lim_ >> norm_plane_lim_ >> zero_H_ >> zero_K_ >> zero_si_;
  is >> min_point_region_ >> approx_tol_ >> anglim_ >> max_nmb_outlier_;
  is >> edge_class_type_ >> classification_type_;
  mainaxis_[0].resize(3);
  mainaxis_[1].resize(3);
  mainaxis_[2].resize(3);
  is >> mainaxis_[0] >> mainaxis_[1] >> mainaxis_[2];
}

 //===========================================================================
void RevEng::writeRegionWithSurf(ostream& of) const
//===========================================================================
{
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      if (regions_[kr]->hasSurface())
	{
	  regions_[kr]->writeRegionPoints(of);
	  regions_[kr]->writeSurface(of);
	}
    }
}
  
 //===========================================================================
void RevEng::writeRegionStage(ostream& of, ostream& ofm, ostream& ofs) const
//===========================================================================
{
  std::cout << "Num regions: " << regions_.size() << ", num surfaces: " << surfaces_.size() << " num edges: " << edges_.size() << std::endl;
  
  vector<Vector3D> small;
  int nmb_one = 0;
  int low = min_point_region_/10;
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      // BoundingBox bbox = regions_[kr]->boundingBox();
      // if (bbox.low().dist(bbox.high()) < 0.1)
      //   std::cout << "Small bounding box" << std::endl;
      // std::set<RevEngPoint*> tmpset(regions_[kr]->pointsBegin(), regions_[kr]->pointsEnd());
      // if (tmpset.size() != regions_[kr]->numPoints())
      // 	std::cout << "Point number mismatch. " << kr << " " << tmpset.size() << " " << regions_[kr]->numPoints() << std::endl;
      int nmb = regions_[kr]->numPoints();
      if (nmb < low && regions_[kr]->hasSurface() == false)
	{
	  for (int ki=0; ki<nmb; ++ki)
	    small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	}
      else if (nmb < min_point_region_ && regions_[kr]->hasSurface() == false)
	{
	  ofm << "400 1 0 0" << std::endl;
	  ofm << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      ofm << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
	    }
	  if (regions_[kr]->hasSurface())
	    regions_[kr]->writeSurface(ofm);
	}
      else
	{
	  if (nmb > 0)
	    {
	      of << "400 1 0 0" << std::endl;
	      of << nmb << std::endl;
	      for (int ki=0; ki<nmb; ++ki)
		{
		  of << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
		}
	    }
	  if (regions_[kr]->hasSurface())
	    regions_[kr]->writeSurface(of);
	}
      if (nmb == 1)
	nmb_one++;
    }
  std::cout << "Number of regions with one point: " << nmb_one << std::endl;
  ofs << "400 1 0 4 0 0 0 255" << std::endl;
  ofs << small.size() << std::endl;
  for (size_t kr=0; kr<small.size(); ++kr)
    ofs << small[kr] << std::endl;
}
