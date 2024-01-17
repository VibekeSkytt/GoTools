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
//#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/SISLconversion.h"
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
//#define DEBUG_EDGE
//#define DEBUG_MONGE
//#define DEBUG_ENHANCE
//#define DEBUG_SEG
#define DEBUG
#define DEBUGONE
#define DEBUG_CHECK
//#define DEBUG_PLANAR
//#define DEBUG_AXIS
//#define DEBUG_GROW
#define DEBUG_VALIDATE

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
  
  int writepoints = 0;
  int nmbpt = tri_sf_->size();
  vector<double> tri_ang(nmbpt);
  for (int ki=0; ki<nmbpt; ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ki]);

      // Compute surface normal from triangulation
      pt->computeTriangNormal(100.0*mean_edge_len_);
      double avlen = pt->getMeanEdgLen();
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

      double avlen = pt->getMeanEdgLen();

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
  int writepoints = 0;
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
  int class_type = getClassificationType();
#ifdef DEBUG_ENHANCE  
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
      double avlen = pt->getMeanEdgLen();
      Point norm = pt->getTriangNormal();
      double ang = pt->getTriangAngle();
      Point Mnorm = pt->getMongeNormal();
      Point Pnorm = pt->getPCANormal();

      double fac = (pt->nmbMonge() == 0) ? 0.0 : 5.0;
#ifdef DEBUG_ENHANCE
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
	  double cfac = getCfac();
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
      else if (edge_class_type == RPFAC_EDGE)
	{
	  double RP_fac = getRPfac();
	  // std::cout << "Edge classification factor with RP: " << RP_fac << std::endl;
	  // std::cout << "New classification factor: " << std::endl;
	  // std::cin >> RP_fac;
	  // setRPfac(RP_fac);
	}
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
      double varh = 0.0;
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
  double zero_si = getShapeIndexZero();
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
  int dix = (int)(dlim*pointdist.size());
  double eps = pointdist[dix];
  if (model_character_ == MEDIUM_ROUGH)
    eps *= 1.5;
  else if (model_character_ == ROUGH)
    eps *= 2.0;

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
  int min_point_in = 50; //10; //20;  // Should be set depending on the total
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
  
}

//===========================================================================
void RevEng::initialSurfaces()
//===========================================================================
{
  // Create surfaces in simple regions (planes and cylinders) and extract
  // deviant points
  int min_point_in = 50; //10; //20;  // Should be set depending on the total
  // number of points. Probably class parameter. Need to look at the use
  double angtol = 5.0*anglim_;
  int regsize = (int)regions_.size();
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
  for (int kr=0; kr<regsize; ++kr)
    {
#ifdef DEBUG
      std::ofstream of0("init_reg.g2");
      regions_[kr]->writeRegionInfo(of0);
#endif
      
      if (regions_[kr]->numPoints() < min_point_region_)
	break;

      vector<vector<RevEngPoint*> > out_groups;
      vector<RevEngPoint*> single;
      vector<shared_ptr<HedgeSurface> > sfs;
      vector<HedgeSurface*> prev_sfs;
      bool repeat = false;
      regions_[kr]->initPlaneCyl(min_point_in, min_point_region_,
				 approx_tol_, angtol, mainaxis_,
				 sfs, out_groups, single, repeat);
      if (single.size() > 0)
	single_points_.insert(single_points_.end(), single.begin(), single.end());
      if (out_groups.size() > 0 || prev_sfs.size() > 0)
	surfaceExtractOutput(kr, out_groups, prev_sfs);

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
      if (tmpset.size() != regions_[ki]->numPoints())
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
      if (tmpset.size() != regions_[ki]->numPoints())
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
      if (tmpset.size() != regions_[ki]->numPoints())
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
	surfaceExtractOutput(ki, added_groups, dummy_surfs);
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
      if (tmpset.size() != regions_[ki]->numPoints())
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
  vector<int> reg_size(surfaces_.size());
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    reg_size[ki] = surfaces_[ki]->numPoints();

  std::sort(reg_size.begin(), reg_size.end());

  Point axis[3];
  int min_num = reg_size[(int)reg_size.size()/4];
  min_num = std::min(min_num, reg_size[reg_size.size()-1]/10);
  min_num = std::max(min_num, reg_size[reg_size.size()-1]/100);
  double max_ang = 0.1*M_PI;
  double angtol = 5.0*anglim_;

  Point plane_axis[3];
  int num_pts1[3];
  computeAxisFromPlane(mainaxis_, min_num, max_ang, plane_axis, num_pts1);

  Point cyl_axis[3];
  int num_pts2[3];
  computeAxisFromCylinder(plane_axis, min_num, max_ang, cyl_axis, num_pts2);

  // Update main axes. Prioritize information from planes
  for (int ka=0; ka<3; ++ka)
    {
      num_pts1[ka] += 2;
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
  mainaxis_[2] = mainaxis_[0].cross(mainaxis_[1]);
  mainaxis_[1] = mainaxis_[2].cross(mainaxis_[0]);

  for (int ka=0; ka<3; ++ka)
    mainaxis_[ka].normalize();

  // Update surfaces
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      bool updated = axisUpdate(surfaces_[ki].get(), max_ang, angtol);
      if (!updated)
	int stop_break0 = 1;
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
	    bool updated = axisUpdate(regions_[ki]->getSurface(0), max_ang,
				      angtol);
	    if (!updated)
	      int stop_break1 = 1;
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
      
  // Update adjacency between regions
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->clearRegionAdjacency();
    }
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      regions_[ki]->setRegionAdjacency();
    }
   int stop_break = 1;
}


//===========================================================================
bool RevEng::axisUpdate(HedgeSurface *hsurf, double max_ang, double angtol)
//===========================================================================
{
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(hsurf->surface());
  int ka = 3;
  if (elem.get())
    {
#ifdef DEBUG_AXIS
      std::ofstream of("axis_adapt.g2");
      elem->writeStandardHeader(of);
      elem->write(of);
      RevEngRegion *reg = hsurf->getRegion(0);
      reg->writeRegionPoints(of);
#endif
      Point vec = elem->direction();
      for (ka=0; ka<3; ++ka)
	{
	  double ang = vec.angle(mainaxis_[ka]);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < max_ang)
	    {
	      bool updated = 
		hsurf->updateSurfaceWithAxis(mainaxis_, ka,
					     approx_tol_, angtol);
	      if (updated)
		break;
	      return false;
	    }
	}
    }
  return (ka < 3) ? true : false;
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
	  segmented = segmentByPlaneGrow(ki, min_point_in, angtol);
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
	      segmented = segmentByAxis(ki, min_point_in);
	    }
	  if (!segmented)
	    segmented = segmentByContext(ki, min_point_in, angtol, true);
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
  for (size_t ki=0; ki<surfaces_.size(); ++ki)
    {
      int sfcode;
      ClassType type = surfaces_[ki]->instanceType(sfcode);
      if (type != Class_Cylinder)
	continue;
      
      shared_ptr<Cylinder> curr =
	dynamic_pointer_cast<Cylinder,ParamSurface>(surfaces_[ki]->surface());
      if (!curr.get())
	continue;

      int num_pts = surfaces_[ki]->numPoints();
      if (num_pts < min_num)
	continue;

      Point vec = curr->direction();
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = vec.angle(initaxis[ka]);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < max_ang)
	    {
	      int num_reg = surfaces_[ki]->numRegions();
	      for (int kb=0; kb<num_reg; ++kb)
		{
		  RevEngRegion *reg = surfaces_[ki]->getRegion(kb);
		  points[ka].push_back(std::make_pair(reg->pointsBegin(),
						      reg->pointsEnd()));
		}
	      num_points[ka] += num_pts;
	    }
	}
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
      
      double fac = 1.0/(double)num_pts;

      int num_reg = surfaces_[ki]->numRegions();
      double avdist = 0.0;
      for (int ka=0; ka<num_reg; ++ka)
	{
	  RevEngRegion *reg = surfaces_[ki]->getRegion(ka);
	  double tmp = reg->getAverageDist();
	  int nmb = reg->numPoints();
	  avdist += nmb*fac*tmp;
	}

      Point normal = curr->direction();
      for (int kb=0; kb<3; ++kb)
	{
	  double ang = normal.angle(initaxis[kb]);
	  ang = std::min(ang, M_PI-ang);
	  if (ang <= max_ang)
	    {
	      planes[kb].push_back(curr);
	      num[kb].push_back(num_pts);
	      avd[kb].push_back(avdist);
	      num_points[kb] += num_pts;
	    }
	}
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
	  if (normal*axis[kb] < 0.0)
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

      if (!found)
	found = regions_[ix]->extractTorus(approx_tol_, min_point_in,
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
      if (tmpset.size() != regions_[ki]->numPoints())
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
      if (tmpset.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, pre. " << ki << " " << tmpset.size() << " " << regions_[ki]->numPoints() << std::endl;
      
 #endif
    }
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
#ifdef DEBUG_CHECK
      std::set<RevEngPoint*> tmpset(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if (tmpset.size() != regions_[ki]->numPoints())
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
      int classtype = regions_[ki]->getClassification();
#ifdef DEBUGONE
      //std::cout << "No " << ki <<  ", classtype: " << classtype << std::endl;
      
      std::ofstream of1("region.g2");
      regions_[ki]->writeRegionInfo(of1);
      std::ofstream of2("unitsphere.g2");
      regions_[ki]->writeUnitSphereInfo(of2);
#endif
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
      if (tmpset2.size() != regions_[kj]->numPoints())
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
      if (tmpset.size() != regions_[ki]->numPoints())
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
      for (size_t ki=0; ki<regions_[ix]->numDivideInfo(); ++ki)
	{
	  segmented = regions_[ix]->divideWithSegInfo((int)ki,
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
      if (tmpset1.size() != regions_[ki]->numPoints())
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
      if (tmpset2.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 2. " << ki << " " << tmpset2.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
    
  for (size_t kr=0; kr<adj_planar.size(); ++kr)
    {
      std::set<RevEngPoint*> tmpset(adj_planar[kr]->pointsBegin(),
				    adj_planar[kr]->pointsEnd());
      if (tmpset.size() != adj_planar[kr]->numPoints())
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
	      if (tmpset4.size() != regions_[ki]->numPoints())
		std::cout << "Point number mismatch, composite 4. " << ki << " " << tmpset4.size() << " " << regions_[ki]->numPoints() << std::endl;
	    }
	  for (size_t kr=0; kr<adj_planar.size(); ++kr)
	    {
	      std::set<RevEngPoint*> tmpset(adj_planar[kr]->pointsBegin(),
					    adj_planar[kr]->pointsEnd());
	      if (tmpset.size() != adj_planar[kr]->numPoints())
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
      if (num > num_points/10)
	repeat = true;
    }
  
#ifdef DEBUG_CHECK
  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      std::set<RevEngPoint*> tmpset5(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd());
      if (tmpset5.size() != regions_[ki]->numPoints())
	std::cout << "Point number mismatch, composite 5. " << ki << " " << tmpset5.size() << " " << regions_[ki]->numPoints() << std::endl;
    }
#endif
#ifdef DEBUG
  if (regions_[ix]->numPoints() < num_points)
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
      if (tmpset6.size() != regions_[ki]->numPoints())
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
      if (tmpset6_2.size() != regions_[ki]->numPoints())
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
      if (tmpset7.size() != regions_[ki]->numPoints())
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
void RevEng::surfaceCreation(int pass)
//===========================================================================
{
  // First pass. Recognize elementary surfaces
  int min_point_in = 50; //10; //20;
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
      if (tmpset.size() != regions_[ka]->numPoints())
	std::cout << "Point number mismatch, pre grow. " << ka << " " << tmpset.size() << " " << regions_[ka]->numPoints() << std::endl;

      bool con = regions_[ka]->isConnected();
      if (!con)
	std::cout << "Disconnected region, ka= " << ka << " " << regions_[ka].get() << std::endl;
    }
#endif

  for (int ki=0; ki<(int)regions_.size(); ++ki)
    {
      if (regions_[ki]->hasSurface())
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
#endif
   
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
  regions_[ix]->growWithSurf(mainaxis_, min_nmb, min_point_region_,
			     approx_tol_, angtol,
			     grown_regions, adj_surfs, (pass>1));
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
	      if (kj == ix)
		{
		  ++kj;
		  continue;
		}

	      if (grown_regions[kr] == regions_[kj].get())
		{
		  regions_.erase(regions_.begin()+kj);
		  if (kj < ix)
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
  vector<shared_ptr<ftSurface> > tmpsfs(surfaces_.begin(), surfaces_.end());
  sfmodel_ = shared_ptr<SurfaceModel>(new SurfaceModel(approx_tol_, 10.0*int_tol_,
						       100*int_tol_, anglim_, 10*anglim_,
						       tmpsfs));
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
  int min_num = std::min(10, tot_num);
  int ixmax = (int)(0.99*num_reg);
  int max_num = std::max(min_num, nmb_pt_reg[ixmax]);
  max_num = std::min(max_num, 10*idel);
  int ixdel = std::max(num_reg/100, 2);
  int prev = nmb_pt_reg[0], prev0 = 0;
  int ix;
  int fac = 2;
  int min_jump = idel; //2;
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
      
  int num = std::max(min_num, std::min(nmb_pt_reg[ix], max_num));
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
    regions_[ki]->store(os);

  os << single_points_.size() << std::endl;
  for (size_t ki=0; ki<single_points_.size(); ++ki)
    os << single_points_[ki]->getIndex() << " ";
  os << std::endl;
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
  std::cout << "Num regions: " << regions_.size() << ", num surfaces: " << surfaces_.size() << std::endl;
  
  vector<Vector3D> small;
  int nmb_one = 0;
  int low = min_point_region_/4;
  for (size_t kr=0; kr<regions_.size(); ++kr)
    {
      // BoundingBox bbox = regions_[kr]->boundingBox();
      // if (bbox.low().dist(bbox.high()) < 0.1)
      //   std::cout << "Small bounding box" << std::endl;
      // std::set<RevEngPoint*> tmpset(regions_[kr]->pointsBegin(), regions_[kr]->pointsEnd());
      // if (tmpset.size() != regions_[kr]->numPoints())
      // 	std::cout << "Point number mismatch. " << kr << " " << tmpset.size() << " " << regions_[kr]->numPoints() << std::endl;
      int nmb = regions_[kr]->numPoints();
      if (nmb < low)
	{
	  for (int ki=0; ki<nmb; ++ki)
	    small.push_back(regions_[kr]->getPoint(ki)->getPoint());
	}
      else if (nmb < min_point_region_)
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
	  of << "400 1 0 0" << std::endl;
	  of << nmb << std::endl;
	  for (int ki=0; ki<nmb; ++ki)
	    {
	      of << regions_[kr]->getPoint(ki)->getPoint() << std::endl;
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
