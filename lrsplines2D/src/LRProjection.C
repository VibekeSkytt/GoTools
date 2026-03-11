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

#include "GoTools/lrsplines2D/LRProjection.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/creators/SmoothSurf.h"
#include "sislP.h"
#include <cmath>
#include <iostream>
#include <fstream>

//#define DEBUG

using std::vector;
using std::set;
using std::map;
using std::array;
using std::cout;
using std::endl;
using namespace Go;

//==============================================================================
void LRProjection::computeCoef(LRSplineSurface *srf, 
			       LRBSpline2D *bspl,
			       int proj_type,
			       Point& coef)
//==============================================================================
{
  int del = (*bspl->supportedElementBegin())->getNmbValPrPoint();
  if (del == 0)
    del = srf->dimension() + 3;  // Parameter pair, point and distance
  if (proj_type == 1)
    {
      // Collect data points
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      for (auto el=bspl->supportedElementBegin(); el!=bspl->supportedElementEnd(); ++el)
	{
	  vector<double> elem_data = (*el)->getDataPoints();
	  data.insert(data.end(), elem_data.begin(), elem_data.end());
	  max_dist = std::max(max_dist, (*el)->getMaxError());
	  av_dist += (*el)->getAccumulatedError();
	  nmb_out += (*el)->getNmbOutsideTol();
	  nmb += (*el)->nmbDataPoints();
	}
      av_dist /= (double)nmb;
  

      int deg1 = bspl->degree(XFIXED);
      int deg2 = bspl->degree(YFIXED);
      if (nmb < deg1*deg2)
	{
	  coef = bspl->Coef();
	}
      else
	{
	  // Perform projection
	  double smoothwgt = 1.0e-3; //1.0e-9; 
	  TPproject(bspl, data, del, smoothwgt, coef);
	  Point coef2 = bspl->Coef();
	  double dist = coef.dist(coef2);
	  if (dist > 10.0*max_dist)
	    coef = coef2;
	}
    }
  else if (proj_type == 2)
    {
      // Radial basis functions with IDW
      int nmb_points = 15; //30;
      double rad;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      IDWproject(bspl, rad, data, del, coef);
    }
  else if (proj_type == 3)
    {
      // Bilinear approximation
      int nmb_points = 8;
      double rad;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      BilinProject(bspl, data, del, rad, coef);
    }
  else if (proj_type == 4)
    {
      // Quadratic approximation
      int nmb_points = 15;
      double rad;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      QuadProject(bspl, data, del, rad, coef);
    }
  else if (proj_type == 5)
    {
      // Cubic approximation
      int nmb_points = 30;
      double rad;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      CubicProject(bspl, data, del, rad, coef);
    }
  int stop_break = 1;
}

//==============================================================================
void LRProjection::RDataSet(LRSplineSurface *srf, LRBSpline2D *bspl, int nmb_pts,
			    double& rad, vector<double>& data,
			    int del, int& nmb, int& nmb_out,
			    double& max_dist, double& av_dist)
//==============================================================================
{
  rad = -1.0;
  Point uv0 = bspl->getGrevilleParameter();
  vector<Element2D*> elems;
  vector<LRBSpline2D*> bsplines;
  bsplines.push_back(bspl);
  vector<double>  uv_dist;
  double umin = uv0[0], umax = uv0[0], vmin = uv0[1], vmax = uv0[1];
  double sf_umin = srf->startparam_u();
  double sf_umax = srf->endparam_u();
  double sf_vmin = srf->startparam_v();
  double sf_vmax = srf->endparam_v();
  size_t start_el = 0;
  int nmb_found = 0;
  for (size_t ki=0; ki<bsplines.size(); ++ki)
    {
      for (auto el=bsplines[ki]->supportedElementBegin();
	   el!=bsplines[ki]->supportedElementEnd(); ++el)
	{
	  auto it = std::find(elems.begin(), elems.end(), *el);
	  if (it == elems.end())
	    elems.push_back(*el);
	}

      for (size_t kj=start_el; kj<elems.size(); ++kj)
	{
	  vector<double> elem_data = elems[kj]->getDataPoints();
	  int nd = (int)elem_data.size()/del;
	  for (int ka=0; ka<nd; ++ka)
	    {
	      double dist = Utils::distance_squared(uv0.begin(), uv0.end(),
						    &elem_data[ka*del]);
	      uv_dist.push_back(sqrt(dist));
	    }
	  umin = std::min(umin, elems[kj]->umin()); 
	  umax = std::max(umax, elems[kj]->umax()); 
	  vmin = std::min(vmin, elems[kj]->vmin()); 
	  vmax = std::max(vmax, elems[kj]->vmax()); 
	}

	  std::sort(uv_dist.begin(), uv_dist.end());
	  if ((int)uv_dist.size() > nmb_pts)
	    {
	      rad = uv_dist[nmb_pts];
	      if ((umin <= uv0[0]-rad || umin == sf_umin) &&
		  (umax >= uv0[0]+rad || umax == sf_umax) &&
		  (vmin <= uv0[1]-rad || vmin == sf_vmin) &&
		  (vmax >= uv0[1]+rad || vmax == sf_vmax))
		break;
	    }

      
	  for (size_t kj=start_el; kj<elems.size(); ++kj)
	    {
	      for (auto b2=elems[kj]->supportBegin(); b2!=elems[kj]->supportEnd(); ++b2)
		{
		  auto itb = std::find(bsplines.begin(), bsplines.end(), *b2);
		  if (itb == bsplines.end())
		    bsplines.push_back(*b2);
		}
	    }
	  start_el = elems.size();
    }
  
  for (size_t ki=0; ki<elems.size(); ++ki)
    {
      vector<double> elem_data = elems[ki]->getDataPoints();
      data.insert(data.end(), elem_data.begin(), elem_data.end());
      max_dist = std::max(max_dist, elems[ki]->getMaxError());
      av_dist += elems[ki]->getAccumulatedError();
      nmb_out += elems[ki]->getNmbOutsideTol();
      nmb += elems[ki]->nmbDataPoints();
    }

  av_dist /= (double)nmb;
}

//==============================================================================
void LRProjection::TPproject(LRBSpline2D *bspl,
			     vector<double>& data,
			     int del,
			     double smoothwgt,
			     Point& coef)
//==============================================================================
{
  // Fetch knot vectors of B-spline
  BSplineUniLR* uni_u = bspl->getUnivariate(XFIXED);
  BSplineUniLR* uni_v = bspl->getUnivariate(YFIXED);
  vector<double> knots1 = uni_u->getKnots();
  vector<double> knots2 = uni_v->getKnots();

  int startmult1 = uni_u->endmult(true);
  int endmult1 = uni_u->endmult(false);
  int startmult2 = uni_v->endmult(true);
  int endmult2 = uni_v->endmult(false);

  int order1 = uni_u->degree() + 1;
  int order2 = uni_v->degree() + 1;

  // Extend knot vectors to full multiplicity in endpoints
  double start1 = knots1[0];
  for (int ka=startmult1; ka<order1; ++ka)
    knots1.insert(knots1.begin(), start1);

  double end1 = knots1[knots1.size()-1];
  for (int ka=endmult1; ka<order1; ++ka)
    knots1.push_back(end1);

  double start2 = knots2[0];
  for (int ka=startmult2; ka<order2; ++ka)
    knots2.insert(knots2.begin(), start2);

  double end2 = knots2[knots2.size()-1];
  for (int ka=endmult2; ka<order2; ++ka)
    knots2.push_back(end2);

  // Define initial spline surface
  int in1 = (int)knots1.size() - order1;
  int in2 = (int)knots2.size() - order2;
  int dim = bspl->dimension();
  vector<double> init_coefs(in1*in2*dim, 0.0);
  shared_ptr<SplineSurface> init_surf(new SplineSurface(in1, in2, order1, order2,
							&knots1[0], &knots2[0],
							&init_coefs[0], dim));

  // Extract point and parameter data
  int num_points = (int)data.size()/del;
  vector<double> param, points;
  param.reserve(2*num_points);
  points.reserve(dim*num_points);
  for (int ka=0; ka<(int)data.size(); ka+=del)
    {
      for (int kb=0; kb<2; ++kb)
	param.push_back(data[ka+kb]);
      for (int kb=0; kb<dim; ++kb)
	points.push_back(data[ka+kb+2]);
    }

  if (num_points > in1*in2)
    smoothwgt *= 0.01;
  
  // Approximation
  SmoothSurf approx;
  int seem[2];
  seem[0] = seem[1] = 0;
  vector<int> coef_known(in1*in2, 0);
  approx.attach(init_surf, seem, &coef_known[0]);

  double wgt1 = 0.0;
  double wgt3 = (std::min(order1,order2) < 4) ? 0.0 : 0.5*smoothwgt;
  double wgt2 = smoothwgt - wgt3 - wgt1;
  approx.setOptimize(wgt1, wgt2, wgt3);

  double approx_wgt = 1.0 - wgt1 - wgt2 - wgt3;
  vector<double> pnt_wgt(num_points, 1.0);
  approx.setLeastSquares(points, param, pnt_wgt, approx_wgt);

  shared_ptr<SplineSurface> approx_surf;
  try {
  approx.equationSolve(approx_surf);
  }
  catch (...)
    {
      coef = bspl->Coef();
      return;
    }
  
  // Fetch coefficient corresponding to input B-spline
  int k1 = order1 - startmult1;
  int k2 = order2 - startmult2;
  int kk = k2*in1 + k1;
  coef = Point(approx_surf->coefs_begin()+kk*dim, approx_surf->coefs_begin()+(kk+1)*dim);
}

//==============================================================================
void LRProjection::IDWproject(LRBSpline2D *bspl, double rad,
			      vector<double>& data, int del, Point& coef)
//==============================================================================
{
  double eps = 1.0e-10;
  Point greville = bspl->getGrevilleParameter();
  int dim = bspl->dimension();
  Point nom(dim);
  nom.setValue(0.0);
  double denom = 0.0;
  int pp = 2;
  int nmbd = (int)data.size()/del;
  for (int ka=0; ka<nmbd; ++ka)
    {
      Point curr(&data[ka*del+2], &data[ka*del+2+dim]);
      double tmp = sqrt(Utils::distance_squared(greville.begin(), greville.end(),
						&data[ka*del]));
      if (tmp < eps)
	{
	  coef = curr;
	  return;
	}
      double tmp2 = std::max(0.0, rad-tmp)/(rad*tmp);
      tmp2 = tmp2*tmp2;

      nom += tmp2*curr;
      denom += tmp2;
    }

  coef = nom/denom;
}

//==============================================================================
void LRProjection::BilinProject(LRBSpline2D *bspl, vector<double>& data, int del,
				double rad, Point& coef)
//==============================================================================
{
  Point par = bspl->getGrevilleParameter();
  double u, v;
  int dim = bspl->dimension();
  vector<double> ls(16, 0.0), rs(4*dim, 0.0);
  vector<int> piv(4);
  double tmp[4];
  for (int ka=0; ka<4; ++ka)
    piv[ka] = ka;
  int nmbd = (int)data.size()/del;
  for (int ka=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      if (dist > rad)
	continue;
      u = data[ka*del];
      v = data[ka*del+1];
      tmp[0] = u;
      tmp[1] = v;
      tmp[2] = u*v;
      tmp[3] = 1.0;
      for (int kb=0; kb<4; ++kb)
	{
	  for (int kc=0; kc<4; ++kc)
	    ls[kc*4+kb] += tmp[kb]*tmp[kc];
	}
      for (int kb=0; kb<dim; ++kb)
	{
	  for (int kc=0; kc<4; ++kc)
	    rs[kc+kb*4] += data[ka*del+2+kb]*tmp[kc];
	}
    }

  coef = Point(dim);
  coef.setValue(0.0);
  int kstat = 0;
  s6lufacp(&ls[0], &piv[0], 4, &kstat);
  if (kstat < 0)
    return;

  for (int kb=0; kb<dim; ++kb)
    {
      s6lusolp(&ls[0], &rs[kb*4], &piv[0], 4, &kstat);
      if (kstat < 0)
	return;
    }

  for (int kb=0; kb<dim; ++kb)
    coef[kb] = rs[kb*4]*par[0] + rs[kb*4+1]*par[1] + rs[kb*4+2]*par[0]*par[1] + rs[kb*4+3];
}

//==============================================================================
void LRProjection::QuadProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef)
//==============================================================================
{
  Point par = bspl->getGrevilleParameter();
  double u, v;
  int dim = bspl->dimension();
  vector<double> ls(36, 0.0), rs(6*dim, 0.0);
  vector<int> piv(6);
  double tmp[6];
  for (int ka=0; ka<6; ++ka)
    piv[ka] = ka;
  int nmbd = (int)data.size()/del;
  for (int ka=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      if (dist > rad)
	continue;
      u = data[ka*del];
      v = data[ka*del+1];
      tmp[0] = u*u;
      tmp[1] = u*v;
      tmp[2] = v*v;
      tmp[3] = u;
      tmp[4] = v;
      tmp[5] = 1.0;
      for (int kb=0; kb<6; ++kb)
	{
	  for (int kc=0; kc<6; ++kc)
	    ls[kc*6+kb] += tmp[kb]*tmp[kc];
	}
      for (int kb=0; kb<dim; ++kb)
	{
	  for (int kc=0; kc<6; ++kc)
	    rs[kc+kb*6] += data[ka*del+2+kb]*tmp[kc];
	}
    }

  coef = Point(dim);
  coef.setValue(0.0);
  int kstat = 0;
  s6lufacp(&ls[0], &piv[0], 6, &kstat);
  if (kstat < 0)
    return;

  for (int kb=0; kb<dim; ++kb)
    {
      s6lusolp(&ls[0], &rs[kb*6], &piv[0], 6, &kstat);
      if (kstat < 0)
	return;
    }

  tmp[0] = par[0]*par[0];
  tmp[1] = par[0]*par[1];
  tmp[2] = par[1]*par[1];
  tmp[3] = par[0];
  tmp[4] = par[1];
  tmp[5] = 1.0;
  for (int kb=0; kb<dim; ++kb)
    for (int kc=0; kc<6; ++kc)
      coef[kb] += rs[kb*6+kc]*tmp[kc];
}

//==============================================================================
void LRProjection::CubicProject(LRBSpline2D *bspl, vector<double>& data, int del,
				double rad, Point& coef)
//==============================================================================
{
  Point par = bspl->getGrevilleParameter();
  double u, v;
  int dim = bspl->dimension();
  vector<double> ls(100, 0.0), rs(10*dim, 0.0);
  vector<int> piv(10);
  double tmp[10];
  for (int ka=0; ka<10; ++ka)
    piv[ka] = ka;
  int nmbd = (int)data.size()/del;
  for (int ka=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      if (dist > rad)
	continue;
      u = data[ka*del];
      v = data[ka*del+1];
      tmp[0] = u*u*u;
      tmp[1] = u*u*v;
      tmp[2] = u*v*v;
      tmp[3] = v*v*v;
      tmp[4] = u*u;
      tmp[5] = u*v;
      tmp[6] = v*v;
      tmp[7] = u;
      tmp[8] = v;
      tmp[9] = 1.0;
      for (int kb=0; kb<10; ++kb)
	{
	  for (int kc=0; kc<10; ++kc)
	    ls[kc*10+kb] += tmp[kb]*tmp[kc];
	}
      for (int kb=0; kb<dim; ++kb)
	{
	  for (int kc=0; kc<10; ++kc)
	    rs[kc+kb*10] += data[ka*del+2+kb]*tmp[kc];
	}
    }

  coef = Point(dim);
  coef.setValue(0.0);
  int kstat = 0;
  s6lufacp(&ls[0], &piv[0], 10, &kstat);
  if (kstat < 0)
    {
      coef = bspl->Coef();
      return;
    }

  for (int kb=0; kb<dim; ++kb)
    {
      s6lusolp(&ls[0], &rs[kb*10], &piv[0], 10, &kstat);
      if (kstat < 0)
	{
	  coef = bspl->Coef();
	  return;
	}
    }

  tmp[0] = par[0]*par[0]*par[0];
  tmp[1] = par[0]*par[0]*par[1];
  tmp[2] = par[0]*par[1]*par[1];
  tmp[3] = par[1]*par[1]*par[1];
  tmp[4] = par[0]*par[0];
  tmp[5] = par[0]*par[1];
  tmp[6] = par[1]*par[1];
  tmp[7] = par[0];
  tmp[8] = par[1];
  tmp[9] = 1.0;
  for (int kb=0; kb<dim; ++kb)
    for (int kc=0; kc<10; ++kc)
      coef[kb] += rs[kb*10+kc]*tmp[kc];
  int stop_break = 1;
}

