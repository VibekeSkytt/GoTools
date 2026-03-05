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
  
  int del = (*bspl->supportedElementBegin())->getNmbValPrPoint();
  if (del == 0)
    del = srf->dimension() + 3;  // Parameter pair, point and distance

  int deg1 = bspl->degree(XFIXED);
  int deg2 = bspl->degree(YFIXED);
  if (proj_type == 1)
    {
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
      if (nmb < deg1*deg2)
	extendedDataSet(bspl, data, nmb, nmb_out, max_dist, av_dist);
      IDWproject(bspl, data, del, coef);
    }
  int stop_break = 1;
}

//==============================================================================
void LRProjection::extendedDataSet(LRBSpline2D *bspl,
				   vector<double>& data,
				   int& nmb, int& nmb_out,
				   double& max_dist, double& av_dist)
//==============================================================================
{
  std::set<Element2D*> elems;
  for (auto el=bspl->supportedElementBegin(); el!=bspl->supportedElementEnd(); ++el)
    {
      elems.insert(*el);
      for (auto b2=(*el)->supportBegin(); b2!=(*el)->supportEnd(); ++b2)
	{
	  for (auto el2=(*b2)->supportedElementBegin(); el2!=(*b2)->supportedElementEnd(); ++el2)
	    elems.insert(*el2);
	}
    }

  for (auto el=elems.begin(); el!=elems.end(); ++el)
    {
      vector<double> elem_data = (*el)->getDataPoints();
      data.insert(data.end(), elem_data.begin(), elem_data.end());
      max_dist = std::max(max_dist, (*el)->getMaxError());
      av_dist += (*el)->getAccumulatedError();
      nmb_out += (*el)->getNmbOutsideTol();
      nmb += (*el)->nmbDataPoints();
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
void LRProjection::IDWproject(LRBSpline2D *bspl,
			      vector<double>& data,
			      int del,
			      Point& coef)
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
      double tmp = Utils::distance_squared(greville.begin(), greville.end(), &data[ka*del]);
      double tmp2 = 1.0/pow(tmp, pp);
      if (tmp < eps)
	{
	  coef = curr;
	  return;
	}

      nom += tmp2*curr;
      denom += tmp2;
    }

  coef = nom/denom;
}

