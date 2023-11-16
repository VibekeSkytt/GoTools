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

#include "GoTools/lrsplines2D/LRContourUtils.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/Integrate.h"

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

//===========================================================================
void 
LRContourUtils::completeContourLoop(shared_ptr<ParamCurve> crv, double isoval,
				    const CurveLoop& parloop, double epsge,
				    vector<pair<vector<shared_ptr<ParamCurve> >, double> >& crv_loops,
				    vector<BoundingBox>& bbox)
//===========================================================================
  {
    // Compute distance between curve endpoints
    Point pos1 = crv->point(crv->startparam());
    Point pos2 = crv->point(crv->endparam());
    
    // Assign bounding boxes
    BoundingBox bb = crv->boundingBox();
    vector<shared_ptr<ParamCurve> > contour_loop;
    contour_loop.push_back(crv);
    double dist0 = pos1.dist(pos2);
    if (dist0 > epsge)
      {
	// Close the curve loop
	int ind1, ind2;
	double par1, par2, dist1, dist2;
	Point pt1, pt2;
	parloop.closestPoint(pos1, ind1, par1, pt1, dist1);
	parloop.closestPoint(pos2, ind2, par2, pt2, dist2);

	if (dist0 < std::min(dist1, dist2))
	    {
	      // Close original loop
	      shared_ptr<ParamCurve> crv2(new SplineCurve(pos1, pos2));
	      BoundingBox bb2 = crv2->boundingBox();
	      contour_loop.push_back(crv2);
	      bb.addUnionWith(bb2);
	      crv_loops.push_back(make_pair(contour_loop,isoval));
	      bbox.push_back(bb);
	    }
	  else
	    {
	      // Include part of surface loop in the curve loop
	      // Compute curve lengths
	      double len0 = crv->estimatedCurveLength();
	      int nmb = parloop.size();
	      if (ind1 > ind2 || (ind1 == ind2 && par1 > par2))
		{
		  std::swap(ind1, ind2);
		  std::swap(par1, par2);
		}
	      double len1 = 
		parloop[ind1]->estimatedCurveLength(par1, 
						    (ind2 != ind1) ? 
						    parloop[ind1]->endparam() : par2);
	      for (int ka=ind1+1; ka<ind2; ++ka)
		len1 += parloop[ka]->estimatedCurveLength();
	      if (ind1 != ind2)
		len1 += parloop[ind2]->estimatedCurveLength(parloop[ind2]->startparam(),
							 par2);

	      double len2 = 
		parloop[ind2]->estimatedCurveLength(par2, parloop[ind2]->endparam());
	      for (int ka=(ind2+1)%nmb; ka!=ind1; ka=(ka+1)%nmb)
		len2 += parloop[ka]->estimatedCurveLength();
	      len2 += parloop[ind1]->estimatedCurveLength(parloop[ind1]->startparam(),
						       par1);
	      
	      vector<shared_ptr<ParamCurve> > contour_loop2;
	      contour_loop2.push_back(shared_ptr<ParamCurve>(crv->clone()));
	      BoundingBox bb3 = bb;
	      if (len2 >= len0 + len1)
		{
		  for (int ka=ind1; ka<=ind2; ++ka)
		    {
		      double t1 = (ka == ind1) ? par1 : parloop[ka]->startparam();
		      double t2 = (ka == ind2) ? par2 : parloop[ka]->endparam();
		      
		      shared_ptr<ParamCurve> crv2;
		      if (ka != ind1 && ka != ind2) 
			crv2 = shared_ptr<ParamCurve>(parloop[ka]->clone());
		      else if (t2 - t1 >= epsge)
			crv2 =
			  shared_ptr<ParamCurve>(parloop[ka]->subCurve(t1, t2));
		      else
			continue;
		      BoundingBox bb2 = crv2->boundingBox();
		      contour_loop.push_back(crv2);
		      bb.addUnionWith(bb2);
		    }
		  crv_loops.push_back(make_pair(contour_loop,isoval));
		  bbox.push_back(bb);
		}

	      if (len1 >= len0 + len2)
		{
		  int ka, kb;
		  int nmb2 = nmb - (ind2-ind1);
		  for (ka=ind2, kb=0; kb<=nmb2; ka=(ka+1)%nmb, ++kb)
		    {
		      double t1 = (ka == ind2) ? par2 : parloop[ka]->startparam();
		      double t2 = (ka == ind1) ? par1 : parloop[ka]->endparam();
		      
		      shared_ptr<ParamCurve> crv2;
		      if (ka != ind1 && ka != ind2) 
			crv2 = shared_ptr<ParamCurve>(parloop[ka]->clone());
		      else if (t2 - t1 >= epsge)
			crv2 =
			  shared_ptr<ParamCurve>(parloop[ka]->subCurve(t1, t2));
		      else
			continue;
		      BoundingBox bb2 = crv2->boundingBox();
		      contour_loop2.push_back(crv2);
		      bb3.addUnionWith(bb2);
		    }
		  crv_loops.push_back(make_pair(contour_loop2,isoval));
		  bbox.push_back(bb3);
		}
	    }
	}
      else
	{
	  crv_loops.push_back(make_pair(contour_loop,isoval));
	  bbox.push_back(bb);
	}  
}

//===========================================================================
double LRContourUtils::computeLoopArea(vector<shared_ptr<SplineCurve> >& curves,
				       double eps)
//===========================================================================
{
  // Check input
  size_t ki, kj;
  for (ki=0; ki<curves.size(); ++ki)
    {
      kj = (ki == curves.size()-1) ? 0 : ki+1;
      if (curves[ki]->dimension() != 2)
	{
	  MESSAGE("Curve dimension different from two");
	  return 0.0;
	}

      Point pt1 = curves[ki]->ParamCurve::point(curves[ki]->endparam());
      Point pt2 = curves[kj]->ParamCurve::point(curves[kj]->startparam());
      if (pt1.dist(pt2) > eps)
	{
	  MESSAGE("Curve sequence not continous/closed");
	  return 0.0;
	}
    }
  
  double area = 0.0;
  for (ki=0; ki<curves.size(); ++ki)
    {
      // Integrate
      vector<double> param;
      vector<double> wgt;
      int order = curves[ki]->order();
      int deg = order*(order-1) - 1;
      vector<double> knots;
      curves[ki]->basis().knotsSimple(knots);
      for (size_t kr=1; kr<knots.size(); ++kr)
	{
	  GaussQuadValues(deg, knots[kr-1], knots[kr],
			  param, wgt);

	  double intval = 0.0;
	  for (size_t kj=0; kj<param.size(); ++kj)
	    {
	      vector<Point> der(2);
	      curves[ki]->point(der, param[kj], 1);
	      double val = -der[0][1]*der[1][0] + der[0][0]*der[1][1];
	      intval += 0.5*wgt[kj]*val;
	    }
	  area += (knots[kr]-knots[kr-1])*intval;
	}
    }

  return area;
}

  


