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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/Integrate.h"
#include "GoTools/geometry/Utils.h"
#include <vector>
#include <fstream>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 2)
    {
      std::cout << "Usage: Input file with closed sequence of 2D spline curves(.g2)" << std::endl;
      return -1;
    }

  std::ifstream is(argv[1]);

  vector<shared_ptr<SplineCurve> > curves;
  while (!is.eof())
    {
      ObjectHeader head;
      shared_ptr<SplineCurve> cv(new SplineCurve());
      is >> head;
      cv->read(is);
      curves.push_back(cv);
      Utils::eatwhite(is);
    }

  if (curves.size() == 0)
    return 0;

  // Check input
  double eps = 1.0e-6;
  size_t ki, kj;
  for (ki=0; ki<curves.size(); ++ki)
    {
      kj = (ki == curves.size()-1) ? 0 : ki+1;
      if (curves[ki]->dimension() != 2)
	{
	  std::cout << "Curve dimension different from two" << std::endl;
	  return -1;
	}

      Point pt1 = curves[ki]->ParamCurve::point(curves[ki]->endparam());
      Point pt2 = curves[kj]->ParamCurve::point(curves[kj]->startparam());
      if (pt1.dist(pt2) > eps)
	{
	  std::cout << "Curve sequence not continous/closed" << std::endl;
	  return -1;
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

  std::cout << "Area: " << area << std::endl;
}

  
      
