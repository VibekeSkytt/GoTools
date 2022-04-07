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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"

#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;

/// Apply a series of refinement to an input SplineSurface or LRSplineSurface.
/// Specification of refinement is performed interactively.

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: spline or lrspline_in (.g2),  lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read input surface
  ObjectHeader header;
  try {
    header.read(filein);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(filein);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no surface" << std::endl;
      exit(-1);
    }

  shared_ptr<LRSplineSurface> lrsf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sf);
  if (!lrsf.get())
    {
      shared_ptr<SplineSurface> splsf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
      if (splsf.get())
	lrsf = shared_ptr<LRSplineSurface>(new LRSplineSurface(splsf.get(), 1.0e-6));
    }
  if (!lrsf.get())
    {
      std::cerr << "Input file contains no spline surface" << std::endl;
      exit(-1);
    }
    
  
  // shared_ptr<LRSplineSurface> tmp2(lrsf->clone());
  // if (tmp2->dimension() == 1)
  //   tmp2->to3D();

  // tmp2->writeStandardHeader(fileout);
  // tmp2->write(fileout);
  // fileout << std::endl;
  // LineCloud lines2 = tmp2->getElementBds();
  // lines2.writeStandardHeader(fileout);
  // lines2.write(fileout);

  std::cout << "Degrees: " << lrsf->degree(XFIXED) << " " << lrsf->degree(YFIXED) << std::endl;
  std::cout << "Knots in the first parameter direction: " << std::endl;
  for (auto it=lrsf->mesh().knotsBegin(XFIXED); it!=lrsf->mesh().knotsEnd(XFIXED); ++it)
    std::cout << *it << " ";
  std::cout << std::endl;
  std::cout << "Knots in the second parameter direction: " << std::endl;
  for (auto it=lrsf->mesh().knotsBegin(YFIXED); it!=lrsf->mesh().knotsEnd(YFIXED); ++it)
    std::cout << *it << " ";
  std::cout << std::endl;
  
  int nmb_refs;
  std::cout << "Number of refinements: ";
  std::cin >> nmb_refs;
  std::cout << "For each refinement specify: direction (0/1), knot value, " << std::endl;
  std::cout << "start and end in the other direction, multiplicity" << std::endl;
  for (int ki=0; ki<nmb_refs; ++ki)
    {
      std::cout << "Refinement no. " << ki+1 << std::endl;
      double parval, start, end;
      int dir;
      int mult;

      std::cin >> dir;
      std::cin >> parval;
      std::cin >> start;
      std::cin >> end;
      std::cin >> mult;

      // check feasibility
      Direction2D d = (dir==0) ? XFIXED : YFIXED;
      if (parval < lrsf->paramMin(d) || parval > lrsf->paramMax(d))
	{
	  std::cout << "Parameter value outside surface domain" << std::endl;
	  continue;
	}
      int count = lrsf->mesh().numKnots(d, parval, start, end);
      if (count <= lrsf->degree(flip(d)) + 1)
	{
	  std::cout << "Refinement not legal" << std::endl;
	  continue;
	}
      //lrsf->refine((dir==0) ? XFIXED : YFIXED, parval, start, end, mult);
      lrsf->refine(d, parval, start, end, mult, true);
    }
  puts("Writing lr-spline to file.");
  // if (lrsf->dimension() == 1)
  // 	lrsf->to3D();
  lrsf->writeStandardHeader(fileout);
  lrsf->write(fileout);
  fileout << std::endl;

  return 0;
}
