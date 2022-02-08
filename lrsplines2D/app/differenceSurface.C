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

#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/LRSurfUtils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: input surface 1(.g2), input surface 2 (.g2), output surface(.g2)" << std::endl;
    return -1;
  }

  std::ifstream is1(argv[1]);
  std::ifstream is2(argv[2]);
  std::ofstream out(argv[3]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

 // Read input surfaces
  ObjectHeader header;
  header.read(is1);
  shared_ptr<GeomObject> geom_obj1(Factory::createObject(header.classType()));
  geom_obj1->read(is1);
  
  shared_ptr<ParamSurface> sf1 = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj1);

  if ((!sf1.get()) || sf1->dimension() != 1)
    {
      std::cout << "First surface not found or dimension not equal to 1" << std::endl;
      return 1;
    }
 
  header.read(is2);
  shared_ptr<GeomObject> geom_obj2(Factory::createObject(header.classType()));
  geom_obj2->read(is2);
  
  shared_ptr<ParamSurface> sf2 = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj2);

  if ((!sf2.get()) || sf2->dimension() != 1)
    {
      std::cout << "Second surface not found or dimension not equal to 1" << std::endl;
      return 1;
    }
 
  shared_ptr<BoundedSurface> bdsf1 = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf1);
  shared_ptr<ParamSurface> psf1 = sf1;
  if (bdsf1.get())
    {
      psf1 = bdsf1->underlyingSurface();
    }

  shared_ptr<LRSplineSurface> lrsf1 = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(psf1);
 
  if (!lrsf1.get())
   {
     std::cout << "Input surface one is not of type LR B-spline:" << std::endl;
     return 1;
   }
 
  shared_ptr<BoundedSurface> bdsf2 = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf2);
  shared_ptr<ParamSurface> psf2 = sf2;
  if (bdsf2.get())
    {
      psf2 = bdsf2->underlyingSurface();
    }

  shared_ptr<LRSplineSurface> lrsf2 = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(psf2);
 
  if (!lrsf2.get())
   {
     std::cout << "Input surface two is not of type LR B-spline:" << std::endl;
     return 1;
   }

  // Ensure corresponding spline spaces
  vector<shared_ptr<LRSplineSurface> > sfs(2);
  sfs[0] = lrsf1;
  sfs[1] = lrsf2;
  LRSurfUtils::defineOnSameMesh(sfs);

  shared_ptr<LRSplineSurface> outlr = lrsf1;
  shared_ptr<LRSplineSurface> zerosf(new LRSplineSurface(*lrsf1));
  LRSplineSurface::BSplineMap::const_iterator it1 = outlr->basisFunctionsBegin();
  LRSplineSurface::BSplineMap::const_iterator it2 = lrsf2->basisFunctionsBegin();
  LRSplineSurface::BSplineMap::const_iterator it3 = zerosf->basisFunctionsBegin();
  Point zero(lrsf1->dimension());  // Initialized to zero
  for (; it1 != outlr->basisFunctionsEnd(); ++it1, ++it2, ++it3) 
    {
      Point c1 = it1->second->Coef();
      Point c2 = it2->second->Coef();
      Point c3 = c1 - c2;
      outlr->setCoef(c3, it1->second.get());
      zerosf->setCoef(zero, it3->second.get());
    }

  std::ofstream ofzero("zerosf.g2");
  zerosf->writeStandardHeader(ofzero);
  zerosf->write(ofzero);

  BoundedSurface *bdsf = NULL;
  if (bdsf1.get())
    bdsf = bdsf1.get();
  else if (bdsf2.get())
    bdsf = bdsf2.get();
  shared_ptr<ParamSurface> outsf;
  if (bdsf)
    {
     // The input surface is trimmed. Trim output surfaces accordingly
      vector<CurveLoop> loops = bdsf->allBoundaryLoops();
      vector<CurveLoop> loops1(loops.size());
      for (size_t ki=0; ki<loops.size(); ++ki)
	{
	  int nmb = loops[ki].size();
	  vector<shared_ptr<ParamCurve> > loop_cvs1(nmb);
	  for (int kj=0; kj<nmb; ++kj)
	    {
	      shared_ptr<ParamCurve> curr = loops[ki][kj];
	      shared_ptr<CurveOnSurface> sfcv = 
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr);
	      if (!sfcv.get())
		{
		  THROW("Trimming curve of wrong type");
		}
	      shared_ptr<ParamCurve> tmp_par1(sfcv->parameterCurve()->clone());
	      shared_ptr<ParamCurve> cv1(new CurveOnSurface(outlr, 
							    tmp_par1,
							    true));
	      loop_cvs1[kj] = cv1;
	    }
	  double eps = std::max(1.0e-12, loops[ki].getSpaceEpsilon());
	  loops1[ki].setCurves(loop_cvs1);
	  loops1[ki].setSpaceEpsilon(eps);
	}
      outsf = shared_ptr<ParamSurface>(new BoundedSurface(outlr, loops1));
      }
  else
    outsf = outlr;

  outsf->writeStandardHeader(out);
  outsf->write(out);
}
