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

#include <fstream>
#include <iostream>
#include <chrono>

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfUtils.h"
#include "GoTools/geometry/SurfaceTools.h"


using namespace std;
using namespace Go;


int main(int varnum, char* vararg[])
{

  if (varnum != 4 && varnum != 5)
    {
      std::cout << "<Surface in (.g2)> <Surfaces out (.g2)> <maximum number of missing knots at finish> (<3D output surfaces>) " << std::endl;
      exit(-1);
    }
  ifstream is(vararg[1]);
  ofstream os(vararg[2]);
  int threshold_missing = atoi(vararg[3]);
  double epsge = 1.0e-6;

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  ObjectHeader header;
  header.read(is);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  shared_ptr<ParamSurface> surf =
    dynamic_pointer_cast<ParamSurface,GeomObject>(geom_obj);
  if (!surf.get())
    {
      std::cout << "Object is not a surface" << std::endl;
      exit(1);
    }
  surf->read(is);
  is.close();
  
  // Split into tensor-produce surfaces
  vector<shared_ptr<ParamSurface> > tpsurfs;
  LRSurfUtils::convert2TPsurfs(surf, epsge, threshold_missing, tpsurfs);

  for (size_t ki=0; ki<tpsurfs.size(); ++ki)
    {
      tpsurfs[ki]->writeStandardHeader(os);
      tpsurfs[ki]->write(os);
    }

  if (tpsurfs.size() > 0 && tpsurfs[0]->dimension() == 1 && varnum == 5)
    {
      ofstream os2(vararg[4]);
      for (size_t ki=0; ki<tpsurfs.size(); ++ki)
	{
	  // Make 3D surface
	  shared_ptr<SplineSurface> ssf =
	    dynamic_pointer_cast<SplineSurface,ParamSurface>(tpsurfs[ki]);
	  shared_ptr<BoundedSurface> bdsf =
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(tpsurfs[ki]);
	  if (bdsf.get())
	    ssf = dynamic_pointer_cast<SplineSurface,ParamSurface>(bdsf->underlyingSurface());
	  if (!ssf.get())
	    continue;

	  int dim = ssf->dimension();
	  int kdim = dim;
	  std::vector<double>::iterator cfstart;
	  vector<double> coefs;
	  BsplineBasis& basisu = ssf->basis_u();
	  BsplineBasis& basisv = ssf->basis_v();
	  int n1 = ssf->numCoefs_u();
	  int n2 = ssf->numCoefs_v();
	  if (ssf->rational())
	    {
	      cfstart = ssf->rcoefs_begin();
	    }
	  else
	    {
	      cfstart = ssf->coefs_begin();
	    }

	  for (int k2=0; k2<n2; ++k2)
	    {
	      double v = basisv.grevilleParameter(k2);
	      for (int k1=0; k1<n1; ++k1)
		{
		  double u = basisu.grevilleParameter(k1);
		  coefs.push_back(u);
		  coefs.push_back(v);
		  for (int kr=0; kr<kdim; ++kr, cfstart++)
		    coefs.push_back(*cfstart);
		}
	    }
	  shared_ptr<SplineSurface> ssf2(new SplineSurface(basisu, basisv, &coefs[0],
							   3, ssf->rational()));

	  shared_ptr<ParamSurface> tmptp;
	  if (!bdsf.get())
	    tmptp = ssf2;
	  else
	    {
	      // Transfer trimming information
	      vector<CurveLoop> loops =
		SurfaceTools::absolutelyAllBoundarySfLoops(bdsf, epsge);
	      vector<vector<shared_ptr<CurveOnSurface> > > loop_cvs2(loops.size());
	      for (size_t kj=0; kj<loops.size(); ++kj)
		{
		  int nmb = loops[kj].size();
		  for (int kr=0; kr<nmb; ++kr)
		    {
		      // Fetch parameter curve
		      shared_ptr<ParamCurve> cv = loops[kj][kr];
		      shared_ptr<CurveOnSurface> sfcv =
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);

		      shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(ssf2, sfcv->parameterCurve(), true));
		      loop_cvs2[kj].push_back(sfcv2);
		    }
		}
	      shared_ptr<BoundedSurface> bdsf2(new BoundedSurface(ssf2, loop_cvs2, bdsf->getEpsGeo()));
	      tmptp = bdsf2;
	    }
	  
	  tmptp->writeStandardHeader(os2);
	  tmptp->write(os2);
	}
    }
}
