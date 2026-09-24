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

#include "GoTools/lrsplines3D/SweepVolumeLR.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using std::vector;

typedef std::map<LRSplineSurface::BSKey, std::unique_ptr<LRBSpline2D> > BSplineMap; // storage of basis functions

//===========================================================================
LRSplineVolume* SweepVolumeLR::linearSweptVolume(const LRSplineSurface &surface,
						 const SplineCurve &curve,
						 const Point &pt)
//===========================================================================
{
  if (curve.rational())
    THROW("linearSweepLR: rational curve not supported");
  if (surface.rational())
    THROW("linearSweepLR: rational surface not supported");
  int dim = surface.dimension();
  if (curve.dimension() != dim || pt.dimension() != dim)
    THROW("linearSweepLR: dimension mismatch");

  
  // Define trivariate mesh
  Mesh2D mesh2d = surface.mesh();
  vector<double> cvknots(curve.knotsBegin(), curve.knotsEnd());
  Mesh3D mesh(mesh2d, cvknots);

  // Define B-splines
  // Define univariate B-splines in the x-direction
  int numuni1 = surface.numUnivariateBSplines(XFIXED);
  vector<std::unique_ptr<BSplineUniLR> > bsplinesuni1(numuni1);
  for (int u_ix=0; u_ix < numuni1; ++u_ix)
    {
      BSplineUniLR *curr = surface.getUnivariateBSpline(XFIXED, u_ix);
      bsplinesuni1[u_ix] =
	std::unique_ptr<BSplineUniLR>(new BSplineUniLR(1, curr->degree(),
						       &curr->kvec()[0],
						       &mesh));
      //bsplinesuni1[u_ix]->setCount(curr->getCount());
      int stop_break = 1;
    }
								 
  // Define univariate B-splines in the y-direction
  int numuni2 = surface.numUnivariateBSplines(YFIXED);
  vector<std::unique_ptr<BSplineUniLR> > bsplinesuni2(numuni2);
  for (int v_ix=0; v_ix < numuni2; ++v_ix)
    {
      BSplineUniLR *curr = surface.getUnivariateBSpline(YFIXED, v_ix);
      bsplinesuni2[v_ix] = 
	std::unique_ptr<BSplineUniLR>(new BSplineUniLR(2, curr->degree(),
								 &curr->kvec()[0],
								 &mesh));
      //bsplinesuni2[v_ix]->setCount(curr->getCount());
    }
								 
  // Create univariate B-splines in the z-direction
  int ncoefw = curve.numCoefs();
  int order = curve.order();
  vector<std::unique_ptr<BSplineUniLR> > bsplinesuni3(ncoefw);
  for (int w_ix = 0; w_ix != ncoefw; ++w_ix)
    {
      bsplinesuni3[w_ix] = 
	std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(3, order-1,
							       cvknots.begin() + w_ix,
							       &mesh)));
      //bsplinesuni3[w_ix]->incrCount();
    }

  // Create trivariate B-splines
  vector<std::unique_ptr<LRBSpline3D> > bsplvol(ncoefw*surface.numBasisFunctions());
  size_t b_ix = 0;
  for (BSplineMap::const_iterator bspl=surface.basisFunctionsBegin();
       bspl!=surface.basisFunctionsEnd(); ++bspl)
    {
      BSplineUniLR* uni1 = bspl->second->getUnivariate(XFIXED);
      int ix_uni1 = 0;
      bool found = BSplineUniUtils::identify_bsplineuni(uni1, bsplinesuni1, ix_uni1);
      if (!found)
	THROW("linearSweepLR: univariate B-spline not found");
	
      BSplineUniLR* uni2 = bspl->second->getUnivariate(YFIXED);
      int ix_uni2 = 0;
      found = BSplineUniUtils::identify_bsplineuni(uni2, bsplinesuni2, ix_uni2);
      if (!found)
	THROW("linearSweepLR: univariate B-spline not found");

      Point sfcoef = bspl->second->Coef();
      double gamma = bspl->second->gamma();
      int ix_uni3 = 0;
      for (auto it=curve.coefs_begin(); it!=curve.coefs_end(); it+=dim, ++ix_uni3)
	{
	  Point cvcoef(it, it+dim);
	  Point volcoef = cvcoef + sfcoef - pt;
	  bsplvol[b_ix++] =
	    std::unique_ptr<LRBSpline3D>(new LRBSpline3D(volcoef, 1.0,
							 bsplinesuni1[ix_uni1].get(),
							 bsplinesuni2[ix_uni2].get(),
							 bsplinesuni3[ix_uni3].get(),
							 gamma));
	}
    }
  
  // Define volume
  LRSplineVolume *sweepvol = new LRSplineVolume(mesh, bsplinesuni1, bsplinesuni2,
						bsplinesuni3, bsplvol,
						surface.getKnotTol());
  return sweepvol;
}



