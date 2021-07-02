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

#include "GoTools/lrsplines2D/LRSurfUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/RectDomain.h"

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

void LRSurfUtils::convert2TPsurfs(shared_ptr<ParamSurface> surf, double epsge,
				  double threshold_missing,
				  vector<shared_ptr<ParamSurface> >& tpsurfs)
{
  shared_ptr<LRSplineSurface> lrsurf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(surf);
  shared_ptr<BoundedSurface> bdsurf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  CurveBoundedDomain bddom;

  // Translate to parameter domain to origo
  RectDomain dom = surf->containingDomain();
  Point vec(-0.5*(dom.umin()+dom.umax()), -0.5*(dom.vmin()+dom.vmax()));
  TrimCrvUtils::translateSurfaceDomain(surf.get(), vec);

  if (bdsurf.get())
    {
      lrsurf = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(bdsurf->underlyingSurface());
      bddom = bdsurf->parameterDomain();
    }

  if (!lrsurf.get())
    {
      THROW("No LR B-spline surface is found");
    }
    
  CurveBoundedDomain *bddomptr = NULL;
  vector<CurveLoop> loops;
  if (bddom.nmbLoops() > 0)
    {
      bddomptr = &bddom;
      loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf, epsge);
    }
      
  //break LR surface up into individual patches
  const vector<pair<shared_ptr<LRSplineSurface>,LRSplineSurface::PatchStatus> > surf_fragments = 
    lrsurf->subdivideIntoSimpler(threshold_missing, epsge, bddomptr);

  // Fetch surfaces andapply trimming if necessary
  for (size_t ki=0; ki<surf_fragments.size(); ++ki)
    {
      LRSplineSurface::PatchStatus stat = surf_fragments[ki].second;
      if (stat == LRSplineSurface::OUTSIDE)
	continue;

      shared_ptr<ParamSurface> curr_sf;
      if (stat == LRSplineSurface::INSIDE)
	{
	  shared_ptr<ParamSurface> curr_sf =
	    shared_ptr<ParamSurface>(surf_fragments[ki].first->asSplineSurface());
	  tpsurfs.push_back(curr_sf);
	}
      else
	{
	  shared_ptr<SplineSurface> curr_tp(surf_fragments[ki].first->asSplineSurface());
	  shared_ptr<BoundedSurface> curr_bd(new BoundedSurface(curr_tp, epsge));
	  
	  // Collect trimming loop pieces which is inside the
	  // current surface
	  vector<shared_ptr<CurveOnSurface> > cv_pieces;
	  for (size_t kj=0; kj<loops.size(); ++kj)
	    {
	      int nmb = loops[kj].size();
	      for (int kr=0; kr<nmb; ++kr)
		{
		  // Fetch parameter curve
		  shared_ptr<ParamCurve> cv = loops[kj][kr];
		  shared_ptr<CurveOnSurface> sfcv =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
		  shared_ptr<ParamCurve> pcv;
		  if (cv->dimension() == 2)
		    pcv = cv;
		  else if (sfcv.get())
		    {
		      if (!sfcv->hasParameterCurve())
			sfcv->ensureParCrvExistence(epsge);
		      pcv = sfcv->parameterCurve();
		    }
		  else
		    continue;  // No parameter curve information
	      
		  // Make CurveOnSurface curve with respect to current
		  // surface
		  shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(curr_tp,
								      pcv, true));
	      
		  // Intersect curve with the trimming loop of the surface
		  vector<shared_ptr<CurveOnSurface> > trim_cvs = 
		    BoundedUtils::intersectWithSurface(*sfcv2, *curr_bd,
						       epsge, true);
		  cv_pieces.insert(cv_pieces.end(), trim_cvs.begin(), 
				   trim_cvs.end());
		}
	    }
	      
	  // Split surface with respect to trim segments
	  vector<shared_ptr<BoundedSurface> > trim_sfs;
	  if (cv_pieces.size() > 0)
	    trim_sfs =
	      BoundedUtils::splitWithTrimSegments(curr_bd, cv_pieces,
						  epsge);
	  
	  // Select inside surface
	  for (size_t kj=0; kj<trim_sfs.size(); ++kj)
	    {
	      double upar, vpar;
	      Point inner = trim_sfs[kj]->getInternalPoint(upar, vpar);
	      Vector2D parpt(upar, vpar);
	      if (bddom.isInDomain(parpt, epsge))
		tpsurfs.push_back(trim_sfs[kj]);
	    }
	}
    }

  // Translate back to original surface domain
  for (size_t ki=0; ki<tpsurfs.size(); ++ki)
    {
      TrimCrvUtils::translateSurfaceDomain(tpsurfs[ki].get(), -vec);
    }
}

