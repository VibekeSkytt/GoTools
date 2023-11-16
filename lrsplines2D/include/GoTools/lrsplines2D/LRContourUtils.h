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

#ifndef LR_CONTOURUTILS_H
#define LR_CONTOURTILS_H

#include <array>
#include <functional>
#include <set>
#include <map>
#include <unordered_map>
#include <tuple>
#include <iostream> // @@ debug

#include "GoTools/utils/checks.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/CurveLoop.h"

namespace Go
{

  /// Utility functionality used in computations involving trimmed
  // LR B-spline surfaces and corresponding contour curves
  namespace LRContourUtils
  {
    void completeContourLoop(shared_ptr<ParamCurve> crv, double isoval,
			     const CurveLoop& loop, double epsge,
			     std::vector<std::pair<std::vector<shared_ptr<ParamCurve> >, double> >& crv_loops,
			     std::vector<BoundingBox>& bbox);

    double computeLoopArea(std::vector<shared_ptr<SplineCurve> >& curves,
			   double eps);

   }; // end namespace LRContourUtils


}; // end namespace Go

#endif
