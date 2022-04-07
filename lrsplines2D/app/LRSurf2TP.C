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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

/// Translates an LRSplineSurface to a SplineSurface by extending knot line
/// segments that do not cover the entire domain of the LR B-spline surface.
/// Writes the LR mesh (mesh1.eps) and the tensor product mesh (mesh2.eps) to file 
/// if visualization by ghostview if requested.
int main( int argc, char* argv[] )
{
  if (argc != 3 && argc != 4) {
    std::cout << "Input parameters : Input file, output file, (mesh to file (0/1), default 0)"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream file2(argv[2]);
  int write_mesh = 0;
  if (argc == 4)
    write_mesh = atoi(argv[3]);

  // Read lrspline surface
  ObjectHeader header;
  header.read(file1);
  shared_ptr<LRSplineSurface> surf(new LRSplineSurface());
  surf->read(file1);

  if (write_mesh == 1)
    {
      std::ofstream ofmesh1("mesh1.eps");
      writePostscriptMesh(*surf, ofmesh1);
    }

  shared_ptr<SplineSurface> splsf(surf->asSplineSurface());
  splsf->writeStandardHeader(file2);
  splsf->write(file2);

  if (write_mesh == 1)
    {
      shared_ptr<LRSplineSurface> surf2(new LRSplineSurface(splsf.get(), 1.0e-8));
      std::ofstream ofmesh2("mesh2.eps");
      writePostscriptMesh(*surf2, ofmesh2);
    }

}

