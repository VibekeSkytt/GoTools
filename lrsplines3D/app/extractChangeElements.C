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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRSpline3DBezierCoefs.h"
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace Go;
using namespace std;

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
      std::cout << "Usage: lr_spline_vol.g2, tri_variate_data.txt, derivative_limit, bezier_binary.bb, bezier_binary_der.bb" << std::endl;
      return 1;
  }

  std::ifstream ifvol(argv[1]);
  std::ifstream ifpts(argv[2]);
  double der_lim = atof(argv[3]);
  std::string of1(argv[4]);
  std::string of2(argv[5]);
  
  ObjectHeader oh;
  oh.read(ifvol);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifvol);
  if (!vol.get())
    {
      std::cout << "Missing volume" << std::endl;
      return 1;
    }

   // Read data points
  int nmb_pts;
  ifpts >> nmb_pts;

  vector<double> pc4d;

  for (int ix=0; ix!=nmb_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifpts >> p0 >> p1 >> p2 >> q0;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
    }

  // Distribute points to elements
  LRSpline3DUtils::distributeDataPoints(vol.get(), pc4d); 
  
  std::cout << "Finished reading" << std::endl;
  LRSpline3DBezierCoefs bez(vol.get());

  bez.getBezierCoefs(0.0, 1, 2, true, der_lim);
  bez.writeToFile(of1, 0);
  bez.writeToFile(of2, 1);
  
  return 0;
}

  

