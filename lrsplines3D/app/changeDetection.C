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

#include "GoTools/lrsplines3D/ChangeDetection.h"
#include "GoTools/lrsplines3D/LRSpline3DBezierCoefs.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Utils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using namespace std;

int main(int argc, char *argv[])
{
  if (argc != 5)
    {
      std::cout << "Usage: list of file names, tol, max nmb iterations, derivative limit" << std::endl;
      return 1;
    }

  std::ifstream ifile(argv[1]);
  double tol = atof(argv[2]);
  int max_iter = atoi(argv[3]);
  double der_lim = atof(argv[4]);

    // Read point clouds 
  vector<string> point_file;
  while (!ifile.eof())
    {
      string file;
      ifile >> file;
      point_file.push_back(file);
      Utils::eatwhite(ifile);
    }

  // Sort point files
  std::sort(point_file.begin(), point_file.end());

  int dim;
  vector<vector<double> > point_seqs(point_file.size());
  for (size_t kj=0; kj<point_file.size(); ++kj)
    {
      std::ifstream inpts(point_file[kj]);
      ObjectHeader header;
      header.read(inpts);
      PointCloud3D cloud;
      cloud.read(inpts);
      int num = cloud.numPoints();
      dim = cloud.dimension();
      vector<double> curr(cloud.rawData(), cloud.rawData()+num*dim);
      point_seqs[kj] = curr;
    }

  ChangeDetection detect(point_seqs, dim);

  double tdel = detect.suggestTimeStep();

  // Possibility to adjust time step

  detect.defineTrivariate(tdel);

  // Volume approximation
  detect.volApprox(tol, max_iter);
  shared_ptr<LRSplineVolume> vol = detect.getApproxVol();
  
  ofstream of1("vol.g2");
  vol->writeStandardHeader(of1);
  vol->write(of1);
  LRSpline3DBezierCoefs bez(vol.get());
  bez.getBezierCoefs(0.0, 1, 2, true, 0.0);
  bez.writeToFile("vol.bb", 0);
  bez.writeToFile("vol_der.bb", 1);

  detect.identifyChanges(der_lim);

  detect.extractChangeData();

  double tol2 = 0.5*tol;
  int max_iter2 = std::max(4, (int)0.75*max_iter);
  detect.surfApprox(tol2, max_iter2);

  detect.differenceSurfaces();
  return 0;
}

  
  
