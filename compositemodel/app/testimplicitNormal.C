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

#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "sisl.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
  if (argc != 4)
    {
    std::cout << "Input parameters : Input file, offset dist, degree"  << std::endl;
    exit(-1);
    }
  
  // Read the points and normals from file
  ifstream input(argv[1]);
  double offdist = atof(argv[2]);
  int degree = atoi(argv[3]);

  if ( !input)
    {
      std::cerr << "Cannot open file for input\n";
      return 1;
    }

  // Read data
  double val;
  vector<double> points_extended;
  while (input >> val)
    {
      points_extended.push_back(val);
    }

  // Extract point data
  vector<double> points3D(points_extended.size()/5);
  vector<double> offpt3D(points_extended.size()/5);
  size_t kj = 0;
  for (size_t ki=0; ki<points_extended.size(); ki+=15)
    {
      Point norm(points_extended[ki+3], points_extended[ki+4], points_extended[ki+5]);
      norm.normalize();
      for (int ka=0; ka<3; ++ka, ++kj)
	{
	  points3D[kj] = points_extended[ki+ka];
	  offpt3D[kj] = points3D[kj] + offdist*norm[ka];
	}
    }  
  PointCloud3D cloud(&points3D[0], points3D.size()/3);
  PointCloud3D offcloud(&offpt3D[0], offpt3D.size()/3);

  std::ofstream of("IPA_cloud.g2");
  of << "400 1 0 4 0 0 255 255" << std::endl;
  cloud.write(of);
  of << "400 1 0 4 0 255 0 255" << std::endl;
  offcloud.write(of);
  
  // Implicitize
  ImplicitizePointCloudAlgo implicitize(cloud, degree);
  implicitize.perform();
  
  // Get result
  BernsteinTetrahedralPoly implicit;
  BaryCoordSystem3D bc;
  double sigma_min;
  implicitize.getResultData(implicit, bc, sigma_min);

    // Differentiate
    Vector4D bdir1(1.0, 0.0, 0.0, 0.0);
    Vector4D bdir2(0.0, 1.0, 0.0, 0.0);
    Vector4D bdir3(0.0, 0.0, 1.0, 0.0);
    Vector4D bdir4(0.0, 0.0, 0.0, 1.0);
    BernsteinTetrahedralPoly deriv1, deriv2, deriv3, deriv4;
    implicit.deriv(1, bdir1, deriv1);
    implicit.deriv(1, bdir2, deriv2);
    implicit.deriv(1, bdir3, deriv3);
    implicit.deriv(1, bdir4, deriv4);

    // Check accuracy
  double avdist = 0.0;
  double avdist2 = 0.0;
  double maxdist = 0.0;
  int dim = 3;
  double avdist_2 = 0.0;
  double avdist2_2 = 0.0;
  double maxdist_2 = 0.0;
  int numpt = cloud.numPoints();
  int nb = (degree+1)*(degree+2)*(degree+3)/6;
  for (int ki=0; ki<numpt; ++ki)
    {
      Vector3D curr = cloud.point(ki);
      Vector4D bary = bc.cartToBary(curr);
      double dist = implicit(bary);

      Vector3D curr2 = offcloud.point(ki);
      Vector4D bary2 = bc.cartToBary(curr2);
      Vector4D vec = bary2 - bary;

      double d1 = deriv1(bary);
      double d2 = deriv2(bary);
      double d3 = deriv3(bary);
      double d4 = deriv4(bary);
      Vector4D grad(d1, d2, d3, d4);
      //Vector4D grad(d2, d3, d4, d1);
      double angle = vec.angle(grad);

      Vector3D grad2 = bc.baryToCart(grad);
      Vector3D vec2 = curr2 - curr;
      double angle2 = vec2.angle(grad2);
      int stop_break = 1;
    }
}
