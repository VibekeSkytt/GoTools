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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRApproxApp.h"

#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
///  
/// This programs reads a point cloud (.txt), parameterizes the points by
/// their x- and y-coordintes and approximates the z-coordinate with an
/// LR B-spline surface. Appropriate parameters are defined.
///
/// The application gotools/lrsplines2D/app/PointCloud2LR performs the
/// corresponding functionality with user defined point cloud and
/// parameters. This application can also support a more sophisticated
/// parameter setting.
///
/// Input/output
/// The input point cloud is read from a hardcoded file. Parameters are
/// hardcoded. The produced surface is written to the file data/surface.g2
///
/// Note that the execution must be started from
/// <build directory>/gotools/lrsplines2D/examples
/// to get access to the data directory.
///
//===========================================================================
int main( int argc, char* argv[] )
{
  // Prepare for input and output file
  std::string infile("data/a_fBn_std0.002.txt");
  std::string outfile("data/surface.g2");

  // Prepare for input data
  std::ifstream input(infile.c_str());

  // Read the point cloud
  // The file type is .txt so the expected format is:
  // Number of points (first line)
  // x y z (all subsequent lines)
  int del = 3;  // Number of entries for each point
  int nmb_pts;  // Number of data points to be read from the file
  vector<double> data;  // Data points to be read from the file
  vector<double> extent(2*del);   // Limits for points in all coordinates
  std::cout << "Read point cloud" << std::endl;
  FileUtils::readTxtPointFile(input, del, data, nmb_pts, extent);

  // Define parameters
  int degree = 2;  // This choice is found to be appropriate for scattered
  // data with noise and unsmooth features, it gives a smooth surface with
  // flexibility
  int order = degree + 1;
  int nmb_coef = 5;  // Gives two inner knots in each parameter direction
  // initally, the minimum number for nmb_coef is order
  bool initmba = false;   // The initial surface will be created with
  // least squares approximation
  double smoothwgt = 1.0e-9; // Weight on smoothing term in least squares
  // approximation. Kept low to emphasize the approximation
  int max_iter = 10;   // The maximum number of steps in the adaptive surface
  // refinement and approximation
  int tomba = 5;       // The iteration level where mba replaces least
  // squares. The switch may also be done earlier by the algorithm
  double epsilon = 0.007;  // Approximation tolerance, used to identify
  // need for refinement. The tolerance must be set according to the
  // magnitude of the points and the noise level. 

  // Define approximation engine
  std::cout << "Define approximatin engine" << std::endl;
  LRSurfApprox approx(nmb_coef, order, nmb_coef, order, data, del-2,
		      epsilon, initmba);

  // Attach parameters
  approx.setSmoothingWeight(smoothwgt);
  approx.setSwitchToMBA(tomba);
  approx.setVerbose(true);   // Writes information about accuracy 
  // during the iterative approximation

  // Define a bound for the surface coefficients based on the range of
  // z-values in the data set
  double zrange = extent[5] - extent[4];
  double zfac = std::max(epsilon, 0.005*zrange);
  approx.addLowerConstraint(extent[4] - zfac);
  approx.addUpperConstraint(extent[5] + zfac);
  approx.setLocalConstraint(zfac);

  // Run approximation and fetch accuracy information
  double maxdist;  // Maximum distance between the surface and the
  // point cloud in z-direction
  double avdist_out;  // Average distance in points with a distance
  // larger than the tolerance (epsilon)
  double avdist_total; // Average distance in all points
  int nmb_out_eps;     // Number of points with a distance larger than
  // epsilon
  shared_ptr<LRSplineSurface> surf;  // Resulting surface (function)
  std::cout << "Perform approximation" << std::endl;
  try {
  surf = approx.getApproxSurf(maxdist, avdist_total,
			      avdist_out, nmb_out_eps, 
			      max_iter);
  }
  catch (...)
    {
      std::cout << "ERROR: Surface approximation failed" << std::endl;
      return 1;
    }

  std::cout << "Approximation completed" << std::endl;
  std::cout << "Total number of points: " << nmb_pts << std::endl;
  std::cout << "Number of elements: " << surf->numElements() << std::endl;
  std::cout << "Number of coefficients: " << surf->numBasisFunctions() << std::endl;
  std::cout << "Maximum distance: " << maxdist << std::endl;
  std::cout << "Average distance: " << avdist_total << std::endl;
  std::cout << "Average distance for points outside of the tolerance: " << avdist_out << std::endl;
  std::cout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl;

  std::cout << "Write surface to file" << std::endl;
  // Write to file
  // The surface can be view with the app goview_vol_and_lr
  std::ofstream sfout(outfile.c_str());
  surf->writeStandardHeader(sfout);
  surf->write(sfout);
  return 0;
}

  
