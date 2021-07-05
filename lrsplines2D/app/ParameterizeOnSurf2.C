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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include <iostream>
#include <fstream>
#include <string.h>

//#define DEBUG

using namespace Go;
using std::vector;
using std::string;

void print_help_text()
{
  std::cout << "Purpose: Parameterize a point cloud by projection onto an LR B-spline surface. \n";
  std::cout << "Mandatory parameters: input surface (.g2), input point cloud (.txt, .xyz or .g2), output parameterized point cloud (.txt) \n";
  std::cout << "It is important that the shape of the surface reflects the terrain, but it does not need to be accurate.\n";
  std::cout << "Optional input parameter: \n";
  std::cout << "-info <filename> : Write accuracy information to file. Otherwise it is written to standard out. \n";
  std::cout << "-h or --help : Write this text\n";
}


int main(int argc, char *argv[])
{
  char* input_type = 0;    // Type of point file
  char *pointfile = 0;     // Input point file
  char *surffile = 0;      // Surface input file
  char *outfile = 0;       // Parameterized points output file
  char *infofile = 0;      // Accuracy information output file

  int ki, kj;
  vector<bool> par_read(argc-1, false);

  // Read optional parameters
  int nmb_par = argc-1;
  for (ki=1; ki<argc; ++ki)
    {
      string arg(argv[ki]);
      if (arg == "-h" || arg == "--help")
	{
	  print_help_text();
	  exit(0);
	}
      else if (arg == "-info")
	{
	  if (ki == argc-1)
	    {
	      std::cout << "ERROR: Missing input" << std::endl;
	      print_help_text();
	      return 1;
	    }
	  infofile = argv[ki+1];
	  par_read[ki-1] = par_read[ki] = true;
	  nmb_par -= 2;
	}
    }

  if (nmb_par != 3) 
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

  for (ki=1; ki<argc; ++ki)
    {
      if (par_read[ki-1])
	continue;
      if (nmb_par == 3)
	{
	  surffile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  pointfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 1)
	{
	  outfile = argv[ki];
	}
    }
  
  time_t tread = time(0);
  tm *ltm = localtime(&tread);
  std::cout << "Reading: Reading input: " << tread;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  // Read surface
  std::ifstream sfin(surffile);
  shared_ptr<LRSplineSurface> sf1;
  ObjectHeader header1;
  try {
    header1.read(sfin);
    sf1 = shared_ptr<LRSplineSurface>(new LRSplineSurface());
    sf1->read(sfin);
  }
  catch (...)
    {
      std::cout << "ERROR: Not a valid surface file" << std::endl;
      return 1;
    }

  if (sf1->dimension() != 3)
    {
      std::cout << "ERROR: 3D surface expected. Dimension = " << sf1->dimension() << std::endl;
      return 1;
    }

  // Read point cloud
  vector<double> data;
  vector<double> extent(6);   // Limits for points in all coordinates
  // Possible types of input files
  char keys[6][8] = {"g2", "txt", "TXT", "xyz", "XYZ", "dat"};
  int ptstype = FileUtils::fileType(pointfile, keys, 6);
  if (ptstype < 0)
    {
      std::cout << "ERROR: File type not recognized" << std::endl;
      return 1;
    }

  int nmb_pts = 0;
  std::ifstream pointsin(pointfile);
  if (ptstype == 0)
    {
     // g2 file
      ObjectHeader header;
      PointCloud3D points;
      try {
	header.read(pointsin);
	points.read(pointsin);
      }
      catch (...)
	{
	  std::cout << "ERROR: Not a valid point file" << std::endl;
	  return -1;
	}
      BoundingBox box = points.boundingBox();
      Point low = box.low();
      Point high = box.high();
      nmb_pts = points.numPoints();
      data.insert(data.end(), points.rawData(), points.rawData()+3*nmb_pts);
      for (int ki=0; ki<3; ++ki)
	{
	  extent[2*ki] = low[ki];
	  extent[2*ki+1] = high[ki];
	}
      std::cout << "Domain: [" << extent[0] << ", " << extent[1] << "] x [" << extent[2];
      std::cout << ", " << extent[3] << "]" << std::endl;
      std::cout << "Range: " << extent[4] << " - " << extent[5] << std::endl;
    }
  else
    FileUtils::readTxtPointFile(pointsin, 3, data, nmb_pts, extent);


  tread = time(0);
  ltm = localtime(&tread);
  std::cout << "End of reading: Points and surface read: " << tread;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  time_t tprocess = time(0);
  ltm = localtime(&tprocess);
  std::cout << "Processing: Parameterizing points: " << tprocess;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  
  // Move point cloud to origo
  Point mid(0.5*(extent[0]+extent[1]), 0.5*(extent[2]+extent[3]), 0.0);
  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=0; kj<3; ++kj)
      data[3*ki+kj] -= mid[kj];
  sf1->translate(-mid);
  
#ifdef DEBUG
  // Write translated surface and points to file in g2 format
  std::ofstream of("translated.g2");
  sf1->writeStandardHeader(of);
  sf1->write(of);
  PointCloud3D points(data.begin(), nmb_pts);
  points.writeStandardHeader(of);
  points.write(of);
#endif

  int dim = sf1->dimension();
  int maxiter = 4;
  double aeps = 0.001;

  double *curr;
  double dist;

  double maxdist = 0.0;
  double avdist = 0.0;


  // Represent the surface as tensor product
  shared_ptr<ParamSurface> tpsf(sf1->asSplineSurface());

  // For each point, project onto surface
  vector<double> param(2*nmb_pts);
  Point seed;
  try {
    for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
      {
	// // Get seed
	Point curr_pt(curr, curr+dim);
	// LRBSpline2D *bspline = LRSplineUtils::mostComparableBspline(sf1.get(), curr_pt);
	// Point seed = bspline->getGrevilleParameter();

	// Perform closest point
	double upar, vpar;
	Point close_pt;
	if (seed.dimension() > 0)
	  tpsf->closestPoint(curr_pt, upar, vpar, close_pt,
			     dist, aeps, maxiter, NULL, seed.begin());
	else
	  tpsf->closestPoint(curr_pt, upar, vpar, close_pt,
			     dist, aeps, maxiter);
	  seed.setValue(upar, vpar);

	maxdist = std::max(maxdist, dist);
	avdist += fabs(dist);

	param[2*ki] = upar;
	param[2*ki+1] = vpar;
      }
  }
  catch(...)
    {
      std::cout << "ERROR: Parameterization computation failed" << std::endl;
      return 1;
    }

  tprocess = time(0);
  ltm = localtime(&tprocess);
  std::cout << "End of processing: Parameterization completed: " << tprocess;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  time_t twrite = time(0);
  ltm = localtime(&twrite);
  std::cout << "Writing: Writing parameterized points to file: " << twrite;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  // Move point cloud to original position
  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=0; kj<3; ++kj)
      data[3*ki+kj] += mid[kj];

  std::ofstream os(outfile);
  (void)os.precision(15);
  for (ki=0; ki<nmb_pts; ++ki)
    {
      os << param[2*ki] << ", " << param[2*ki+1] << ", ";
      os << data[3*ki] << ", " << data[3*ki+1] << ", " << data[3*ki+2] << "\n";
    }

  avdist /= (double)nmb_pts;
  if (infofile)
    {
      std::ofstream infoout(infofile);   // Accuracy information output stream
      infoout << "Number of points: " << nmb_pts << std::endl;
      infoout << "Maximum distance between points and surface: " << maxdist << std::endl;
      infoout << "Average distance between points and surface: " << avdist << std::endl;
    }
  else
    {
      std::cout << "INFO: Number of points: " << nmb_pts << std::endl;
      std::cout << "INFO: Maximum distance between points and surface: " << maxdist << std::endl;
      std::cout << "INFO: Average distance between points and surface: " << avdist << std::endl;
    }

  twrite = time(0);
  ltm = localtime(&twrite);
  std::cout << "End of writing: Finished writing: " << twrite;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  return 0;
}

      




 
