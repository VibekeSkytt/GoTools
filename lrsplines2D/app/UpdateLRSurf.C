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

#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <getopt.h>

//#define DEBUG

using namespace Go;
using std::vector;
using std::string;

void print_help_text()
{
  std::cout << "Purpose: Update an LR B-spline surface with a new point cloud. \n";
  std::cout << "Mandatory parameters: input surface (.g2), input point cloud (.txt, .xyz or .g2), output surface (.g2), tolerance, number of iterations. \n";
  std::cout << "The output surface sill be non-trimmed even if the input surface is trimmed \n";
  std::cout << "Note that the number of iterations must take into consideration ";
  std::cout << "that the surface already has a certain size.\n";
  std::cout << "For a large surface, the number of iterations should be low (1-3) \n";
  std::cout << "0 iterations will only provide information about the distance, not update the surface.\n";
  std::cout << "1 iteration will update the surface within the current spline space.\n";
  std::cout << "An adaptive approximation procedure is applied which for the";
  std::cout << " specified number of iterations: \n";
  std::cout << " - Approximates the points with a surface in the current spline space \n";
  std::cout << " - Computes the approximation accuracy \n";
  std::cout << " - Refines the surfaces in areas where the tolerance is not met \n";
  std::cout << "The approximation is completed when the tolerance is met or the";
  std::cout << " specified number of iterations is exceeded.\n";
  std::cout << "For a 1D surface, the points are expected to be given as x, y, z and be parameterized on x and y. \n";
  std::cout << "For a 3D surface, the points must be given as u, v, x, y, z. ";
  std::cout << "In that case, only a .txt file is accepted. \n";
  std::cout << "Optional input parameters: \n";
  std::cout << "-dist <filename (.txt)> : Write distance field to file (x, y, z, distance) \n";
  std::cout << "-info <filename> : Write accuracy information to file \n";
  std::cout << "-h or --help : Write this text\n";
}

int compareext(const char *str1, char str2[][8], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}

int main(int argc, char *argv[])
{
  char* input_type = 0;    // Type of point file
  char *pointfile = 0;     // Input point file
  char *surffile = 0;       // Surface output file
  char *outfile = 0;       // Surface output file
  char *infofile = 0;      // Accuracy information output file
  int del = 3;             // Number of entries for each point
  double AEPSGE = 0.5;     // Requested accuracy
  int max_iter = 6;        // Maximum number of iteations in adaptive alogrithm
  char *field_out = 0;     // Distance field output file

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
      else if (arg == "-dist")
	{
	  if (ki == argc-1)
	    {
	      std::cout << "ERROR: Missing input" << std::endl;
	      print_help_text();
	      return 1;
	    }
	  field_out = argv[ki+1];
	  par_read[ki-1] = par_read[ki] = true;
	  nmb_par -= 2;
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

  // Read remaining parameters
  if (nmb_par != 5)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }


  for (ki=1; ki<argc; ++ki)
    {
      if (par_read[ki-1])
	continue;
      if (nmb_par == 5)
	{
	  surffile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 4)
	{
	  pointfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 3)
	{
	  outfile = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  AEPSGE = atof(argv[ki]);
	  nmb_par--;
	}
      else
	{
	  max_iter = atoi(argv[ki]);
	}
    }
      
   time_t tread = time(0);
  tm *ltm = localtime(&tread);
  std::cout << "Reading: Reading input: " << tread;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  // Read input surface
  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  std::ifstream is3(surffile);
  ObjectHeader header;
  shared_ptr<GeomObject> geom_obj;
  try {
  header.read(is3);
  geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
  geom_obj->read(is3);
  }
  catch(...)
    {
      std::cout << "WARNING: No surface found" << std::endl;
      exit(0);
    }
  
  shared_ptr<ParamSurface> sf0 = 
    dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf0.get())
    {
      std::cout << "ERROR: Failed reading surface: " << std::endl;
      return 1;
    }
  shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface,
							  ParamSurface>(sf0);
  if (bd_sf.get())
    sf0 = bd_sf->underlyingSurface();
  shared_ptr<LRSplineSurface> sf1 = dynamic_pointer_cast<LRSplineSurface,
							  ParamSurface>(sf0);
  if (!sf1.get())
   {
     std::cout << "ERROR: Input surface is not of type LR B-spline:" << std::endl;
     return 1;
   }
 

  if (!sf1.get() || (sf1->dimension() != 1 && sf1->dimension() != 3))
    {
      std::cout << "ERROR: Not a valid surface" << std::endl;
      return -1; 
    }

  std::ofstream ofsf("sfout.g2");
  sf1->writeStandardHeader(ofsf);
  sf1->write(ofsf);

  if (sf1->dimension() == 3)
    del = 5;

  double smoothwg = 1.0e-10; 
  int mba = 1;      // Use LR-MBA approximation
  int tomba = 0; // Already MBA

  // Read point cloud
  vector<double> data;
  vector<double> extent(2*del);   // Limits for points in all coordinates
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
    FileUtils::readTxtPointFile(pointsin, del, data, nmb_pts, extent);


  tread = time(0);
  ltm = localtime(&tread);
  std::cout << "End of reading: Point cloud read: " << tread;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  time_t tprocess = time(0);
  ltm = localtime(&tprocess);
  std::cout << "Processing:  Surface update starts: " << tprocess;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

 // Move point cloud to origo
  double umin = sf1->paramMin(XFIXED);
  double umax = sf1->paramMax(XFIXED);
  double vmin = sf1->paramMin(YFIXED);
  double vmax = sf1->paramMax(YFIXED);
  Point mid;
  if (sf1->dimension() == 1)
    mid.setValue(0.5*(umin+umax), 0.5*(vmin+vmax), 0.0);
  else
    {
      mid = Point(0.5*(extent[2*(del-3) + 2*(del-3)+1]),
		  0.5*(extent[2*(del-2) + 2*(del-2)+1]), 0.0);
    }
  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      {
	data[del*ki+kj] -= mid[kj-del+3];
      }

  // Move surface to origo
  try {
    if (sf1->dimension() == 1)
      {
	sf1->setParameterDomain(umin - mid[0], umax - mid[0],
				vmin - mid[1], vmax - mid[1]);
      }
    else
      sf1->translate(-mid);
  }
  catch (...)
    {
      std::cout << "ERROR: Translation to origo failed" << std::endl;
      return 1;
    }
      

#ifdef DEBUG
  // Write translated surface and points to g2 format
  vector<double> data2;
  data2.reserve(nmb_pts*3);
  for (ki=0, kj=0; ki<nmb_pts; ++ki, kj+=del)
    data2.insert(data2.end(), data.begin()+kj, data.begin()+kj+3);
  PointCloud3D cloud(data2.begin(), nmb_pts);

  std::ofstream of1("translated_sf.g2");
  std::ofstream of2("translated_points.g2");
  cloud.writeStandardHeader(of2);
  cloud.write(of2);
#endif
  
  
  bool repar = true;
  LRSurfApprox approx(sf1, data, AEPSGE, true, repar, true);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(true);
  if (mba)
    approx.setUseMBA(true);
  else
    {
      approx.setSwitchToMBA(tomba);
    }
  approx.setVerbose(true);

  if (del == 3)
    {
      double zrange = extent[5] - extent[4];
      approx.addLowerConstraint(extent[4] - 0.1*(zrange));
      approx.addUpperConstraint(extent[5] + 0.1*(zrange));
    }

  double maxdist, avdist, avdist_total; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<LRSplineSurface> surf = approx.getApproxSurf(maxdist, avdist_total,
							  avdist, nmb_out_eps, 
							  max_iter);

 
  tprocess = time(0);
  ltm = localtime(&tprocess);
  std::cout << "End of processing: Approximation completed: " << tprocess;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  time_t twrite = time(0);
  ltm = localtime(&twrite);
  std::cout << "Writing: Write updated surface: " << twrite;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  std::ofstream sfout(outfile);     // Surface output stream

  if (infofile)
    {
      std::ofstream infoout(infofile);   // Accuracy information output stream

      infoout << "Total number of points: " << nmb_pts << std::endl;
      infoout << "Number of elements: " << surf->numElements() << std::endl;
      infoout << "Maximum distance: " << maxdist << std::endl;
      infoout << "Average distance: " << avdist_total << std::endl;
      infoout << "Average distance for points outside of the tolerance: " << avdist << std::endl;
      infoout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl;
    } 
  else
    {
      std::cout << "INFO: Total number of points: " << nmb_pts << std::endl;
      std::cout << "INFO: Number of elements: " << surf->numElements() << std::endl;
      std::cout << "INFO: Maximum distance: " << maxdist << std::endl;
      std::cout << "INFO: Average distance: " << avdist_total << std::endl;
      std::cout << "INFO: Average distance for points outside of the tolerance: " << avdist << std::endl;
      std::cout << "INFO: Number of points outside the tolerance: " << nmb_out_eps << std::endl;
    } 


  if (surf.get())
    {
#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      if (surf->dimension() == 3)
	{
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
	}
      else
	{
	  shared_ptr<LRSplineSurface> surf2(surf->clone());
	  surf2->to3D();
	  surf2->writeStandardHeader(of1);
	  surf2->write(of1);
	}
#endif
	  
      // Translate
      try {
	if (surf->dimension() == 3)
	  {
	    surf->translate(mid);
	  }
	else
	  {
	    // Update parameter domain
	    double umin = surf->paramMin(XFIXED);
	    double umax = surf->paramMax(XFIXED);
	    double vmin = surf->paramMin(YFIXED);
	    double vmax = surf->paramMax(YFIXED);
	  
	    surf->setParameterDomain(umin + mid[0], umax + mid[0],
				     vmin + mid[1], vmax + mid[1]);
	  }
      }
      catch (...)
	{
	  std::cout << "ERROR: Translation to origo failed" << std::endl;
	  return 1;
	}
      
      try {
	surf->writeStandardHeader(sfout);
	surf->write(sfout);
      }
      catch (...)
	{
	  std::cout << "ERROR: Failed writing surface to file" << std::endl;
	  return 1;
	}

      twrite = time(0);
      ltm = localtime(&twrite);
      std::cout << "End of writing: Finished writing surface to file and accuracy info: " << twrite;
      std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
      std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
      std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

       if (field_out)
	{
	  twrite = time(0);
	  ltm = localtime(&twrite);
	  std::cout << "Writing: Write distance field: " << twrite;
	  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
	  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
	  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

	  // Fetch data points with distance information
	  vector<double> pnts_dist;
	  pnts_dist.reserve(4*nmb_pts);
	  try {
	    LRSplineSurface::ElementMap::const_iterator elem = surf->elementsBegin();
	    LRSplineSurface::ElementMap::const_iterator last = surf->elementsEnd();
	    for (; elem != last; ++elem)
	      {
		if (!elem->second->hasDataPoints())
		  continue;
		vector<double>& points = elem->second->getDataPoints();
		pnts_dist.insert(pnts_dist.end(), points.begin(), points.end());
	      }
	  }
	  catch (...)
	    {
	      std::cout << "ERROR: Extraction of trimming information failed" << std::endl;
	      return 1;
	    }

	  // Translate to initial domain
	  for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
	    {
	      pnts_dist[kj] += mid[0];
	      pnts_dist[kj+1] += mid[1];
	    }

	  // Write to file
	  // Find file extension
	  char *loc2 = strchr(field_out, '.');
	  char *last2 = 0;
	  while (loc2 != NULL)
	    {
	      last2 = loc2;
	      loc2 = strchr(loc2+1, '.');
	    }
	  if (last2 == NULL)
	    {
	      std::cout << "ERROR: Missing file extension of output file" << std::endl;
	      return 1;
	    }
	  char *out_type = last2+1;
	  int dist_type;
	  try {
	    dist_type = compareext(out_type, keys, 9);
	  }
	  catch (...)
	    {
	      std::cout << "ERROR: Invalide file extension of output file" << std::endl;
	      return 1;
	    }

	  if (dist_type <= 1)
	    {
	      std::ofstream field_info(field_out);
	      (void)field_info.precision(15);
	      for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
		{
		  for (ki=0; ki<4; ++ki)
		    field_info << pnts_dist[kj+ki] << " ";
		  field_info << std::endl;
		}
	    }
	  else
	    {
	      std::cout << "INFO: Distance output file not recognized, no distance output created" << std::endl;
	    }
	  twrite = time(0);
	  ltm = localtime(&twrite);
	  std::cout << "End of writing: Finished writing distance field to file: " << twrite;
	  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
	  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
	  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;
	}

#ifdef DEBUG
      if (surf->dimension() == 1)
	{
	  std::ofstream of2("surf_3D.g2");
	  surf->to3D();
	  surf->writeStandardHeader(of2);
	  surf->write(of2);
	}
#endif
    }
  return 0;
}

