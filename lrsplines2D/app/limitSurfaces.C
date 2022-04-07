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
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRApproxApp.h"
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
  std::cout << "Purpose: Compute limit surfaces from an LR B-spline surface with respect to a corresponding point cloud. \n";
  std::cout << "Mandatory parameters: input surface (.g2), input point cloud (.txt, .xyz or .g2), number of MBA iterations, output surface1 (.g2), output surface2 (.g2) \n";
  std::cout << "-h or --help : Write this text\n";
}

int compareext(const char *str1, char str2[][8], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}

/// Compute surfaces bounding a point cloud on both sides of a LRSplineSurface approximating
/// the point cloud. The points are expected to be defined by (x,y,z) and be parameterized
/// in x and y in the surface approximation. The surface is, thus, a function f(x,y)
/// approximating the z-coordinates of the points. The bounding surfaces are represented
/// in the same spline space as the input surface.

int main(int argc, char *argv[])
{
  char* input_type = 0;    // Type of point file
  char *pointfile = 0;     // Input point file
  char *surffile = 0;       // Surface output file
  char *outfile1 = 0;       // Surface output file
  char *outfile2 = 0;       // Surface output file
  int nmb_MBA;

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
	  nmb_MBA = atoi(argv[ki]);
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  outfile1 = argv[ki];
	  nmb_par--;
	}
      else if (nmb_par == 1)
	{
	  outfile2 = argv[ki];
	  nmb_par--;
	}
    }
      
   time_t tread = time(0);
  tm *ltm = localtime(&tread);
  std::cout << "Reading: Reading input: " << tread;
  std::cout << " | " << 1900+ltm->tm_year << "-" << 1+ltm->tm_mon << ":";
  std::cout << ltm->tm_mday << " " << ltm->tm_hour << ":";
  std::cout << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read input surface
  std::ifstream insurf(surffile);
  shared_ptr<GeomObject> geom_obj;
  ObjectHeader header;
  try {
    header.read(insurf);
  }
  catch (...)
    {
      std::cout << "ERROR: Input object not recognized. Exiting" << std::endl;
      return 1;
    }
  try {
    geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
    geom_obj->read(insurf);
  }
  catch (...)
    {
      std::cout << "ERROR: Input surface could not be read. Exiting" << std::endl;
      return 1;
    }
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cout << "ERROR: Input file contains no surface" << std::endl;
      return 1;
    }
 

  if (!sf.get() || sf->dimension() != 1)
    {
      std::cout << "ERROR: Not a valid surface" << std::endl;
      return -1; 
    }
  int del = 2 + sf->dimension();

  // Find file extension
  char *loc;
  char *last = 0;
  loc = strchr(pointfile, '.');
  while (loc != NULL)
    {
      last = loc;
      loc = strchr(loc+1, '.');
    }
  if (last == NULL)
    {
      std::cout << "ERROR: Missing file extension of input file" << std::endl;
      return 1;
    }
  input_type = last+1;

  // Check type of input file
  char keys[5][8] = {"txt", "TXT", "xyz", "XYZ", "g2"};
  int type;
  try {
    type = compareext(input_type, keys, 5);
  }
  catch (...)
    {
      std::cout << "ERROR: Invalide file extension of input file" << std::endl;
      return 1;
    }
  if (type < 0)
    {
      std::cout << "ERROR: Not a valid point file type" << std::endl;
      return 1;
    }
  if (type > 3 && del == 5)
    {
      std::cout << "ERROR: Parameterized points expected, pointfile type not allowed" << std::endl;
      return 1;
    }
  
  // Read points
  int nmb_pts = 0;
  vector<double> data;
  vector<double> extent(2*del);   // Limits for points in all coordinates
  std::ifstream pointsin(pointfile);
  if (type <= 3)
    {
      // Read xyz (separated by , or space)
      try {
      FileUtils::readTxtPointFile(pointsin, del, data, nmb_pts, extent);
      }
      catch (...)
	{
	  std::cout << "ERROR: Not a valid point file" << std::endl;
	  return 1;
	}
    }
  else 
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
	  return 1;
	}
      BoundingBox box = points.boundingBox();
      Point low = box.low();
      Point high = box.high();
      nmb_pts = points.numPoints();
      data.insert(data.end(), points.rawData(), points.rawData()+3*nmb_pts);
      for (ki=0; ki<3; ++ki)
	{
	  extent[2*ki] = low[ki];
	  extent[2*ki+1] = high[ki];
	}
    }

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
  shared_ptr<BoundedSurface> bd_sf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
  RectDomain dom = (bd_sf.get()) ?
    bd_sf->underlyingSurface()->containingDomain() :
    sf->containingDomain();
  double umin = dom.umin();
  double umax = dom.umax();
  double vmin = dom.vmin();
  double vmax = dom.vmax();
  Point mid;
  if (sf->dimension() == 1)
    mid.setValue(0.5*(umin+umax), 0.5*(vmin+vmax), 0.0);


  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      {
	data[del*ki+kj] -= mid[kj-del+3];
      }

  // Move surface to origo
  try {
    if (sf->dimension() == 1)
      {
	sf->setParameterDomain(umin - mid[0], umax - mid[0],
				vmin - mid[1], vmax - mid[1]);
      }
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
  

  shared_ptr<ParamSurface> limsf1;
  shared_ptr<ParamSurface> limsf2;

  LRApproxApp::limitingSurfs(data, sf, nmb_MBA, limsf1, limsf2);
 
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

  std::ofstream sfout1(outfile1);     // Surface output stream
  std::ofstream sfout2(outfile2);

  // Translate
  try {
    // Update parameter domain
    limsf1->setParameterDomain(umin, umax, vmin, vmax);
    limsf2->setParameterDomain(umin, umax, vmin, vmax);
  }
  catch (...)
    {
      std::cout << "ERROR: Translation to origo failed" << std::endl;
      return 1;
    }
  
  try {
    limsf1->writeStandardHeader(sfout1);
    limsf1->write(sfout1);
    limsf2->writeStandardHeader(sfout2);
    limsf2->write(sfout2);
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

  return 0;
}

