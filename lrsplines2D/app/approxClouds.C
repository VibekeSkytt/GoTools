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
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRSurfUtils.h"
#include "GoTools/utils/timeutils.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/stat.h>
#include <boost/timer.hpp>
#include <time.h>
#include <dirent.h>

//#define DEBUG
//#define DEBUG_EL
//#define DEBUG2


using namespace Go;
using std::vector;
using std::string;

void print_help_text()
{
  std::cout << "Purpose: Approximate a number of point clouds by LR B-spline surfaces defined on the same mesh. \n";
  std::cout << "Mandatory parameters: input folder (point clouds .g2), output folder (surfaces .g2), tolerance, number of iterations. \n";
  std::cout << "An adaptive approximation procedure is applied which for the";
  std::cout << " specified number of iterations: \n";
  std::cout << " - Approximates each point set with a surface in the current spline space \n";
  std::cout << " - Computes the approximation accuracy \n";
  std::cout << " - Refines the surfaces in areas where the tolerance is not met \n";
  std::cout << "The approximation is completed when the tolerance is met or the";
  std::cout << " specified number of iterations is exceeded.\n";
  std::cout << "The surfaces are defined on a union of the mesh corresponding to the individual surfaces, ";
  std::cout << "and approximated on the refined mesh. \n";
  std::cout << "The number of iterations is recommended to lie in the interval [4:7]. \n";
  std::cout << "The points are expected to be given as x, y, z and be parameterized on x and y, but 3D parameterized points are accepted.\n";
  std::cout << "Then the points are given as u, v, x, y, z. \n";
  std::cout << "Optional input parameters: \n";
  std::cout << "-dist <filename (.txt)> : Write distance field to file (x, y, z, distance) \n";
  std::cout << "-info <filename> : Write accuracy information to file \n";
  std::cout << "-smooth <weight> : Overrule default smoothing weight (1.0e-10) in least squares approximation \n";
  std::cout << "-mba <0/1/n/-1>: 0 = use only least squares approximation \n";
  std::cout << "                 1 = use only multilevel B-spline approximation (MBA) \n";
  std::cout << "                 n = start with least squares, turn to MBA after n iterations \n";
  std::cout << "                -1 = initiate computation using MBA \n";
  std::cout << "Default setting is start with least squares, turn to MBA for the last iterations \n";
  std::cout << "-degree <polynomial degree> : 2 or 3 recommended \n";
  std::cout << "-nmb_coef <initial value> : Initial number of coefficients in each parameter direction \n";
  std::cout << "-distributecf <0/1> : Modify initial number of coefficients according to relative size of parameter domain \n";
  std::cout << "-minsize <size> : Minimum element size, all directions \n";
  std::cout << "-reltol <0/1>: Apply relative tolerance flag. Default false \n";
  std::cout << "-tolfac1: Factor for modification of tolerance, positive heights. Default 0.0 \n";
  std::cout << "-tolfac2: Factor for modification of tolerance, negative heights. Default 0.0 \n";
  std::cout << "-verbose <0/1>: Write accuracy information at each iteration level. Default = 0 \n";
  std::cout << "-refstrat: Print parameter information to define refinement strategy. \n";
  std::cout << "-h or --help : Write this text\n";
}

void print_refine_info()
{
  std::cout << "Define refinement strategy from command line with the following parameters: \n";
  std::cout << "-refcat <code>: F = full span, Ml/Mu/Mc = minimum span, S = structured mesh, \n";
  std::cout << "R = restricted mesh, RL = restricted mesh with element extension. Default = F. \n";
  std::cout << "-alterdir <code>: B = refine in both parameter directions, A = refine in alternating directions \n";
  std::cout << "-threshold <code>: no = none, td = with respect to distance, tn = with respect to number, \n";
  std::cout << "tdk = distance and number for category R and S. Default = no \n";
  std::cout << "-swapcat <size>: Swap strategy when approximation efficience drops below size. Default = -100 \n";
  std::cout << "-refcat2 <code>: After swap, as -refcat \n";
  std::cout << "-threshold2 <code>: After swap, as -threshold \n";
}



int fetchIntParameter(int argc, char *argv[], int ki, int& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = atoi(argv[ki+1]);
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int fetchDoubleParameter(int argc, char *argv[], int ki, double& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = atof(argv[ki+1]);
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

int fetchCharParameter(int argc, char *argv[], int ki, char*& parameter, 
		      int& nmb_par, vector<bool>& par_read)
{
  if (ki == argc-1)
    {
      std::cout << "ERROR: Missing input" << std::endl;
      print_help_text();
      return -1;
    }
  parameter = argv[ki+1];
  par_read[ki-1] = par_read[ki] = true;
  nmb_par -= 2;
  return 0;
}

void combinePathFile(const char *path, string& file, string& file2)
{
  char outfile[160];
  strcpy(outfile, path);
  char tmp1[2];
  sprintf(tmp1,"/");
  strncat(outfile, tmp1, 1);
  strcat(outfile, file.c_str());
  file2 = std::string(outfile);
}

int main(int argc, char *argv[])
{
  char* input_type = 0;    // Type of point file
  char *folderin = 0;     // Input point file
  char *folderout = 0;       // Surface output file
  char *infofile = 0;      // Accuracy information output file
  int del = 3;             // Number of entries for each point
  double AEPSGE = 0.5;     // Requested accuracy
  int max_iter = 6;        // Maximum number of iterations in adaptive alogrithm
  char *field_out = 0;     // Distance field output file
  double smoothwg = 1.0e-9; 
  int initmba = 0; //1;  // Initiate surface using the mba method
  int mba = 0;      // Use least squares approximation
  int tomba = std::min(5, max_iter-1);    // Turn to the mba method at 
  // iteration level 5 or in the last iteration
  int degree = 2;
  int outlierflag = 0;
  int reltol = 0;
  double tolfac1 = 0.0, tolfac2 = 0.0;
  double minsize = -1.0;
  int verbose = 0;
  
  int initncoef = 10;
  int distribute_ncoef = 0;

  int refcat1=1, refcat2=0, threshold1=-1, threshold2=-1, alter=0;
  double swap = -100.0;

  int ka;
  vector<bool> par_read(argc-1, false);

  // Read optional parameters
  int nmb_par = argc-1;
  for (ka=1; ka<argc; ++ka)
    {
      string arg(argv[ka]);
      if (arg == "-h" || arg == "--help")
	{
	  print_help_text();
	  exit(0);
	}
      else if (arg == "-refstrat")
	{
	  print_refine_info();
	  exit(0);
	}
      else if (arg == "-dist")
	{
	  int stat = fetchCharParameter(argc, argv, ka, field_out, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-info")
	{
	  int stat = fetchCharParameter(argc, argv, ka, infofile, 
					nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-smooth")
	{
	  int stat = fetchDoubleParameter(argc, argv, ka, smoothwg, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-mba")
	{
	  int mm;
	  int stat = fetchIntParameter(argc, argv, ka, mm, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  if (mm == 0)
	    tomba = 100;
	  else if (mm == 1)
	    mba = 1;
	  else if (mm < 0)
	    initmba = 1;
	  else
	    tomba = mm;
	}
      else if (arg == "-degree")
	{
	  int stat = fetchIntParameter(argc, argv, ka, degree, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-nmb_coef")
	{
	  int stat = fetchIntParameter(argc, argv, ka, initncoef, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-distributecf")
	{
	  int stat = fetchIntParameter(argc, argv, ka, distribute_ncoef, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
       else if (arg == "-outlier")
	{
	  int stat = fetchIntParameter(argc, argv, ka, outlierflag, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-minsize")
	{
	  int stat = fetchDoubleParameter(argc, argv, ka, minsize, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-reltol")
	{
	  int stat = fetchIntParameter(argc, argv, ka, reltol, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-tolfac1")
	{
	  int stat = fetchDoubleParameter(argc, argv, ka, tolfac1, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-tolfac2")
	{
	  int stat = fetchDoubleParameter(argc, argv, ka, tolfac2, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-verbose")
	{
	  int stat = fetchIntParameter(argc, argv, ka, verbose, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
      else if (arg == "-refcat")
	{
	  char* cat0;
	  int stat = fetchCharParameter(argc, argv, ka, cat0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string cat(cat0);
	  if (cat == "F")
	    refcat1 = 1;
	  else if (cat == "Ml")
	    refcat1 = 2;
	  else if (cat == "Mu")
	    refcat1 = 3;
	  else if (cat == "Mc")
	    refcat1 = 4;
	  else if (cat == "S")
	    refcat1 = 5;
	  else if (cat == "R")
	    refcat1 = 6;
	  else if (cat == "RL")
	    refcat1 = 7;
	  if (refcat2 < 1 || refcat2 > 7)
	    refcat2 = refcat1;
	}
      else if (arg == "-refcat2")
	{
	  char* cat0;
	  int stat = fetchCharParameter(argc, argv, ka, cat0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string cat(cat0);
	  if (cat == "F")
	    refcat2 = 1;
	  else if (cat == "Ml")
	    refcat2 = 2;
	  else if (cat == "Mu")
	    refcat2 = 3;
	  else if (cat == "Mc")
	    refcat2 = 4;
	  else if (cat == "S")
	    refcat2 = 5;
	  else if (cat == "R")
	    refcat2 = 6;
	  else if (cat == "RL")
	    refcat2 = 7;
	}
      else if (arg == "-threshold")
	{
	  char* thresh0;
	  int stat = fetchCharParameter(argc, argv, ka, thresh0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string thresh(thresh0);
	  if (thresh == "no")
	    threshold1 = 0;
	  else if (thresh == "td")
	    threshold1 = 1;
	  else if (thresh == "tn")
	    threshold1 = (refcat1 >= 6) ? 3 : 2;
	  else if (thresh == "tdk")
	    threshold1 = 4;
	  if (threshold2 < 0)
	    threshold2 = threshold1;
	}
      else if (arg == "-threshold2")
	{
	  char* thresh0;
	  int stat = fetchCharParameter(argc, argv, ka, thresh0, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;

	  string thresh(thresh0);
	  if (thresh == "no")
	    threshold2 = 0;
	  else if (thresh == "td")
	    threshold2 = 1;
	  else if (thresh == "tn")
	    threshold2 = (refcat1 >= 6) ? 3 : 2;
	  else if (thresh == "tdk")
	    threshold2 = 4;
	}
      else if (arg == "-alterdir")
	{
	  char *altdir;
	  int stat = fetchCharParameter(argc, argv, ka, altdir, 
				       nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	  alter = (altdir[0] == 'A') ? 1 : 0;
	}
      else if (arg == "-swap")
	{
	  int stat = fetchDoubleParameter(argc, argv, ka, swap, 
					  nmb_par, par_read);
	  if (stat < 0)
	    return 1;
	}
	
    }

  // Read remaining parameters
  if (nmb_par != 4)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

  for (ka=1; ka<argc; ++ka)
    {
      if (par_read[ka-1])
	continue;
      if (nmb_par == 4)
	{
	  folderin = argv[ka];
	  nmb_par--;
	}
      else if (nmb_par == 3)
	{
	  folderout = argv[ka];
	  nmb_par--;
	}
      else if (nmb_par == 2)
	{
	  AEPSGE = atof(argv[ka]);
	  nmb_par--;
	}
      else
	{
	  max_iter = atoi(argv[ka]);
	}
    }


  // Read point files
  DIR *dir;
  struct dirent *ent;
  vector<string> infile;
  vector<string> pointfile;
  string ext("_sf.g2");
  string delimiter("/");
  if ((dir = opendir (folderin)) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      string data_file(ent->d_name);
      if (data_file == "." || data_file == "..")
	continue;
      pointfile.push_back(data_file);
      string data_file2;
      combinePathFile(folderin, data_file, data_file2);
      infile.push_back(data_file2);
    }
  } else {
    /* could not open directory */
    std::cout << "ERROR: Could not open file directory" << std::endl;
    return 1;
  }

  std::ofstream infoout(infofile);
  for (size_t kj=0; kj<infile.size(); ++kj)
    {
      // read point cloud      
      std::ifstream pointsin(infile[kj]);
      vector<double> extent;   // Limits for points in all coordinates
      int nmb_pts = 0;
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
      vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);
      for (ka=0; ka<3; ++ka)
	{
	  extent.push_back(low[ka]);
	  extent.push_back(high[ka]);
	}
      std::cout << "Domain: [" << extent[0] << ", " << extent[1] << "] x [";
      std::cout << extent[2] << ", " << extent[3] << "]" << std::endl;
      std::cout << "Range: " << extent[4] << " - " << extent[5] << std::endl;


      // Move point cloud to origo
      Point mid(0.5*(extent[0] + extent[1]), 0.5*(extent[2] + extent[3]), 0.0);
      for (int ka=0; ka<nmb_pts; ++ka)
	for (int kb=del-3; kb<del-1; ++kb)
	  {
	    data[del*ka+kb] -= mid[kb-del+3];
	  }
      
      int order = degree + 1; 
      int nmb_coef = std::max(order, initncoef); //std::max(order, 14);
      int nm = nmb_coef*nmb_coef;
      double dom = (extent[1]-extent[0])*(extent[3]-extent[2]);
      double c1 = std::pow((double)nm/dom, 1.0/2.0);
      int nc[2];
      for (int kb=0; kb<2; ++kb)
	{
	  double tmp = c1*(extent[2*kb+1]-extent[2*kb]);
	  nc[kb] = (int)tmp;
	  //if (tmp - (double)nc[kj] > 0.5)
	  ++nc[kb];
	  nc[kb] = std::max(nc[kb], order);
	}
      double mba_coef = 0.0;
      if (initmba)
	mba_coef = 0.5*(extent[2*(del-1)] + extent[2*(del-1)+1]);
      vector<shared_ptr<LRSplineSurface> > sfs;
      for (size_t ki=0; ki<data.size(); ++ki)
	{
	  shared_ptr<LRSurfApprox> approx;
	  if (distribute_ncoef)
	    approx = shared_ptr<LRSurfApprox>(new LRSurfApprox(nc[0], order, nc[1], order, data, del-2, 
							       AEPSGE, initmba ? true : false, 
							       mba_coef, true, true));
	  else
	    approx = shared_ptr<LRSurfApprox>(new LRSurfApprox(nmb_coef, order, nmb_coef, order, data, del-2, 
							       AEPSGE, initmba ? true : false, 
							       mba_coef, true, true));
	  approx->setSmoothingWeight(smoothwg);
	  approx->setSmoothBoundary(true);
	  if (mba)
	    approx->setUseMBA(true);
	  else
	    {
	      // if (initmba)
	      // 	approx->setInitMBA(initmba, 0.5*(low[2]+high[2]));
	      //if (del == 3)
	      approx->setSwitchToMBA(tomba);
	      approx->setMakeGhostPoints(false /*true*/);
	    }

	  if (minsize > 0.0)
	    approx->setMinimumElementSize(minsize, minsize);
      
	  approx->setVerbose(verbose);


	  // Refine strategy
	  approx->setRefinementStrategy(refcat1, alter, threshold1, swap, refcat2, threshold2);
      
	  if (del == 3)
	    {
	      double zrange = extent[5] - extent[4];
	      double zfac = std::max(AEPSGE, 0.005*zrange);
	      approx->addLowerConstraint(extent[4] - zfac);
	      approx->addUpperConstraint(extent[5] + zfac);
	      approx->setLocalConstraint(zfac);
	    }

	  double maxdist, avdist, avdist_total; // will be set below
	  int nmb_out_eps;        // will be set below
	  double maxout, avout;
	  shared_ptr<LRSplineSurface> surf;
	  try {
	    surf = approx->getApproxSurf(maxdist, avdist_total,
					 avdist, nmb_out_eps, 
					 max_iter);
	  }
	  catch (...)
	    {
	      std::cout << "ERROR: Surface approximation failed" << std::endl;
	      return 1;
	    }


	  if (infofile)
	    {
	      infoout << "Data set number " << kj+1 << std::endl;
	      infoout << pointfile[kj] << std::endl;
	      infoout << "Total number of points: " << nmb_pts << std::endl;
	      infoout << "Number of elements: " << surf->numElements() << std::endl;
	      infoout << "Number of coefficients: " << surf->numBasisFunctions() << std::endl;
	      infoout << "Maximum distance: " << maxdist << std::endl;
	      infoout << "Average distance: " << avdist_total << std::endl;
	      infoout << "Average distance for points outside of the tolerance: " << avdist << std::endl;
	      infoout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl << std::endl;
	    }
	  else
	    {
	      std::cout << "Data set number " << kj+1 << std::endl;
	      std::cout << pointfile[kj] << std::endl;
	      std::cout << "INFO: Total number of points: " << nmb_pts << std::endl;
	      std::cout << "INFO: Number of elements: " << surf->numElements() << std::endl;
	      std::cout << "INFO: Number of coefficients: " << surf->numBasisFunctions() << std::endl;
	      std::cout << "INFO: Maximum distance: " << maxdist << std::endl;
	      std::cout << "INFO: Average distance: " << avdist_total << std::endl;
	      std::cout << "INFO: Average distance for points outside of the tolerance: " << avdist << std::endl;
	      std::cout << "INFO: Number of points outside the tolerance: " << nmb_out_eps << std::endl << std::endl;
	    }



	  // Translate
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

	  // Construct output file name
	  string name;
	  FileUtils::extractPathName(pointfile[kj], name);
	  string base, base0, fileout;
	  FileUtils::extendName(folderout, delimiter.c_str(), base0);
	  FileUtils::extendName(base0.c_str(), name.c_str(), base);
	  FileUtils::extendName(base.c_str(), ext.c_str(), fileout);
	  std::cout << fileout << std::endl;
	  std::ofstream sfout(fileout);
	  surf->writeStandardHeader(sfout);
	  surf->write(sfout);

	}
    }
  return 0;
}

