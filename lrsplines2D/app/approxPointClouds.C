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

//#define DEBUG
//#define DEBUG_EL
//#define DEBUG2


using namespace Go;
using std::vector;
using std::string;

void print_help_text()
{
  std::cout << "Purpose: Approximate a number of point clouds by LR B-spline surfaces defined on the same mesh. \n";
  std::cout << "Mandatory parameters: input point clouds (.g2), output surfaces (.g2), tolerance, number of iterations. \n";
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
  std::cout << "-tolfile: File specifying domains with specific tolerances, global tolerance apply outside domains. PointCloud2LR -tolfile for file format \n";
  std::cout << "-toldoc: Documentation on file format for tolerance domains. \n";
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


void print_tol_file_format()
{
  std::cout << "File specifying domains/boxes with different tolerances. \n";
  std::cout << "If not otherwise stated, the global tolerance will applies outside the boxes \n";
  std::cout << "Numbers are given as floats and separated by space \n";
  std::cout << "Line1: Number of boxes (integer), whether or not the tolerance is specified outside the boxes (0/1) \n";
  std::cout << "Following lines: \n";
  std::cout << "xmin ymin xmax ymax tolerance \n";
  std::cout << "Last line (if given): Tolerance for remaining area, overrules tolerance given in parameter line. \n";
  std::cout << "Ensure non-overlapping boxes. No test applied. \n";
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

int main(int argc, char *argv[])
{
  char* input_type = 0;    // Type of point file
  char *pointfile = 0;     // Input point file
  char *surffile = 0;       // Surface output file
  char *infofile = 0;      // Accuracy information output file
  char *tolfile = 0;       // File specifying varying tolerances
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
      else if (arg == "-toldoc")
	{
	  print_tol_file_format();
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
      else if (arg == "-tolfile")
	{
	  int stat = fetchCharParameter(argc, argv, ka, tolfile, 
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
	  pointfile = argv[ka];
	  nmb_par--;
	}
      else if (nmb_par == 3)
	{
	  surffile = argv[ka];
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

  // Read point clouds
  vector<vector<double> > data;
  vector<double> extent;   // Limits for points in all coordinates

  std::ifstream pointsin(pointfile);
  int kb = 0;
  while (!pointsin.eof())
    {
      // read point cloud      
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
      Utils::eatwhite(pointsin);
      
      BoundingBox box = points.boundingBox();
      Point low = box.low();
      Point high = box.high();
      nmb_pts = points.numPoints();
      vector<double> data2(points.rawData(), points.rawData()+3*nmb_pts);
      data.push_back(data2);
      for (ka=0; ka<3; ++ka)
	{
	  extent.push_back(low[ka]);
	  extent.push_back(high[ka]);
	}
      std::cout << "Domain " << kb+1 << ": [" << extent[6*kb] << ", " << extent[6*kb+1] << "] x [";
      std::cout << extent[6*kb+2] << ", " << extent[6*kb+3] << "]" << std::endl;
      std::cout << "Range: " << extent[6*kb+4] << " - " << extent[6*kb+5] << std::endl;
      ++kb;
    }


  // Modify data extent if necessary
  double domain[6];
  domain[0] = domain[2] = domain[4] = std::numeric_limits<double>::max();
  domain[1] = domain[3] = domain[5] = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<extent.size(); ki+=6)
    for (int ka=0; ka<3; ++ka)
	{
	  domain[2*ka] = std::min(domain[2*ka], extent[ki+2*ka]);
	  domain[2*ka+1] = std::max(domain[2*ka+1], extent[ki+2*ka+1]);
	}
  
  bool use_stdd = false;
  vector<LRSurfApprox::TolBox> tolerances;
  if (tolfile != 0)
    {
      std::ifstream tolin(tolfile);
      int nmb_box, last;
      double umin, umax, vmin, vmax, tol;
      tolin >> nmb_box >> last;
      tolerances.resize(nmb_box);
      for (int ka=0; ka<nmb_box; ++ka)
	{
	  tolin >> umin >> umax >> vmin >> vmax >> tol;
	  tolerances[ka].setVal(std::max(umin,extent[0]), std::min(umax,extent[1]),
				std::max(vmin,extent[2]), std::min(vmax,extent[3]), tol);
	  if (tol < 0.0)
	    use_stdd = true;
	}
      if (last > 0)
	{
	  tolin >> AEPSGE;
	  if (AEPSGE < 0.0)
	    use_stdd = true;
	}
    }
  if (AEPSGE < 0.0)
    use_stdd = true;
	  

  // Move point clouds to origo
  //Point mid(0.5*(domain[0] + domain[1]), 0.5*(domain[2] + domain[3]), 0.0);
  Point mid(0.0, 0.0, 0.0); //0.5*(domain[0] + domain[1]), 0.5*(domain[2] + domain[3]), 0.0);
  for (size_t ki=0; ki<data.size(); ++ki)
    {
      int nmb_pts = (int)data[ki].size()/3;
      for (int ka=0; ka<nmb_pts; ++ka)
	for (int kb=del-3; kb<del-1; ++kb)
	  {
	    data[ki][del*ka+kb] -= mid[kb-del+3];
	  }
    }
  double domain2[4];
  for (int ka=0; ka<2; ++ka)
    for (int kb=0; kb<2; ++kb)
      domain2[2*ka+kb] = domain[2*ka+kb]-mid[ka];

  for (size_t kj=0; kj<tolerances.size(); ++kj)
    tolerances[kj].translateBox(-mid[0], -mid[1]);

  if (/*true)/*/use_stdd)
    {
      for (size_t ki=0; ki<data.size(); ++ki)
	{
	  int nmb_pts = (int)data[ki].size()/3;
	  double avheight = 0.0;
	  for (int ka=0; ka<nmb_pts; ++ka)
	    {
	      double height = data[ki][del*ka+del-1];
	      avheight += (height/(double)nmb_pts);
	    }
	  double stdd = 0.0;
	  for (int ka=0; ka<nmb_pts; ++ka)
	    {
	      double height = data[ki][del*ka+del-1];
	      stdd += (pow(avheight-height,2)/(double)nmb_pts);
	    }
	  stdd = sqrt(stdd);
	  std::cout << "Standard deviation, " << ki+1 <<": " << stdd << std::endl;
	}
    }
     

  time_t start = time(NULL);



 boost::timer t;
  double duration;

  t.restart();

  // if (del > 3)
  //   {
  //     initmba = 0;
  //     mba = 0;
  //   }
  int order = degree + 1; 
  int nmb_coef = std::max(order, initncoef); //std::max(order, 14);
  int nm = nmb_coef*nmb_coef;
  double dom = (domain[1]-domain[0])*(domain[3]-domain[2]);
  double c1 = std::pow((double)nm/dom, 1.0/2.0);
  int nc[2];
  for (int kb=0; kb<2; ++kb)
    {
      double tmp = c1*(domain[2*kb+1]-domain[2*kb]);
      nc[kb] = (int)tmp;
      //if (tmp - (double)nc[kj] > 0.5)
      ++nc[kb];
      nc[kb] = std::max(nc[kb], order);
    }
  double mba_coef = 0.0;
  if (initmba)
    mba_coef = 0.5*(domain[2*(del-1)] + domain[2*(del-1)+1]);
  vector<shared_ptr<LRSplineSurface> > sfs;
  for (size_t ki=0; ki<data.size(); ++ki)
    {
      shared_ptr<LRSurfApprox> approx;
      if (distribute_ncoef)
	approx = shared_ptr<LRSurfApprox>(new LRSurfApprox(nc[0], order, nc[1], order, data[ki], del-2, 
							   domain2, AEPSGE, initmba ? true : false, 
							   mba_coef, true, true));
      else
	approx = shared_ptr<LRSurfApprox>(new LRSurfApprox(nmb_coef, order, nmb_coef, order, data[ki], del-2, 
							   domain2, AEPSGE, initmba ? true : false, 
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
      if (reltol > 0)
	{
	  double tol1 = domain[4]<0.0 ? AEPSGE - tolfac2*domain[4] : AEPSGE + tolfac1*domain[4];
	  double tol2 = domain[5]<0.0 ? AEPSGE - tolfac2*domain[5] : AEPSGE + tolfac1*domain[5];
	  std::cout << "Variable tolerance: " << tol1 << " - " << tol2 << std::endl;
	  approx->setVarTol(tolfac1, tolfac2);
	}

      if (tolerances.size() > 0)
	approx->setVarTolBox(tolerances);
      
      approx->setVerbose(verbose);


      // Refine strategy
      approx->setRefinementStrategy(refcat1, alter, threshold1, swap, refcat2, threshold2);
      
      if (del == 3)
	{
	  double zrange = domain[5] - domain[4];
	  double zfac = std::max(AEPSGE, 0.005*zrange);
	  approx->addLowerConstraint(domain[4] - zfac);
	  approx->addUpperConstraint(domain[5] + zfac);
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

      sfs.push_back(surf);
    }

  std::ofstream of("sfs.g2");
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      sfs[ki]->writeStandardHeader(of);
      sfs[ki]->write(of);
    }

  // Define surfaces on the same mesh
  LRSurfUtils::defineOnSameMesh(sfs);

  // Update surface approximation with extended mesh
  vector<double> maxdist(data.size()), avdist(data.size()), avdist_total(data.size());
  vector<int> nmb_out_eps(data.size());
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      LRSurfApprox approx(sfs[ki], data[ki], AEPSGE, true, true, true);
      approx.setSmoothingWeight(smoothwg);
      approx.setSmoothBoundary(true);
      if (mba || max_iter >= tomba)
	approx.setUseMBA(true);
      else
	{
	  approx.setSwitchToMBA(tomba-max_iter);
	}
      approx.setVerbose(true);

      if (del == 3)
	{
	  double zrange = extent[5] - extent[4];
	  approx.addLowerConstraint(extent[4] - 0.1*(zrange));
	  approx.addUpperConstraint(extent[5] + 0.1*(zrange));
	}

      shared_ptr<LRSplineSurface> surf = approx.getApproxSurf(maxdist[ki], avdist_total[ki],
							      avdist[ki], nmb_out_eps[ki], 
							      0);
      sfs[ki] = surf;
     }
  
  duration = t.elapsed();
  std::cout << "Duration: " << duration << std::endl;
  double min = floor(duration/60);
  double sec = duration - 60*min;
  std::cout << min << "m" << sec << "s" << std::endl;
  time_t end = time(NULL);
  std::cout<<"Execution Time: "<< (double)(end-start)<<" Seconds"<<std::endl;

  std::ofstream sfout(surffile);     // Surface output stream
  if (infofile)
    {
      std::ofstream infoout(infofile);   // Accuracy information output stream
      for (size_t ki=0; ki<data.size(); ++ki)
	{
	  infoout << "Data set number " << ki+1 << std::endl;
	  infoout << "Total number of points: " << data[ki].size()/3 << std::endl;
	  infoout << "Number of elements: " << sfs[ki]->numElements() << std::endl;
	  infoout << "Number of coefficients: " << sfs[ki]->numBasisFunctions() << std::endl;
	  infoout << "Maximum distance: " << maxdist[ki] << std::endl;
	  infoout << "Average distance: " << avdist_total[ki] << std::endl;
	  infoout << "Average distance for points outside of the tolerance: " << avdist[ki] << std::endl;
	  infoout << "Number of points outside the tolerance: " << nmb_out_eps[ki] << std::endl;
	}
    }
  else
    {
      for (size_t ki=0; ki<data.size(); ++ki)
	{
	  std::cout << "INFO: Total number of points: " << data[ki].size()/3 << std::endl;
	  std::cout << "INFO: Number of elements: " << sfs[ki]->numElements() << std::endl;
	  std::cout << "INFO: Number of coefficients: " << sfs[ki]->numBasisFunctions() << std::endl;
	  std::cout << "INFO: Maximum distance: " << maxdist[ki] << std::endl;
	  std::cout << "INFO: Average distance: " << avdist_total[ki] << std::endl;
	  std::cout << "INFO: Average distance for points outside of the tolerance: " << avdist[ki] << std::endl;
	  std::cout << "INFO: Number of points outside the tolerance: " << nmb_out_eps[ki] << std::endl;
	}
    } 


  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      // Translate
      if (sfs[ki]->dimension() == 3)
	{
	  sfs[ki]->translate(mid);
	}
      else
	{
	  // Update parameter domain
	  double umin = sfs[ki]->paramMin(XFIXED);
	  double umax = sfs[ki]->paramMax(XFIXED);
	  double vmin = sfs[ki]->paramMin(YFIXED);
	  double vmax = sfs[ki]->paramMax(YFIXED);

	  sfs[ki]->setParameterDomain(umin + mid[0], umax + mid[0],
				   vmin + mid[1], vmax + mid[1]);
	}

      sfs[ki]->writeStandardHeader(sfout);
      sfs[ki]->write(sfout);

    }
  return 0;
}

