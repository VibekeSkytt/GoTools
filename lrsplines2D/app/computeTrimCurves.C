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

#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include "GoTools/lrsplines2D/TrimUtils.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

using namespace Go;

using std::cout;
using std::endl;
using std::vector;
using std::string;

//#define DEBUG

void print_help_text()
{
  std::cout << "computeTrimCurves computes 2D trimming curves with respect to a point cloud. \n";
  std::cout << "A tightness parameter governs how close to the point cloud the trimming curve will pass. \n";
  std::cout << "The density and uniformity of the cloud should be reflected in this parameter. \n";
  std::cout << "The parameter should be in the range [1:8] and higher for denser point clouds. \n";
  std::cout << "The smaller number should be used only when the point cloud is very ragged towards the boundary. \n";
  std::cout << "Usage: Input cloud (.txt, .g2), tightness parameter, output curves (.g2) \n";
  std::cout << "-h or --help for help text" << std::endl;
}


int main(int argc, char* argv[])
{
  for (int kii=1; kii<argc; ++kii)
    {
      string arg(argv[kii]);
      if (arg == "-h" || arg == "--help")
	{
	  print_help_text();
	  exit(0);
	}
    }

  if (argc != 4 && argc != 5)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

    
  char *pointfile(argv[1]);
  int tightness = atoi(argv[2]);
  std::ofstream fileout(argv[3]);
  int nmb_div = -1;
  if (argc == 5)
    nmb_div = atoi(argv[4]);


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

  // Set parameters for computations of trimming sequence
  int max_rec;
  if (tightness <= 2)
    {
      max_rec = 1;
      if (nmb_div < 0)
	nmb_div = (tightness == 2) ? 20 : 15;
    }
  else if (tightness <= 5)
    {
      max_rec = 2;
      if (nmb_div < 0)
	nmb_div = (tightness == 3) ? 8 : ((tightness == 4) ? 12 : 15);
    }
  else
    {
      max_rec = 3;
      if (nmb_div < 0)
	nmb_div = (tightness == 6) ? 10 : ((tightness == 7) ? 12 : 15);
    }


  // Compute trimming seqence
  bool only_outer = true;
  vector<vector<double> > seqs;
  //TrimUtils trimutil(points2, 1, domain);
  TrimUtils trimutil(&data[0], nmb_pts, 1, &extent[0]);
  trimutil.computeTrimSeqs(max_rec, nmb_div, seqs, only_outer);
  
  double udel, vdel;
  trimutil.getDomainLengths(udel, vdel);

  // @@@ VSK 0115. Use only the first sequence. If there are more, they
  // are likely to be incomplete.
  // if (only_outer && seqs.size() > 1)
  //   seqs.erase(seqs.begin()+1, seqs.end());
  int nmb_loops = only_outer ? 1 : (int)seqs.size();
  vector<vector<vector<double> > >  all_seqs(nmb_loops);
  for (int kh=0; kh<nmb_loops; ++kh)
    all_seqs[kh].push_back(seqs[kh]);


  // Compute trimming loop
  // First extract parts of the trimming sequences following iso trim curves
  int nmb_match = 4;
  double tol = 1.0e-5;
  double eps = std::max(udel, vdel);
  vector<vector<shared_ptr<SplineCurve> > > loop(nmb_loops);
  for (int kh=0; kh<nmb_loops; ++kh)
    {
      for (int kj=0; kj<4; ++kj)
      	{
      	  int ix = (kj < 2) ? 0 : 1;
      	  for (int ki=0; ki<(int)all_seqs[kh].size();)
      	    {
      	      vector<vector<double> > split_seqs1 = 
      		TrimCrvUtils::extractConstParSeqs(all_seqs[kh][ki], ix, 
      						  extent[kj], nmb_match, 
      						  tol, eps);
      	      if (split_seqs1.size() > 1 || split_seqs1[0].size() == 4)
      		{
      		  all_seqs[kh].erase(all_seqs[kh].begin()+ki);
      		  all_seqs[kh].insert(all_seqs[kh].begin()+ki, 
      				      split_seqs1.begin(), split_seqs1.end());
      		  ki += (int)split_seqs1.size();
      		}
      	      else
      		++ki;
      	    }
      	}

      // Split the remaining sequences from the outer loop in kinks
      double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
      vector<vector<double> > split_seqs;
      for (size_t ki=0; ki<all_seqs[kh].size(); ++ki)
	{
	  vector<vector<double> > curr_seqs = 
	    TrimCrvUtils::splitCurvePointsInKinks(all_seqs[kh][ki], kink_tol);
	  split_seqs.insert(split_seqs.end(), curr_seqs.begin(), curr_seqs.end());
	}

      // Ensure a closed trimming loop
      TrimCrvUtils::makeConnectedLoop(split_seqs, tol);
  
      // Create trimming curves
      const int par_dim = 2;
      const int max_iter = 5;
      vector<shared_ptr<SplineCurve> > par_cvs;
      for (size_t ki = 0; ki < split_seqs.size(); ++ki)
	{
	  shared_ptr<SplineCurve> spline_cv_appr_2d
	    (TrimCrvUtils::approximateTrimPts(split_seqs[ki], par_dim, eps, 
					      max_iter));
	  par_cvs.push_back(spline_cv_appr_2d);
	}
  
      // The curve should be CCW.
      // @@@ VSK 150310. Only outer curve loops
      // Assume one outer and the rest inner (this is not necessarily true)
      const double int_tol = 1e-06;
      vector<shared_ptr<ParamCurve> > par_cvs2(par_cvs.begin(), par_cvs.end());
      bool loop_is_ccw = LoopUtils::loopIsCCW(par_cvs2, eps, int_tol);
      if ((kh==0 && !loop_is_ccw) || (kh>0 && loop_is_ccw))
	{
	  //MESSAGE("We should change direction of the loop cv!");
	  for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	    {
	      par_cvs[ki]->reverseParameterDirection();
	    }
	  reverse(par_cvs.begin(), par_cvs.end());
	}
      loop[kh] = par_cvs;
    }
      
   // Write output curves
  for (int kh=0; kh<nmb_loops; ++kh)
    {
      for (size_t kr=0; kr<loop[kh].size(); ++kr)
	{
	  loop[kh][kr]->writeStandardHeader(fileout);
	  loop[kh][kr]->write(fileout);
	}
    }

   return 0;
}
