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
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>

using namespace Go;
using std::vector;
using std::string;

void combinePathFile(string& path, string& file, string& file2)
{
  char outfile[160];
  strcpy(outfile, path.c_str());
  char tmp1[2];
  sprintf(tmp1,"/");
  strncat(outfile, tmp1, 1);
  strcat(outfile, file.c_str());
  file2 = std::string(outfile);
}


int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: foldef in (.g2 files), trivariate point + hight out, time step" << std::endl;
    return -1;
  }

  std::string inpath(argv[1]);
  std::ofstream outfile(argv[2]);
  double zdel = atof(argv[3]);
  vector<string> infile;
  std::string separator = "/";
  std::string extension = ".g2";

  // Read point files
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir (inpath.c_str())) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      string data_file(ent->d_name);
      if (data_file == "." || data_file == "..")
	continue;
      string data_file2;
      combinePathFile(inpath, data_file, data_file2);
      infile.push_back(data_file2);
    }
  } else {
    /* could not open directory */
    std::cout << "ERROR: Could not open file directory" << std::endl;
    return 1;
  }

  // Sort point files
  std::sort(infile.begin(), infile.end());

  // Read point data
  vector<PointCloud3D> clouds(infile.size());
  int nmb_pts = 0;
  double domain[4];
  domain[0] = domain[2] = std::numeric_limits<double>::max();
  domain[1] = domain[3] = std::numeric_limits<double>::lowest();
  for (size_t kj=0; kj<infile.size(); ++kj)
    {
      std::ifstream inpts(infile[kj]);
      ObjectHeader header;
      header.read(inpts);
      clouds[kj].read(inpts);

      nmb_pts += clouds[kj].numPoints();
      BoundingBox bb = clouds[kj].boundingBox();

      domain[0] = std::min(domain[0], bb.low()[0]);
      domain[1] = std::max(domain[1], bb.high()[0]);
      domain[2] = std::min(domain[2], bb.low()[1]);
      domain[3] = std::max(domain[3], bb.high()[1]);
    }

  //double zrange = 2.0*(double)(infile.size()-1); 
  //double zrange = 0.5*(domain[1]+domain[3]-domain[0]-domain[2]);
  //double zrange = 0.25*(domain[1]+domain[3]-domain[0]-domain[2]);
  //double zdel = zrange/(double)(clouds.size()-1);

  // Write 4D point cloud
  outfile << nmb_pts << std::endl;
  for (size_t kj=0; kj<clouds.size(); ++kj)
    {
      double *data = clouds[kj].rawData();
      int npt3 = 3*clouds[kj].numPoints();
      for (int ka=0; ka<npt3; ka += 3)
	{
	  int kb;
	  for (kb=0; kb<2; ++kb)
	    outfile << data[ka+kb] << " ";
	  outfile << kj*zdel << " ";
	  for (; kb<3; ++kb)
	    outfile << data[ka+kb] << " ";
	  outfile << std::endl;
	}
    }
}

  
