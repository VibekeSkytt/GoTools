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
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/FileUtils.h"
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

int compare_extension(const char *str1, char str2[][8], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}


int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: foldef in, folder out" << std::endl;
    return -1;
  }

  std::string inpath(argv[1]);
  std::string outpath(argv[2]); 
  vector<string> infile;
  vector<string> outfile;
  std::string separator = "/";
  std::string extension = ".g2";

  // Possible types of input files
  char keys[12][8] = {"txt", "TXT", "xyz", "XYZ", "dat", "asc"};

  // Read data files
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir (inpath.c_str())) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      //  Find file extension
      char *loc;
      char *last = 0;
      loc = strchr(ent->d_name, '.');
      while (loc != NULL)
	{
	  last = loc;
	  loc = strchr(loc+1, '.');
	}
      if (last == NULL)
	continue;

      char *input_type = last+1;

      // Check type
      int type;
      try {
	type = FileUtils::fileType(ent->d_name, keys, 6);
      }
      catch (...)
	{
	  continue;
	}

      if (type >=0)
	{
	  string data_file(ent->d_name);
	  string data_file2;
	  combinePathFile(inpath, data_file, data_file2);
	  infile.push_back(data_file2);
	  string data_file3 = outpath + separator + data_file + extension;
	  outfile.push_back(data_file3);
	}
    }
    closedir (dir);
  } else {
    /* could not open directory */
    std::cout << "ERROR: Could not open file directory" << std::endl;
    return 1;
  }

  vector<int> nmb_points;
  vector<double> bb;
  vector<string> file_name;

  double u1, u2, v1, v2;
  PointCloud3D points;
  for (size_t kj=0; kj<infile.size(); ++kj)
    {
      // Read point file
      vector<double> data;
      vector<double> extent(6);
      int nmb_pts;
      try {
	std::ifstream pointsin(infile[kj].c_str());
	FileUtils::readTxtPointFile(pointsin, 3, data, nmb_pts, extent);
      }
      catch (...)
	{
	  std::cout << "ERROR: Failed reading point file: " << infile[kj] << std::endl;
	  continue;
	}

   
      points = PointCloud3D(data.begin(), data.size()/3);


      // Write to g2 format
      std::ofstream ptsout(outfile[kj].c_str());
      points.writeStandardHeader(ptsout);
      points.write(ptsout);
    }
}

