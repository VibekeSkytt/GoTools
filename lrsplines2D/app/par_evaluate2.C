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
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

void polynomialTerms(double u, double v, int num, int max_num,
		     vector<double>& terms)
{
  vector<double> u_terms(num+1), v_terms(num+1);
  u_terms[0] = v_terms[0] = 1.0;

  for (int ki=1; ki<=num; ++ki)
    {
      u_terms[ki] = u_terms[ki-1]*u;
      v_terms[ki] = v_terms[ki-1]*v;
    }

  for (int ki=0, kr=0; ki<=num; ++ki)
    for (int kj=0; kj<=num; ++kj)
      {
	if (ki+kj <= max_num)
	  terms[kr++] = u_terms[ki]*v_terms[kj];
      }
}

int main(int argc, char *argv[])
{
  if (argc < 7) {
    std::cout << "Parameters : Input points (.txt), no. of regular points, output points (.txt), degree, total degree, polynomial factors"  << std::endl;
    exit(-1);
  }
  std::string pointfile(argv[1]);
  int num = atoi(argv[2]);
  std::ofstream pointout(argv[3]);
  int degree = atoi(argv[4]);
  int tot_degree = atoi(argv[5]);
  int num_terms = (tot_degree == degree) ? (degree+1)*(degree+2)/2 :
    (degree+1)*(degree+1);
  if (argc != num_terms + 6)
    {
      std::cout << "Expecting " << num_terms << " polynomial factors" << std::endl;
      exit(-1);
    }

  pointout.precision(18);
  
  num_terms = std::min(num_terms, 16);
  vector<double> pol(num_terms);
  for (int ka=0; ka<num_terms; ++ka)
    pol[ka] = atof(argv[6+ka]);

		   
  // Read parameter values (and points)
  int del = 5;
  int nmb_pts = 0;
  vector<double> data;
  vector<double> extent(2*del);   // Limits for points in all coordinates
  if (pointfile != "no")
    {
      std::ifstream is(pointfile.c_str());
      FileUtils::readTxtPointFile(is, del, data, nmb_pts, extent);
    }

  vector<double> term(num_terms);
  
  // Evaluate
  if (nmb_pts > 0)
    {
      for (int ki=0; ki<nmb_pts; ++ki)
	{
	  polynomialTerms(data[ki*del], data[ki*del+1], degree, tot_degree, term);
	  double val = 0.0;
	  for (int ka=0; ka<num_terms; ++ka)
	    val += pol[ka]*term[ka];
	  pointout << data[ki*del] << " " <<  data[ki*del+1] << " " << val << std::endl;
	}
    }
  else
    {
      extent[0] = extent[2] = 0.0;
      extent[1] = extent[3] = 1.0;
    }
  if (num > 1)
    {
      double udel = (extent[1] - extent[0])/(double)(num - 1);
      double vdel = (extent[3] - extent[2])/(double)(num - 1);

      double u=extent[0], v=extent[2];
      int ka, kb;
      for (kb=0; kb<num; ++kb, v+=vdel)
	for (ka=0, u=extent[0]; ka<num; ++ka, u+=udel)
	  {
	    polynomialTerms(u, v, degree, tot_degree, term);
	    double val = 0.0;
	    for (int ka=0; ka<num_terms; ++ka)
	      val += pol[ka]*term[ka];
	    pointout << u << " " << v << " " << val << std::endl;
	  }
    }
}

