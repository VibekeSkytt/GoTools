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

void polynomialTerms(double u, double v, double w,
		     int num1, int num2, int num3, int max_num,
		     vector<double>& terms)
{
  vector<double> u_terms(num1+1), v_terms(num2+1), w_terms(num3+1);
  u_terms[0] = v_terms[0] = w_terms[0] = 1.0;

  for (int ki=1; ki<=num1; ++ki)
    u_terms[ki] = u_terms[ki-1]*u;
  for (int ki=1; ki<=num2; ++ki)
     v_terms[ki] = v_terms[ki-1]*v;
  for (int ki=1; ki<=num3; ++ki)
     w_terms[ki] = w_terms[ki-1]*w;

  for (int kh=0, kr=0; kh<=num3; ++kh)
    for (int kj=0; kj<=num2; ++kj)
      for (int ki=0; ki<=num1; ++ki)
	{
	  if (ki+kj+kh <= max_num)
	    terms[kr++] = u_terms[ki]*v_terms[kj]*w_terms[kh];
	}
}

void Chebyshev(int degree1, int degree2, int degree3, vector<double>& pol)
{
  int ki, kj, kh, kr;
  vector<vector<double> > Tc(6);
  Tc[0].push_back(1.0);
  Tc[1].insert(Tc[1].end(), {0.0, 1.0});
  Tc[2].insert(Tc[2].end(), {-1.0, 0.0, 2.0});
  Tc[3].insert(Tc[3].end(), {0.0, -3.0, 0.0, 4.0});
  Tc[4].insert(Tc[4].end(), {1.0, 0.0, -8.0, 0.0, 8.0});
  Tc[5].insert(Tc[5].end(), {0.0, 5.0, 0.0, -20.0, 0.0, 16.0});
  int deg1 = std::min(degree1, 5);
  int deg2 = std::min(degree2, 5);
  int deg3 = std::min(degree3, 5);
  for (kh=0, kr=0; kh<=deg3; ++kh)
    for (kj=0; kj<=deg2; ++kj)
      for (ki=0; ki<=deg1; ++ki)
	{
	  pol[kr++] = Tc[deg1][ki]*Tc[deg2][kj]*Tc[deg3][kh];
	}
}

int countTerms(int deg1, int deg2, int deg3, int tot)
{
  int num = 0;
  for (int kh=0, kr=0; kh<=deg3; ++kh)
    for (int kj=0; kj<=deg2; ++kj)
      for (int ki=0; ki<=deg1; ++ki)
	{
	  if (ki+kj+kh <= tot)
	    ++num;
	}
  return num;
}

int main(int argc, char *argv[])
{
  if (argc < 7) {
    std::cout << "Parameters : No. of regular points, output points (.txt), extended total degree (0/1/2 (Chebyshev polynomial)), degree1, degree2, degree3, polynomial factors"  << std::endl;
    exit(-1);
  }
  int num = atoi(argv[1]);
  std::ofstream pointout(argv[2]);
  int total = atoi(argv[3]);
  int degree1 = atoi(argv[4]);
  int degree2 = atoi(argv[5]);
  int degree3 = atoi(argv[6]);
  int maxdeg = std::max(degree1,std::max(degree2,degree3));
  int mindeg = std::min(degree1,std::min(degree2,degree3));
  int tot_degree = (total) ? degree1*degree2*degree3 : maxdeg;
  int num_terms = (tot_degree == maxdeg) ?
    countTerms(degree1, degree2, degree3, maxdeg) : (degree1+1)*(degree2+1)*(degree3+1);
  if (total != 2 && argc != num_terms + 7)
    {
      std::cout << "Expecting " << num_terms << " polynomial factors" << std::endl;
      exit(-1);
    }

  pointout.precision(18);
  
  vector<double> pol(num_terms, 0.0);
  if (total == 2)
    Chebyshev(degree1, degree2, degree3, pol);
  else
    {
      for (int ka=0; ka<num_terms; ++ka)
	pol[ka] = atof(argv[7+ka]);
    }

		   
  vector<double> range(6);
  range[0] = range[2] = range[4] = -1;
  range[1] = range[3] = range[5] = 1;
  
  vector<double> term(num_terms, 0.0);
  
  // Evaluate
  if (num > 1)
    {
      double udel = (range[1] - range[0])/(double)(num - 1);
      double vdel = (range[3] - range[2])/(double)(num - 1);
      double wdel = (range[5] - range[4])/(double)(num - 1);

      double u=range[0], v=range[2], w=range[4];
      int ka, kb, kc, kd;
      for (kc=0; kc<num; ++kc, w+=wdel)
	for (kb=0, v=range[2]; kb<num; ++kb, v+=vdel)
	  for (ka=0, u=range[0]; ka<num; ++ka, u+=udel)
	  {
	    polynomialTerms(u, v, w, degree1, degree2, degree3, tot_degree, term);
	    double val = 0.0;
	    for (kd=0; kd<num_terms; ++kd)
	      val += pol[kd]*term[kd];
	    pointout << u << " " << v << " " << w << " " << val << std::endl;
	  }
    }
}

