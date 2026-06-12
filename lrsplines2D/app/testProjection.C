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
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRProjection.h"
#include "sislP.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

void polynomialTerms(double u, double v, int num1, int num2, int max_num,
		     double *terms)
{
  vector<double> u_terms(num1+1), v_terms(num2+1);
  u_terms[0] = v_terms[0] = 1.0;

  for (int ki=1; ki<=num1; ++ki)
    u_terms[ki] = u_terms[ki-1]*u;
  for (int ki=1; ki<=num2; ++ki)
     v_terms[ki] = v_terms[ki-1]*v;

  for (int kj=0, kr=0; kj<=num2; ++kj)
    for (int ki=0; ki<=num1; ++ki)
      {
	if (ki+kj <= max_num)
	  terms[kr++] = u_terms[ki]*v_terms[kj];
      }
}

void Chebyshev(int degree1, int degree2, vector<double>& pol)
{
  int ki, kj, kr;
  vector<vector<double> > Tc(6);
  Tc[0].push_back(1.0);
  Tc[1].insert(Tc[1].end(), {0.0, 1.0});
  Tc[2].insert(Tc[2].end(), {-1.0, 0.0, 2.0});
  Tc[3].insert(Tc[3].end(), {0.0, -3.0, 0.0, 4.0});
  Tc[4].insert(Tc[4].end(), {1.0, 0.0, -8.0, 0.0, 8.0});
  Tc[5].insert(Tc[5].end(), {0.0, 5.0, 0.0, -20.0, 0.0, 16.0});
  int deg1 = std::min(degree1, 6);
  int deg2 = std::min(degree2, 6);
  for (kj=0, kr=0; kj<=deg2; ++kj)
    for (ki=0; ki<=deg1; ++ki)
      {
	pol[kr++] = Tc[deg1][ki]*Tc[deg2][kj];
      }

  int num = 50;
  std::ofstream ofx("Cpolx.g2");
  std::ofstream ofy("Cpoly.g2");
  ofx << "410 1 0 0" << std::endl;
  ofx << num-1 << std::endl;
  double x = -1;
  double x2 = -1;
  double del = 2.0/(double)(num-1);
  double val = Tc[deg1][0];
  for (ki=1; ki<=deg1; ++ki, x2*=x)
    val += Tc[deg1][ki]*x2;
  Point prev(x, val, 0.0);
  for (kj=1, x+=del; kj<num; ++kj, x+=del)
    {
      val = Tc[deg1][0];
      x2 = x;
      for (ki=1; ki<=deg1; ++ki, x2*=x)
	val += Tc[deg1][ki]*x2;
      Point curr(x, val, 0.0);
      ofx << prev << " " << curr << std::endl;
      prev = curr;
    }
  
  ofy << "410 1 0 0" << std::endl;
  ofy << num-1 << std::endl;
  double y = -1;
  double y2 = -1;
  val = Tc[deg2][0];
  for (ki=1; ki<=deg2; ++ki, y2*=y)
    val += Tc[deg2][ki]*y2;
   prev.setValue(y, val, 0.0);
  for (kj=1, y+=del; kj<num; ++kj, y+=del)
    {
      val = Tc[deg2][0];
      y2 = y;
      for (ki=1; ki<=deg2; ++ki, y2*=y)
	val += Tc[deg2][ki]*y2;
      Point curr(y, val, 0.0);
      ofy << prev << " " << curr << std::endl;
      prev = curr;
    }
}

int main(int argc, char *argv[])
{
  if (argc < 8) {
    std::cout << "Parameters: lr func in (.g2), lr func out (.g2), use computed polynomial (0/1), num_sample, degree1, degree2, extended total degree (0/1/2 (Chebyshev polynomial)), polynomial factors"  << std::endl;
    exit(-1);
  }
  std::ifstream sf_in(argv[1]);
  std::ofstream sf_out(argv[2]);
  bool use_comp = atoi(argv[3]);
  int num_sample = atoi(argv[4]);
  int degree1 = atoi(argv[5]);
  int degree2 = atoi(argv[6]);
  int total = atoi(argv[7]);
  int maxdeg = std::max(degree1,degree2);
  int mindeg = std::min(degree1,degree2);
  int tot_degree = (total) ? degree1*degree2 : maxdeg;
  int num_terms = (tot_degree == maxdeg) ?
    maxdeg*(mindeg+1) - mindeg*(mindeg-1)/2 + 1 : (degree1+1)*(degree2+1);
  if (total != 2 && argc != num_terms + 8)
    {
      std::cout << "Expecting " << num_terms << " polynomial factors" << std::endl;
      exit(-1);
    }

  ObjectHeader header;
  header.read(sf_in);
  shared_ptr<LRSplineSurface> lrsf(new LRSplineSurface());
  lrsf->read(sf_in);
  
  //num_terms = std::min(num_terms, 16);
  vector<double> pol(num_terms, 0.0);
  if (total == 2)
    Chebyshev(degree1, degree2, pol);
  else
    {
      for (int ka=0; ka<num_terms; ++ka)
	pol[ka] = atof(argv[7+ka]);
    }

  double u1 = -1.0, u2 = 1.0, v1 = -1.0, v2 = 1.0;
  lrsf->setParameterDomain(u1, u2, v1, v2);
  
  vector<Point> par(num_terms);
  if (num_terms == 1)
    par[0] = Point(0.0, 0.0);
  else if (total)
    {
      double udel = (u2 - u1)/(double)degree1;
      double vdel = (v2 - v1)/(double)degree2;
      double u, v;
      int ki, kj, kr;
      for (kj=0, kr=0, v=v1; kj<=degree2; ++kj, v+=vdel)
	{
	  v = std::min(v, v2);
	  for (ki=0, u=u1; ki<=degree1; ++ki, ++kr, u+=udel)
	    {
	      u = std::min(u, u2);
	      par[kr] = Point(u, v);
	    }
	}
    }
  else
    {
      int tot_num = (maxdeg+1)*(maxdeg+2)/2;
      vector<Point> tmp_par;
      tmp_par.reserve(tot_num);
      double tdel = 1.0/(double)tot_num;
      int ki;
      vector<double> u(maxdeg+1), v(maxdeg+2);
      for (ki=0; ki<=maxdeg; ++ki)
	u[ki] = cos(ki*M_PI/(double)maxdeg);
      for (ki=0; ki<=maxdeg+1; ++ki)
	v[ki] = cos(ki*M_PI/(double)(maxdeg+1));
      for (ki=0; ki<=maxdeg; ++ki)
	{
	  if (ki%2 == 0)
	    {
	      for (int kj=1; kj<(int)v.size(); kj+=2)
		tmp_par.push_back(Point(u[ki], v[kj]));
	    }
	  else
	    {
	      for (int kj=0; kj<(int)v.size(); kj+=2)
		tmp_par.push_back(Point(u[ki], v[kj]));
	    }
	}

      // for (size_t kr=0; kr<tmp_par.size(); ++kr)
      // 	tmp_par[kr] = Point((tmp_par[kr][0]+1.0)/2.0,(tmp_par[kr][1]+1.0)/2.0); 
				  
      // double u, v, t;
      // for (ki=0, t=0; ki<tot_num; ++ki, t+=tdel)
      // 	{
      // 	  u = cos(ki*M_PI/maxdeg); //0.5*(1.0 - cos(maxdeg*t));
      // 	  v = cos( 0.5*(1.0 - cos((maxdeg+1)*t));
      // 	  tmp_par.push_back(Point(u, v));
      // 	}

      if (tot_num > num_terms)
	{
	  int diff = tot_num - num_terms;
	  int kdel = tot_num/diff;
	  int kcurr = tot_num - kdel/2 - 1;
	  for (ki=0; ki<diff; ++ki, kcurr-=kdel)
	    tmp_par.erase(tmp_par.begin()+kcurr);
	}
      for (ki=0; ki<num_terms; ++ki)
	par[ki] = tmp_par[ki];

      int stop_break = 1;
    }
  
  // par[1].push_back(Point(0.25, 0.25));
  // par[1].push_back(Point(0.75, 0.25));
  // par[1].push_back(Point(0.5, 0.75));
  
  // par[2].push_back(Point(0.25, 0.25));
  // par[2].push_back(Point(0.75, 0.25));
  // par[2].push_back(Point(0.25, 0.75));
  // par[2].push_back(Point(0.75, 0.75));
  
  // par[3].push_back(Point(0.0, 0.0));
  // par[3].push_back(Point(0.75, 0.0));
  // par[3].push_back(Point(0.25, 0.5));
  // par[3].push_back(Point(1.0, 0.5));
  // par[3].push_back(Point(0.0, 1.0));
  // par[3].push_back(Point(1.0, 1.0));
  
  // par[4].push_back(Point(0.0, 0.0));
  // par[4].push_back(Point(0.5, 0.0));
  // par[4].push_back(Point(1.0, 0.0));
  // par[4].push_back(Point(0.0, 0.5));
  // par[4].push_back(Point(0.5, 0.5));
  // par[4].push_back(Point(1.0, 0.5));
  // par[4].push_back(Point(0.0, 1.0));
  // par[4].push_back(Point(0.5, 1.0));
  // par[4].push_back(Point(1.0, 1.0));
  
  // par[5].push_back(Point(0.0, 0.0));
  // par[5].push_back(Point(0.3, 0.0));
  // par[5].push_back(Point(0.6, 0.0));
  // par[5].push_back(Point(0.3, 0.3));
  // par[5].push_back(Point(0.8, 0.3));
  // par[5].push_back(Point(0.0, 0.7));
  // par[5].push_back(Point(0.5, 0.7));
  // par[5].push_back(Point(1.0, 0.7));
  // par[5].push_back(Point(0.0, 1.0));
  // par[5].push_back(Point(0.5, 1.0));
  // par[5].push_back(Point(1.0, 1.0));
  
  // par[6].push_back(Point(0.0, 0.0));
  // par[6].push_back(Point(0.25, 0.0));
  // par[6].push_back(Point(0.75, 0.0));
  // par[6].push_back(Point(1.0, 0.0));
  // par[6].push_back(Point(0.0, 0.25));
  // par[6].push_back(Point(0.25, 0.25));
  // par[6].push_back(Point(0.75, 0.25));
  // par[6].push_back(Point(1.0, 0.25));
  // par[6].push_back(Point(0.0, 0.75));
  // par[6].push_back(Point(0.25, 0.75));
  // par[6].push_back(Point(0.75, 0.75));
  // par[6].push_back(Point(1.0, 0.75));
  // par[6].push_back(Point(0.0, 1.0));
  // par[6].push_back(Point(0.25, 1.0));
  // par[6].push_back(Point(0.75, 1.0));
  // par[6].push_back(Point(1.0, 1.0));

  vector<double> ls(num_terms*num_terms, 0.0), rs(num_terms, 0.0);
  vector<int> piv(num_terms);
  for (int ka=0; ka<num_terms; ++ka)
    piv[ka] = ka;

  for (int ka=0; ka<num_terms; ++ka)
    {
      polynomialTerms(par[ka][0], par[ka][1], degree1, degree2,
		      tot_degree, &ls[ka*num_terms]);
      double val = 0.0;
      for (int kb=0; kb<num_terms; ++kb)
	    val += pol[kb]*ls[ka*num_terms+kb];
      rs[ka] = val;
    }

    // Solve
  int kstat = 0;
  s6lufacp(&ls[0], &piv[0], num_terms, &kstat);
  if (kstat < 0)
    {
      std::cout << "Singular equation system" << std::endl;
      exit(-1);
    }

  s6lusolp(&ls[0], &rs[0], &piv[0], num_terms, &kstat);
  if (kstat < 0)
    {
      std::cout << "Failed solving equation system" << std::endl;
      exit(-1);
    }

  std::cout << "Point distances: " << std::endl;
  vector<double> tmp(num_terms);
  for (int ka=0; ka<num_terms; ++ka)
    {
      polynomialTerms(par[ka][0], par[ka][1], degree1, degree2,
		      tot_degree, &tmp[0]);
      double val1 = 0.0, val2 = 0.0;
      for (int kb=0; kb<num_terms; ++kb)
	val1 += pol[kb]*tmp[kb];
      for (int kb=0; kb<num_terms; ++kb)
	val2 += rs[kb]*tmp[kb];
      std::cout << val1 << " " << val2 << " " << val1-val2 << std::endl;
    }
  
  std::cout << "Differences: " << std::endl;
  for (int ka=0; ka<num_terms; ++ka)
    {
      double diff = rs[ka] - pol[ka];
      std::cout << pol[ka] << " " << rs[ka] << " " << diff << std::endl;
    }

    int max_level = 10;  // Should always be enough
    for (int level=0; level<max_level; ++level)
      {
	int num_update = 0;
	for (auto bspl=lrsf->basisFunctionsBegin();
	     bspl!=lrsf->basisFunctionsEnd(); ++bspl)
	  {
	    int blevel = bspl->second->getNestLevel();
	    if (blevel != level)
	      continue;

	    Point coef = LRProjection::Polynomial2Coef((use_comp) ? rs : pol,
						       degree1, degree2, tot_degree,
						       u1, u2, v1, v2, bspl->second.get());
	    if (blevel > 0)
	      {
		std::cout << "Level: " << blevel << std::endl;
		
		bool OK = bspl->second->adaptProjCoef(coef);
	      }
	    lrsf->setCoef(coef, bspl->second.get());
	    num_update++;
	  }
	if (num_update == 0)
	  break;
      }

    std::ofstream of("sample_pts.g2");
    of.precision(15);
    of << "400 1 0 4 255 0 0 255" << std::endl;
    of << num_sample*num_sample << std::endl;
    double udel = (u2 - u1)/(double)(num_sample-1);
    double vdel = (v2 - v1)/(double)(num_sample-1);
    double u, v;
    int ki, kj;
    double maxdist = 0.0, avdist = 0.0;
    double fac = 1.0/(double)(num_sample*num_sample);
    for (kj=0, v=v1; kj<num_sample; ++kj, v+=vdel)
      {
	v = std::min(v, v2);
	for (ki=0, u=u1; ki<num_sample; ++ki, u+=udel)
	  {
	    u = std::min(u, u2);
	    polynomialTerms(u, v, degree1, degree2, tot_degree, &tmp[0]);
	    double val = 0.0;
	    for (int kb=0; kb<num_terms; ++kb)
	      val += pol[kb]*tmp[kb];

	    of << u << " " << v << " " << val << std::endl;

	    Point pos;
	    lrsf->point(pos, u, v);
	    double dd = fabs(val-pos[0]);
	    maxdist = std::max(maxdist, dd);
	    avdist += fac*dd;
	  }
      }
    std::cout << "Maxdist: " << maxdist << ", average: " << avdist << std::endl;
    
    lrsf->writeStandardHeader(sf_out);
    lrsf->write(sf_out);
}



     
