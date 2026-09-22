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
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRProjection.h"
#include "GoTools/utils/QRFactorization.h"
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

void readParVal(std::string& pointfile, int num_sample, int degree1,
		int degree2, int tot_degree, vector<double>& pol,
		vector<double>& range, vector<double>& Aqr,
		vector<double>& bqr, vector<double>& Als,
		vector<double>& bls)
{
  int del = 5;
  int nmb_pts = 0;
  vector<double> data;
  vector<double> extent(2*del, 0.0);   // Limits for points in all coordinates
  if (pointfile != "no")
    {
      std::ifstream is(pointfile.c_str());
      FileUtils::readTxtPointFile(is, del, data, nmb_pts, extent);
    }

  std::cout << "nmb_pts: " << nmb_pts << std::endl;
  if (nmb_pts == 0 && num_sample <= 1)
    return;
  std::cout << "num_sample: " << num_sample << std::endl;
  
  // Evaluate
  int num_terms = (int)pol.size();
  int num_point = nmb_pts + num_sample*num_sample;
  Aqr.resize(num_point*num_terms);
  bqr.resize(num_point);
  Als.resize(num_terms*num_terms, 0.0);
  bls.resize(num_terms, 0.0);
  vector<double> term(num_terms, 0.0);
  int kr=0;
  
  std::ofstream of("scattered_pts.g2");
  of << "400 1 0 4 55 100 100 255" << std::endl;
  of << num_point << std::endl;
  if (nmb_pts > 0)
    {
      double facu = (range[1] - range[0])/(extent[1] - extent[0]);
      double facv = (range[3] - range[2])/(extent[3] - extent[2]);
      for (int ki=0; ki<nmb_pts; ++ki)
	{
	  double u = range[0] + (data[ki*del] - extent[0])*facu;
	  double v = range[2] + (data[ki*del+1] - extent[2])*facv;
	  polynomialTerms(u, v, degree1, degree2, tot_degree, &term[0]);
	  double val = 0.0;
	  for (int ka=0; ka<num_terms; ++ka)
	    {
	      val += pol[ka]*term[ka];
	      Aqr[ki*num_terms+ka] = term[ka];
	      for (int kb=0; kb<num_terms; ++kb)
		Als[kb*num_terms+ka] += term[ka]*term[kb];
	    }
	  bqr[kr++] = val;

	  for (int kb=0; kb<num_terms; ++kb)
	    bls[kb] += val*term[kb];

	  of << u << " " << v << " " << val << std::endl;
	}
    }
  
  std::cout << "file read "  << std::endl;
  if (num_sample > 1)
    {
      double udel = (range[1] - range[0])/(double)(num_sample - 1);
      double vdel = (range[3] - range[2])/(double)(num_sample - 1);

      double u=range[0], v=range[2];
      int ka, kb, kc;
      for (kb=0, kc=0; kb<num_sample; ++kb, v+=vdel)
	for (ka=0, u=range[0]; ka<num_sample; ++ka, ++kc, u+=udel)
	  {
	    polynomialTerms(u, v, degree1, degree2, tot_degree, &term[0]);
	    double val = 0.0;
	    for (int ki=0; ki<num_terms; ++ki)
	      {
		val += pol[ki]*term[ki];
		Aqr[(nmb_pts+kc)*num_terms+ki] = term[ki];
		for (int kj=0; kj<num_terms; ++kj)
		  Als[kj*num_terms+ki] += term[ki]*term[kj];
	      }
	    
	    bqr[kr++] = val;
	    for (int kj=0; kj<num_terms; ++kj)
	      bls[kj] += val*term[kj];
	    
	    of << u << " " << v << " " << val << std::endl;
	  }
    }
}

int main(int argc, char *argv[])
{
  if (argc < 9) {
    std::cout << "Parameters: lr func in (.g2), lr func out (.g2), input points (.txt), no. of regular points, num_sample, degree1, degree2, extended total degree (0/1/2 (Chebyshev polynomial)), polynomial factors"  << std::endl;
    exit(-1);
  }
  std::ifstream sf_in(argv[1]);
  std::ofstream sf_out(argv[2]);
  std::string pointfile(argv[3]);
  int num_reg = atoi(argv[4]);
  int num_sample = atoi(argv[5]);
  int degree1 = atoi(argv[6]);
  int degree2 = atoi(argv[7]);
  int total = atoi(argv[8]);
  int maxdeg = std::max(degree1,degree2);
  int mindeg = std::min(degree1,degree2);
  int tot_degree = (total) ? degree1+degree2 : maxdeg;
  int num_terms = (tot_degree == maxdeg) ?
    maxdeg*(mindeg+1) - mindeg*(mindeg-1)/2 + 1 : (degree1+1)*(degree2+1);
  std::cout << degree1 << " " << degree2 << " " << tot_degree << std::endl;
  std::cout << "Total degree: " << total << ", num terms: " << num_terms << std::endl;
  if (total != 2 && argc != num_terms + 9)
    {
      std::cout << "Expecting " << num_terms << " polynomial factors" << std::endl;
      exit(-1);
    }

  shared_ptr<ParamSurface> dummy;
  vector<double> range(4);
  range[0] = range[2] = -1;
  range[1] = range[3] = 1;
  
  ObjectHeader header;
  header.read(sf_in);
  vector<shared_ptr<LRSplineSurface> > lrsf(8);
  lrsf[0] = shared_ptr<LRSplineSurface>(new LRSplineSurface());
  lrsf[0]->read(sf_in);
  lrsf[0]->setParameterDomain(range[0], range[1], range[2], range[3]);
  
  for (int ka=1; ka<8; ++ka)
    lrsf[ka] = shared_ptr<LRSplineSurface>(lrsf[0]->clone());
  
  vector<double> pol(num_terms, 0.0);
  if (total == 2)
    Chebyshev(degree1, degree2, pol);
  else
    {
      for (int ka=0; ka<num_terms; ++ka)
	pol[ka] = atof(argv[9+ka]);
    }
  for (int ka=0; ka<num_terms; ++ka)
    std::cout << pol[ka] << " ";
  std::cout << std::endl;

  vector<double> Aqr, bqr, Als, bls;
  readParVal(pointfile, num_reg, degree1, degree2, tot_degree,
	     pol, range, Aqr, bqr, Als, bls);

  std::cout << "Scattered data sampled" << std::endl;
  vector<Point> par(num_terms);
  if (num_terms == 1)
    par[0] = Point(0.0, 0.0);
  else if (total)
    {
      double udel = (range[1] - range[0])/(double)degree1;
      double vdel = (range[3] - range[2])/(double)degree2;
      double u, v;
      int ki, kj, kr;
      for (kj=0, kr=0, v=range[2]; kj<=degree2; ++kj, v+=vdel)
	{
	  v = std::min(v, range[3]);
	  for (ki=0, u=range[0]; ki<=degree1; ++ki, ++kr, u+=udel)
	    {
	      u = std::min(u, range[1]);
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
  

  vector<double> ls(num_terms*num_terms, 0.0), rs(num_terms, 0.0);
  vector<int> piv(num_terms);
  for (int ka=0; ka<num_terms; ++ka)
    piv[ka] = ka;

  std::cout << "num_terms: " << num_terms << std::endl;
  
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

  std::cout << "Ready for QR" << std::endl;
  vector<double> Q, R, x;
  int num_pt = (int)bqr.size();
  int num_sf = (num_pt > 1) ? 8 : 4;
  if (num_pt > 1)
    {
      QRFactorization::QRDecomp(Aqr, num_terms, num_pt, Q, R);
      QRFactorization::QRSolve(Q, R, num_terms, num_pt, bqr, 1, x);
  
      std::cout << "Ready for ls" << std::endl;
      for (int ka=0; ka<num_terms; ++ka)
	piv[ka] = ka;
      s6lufacp(&Als[0], &piv[0], num_terms, &kstat);
      if (kstat < 0)
	return false;

      s6lusolp(&Als[0], &bls[0], &piv[0], num_terms, &kstat);
      if (kstat < 0)
	return false;
    }
  
  std::cout << std::endl << "Interpolation point distances: " << std::endl;
  vector<double> tmp(num_terms);
  double max_int = 0.0, av_int = 0.0, max_qr = 0.0, av_qr = 0.0, max_ls = 0.0, av_ls = 0.0;
  double fac_term = 1.0/(double)num_terms;
  for (int ka=0; ka<num_terms; ++ka)
    {
      polynomialTerms(par[ka][0], par[ka][1], degree1, degree2,
		      tot_degree, &tmp[0]);
      double val1 = 0.0, val2 = 0.0, val3 = 0.0, val4 = 0.0;
      for (int kb=0; kb<num_terms; ++kb)
	{
	  val1 += pol[kb]*tmp[kb];
	  val2 += rs[kb]*tmp[kb];
	  if (num_sf > 2)
	    {
	      val3 += x[kb]*tmp[kb];
	      val4 += bls[kb]*tmp[kb];
	    }
	}
      max_int = std::max(max_int, fabs(val2-val1));
      av_int += fac_term*fabs(val2-val1);
      //std::cout << "Interpolation:" << val1 << " " << val2 << " " << val1-val2 << std::endl;
	  if (num_sf > 2)
	    {
	      max_qr = std::max(max_qr, fabs(val3-val1));
	      av_qr += fac_term*fabs(val3-val1);
	      max_ls = std::max(max_ls, fabs(val4-val1));
	      av_ls += fac_term*fabs(val4-val1);
	      
	      //std::cout << "QR: " << val1 << " " << val3 << " " << val1-val3 << std::endl;
	      //std::cout << "ls: " << val1 << " " << val4 << " " << val1-val4 << std::endl;
	    }
    }
  std::cout << "Interpolation, maxdist: " << max_int << ", average: " << av_int << std::endl;
  std::cout << "QR, maxdist: " << max_qr << ", average: " << av_qr << std::endl;
  std::cout << "ls, maxdist: " << max_ls << ", average: " << av_ls << std::endl;
  
  std::cout << std::endl << "Differences polynom: " << std::endl;
  max_int = av_int = max_qr = av_qr = max_ls = av_ls = 0.0;
  for (int ka=0; ka<num_terms; ++ka)
    {
      double diff1 = rs[ka] - pol[ka];
      max_int = std::max(max_int, fabs(diff1));
      av_int += fac_term*fabs(diff1);
      //std::cout << "Interpolation:" << pol[ka] << " " << rs[ka] << " " << diff1 << std::endl;
	  if (num_sf > 2)
	    {
	      double diff2 = x[ka] - pol[ka];
	      double diff3 = bls[ka] - pol[ka];
	      max_qr = std::max(max_qr, fabs(diff2));
	      av_qr += fac_term*fabs(diff2);
	      max_ls = std::max(max_ls, fabs(diff3));
	      av_ls += fac_term*fabs(diff3);
	      //std::cout << "QR: " << pol[ka] << " " << x[ka] << " " << diff2 << std::endl;
	      //std::cout << "ls: " << pol[ka] << " " << bls[ka] << " " << diff3 << std::endl;
	    }
    }

  std::cout << "Interpolation, maxdist: " << max_int << ", average: " << av_int << std::endl;
  std::cout << "QR, maxdist: " << max_qr << ", average: " << av_qr << std::endl;
  std::cout << "ls, maxdist: " << max_ls << ", average: " << av_ls << std::endl << std::endl;
  
    int max_level = 10;  // Should always be enough
    for (int ka=0; ka<num_sf/2; ++ka)
      {
	vector<double> curr_pol;
	if (ka == 3)
	  curr_pol = bls;
	else if (ka == 2)
	  curr_pol = x;
	else if (ka == 1)
	  curr_pol = rs;
	else
	  curr_pol = pol;
	for (int kb=0; kb<2; ++kb)
	  {
	    for (int level=0; level<max_level; ++level)
	      {
		int num_update = 0;
		for (auto bspl=lrsf[2*ka+kb]->basisFunctionsBegin();
		     bspl!=lrsf[2*ka+kb]->basisFunctionsEnd(); ++bspl)
		  {
		    int blevel = bspl->second->getNestLevel();
		    if (blevel != level)
		      continue;

		    Point coef;
		    if (kb == 0)
		      coef =
		      LRProjection::Polynomial2Coef(curr_pol, degree1, degree2,
						    tot_degree,
						    range[0], range[1],
						    range[2], range[3],
						    bspl->second.get(), dummy);
		    else
		      coef =
			LRProjection::Polynomial2Coef2(curr_pol, degree1,
						       degree2, tot_degree,
						       bspl->second.get(), dummy);
		    if (blevel > 0)
		      {
			std::cout << "Level: " << blevel << std::endl;
		
			bool OK = bspl->second->adaptProjCoef(coef);
		      }
		    lrsf[2*ka+kb]->setCoef(coef, bspl->second.get());
		    num_update++;
		  }
		if (num_update == 0)
		  break;
	      }
	  }
      }
    std::ofstream of("sample_pts.g2");
    of.precision(15);
    double udel = (range[1] - range[0])/(double)(num_sample-1);
    double vdel = (range[3] - range[2])/(double)(num_sample-1);
    double u, v;
    int ki, kj;
    double fac = 1.0/(double)(num_sample*num_sample);
    for (int ka=0; ka<num_sf/2; ++ka)
      {
	of << "400 1 0 4 255 0 0 255" << std::endl;
	of << num_sample*num_sample << std::endl;
	double maxdist = 0.0, avdist = 0.0;
	double maxdist2 = 0.0, avdist2 = 0.0;
	vector<double> curr_pol;
	if (ka == 3)
	  curr_pol = bls;
	else if (ka == 2)
	  curr_pol = x;
	else if (ka == 1)
	  curr_pol = rs;
	else
	  curr_pol = pol;
	for (int kc=0; kc<2; ++kc)
	  {
	    for (kj=0, v=range[2]; kj<num_sample; ++kj, v+=vdel)
	      {
		v = std::min(v, range[3]);
		for (ki=0, u=range[0]; ki<num_sample; ++ki, u+=udel)
		  {
		    u = std::min(u, range[1]);
		    polynomialTerms(u, v, degree1, degree2, tot_degree, &tmp[0]);
		    double val = 0.0, val2 = 0.0;
		    for (int kb=0; kb<num_terms; ++kb)
		      {
			val += pol[kb]*tmp[kb];
			val2 += curr_pol[kb]*tmp[kb];
		      }
		    of << u << " " << v << " " << val << std::endl;

		    Point pos;
		    lrsf[2*ka+kc]->point(pos, u, v);
		    double dd = fabs(val-pos[0]);
		    double dd2 = fabs(val2-pos[0]);
		    //std::cout << u << " " << v << " " << val << " " << pos[0] << " " << dd << std::endl;
		    maxdist = std::max(maxdist, dd);
		    avdist += fac*dd;
		    maxdist2 = std::max(maxdist2, dd2);
		    avdist2 += fac*dd2;
		  }
	      }
	  
	    std::cout << "Surf: " << ka << ", maxdist: " << maxdist << ", average: " << avdist << std::endl;
    
	    std::cout << "maxdist2: " << maxdist2 << ", average2: " << avdist2 << std::endl;
	    lrsf[ka]->writeStandardHeader(sf_out);
	    lrsf[ka]->write(sf_out);
	  }
      }
}



     
