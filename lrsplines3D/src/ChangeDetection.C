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

#include "GoTools/lrsplines3D/ChangeDetection.h"
#include "GoTools/lrsplines3D/IdentifyChanges.h"
#include "GoTools/lrsplines3D/LRVolApprox.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRSurfUtils.h"

using namespace Go;
using std::vector;

//==============================================================================
ChangeDetection::ChangeDetection(vector<vector<double> >& point_clouds, int del)
//==============================================================================
  : del_(del), init_num_coef_(6), degree_(2), mba_level_(0.0), der_lim_(0.0)
{
  point_seqs_.resize(point_clouds.size());
  seq_box_.resize(point_clouds.size());
  Point low(del_), high(del_);
  double div = 1.0/(double)point_clouds.size();
  for (size_t ki=0; ki<point_clouds.size(); ++ki)
    {
      point_seqs_[ki].insert(point_seqs_[ki].begin(), point_clouds[ki].begin(),
			     point_clouds[ki].end());  // Copying
      seq_box_[ki].setFromArray(&point_seqs_[ki][0],
				&point_seqs_[ki][point_seqs_[ki].size()], del_);
      low += div*seq_box_[ki].low();
      high += div*seq_box_[ki].high();
    }
  av_box_ = BoundingBox(low, high);
}

//==============================================================================
double ChangeDetection::suggestTimeStep()
//==============================================================================
{
  if (del_ < 2)
    return -1;
  double tdel = 0.0;
  double area = 0.0;
  for (size_t ki=0; ki<point_seqs_.size(); ++ki)
    {
      size_t npt = point_seqs_[ki].size()/del_;
      double area0 = (seq_box_[ki].high()[0]-seq_box_[ki].low()[0])*
	(seq_box_[ki].high()[1] - seq_box_[ki].low()[1]);
      double tdel0 = sqrt(area0/(double)npt);
      area += area0;
      tdel += tdel0;
    }
  area /= (double)point_seqs_.size();
  double a2 = sqrt(area);

  double fac1 = 2.0;
  double fac2 = 100.0;
  tdel = std::min(tdel, fac1*a2);
  tdel = std::max(tdel, a2/fac2);
  
  tdel /= (double)point_seqs_.size();

  time_delta_ = tdel;
  return time_delta_;
}

//==============================================================================
void ChangeDetection::defineTrivariate(double time_step)
//==============================================================================
{
  time_delta_ = time_step;
  size_t all_num = 0;
  size_t ki;
  for (ki=0; ki<point_seqs_.size(); ++ki)
    all_num += point_seqs_[ki].size();

  all_num /= del_;
  tri_pts_.resize(all_num*(del_ + 1));

  time_.resize(point_seqs_.size());
  
  double tpar = 0.0;
  int kb=0;
  mba_level_ = 0.0;
  domain_[0] = domain_[2] = domain_[4] = std::numeric_limits<double>::max();
  domain_[1] = domain_[3] = domain_[5] = std::numeric_limits<double>::lowest();
  for (ki=0; ki<point_seqs_.size(); ++ki, tpar+=time_step)
    {
      time_[ki] = tpar;
      for (int kj=0; kj<(int)point_seqs_[ki].size(); kj+=del_)
	{
	  for (int ka=0; ka<del_-1; ++ka)
	    {
	      tri_pts_[kb++] = point_seqs_[ki][kj+ka];
	      domain_[2*ka] = std::min(domain_[2*ka], point_seqs_[ki][kj+ka]);
	      domain_[2*ka+1] = std::max(domain_[2*ka+1], point_seqs_[ki][kj+ka]);
	    }
	  tri_pts_[kb++] = tpar;
	  tri_pts_[kb++] = point_seqs_[ki][kj+del_-1];
	  domain_[4] = std::min(domain_[4], tpar);
	  domain_[5] = std::max(domain_[5], tpar);

	  mba_level_ += point_seqs_[ki][kj+del_-1]/(double)all_num;
	}
    }
  
}

//==============================================================================
void ChangeDetection::volApprox(double tol, int max_iter)
//==============================================================================
{
  // Define initial number of coefficients in all parameter directions
  int order = degree_ + 1;
  int nm = init_num_coef_*init_num_coef_*init_num_coef_;
  double dom = (domain_[1]-domain_[0])*(domain_[3]-domain_[2])*
    (domain_[5]-domain_[4]);
  double c1 = std::pow((double)nm/dom, 1.0/3.0);
  int nc[3];
  for (int kj=0; kj<3; ++kj)
    {
      double tmp = c1*(domain_[2*kj+1]-domain_[2*kj]);
      nc[kj] = (int)tmp;
	++nc[kj];
      nc[kj] = std::max(nc[kj], order);
    }
  nc[2] = std::max(nc[2], order+1);  // At least one inner knot in the
  // time direction
  std::cout << "Number of coefficients: " << nc[0] << ", " << nc[1] << ", " << nc[2] << std::endl;

  // Approximation engine
  int dim = del_ - 2;
  LRVolApprox vol_approx(nc[0], order, nc[1], order, nc[2], order,
  			 tri_pts_, dim, domain_, tol, mba_level_);
  vol_approx.setInitMBA(true);
  vol_approx.setVerbose(true);  // Currently

  // Approximate
  std::cout << "Starting approximation..." << std::endl;

  vol_ = vol_approx.getApproxVol(max_dist_, av_dist_all_, av_dist_out_,
				 num_out_, max_iter);
  std::cout << "Resulting number of mesh positions: " << vol_->mesh().numDistinctKnots(1) << ", " << vol_->mesh().numDistinctKnots(2) << ", " << vol_->mesh().numDistinctKnots(3) << ", " << std::endl;

  
}


//==============================================================================
void ChangeDetection::identifyChanges(double limit)
//==============================================================================
{
  der_lim_ = limit;
  IdentifyChanges change_engine(vol_.get(), tri_pts_, true);
  change_engine.getChangeDom(der_lim_, change_box_);

  // Extend boxes with two or less acquisitions in the time direction
       for (size_t ki=0; ki<change_box_.size(); ++ki)
	 {
	   Point low = change_box_[ki].low();
	   Point high = change_box_[ki].high();
	   double t1 = low[2];
	   double t2 = high[2];
	   size_t kr, kh;
	   for (kr=0; kr<time_.size() && time_[kr]<t1; ++kr);
	   for (kh=kr+1; kh<time_.size() && time_[kh]<t2; ++kh);
	   if (kh - kr <= 2)
	     {
	       if (kr > 0)
		 low[2] = time_[kr-1];
	       if (kh < time_.size())
		 high[2] = time_[kh];
	       change_box_[ki].unset();
	       change_box_[ki].setFromPoints(low, high);
	     }
	 }
  int stop_break = 1;
}

//==============================================================================
void ChangeDetection::extractChangeData()
//==============================================================================
{
  change_pts_.resize(change_box_.size());
  for (size_t ki=0; ki<change_box_.size(); ++ki)
    {
      double x1 = change_box_[ki].low()[0];
      double x2 = change_box_[ki].high()[0];
      double y1 = change_box_[ki].low()[1];
      double y2 = change_box_[ki].high()[1];
      double t1 = change_box_[ki].low()[2];
      double t2 = change_box_[ki].high()[2];
      for (size_t kj=0; kj<point_seqs_.size(); ++kj)
	{
	  if (time_[kj] < t1 || time_[kj] > t2)
	    continue;
	  vector<double> sub_data;
	  for (size_t kr=0; kr<point_seqs_[kj].size(); kr+=del_)
	    {
	      if (point_seqs_[kj][kr] >= x1 && point_seqs_[kj][kr] <= x2 &&
		  point_seqs_[kj][kr+1] >= y1 && point_seqs_[kj][kr+1] <= y2)
		sub_data.insert(sub_data.end(), &point_seqs_[kj][kr],
				&point_seqs_[kj][kr+del_]);
	    }
	  if (sub_data.size() > 0)
	    {
	      change_pts_[ki].subseq_pts_.push_back(sub_data);
	      change_pts_[ki].pts_seq_ix_.push_back((int)kj);
	    }
	}
    }
}

void ChangeDetection::surfApprox(double tol, int max_iter)
{
  int order = degree_ + 1;
  int nmb_coef = order;

  for (size_t ki=0; ki<change_pts_.size(); ++ki)
    {
      double domain[4];
      domain[0] = change_box_[ki].low()[0];
      domain[1] = change_box_[ki].high()[0];
      domain[2] = change_box_[ki].low()[1];
      domain[3] = change_box_[ki].high()[1];
      int nm = nmb_coef*nmb_coef;
      double dom = (domain[1]-domain[0])*(domain[3]-domain[2]);
      double c1 = std::pow((double)nm/dom, 1.0/2.0);
      int nc[2];
      for (int kb=0; kb<2; ++kb)
	{
	  double tmp = c1*(domain[2*kb+1]-domain[2*kb]);
	  nc[kb] = (int)tmp;
	  ++nc[kb];
	  nc[kb] = std::max(nc[kb], order);
	}

     for (size_t kj=0; kj<change_pts_[ki].subseq_pts_.size(); ++kj)
	{
	  LRSurfApprox approx(nc[0], order, nc[1], order,
			      change_pts_[ki].subseq_pts_[kj],
			      del_-2, 
			      domain, tol, false, 0.0, true, true);
	  approx.setUseMBA(true);
	  if (del_ == 3)
	    {
	      double zmin = std::numeric_limits<double>::max();
	      double zmax = std::numeric_limits<double>::lowest();
	      for (size_t kr=0; kr<change_pts_[ki].subseq_pts_[kj].size(); kr+=3)
		{
		  double zval = change_pts_[ki].subseq_pts_[kj][kr+2];
		  zmin = std::min(zmin, zval);
		  zmax = std::max(zmax, zval);
		}
	      double zfac = std::max(tol, 0.005*(zmax - zmin));
	      approx.addLowerConstraint(domain[4] - zfac);
	      approx.addUpperConstraint(domain[5] + zfac);
	      approx.setLocalConstraint(zfac);
	    }
	  // Refinement strategy: Full span, both parameter directions, no threshold
	  approx.setRefinementStrategy(1, 0, 0, -100, 0, 0);
	  approx.setVerbose(true);
	  
	  double maxdist, avdist_out, avdist_all; // will be set below
	  int nmb_out_eps;        // will be set below
	  shared_ptr<LRSplineSurface> surf;
	  try {
	    surf = approx.getApproxSurf(maxdist, avdist_all,
					avdist_out, nmb_out_eps, 
					max_iter);
	  }
	  catch (...)
	    {
	      std::cout << "ERROR: Surface approximation failed" << std::endl;
	      break;
	    }

	  change_pts_[ki].subseq_sfs_.push_back(surf);
	  change_pts_[ki].maxdist_.push_back(maxdist);
	  change_pts_[ki].avdist_.push_back(avdist_all);
	  change_pts_[ki].avdistout_.push_back(avdist_out);
	  change_pts_[ki].nmb_out_.push_back(nmb_out_eps);
	}
    }
}

void ChangeDetection::differenceSurfaces()
{
  for (size_t ki=0; ki<change_pts_.size(); ++ki)
    {
      LRSurfUtils::defineOnSameMesh(change_pts_[ki].subseq_sfs_);
      for (size_t kj=1; kj<change_pts_[ki].subseq_sfs_.size(); ++kj)
	{
	  shared_ptr<LRSplineSurface> diffsf(new LRSplineSurface(*change_pts_[ki].subseq_sfs_[kj-1]));
	  LRSplineSurface::BSplineMap::const_iterator it1 = change_pts_[ki].subseq_sfs_[kj-1]->basisFunctionsBegin();
	  LRSplineSurface::BSplineMap::const_iterator it2 = change_pts_[ki].subseq_sfs_[kj]->basisFunctionsBegin();
	  LRSplineSurface::BSplineMap::const_iterator it3 = diffsf->basisFunctionsBegin();
	  double min_diff = std::numeric_limits<double>::max();;
	  double max_diff = std::numeric_limits<double>::lowest();;
	  double av_diff = 0.0;
	  for (; it1 != change_pts_[ki].subseq_sfs_[kj-1]->basisFunctionsEnd();
	       ++it1, ++it2, ++it3) 
	    {
	      Point c1 = it1->second->Coef();
	      Point c2 = it2->second->Coef();
	      Point c3 = c2 - c1;
	      min_diff = std::min(min_diff, c3[0]);
	      max_diff = std::max(max_diff, c3[0]);
	      av_diff += c3[0];
	      diffsf->setCoef(c3, it3->second.get());
	    }
	  av_diff /= change_pts_[ki].subseq_sfs_[kj-1]->numBasisFunctions();
	  change_pts_[ki].subseq_diffsfs_.push_back(diffsf);
	  change_pts_[ki].min_diff_.push_back(min_diff);
	  change_pts_[ki].max_diff_.push_back(max_diff);
	  change_pts_[ki].av_diff_.push_back(av_diff);
	}
    }
  }
