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

#include "GoTools/lrsplines3D/IdentifyChanges.h"
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/lrsplines3D/LRSpline3DEvalGrid.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/creators/ApproxCurve.h"
#include "sislP.h"
#include <fstream>

using namespace Go;
using std::vector;

#define MAX_COLORS 12
int colours[MAX_COLORS][3] = {
  {255, 0, 0},
  {0, 255, 0},
  {0, 0, 255},
  {255, 255, 0},
  {255, 0, 255},
  {0, 255, 255},
  {128, 255, 0},
  {255, 128, 0},
  {128, 0, 255},
  {255, 0, 128},
  {0, 128, 255},
  {0, 255, 128},
};

// =============================================================================
IdentifyChanges::IdentifyChanges(LRSplineVolume* vol, vector<double>& points,
				 bool distribute)
// =============================================================================
  : points_(points), cdir_(0), min_der_(0.0), max_der_(0.0), av_der_(0.0)
{
  lrvol_ = shared_ptr<LRSplineVolume>(vol->clone());
  order_u_ = lrvol_->degree(XDIR) + 1;
  order_v_ = lrvol_->degree(YDIR) + 1;
  order_w_ = lrvol_->degree(ZDIR) + 1;
  dim_ = lrvol_->dimension();
  const Mesh3D mesh = lrvol_->mesh();
  knots_u_.insert(knots_u_.end(), mesh.knotsBegin(XDIR), mesh.knotsEnd(XDIR)); 
  knots_v_.insert(knots_v_.end(), mesh.knotsBegin(YDIR), mesh.knotsEnd(YDIR)); 
  knots_w_.insert(knots_w_.end(), mesh.knotsBegin(ZDIR), mesh.knotsEnd(ZDIR)); 

  nmb_pts_ = points_.size()/dim_;

  // Distribute points to volume
  if (distribute)
    LRSpline3DUtils::distributeDataPoints(lrvol_.get(), points_,
					  false, true, true); 

}

// =============================================================================
void IdentifyChanges::getChangeDom(double limit,
				   vector<BoundingBox>& tri_dom)
// =============================================================================
{
  computeElemInfo();

  std::ofstream of("eleminfo.g2");
  of << "410 1 0 4 0 125 120 255" << std::endl;
  of << 12*eleminfo_.size() << std::endl;
  for (size_t ki=0; ki<eleminfo_.size(); ++ki)
    {
      Point low = eleminfo_[ki].tri_box_.low();
      Point high = eleminfo_[ki].tri_box_.high();
      Point pt2(high[0], low[1], low[2]);
      Point pt3(low[0], high[1], low[2]);
      Point pt4(high[0], high[1], low[2]);
      Point pt5(low[0], low[1], high[2]);
      Point pt6(high[0], low[1], high[2]);
      Point pt7(low[0], high[1], high[2]);
      of << low << " " << pt2 << std::endl;
      of << low << " " << pt3 << std::endl;
      of << pt2 << " " << pt4 << std::endl;
      of << pt3 << " " << pt4 << std::endl;
      of << high << " " << pt6 << std::endl;
      of << high << " " << pt7 << std::endl;
      of << pt6 << " " << pt5 << std::endl;
      of << pt7 << " " << pt5 << std::endl;
      of << low << " " << pt5 << std::endl;
      of << pt2 << " " << pt6 << std::endl;
      of << pt3 << " " << pt7 << std::endl;
      of << pt4 << " " << high << std::endl;
    }

  
  double limit2 = limit; //fabs(av_der_);
  // if (limit2 < 1.0e-4)
  //   limit2 = 1.0e-4;

  vector<double> maxder(eleminfo_.size());
  vector<double> parder(eleminfo_.size());
  vector<int> type(eleminfo_.size(), 1);
  for (size_t ki=0; ki<eleminfo_.size(); ++ki)
    {
      maxder[ki] = std::max(fabs(eleminfo_[ki].max_val_),
			    fabs(eleminfo_[ki].min_val_));
      parder[ki] = 0.01*ki;
    }

  std::sort(maxder.begin(), maxder.end());
  vector<double> diff(eleminfo_.size()-1);
  for (size_t ki=1; ki<eleminfo_.size(); ++ki)
    diff[ki-1] = maxder[ki] - maxder[ki-1];

  vector<double> diff2(diff.begin(), diff.end());
  std::sort(diff2.begin(), diff2.end());

  double eps = 0.1*std::max(fabs(max_der_), fabs(min_der_));
  ApproxCurve approx(maxder, parder, 1, eps, 3, 3);
  approx.setSmooth(0.1);
  double maxd, avd;
  shared_ptr<SplineCurve> cv = approx.getApproxCurve(maxd, avd);
  shared_ptr<SplineCurve> cv2(cv->derivCurve(1));
  int in = cv2->numCoefs();
  int ik = cv2->order();
  std::ofstream of0("dercv0.g2");
  of0 << "100 1 0 0" << std::endl;
  of0 << "3 0" << std::endl;
  of0 << in << " " << ik << std::endl;
  for (auto it=cv2->basis().begin(); it!=cv2->basis().end(); ++it)
    of0 << *it << " ";
  of0 << std::endl;
  int ka = 0;
  for (auto it2=cv2->coefs_begin(); ka<in; ++it2, ++ka)
    of0 << cv2->basis().grevilleParameter(ka) << " " << *it2 << " 0.0" << std::endl;

  SISLCurve *qc = 0;
  double *gpar = 0;
  double endpar;
  int npar;
  int kstat = 0;
  s1357(&maxder[0], (int)maxder.size(), 1, &type[0], &parder[0],
	0, 0, 1, 3, parder[0], &endpar, &qc, &gpar, &npar, &kstat);
  if (kstat >= 0 && qc != 0)
    {
      SISLCurve *qc2 = 0;
      s1720(qc, 1, &qc2, &kstat);
      
      std::ofstream of("dercv.g2");
      of << "100 1 0 0" << std::endl;
      of << "3 0" << std::endl;
      of << qc2->in << " " << qc2->ik << std::endl;
      for (int ka=0; ka<qc2->in+qc2->ik; ++ka)
	of << qc2->et[ka] << " ";
      of << std::endl;
      for (int ka=0; ka<qc2->in; ++ka)
	{
	  double tmp = 0.0;
	  for (int kb=1; kb<qc2->ik; ++kb)
	    tmp += qc2->et[ka+kb];
	  tmp /= (double)(qc2->ik-1);
	  of << tmp << " " << qc2->ecoef[ka] << " 0.0" << std::endl;
	}
      // shared_ptr<SplineCurve> dercv(SISLCurve2Go(qc));
      // dercv->writeStandardHeader(of);
      // dercv->write(of);
      if (qc2) freeCurve(qc2);
    }
  if (qc) freeCurve(qc);
  if (gpar) free(gpar);
  
  setCellFlags(limit2);

  analyzeCellFlag();

  std::ofstream of2("change_box.g2");
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      Point low(knots_u_[change_dom_[ki].x1_],
		knots_v_[change_dom_[ki].y1_],
		knots_w_[change_dom_[ki].z1_]);
      Point high(knots_u_[change_dom_[ki].x2_+1],
		 knots_v_[change_dom_[ki].y2_+1],
		 knots_w_[change_dom_[ki].z2_+1]);
      tri_dom.push_back(BoundingBox(low, high));

      of2 << "410 1 0 4 0 0 0 255" << std::endl;
      of2 << "12" << std::endl;
      Point pt2(high[0], low[1], low[2]);
      Point pt3(low[0], high[1], low[2]);
      Point pt4(high[0], high[1], low[2]);
      Point pt5(low[0], low[1], high[2]);
      Point pt6(high[0], low[1], high[2]);
      Point pt7(low[0], high[1], high[2]);
      of2 << low << " " << pt2 << std::endl;
      of2 << low << " " << pt3 << std::endl;
      of2 << pt2 << " " << pt4 << std::endl;
      of2 << pt3 << " " << pt4 << std::endl;
      of2 << high << " " << pt6 << std::endl;
      of2 << high << " " << pt7 << std::endl;
      of2 << pt6 << " " << pt5 << std::endl;
      of2 << pt7 << " " << pt5 << std::endl;
      of2 << low << " " << pt5 << std::endl;
      of2 << pt2 << " " << pt6 << std::endl;
      of2 << pt3 << " " << pt7 << std::endl;
      of2 << pt4 << " " << high << std::endl;
    }
  int stop_break = 1;
}

// =============================================================================
void IdentifyChanges::computeElemInfo()
// =============================================================================
{
  // Compute derivative in time direction defined as Bezier patches
  computeTimeDeriv();
  
  // Traverse the Bezier representation of the time derivatires and
  // identify patches with high variability that are candidate for
  // split
  const Mesh3D& mesh = lrvol_->mesh();
  vector<LRSplineVolume::Refinement3D> mesh_refs;
  min_der_ =  std::numeric_limits<double>::max();
  max_der_ = std::numeric_limits<double>::lowest();
  av_der_ = 0.0;
  double fac = 1.0/(double)eleminfo_.size();
  int xnmb = order_u_;
  int ynmb = order_v_;
  int tnmb = order_w_;  for (size_t kj=0; kj<eleminfo_.size(); ++kj)
    {
      if (eleminfo_[kj].bez2_coefs_.size() == 0)
	continue;   // Element without points

      int ncoef = (int)eleminfo_[kj].bez2_coefs_.size()/dim_;
      double min_val =  std::numeric_limits<double>::max();
      double max_val = std::numeric_limits<double>::lowest();
      double av_val = 0.0;
      int kr, kx, ky, kt;
      for (kr=0, kt=0; kt<tnmb; ++kt)
	for (ky=0; ky<ynmb; ++ky)
	  for (kx=0; kx<xnmb; ++kx, kr+=dim_)
	    {
	      double curr_coef = eleminfo_[kj].bez2_coefs_[kr+cdir_];
	      min_val = std::min(min_val, curr_coef);
	      max_val = std::max(max_val, curr_coef);
	      av_val += curr_coef;
	      if (curr_coef < 0.0)
		eleminfo_[kj].acc_neg_ += curr_coef;
	      else
		eleminfo_[kj].acc_pos_ += curr_coef;

	      // Check sign of neighbours
	      int nn = 0;
	      int nsame1 = 0, nsame2 = 0;
	      int i1, i2, i3;
	      for (i3=std::max(0,kt-1); i3<std::min(kt+2,tnmb); ++i3)
		for (i2=std::max(0,ky-1); i2<std::min(ky+2,ynmb); ++i2)
		  for (i1=std::max(0,kx-1); i1<std::min(kx+2,xnmb); ++i1)
		    {
		      if (i1==kx && i2==ky && i3==kt)
			continue;
		      ++nn;
		      double next_coef =
			eleminfo_[kj].bez2_coefs_[(i3*ynmb+i2)*xnmb*dim_+cdir_];
		      if (curr_coef*next_coef > 0.0)
			{
			  if (curr_coef > 0.0)
			    ++nsame1;
			  else
			    ++nsame2;
			}
		    }
	      if (nsame1 == nn)
		eleminfo_[kj].same_pos_++;
	      else if (nsame2 == nn)
		eleminfo_[kj].same_neg_++;
	    }
      av_val /= (double)ncoef;

      eleminfo_[kj].min_val_ = min_val;
      eleminfo_[kj].max_val_ = max_val;
      eleminfo_[kj].av_val_ = av_val;
      min_der_ = std::min(min_der_, min_val);
      max_der_ = std::max(max_der_, max_val);
      av_der_ += fac*av_val;

      // Count number of coefficients with the same sign from start and end
      int nstart=1, nend=1;
      double cfstart = eleminfo_[kj].bez2_coefs_[cdir_];
      double cfend = eleminfo_[kj].bez2_coefs_[eleminfo_[kj].bez2_coefs_.size() - dim_ + cdir_];
      if (cfstart*cfend < 0.0)
	{
	  for (int ka=dim_; ka<(int)eleminfo_[kj].bez2_coefs_.size(); ka+=dim_)
	    {
	      double curr_coef = eleminfo_[kj].bez2_coefs_[ka+cdir_];
	      if ((cfstart > 0.0 && curr_coef > 0.0) ||
		  (cfstart < 0.0 && curr_coef < 0.0))
		nstart++;
	      else
		break;
	    }
	  for (int ka=(int)eleminfo_[kj].bez2_coefs_.size()-2*dim_; ka>=0; ka-=dim_)
	    {
	      double curr_coef = eleminfo_[kj].bez2_coefs_[ka+cdir_];
	      if ((cfend > 0.0 && curr_coef > 0.0) ||
		  (cfend < 0.0 && curr_coef < 0.0))
		nend++;
	      else
		break;
	    }
	}

      if (nstart >= xnmb*ynmb && nend >= xnmb*ynmb)
	eleminfo_[kj].time_change_ = true;
    }
}

// =============================================================================
void IdentifyChanges::setCellFlags(double limit)
// =============================================================================
{
  double eps = 1.0e-8;
  double limit2 = 0.5*limit;
  double frac = 0.5;

  const Mesh3D& mesh = lrvol_->mesh();
  vector<LRSplineVolume::Refinement3D> mesh_refs;
  for (size_t kj=0; kj<eleminfo_.size(); ++kj)
    {
     if (std::max(fabs(eleminfo_[kj].min_val_),
		  fabs(eleminfo_[kj].max_val_)) < limit)
	continue;  // No change identified

      if (eleminfo_[kj].min_val_*eleminfo_[kj].max_val_ >= 0.0)
	continue;   // Consistent direction

      if (std::min(fabs(eleminfo_[kj].min_val_),
		   fabs(eleminfo_[kj].max_val_)) < 0.1*limit)
	continue;  // One dominant direction
      
      int ka;
      for (ka=0; ka<3; ++ka)
	{
	  int ix1 = mesh.getKnotIdx(ka, eleminfo_[kj].tri_box_.low()[ka], eps);
	  int ix2 = mesh.getKnotIdx(ka, eleminfo_[kj].tri_box_.high()[ka], eps);
	  if (ix2 - ix1 > 1)
	    break;
	}
      if (ka == 3)
	continue;   // No internal candidate mesh line

      int stop_break = 1;
    }

  if (mesh_refs.size() > 0)
    {
      // Refine volume and recompute element information
    }

  vector<vector<Point> > midpt(10);
  int num_cell1 = mesh.numDistinctKnots(XDIR) - 1;
  int num_cell2 = mesh.numDistinctKnots(YDIR) - 1;
  int num_cell3 = mesh.numDistinctKnots(ZDIR) - 1;
  cell_flag_.resize(num_cell1*num_cell2*num_cell3);
  std::fill(cell_flag_.begin(), cell_flag_.end(), FLAT);
  
  for (size_t kj=0; kj<eleminfo_.size(); ++kj)
    {
      if (eleminfo_[kj].bez2_coefs_.size() == 0)
	continue;   // Element without points

      int flag = SCATTERED;
      int size = eleminfo_[kj].bez2_coefs_.size();
      int sgn1 = (eleminfo_[kj].bez2_coefs_[cdir_] < 0.0) ? -1 : 1;
      int sgn2 = (eleminfo_[kj].bez2_coefs_[size - dim_ + cdir_] < 0.0) ? -1 : 1;
      bool pos = (eleminfo_[kj].same_pos_ >= order_u_*order_v_);
      bool neg = (eleminfo_[kj].same_neg_ >= order_u_*order_v_);
      double lim = (pos || neg) ? limit2 : limit;
      if (eleminfo_[kj].time_change_)
	{
	  int first = FLAT, second = FLAT;
	  if (sgn1 > 0 && (eleminfo_[kj].max_val_ > lim ||
			   frac*eleminfo_[kj].max_val_ > fabs(eleminfo_[kj].min_val_)))
	    first = UP;
	  else if (sgn1 < 0 && (eleminfo_[kj].min_val_ < -lim ||
				frac*fabs(eleminfo_[kj].min_val_) > eleminfo_[kj].max_val_))
	    first = DOWN;
	  if (sgn2 > 0 && (eleminfo_[kj].max_val_ > lim ||
			   frac*eleminfo_[kj].max_val_ > fabs(eleminfo_[kj].min_val_)))
	    second = UP;
	  else if (sgn2 < 0 && (eleminfo_[kj].min_val_ < -lim ||
				frac*fabs(eleminfo_[kj].min_val_) > eleminfo_[kj].max_val_))
	    second = DOWN;
	  if (first == UP)
	    flag = (second == UP) ? UP : ((second == FLAT) ? UP_FLAT : UP_DOWN);
	  else if (first == FLAT)
	    flag = (second == UP) ? FLAT_UP : ((second == FLAT) ? FLAT : FLAT_DOWN);
	  else if (first == DOWN)
	    flag = (second == UP) ? DOWN_UP : ((second == FLAT) ? DOWN_FLAT : DOWN);
	}
      else
	{
	  if (std::max(fabs(eleminfo_[kj].min_val_),
		       fabs(eleminfo_[kj].max_val_)) < lim)
	    continue;  // No change identified

	  // Define flag
	  if (fabs(eleminfo_[kj].min_val_) < lim ||
	      eleminfo_[kj].min_val_ > 0.0 /*|| eleminfo_[kj].av_val_ > 0.1*limit*/)
	    flag = UP;
	  else if (fabs(eleminfo_[kj].max_val_) < lim ||
		   eleminfo_[kj].max_val_ < 0.0 /*|| eleminfo_[kj].av_val_ < -0.1*limit*/)
	    flag = DOWN;
	}

      // Find associated cells
      int ix1[3], ix2[3];
      for (int ka=0; ka<3; ++ka)
	{
	  ix1[ka] = mesh.getKnotIdx(ka, eleminfo_[kj].tri_box_.low()[ka], eps);
	  ix2[ka] = mesh.getKnotIdx(ka, eleminfo_[kj].tri_box_.high()[ka], eps);
	}

      // if (flag == 10)
      // 	flag = 0;  // TEST
      int xydel = num_cell1*num_cell2;
      for (int kc=ix1[2]; kc<ix2[2]; ++kc)
	for (int kb=ix1[1]; kb<ix2[1]; ++kb)
	  for (int ka=ix1[0]; ka<ix2[0]; ++ka)
	    {
	      cell_flag_[kc*xydel+kb*num_cell1+ka] = flag;
	      Point mid(0.5*(mesh.kval(1,ka)+mesh.kval(1,ka+1)),
			0.5*(mesh.kval(2,kb)+mesh.kval(2,kb+1)),
			0.5*(mesh.kval(3,kc)+mesh.kval(3,kc+1)));
	      midpt[flag].push_back(mid);
	    }
  
      int stop_break = 1;
    }

  std::ofstream of("elemmid.g2");
  for (int ka=1; ka<9; ++ka)
    {
      if (midpt[ka].size() > 0)
	{
	  of << "400 1 0 4 " << colours[ka][0] << " " << colours[ka][1] << " " << colours[ka][2] << " 255" << std::endl;
	  of << midpt[ka].size() << std::endl;
	  for (size_t ki=0; ki<midpt[ka].size(); ++ki)
	    of << midpt[ka][ki] << std::endl;
	}
    }

  splitCells();

  int stop_break2 = 1;
 }

// =============================================================================
void IdentifyChanges::splitCells()
// =============================================================================
{
  size_t no_level = (knots_u_.size()-1)*(knots_v_.size()-1);
  for (size_t ki=0; ki<knots_w_.size()-1; ++ki)
    {
      // Check if the cell layer should be split
      size_t kj;
      for (kj=0; kj<no_level; ++kj)
	{
	  int flag = cell_flag_[ki*no_level+kj];
	  if (!(flag <= UP || flag == SCATTERED))
	    break;
	}
      if (kj < no_level)
	{
	  // Split cell layer. Define new knot value in z
	  double zknot = 0.5*(knots_w_[ki] + knots_w_[ki+1]);
	  knots_w_.insert(knots_w_.begin()+ki+1, zknot);

	  vector<double> new_flag(no_level, FLAT);
	  for (size_t kr=0; kr<no_level; ++kr)
	    {
	      int flag = cell_flag_[ki*no_level+kr];
	      if (flag == DOWN_FLAT || flag == DOWN_UP)
		cell_flag_[ki*no_level+kr] = DOWN;
	      else if (flag == FLAT_DOWN || flag == FLAT_UP)
		cell_flag_[ki*no_level+kr] = FLAT;
	      else if (flag == UP_DOWN || flag == UP_FLAT)
		cell_flag_[ki*no_level+kr] = UP;

	      if (flag == FLAT_UP || flag == DOWN_UP)
		new_flag[kr] = UP;
	      else if (flag == FLAT_DOWN || flag == UP_DOWN)
		new_flag[kr] = DOWN;
	      else if (flag == DOWN_FLAT || flag == UP_FLAT)
		new_flag[kr] = FLAT;
	      else
		new_flag[kr] = flag;
	    }
	  cell_flag_.insert(cell_flag_.begin()+(ki+1)*no_level,
			    new_flag.begin(), new_flag.end());
	  ++ki;
	}
    }

  vector<vector<Point> > midpt(10);
  int num_cell1 = (int)knots_u_.size() - 1;
  int num_cell2 = (int)knots_v_.size() - 1;
  int num_cell3 = (int)knots_w_.size() - 1;
  for (int kz=0, ka=0; kz<num_cell3; ++kz)
    {
       for (int ky=0; ky<num_cell2; ++ky)
	{
	  for (int kx=0; kx<num_cell1; ++kx, ++ka)
	    {
	      int flag = cell_flag_[ka];
	      Point mid(0.5*(knots_u_[kx]+knots_u_[kx+1]),
			0.5*(knots_v_[ky]+knots_v_[ky+1]),
			0.5*(knots_w_[kz]+knots_w_[kz+1]));
	      midpt[flag].push_back(mid);
	    }
	}
    }
 
  std::ofstream of("cellmid.g2");
  for (int ka=1; ka<9; ++ka)
    {
      if (midpt[ka].size() > 0)
	{
	  of << "400 1 0 4 " << colours[ka][0] << " " << colours[ka][1] << " " << colours[ka][2] << " 255" << std::endl;
	  of << midpt[ka].size() << std::endl;
	  for (size_t ki=0; ki<midpt[ka].size(); ++ki)
	    of << midpt[ka][ki] << std::endl;
	}
    }
  int stop_break = 1;
}

// =============================================================================
void IdentifyChanges::analyzeCellFlag()
// =============================================================================
{
  int num_cell1 = (int)knots_u_.size() - 1;
  int num_cell2 = (int)knots_v_.size() - 1;
  int num_cell3 = (int)knots_w_.size() - 1;

  vector<vector<ChangeRectangle> > rectxy(num_cell3);
  for (int kz=0; kz<num_cell3; ++kz)
    {
      // Intervals in x
      vector<vector<ChangeRectangle > > rectx(num_cell2);
      for (int ky=0; ky<num_cell2; ++ky)
	{
	  int kx=0, kx2=0;
	  int kstart = (num_cell2*kz + ky)*num_cell1;
	  for (kx=0; kx<num_cell1; kx=kx2)
	    {
	      for (; kx<num_cell1 && cell_flag_[kstart+kx] == FLAT; ++kx);
	      int flag = cell_flag_[kstart+kx];
	      int no_pos = (flag == UP);
	      int no_neg = (flag == DOWN);
	      int no_scattered = (flag == SCATTERED);
	      for (kx2=kx+1; kx2<num_cell1 &&
		     (cell_flag_[kstart+kx2] == cell_flag_[kstart+kx] ||
		      cell_flag_[kstart+kx2] == SCATTERED); ++kx2)
		{
		  if (cell_flag_[kstart+kx2] == UP)
		    no_pos++;
		  else if (cell_flag_[kstart+kx2] == DOWN)
		    no_neg++;
		  else
		    no_scattered++;
		}
	      if (kx < num_cell1)
		rectx[ky].push_back(ChangeRectangle(kx, kx2-1, ky, ky, kz, kz,
						    no_pos, no_neg, no_scattered));
	    
	  // int no_pos2=0, no_neg2=0, no_scattered2 = 0;
	  // int no_x = (int)knots_u_.size() - 1;
	  // int no_y = (int)knots_v_.size() - 1;
	  // for (int ka=kx; ka<kx2; ++ka)
	  // {
	  //   int flag = cell_flag_[(kz*no_y+ky)*no_x+ka];
	  //   if (flag == UP)
	  //     no_pos2++;
	  //   else if (flag == DOWN)
	  //     no_neg2++;
	  //   else if (flag == SCATTERED)
	  //     no_scattered2++;
	  // }
	  // if (no_pos2 != no_pos || no_neg2 != no_neg ||
	  //     no_scattered2 != no_scattered)
	  //   std::cout << "Sign mismatchx: " << kz << " " << ky << std::endl;
	    }
	  // Merge across small breaks
	  bool domerge = true; //false;
	  if (domerge)
	    {
	  int stop_break1 = 1;
	  for (size_t ki=1; ki<rectx[ky].size(); )
	    {
	      int flag1 = rectx[ky][ki-1].getMainFlag();
	      int flag2 = rectx[ky][ki].getMainFlag();
	      if (flag1 != flag2 && !(flag1 == SCATTERED || flag2 == SCATTERED))
		{
		  ++ki;
		  continue;
		}
	      int len1 = rectx[ky][ki-1].x2_ - rectx[ky][ki-1].x1_ + 1;
	      int len2 = rectx[ky][ki].x2_ - rectx[ky][ki].x1_ + 1;
	      int lenm = rectx[ky][ki].x1_ - rectx[ky][ki-1].x2_ - 1;
	      int lenlim = std::min(num_cell1/10, std::max(len1,len2)/2);
	      if (lenm < lenlim)
		{
		  rectx[ky][ki-1].x2_ = rectx[ky][ki].x2_;
		  rectx[ky][ki-1].no_pos_ += rectx[ky][ki].no_pos_;
		  rectx[ky][ki-1].no_neg_ += rectx[ky][ki].no_neg_;
		  rectx[ky][ki-1].no_scattered_ += rectx[ky][ki].no_scattered_;
		  rectx[ky].erase(rectx[ky].begin()+ki);
		}
	      else
		++ki;
	    }
	    }
	  int stop_break3 = 1;
	}

      // Rectangles in xy
      for (int ky=0; ky<num_cell2; ++ky)
	{
	  if (rectx[ky].size() == 0)
	    continue;

	  for (size_t ki=0; ki<rectx[ky].size(); ++ki)
	    {
	      // Check if the interval is part of an existing (possibly
	      // modified change rectangle
	      size_t kj;
	      for (kj=0; kj<rectxy[kz].size(); ++kj)
		{
		  if (rectxy[kz][kj].y2_ < ky-1)
		    continue;  // Not continous

		  // Check flag
		  int flag1 = rectx[ky][ki].getMainFlag();
		  int flag2 = rectxy[kz][kj].getMainFlag();
		  if (flag1 != flag2 && !(flag1 == SCATTERED || flag2 == SCATTERED))
		    continue;
		  
		  // Check for overlap 
		  if (rectx[ky][ki].x2_ < rectxy[kz][kj].x1_ ||
		      rectx[ky][ki].x1_ > rectxy[kz][kj].x2_)
		    continue;
		  
		  rectxy[kz][kj].x1_  =
		    std::min(rectxy[kz][kj].x1_, rectx[ky][ki].x1_);
		  rectxy[kz][kj].x2_  =
		    std::max(rectxy[kz][kj].x2_, rectx[ky][ki].x2_);
		  rectxy[kz][kj].y2_ = ky;
		  for (int ka=rectxy[kz][kj].x1_; ka<=rectxy[kz][kj].x2_; ++ka)
		    {
		      int flag = cell_flag_[(kz*num_cell2+ky)*num_cell1+ka];
		      if (flag == UP)
			rectxy[kz][kj].no_pos_++;
		      else if (flag == DOWN)
			rectxy[kz][kj].no_neg_++;
		      else if (flag == SCATTERED)
			rectxy[kz][kj].no_scattered_++;
		    }
		  break;
		}
	      if (kj == rectxy[kz].size())
		{
		  // Define a new change rectangle
		  rectxy[kz].push_back(ChangeRectangle(rectx[ky][ki].x1_,
						       rectx[ky][ki].x2_,
						       ky, ky, kz, kz,
						       rectx[ky][ki].no_pos_,
						       rectx[ky][ki].no_neg_,
						       rectx[ky][ki].no_scattered_));
		}
      // int no_pos=0, no_neg=0, no_scattered = 0;
      // int no_x = (int)knots_u_.size() - 1;
      // int no_y = (int)knots_v_.size() - 1;
      // for (int kb=rectxy[kz][kj].y1_; kb<=rectxy[kz][kj].y2_; ++kb)
      // 	for (int ka=rectxy[kz][kj].x1_; ka<=rectxy[kz][kj].x2_; ++ka)
      // 	  {
      // 	    int flag = cell_flag_[(kz*no_y+kb)*no_x+ka];
      // 	    if (flag == UP)
      // 	      no_pos++;
      // 	    else if (flag == DOWN)
      // 	      no_neg++;
      // 	    else if (flag == SCATTERED)
      // 	      no_scattered++;
      // 	  }
      // if (no_pos != rectxy[kz][kj].no_pos_ || no_neg != rectxy[kz][kj].no_neg_ ||
      // 	  no_scattered != rectxy[kz][kj].no_scattered_)
      // 	std::cout << "Sign mismatch0: " << kz << " " << kj << std::endl;
      	    }
      	  int stop_break4 = 1;
	}
      
      // Include rectangles surrounded by larger rectangles into the large one
      for (size_t ki=0; ki<rectxy[kz].size(); ++ki)
	{
	  int x11 = rectxy[kz][ki].x1_;
	  int x12 = rectxy[kz][ki].x2_;
	  int y11 = rectxy[kz][ki].y1_;
	  int y12 = rectxy[kz][ki].y2_;
	  int z11 = rectxy[kz][ki].z1_;
	  int z12 = rectxy[kz][ki].z2_;
	  for (size_t kj=ki+1; kj<rectxy[kz].size(); )
	    {
	      int x21 = rectxy[kz][kj].x1_;
	      int x22 = rectxy[kz][kj].x2_;
	      int y21 = rectxy[kz][kj].y1_;
	      int y22 = rectxy[kz][kj].y2_;
	      int z21 = rectxy[kz][kj].z1_;
	      int z22 = rectxy[kz][kj].z2_;
	      if (x11 >= x21 && x12 <= x22 && y11 >= y21 && y12 <= y22 &&
		  z11 >= z21 && z12 <= z22)
		{
		  std::swap(rectxy[kz][ki], rectxy[kz][kj]);
		  rectxy[kz][ki].no_pos_ += rectxy[kz][kj].no_pos_;
		  rectxy[kz][ki].no_neg_ += rectxy[kz][kj].no_neg_;
		  rectxy[kz][ki].no_scattered_ += rectxy[kz][kj].no_scattered_;
		  rectxy[kz].erase(rectxy[kz].begin()+kj);
		}
	      else if  (x11 <= x21 && x12 >= x22 && y11 <= y21 && y12 >= y22 &&
			z11 <= z21 && z12 >= z22)
		{
		  rectxy[kz][ki].no_pos_ += rectxy[kz][kj].no_pos_;
		  rectxy[kz][ki].no_neg_ += rectxy[kz][kj].no_neg_;
		  rectxy[kz][ki].no_scattered_ += rectxy[kz][kj].no_scattered_;
		  rectxy[kz].erase(rectxy[kz].begin()+kj);
		}
	      else
		++kj;
	    }
	}
  std::ofstream of2("xy_box.g2");
  for (size_t ki=0; ki<rectxy[kz].size(); ++ki)
    {
      Point low(knots_u_[rectxy[kz][ki].x1_],
		knots_v_[rectxy[kz][ki].y1_],
		knots_w_[rectxy[kz][ki].z1_]);
      Point high(knots_u_[rectxy[kz][ki].x2_+1],
		 knots_v_[rectxy[kz][ki].y2_+1],
		 knots_w_[rectxy[kz][ki].z2_+1]);

      of2 << "410 1 0 4 0 0 0 255" << std::endl;
      of2 << "12" << std::endl;
      Point pt2(high[0], low[1], low[2]);
      Point pt3(low[0], high[1], low[2]);
      Point pt4(high[0], high[1], low[2]);
      Point pt5(low[0], low[1], high[2]);
      Point pt6(high[0], low[1], high[2]);
      Point pt7(low[0], high[1], high[2]);
      of2 << low << " " << pt2 << std::endl;
      of2 << low << " " << pt3 << std::endl;
      of2 << pt2 << " " << pt4 << std::endl;
      of2 << pt3 << " " << pt4 << std::endl;
      of2 << high << " " << pt6 << std::endl;
      of2 << high << " " << pt7 << std::endl;
      of2 << pt6 << " " << pt5 << std::endl;
      of2 << pt7 << " " << pt5 << std::endl;
      of2 << low << " " << pt5 << std::endl;
      of2 << pt2 << " " << pt6 << std::endl;
      of2 << pt3 << " " << pt7 << std::endl;
      of2 << pt4 << " " << high << std::endl;
      // int no_pos=0, no_neg=0, no_scattered = 0;
      // int no_x = (int)knots_u_.size() - 1;
      // int no_y = (int)knots_v_.size() - 1;
      // for (int kb=rectxy[kz][ki].y1_; kb<=rectxy[kz][ki].y2_; ++kb)
      // 	for (int ka=rectxy[kz][ki].x1_; ka<=rectxy[kz][ki].x2_; ++ka)
      // 	  {
      // 	    int flag = cell_flag_[(kz*no_y+kb)*no_x+ka];
      // 	    if (flag == UP)
      // 	      no_pos++;
      // 	    else if (flag == DOWN)
      // 	      no_neg++;
      // 	    else if (flag == SCATTERED)
      // 	      no_scattered++;
      // 	  }
      // if (no_pos != rectxy[kz][ki].no_pos_ || no_neg != rectxy[kz][ki].no_neg_ ||
      // 	  no_scattered != rectxy[kz][ki].no_scattered_)
      // 	std::cout << "Sign mismatch: " << kz << " " << ki << std::endl;
	}
      int stop_break2 = 1;
    }
      
  // Boxes in xyz
  for (int kz=0; kz<num_cell3; ++kz)
    {
      if (rectxy[kz].size() == 0)
	continue;
      
      for (size_t ki=0; ki<rectxy[kz].size(); ++ki)
	{
	  // Check if the interval is part of an existing (possibly
	  // modified change rectangle
	  size_t kj;
	  for (kj=0; kj<change_dom_.size(); ++kj)
	    {
	      if (ki == kj)
		continue;
	      if (change_dom_[kj].z2_ < kz-1)
		continue;  // Not contininous

	      // Check flag
	      int flag1 = rectxy[kz][ki].getMainFlag();
		int flag2 = change_dom_[kj].getMainFlag();
	      if (flag1 != flag2 && !(flag1 == SCATTERED || flag2 == SCATTERED))
		continue;
		  
	      // Check for overlap 
	      if (rectxy[kz][ki].x2_ < change_dom_[kj].x1_ ||
		  rectxy[kz][ki].x1_ > change_dom_[kj].x2_ ||
		  rectxy[kz][ki].y2_ < change_dom_[kj].y1_ ||
		  rectxy[kz][ki].y1_ > change_dom_[kj].y2_)
		continue;
		  
	      change_dom_[kj].x1_  =
		std::min(change_dom_[kj].x1_, rectxy[kz][ki].x1_);
	      change_dom_[kj].x2_  =
		std::max(change_dom_[kj].x2_, rectxy[kz][ki].x2_);
	      change_dom_[kj].y1_  =
		std::min(change_dom_[kj].y1_, rectxy[kz][ki].y1_);
	      change_dom_[kj].y2_  =
		std::max(change_dom_[kj].y2_, rectxy[kz][ki].y2_);
	      change_dom_[kj].z2_ = kz;
	      if (change_dom_[kj].z1_ > change_dom_[kj].z2_)
		std::swap(change_dom_[kj].z1_, change_dom_[kj].z2_);
	      for (int kb=change_dom_[kj].y1_; kb<=change_dom_[kj].y2_; ++kb)
		for (int ka=change_dom_[kj].x1_; ka<=change_dom_[kj].x2_; ++ka)
		  {
		    int flag = cell_flag_[(kz*num_cell2+kb)*num_cell1+ka];
		    if (flag == UP)
		      change_dom_[kj].no_pos_++;
		    else if (flag == DOWN)
		      change_dom_[kj].no_neg_++;
		    else if (flag == SCATTERED)
		      change_dom_[kj].no_scattered_++;
		  }
	      break;
	    }
	  if (kj == change_dom_.size())
	    {
	      // Define a new change rectangle
	      change_dom_.push_back(ChangeRectangle(rectxy[kz][ki].x1_,
						    rectxy[kz][ki].x2_,
						    rectxy[kz][ki].y1_,
						    rectxy[kz][ki].y2_,
						    kz, kz,
						    rectxy[kz][ki].no_pos_,
						    rectxy[kz][ki].no_neg_,
						    rectxy[kz][ki].no_scattered_));
	    }
	}
      int stop_break4 = 1;
    }
  
  // Sort
  vector<double> dom_size(change_dom_.size());
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      dom_size[ki] = (change_dom_[ki].x2_ - change_dom_[ki].x1_)*
	(change_dom_[ki].y2_ - change_dom_[ki].y1_)*
	(change_dom_[ki].z2_ - change_dom_[ki].z1_);
    }
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    for (size_t kj=ki+1; kj<change_dom_.size(); ++kj)
      if (dom_size[kj] > dom_size[ki])
	{
	  std::swap(dom_size[ki], dom_size[kj]);
	  std::swap(change_dom_[ki], change_dom_[kj]);
	}
  
  // Merge across small breaks
  int gap_lim = 2;
  double lapfrac = 0.6;
  size_t num_dom = change_dom_.size();
  while (true)
    {
      for (size_t ki=0; ki<change_dom_.size(); )
	{
	  int x11 = change_dom_[ki].x1_;
	  int x12 = change_dom_[ki].x2_;
	  int y11 = change_dom_[ki].y1_;
	  int y12 = change_dom_[ki].y2_;
	  int z11 = change_dom_[ki].z1_;
	  int z12 = change_dom_[ki].z2_;
	  int flag1 = change_dom_[ki].getMainFlag();
	  size_t kj;
	  for (kj=ki+1; kj<change_dom_.size(); ++kj)
	    {
	      int x21 = change_dom_[kj].x1_;
	      int x22 = change_dom_[kj].x2_;
	      int y21 = change_dom_[kj].y1_;
	      int y22 = change_dom_[kj].y2_;
	      int z21 = change_dom_[kj].z1_;
	      int z22 = change_dom_[kj].z2_;
	      int flag2 = change_dom_[kj].getMainFlag();

	      // Check overlap
	      int lapx = std::min(x12,x22) - std::max(x11,x21);
	      int lapy = std::min(y12,y22) - std::max(y11,y21);
	      double lenx = 0.5*((double)(x12-x11) + (double)(x22-x21));
	      double leny = 0.5*((double)(y12-y11) + (double)(y22-y21));
	      if (lapx < 0 || (double)lapx < lapfrac*lenx ||
		  lapy < 0 || (double)lapy < lapfrac*leny)
		continue;

	      // Check gap
	      if (z11 > z22+gap_lim || z21 > z12+gap_lim)
		continue;

	      if (flag1 != flag2 && !(flag1 == SCATTERED || flag2 == SCATTERED))
		continue;

	      // Combine
	      change_dom_[ki].x1_ = std::min(x11, x21);
	      change_dom_[ki].x2_ = std::max(x12, x22);
	      change_dom_[ki].y1_ = std::min(y11, y21);
	      change_dom_[ki].y2_ = std::max(y12, y22);
	      change_dom_[ki].z1_ = std::min(z11, z21);
	      change_dom_[ki].z2_ = std::max(z12, z22);
	      change_dom_[ki].no_pos_ += change_dom_[kj].no_pos_;
	      change_dom_[ki].no_neg_ += change_dom_[kj].no_neg_;
	      change_dom_[ki].no_scattered_ += change_dom_[kj].no_scattered_;
	      change_dom_.erase(change_dom_.begin()+kj);
	      break;
	    }
	  if (kj == change_dom_.size())
	    ++ki;
	}
      if (change_dom_.size() == num_dom)
	break;
      num_dom = change_dom_.size();
    }

  // Include rectangles surrounded by larger rectangles into the large one
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      int x11 = change_dom_[ki].x1_;
      int x12 = change_dom_[ki].x2_;
      int y11 = change_dom_[ki].y1_;
      int y12 = change_dom_[ki].y2_;
      int z11 = change_dom_[ki].z1_;
      int z12 = change_dom_[ki].z2_;
      for (size_t kj=ki+1; kj<change_dom_.size(); )
	{
	  int x21 = change_dom_[kj].x1_;
	  int x22 = change_dom_[kj].x2_;
	  int y21 = change_dom_[kj].y1_;
	  int y22 = change_dom_[kj].y2_;
	  int z21 = change_dom_[kj].z1_;
	  int z22 = change_dom_[kj].z2_;
	  if (x11 >= x21 && x12 <= x22 && y11 >= y21 && y12 <= y22 &&
	      z11 >= z21 && z12 <= z22)
	    {
	      std::swap(change_dom_[ki], change_dom_[kj]);
	      change_dom_[ki].no_pos_ += change_dom_[kj].no_pos_;
	      change_dom_[ki].no_neg_ += change_dom_[kj].no_neg_;
	      change_dom_[ki].no_scattered_ += change_dom_[kj].no_scattered_;
	      change_dom_.erase(change_dom_.begin()+kj);
	    }
	  else if  (x11 <= x21 && x12 >= x22 && y11 <= y21 && y12 >= y22 &&
		    z11 <= z21 && z12 >= z22)
	    {
	      change_dom_[ki].no_pos_ += change_dom_[kj].no_pos_;
	      change_dom_[ki].no_neg_ += change_dom_[kj].no_neg_;
	      change_dom_[ki].no_scattered_ += change_dom_[kj].no_scattered_;
	      change_dom_.erase(change_dom_.begin()+kj);
	    }
	  else
	    ++kj;
	}
    }

  // Sort
  dom_size.resize(change_dom_.size());
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      dom_size[ki] = (change_dom_[ki].x2_ - change_dom_[ki].x1_)*
	(change_dom_[ki].y2_ - change_dom_[ki].y1_)*
	(change_dom_[ki].z2_ - change_dom_[ki].z1_);
    }
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    for (size_t kj=ki+1; kj<change_dom_.size(); ++kj)
      if (dom_size[kj] > dom_size[ki])
	{
	  std::swap(dom_size[ki], dom_size[kj]);
	  std::swap(change_dom_[ki], change_dom_[kj]);
	}
  
  std::ofstream of2("change_box0.g2");
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      Point low(knots_u_[change_dom_[ki].x1_],
		knots_v_[change_dom_[ki].y1_],
		knots_w_[change_dom_[ki].z1_]);
      Point high(knots_u_[change_dom_[ki].x2_+1],
		 knots_v_[change_dom_[ki].y2_+1],
		 knots_w_[change_dom_[ki].z2_+1]);

      of2 << "410 1 0 4 0 0 0 255" << std::endl;
      of2 << "12" << std::endl;
      Point pt2(high[0], low[1], low[2]);
      Point pt3(low[0], high[1], low[2]);
      Point pt4(high[0], high[1], low[2]);
      Point pt5(low[0], low[1], high[2]);
      Point pt6(high[0], low[1], high[2]);
      Point pt7(low[0], high[1], high[2]);
      of2 << low << " " << pt2 << std::endl;
      of2 << low << " " << pt3 << std::endl;
      of2 << pt2 << " " << pt4 << std::endl;
      of2 << pt3 << " " << pt4 << std::endl;
      of2 << high << " " << pt6 << std::endl;
      of2 << high << " " << pt7 << std::endl;
      of2 << pt6 << " " << pt5 << std::endl;
      of2 << pt7 << " " << pt5 << std::endl;
      of2 << low << " " << pt5 << std::endl;
      of2 << pt2 << " " << pt6 << std::endl;
      of2 << pt3 << " " << pt7 << std::endl;
      of2 << pt4 << " " << high << std::endl;
    }
  
  // Check for overlap and adjacency
  num_dom = change_dom_.size();
  while (true)
    {
      for (size_t ki=0; ki<change_dom_.size(); ++ki)
	{
	  int x11 = change_dom_[ki].x1_;
	  int x12 = change_dom_[ki].x2_;
	  int y11 = change_dom_[ki].y1_;
	  int y12 = change_dom_[ki].y2_;
	  int z11 = change_dom_[ki].z1_;
	  int z12 = change_dom_[ki].z2_;
	  int flag1 = change_dom_[ki].getMainFlag();
	  for (size_t kj=ki+1; kj<change_dom_.size(); ++kj)
	    {
	      int x21 = change_dom_[kj].x1_;
	      int x22 = change_dom_[kj].x2_;
	      int y21 = change_dom_[kj].y1_;
	      int y22 = change_dom_[kj].y2_;
	      int z21 = change_dom_[kj].z1_;
	      int z22 = change_dom_[kj].z2_;
	      int flag2 = change_dom_[kj].getMainFlag();
	  
	      // Check overlap
	      int lapx = std::min(x12,x22) - std::max(x11,x21);
	      int lapy = std::min(y12,y22) - std::max(y11,y21);
	      int lapz = std::min(z12,z22) - std::max(z11,z21);
	      if (lapx >= 0 && lapy >= 0 && lapz >=0)
		{
		  // Overlap
		  resolveOverlap(ki, kj);
		  int stop_lap = 1;
		}
	      else if (lapx >= -1 && lapy >= -1 && lapz >=-1)
		{
		  // Adjacency
		  resolveAdjacency(ki, kj);
		  int stop_adj = 1;
		}
	    }
	}
      if (change_dom_.size() == num_dom)
	break;
      num_dom = change_dom_.size();
    }


  int stop_break5 = 1;
}

// =============================================================================
bool IdentifyChanges::hasChangeFlag(int dom[6], int val, int dir)
// =============================================================================
{
  int dom2[6];
  for (int ki=0; ki<6; ++ki)
    dom2[ki] = dom[ki];
  dom2[2*dir] = dom2[2*dir+1] = val;
  int no_x = (int)knots_u_.size() - 1;
  int no_y = (int)knots_v_.size() - 1;
  int no_z = (int)knots_w_.size() - 1;
  for (int kc=dom2[4]; kc<=dom2[5]; ++kc)
    for (int kb=dom2[2]; kb<=dom2[3]; ++kb)
      for (int ka=dom2[0]; ka<=dom2[1]; ++ka)
	{
	  int flag = cell_flag_[(kc*no_y+kb)*no_x+ka];
	  if (flag != FLAT)
	    return true;
	}
  
  return false;
}

// =============================================================================
void IdentifyChanges::updateFlagNumbers(vector<int>& cell_ix)
// =============================================================================
{
  int no_x = (int)knots_u_.size() - 1;
  int no_y = (int)knots_v_.size() - 1;
  int no_z = (int)knots_w_.size() - 1;
  for (size_t kr=0; kr<cell_ix.size(); ++kr)
    {
      int no1 = 0, no2 = 0, no3 = 0;
      for (int kc=change_dom_[kr].z1_; kc<=change_dom_[kr].z2_; ++kc)
	for (int kb=change_dom_[kr].y1_; kb<=change_dom_[kr].y2_; ++kb)
	  for (int ka=change_dom_[kr].x1_; ka<=change_dom_[kr].x2_; ++ka)
	    {
	      int flag = cell_flag_[(kc*no_y+kb)*no_x+ka];
	      if (flag == UP)
		no1++;
	      else if (flag == DOWN)
		no2++;
	      else if (flag == SCATTERED)
		no3++;
	    }
      change_dom_[kr].setFlagCount(no1, no2, no3);
    }
 }

// =============================================================================
void IdentifyChanges::countFlags(int dom[6], int& no_pos, int& no_neg,
				 int& no_scattered)
// =============================================================================
{
  int no_x = (int)knots_u_.size() - 1;
  int no_y = (int)knots_v_.size() - 1;
  int no_z = (int)knots_w_.size() - 1;
  no_pos = no_neg = no_scattered = 0;
  for (int kc=dom[4]; kc<=dom[5]; ++kc)
    for (int kb=dom[2]; kb<=dom[3]; ++kb)
      for (int ka=dom[0]; ka<=dom[1]; ++ka)
	{
	  int flag = cell_flag_[(kc*no_y+kb)*no_x+ka];
	  if (flag == UP)
	    no_pos++;
	  else if (flag == DOWN)
	    no_neg++;
	  else if (flag == SCATTERED)
	    no_scattered++;
	}
  
  }

// =============================================================================
void IdentifyChanges::collectInvolved(int dom[6], size_t ix, vector<size_t>& involved,
				      int dom2[6])
// =============================================================================
{
  for (int ka=0; ka<6; ++ka)
    dom2[ka] = dom[ka];
  
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      if (ki == ix)
	continue;
      int x3 = change_dom_[ki].x1_;
      int x4 = change_dom_[ki].x2_;
      int y3 = change_dom_[ki].y1_;
      int y4 = change_dom_[ki].y2_;
      int z3 = change_dom_[ki].z1_;
      int z4 = change_dom_[ki].z2_;
      if (((x3 >= dom[0] && x3 <= dom[1]) || (x4 >= dom[0] && x4 <= dom[1])) &&
	  ((y3 >= dom[2] && y3 <= dom[3]) || (y4 >= dom[2] && y4 <= dom[3])) &&
	  ((z3 >= dom[4] && z3 <= dom[5]) || (z4 >= dom[4] && z4 <= dom[5])))
	{
	  involved.push_back(ki);
	  dom2[0] = std::min(dom2[0], x3);
	  dom2[1] = std::max(dom2[1], x4);
	  dom2[2] = std::min(dom2[2], y3);
	  dom2[3] = std::max(dom2[3], y4);
	  dom2[4] = std::min(dom2[4], z3);
	  dom2[5] = std::max(dom2[5], z4);
	}
    }
}

// =============================================================================
void IdentifyChanges::resolveOverlap(size_t ix1, size_t ix2)
// =============================================================================
{
  int dlim[12];
  int ix=0;
  for (size_t kr=ix1; ix<2; ++ix, kr=ix2)
    {
      dlim[6*ix] = change_dom_[kr].x1_;
      dlim[6*ix+1] = change_dom_[kr].x2_;
      dlim[6*ix+2] = change_dom_[kr].y1_;
      dlim[6*ix+3] = change_dom_[kr].y2_;
      dlim[6*ix+4] = change_dom_[kr].z1_;
      dlim[6*ix+5] = change_dom_[kr].z2_;
    }
  
  // Potential extended change box
  int elim[6];
  for (ix=0; ix<6; ix+=2)
    {
      elim[ix] = std::min(dlim[ix], dlim[6+ix]);
      elim[ix+1] = std::max(dlim[ix+1], dlim[6+ix+1]);
    }

  // Overlap area
  int alim[6];
  for (ix=0; ix<6; ix+=2)
    {
      alim[ix] = std::max(dlim[ix], dlim[6+ix]);
      alim[ix+1] = std::min(dlim[ix+1], dlim[6+ix+1]);
    }

  // The flag numbers may be inaccurate. Update
  vector<int> cell_ix(2);
  cell_ix[0] = ix1;
  cell_ix[1] = ix2;
  
  // Count flags in overlap area
  int no_pos = 0, no_neg = 0, no_scattered = 0;
  countFlags(alim, no_pos, no_neg, no_scattered);
  
  // Count flags in extended area
  int no_pos2 = 0, no_neg2 = 0, no_scattered2 = 0;
  countFlags(elim, no_pos2, no_neg2, no_scattered2);
  
  // Collect all other involved boxes
  vector<size_t> cand;
  int elim2[6];
  collectInvolved(elim, ix1, cand, elim2);

  // Compute domain sizes
  int extended=1, extended2=1, overlap=1, dom1=1, dom2=1;
  for (ix=0; ix<6; ix+=2)
    {
      extended *= (elim[ix+1]-elim[ix]+1);
      extended2 *= (elim2[ix+1]-elim2[ix]+1);
      overlap *= (alim[ix+1]-alim[ix]+1);
      dom1 *= (dlim[ix+1]-dlim[ix]+1);
      dom2 *= (dlim[6+ix+1]-dlim[6+ix]+1);
    }
  double domain_frac = (double)(dom1+dom2-overlap)/(double)extended;
  double extended_frac = (double)extended/(double)extended2;

  double frac1 = (double)std::min(change_dom_[ix1].no_pos_, change_dom_[ix1].no_neg_)/
    (double)std::max(change_dom_[ix1].no_pos_, change_dom_[ix1].no_neg_);
  double frac2 = (double)std::min(change_dom_[ix2].no_pos_, change_dom_[ix2].no_neg_)/
    (double)std::max(change_dom_[ix2].no_pos_, change_dom_[ix2].no_neg_);
  double frac3 = (double)std::min(no_pos2, no_neg2)/(double)std::max(no_pos2, no_neg2);
  double fac = 1.5;
  bool split = false;
  if (frac3 > fac*std::max(frac1, frac2))
    split = true;

  if (split)
    {
      int dir = (alim[1]-alim[0] < std::min(alim[3]-alim[2], alim[5]-alim[4])) ? 0 :
	((alim[3]-alim[2] < alim[5]-alim[4]) ? 1 : 2);
      int sz1=1, sz2=1;
      for (ix=0; ix<3; ++ix)
	{
	  if (dir == ix)
	    continue;
	  sz1 *= (dlim[2*ix+1] - dlim[2*ix] + 1);
	  sz2 *= (dlim[6+2*ix+1] - dlim[6+2*ix] + 1);
	}

      int wf;
      bool first = ((no_pos > no_neg && change_dom_[ix1].getMainFlag() == UP) ||
		    (no_pos < no_neg && change_dom_[ix1].getMainFlag() == DOWN));
      int ixk;
      if (sz1 < sz2 || (sz1 == sz2 && first))
	{
	  wf = (dlim[2*dir] < dlim[6+2*dir]) ? dlim[6+2*dir] : dlim[6+2*dir+1];
	  ixk = 6;
	}
      else
	{
	  wf = (dlim[2*dir] < dlim[6+2*dir]) ? dlim[2*dir+1] : dlim[2*dir];
	  ixk = 0;
	}

      // Kept domain
      int ndom[6], ndom3[6];
      int ixk2 = 6 - ixk;
      int kx2 = (dlim[ixk+2*dir] < dlim[ixk2+2*dir]) ? 1 : 0;
      for (ix=0; ix<6; ++ix)
	ndom3[ix] = ndom[ix] = elim[ix];
      ndom3[2*dir+kx2] = ndom[2*dir+kx2] = wf;

      // Check if the kept domain can be reduced
      for (ix=0; ix<3; ++ix)
	{
	  if (ix == dir)
	    continue;
	  int ka;
	  for (ka=ndom[2*ix]; ka<dlim[ixk+2*ix]; ++ka)
	    {
	      bool change = hasChangeFlag(ndom, ka, ix);
	      if (change)
		break;
	    }
	  ndom[2*ix] = ka;
	  for (ka=ndom[2*ix+1]; ka>dlim[ixk+2*ix+1]; --ka)
	    {
	      bool change = hasChangeFlag(ndom, ka, ix);
	      if (change)
		break;
	    }
	  ndom[2*ix+1] = ka;
	}
      int no_pos3 = 0, no_neg3 = 0, no_scattered3 = 0;
      countFlags(ndom, no_pos3, no_neg3, no_scattered3);
      change_dom_[(ixk==0) ? ix1 : ix2] =
	ChangeRectangle(ndom[0], ndom[1], ndom[2], ndom[3], ndom[4], ndom[5],
			no_pos3, no_neg3, no_scattered3);

      // The other
      int ndom2[6];
      for (ix=0; ix<6; ++ix)
	ndom2[ix] = dlim[ixk2+ix];
      ndom2[2*dir+1-kx2] = (kx2 == 0) ? wf-1 : wf+1;
      no_pos3 = no_neg3 = no_scattered3 = 0;
      countFlags(ndom2, no_pos3, no_neg3, no_scattered3);
      change_dom_[(ixk==0) ? ix2 : ix1] =
	ChangeRectangle(ndom2[0], ndom2[1], ndom2[2], ndom2[3], ndom2[4], ndom2[5],
			no_pos3, no_neg3, no_scattered3);
      
	
      int stop_split = 1;
    }
  else
    {
      change_dom_[ix1] = ChangeRectangle(elim[0], elim[1], elim[2], elim[3],
					 elim[4], elim[5], no_pos2,
					 no_neg2, no_scattered2);

      for (int kr=(int)cand.size()-1; kr>=0; --kr)
	{
	  int stat = change_dom_[cand[kr]].insideStatus(elim);
	  if (stat == 1)
	    change_dom_.erase(change_dom_.begin()+cand[kr]);
	}
    }
  int stop_break = 1;
}

// =============================================================================
void IdentifyChanges::resolveAdjacency(size_t ix1, size_t ix2)
// =============================================================================
{
  int dlim[12];
  int ix=0;
  for (size_t kr=ix1; ix<2; ++ix, kr=ix2)
    {
      dlim[6*ix] = change_dom_[kr].x1_;
      dlim[6*ix+1] = change_dom_[kr].x2_;
      dlim[6*ix+2] = change_dom_[kr].y1_;
      dlim[6*ix+3] = change_dom_[kr].y2_;
      dlim[6*ix+4] = change_dom_[kr].z1_;
      dlim[6*ix+5] = change_dom_[kr].z2_;
    }
  
  // Potential extended change box
  int elim[6];
  for (ix=0; ix<6; ix+=2)
    {
      elim[ix] = std::min(dlim[ix], dlim[6+ix]);
      elim[ix+1] = std::max(dlim[ix+1], dlim[6+ix+1]);
    }

  // Update flag numbers
  vector<int> cell_ix(2);
  cell_ix[0] = ix1;
  cell_ix[1] = ix2;
  updateFlagNumbers(cell_ix);

  // Count flags in extended area
  int no_pos = 0, no_neg = 0, no_scattered = 0;
  countFlags(elim, no_pos, no_neg, no_scattered);
  
  // Collect all other involved boxes
  vector<size_t> cand;
  int elim2[6];
  collectInvolved(elim, ix1, cand, elim2);

  // Compute domain sizes
  int extended=1, extended2=1, dom1=1, dom2=1;
  for (ix=0; ix<6; ix+=2)
    {
      extended *= (elim[ix+1]-elim[ix]+1);
      extended2 *= (elim2[ix+1]-elim2[ix]+1);
      dom1 *= (dlim[ix+1]-dlim[ix]+1);
      dom2 *= (dlim[6+ix+1]-dlim[6+ix]+1);
    }
  double domain_frac = (double)(dom1+dom2)/(double)extended;
  double extended_frac = (double)extended/(double)extended2;

  double frac1 = (double)std::min(change_dom_[ix1].no_pos_, change_dom_[ix1].no_neg_)/
    (double)std::max(change_dom_[ix1].no_pos_, change_dom_[ix1].no_neg_);
  double frac2 = (double)std::min(change_dom_[ix2].no_pos_, change_dom_[ix2].no_neg_)/
    (double)std::max(change_dom_[ix2].no_pos_, change_dom_[ix2].no_neg_);
  double frac3 = (double)std::min(no_pos, no_neg)/(double)std::max(no_pos, no_neg);
  double fac = 1.5;
  double ext_fac = 0.8;
  
  if (frac3 <= fac*std::max(frac1, frac2) && extended_frac > ext_fac &&
      domain_frac > ext_fac)
    {
      // Extend box to include adjacent domain
      change_dom_[ix1] = ChangeRectangle(elim[0], elim[1], elim[2], elim[3],
					 elim[4], elim[5], no_pos,
					 no_neg, no_scattered);

      for (int kr=(int)cand.size()-1; kr>=0; --kr)
	{
	  int stat = change_dom_[cand[kr]].insideStatus(elim);
	  if (stat == 1)
	    change_dom_.erase(change_dom_.begin()+cand[kr]);
	}
    }

  int stop_break = 1;
}

// =============================================================================
void IdentifyChanges::computeTimeDeriv()
// =============================================================================
{
  eleminfo_.clear();
  eleminfo_.resize(lrvol_->numElements());
  LRSpline3DEvalGrid eval_grid(lrvol_.get());

  // Loop through elements, sample and compute Bezier coefficients from samples
  int i = 0;
  int dir = 2;
    
  // Storage for points and coefs
  std::vector<double> samplepts(dim_*order_u_*order_v_*order_w_);
  for(auto it=eval_grid.elements_begin(); it!=eval_grid.elements_end();
      it++, i++)
    {
      double ll_x, ll_y, ll_z, ur_x, ur_y, ur_z;
      eval_grid.origLow(*it, ll_x, ll_y, ll_z);
      eval_grid.origHigh(*it, ur_x, ur_y, ur_z);
      Point low(ll_x, ll_y, ll_z);
      Point high(ur_x, ur_y, ur_z);
      eleminfo_[i].tri_box_ = BoundingBox(low, high);

      eval_grid.evaluateGrid(*it, 1, dir, &samplepts[0]);
	
      // Check if the element contains data points
      if ((*it)->hasDataPoints())
	// Compute Bezier coefficients
	LRSpline3DUtils::computeCoefsFromPts(samplepts, order_u_,
					     order_v_, order_w_,
					     dim_, eleminfo_[i].bez2_coefs_);
      int stop_break = 1;
    }
}

   

  

