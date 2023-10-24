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

using namespace Go;
using std::vector;

// =============================================================================
IdentifyChanges::IdentifyChanges(LRSplineVolume* vol, vector<double>& points,
				 bool distribute)
// =============================================================================
  : points_(points), cdir_(0)
{
  lrvol_ = shared_ptr<LRSplineVolume>(vol->clone());
  order_u_ = lrvol_->degree(XDIR) + 1;
  order_v_ = lrvol_->degree(YDIR) + 1;
  order_w_ = lrvol_->degree(ZDIR) + 1;
  dim_ = lrvol_->dimension();

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
  setCellFlags(limit);

  analyzeCellFlag();

  const Mesh3D& mesh = lrvol_->mesh();
  for (size_t ki=0; ki<change_dom_.size(); ++ki)
    {
      Point low(mesh.kval(XDIR, change_dom_[ki].x1_),
		mesh.kval(YDIR, change_dom_[ki].y1_),
		mesh.kval(ZDIR, change_dom_[ki].z1_));
      Point high(mesh.kval(XDIR, change_dom_[ki].x2_+1),
		 mesh.kval(YDIR, change_dom_[ki].y2_+1),
		 mesh.kval(ZDIR, change_dom_[ki].z2_+1));
      tri_dom.push_back(BoundingBox(low, high));
    }
  
}

// =============================================================================
void IdentifyChanges::setCellFlags(double limit)
// =============================================================================
{
  double eps = 1.0e-8;
  
  // Compute derivative in time direction defined as Bezier patches
  computeTimeDeriv();
  
  // Traverse the Bezier representation of the time derivatires and
  // identify patches with high variability that are candidate for
  // split
  const Mesh3D& mesh = lrvol_->mesh();
  vector<LRSplineVolume::Refinement3D> mesh_refs;
  for (size_t kj=0; kj<eleminfo_.size(); ++kj)
    {
      if (eleminfo_[kj].bez2_coefs_.size() == 0)
	continue;   // Element without points
      
      double min_val =  std::numeric_limits<double>::max();
      double max_val = std::numeric_limits<double>::lowest();
      double av_val = 0.0;
      for (size_t kr=0; kr<eleminfo_[kj].bez2_coefs_.size(); kr+=dim_)
	{
	  double curr_coef = eleminfo_[kj].bez2_coefs_[kr+cdir_];
	  min_val = std::min(min_val, curr_coef);
	  max_val = std::max(max_val, curr_coef);
	  av_val += curr_coef;
	}
      av_val /= eleminfo_[kj].bez2_coefs_.size()/dim_;

      eleminfo_[kj].min_val_ = min_val;
      eleminfo_[kj].max_val_ = max_val;
      eleminfo_[kj].av_val_ = av_val;

      if (std::max(fabs(min_val), fabs(max_val)) < limit)
	continue;  // No change identified

      if (min_val*max_val >= 0.0)
	continue;   // Consistent direction

      if (std::min(fabs(min_val), fabs(max_val)) < 0.1*limit)
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
  
  int num_cell1 = mesh.numDistinctKnots(XDIR) - 1;
  int num_cell2 = mesh.numDistinctKnots(YDIR) - 1;
  int num_cell3 = mesh.numDistinctKnots(ZDIR) - 1;
  cell_flag_.resize(num_cell1*num_cell2*num_cell3);
  std::fill(cell_flag_.begin(), cell_flag_.end(), 0);
  
  for (size_t kj=0; kj<eleminfo_.size(); ++kj)
    {
      if (eleminfo_[kj].bez2_coefs_.size() == 0)
	continue;   // Element without points
      if (std::max(fabs(eleminfo_[kj].min_val_),
		   fabs(eleminfo_[kj].max_val_)) < limit)
	continue;  // No change identified

      // Define flag
      int flag = 10;  // Indicates both directions
      if (fabs(eleminfo_[kj].min_val_) < limit ||
	  eleminfo_[kj].min_val_ > 0.0 /*|| eleminfo_[kj].av_val_ > 0.1*limit*/)
	flag = 1;
      else if (fabs(eleminfo_[kj].max_val_) < limit ||
	       eleminfo_[kj].max_val_ < 0.0 /*|| eleminfo_[kj].av_val_ < -0.1*limit*/)
	flag = -1;
      
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
	    cell_flag_[kc*xydel+kb*num_cell1+ka] = flag;

      int stop_break = 1;
    }

  analyzeCellFlag();
  int stop_break2 = 1;
 }

// =============================================================================
void IdentifyChanges::analyzeCellFlag()
// =============================================================================
{
  const Mesh3D& mesh = lrvol_->mesh();
  int num_cell1 = mesh.numDistinctKnots(XDIR) - 1;
  int num_cell2 = mesh.numDistinctKnots(YDIR) - 1;
  int num_cell3 = mesh.numDistinctKnots(ZDIR) - 1;

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
	      for (; kx<num_cell1 && cell_flag_[kstart+kx] == 0; ++kx);
	      int flag = cell_flag_[kstart+kx];
	      for (kx2=kx+1; kx2<num_cell1 &&
		     (cell_flag_[kstart+kx2] == cell_flag_[kstart+kx] ||
		      cell_flag_[kstart+kx2] == 10); ++kx2)
		if (flag == 10)
		  flag =  cell_flag_[kstart+kx2];
	      if (kx < num_cell1)
		rectx[ky].push_back(ChangeRectangle(kx, kx2-1, ky, ky, kz, kz, flag));
	    }

	  // Merge across small breaks
	  bool domerge = true; //false;
	  if (domerge)
	    {
	  int stop_break1 = 1;
	  for (size_t ki=1; ki<rectx[ky].size(); )
	    {
	      int flag1 = rectx[ky][ki-1].flag_;
	      int flag2 = rectx[ky][ki].flag_;
	      if (abs(flag1) == 1 && abs(flag2) == 1 && flag1*flag2 < 0)
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
		  rectx[ky][ki-1].flag_ = (flag1 == 10) ? flag2 : flag1;
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
		    continue;  // Not contininous

		  // Check flag
		  int flag1 = rectx[ky][ki].flag_;
		  int flag2 = rectxy[kz][kj].flag_;
		  if (abs(flag1) == 1 && abs(flag2) == 1 && flag1*flag2 < 0)
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
		  rectxy[kz][kj].flag_ = (flag1 == 10) ? flag2 : flag1;
		  break;
		}
	      if (kj == rectxy[kz].size())
		{
		  // Define a new change rectangle
		  rectxy[kz].push_back(ChangeRectangle(rectx[ky][ki].x1_,
						       rectx[ky][ki].x2_,
						       ky, ky, kz, kz,
						       rectx[ky][ki].flag_));
		}
	    }
	  int stop_break4 = 1;
	}
      // Include rectangles surrounded by larger rectangles into the large one
      for (size_t ki=1; ki<rectxy[kz].size(); )
	{
	  int x11 = rectxy[kz][ki-1].x1_;
	  int x12 = rectxy[kz][ki-1].x2_;
	  int y11 = rectxy[kz][ki-1].y1_;
	  int y12 = rectxy[kz][ki-1].y2_;
	  int x21 = rectxy[kz][ki].x1_;
	  int x22 = rectxy[kz][ki].x2_;
	  int y21 = rectxy[kz][ki].y1_;
	  int y22 = rectxy[kz][ki].y2_;
	  if (x11 >= x21 && x12 <= x22 && y11 >= y21 && y12 <= y22)
	    {
	      rectxy[kz][ki].flag_ =
		(rectxy[kz][ki-1].flag_ == rectxy[kz][ki].flag_) ?
		rectxy[kz][ki].flag_ : 20;
	      rectxy[kz].erase(rectxy[kz].begin()+ki-1);
	    }
	  else if  (x11 <= x21 && x12 >= x22 && y11 <= y21 && y12 >= y22)
	    {
	      rectxy[kz][ki-1].flag_ = 20;
	      rectxy[kz].erase(rectxy[kz].begin()+ki);
	    }
	  else
	    ++ki;
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
	      if (change_dom_[kj].z2_ < kz-1)
		continue;  // Not contininous

	      // Check flag
	      int flag1 = rectxy[kz][ki].flag_;
	      int flag2 = change_dom_[kj].flag_;
	      if (abs(flag1) == 1 && abs(flag2) == 1 && flag1*flag2 < 0)
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
	      change_dom_[kj].flag_ = (flag1 == 10) ? flag2 : flag1;
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
						    rectxy[kz][ki].flag_));
	    }
	}
      int stop_break4 = 1;
    }
  
  // Include rectangles surrounded by larger rectangles into the large one
  for (size_t ki=1; ki<change_dom_.size(); )
    {
      int x11 = change_dom_[ki-1].x1_;
      int x12 = change_dom_[ki-1].x2_;
      int y11 = change_dom_[ki-1].y1_;
      int y12 = change_dom_[ki-1].y2_;
      int z11 = change_dom_[ki-1].z1_;
      int z12 = change_dom_[ki-1].z2_;
      int x21 = change_dom_[ki].x1_;
      int x22 = change_dom_[ki].x2_;
      int y21 = change_dom_[ki].y1_;
      int y22 = change_dom_[ki].y2_;
      int z21 = change_dom_[ki].z1_;
      int z22 = change_dom_[ki].z2_;
      if (x11 >= x21 && x12 <= x22 && y11 >= y21 && y12 <= y22 &&
	  z11 >= z21 && z12 <= z22)
	{
	  change_dom_[ki].flag_ = 20;
	  change_dom_.erase(change_dom_.begin()+ki-1);
	}
      else if  (x11 <= x21 && x12 >= x22 && y11 <= y21 && y12 >= y22 &&
		z11 <= z21 && z12 >= z22)
	{
	  change_dom_[ki-1].flag_ = 20;
	  change_dom_.erase(change_dom_.begin()+ki);
	}
      else
	++ki;
    }
  
  int stop_break5 = 1;
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

   

  

