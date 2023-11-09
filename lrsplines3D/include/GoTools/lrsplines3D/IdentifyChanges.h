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

#ifndef _IDENTIFYCHANGES_H
#define _IDENTIFYCHANGES_H

#include "GoTools/utils/BoundingBox.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"

namespace Go
{
// =============================================================================
  class IdentifyChanges
// =============================================================================
  {
  public:

    enum
      {
       FLAT, DOWN, UP, DOWN_FLAT, DOWN_UP, FLAT_DOWN, FLAT_UP, UP_DOWN, UP_FLAT, SCATTERED
      };

    struct ElementInfo {
      BoundingBox tri_box_;
      std::vector<double> bez2_coefs_;
      double min_val_;
      double max_val_;
      double av_val_;
      double acc_neg_;
      double acc_pos_;
      int same_neg_, same_pos_;
      bool time_change_;

      ElementInfo()
      {
	min_val_ = max_val_ = av_val_ = 0.0;
	acc_neg_ = acc_pos_ = 0.0;
	same_neg_ = same_pos_  = 0;
	time_change_ = false;
      }
    };

    struct ChangeRectangle {
      int x1_, x2_, y1_, y2_, z1_, z2_;
      int no_pos_, no_neg_, no_scattered_;

      ChangeRectangle(int x1, int x2, int y1, int y2, int z1, int z2,
		      int no_pos, int no_neg, int no_scattered)
      {
	x1_ = x1;
	x2_ = x2;
	y1_ = y1;
	y2_ = y2;
	z1_ = z1;
	z2_ = z2;
	no_pos_ = no_pos;
	no_neg_ = no_neg;
	no_scattered_ = no_scattered;
      }

      void setFlagCount(int no_pos, int no_neg, int no_scattered)
      {
	no_pos_ = no_pos;
	no_neg_ = no_neg;
	no_scattered_ = no_scattered;
      }
      
      int getMainFlag()
      {
	int flag = (no_pos_ > std::max(no_neg_, no_scattered_) ? UP :
		    ((no_neg_ > no_scattered_) ? DOWN : SCATTERED));
	return flag;
      }

      int insideStatus(int dom[6])
      {
	int status;   // 0 = out, 1 = in, 2 = crosses boundary
	if (x1_>=dom[0] && x2_<=dom[1] && y1_>=dom[2] && y2_<=dom[3] &&
	    z1_>=dom[4] && z2_<=dom[5])
	  status = 1;
	else if (dom[0]>x2_ || dom[1]<x1_ || dom[2]>y2_ || dom[3]<y1_ ||
		 dom[4]>z2_ || dom[5]<z1_)
	  status = 0;
	else
	  status = 2;
	return status;
      }
    };
    
    IdentifyChanges(LRSplineVolume* vol, std::vector<double>& points, bool distribute);

    void getChangeDom(double limit, std::vector<BoundingBox>& tri_dom);

    void getDerSize(double& max_der, double& min_der, double& av_der)
    {
      max_der = max_der_;
      min_der = min_der_;
      av_der = av_der_;
    }
    
  private:
    shared_ptr<LRSplineVolume> lrvol_;
    int order_u_, order_v_, order_w_;
    std::vector<double> knots_u_, knots_v_, knots_w_;
    int dim_;
    int nmb_pts_;
    std::vector<double>& points_;  // Reference to input points and parameter values
    int cdir_;

    std::vector<int> cell_flag_;
    std::vector<ElementInfo> eleminfo_;
    std::vector<ChangeRectangle> change_dom_;
    double min_der_, max_der_, av_der_;
    
    void computeTimeDeriv();

    void computeElemInfo();

    void setCellFlags(double limit);

    void analyzeCellFlag();

    void splitCells();

    void resolveOverlap(size_t ix1, size_t ix2);

    void resolveAdjacency(size_t ix1, size_t ix2);

    bool hasChangeFlag(int dom[6], int val, int dir);

    void updateFlagNumbers(std::vector<int>& cell_ix);

    void countFlags(int dom[6], int& no_pos, int& no_neg,  int& no_scattered);

    void collectInvolved(int dom[6], size_t ix, std::vector<size_t>& involved,
			 int dom2[6]);
  };
  
}

#endif
