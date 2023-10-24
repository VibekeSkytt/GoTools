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

    struct ElementInfo {
      BoundingBox tri_box_;
      std::vector<double> bez2_coefs_;
      double min_val_;
      double max_val_;
      double av_val_;

      ElementInfo()
      {
	min_val_ = max_val_ = av_val_ = 0.0;
      }
    };

    struct ChangeRectangle {
      int x1_, x2_, y1_, y2_, z1_, z2_;
      int flag_;

      ChangeRectangle(int x1, int x2, int y1, int y2, int z1, int z2, int flag)
      {
	x1_ = x1;
	x2_ = x2;
	y1_ = y1;
	y2_ = y2;
	z1_ = z1;
	z2_ = z2;
	flag_ = flag;
      }
    };
    
    IdentifyChanges(LRSplineVolume* vol, std::vector<double>& points, bool distribute);

    void getChangeDom(double limit, std::vector<BoundingBox>& tri_dom);
    
  private:
    shared_ptr<LRSplineVolume> lrvol_;
    int order_u_, order_v_, order_w_;
    int dim_;
    int nmb_pts_;
    std::vector<double>& points_;  // Reference to input points and parameter values
    int cdir_;

    std::vector<int> cell_flag_;
    std::vector<ElementInfo> eleminfo_;
    std::vector<ChangeRectangle> change_dom_;
    
    void computeTimeDeriv();

    void setCellFlags(double limit);

    void analyzeCellFlag();
  };
  
}

#endif
