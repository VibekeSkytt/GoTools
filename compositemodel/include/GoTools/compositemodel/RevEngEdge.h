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

#ifndef _REVENGEDGE_H
#define _REVENGEDGE_H

#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/geometry/CurveOnSurface.h"

namespace Go
{
  enum
  {
   BLEND_NOT_SET, NOT_BLEND, STRAIGHT_BLEND, CIRCULAR_BLEND, OTHER_BLEND
  };

  class RevEngEdge
  {
  public:
    // Constructor
    RevEngEdge();

    RevEngEdge(RevEngRegion* reg1, RevEngRegion* reg2);

    RevEngEdge(int type, RevEngRegion* reg1, double dist1,
	       std::vector<shared_ptr<CurveOnSurface> > cvs1,
	       RevEngRegion* reg2, double dist2,
	       std::vector<shared_ptr<CurveOnSurface> > cvs2);

    // Destructor
    ~RevEngEdge();

    void store(std::ostream& os);

    void read(std::istream& is, int& reg_id1, int& reg_id2,
	      int& reg_id3, std::vector<int>& blend_id);

    int getId()
    {
      return Id_;
    }

    void setId(int Id)
    {
      Id_ = Id;
    }

    void setReg1(RevEngRegion *reg);
    
    void setReg2(RevEngRegion *reg);
    
    void addBlendRegion(RevEngRegion* reg)
    {
      blend_regs_.push_back(reg);
    }

    void addBlendRegions(std::vector<RevEngRegion*>& regs)
    {
      blend_regs_.insert(blend_regs_.end(), regs.begin(), regs.end());
    }

    void clearBlendRegions()
    {
      blend_regs_.clear();
    }
    
    std::vector<shared_ptr<ParamCurve> > getSpaceCurves();

    int numBlendRegs()
    {
      return (int)blend_regs_.size();
    }

    RevEngRegion* getBlendReg(int ix)
    {
      if (ix < 0 || ix >= (int)blend_regs_.size())
	return 0;
      else
	return blend_regs_[ix];
    }

    void getAdjacent(RevEngRegion*& reg1, RevEngRegion*& reg2)
    {
      reg1 = adjacent1_;
      reg2 = adjacent2_;
    }

    void getAllBlendRegs(std::vector<RevEngRegion*>& blend_regs)
    {
      if (blend_regs_.size() > 0)
	blend_regs.insert(blend_regs.end(), blend_regs_.begin(),
			  blend_regs_.end());
    }

    void getCurve(std::vector<shared_ptr<CurveOnSurface> >& cvs, bool first=true)
    {
      if (first)
	cvs.insert(cvs.end(), cvs1_.begin(), cvs1_.end());
      else
	cvs.insert(cvs.end(), cvs2_.begin(), cvs2_.end());
    }

    void getDistances(double& dist1, double& dist2)
    {
      dist1 = distance1_;
      dist2 = distance2_;
    }

    void setBlendRegSurf(RevEngRegion* blend)
    {
      defined_blend_ = blend;
    }

    RevEngRegion* getBlendRegSurf()
    {
      return defined_blend_;
    }

    void setAltRadius(double radius)
    {
      alt_rad_ =  radius;
    }

    double getAltRadius()
    {
      return alt_rad_;
    }
	
  private:
    int Id_;
    RevEngRegion* adjacent1_;
    RevEngRegion* adjacent2_;
    RevEngRegion* defined_blend_;
    std::vector<shared_ptr<CurveOnSurface> > cvs1_;
    std::vector<shared_ptr<CurveOnSurface> > cvs2_;
    std::vector<RevEngRegion*> blend_regs_;
    int blend_type_;
    double distance1_;
    double distance2_;
    double alt_rad_;
  };
}

#endif
