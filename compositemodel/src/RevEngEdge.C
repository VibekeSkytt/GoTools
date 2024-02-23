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

#include "GoTools/compositemodel/RevEngEdge.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::istream;
using std::ostream;

//===========================================================================
RevEngEdge::RevEngEdge()
//===========================================================================
  : adjacent1_(0), adjacent2_(0), defined_blend_(0),
    blend_type_(BLEND_NOT_SET), distance1_(0.0), distance2_(0.0),
    alt_rad_(-1.0)
{
}

//===========================================================================
RevEngEdge::RevEngEdge(RevEngRegion* reg1, RevEngRegion* reg2)
//===========================================================================
  : adjacent1_(reg1), adjacent2_(reg2), defined_blend_(0),
    blend_type_(BLEND_NOT_SET), distance1_(0.0), distance2_(0.0),
    alt_rad_(-1.0)
{
}

//===========================================================================
RevEngEdge::RevEngEdge(int type, RevEngRegion* reg1, double dist1,
		       vector<shared_ptr<CurveOnSurface> > cvs1,
		       RevEngRegion* reg2, double dist2,
		       vector<shared_ptr<CurveOnSurface> > cvs2)
//===========================================================================
  : adjacent1_(reg1), adjacent2_(reg2), defined_blend_(0), blend_type_(type),
    distance1_(dist1), distance2_(dist2),
    alt_rad_(-1.0)
{
  cvs1_.insert(cvs1_.end(), cvs1.begin(), cvs1.end());
  cvs2_.insert(cvs2_.end(), cvs2.begin(), cvs2.end());
}

//===========================================================================
RevEngEdge::~RevEngEdge()
//===========================================================================
{
}


//===========================================================================
void RevEngEdge::setReg1(RevEngRegion *reg)
//===========================================================================
{
  adjacent1_ = reg;
  shared_ptr<ParamSurface> surf;
  if (reg->hasSurface())
    surf = reg->getSurface(0)->surface();
  else if (reg->hasBaseSf())
    surf = reg->getBase();

  if (surf.get())
    {
      for (size_t ki=0; ki<cvs1_.size(); ++ki)
	cvs1_[ki]->setUnderlyingSurface(surf);
    }
}

//===========================================================================
void RevEngEdge::setReg2(RevEngRegion *reg)
//===========================================================================
{
  adjacent2_ = reg;
  shared_ptr<ParamSurface> surf;
  if (reg->hasSurface())
    surf = reg->getSurface(0)->surface();
  else if (reg->hasBaseSf())
    surf = reg->getBase();

  if (surf.get())
    {
      for (size_t ki=0; ki<cvs2_.size(); ++ki)
	cvs2_[ki]->setUnderlyingSurface(surf);
    }
}

//===========================================================================
void RevEngEdge::store(ostream& os)
//===========================================================================
{
  os << Id_ << std::endl;
  int id1 = adjacent1_->getId();
  int id2 = adjacent2_->getId();
  int id3 = (defined_blend_ == 0) ? -1 : defined_blend_->getId();
  os << id1 << " " << id2 << " " << id3 << std::endl;
  os << cvs1_.size() << " " << cvs2_.size() << std::endl;
  
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      bool space = cvs1_[ki]->hasSpaceCurve();
      bool param = cvs1_[ki]->hasParameterCurve();
      os << space << " " << param << std::endl;
      if (space)
	{
	  shared_ptr<ParamCurve> spacecv = cvs1_[ki]->spaceCurve();
	  shared_ptr<ParamCurve> parcv = cvs1_[ki]->parameterCurve();
	  spacecv->writeStandardHeader(os);
	  spacecv->write(os);
	  parcv->writeStandardHeader(os);
	  parcv->write(os);
	}
    }
      
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    {
      bool space = cvs2_[ki]->hasSpaceCurve();
      bool param = cvs2_[ki]->hasParameterCurve();
      os << space << " " << param << std::endl;
      if (space)
	{
	  shared_ptr<ParamCurve> spacecv = cvs2_[ki]->spaceCurve();
	  shared_ptr<ParamCurve> parcv = cvs2_[ki]->parameterCurve();
	  spacecv->writeStandardHeader(os);
	  spacecv->write(os);
	  parcv->writeStandardHeader(os);
	  parcv->write(os);
	}
    }

  os << blend_regs_.size() << std::endl;
  for (size_t ki=0; ki < blend_regs_.size(); ++ki)
    os << blend_regs_[ki]->getId() << " ";
  os << std::endl;

  os << blend_type_ << " " << distance1_ << " " << distance2_ << std::endl;
      
}


//===========================================================================
void RevEngEdge::read(istream& is, int& reg_id1, int& reg_id2, int& reg_id3,
		      vector<int>& blend_id)
//===========================================================================
{
  is >> Id_;
  is >> reg_id1 >> reg_id2 >> reg_id3;
  int num_cv1, num_cv2;
  is >> num_cv1 >> num_cv2;
  if (num_cv1 > 0)
    cvs1_.resize(num_cv1);
  if (num_cv2 > 0)
    cvs2_.resize(num_cv2);
  for (int ka=0; ka<num_cv1; ++ka)
    {
      int space, param;
      is >> space >> param;
      shared_ptr<ParamCurve> spacecv, parcv;
      shared_ptr<ParamSurface> dummy;
      if (space)
	{
	  ObjectHeader header;
	  header.read(is);
	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
	  obj->read(is);
	  spacecv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
	}

      if (param)
	{
	  ObjectHeader header;
	  header.read(is);
	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
	  obj->read(is);
	  parcv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
	}
      
      cvs1_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface(dummy, parcv,
								spacecv, false,
								-1, -1, 0.0,
								-1, true));
    }

  for (int ka=0; ka<num_cv2; ++ka)
    {
      int space, param;
      is >> space >> param;
      shared_ptr<ParamCurve> spacecv, parcv;
      shared_ptr<ParamSurface> dummy;
      if (space)
	{
	  ObjectHeader header;
	  header.read(is);
	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
	  obj->read(is);
	  spacecv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
	}

      if (param)
	{
	  ObjectHeader header;
	  header.read(is);
	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
	  obj->read(is);
	  parcv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
	}
      
      cvs2_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface(dummy, parcv,
								spacecv, false,
								-1, -1, 0.0,
								-1, true));
    }

  int num_blend_reg;
  is >> num_blend_reg;
  if (num_blend_reg > 0)
    blend_id.resize(num_blend_reg);
  for (int ka=0; ka<num_blend_reg; ++ka)
    is >> blend_id[ka];
  
  is >> blend_type_ >> distance1_ >> distance2_;
}


//===========================================================================
vector<shared_ptr<ParamCurve> > RevEngEdge::getSpaceCurves()
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > curves;
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    curves.push_back(cvs1_[ki]->spaceCurve());

  return curves;
}
