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
    blend_type_(BLEND_NOT_SET), outer1_(false), outer2_(false),
    distance_(0.0), radius_(0.0)
{
}

//===========================================================================
RevEngEdge::RevEngEdge(RevEngRegion* reg1, RevEngRegion* reg2)
//===========================================================================
  : adjacent1_(reg1), adjacent2_(reg2), defined_blend_(0),
    blend_type_(BLEND_NOT_SET), outer1_(false), outer2_(false),
    distance_(0.0), radius_(0.0)
{
}

//===========================================================================
RevEngEdge::RevEngEdge(int type, RevEngRegion* reg1, 
		       vector<shared_ptr<CurveOnSurface> > cvs1,
		       bool out1, RevEngRegion* reg2,
		       vector<shared_ptr<CurveOnSurface> > cvs2,
		       bool out2, double distance, double radius)
//===========================================================================
  : adjacent1_(reg1), adjacent2_(reg2), defined_blend_(0), blend_type_(type),
    outer1_(out1), outer2_(out2), distance_(distance), radius_(radius)
{
  cvs1_.insert(cvs1_.end(), cvs1.begin(), cvs1.end());
  cvs2_.insert(cvs2_.end(), cvs2.begin(), cvs2.end());
}

//===========================================================================
RevEngEdge::~RevEngEdge()
//===========================================================================
{
  if (adjacent1_)
    adjacent1_->removeRevEngEdge(this);
  if (adjacent2_)
    adjacent2_->removeRevEngEdge(this);
  for (size_t ki=0; ki<blend_regs_.size(); ++ki)
    blend_regs_[ki]->removeAssociatedBlend();
  if (defined_blend_)
    defined_blend_->removeBlendEdge();
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
void RevEngEdge::fixMismatchCurves(double tol)
//===========================================================================
{
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    cvs1_[ki]->fixMismatchCurves(tol);
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    cvs2_[ki]->fixMismatchCurves(tol);
}

//===========================================================================
void RevEngEdge::replaceSurf(RevEngRegion* reg,
			     shared_ptr<ParamSurface>& new_surf, double tol)
//===========================================================================
{
  if (reg == adjacent1_)
    {
      for (size_t ki=0; ki<cvs1_.size(); ++ki)
	{
	  cvs1_[ki]->setUnderlyingSurface(new_surf);
	  cvs1_[ki]->unsetParameterCurve();
	  cvs1_[ki]->ensureParCrvExistence(tol);
	}
    }
  else if (reg == adjacent2_)
    {
      for (size_t ki=0; ki<cvs2_.size(); ++ki)
	{
	  cvs2_[ki]->setUnderlyingSurface(new_surf);
	  cvs2_[ki]->unsetParameterCurve();
	  cvs2_[ki]->ensureParCrvExistence(tol);
	}
     }
      
}

//===========================================================================
void RevEngEdge::store(ostream& os)
//===========================================================================
{
  os << Id_ << std::endl;
  int id1 = (adjacent1_) ? adjacent1_->getId() : -1;
  int id2 = (adjacent2_) ? adjacent2_->getId() : -1;
  int id3 = (defined_blend_ == 0) ? -1 : defined_blend_->getId();
  os << id1 << " " << id2 << " " << id3 << std::endl;
  os << cvs1_.size() << " " << cvs2_.size() << std::endl;
  
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      cvs1_[ki]->write(os);
      // bool space = cvs1_[ki]->hasSpaceCurve();
      // bool param = cvs1_[ki]->hasParameterCurve();
      // os << space << " " << param << std::endl;
      // if (space)
      // 	{
      // 	  shared_ptr<ParamCurve> spacecv = cvs1_[ki]->spaceCurve();
      // 	  shared_ptr<ParamCurve> parcv = cvs1_[ki]->parameterCurve();
      // 	  spacecv->writeStandardHeader(os);
      // 	  spacecv->write(os);
      // 	  parcv->writeStandardHeader(os);
      // 	  parcv->write(os);
      // 	}
    }
      
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    {
      cvs2_[ki]->write(os);
      // bool space = cvs2_[ki]->hasSpaceCurve();
      // bool param = cvs2_[ki]->hasParameterCurve();
      // os << space << " " << param << std::endl;
      // if (space)
      // 	{
      // 	  shared_ptr<ParamCurve> spacecv = cvs2_[ki]->spaceCurve();
      // 	  shared_ptr<ParamCurve> parcv = cvs2_[ki]->parameterCurve();
      // 	  spacecv->writeStandardHeader(os);
      // 	  spacecv->write(os);
      // 	  parcv->writeStandardHeader(os);
      // 	  parcv->write(os);
      // 	}
    }

  os << blend_regs_.size() << std::endl;
  for (size_t ki=0; ki < blend_regs_.size(); ++ki)
    os << blend_regs_[ki]->getId() << " ";
  os << std::endl;

  os << blend_type_ << " " << distance_ << " " << radius_ << " ";
  os << outer1_ << " " << outer2_ << std::endl;
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
      cvs1_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface());
      cvs1_[ka]->read(is);
      // int space, param;
      // is >> space >> param;
      // shared_ptr<ParamCurve> spacecv, parcv;
      // shared_ptr<ParamSurface> dummy;
      // if (space)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  spacecv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}

      // if (param)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  parcv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}
      
      // cvs1_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface(dummy, parcv,
      // 								spacecv, false,
      // 								-1, -1, 0.0,
      // 								-1, true));
    }

  for (int ka=0; ka<num_cv2; ++ka)
    {
      cvs2_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface());
      cvs2_[ka]->read(is);
      // int space, param;
      // is >> space >> param;
      // shared_ptr<ParamCurve> spacecv, parcv;
      // shared_ptr<ParamSurface> dummy;
      // if (space)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  spacecv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}

      // if (param)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  parcv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}
      
      // cvs2_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface(dummy, parcv,
      // 								spacecv, false,
      // 								-1, -1, 0.0,
      // 								-1, true));
    }

  int num_blend_reg;
  is >> num_blend_reg;
  if (num_blend_reg > 0)
    blend_id.resize(num_blend_reg);
  for (int ka=0; ka<num_blend_reg; ++ka)
    is >> blend_id[ka];
  
  is >> blend_type_ >> distance_ >> radius_ >> outer1_ >> outer2_;
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

//===========================================================================
void RevEngEdge::getCrvEndPoints(Point& pos1, Point& pos2)
//===========================================================================
{
  if (cvs1_.size() > 0)
    {
      cvs1_[0]->point(pos1, cvs1_[0]->startparam());
      cvs1_[cvs1_.size()-1]->point(pos2, cvs1_[cvs1_.size()-1]->endparam());
    }
}

//===========================================================================
void RevEngEdge::closestPoint(const Point& pos, double& par, Point& close,
			      double& dist)
//===========================================================================
{
  dist = std::numeric_limits<double>::max();
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      double par1, dist1;
      Point close1;
      cvs1_[ki]->closestPoint(pos, cvs1_[ki]->startparam(),
			      cvs1_[ki]->endparam(), par1, close1, dist1);
      if (dist1 < dist)
	{
	  par = par1;
	  close = close1;
	  dist = dist1;
	}
    }
}

//===========================================================================
int RevEngEdge::closedSfAtEnd(double tol, double& par, Point& pos, bool at_start)
//===========================================================================
{
  shared_ptr<ParamSurface> surf1 = cvs1_[0]->underlyingSurface();
  shared_ptr<ParamSurface> surf2 = cvs2_[0]->underlyingSurface();

  // The parameter interval and the space curve of the two curve sequences
  // correspond
  int ix = (at_start) ? 0 : (int)cvs1_.size()-1;
  par = (at_start) ? cvs1_[ix]->startparam() : cvs1_[ix]->endparam();
  pos = cvs1_[ix]->ParamCurve::point(par);

  double eps = 1.0e-9;
  double u1, u2, v1, v2, d1, d2;
  Point close1, close2;
  surf1->closestBoundaryPoint(pos, u1, v1, close1, d1, eps);
  surf2->closestBoundaryPoint(pos, u2, v2, close2, d2, eps);
  if (d1 > tol && d2 > tol)
    return 0;

  if (d1 <= tol)
    {
      RectDomain dom = surf1->containingDomain();
      if (fabs(u1-dom.umin()) <= eps)
	{
	  Point pos2 = surf1->point(dom.umax(), v1);
	  if (pos.dist(pos2))
	    return 1;
	}
      else if (fabs(dom.umax()-u1) <= eps)
	{
	  Point pos2 = surf1->point(dom.umin(), v1);
	  if (pos.dist(pos2))
	    return 1;
	}
      if (fabs(v1-dom.vmin()) <= eps)
	{
	  Point pos2 = surf1->point(u1, dom.vmax());
	  if (pos.dist(pos2))
	    return 1;
	}
      else if (fabs(dom.vmax()-v1) <= eps)
	{
	  Point pos2 = surf1->point(u1, dom.vmin());
	  if (pos.dist(pos2))
	    return 1;
	}
    }

  if (d2 <= tol)
    {
      RectDomain dom = surf2->containingDomain();
      if (fabs(u2-dom.umin()) <= eps)
	{
	  Point pos2 = surf2->point(dom.umax(), v2);
	  if (pos.dist(pos2))
	    return 2;
	}
      else if (fabs(dom.umax()-u2) <= eps)
	{
	  Point pos2 = surf2->point(dom.umin(), v2);
	  if (pos.dist(pos2))
	    return 2;
	}
      if (fabs(v2-dom.vmin()) <= eps)
	{
	  Point pos2 = surf2->point(u2, dom.vmax());
	  if (pos.dist(pos2))
	    return 2;
	}
      else if (fabs(dom.vmax()-v2) <= eps)
	{
	  Point pos2 = surf2->point(u2, dom.vmin());
	  if (pos.dist(pos2))
	    return 2;
	}
    }

  return 0;
}

//===========================================================================
bool RevEngEdge::isClosed(double tol)
//===========================================================================
{
  double par1 = cvs1_[0]->startparam();
  double par2 = cvs1_[cvs1_.size()-1]->endparam();
  Point pos1 = cvs1_[0]->ParamCurve::point(par1);
  Point pos2 = cvs1_[cvs1_.size()-1]->ParamCurve::point(par2);
  double dist = pos1.dist(pos2);
  return (dist <= tol);
}

