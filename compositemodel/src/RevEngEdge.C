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
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ClosestPoint.h"
#include <fstream>

#define DEBUG

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
	  if (!cvs1_[ki]->hasParameterCurve())
	    std::cout << "RevEngEdge::replaceSurf: No parameter curve" << std::endl;
	}
    }
  else if (reg == adjacent2_)
    {
      for (size_t ki=0; ki<cvs2_.size(); ++ki)
	{
	  cvs2_[ki]->setUnderlyingSurface(new_surf);
	  cvs2_[ki]->unsetParameterCurve();
	  cvs2_[ki]->ensureParCrvExistence(tol);
	  if (!cvs2_[ki]->hasParameterCurve())
	    std::cout << "RevEngEdge::replaceSurf: No parameter curve" << std::endl;
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
  int ix;
  closestPoint(pos, par, close, dist, ix);
}

//===========================================================================
void RevEngEdge::closestPoint(const Point& pos, double& par, Point& close,
			      double& dist, int& ix)
//===========================================================================
{
  dist = std::numeric_limits<double>::max();
  ix = -1;
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
	  ix = (int)ki;
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

//===========================================================================
bool RevEngEdge::isAdjacent(RevEngEdge* other, double tol, double& par1, double& par2)
//===========================================================================
{
  double tpar1 = cvs1_[0]->startparam();
  double tpar2 = cvs1_[cvs1_.size()-1]->endparam();
  Point pos1 = cvs1_[0]->ParamCurve::point(tpar1);
  Point pos2 = cvs1_[cvs1_.size()-1]->ParamCurve::point(tpar2);
  double tpar3 = other->startparam();
  double tpar4 = other->endparam();
  Point pos3 = other->point(tpar3);
  Point pos4 = other->point(tpar4);

  double dd1 = pos1.dist(pos3);
  double dd2 = pos1.dist(pos4);
  double dd3 = pos2.dist(pos3);
  double dd4 = pos2.dist(pos4);
  bool adjacent = false;
  if (dd1 <= tol && dd1 <= std::min(dd2, std::min(dd3,dd4)))
    {
      par1 = tpar1;
      par2 = tpar3;
      adjacent = true;
    }
  else if (dd2 <= tol && dd2 <= std::min(dd3, dd4))
    {
      par1 = tpar1;
      par2 = tpar4;
      adjacent = true;
    }
    else if (dd3 <= tol && dd3 <= dd4)
    {
      par1 = tpar2;
      par2 = tpar3;
      adjacent = true;
    }
    else if (dd4 <= tol)
     {
      par1 = tpar2;
      par2 = tpar4;
      adjacent = true;
    }
    return adjacent;
}

//===========================================================================
Point RevEngEdge::point(double par)
//===========================================================================
{
  size_t ix = 0;
  for (ix; ix<cvs1_.size() && cvs1_[ix]->endparam() <= par; ++ix);
  if (ix == cvs1_.size())
    --ix;

  return cvs1_[ix]->ParamCurve::point(par);
}

//===========================================================================
bool RevEngEdge::append(RevEngEdge* other, double tol)
//===========================================================================
{
  double tpar1, tpar2;
  bool adjacent = isAdjacent(other, tol, tpar1, tpar2);
  if (!adjacent)
    return false;   // Cannot append
  
  RevEngRegion *adj3, *adj4;
  other->getAdjacent(adj3, adj4);
  bool first;
  if (adjacent1_ == adj3 && adjacent2_ == adj4)
    first = true;
  else if (adjacent1_ == adj4 && adjacent2_ == adj3)
    first = false;
  else
    return false;  // Not an appendable configuration

  size_t ncv1 = cvs1_.size();
  size_t ix1 = 0;
  for (; ix1<ncv1 && cvs1_[ix1]->endparam() <= tpar1; ++ix1);
  if (ix1 == ncv1)
    --ix1;
  Point ppos1 = cvs1_[ix1]->faceParameter(tpar1);
  Point ppos2 = cvs2_[ix1]->faceParameter(tpar1);
  
  vector<shared_ptr<CurveOnSurface> > cvs3, cvs4;
  other->getCurve(cvs3, true);
  other->getCurve(cvs4, false);
  size_t ncv3 = cvs3.size();
  size_t ix2 = 0;
  for (; ix2<ncv3 && cvs3[ix2]->endparam() <= tpar2; ++ix2);
  if (ix2 == ncv3)
    --ix2;
  Point ppos3 = cvs3[ix2]->faceParameter(tpar2);
  Point ppos4 = cvs4[ix2]->faceParameter(tpar2);

  if ((!adjacent1_->hasSurface()) || (!adjacent2_->hasSurface()))
    return false;  // Unexpected
  
  shared_ptr<ParamSurface> surf1 = adjacent1_->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adjacent2_->getSurface(0)->surface();
  Point par_eps1 = SurfaceTools::getParEpsilon(*surf1, tol);
  Point par_eps2 = SurfaceTools::getParEpsilon(*surf2, tol);
  double epspar1 = 0.5*(par_eps1[0], par_eps1[1]);
  double epspar2 = 0.5*(par_eps2[0], par_eps2[1]);
  double dd1 = ppos1.dist(first ? ppos3 : ppos4);
  double dd2 = ppos2.dist(first ? ppos4 : ppos3);

  if (dd1 > epspar1 || dd2 > epspar2)
    return false;

  // Try to append
  // Copy curves to keep originals in case of failure
  vector<shared_ptr<CurveOnSurface> > cvs1_2(ncv1);
  for (size_t ki=0; ki<ncv1; ++ki)
    cvs1_2[ki] = shared_ptr<CurveOnSurface>(cvs1_[ki]->clone());
  vector<shared_ptr<CurveOnSurface> > cvs2_2(ncv1);
  for (size_t ki=0; ki<ncv1; ++ki)
    cvs2_2[ki] = shared_ptr<CurveOnSurface>(cvs2_[ki]->clone());
  vector<shared_ptr<CurveOnSurface> > cvs3_2(ncv3);
  for (size_t ki=0; ki<ncv3; ++ki)
    cvs3_2[ki] = shared_ptr<CurveOnSurface>(cvs3[ki]->clone());
  vector<shared_ptr<CurveOnSurface> > cvs4_2(ncv3);
  for (size_t ki=0; ki<ncv3; ++ki)
    cvs4_2[ki] = shared_ptr<CurveOnSurface>(cvs4[ki]->clone());

  if (fabs(tpar1-cvs1_2[0]->startparam()) < fabs(cvs1_2[ncv1-1]->endparam()-tpar1))
    {
      for (size_t ki=0; ki<ncv1; ++ki)
	{
	  cvs1_2[ki]->reverseParameterDirection();
	  cvs2_2[ki]->reverseParameterDirection();
	}
      for (size_t ki=0; ki<ncv1/2; ++ki)
	{
	  std::swap(cvs1_2[ki], cvs1_2[ncv1-1-ki]);
	  std::swap(cvs2_2[ki], cvs2_2[ncv1-1-ki]);
	}
    }

  if (fabs(tpar2-cvs3_2[0]->startparam()) > fabs(cvs3_2[ncv3-1]->endparam()-tpar2))
    {
      for (size_t ki=0; ki<ncv3; ++ki)
	{
	  cvs3_2[ki]->reverseParameterDirection();
	  cvs4_2[ki]->reverseParameterDirection();
	}
      for (size_t ki=0; ki<ncv3/2; ++ki)
	{
	  std::swap(cvs3_2[ki], cvs3_2[ncv3-1-ki]);
	  std::swap(cvs4_2[ki], cvs4_2[ncv3-1-ki]);
	}
    }

  // Do append
  double dist1, dist2;
  if (first)
    {
      cvs1_2[ncv1-1]->appendCurve(cvs3_2[0].get(), 1, dist1, false, epspar1);
      cvs2_2[ncv1-1]->appendCurve(cvs4_2[0].get(), 1, dist2, false, epspar2);
    }
  else
    {
      cvs1_2[ncv1-1]->appendCurve(cvs4_2[0].get(), 1, dist1, false, epspar1);
      cvs2_2[ncv1-1]->appendCurve(cvs2_2[0].get(), 1, dist2, false, epspar2);
    }

  if (dist1 > tol || dist2 > tol)
    return false;

  // Check that the parameter curves of the joined curves exists
  if ((!cvs1_2[ncv1-1]->hasParameterCurve()) ||
      (!cvs2_2[ncv1-1]->hasParameterCurve()))
    return false;

  // Replace curves
  cvs1_.clear();
  cvs1_.insert(cvs1_.end(), cvs1_2.begin(), cvs1_2.end());
  if (first && cvs3_2.size() > 1)
    cvs1_.insert(cvs1_.end(), cvs3_2.begin()+1, cvs3_2.end());
  else if ((!first) && cvs4_2.size() > 1)
    cvs1_.insert(cvs1_.end(), cvs4_2.begin()+1, cvs4_2.end());

  cvs2_.clear();
  cvs2_.insert(cvs2_.end(), cvs2_2.begin(), cvs2_2.end());
  if (first && cvs4_2.size() > 1)
    cvs2_.insert(cvs2_.end(), cvs4_2.begin()+1, cvs4_2.end());
  else if ((!first) && cvs4_2.size() > 1)
    cvs2_.insert(cvs2_.end(), cvs3_2.begin()+1, cvs3_2.end());

  blend_regs_.insert(blend_regs_.end(), other->blend_regs_.begin(),
		     other->blend_regs_.end());

  distance_ = 0.5*(distance_ + other->distance_);
  radius_ = 0.5*(radius_ + other->radius_);

  return true;
}

//===========================================================================
int RevEngEdge::missingParCrv()
//===========================================================================
{
  int missing = 0;
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    if (!cvs1_[ki]->hasParameterCurve())
      {
	missing += 1;
	break;
      }
  
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    if (!cvs2_[ki]->hasParameterCurve())
      {
	missing += 2;
	break;
      }
  return missing;
}

//===========================================================================
void RevEngEdge::splitAtSeam(double tol,
			     vector<shared_ptr<RevEngEdge> >& added_edgs,
			     vector<shared_ptr<RevEngRegion> >& added_regs,
			     vector<shared_ptr<HedgeSurface> >& added_sfs)
//===========================================================================
{
  if ((!adjacent1_->hasSurface()) || (!adjacent2_->hasSurface()))
    return;  // Something is wrong
#ifdef DEBUG
  std::ofstream of("edge_split.g2");
  adjacent1_->writeRegionPoints(of);
  adjacent2_->writeRegionPoints(of);
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp = cvs1_[ki]->spaceCurve();
      tmp->writeStandardHeader(of);
      tmp->write(of);
    }
#endif
  shared_ptr<ParamSurface> surf1 = adjacent1_->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adjacent2_->getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);

  bool close_u1=false, close_u2=false, close_v1=false, close_v2=false;
  if (elem1.get())
    elem1->isClosed(close_u1, close_v1);
  if (elem2.get())
    elem2->isClosed(close_u2, close_v2);

  RectDomain dom1 = surf1->containingDomain();
  RectDomain dom2 = surf2->containingDomain();
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      if (!cvs1_[ki]->hasParameterCurve())
	{
	  // Try to split curve
	  if (close_u1)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf1->constParamCurves(dom1.umin(), false);

	      shared_ptr<ParamCurve> space = cvs1_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(), 
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 1;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	  
	  if (close_v1)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf1->constParamCurves(dom1.vmin(), true);

	      shared_ptr<ParamCurve> space = cvs1_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(), 
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 1;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	}
    }

  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    {
      if (!cvs2_[ki]->hasParameterCurve())
	{
	  // Try to split curve
	  if (close_u2)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf2->constParamCurves(dom2.umin(), false);

	      shared_ptr<ParamCurve> space = cvs2_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(), 
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 2;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	  
	  if (close_v1)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf1->constParamCurves(dom1.vmin(), true);

	      shared_ptr<ParamCurve> space = cvs1_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(),
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 2;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	}
    }

  
}


//===========================================================================
shared_ptr<RevEngEdge>
RevEngEdge::doSplit(size_t ix, int side, double par, double tol,
		    vector<shared_ptr<RevEngRegion> >& added_regs,
		    vector<shared_ptr<HedgeSurface> >& added_sfs)
//===========================================================================
{
  shared_ptr<RevEngEdge> new_edg;

  if (par <= cvs1_[ix]->startparam() || par >= cvs1_[ix]->endparam())
    return new_edg;

  // Split curves
  vector<shared_ptr<ParamCurve> > sub1 = cvs1_[ix]->split(par);
  vector<shared_ptr<ParamCurve> > sub2 = cvs2_[ix]->split(par);
  if (sub1.size() != 2 || sub2.size() != 2)
    return new_edg;
  
  // Split blend regions
  vector<RevEngRegion*> move_reg;
  for (size_t ki=0; ki<blend_regs_.size(); )
    {
      vector<RevEngPoint*> points = blend_regs_[ki]->getPoints();
      vector<RevEngPoint*> keep, move;
      
      for (size_t kj=0; kj<points.size(); ++kj)
	{
	  Vector3D xyz = points[kj]->getPoint();
	  Point pnt(xyz[0], xyz[1], xyz[2]);
	  double par2, dist;
	  Point close;
	  int ix2;
	  closestPoint(pnt, par2, close, dist, ix2);
	  if (ix2 < ix || par2 <= par)
	    keep.push_back(points[kj]);
	  else
	    move.push_back(points[kj]);
	}

      if (move.size() == 0)
	++ki; // Do nothing
      else if (keep.size() == 0)
	{
	  // Move regions to new edge
	  move_reg.push_back(blend_regs_[ki]);
	  blend_regs_.erase(blend_regs_.begin()+ki);
	}
      else
	{
	  // Split region
	  blend_regs_[ki]->removePoints(move);  // This is not the most effective
	  // method, but the simplest to implement
	  blend_regs_[ki]->updateInfo();

	  shared_ptr<RevEngRegion> new_reg(new RevEngRegion(blend_regs_[ki]->getClassificationType(),
							    blend_regs_[ki]->getEdgeClassificationType(),
							    move));
	  added_regs.push_back(new_reg);
	  move_reg.push_back(new_reg.get());
	  ++ki;
	}
    }

  // Distribute curves
  vector<shared_ptr<CurveOnSurface> > cvs1_2, cvs2_2;
  cvs1_2.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub1[1]));
  cvs1_2[0]->ensureParCrvExistence(tol);
  for (size_t ki=ix+1; ki<cvs1_.size(); ++ki)
    cvs1_2.push_back(cvs1_[ki]);
  cvs2_2.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub2[1]));
  cvs2_2[0]->ensureParCrvExistence(tol);
  for (size_t ki=ix+1; ki<cvs2_.size(); ++ki)
    cvs2_2.push_back(cvs2_[ki]);

  cvs1_[ix] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub1[0]);
  cvs1_[ix]->ensureParCrvExistence(tol);
  cvs2_[ix] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub2[0]);
  cvs2_[ix]->ensureParCrvExistence(tol);
  if (ix < cvs1_.size()-1)
    {
      cvs1_.erase(cvs1_.begin()+ix+1, cvs1_.end());
      cvs2_.erase(cvs1_.begin()+ix+1, cvs1_.end());
    }

  if (defined_blend_)
    {
      // Split associated surface
      int stop_blend = 1;
    }

  new_edg = shared_ptr<RevEngEdge>(new RevEngEdge(blend_type_, adjacent1_,
						  cvs1_2, outer1_, adjacent2_,
						  cvs2_2, outer2_, radius_,
						  distance_));
  if (move_reg.size() > 0)
    new_edg->addBlendRegions(move_reg);
  
  return new_edg;
}

