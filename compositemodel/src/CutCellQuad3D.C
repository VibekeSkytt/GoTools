/* (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
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

#include "GoTools/compositemodel/CutCellQuad3D.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ModifyFaceSet.h"
#include "GoTools/creators/CutCellQuad.h"
#include <fstream>

#define DEBUG

using std::vector;
using std::set;
using std::pair;
using std::make_pair;
using namespace Go;

//==============================================================================
CutCellQuad3D::CutCellQuad3D(shared_ptr<SurfaceModel> sfmodel, double tol)
//==============================================================================
  : tol_(tol)
{
  body_ = shared_ptr<Body>(new Body(sfmodel));
}

//==============================================================================
int CutCellQuad3D::cellStat(const Point& ll, const Point& ur, int& coinc)
//==============================================================================
{
  coinc = 0;
  vector<shared_ptr<SplineSurface> > cell_sfs;
  createCellSfs(ll, ur, cell_sfs);

  vector<BoundingBox> sidebb(cell_sfs.size());
  for (size_t kr=0; kr<cell_sfs.size(); ++kr)
    sidebb[kr] = cell_sfs[kr]->boundingBox();
  
  tpTolerances top = body_->getTolerances();
  vector<shared_ptr<ParamSurface> > cell_sfs2(cell_sfs.begin(), cell_sfs.end());
  shared_ptr<SurfaceModel> cell_mod(new SurfaceModel(top.gap, top.gap, top.neighbour,
						     top.kink, top.bend, cell_sfs2));
  shared_ptr<Body> cell_body(new Body(cell_mod));
  BoundingBox cellbb = cell_mod->boundingBox();

  // Check for intersection and coincidence
  double eps = 1.0e-6; 
  int nmbshell = body_->nmbOfShells();
  int ki, kj;
  Identity ident;
  int intern = -1; 
  for (ki=0; ki<nmbshell; ++ki)
    {
      shared_ptr<SurfaceModel> shell = body_->getShell(ki);
      int nmbface = shell->nmbEntities();
      for (kj=0; kj<nmbface; ++kj)
	{
	  shared_ptr<ParamSurface> sf = shell->getSurface(kj);
	  double u, v;
	  Point sfmid = sf->getInternalPoint(u, v);
	  if (((fabs(ll[0]-sfmid[0]) < tol_ || fabs(ur[0]-sfmid[0]) < tol_) &&
	       sfmid[1] >= ll[1] && sfmid[1] <= ur[1] && sfmid[2] >= ll[2] && sfmid[2] <= ur[2]) ||
	      ((fabs(ll[1]-sfmid[1]) < tol_ || fabs(ur[1]-sfmid[1]) < tol_) &&
		sfmid[0] >= ll[0] && sfmid[0] <= ur[0] && sfmid[2] >= ll[2] && sfmid[2] <= ur[2]) ||
	      ((fabs(ll[2]-sfmid[2]) < tol_ || fabs(ur[2]-sfmid[2]) < tol_) &&
		sfmid[0] >= ll[0] && sfmid[0] <= ur[0] && sfmid[1] >= ll[1] && sfmid[1] <= ur[1]))
	    {
	      // Potentially coincident surfaces
	      ;
	    }
	  else
	    {
	      bool sfin = cell_body->isInside(sfmid);
	      if (intern == -1)
		intern = (sfin) ? 1 : 0;
	      else if ((sfin && (!intern)) || ((!sfin) && intern))
		return 2;
	    }
	  BoundingBox sfbb = sf->boundingBox();
	  if (!sfbb.overlaps(cellbb))
	    continue;
	  
	  for (size_t kr=0; kr<cell_sfs.size(); ++kr)
	    {
	      if (!sfbb.overlaps(sidebb[kr]))
		continue;
#ifdef DEBUG
	      std::ofstream of("curr_sfs.g2");
	      sf->writeStandardHeader(of);
	      sf->write(of);
	      cell_sfs[kr]->writeStandardHeader(of);
	      cell_sfs[kr]->write(of);
#endif
	      // Possible intersection. Check for coincidence
	      int stat = ident.identicalSfs(sf, cell_sfs[kr], tol_);
	      if (stat != 0)
		{
		  coinc = 1;
		  continue;   // Avoid performing surface-surface intersection
		}

	      // Intersect
	      shared_ptr<BoundedSurface> bd1, bd2;
	      vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
	      BoundedUtils::getSurfaceIntersections(sf, cell_sfs[kr], eps,
						    int_cv1, bd1,
						    int_cv2, bd2, true);
	      if (int_cv1.size() > 0 || int_cv2.size() > 0)
		return 2;
	    }
	}
    }
  // Not a cut cell (except possible in connection with coincident boundary surfaces)
  // Check if cell midpoint is inside the shells
  Point mid = 0.5*(ll + ur);
  bool inside = body_->isInside(mid);
  return (inside) ? 1 : 0;
}

//==============================================================================
void
CutCellQuad3D::quadrature(const Point& ll, const Point& ur,
			  vector<double>& quadraturepoints,
			  vector<double>& pointsweights,
			  vector<vector<shared_ptr<ParamSurface> > >& unresolved_cells,
			  vector<double>& surfquads,
			  vector<double>& surfnorms,
			  vector<double>& sfptweights,
			  vector<vector<shared_ptr<ParamSurface> > >& small_sfs,
			  int stat, int coinc)
//==============================================================================
{
  if (stat < 0 || stat > 2)
    stat = cellStat(ll, ur, coinc);

    if (stat == 0)
    return;  // Outside, no quadrature points

  else if (stat == 1)
    {
      // Not a cut cell. Adapt quadrature points to cell size
      vector<double> quadpar1(quadpar_.size());
      vector<double> quadpar2(quadpar_.size());
      vector<double> quadpar3(quadpar_.size());
      vector<double> wgt1(quadpar_.size());
      vector<double> wgt2(quadpar_.size());
      vector<double> wgt3(quadpar_.size());
      double del1 = ur[0] - ll[0];
      double del2 = ur[1] - ll[1];
      double del3 = ur[2] - ll[2];
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar1[ki] = ll[0] + quadpar_[ki]*del1;
	  wgt1[ki] = weights_[ki]*del1;
	}
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar2[ki] = ll[1] + quadpar_[ki]*del2;
	  wgt2[ki] = weights_[ki]*del2;
	}
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar3[ki] = ll[2] + quadpar_[ki]*del3;
	  wgt3[ki] = weights_[ki]*del3;
	}
      quadraturepoints.resize(3*quadpar_.size()*quadpar_.size()*quadpar_.size());
      pointsweights.resize(quadpar_.size()*quadpar_.size()*quadpar_.size());
      for (size_t kk=0; kk<quadpar_.size(); ++kk)
	for (size_t kj=0; kj<quadpar_.size(); ++kj)
	  for (size_t ki=0; ki<quadpar_.size(); ++ki)
	    {
	      quadraturepoints[3*((kk*quadpar_.size()+kj)*quadpar_.size()+ki)] = quadpar1[ki];
	      quadraturepoints[3*((kk*quadpar_.size()+kj)*quadpar_.size()+ki)+1] = quadpar2[kj];
	      quadraturepoints[3*((kk*quadpar_.size()+kj)*quadpar_.size()+ki)+2] = quadpar3[kk];
	      pointsweights[(kk*quadpar_.size()+kj)*quadpar_.size()+ki] = wgt1[ki]*wgt2[kj]*wgt3[kk];
	    }

      // Check for coincidence between a piece of the boundary loop and the
      // cell boundaries
      if (coinc != 0)
	{
	  // Coincidence possible or identified
	}
    }
  else
    {
      // Cut cell
      std::cout << "Cut cell" << std::endl;
      vector<shared_ptr<Body> > cutcells;
      createCutCell(ll, ur, cutcells);

      for (size_t ki=0; ki<cutcells.size(); ++ki)
	{
	  quadraturePoints(ll, ur, cutcells[ki], quadraturepoints, pointsweights, unresolved_cells,
			   surfquads, surfnorms, sfptweights, small_sfs);
	}

      int stop_break0 = 1;
    }

    int stop_break = 1; 
}

//==============================================================================
void CutCellQuad3D::splitCell(shared_ptr<Body> body,
			      shared_ptr<ParamSurface> splitsf,
			      vector<shared_ptr<Body> >& subcell)
//==============================================================================
{
  tpTolerances top = body->getTolerances();
  vector<shared_ptr<ParamSurface> > split_sfs;
  split_sfs.push_back(splitsf);
  shared_ptr<SurfaceModel> other_mod(new SurfaceModel(top.gap, top.gap, top.neighbour,
						     top.kink, top.bend, split_sfs));

  // Intersect CAD model (exclusive voids) and element and sort pieces according to
  // inside and outside
  shared_ptr<SurfaceModel> outer = body->getShell(0);
  vector<shared_ptr<SurfaceModel> > split_mod = outer->splitSurfaceModels(other_mod);
  
  // Check for identical faces between the surfaces of the initial face
  // internal to the initial volume that are coincident with a surface in any
  // of the new volume pieces. These surfaces must be removed
  removeCoincFaces(split_mod[0], split_mod[1], split_mod[2], top.gap);

  int nmb = split_mod[2].get() ? split_mod[2]->nmbEntities() : 0;
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ftSurface> face1 = split_mod[2]->getFace(ki);
      shared_ptr<ParamSurface> surf1 = split_mod[2]->getSurface(ki);
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();
      
      shared_ptr<ftSurface> face1_2(new ftSurface(surf1, -1));
      if (split_mod[0].get())
	split_mod[0]->append(face1_2, true, false, true);
      else
	{
	  vector<shared_ptr<ftSurface> > faces;
	  faces.push_back(face1_2);
	  split_mod[0] = shared_ptr<SurfaceModel>(new SurfaceModel(top.gap, top.gap,
								   top.neighbour,
								   top.kink, top.bend,
								   faces));
	}

      shared_ptr<ftSurface> face2(new ftSurface(surf2, -1));
      if (split_mod[1].get())
	split_mod[1]->append(face2, true, false, true);
      else
	{
	  vector<shared_ptr<ftSurface> > faces;
	  faces.push_back(face2);
	  split_mod[1] = shared_ptr<SurfaceModel>(new SurfaceModel(top.gap, top.gap,
								   top.neighbour,
								   top.kink, top.bend,
								   faces));
	}
    }

  if (body->nmbOfShells() > 1)
    MESSAGE("Voids not handled");

  std::ofstream of1("submod1.g2");
  std::ofstream of2("submod2.g2");
  int nmb1 = split_mod[0].get() ? split_mod[0]->nmbEntities() : 0;
  int nmb2 = split_mod[1].get() ? split_mod[1]->nmbEntities() : 0;
  for (int ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> sf = split_mod[0]->getSurface(ki);
      sf->writeStandardHeader(of1);
      sf->write(of1);
    }
  for (int ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamSurface> sf = split_mod[1]->getSurface(ki);
      sf->writeStandardHeader(of2);
      sf->write(of2);
    }
  
  for (int ki=0; ki<2; ++ki)
    {
      if (!split_mod[ki].get())
	continue;
      
        int nmbbd = split_mod[ki]->nmbBoundaries();
	if (nmbbd > 0)
	  {
	    vector<shared_ptr<ftEdge> > bdedg = split_mod[ki]->getBoundaryEdges();
#ifdef DEBUG
	    std::ofstream of2("ededg.g2");
	    for (size_t kr=0; kr<bdedg.size(); ++kr)
	      {
		shared_ptr<ParamCurve> cv = bdedg[kr]->geomCurve();
		shared_ptr<SplineCurve> splcv(cv->geometryCurve());
		splcv->writeStandardHeader(of2);
		splcv->write(of2);
	      }
#endif
      
	    THROW("Cut cell not closed");
	  }

	  vector<shared_ptr<SurfaceModel> > connected_mod = split_mod[ki]->getConnectedModels();

	  for (size_t kj=0; kj<connected_mod.size(); ++kj)
	    subcell.push_back(shared_ptr<Body>(new Body(connected_mod[kj])));


    }
  int stop_break = 1;
}

//==============================================================================
void CutCellQuad3D::createCutCell(const Point& ll, const Point& ur,
				  vector<shared_ptr<Body> >& cutcell)
//==============================================================================
{
  vector<shared_ptr<SplineSurface> > cell_sfs;
  createCellSfs(ll, ur, cell_sfs);

  // Make surface model
  tpTolerances top = body_->getTolerances();
  vector<shared_ptr<ParamSurface> > cell_sfs2(cell_sfs.begin(), cell_sfs.end());
  shared_ptr<SurfaceModel> cell_mod(new SurfaceModel(top.gap, top.gap, top.neighbour,
						     top.kink, top.bend, cell_sfs2));
  shared_ptr<Body> cell_body(new Body(cell_mod));

  // Intersect CAD model (exclusive voids) and element and sort pieces according to
  // inside and outside
  shared_ptr<SurfaceModel> outer = body_->getShell(0);
  vector<shared_ptr<SurfaceModel> > sub_mod = cell_mod->splitSurfaceModels(outer);

  // Set code for surface origin
  int nmb1 = (sub_mod[0].get()) ? sub_mod[0]->nmbEntities() : 0;
  for (int ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> sf = sub_mod[0]->getSurface(ki);
      if (sf->getFlagCode() != 3)
	sf->setFlagCode(1);
    }
  // Collect faces being internal to both surface models
  Identity identity;
  int nmb2 = (sub_mod[2].get()) ? sub_mod[2]->nmbEntities() : 0;
  if (nmb2 > 0)
    {
      vector<shared_ptr<ftSurface> > addfaces;
      int nmb2 = sub_mod[2]->nmbEntities();
      for (int ki=0; ki<nmb2; ++ki)
	{
	  shared_ptr<ParamSurface> sf = sub_mod[2]->getSurface(ki);

	  // Check for coincidence
	  bool coinc = false;
	  for (int kj=0; kj<nmb1; ++kj)
	    {
	      shared_ptr<ParamSurface> sf1 = sub_mod[0]->getSurface(kj);
	      int stat = identity.identicalSfs(sf1, sf, tol_);
	      if (stat > 0)
		{
		  sf1->setFlagCode(3);
		  coinc = true;
		  break;
		}
	    }
	  if (!coinc)
	    {
	      sf->setFlagCode(2);
	      addfaces.push_back(sub_mod[2]->getFace(ki));
	    }
	}

      if (addfaces.size() > 0)
	{
	  for (size_t kr=0; kr<addfaces.size(); ++kr)
	    sub_mod[2]->removeFace(addfaces[kr]);
	  if (nmb1 > 0)
	    sub_mod[0]->append(addfaces);
	  else
	    sub_mod[0] = shared_ptr<SurfaceModel>(new SurfaceModel(top.gap, top.gap,
								   top.neighbour,
								   top.kink, top.bend,
								   addfaces));
	}
    }

  // Voids
#ifdef DEBUG
  std::ofstream of("cutcell.g2");
  int nmb = (sub_mod[0].get()) ? sub_mod[0]->nmbEntities() : 0;
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> sf = sub_mod[0]->getSurface(ki);
      sf->writeStandardHeader(of);
      sf->write(of);
    }
#endif
  int nmbbd = sub_mod[0]->nmbBoundaries();
  if (nmbbd > 0)
    {
      vector<shared_ptr<ftEdge> > bdedg = sub_mod[0]->getBoundaryEdges();
      std::ofstream of2("ededg.g2");
      for (size_t kr=0; kr<bdedg.size(); ++kr)
	{
	  shared_ptr<ParamCurve> cv = bdedg[kr]->geomCurve();
	  shared_ptr<SplineCurve> splcv(cv->geometryCurve());
	  splcv->writeStandardHeader(of2);
	  splcv->write(of2);
	}
      
      THROW("Cut cell not closed");
    }

  vector<shared_ptr<SurfaceModel> > connected_mod = sub_mod[0]->getConnectedModels();

  for (size_t kj=0; kj<connected_mod.size(); ++kj)
    cutcell.push_back(shared_ptr<Body>(new Body(connected_mod[kj])));

  int stop_break = 1;
}

//===========================================================================
// 
// 
void
CutCellQuad3D::removeCoincFaces(shared_ptr<SurfaceModel>& mod1,
				shared_ptr<SurfaceModel>& mod2,
				shared_ptr<SurfaceModel>& mod3,
				double tol)

//===========================================================================
{
  // Check for surfaces in mod3 that are identical with surfaces in
  // mod1 or mod2. Remove those surfaces
  int nmb = (!mod3.get()) ? 0 : mod3->nmbEntities();
  for (int kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ParamSurface> surf1 = mod3->getSurface(kj);
      int nmb2 = (!mod1.get()) ? 0 : mod1->nmbEntities();
      int kr;
      for (kr=0; kr<nmb2; ++kr)
	{
	  shared_ptr<ParamSurface> surf2 = mod1->getSurface(kr);
	  Identity ident;
	  int coinc = ident.identicalSfs(surf1, surf2, tol);
	  if (coinc > 0)
	    break;
	}
      if (kr < nmb2)
	{
	  mod3->removeFace(mod3->getFace(kj));
	  --kj;
	  --nmb;
	}
      else
	{
	  nmb2 = (!mod2.get()) ? 0 : mod2->nmbEntities();
	  for (kr=0; kr<nmb2; ++kr)
	    {
	      shared_ptr<ParamSurface> surf2 = mod2->getSurface(kr);
	      Identity ident;
	      int coinc = ident.identicalSfs(surf1, surf2, tol);
	      if (coinc > 0)
		break;
	    }
	  if (kr < nmb2)
	    {
	      mod3->removeFace(mod3->getFace(kj));
	      --kj;
	      --nmb;
	    }
	}
    }
}

//==============================================================================
void CutCellQuad3D::createCellSfs(const Point& ll, const Point& ur,
				  vector<shared_ptr<SplineSurface> >& cell_sfs)
//==============================================================================
{
  cell_sfs.resize(6);

  double coefs1[] = {ll[0], ll[1], ll[2], ll[0], ur[1], ll[2], ur[0], ll[1], ll[2], ur[0], ur[1], ll[2]};
  double coefs2[] = {ll[0], ll[1], ur[2], ur[0], ll[1], ur[2], ll[0], ur[1], ur[2], ur[0], ur[1], ur[2]};
  double coefs3[] = {ll[0], ll[1], ll[2], ll[0], ll[1], ur[2], ll[0], ur[1], ll[2], ll[0], ur[1], ur[2]};
  double coefs4[] = {ur[0], ll[1], ll[2], ur[0], ur[1], ll[2], ur[0], ll[1], ur[2], ur[0], ur[1], ur[2]};
  double coefs5[] = {ll[0], ll[1], ll[2], ur[0], ll[1], ll[2], ll[0], ll[1], ur[2], ur[0], ll[1], ur[2]};
  double coefs6[] = {ll[0], ur[1], ll[2], ll[0], ur[1], ur[2], ur[0], ur[1], ll[2], ur[0], ur[1], ur[2]};

  double knots1[] = {ll[0], ll[0], ur[0], ur[0]};
  double knots2[] = {ll[1], ll[1], ur[1], ur[1]};
  double knots3[] = {ll[2], ll[2], ur[2], ur[2]};

  cell_sfs[0] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots2[0],
							    &coefs1[0], 3, false));
  cell_sfs[1] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots2[0],
							    &coefs2[0], 3, false));
  cell_sfs[2] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots2[0], &knots3[0],
							    &coefs3[0], 3, false));
  cell_sfs[3] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots2[0], &knots3[0],
							    &coefs4[0], 3, false));
  cell_sfs[4] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots3[0],
							    &coefs5[0], 3, false));
  cell_sfs[5] = shared_ptr<SplineSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0], &knots3[0],
							    &coefs6[0], 3, false));
#ifdef DEBUG
  std::ofstream of("cell_sfs.g2");
  for (size_t ki=0; ki<cell_sfs.size(); ++ki)
    {
      cell_sfs[ki]->writeStandardHeader(of);
      cell_sfs[ki]->write(of);
    }
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << ll << std::endl;
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << ur << std::endl;

  int stop_break = 1;
#endif
}

//==============================================================================
void
CutCellQuad3D::quadraturePoints(const Point& ll, const Point& ur,
				shared_ptr<Body> body,
				vector<double>& quadraturepoints,
				vector<double>& pointsweights,
				vector<vector<shared_ptr<ParamSurface> > >& unresolved_cells,
				vector<double>& surfquads,
				vector<double>& surfnorms,
				vector<double>& sfptweights,
				vector<vector<shared_ptr<ParamSurface> > >& small_sfs)
//==============================================================================
{
  int splitdir;
  double splitval;
  BoundingBox bbox = body->boundingBox();
  Point ll3 = ll;
  Point ur3 = ur;
  shared_ptr<SplineCurve> rulecv;
  Point ruledir;
  int basedir = selectBaseDir(ll3, ur3, body, splitdir, splitval, rulecv, ruledir);
  if (basedir < 0)
    {
      if (splitdir >= 0)
	{
	  // Split cut cell model and treat sub models recursively
	  Point pos = 0.5*(bbox.low() + bbox.high());
	  pos[splitdir] = splitval;
	  Point norm(0.0, 0.0, 0.0);
	  norm[splitdir] = 1.0;

	  shared_ptr<Plane> splitplane(new Plane(pos, norm));
	  int udir = (splitdir == 0) ? 1 : 0;
	  int vdir = (splitdir == 2) ? 1 : 2;
	  double udel = bbox.high()[udir] - bbox.low()[udir];
	  double vdel = bbox.high()[vdir] - bbox.low()[vdir];
	  splitplane->setParameterBounds(-udel, -vdel, udel, vdel);

	  vector<shared_ptr<Body> > subcell;
	  splitCell(body, splitplane, subcell);

	  for (size_t ki=0; ki<subcell.size(); ++ki)
	    quadraturePoints(ll, ur, subcell[ki], quadraturepoints, pointsweights,
			     unresolved_cells, surfquads, surfnorms, sfptweights, small_sfs);

	  int stop_break = 1;
	}
      else if (rulecv.get())
	{
	  // Linearily extend curve to make sure it covers the entire model
	  double len = bbox.low().dist(bbox.high());
	  rulecv->enlarge(0.5*len, true);
	  rulecv->enlarge(0.5*len, false);

	  // Create ruled surface
	  shared_ptr<SplineCurve> cv2(rulecv->clone());
	  Point vec = len*ruledir;
	  cv2->translateCurve(vec);

	  vector<double> coefs(rulecv->coefs_begin(), rulecv->coefs_end());
	  coefs.insert(coefs.end(), cv2->coefs_begin(), cv2->coefs_end());
	  vector<double> knots2(4);
	  knots2[0] = knots2[1] = 0.0;
	  knots2[2] = knots2[3] = len;
	  shared_ptr<SplineSurface> ruledsf(new SplineSurface(rulecv->numCoefs(), 2,
							      rulecv->order(), 2,
							      rulecv->basis().begin(),
							      &knots2[0], &coefs[0],
							      rulecv->dimension(), rulecv->rational()));
	  std::ofstream ofrule("ruledsf.g2");
	  ruledsf->writeStandardHeader(ofrule);
	  ruledsf->write(ofrule);
	  
	  vector<shared_ptr<Body> > subcell;
	  splitCell(body, ruledsf, subcell);

	  for (size_t ki=0; ki<subcell.size(); ++ki)
	    quadraturePoints(ll, ur, subcell[ki], quadraturepoints, pointsweights,
			     unresolved_cells, surfquads, surfnorms, sfptweights, small_sfs);

	  int stop_rule = 1;
							      
	}
      else
	std::cout << "No suitable base direction found" << std::endl;
    }
  else
    {
#ifdef DEBUG
      std::ofstream ofp("currquads.g2");
#endif
      int numshell = body->nmbOfShells();

      // 2D cell
      Point ll2(2), ur2(2);
      int ix = 0;
      for (int kc=0; kc<3; ++kc)
	{
	  if (kc == basedir)
	    continue;
	  ll2[ix] = ll[kc];
	  ur2[ix++] = ur[kc];
	}
      
      double tmin = ll3[basedir];
      double tmax = ur3[basedir];
      double del = tmax - tmin;

      // Create 2D models
      double quadpar;
      double wgt;
      Point pos = 0.5*(bbox.low() + bbox.high());
      Point norm(0.0, 0.0, 0.0);
      norm[basedir] = 1.0;
      for (size_t ki=0; ki<quadpar_.size(); ++ki)
	{
	  quadpar = tmin + quadpar_[ki]*del;
	  wgt = weights_[ki]*del;

	  // Intersect with plane through quadrature parameter
	  pos[basedir] = quadpar;
	  ftPlane plane(norm, pos);
	  
#ifdef DEBUG
	  std::ofstream of("intcvs.g2");
#endif
	  vector<shared_ptr<ParamCurve> > bd_curves;
	  for (int kd=0; kd<numshell; ++kd)
	    {
	      ftCurve intres = body->getShell(kd)->intersect(plane);

	      // Collect intersection curves and remove basedir coordinate
	      int num = intres.numSegments();
	      for (int ka=0; ka<num; ++ka)
		{
		  shared_ptr<ParamCurve> cv = intres.segment(ka).spaceCurve();
		  shared_ptr<CurveOnSurface> sfcv =
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
		  if (sfcv)
		    cv = sfcv->spaceCurve();
#ifdef DEBUG
		  if (cv.get())
		    {
		      cv->writeStandardHeader(of);
		      cv->write(of);
		    }
#endif
		  
		  shared_ptr<SplineCurve> spl2;
		  SplineCurve *spl = cv->getSplineCurve();
		  if (!spl)
		    {
		      spl2 = shared_ptr<SplineCurve>(cv->geometryCurve());
		      spl = spl2.get();
		    }

		  int dim = spl->dimension();
		  bool rat = spl->rational();
		  int numc = spl->numCoefs();
		  vector<double> coefs;
		  vector<double>::iterator cf = spl->coefs_begin();
		  for (int kb=0; kb<numc; ++kb)
		    {
		      for (int kc=0; kc<dim; ++kc, cf++)
			{
			  if (kc == basedir)
			    continue;
			  coefs.push_back(*cf);
			  if (rat)
			    {
			      cf++;
			      coefs.push_back(*cf);
			    }
			}
		    }
		  shared_ptr<SplineCurve> bdspl(new SplineCurve(spl->basis(), &coefs[0],
								dim-1, rat));
		  bd_curves.push_back(bdspl);
		}
	    }

	  // Compute 2D quadrature points
	  // NB! What about disjunct curve loops? Are they handled in 2D code?
	  CutCellQuad quad2D(bd_curves, 2.0*tol_);
	  quad2D.setQuadratureInfo(quadpar_, weights_, min_cell_size_);
	  vector<double> quadpts;
	  vector<double> wgts;
	  vector<vector<shared_ptr<ParamCurve> > > unresolved;
	  vector<double> cvquads;
	  vector<double> cvnorms;
	  vector<double> cvwgts;
	  vector<vector<shared_ptr<ParamCurve> > > shortcvs;
	  quad2D.quadrature(ll2, ur2, quadpts, wgts, unresolved, cvquads, cvnorms,
			    cvwgts, shortcvs);

	  // Expand quadrature results with basedir information
	  vector<double> currquads;
	  vector<double> currwgts;
	  size_t kr, kh;
	  for (kr=0, kh=0; kr<quadpts.size(); kr+=2, ++kh)
	    {
	      currwgts.push_back(wgts[kh]*wgt);
	      size_t kw = 0;
	      for (int ka=0; ka<3; ++ka)
		{
		  if (ka == basedir)
		    currquads.push_back(quadpar);
		  else
		    {
		      currquads.push_back(quadpts[kr+kw]);
		      ++kw;
		    }
		}
	    }
	  quadraturepoints.insert(quadraturepoints.end(), currquads.begin(),
				  currquads.end());
	  pointsweights.insert(pointsweights.end(), currwgts.begin(), currwgts.end());

#ifdef DEBUG
	  ofp << "400 1 0 4 255 0 0 255" << std::endl;
	  ofp << currwgts.size() << std::endl;
	  for (kr=0; kr<currquads.size(); kr+=3)
	    {
	      Point pt(currquads[kr], currquads[kr+1], currquads[kr+2]);
	      ofp << pt << std::endl;
	    }
#endif
	}
    
      // Compute surface curvature
#ifdef DEBUG
      std::ofstream sfof("sfquads.g2");
#endif
      for (int ka=0; ka<numshell; ++ka)
	{
	  shared_ptr<SurfaceModel> shell = body->getShell(ka);
	  int num = shell->nmbEntities();
	  for (int kb=0; kb<num; ++kb)
	    {
	      shared_ptr<ParamSurface> surf = shell->getSurface(kb);
	      int code = surf->getFlagCode();
	      if (code != 2)
		continue;   // Not a boundary surface

	      // Fetch 2D boundary loop
	      vector<shared_ptr<ParamCurve> > bdcvs;
	      vector<CurveLoop> loop = SurfaceTools::absolutelyAllBoundarySfLoops(surf, tol_);
	      int nmbcvs = loop[0].size();
	      for (int kc=0; kc<nmbcvs; ++kc)
		{
		  shared_ptr<ParamCurve> parcv = loop[0][kc];
		  shared_ptr<CurveOnSurface> sfcv =
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(parcv);
		  if (!sfcv.get())
		    THROW("Not a curve-on-surface curve");
		  sfcv->ensureParCrvExistence(tol_);
		  bdcvs.push_back(sfcv->parameterCurve());
		}

	      if (bdcvs.size() == 0)
		THROW("Missing 2D curves");
	      BoundingBox bdbox = bdcvs[0]->boundingBox();
	      for (size_t ki=1; ki<bdcvs.size(); ++ki)
		{
		  BoundingBox bdbox2 = bdcvs[ki]->boundingBox();
		  bdbox.addUnionWith(bdbox2);
		}

	      // Define 2D quadrature points
	      CutCellQuad quad2D(bdcvs, 2.0*tol_);
	      quad2D.setQuadratureInfo(quadpar_, weights_, min_cell_size_);
	      vector<double> quadpts;
	      vector<double> wgts;
	      vector<vector<shared_ptr<ParamCurve> > > unresolved;
	      vector<double> cvquads;
	      vector<double> cvnorms;
	      vector<double> cvwgts;
	      vector<vector<shared_ptr<ParamCurve> > > shortcvs;
	      quad2D.quadrature(bdbox.low(), bdbox.high(), quadpts, wgts, unresolved, cvquads, 
				cvnorms, cvwgts, shortcvs);

	      for (size_t ki=0; ki<wgts.size(); ++ki)
		{
		  vector<Point> der(3);
		  surf->point(der, quadpts[2*ki], quadpts[2*ki+1], 1);
		  surfquads.insert(surfquads.end(), der[0].begin(), der[0].end());
		  double det = der[1].length()*der[2].length();
		  sfptweights.push_back(wgts[ki]*det);
		  Point norm = der[1].cross(der[2]);
		  (void)norm.normalize_checked();
		  surfnorms.insert(surfnorms.end(), norm.begin(), norm.end());
		}

#ifdef DEBUG
	      for (size_t ki=0; ki<wgts.size(); ++ki)
		{
		  sfof << "400 1 0 4 155 0 100 255" << std::endl;
		  sfof << wgts.size() << std::endl;
		  for (size_t ki=0; ki<wgts.size(); ++ki)
		    {
		      Point pos = surf->point(quadpts[2*ki], quadpts[2*ki+1]);
		      sfof << pos << std::endl;
		    }
		}
		  
	      int stop_bread_bd = 1;
#endif
	    }
	}
      
      int stop_break2 = 1;
    }
}

struct edgeInfo
{
  int splitdir_;
  double splitval_;
  int prio_;

  edgeInfo()
  {
    splitdir_ = -1;
    prio_ = -1;
  }

  edgeInfo(int splitdir, double splitval, int prio)
  {
    splitdir_ = splitdir;
    splitval_ = splitval;
    prio_ = prio;
  }
  
};

struct ruledInfo
{
  int ruleddir_;
  int sgn_;
  int sfix1_;
  int sfix2_;
  shared_ptr<SplineCurve> crv_;

  ruledInfo()
  {
    ruleddir_ = -1;
    sgn_ = 1;
    sfix1_ = sfix2_ = -1;
  }
  
  ruledInfo(int ruleddir, int sgn, int sfix1, int sfix2, shared_ptr<SplineCurve> crv)
  {
    ruleddir_ = ruleddir;
    sgn_ = sgn;
    sfix1_ = sfix1;
    sfix2_ = sfix2;
    crv_ = crv;
  }

};
  
//==============================================================================
int
CutCellQuad3D::selectBaseDir(Point& ll, Point& ur,
			     shared_ptr<Body> body,
			     int& splitdir, double& splitval,
			     shared_ptr<SplineCurve>& rulecv,
			     Point& ruledir)
//==============================================================================
{
  splitdir = -1;
  splitval = 0.0;
  
  // Collect surfaces and generate normal cones corresponding to all surfaces
  double numtol = 1.0e-10;
  double angtol = 0.001;
  double anglim = M_PI/6.0;
  int numshells = body->nmbOfShells();
  int numfaces = body->nmbOfFaces();
  vector<shared_ptr<ParamSurface> > sfs(numfaces);
  vector<DirectionCone> ncones(numfaces);
  int curved[3];
  vector<pair<double,double> > minmaxang(3*numfaces);
  for (int kc=0; kc<3; ++kc)
    {
      curved[kc] = 0;
    }
  int ix = 0;
  for (int ka=0; ka<numshells; ++ka)
    {
      shared_ptr<SurfaceModel> shell = body->getShell(ka);
      int num = shell->nmbEntities();
      for (int kb=0; kb<num; ++kb)
	{
	  sfs[ix] = shell->getSurface(kb);
	  ncones[ix] = sfs[ix]->normalCone();
	  pair<double,double> dirang[3];
	  sfs[ix]->dirNormAngles(dirang);
	  for (int kc=0; kc<3; ++kc)
	    {
	      minmaxang[3*ix+kc] = dirang[kc];
	      if (dirang[kc].second-dirang[kc].first > angtol)
		curved[kc]++;
	    }
	  ++ix;
	}
    }
#ifdef DEBUG
  std::ofstream of("selsfs.g2");
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      sfs[ki]->writeStandardHeader(of);
      sfs[ki]->write(of);
    }
#endif
  
  // Check the coordinate directions for a suitable base direction
  // Identify surfaces that overlap with each coordinate direction
  vector<size_t> cand[3];
  int nonlin[3];
  for (ix=0; ix<3; ++ix)
    {
      nonlin[ix] = 0;
      Point dir(0.0, 0.0, 0.0);
      dir[ix] = 1.0;
      Point dir2 = -dir;

      for (size_t ki=0; ki<ncones.size(); ++ki)
	{
	  // Check normal cone overlap
	  // bool overlap1 = ncones[ki].containsDirection(dir);
	  // bool overlap2 = ncones[ki].containsDirection(dir2);
	  bool overlap = (fabs(ncones[ki].centre()*dir) > numtol);
	  if (overlap)
	    {
	      cand[ix].push_back(ki);
	      if (ncones[ki].angle() > angtol)
		nonlin[ix]++;
	    }
	}

    }

 BoundingBox bbox = body->boundingBox();
 for (int ix=0; ix<3; ++ix)
   {
     ll[ix] = std::max(ll[ix], bbox.low()[ix]);
     ur[ix] = std::min(ur[ix], bbox.high()[ix]);
   }

  // Fetch sharp edges
 vector<ftEdge*> convex, concave, all_cvs;
 fetchSharpEdges(body, convex, concave);
 all_cvs.insert(all_cvs.end(), convex.begin(), convex.end());
 all_cvs.insert(all_cvs.end(), concave.begin(), concave.end());
 
#ifdef DEBUG
  std::ofstream ofconv("convex.g2");
  std::ofstream ofconc("concave.g2");
  for (size_t ki=0; ki<convex.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = convex[ki]->geomCurve();
      SplineCurve *spl = cv->getSplineCurve();
      if (spl)
	{
	  spl->writeStandardHeader(ofconv);
	  spl->write(ofconv);
	}
    }
  for (size_t ki=0; ki<concave.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = concave[ki]->geomCurve();
      SplineCurve *spl = cv->getSplineCurve();
      if (spl)
	{
	  spl->writeStandardHeader(ofconc);
	  spl->write(ofconc);
	}
    }
#endif
  
  if (nonlin[0] == 0 && nonlin[1] == 0 && nonlin[2] == 0 &&
      concave.size() == 0)
   {
     // All boundary surfaces are planar
     for (ix=0; ix<3; ++ix)
       {
  	 Point dir(0.0, 0.0, 0.0);
  	 dir[ix] = 1.0;
  	 Point dir2 = -dir;
  	 int nmb_dir = 0;
  	 for (size_t ki=0; ki<cand[ix].size(); ++ki)
  	   {
  	     if (ncones[cand[ix][ki]].containsDirection(dir) ||
  		 ncones[cand[ix][ki]].containsDirection(dir2))
  	       nmb_dir++;
  	   }
       if (nmb_dir == 2)
  	 return ix;  
       }
   }

  // Check for a simple configuration that can be handled by the 2D code
  for (ix=0; ix<3; ++ix)
    {
      if (nonlin[ix] > 0)
  	continue;
      
      // Check if the side surface normal coincides with the coordinate directions
      Point dir(0.0, 0.0, 0.0);
      dir[ix] = 1.0;
      Point dir2 = -dir;
      size_t ki;
      for (ki=0; ki<cand[ix].size(); ++ki)
  	{
  	  if (dir.angle(ncones[cand[ix][ki]].centre()) > angtol &&
  	      dir2.angle(ncones[cand[ix][ki]].centre()) > angtol)
  	    break;
  	}
      if (ki == cand[ix].size())
  	return ix;   // The configuration corresponds to a linear sweep
    }

  // A further search for a simple parameter direction or planar splits
  int numplane[3];
  vector<vector<ftEdge*> > planar(3);
  vector<vector<double> > planar_par(3);
  vector<vector<double> > sharp_bounds(3);
  double acc_coneangle[6];
  vector<int> cand_sidesfs;
  for (ix=0; ix<3; ++ix)
    {
      // Count number of planar surfaces with a surface normal coinciding with the
      // coordinate axis
      numplane[ix] = 0;
      Point dir(0.0, 0.0, 0.0);
      dir[ix] = 1.0;
      Point dir2 = -dir;
      for (size_t ki=0; ki<cand[ix].size(); ++ki)
	{
	  if (ncones[cand[ix][ki]].angle() <= angtol &&
	      (dir.angle(ncones[cand[ix][ki]].centre()) <= angtol ||
	       dir2.angle(ncones[cand[ix][ki]].centre()) <= angtol))
	    {
	      ++numplane[ix];
	      cand_sidesfs.push_back(cand[ix][ki]);
	    }
	}

      // Identify sharp curves lying in an internal plane in this parameter direction
      double low = ll[ix];
      double high = ur[ix];
      planar[ix].insert(planar[ix].end(), convex.begin(), convex.end());
      planar[ix].insert(planar[ix].end(), concave.begin(), concave.end());
      for (size_t ki=0; ki<planar[ix].size(); )
	{
	  shared_ptr<ParamCurve> crv = planar[ix][ki]->geomCurve();
	  BoundingBox bbcv = crv->boundingBox();
	  double low2 = bbcv.low()[ix];
	  double high2 = bbcv.high()[ix];
	  if (high2 - low2 > tol_)
	    {
	      // Not planar, remember end positions of curve
	      vector<Point> der1(2), der2(2);
	      crv->point(der1, planar[ix][ki]->tMin(), 1);
	      crv->point(der2, planar[ix][ki]->tMax(), 1);
	      if (der1[0][ix] > low+tol_ && der1[0][ix] < high-tol_)
		{
		  // Check for a smooth transisition
		  size_t kj;
		  for (kj=0; kj<all_cvs.size(); ++kj)
		    {
		      if (all_cvs[kj] == planar[ix][ki])
			continue;
		      shared_ptr<ParamCurve> crv2 = all_cvs[kj]->geomCurve();
		      vector<Point> der3(2), der4(2);
		      crv2->point(der3, all_cvs[kj]->tMin(), 1);
		      crv2->point(der4, all_cvs[kj]->tMax(), 1);

		      double d1 = der1[0].dist(der3[0]);
		      double d2 = der1[0].dist(der4[0]);
		      double a1 = der1[1].angle(der3[1]);
		      a1 = std::min(a1, M_PI-a1);
		      double a2 = der1[1].angle(der4[1]);
		      a2 = std::min(a2, M_PI-a2);
		      if (d1 < tol_ && a1 < angtol)
			break;
		      if (d2 < tol_ && a2 < angtol)
			break;
		    }
		  if (kj == all_cvs.size())
		    sharp_bounds[ix].push_back(der1[0][ix]);
		}
	      if (der2[0][ix] > low+tol_ && der2[0][ix] < high-tol_)
		{
		  // Check for a smooth transisition
		  size_t kj;
		  for (kj=0; kj<all_cvs.size(); ++kj)
		    {
		      if (all_cvs[kj] == planar[ix][ki])
			continue;
		      shared_ptr<ParamCurve> crv2 = all_cvs[kj]->geomCurve();
		      vector<Point> der3(2), der4(2);
		      crv2->point(der3, all_cvs[kj]->tMin(), 1);
		      crv2->point(der4, all_cvs[kj]->tMax(), 1);

		      double d1 = der2[0].dist(der3[0]);
		      double d2 = der2[0].dist(der4[0]);
		      double a1 = der2[1].angle(der3[1]);
		      a1 = std::min(a1, M_PI-a1);
		      double a2 = der2[1].angle(der4[1]);
		      a2 = std::min(a2, M_PI-a2);
		      if (d1 < tol_ && a1 < angtol)
			break;
		      if (d2 < tol_ && a2 < angtol)
			break;
		    }
		  if (kj == all_cvs.size())
		    sharp_bounds[ix].push_back(der2[0][ix]);
		}

	      planar[ix].erase(planar[ix].begin() + ki);
	    }
	  else if (fabs(high2-low) < tol_ || fabs(high-low2) < tol_)
	    planar[ix].erase(planar[ix].begin() + ki);
	  else
	    {
	      // Internal planar curve
	      planar_par[ix].push_back(0.5*(low2+high2));
	      ++ki;
	    }
	}

      // Compute accumulated normal cone angle along candidate direction of remaining surfaces
      acc_coneangle[2*ix] = std::numeric_limits<double>::max();
      acc_coneangle[2*ix+1] = std::numeric_limits<double>::lowest();
      for (int ix2=0; ix2<3; ++ix2)
	{
	  if (ix2 == ix)
	    continue;

	  for (size_t ki=0; ki<cand[ix2].size(); ++ki)
	    {
	      int sfix = cand[ix2][ki];
	      acc_coneangle[2*ix] = std::min(acc_coneangle[2*ix],
					     minmaxang[3*sfix+ix].first);
	      acc_coneangle[2*ix+1] = std::max(acc_coneangle[2*ix+1],
					       minmaxang[3*sfix+ix].second);
	    }
	}
    }

  // Clean up in sharp_bounds
  for (ix=0; ix<3; ++ix)
    {
      if (sharp_bounds[ix].size() == 0)
	continue;

      std::sort(sharp_bounds[ix].begin(), sharp_bounds[ix].end());
      size_t ki, kj;
      for (ki=0; ki<sharp_bounds[ix].size(); ki=kj)
	{
	  double mean = sharp_bounds[ix][ki];
	  for (kj=ki+1; kj<sharp_bounds[ix].size(); ++kj)
	    {
	      if (sharp_bounds[ix][kj] > sharp_bounds[ix][ki]+tol_)
		break;
	      mean += sharp_bounds[ix][kj];
	    }
	  if (kj - ki > 1)
	    {
	      mean /= (double)(kj-ki);
	      sharp_bounds[ix][ki] = mean;
	      sharp_bounds[ix].erase(sharp_bounds[ix].begin()+ki+1,
				     sharp_bounds[ix].begin()+kj);
	      kj = ki+1;
	    }
	}
      
      for (ki=0; ki<sharp_bounds[ix].size(); )
	{
	  for (kj=0; kj<planar_par[ix].size(); ++kj)
	    {
	      if (fabs(sharp_bounds[ix][ki]-planar_par[ix][kj]) < tol_)
		{
		  sharp_bounds[ix].erase(sharp_bounds[ix].begin()+ki,
					 sharp_bounds[ix].begin()+ki+1);
		  break;
		}
	    }
	  if (kj == planar_par[ix].size())
	    ++ki;
	}
    }
  
#ifdef DEBUG
 std::ofstream edgof("planar_edgs.g2");
 for (int ix=0; ix<3; ++ix)
   {
     for (size_t ki=0; ki<planar[ix].size(); ++ki)
       {
	 shared_ptr<ParamCurve> crv = planar[ix][ki]->geomCurve();
	 if (crv.get())
	   {
	     shared_ptr<SplineCurve> crv2(crv->geometryCurve());
	     crv2->writeStandardHeader(edgof);
	     crv2->write(edgof);
	   }
       }
   }
#endif

 // Check for a base direction
  for (ix=0; ix<3; ++ix)
    {
      if (planar_par[ix].size() + sharp_bounds[ix].size() > 0)
	continue;   // Edges in potential direction

      if (acc_coneangle[2*ix+1] - acc_coneangle[2*ix] > anglim)
	continue; // Too high curvature

      return ix;  // Direction found
    }
  

  // Collect splitting information
  vector<edgeInfo> split_info;
  vector<ruledInfo> ruled_info;
  for (size_t ki=0; ki<concave.size(); ++ki)
    {
      // Identify associated surfaces
      ftEdgeBase *twin = concave[ki]->twin();
      if (!twin)
	continue;
     shared_ptr<ParamSurface> surf1 = concave[ki]->face()->asFtSurface()->surface();
      shared_ptr<ParamSurface> surf2 = twin->face()->asFtSurface()->surface();
      int sfix1=-1, sfix2=-1;
      for (size_t kj=0; kj<sfs.size(); ++kj)
	{
	  if (sfs[kj].get() == surf1.get())
	    sfix1 = (int)kj;
	  else if (sfs[kj].get() == surf2.get())
	    sfix2 = (int)kj;
	}

      // Find ruled direction
      double minang = M_PI;
      int minix = -1;
      int sgn = 0;
      for (int ix2=0; ix2<3; ++ix2)
	{
	  Point dir(0.0, 0.0, 0.0);
	  dir[ix2] = 1.0;
	  double ang1 = dir.angle(ncones[sfix1].centre());
	  double ang2 = dir.angle(ncones[sfix2].centre());
	  if (ang1 < minang || M_PI-ang1 < minang)
	    {
	      minang = ang1;
	      minix = ix2;
	      sgn = (ang1 < M_PI-ang1) ? -1 : 1;
	    }
	  if (ang2 < minang || M_PI-ang2 < minang)
	    {
	      minang = ang2;
	      minix = ix2;
	      sgn = (ang2 < M_PI-ang2) ? -1 : 1;
	    }
	}
#ifdef DEBUG
      std::cout << "Ruled direction: " << minix << std::endl;
#endif
      shared_ptr<ParamCurve> crv = concave[ki]->geomCurve();
      shared_ptr<SplineCurve> rulecv(crv->geometryCurve());
      shared_ptr<SplineCurve> subrule(rulecv->subCurve(concave[ki]->tMin(),
						       concave[ki]->tMax()));
      ruledInfo rinfo(minix, sgn, sfix1, sfix2, subrule);
      ruled_info.push_back(rinfo);
    }

  for (ix=0; ix<3; ++ix)
    {
      for (size_t ki=0; ki<planar_par[ix].size(); ++ki)
	{
	  // Check if the edge is associated a candidate side surface
	  // Identify associated surfaces
	  ftEdgeBase *twin = planar[ix][ki]->twin();
	  if (!twin)
	    continue;
	  shared_ptr<ParamSurface> surf1 = planar[ix][ki]->face()->asFtSurface()->surface();
	  shared_ptr<ParamSurface> surf2 = twin->face()->asFtSurface()->surface();
	  int sfix1=-1, sfix2=-1;
	  size_t kj;
	  for (kj=0; kj<sfs.size(); ++kj)
	    {
	      if (sfs[kj].get() == surf1.get())
		sfix1 = (int)kj;
	      else if (sfs[kj].get() == surf2.get())
		sfix2 = (int)kj;
	    }

	  for (kj=0; kj<cand_sidesfs.size(); ++kj)
	    if (sfix1 == cand_sidesfs[kj] || sfix2 == cand_sidesfs[kj])
	      break;
	  
	  int prio = (kj == cand_sidesfs.size()) ? 2 : 0;
	  

	  edgeInfo info(ix, planar_par[ix][ki], prio);
	  split_info.push_back(info);
	}

      for (size_t ki=0; ki<sharp_bounds[ix].size(); ++ki)
	{
	  edgeInfo info(ix, sharp_bounds[ix][ki], 1);
	  split_info.push_back(info);
	}
    }
  
 if (ruled_info.size() > 0)
   {
     // Should merge curves if possible or select best curve if more than one exists.
     // Simple first approach
     rulecv = ruled_info[0].crv_;
     ruledir = Point(0.0, 0.0, 0.0);
     ruledir[ruled_info[0].ruleddir_] = 1.0;
     ruledir *= ruled_info[0].sgn_;
     return -1;
   }
 
  
 if (split_info.size() > 0)
   {
     // Select the middlemost candidate, prefer planar curves
     double maxfrac = 0.0;
     size_t maxix = -1;
     int maxprio = -1;
     for (size_t ki=0; ki<split_info.size(); ++ki)
       {
	 int dir = split_info[ki].splitdir_;
	 double val = split_info[ki].splitval_;
	 double frac = std::min(val-ll[dir], ur[dir]-val)/(ur[dir]-ll[dir]);
	 if (split_info[ki].prio_ > maxprio ||
	     (split_info[ki].prio_ == maxprio && frac > maxfrac))
	 {
	   maxprio = split_info[ki].prio_;
	   maxfrac = frac;
	   maxix = ki;
	 }
       }

     splitdir = split_info[maxix].splitdir_;
     splitval = split_info[maxix].splitval_;
   }

  return -1;  // No suitable direction found
}

//==============================================================================
void
CutCellQuad3D::fetchSharpEdges(shared_ptr<Body> body, vector<ftEdge*>& convex,
			       vector<ftEdge*>& concave)
//==============================================================================
{
  int numshells = body->nmbOfShells();
  int nmb_samples = 10;
  for (int ka=0; ka<numshells; ++ka)
    {
      shared_ptr<SurfaceModel> shell = body->getShell(ka);

      // Collect all sharp edges
      vector<ftEdge*> sharp_edges;
      shell->getCorners(sharp_edges);

      // Sort into convex and concave edges (majority vote)
      for (size_t ki=0; ki<sharp_edges.size(); ++ki)
	{
	  if (!sharp_edges[ki]->twin())
	    {
	      sharp_edges.erase(sharp_edges.begin()+ki);
	      continue;  // Gaps are not expected and not treated
	    }

	  int nmb_concave = 0;
	  int nmb_convex = 0;
	  double t1 = sharp_edges[ki]->tMin();
	  double t2 = sharp_edges[ki]->tMax();
	  double del = (t2 - t1)/(double)(nmb_samples-1);
	  double par;
	  int kr;
	  for (kr=0, par=t1; kr<nmb_samples; ++kr, par+=del)
	    {
	      Point pos1 = sharp_edges[ki]->point(par);
	      Point norm1 = sharp_edges[ki]->normal(par);
	      Point tan1 = sharp_edges[ki]->tangent(par);
	      Point vec1 = tan1.cross(norm1);
	      
	      double clo_t, clo_dist;
	      Point clo_pt;
	      sharp_edges[ki]->twin()->closestPoint(pos1, clo_t, clo_pt,
						    clo_dist);
	      Point norm2 = sharp_edges[ki]->twin()->normal(clo_t);
	      Point tan2 = sharp_edges[ki]->twin()->tangent(clo_t);
	      Point vec2 = tan2.cross(norm2);
	      double ang2 = vec1.angle(norm2);
	      double ang3 = vec2.angle(norm1);
	      if (ang2 + ang3 > M_PI)
		nmb_concave++;
	      else
		nmb_convex++;
	      int stop_break;
	    }
	  if (nmb_convex >= nmb_concave)
	    convex.push_back(sharp_edges[ki]);
	  else
	    concave.push_back(sharp_edges[ki]);
	}
     
    }
}  



