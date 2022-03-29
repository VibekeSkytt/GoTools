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

#include "GoTools/creators/CutCellQuad.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/GoIntersections.h"
#include <iostream>
#include <fstream>

//#define DEBUG

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

//==============================================================================
CutCellQuad::CutCellQuad(vector<shared_ptr<ParamCurve> >& bd_curves,
			 double tol)
//==============================================================================
  : tol_(tol), angtol_(1.0e-2)
{
  //vector<shared_ptr<CurveLoop> > loops;
  sortBoundary(bd_curves, loops_);

  //domain_ = shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(loops));
}

//==============================================================================
int CutCellQuad::cellStat(const Point& ll, const Point& ur)
//==============================================================================
{
  // Create cell boundaries
  vector<shared_ptr<SplineCurve> > cell_bd(4);
  cell_bd[0] = shared_ptr<SplineCurve>(new SplineCurve(ll, Point(ur[0],ll[1])));
  cell_bd[1] = shared_ptr<SplineCurve>(new SplineCurve(Point(ur[0],ll[1]), ur));
  cell_bd[2] = shared_ptr<SplineCurve>(new SplineCurve(ur, Point(ll[0],ur[1])));
  cell_bd[3] = shared_ptr<SplineCurve>(new SplineCurve(Point(ll[0],ur[1]), ll));

  // Check for intersections between the boundary curves and the curve model
  double eps = 1.0e-10;
  shared_ptr<CurveBoundedDomain> domain(new CurveBoundedDomain(loops_[0]));
  int nmb_outside = 0;
  for (size_t ki=0; ki<4; ++ki)
    {
      vector<double> start_end_par;
      domain->findPcurveInsideSegments(*cell_bd[ki], tol_, start_end_par, false);
      if (start_end_par.size() > 2)
	return 2;

      if (start_end_par.size() == 0)
	nmb_outside++;

      double t1 = cell_bd[ki]->startparam();
      double t2 = cell_bd[ki]->endparam();
      for (size_t kj=0; kj<start_end_par.size(); ++kj)
	if (start_end_par[kj] > t1 + eps && start_end_par[kj] < t2 - eps)
	  return 2;
    }
  if (nmb_outside > 0 && nmb_outside < 4)
    return 2;

  // Inside/outside test
  Vector2D mid(0.5*(ll[0]+ur[0]), 0.5*(ll[1]+ur[1]));
  bool inside = domain->isInDomain(mid, tol_);

  if (inside)
    {
      RectDomain celldom(Vector2D(ll[0], ll[1]), Vector2D(ur[0], ur[1]));

      // Check for intersection with trimming loop
      for (size_t ki=1; ki<loops_.size(); ++ki)
	{
	  shared_ptr<CurveBoundedDomain> domain2(new CurveBoundedDomain(loops_[ki]));
	  for (size_t ki=0; ki<4; ++ki)
	    {
	      vector<double> start_end_par;
	      domain2->findPcurveInsideSegments(*cell_bd[ki], tol_, start_end_par, false);
	      if (start_end_par.size() > 2)
		return 2;
	      
	      if (start_end_par.size() == 0)
		nmb_outside++;

	      double t1 = cell_bd[ki]->startparam();
	      double t2 = cell_bd[ki]->endparam();
	      for (size_t kj=0; kj<start_end_par.size(); ++kj)
		if (start_end_par[kj] > t1 + eps && start_end_par[kj] < t2 - eps)
		  return 2;
	    }
	  if (nmb_outside > 0 && nmb_outside < 4)
	    return 2;
	}

      // Check for totally included trimming loop
     for (size_t ki=1; ki<loops_.size(); ++ki)
	{
	  shared_ptr<ParamCurve> cv = (*loops_[ki])[0];
	  Point pos = cv->point(cv->startparam());
	  if (celldom.isInDomain(Vector2D(pos[0],pos[1]), tol_))
	    return 3;
	}
    }
  
  return (inside) ? 1 : 0;
}

//==============================================================================
void
CutCellQuad::quadrature(const Point& ll, const Point& ur,
			vector<double>& quadraturepoints,
			vector<double>& pointsweights,
			vector<vector<shared_ptr<ParamCurve> > >& unresolved_cells,
			vector<double>& curvequads,
			vector<double>& curvenorms,
			vector<double>& crvptweights,
			vector<vector<shared_ptr<ParamCurve> > >& short_curves,
			int stat)
//==============================================================================
{
  // Create cell boundaries
  vector<shared_ptr<SplineCurve> > cell_bd(4);
  cell_bd[0] = shared_ptr<SplineCurve>(new SplineCurve(ll, Point(ur[0],ll[1])));
  cell_bd[1] = shared_ptr<SplineCurve>(new SplineCurve(Point(ur[0],ll[1]), ur));
  cell_bd[2] = shared_ptr<SplineCurve>(new SplineCurve(ur, Point(ll[0],ur[1])));
  cell_bd[3] = shared_ptr<SplineCurve>(new SplineCurve(Point(ll[0],ur[1]), ll));

#ifdef DEBUG
  std::ofstream ofbd("cell_bd.g2");
  for (int kk=0; kk<4; ++kk)
    {
      cell_bd[kk]->writeStandardHeader(ofbd);
      cell_bd[kk]->write(ofbd);
    }
#endif
  
  shared_ptr<CurveBoundedDomain> domain(new CurveBoundedDomain(loops_[0]));
  if (stat < 0 || stat > 3)
    {
      double eps = 1.0e-10;
      int nmb_outside = 0;
      for (size_t ki=0; ki<4; ++ki)
	{
	  vector<double> start_end_par;
	  domain->findPcurveInsideSegments(*cell_bd[ki], tol_, start_end_par, false);
	  if (start_end_par.size() > 2)
	    {
	      stat = 2;
	      break;
	    }

	  if (start_end_par.size() == 0)
	    nmb_outside++;

	  double t1 = cell_bd[ki]->startparam();
	  double t2 = cell_bd[ki]->endparam();
	  for (size_t kj=0; kj<start_end_par.size(); ++kj)
	    if (start_end_par[kj] > t1 + eps && start_end_par[kj] < t2 - eps)
	      {
		stat = 2;
		break;
	      }
	  if (stat == 2)
	    break;
	}
      if (nmb_outside > 0 && nmb_outside < 4)
	stat = 2;
      
      if (stat < 3)
	{
	  Vector2D mid(0.5*(ll[0]+ur[0]), 0.5*(ll[1]+ur[1]));
	  bool inside = domain->isInDomain(mid, tol_);
	  if (inside)
	    {
	      RectDomain celldom(Vector2D(ll[0], ll[1]), Vector2D(ur[0], ur[1]));
		
	      // Check for intersection with trimming loop
	      for (size_t ki=1; ki<loops_.size(); ++ki)
		{
		  shared_ptr<CurveBoundedDomain> domain2(new CurveBoundedDomain(loops_[ki]));
		  for (size_t ki=0; ki<4; ++ki)
		    {
		      vector<double> start_end_par;
		      domain2->findPcurveInsideSegments(*cell_bd[ki], tol_, start_end_par, false);
		      if (start_end_par.size() > 2)
			stat = 2;
	      
		      if (start_end_par.size() == 0)
			nmb_outside++;
		      
		      double t1 = cell_bd[ki]->startparam();
		      double t2 = cell_bd[ki]->endparam();
		      for (size_t kj=0; kj<start_end_par.size(); ++kj)
			if (start_end_par[kj] > t1 + eps && start_end_par[kj] < t2 - eps)
			  stat = 2;
		    }
		  if (nmb_outside > 0 && nmb_outside < 4)
		    stat = 2;

		  if (stat == 2)
		    break;
		}
	      
	      if (stat != 2)
		{
		  // Check for totally included trimming loop
		  for (size_t ki=1; ki<loops_.size(); ++ki)
		    {
		      shared_ptr<ParamCurve> cv = (*loops_[ki])[0];
		      Point pos = cv->point(cv->startparam());
		      if (celldom.isInDomain(Vector2D(pos[0],pos[1]), tol_))
			stat = 3;
		    }
		}
	    }
	  if (stat < 2)
	    stat = (inside) ? 1 : 0;
	}
    }

  if (stat == 0)
    return;  // Outside, no quadrature points

  else if (stat == 1)
    {
      // Not a cut cell. Adapt quadrature points to cell size
      vector<double> quadpar1(quadpar_.size());
      vector<double> quadpar2(quadpar_.size());
      vector<double> wgt1(quadpar_.size());
      vector<double> wgt2(quadpar_.size());
      double del1 = ur[0] - ll[0];
      double del2 = ur[1] - ll[1];
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
      quadraturepoints.resize(2*quadpar_.size()*quadpar_.size());
      pointsweights.resize(quadpar_.size()*quadpar_.size());
      for (size_t kj=0; kj<quadpar_.size(); ++kj)
	for (size_t ki=0; ki<quadpar_.size(); ++ki)
	  {
	    quadraturepoints[2*(kj*quadpar_.size()+ki)] = quadpar1[ki];
	    quadraturepoints[2*(kj*quadpar_.size()+ki)+1] = quadpar2[kj];
	    pointsweights[kj*quadpar_.size()+ki] = wgt1[ki]*wgt2[kj];
	  }

      // Check for coincidence between a piece of the boundary loop and the
      // cell boundaries
      vector<shared_ptr<ParamCurve> > cell_bd2(cell_bd.begin(), cell_bd.end());
      shared_ptr<CurveLoop> cell_loop(new CurveLoop(cell_bd2, tol_));
      shared_ptr<CurveBoundedDomain> cvdom(new CurveBoundedDomain(cell_loop));
      vector<vector<shared_ptr<ParamCurve> > > bd_segs;
      identifyBdCvs(loops_, cvdom, bd_segs);
      for (size_t ki=0; ki<bd_segs.size(); ++ki)
	{
	  // Compute curve length
	  double len = 0.0;
	  for (size_t kj=0; kj<bd_segs[ki].size(); ++kj)
	    len += bd_segs[ki][kj]->length(tol_);
	  if (len < min_cell_size_)
	    short_curves.push_back(bd_segs[ki]);
	  else
	    {
	      computeQuadraturePoints(bd_segs[ki], curvequads, curvenorms, crvptweights);
	    }
	}
    }
  else
    {
      // Cut cell
#ifdef DEBUG
      std::cout << "Cut cell" << std::endl;
#endif
      vector<vector<shared_ptr<CurveLoop> > > cut_bd;
      vector<shared_ptr<ParamCurve> > cell_bd2(cell_bd.begin(), cell_bd.end());
      shared_ptr<CurveLoop> cell_loop(new CurveLoop(cell_bd2, tol_));
      cutcelldom_ = shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(cell_loop));
	  
      if (stat == 2)
	{
	  // Create cut cell representation
	  createCutCell(loops_, cell_bd, cut_bd, true);
	}
      else
	{
	  cut_bd.resize(1);
	  cut_bd[0].push_back(cell_loop);

	  // Identify internal trimming loops inside the cell
	  RectDomain celldom(Vector2D(ll[0], ll[1]), Vector2D(ur[0], ur[1]));
	  for (size_t ki=1; ki<loops_.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = (*loops_[ki])[0];
	      Point pos = cv->point(cv->startparam());
	      if (celldom.isInDomain(Vector2D(pos[0],pos[1]), tol_))
		cut_bd[0].push_back(loops_[ki]);
	    }
	}

      for (size_t ki=0; ki<cut_bd.size(); ++ki)
	{
	  quadraturePoints(cut_bd[ki], quadraturepoints, pointsweights, unresolved_cells,
			   curvequads, curvenorms, crvptweights, short_curves);
	}
    }
  
#ifdef DEBUG
  std::ofstream ofq("quadpts.g2");
  ofq << "400 1 0 4 100 100 55 255" << std::endl;
  ofq << pointsweights.size() << std::endl;
  for (size_t kj=0; kj<quadraturepoints.size(); kj+=2)
    {
      Point pt(quadraturepoints[kj], quadraturepoints[kj+1], 0.0);
      ofq << pt << std::endl;
    }
  int stop_break = 1;
#endif
}


//==============================================================================
void CutCellQuad::sortBoundary(vector<shared_ptr<ParamCurve> >& bd_curves,
			       vector<shared_ptr<CurveLoop> >& bd_loops)
//==============================================================================
{
  if (bd_curves.size() == 0)
    return;
  
  // Compute end points of curves and ensure that the curves are in 2D
  vector<Point> endpt(2*bd_curves.size());
  for (size_t ki=0; ki<bd_curves.size(); ++ki)
    {
      if (bd_curves[ki]->dimension() != 2)
	THROW("Curve dimension different from 2");

      endpt[2*ki] = bd_curves[ki]->point(bd_curves[ki]->startparam());
      endpt[2*ki+1] = bd_curves[ki]->point(bd_curves[ki]->endparam());
      Point mid = bd_curves[ki]->point(0.5*(bd_curves[ki]->startparam()+
					    bd_curves[ki]->endparam()));
      if (endpt[2*ki].dist(mid) + endpt[2*ki+1].dist(mid) <= tol_)
	THROW("Curve too short");
    }

  // Sort curves into loops and turn direction of curves if necessary.
  // Two curves should meet at every end point (possibly the two ends of
  // a closed curve).
  double tol = tol_; //1.5*tol_;
  size_t ix = 0;
  vector<vector<shared_ptr<ParamCurve> > > loop_cvs;
  vector<Point> endloop;
  for (size_t ki=0; ki<bd_curves.size(); ++ki)
    {
      size_t kj;
      if (endpt[2*ki].dist(endpt[2*ki+1]) <= tol)
	kj = loop_cvs.size();  // Closed curve
      else
	{
	  for (kj=0; kj<loop_cvs.size(); ++kj)
	    {
	      if (endloop[2*kj].dist(endloop[2*kj+1]) <= tol)
		continue;  // Loop already closed
	  
	      double d1 = endloop[2*kj].dist(endpt[2*ki]);
	      double d2 = endloop[2*kj].dist(endpt[2*ki+1]);
	      double d3 = endloop[2*kj+1].dist(endpt[2*ki]);
	      double d4 = endloop[2*kj+1].dist(endpt[2*ki+1]);
	      if (d1 <= tol || d2 <= tol)
		{
		  if (d1 < d2)
		    {
		      bd_curves[ki]->reverseParameterDirection();
		      endloop[2*kj] = endpt[2*ki+1];
		    }
		  else
		    endloop[2*kj] = endpt[2*ki];
		  loop_cvs[kj].insert(loop_cvs[kj].begin(), bd_curves[ki]);
		  break;
		}
	      else if (d3 <= tol || d4 <= tol)
		{
		  if (d4 < d3)
		    {
		      bd_curves[ki]->reverseParameterDirection();
		      endloop[2*kj+1] = endpt[2*ki];
		    }
		  else
		    endloop[2*kj+1] = endpt[2*ki+1];
		  loop_cvs[kj].push_back(bd_curves[ki]);
		  break;
		}
	    }
		
	}
      if (kj == loop_cvs.size())
	{
	  vector<shared_ptr<ParamCurve> > curr_loop;
	  curr_loop.push_back(bd_curves[ki]);
	  loop_cvs.push_back(curr_loop);
	  endloop.push_back(endpt[2*ki]);
	  endloop.push_back(endpt[2*ki+1]);
	}
    }

  // Check if all loops are closed and make sure that each loop contains at least
  // two curves. Closed single curves loops can mess up the topology computations
  for (size_t kj=0; kj<loop_cvs.size(); ++kj)
    {
      if (endloop[2*kj].dist(endloop[2*kj+1]) > tol)
	THROW("Boundary loop not closed");

      if (loop_cvs[kj].size() == 1)
	{
	  double par = 0.5*(loop_cvs[kj][0]->startparam()+loop_cvs[kj][0]->endparam());
	  double partol = 0.1*(loop_cvs[kj][0]->startparam()+loop_cvs[kj][0]->endparam());
	  vector<shared_ptr<ParamCurve> > sub_cvs =
	    loop_cvs[kj][0]->split(par, partol);
	  loop_cvs[kj][0] = sub_cvs[0]; 
	  loop_cvs[kj].push_back(sub_cvs[1]);
	}
    }
 
  // Ensure correct orientation of loops. First check orientation
  // Create curve loops and domains
  double int_tol = tol_; //std::max(0.1*tol_, std::min(1.0e-4, tol_));
  vector<bool> CCW(loop_cvs.size());
  bd_loops.resize(loop_cvs.size());
  vector<shared_ptr<CurveBoundedDomain> > cvdom(loop_cvs.size());
  vector<Vector2D> internalPt(loop_cvs.size());
  for (size_t kj=0; kj<loop_cvs.size(); ++kj)
    {
      CCW[kj] = LoopUtils::loopIsCCW(loop_cvs[kj], tol_, int_tol);
      bd_loops[kj] = shared_ptr<CurveLoop>(new CurveLoop(loop_cvs[kj], tol_));
      cvdom[kj] =
	shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(bd_loops[kj]));

      double u, v;
      cvdom[kj]->getInternalPoint(u, v);
      internalPt[kj] = Vector2D(u, v);
    }

  // Find outer loop
  int ixout = 0;
  for (size_t ki=1; ki<loop_cvs.size(); ++ki)
    {
      bool in1 = cvdom[ixout]->isInDomain(internalPt[ki], int_tol);
      bool in2 = cvdom[ki]->isInDomain(internalPt[ixout], int_tol);
      if (in2 && (!in1))
	ixout = (int)ki;
    }

  // Make sure that outer loop has a counter clockwise orientation while
  // the remaining have a clockwise orientation
  // NB! This will fail if the input curves has several levels of curve
  // inside curve
  if (!CCW[ixout])
    bd_loops[ixout]->turnOrientation();
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      if ((int)ki == ixout)
	continue;
  
      if (CCW[ki])
	bd_loops[ki]->turnOrientation();
    }

  if (ixout != 0)
    std::swap(bd_loops[ixout], bd_loops[0]);  // Largest loop in first
}

//==============================================================================
void
CutCellQuad::identifyBdCvs(vector<shared_ptr<CurveLoop> >& bd_loops,
			   shared_ptr<CurveBoundedDomain> cell_domain,
			   vector<vector<shared_ptr<ParamCurve> > >& bd_seg)
//==============================================================================
{
  // First intersect fetch the part of the curve model boundares lying inside the cell
  vector<shared_ptr<ParamCurve> > bd_cvs;
  for (size_t ki=0; ki<bd_loops.size(); ++ki)
    {
      int nmb = bd_loops[ki]->size();
      for (int ka=0; ka<nmb; ++ka)
	{
	  shared_ptr<SplineCurve> cv((*bd_loops[ki])[ka]->geometryCurve());
	  vector<double> par_start_end;
	  cell_domain->findPcurveInsideSegments(*cv, tol_, par_start_end, true, true);
	  for (size_t kj=0; kj<par_start_end.size(); kj+=2)
	    bd_cvs.push_back(shared_ptr<ParamCurve>(cv->subCurve(par_start_end[kj],
								 par_start_end[kj+1])));
	  if (par_start_end.size() == 0)
	    {
	      // Check internal point
	      Point pos = cv->ParamCurve::point(0.5*(cv->startparam()+cv->endparam()));
	      if (cell_domain->isInDomain(Vector2D(pos[0],pos[1]), tol_))
		bd_cvs.push_back(cv);
	    }
	}
    }

  if (bd_cvs.size() == 0)
    return;   // No boundary curves inside cell
  if (bd_cvs.size() == 1)
    {
      bd_seg.push_back(bd_cvs);
      return;   // No sorting and grouping is required
    }
  
  // Group curves with G1 continuity
  // As the curves originate from a set of boundary loops, we can assume
  // consistent sequence and orientation
  vector<Point> endpts;
  vector<Point> der1(2), der2(2);
  bd_cvs[0]->point(der1, bd_cvs[0]->startparam(), 1);
  bd_cvs[0]->point(der2, bd_cvs[0]->endparam(), 1);
  vector<shared_ptr<ParamCurve> > curr_seg;
  curr_seg.push_back(bd_cvs[0]);
  bd_seg.push_back(curr_seg);
  size_t ix = 0;
  for (size_t kj=ix+1; kj<bd_cvs.size(); ++kj)
    {
      vector<Point> der3(2), der4(2);
      bd_cvs[kj]->point(der3, bd_cvs[kj]->startparam(), 1);
      bd_cvs[kj]->point(der4, bd_cvs[kj]->endparam(), 1);
      bool match = false;
      if (der2[0].dist(der3[0]) < tol_)
	{
	  // Check tangent continuity
	  double ang = der2[1].angle(der3[1]);
	  if (ang < angtol_)
	    {
	      bd_seg[ix].push_back(bd_cvs[kj]);
	      der2 = der4;
	      match = true;
	    }
	}
      if (!match)
	{
	  vector<shared_ptr<ParamCurve> > curr_seg2;
	  curr_seg2.push_back(bd_cvs[kj]);
	  bd_seg.push_back(curr_seg2);
	  der1 = der3;
	  der2 = der4;
	  ix++;
	}
    }

  if (bd_seg.size() > 1)
    {
      // Check for match across seam
      size_t ll1 = bd_seg.size()-1;
      size_t ll2 = bd_seg[ll1].size()-1;
      bd_seg[ll1][ll2]->point(der1, bd_seg[ll1][ll2]->endparam(), 1);
      bd_seg[0][0]->point(der2, bd_seg[0][0]->startparam(), 1);
      if (der1[0].dist(der2[0]) < tol_)
	{
	  // Check tangent continuity
	  double ang = der1[1].angle(der2[1]);
	  if (ang < angtol_)
	    {
	      bd_seg[0].insert(bd_seg[0].begin(), bd_seg[ll1].begin(), bd_seg[ll1].end());
	      bd_seg.pop_back();
	    }
	}
    }
}

//==============================================================================
void
CutCellQuad::createCutCell(vector<shared_ptr<CurveLoop> >& bd_loops,
			   vector<shared_ptr<SplineCurve> >& cell_bd,
			   vector<vector<shared_ptr<CurveLoop> > >& trim_loops,
			   bool test_inside)
//==============================================================================
{
  // First intersect fetch the part of the cell boundares lying inside the curve model
  shared_ptr<CurveBoundedDomain> domain(new CurveBoundedDomain(bd_loops));
  vector<shared_ptr<ParamCurve> > split_cvs;
  for (size_t ki=0; ki<cell_bd.size(); ++ki)
    {
      vector<double> par_start_end;
      domain->findPcurveInsideSegments(*cell_bd[ki], tol_, par_start_end, false, true);
      for (size_t kj=0; kj<par_start_end.size(); kj+=2)
	split_cvs.push_back(shared_ptr<ParamCurve>(cell_bd[ki]->subCurve(par_start_end[kj],
									 par_start_end[kj+1])));
    }

  // Sort curves according to decreasing length to avoid starting the boundary loop
  // recognition from a short curve
  for (size_t ki=0; ki<split_cvs.size(); ++ki)
    for (size_t kj=ki+1; kj<split_cvs.size(); ++kj)
      {
	double len1 = split_cvs[ki]->estimatedCurveLength(2);
	double len2 = split_cvs[kj]->estimatedCurveLength(2);
	if (len2 > len1)
	  std::swap(split_cvs[ki], split_cvs[kj]);
      }
#ifdef DEBUG  
  std::ofstream ofcv("split_cvs.g2");
  for (size_t ki=0; ki<split_cvs.size(); ++ki)
    {
      split_cvs[ki]->writeStandardHeader(ofcv);
      split_cvs[ki]->write(ofcv);
    }
#endif
  
  if (split_cvs.size() == 0)
    {
      // Touch situation
      return;
    }
      
  // Prepare for boundary loop construction of cut cell(s)
  size_t nmb_cvs = split_cvs.size();
  for (size_t ki=0; ki<nmb_cvs; ++ki)
    {
      shared_ptr<ParamCurve> tmp_crv = 
	shared_ptr<ParamCurve>(split_cvs[ki]->clone());
      tmp_crv->reverseParameterDirection();
      split_cvs.push_back(tmp_crv);
    }

  // Define tolerances
  double min_loop_tol = 1.1*loops_[0]->getMaxCurveDist();
  min_loop_tol = std::max(min_loop_tol, tol_);
  double scaling = 0.5;
  Point par_eps(scaling*min_loop_tol, scaling*min_loop_tol);
  
  CurveBoundedDomain cvdom(bd_loops);
  RectDomain dom = cvdom.containingDomain();

  // Compute trimmed loops
  vector<vector<shared_ptr<ParamCurve> > > trim_cvs =
    BoundedUtils::getBoundaryLoops(bd_loops, split_cvs, tol_, min_loop_tol,
				   min_loop_tol, par_eps);
#ifdef DEBUG
  std::ofstream of("cut_bd.g2");
  for (size_t ki=0; ki<trim_cvs.size(); ++ki)
    for (size_t kj=0; kj<trim_cvs[ki].size(); ++kj)
      {
	trim_cvs[ki][kj]->writeStandardHeader(of);
	trim_cvs[ki][kj]->write(of);
      }
#endif

  vector<shared_ptr<ParamCurve> > cell_bd2(cell_bd.begin(), cell_bd.end());
  sortLoops(trim_cvs, cell_bd2, trim_loops, test_inside);

#ifdef DEBUG
  std::cout << "Number of cell pieces: " << trim_loops.size() << std::endl;

  
  std::ofstream of2("cut_bd2.g2");
  for (size_t ki=0; ki<trim_loops.size(); ++ki)
    for (size_t kj=0; kj<trim_loops[ki].size(); ++kj)
      {
	int nmb = trim_loops[ki][kj]->size();
	for (int ka=0; ka<nmb; ++ka)
	  {
	    (*trim_loops[ki][kj])[ka]->writeStandardHeader(of2);
	    (*trim_loops[ki][kj])[ka]->write(of2);
	  }
      }

  int stop_break = 1;
#endif
}


//==============================================================================
void
CutCellQuad::sortLoops(vector<vector<shared_ptr<ParamCurve> > >& loop_cvs,
		       vector<shared_ptr<ParamCurve> >& cell_cvs,
		       vector<vector<shared_ptr<CurveLoop> > >& trim_loops,
		       bool test_inside)
//==============================================================================
{
  // Collect the loops lying inside the cell
  double int_tol = tol_; //std::max(1.0e-6, 0.1*tol_);
  shared_ptr<CurveLoop> cell_loop;
  shared_ptr<CurveBoundedDomain> cell_domain;
  if (test_inside)
    {
      cell_loop = shared_ptr<CurveLoop>(new CurveLoop(cell_cvs, tol_));
      cell_domain = shared_ptr<CurveBoundedDomain>(new CurveBoundedDomain(cell_loop));
    }
  vector<shared_ptr<CurveLoop> > in_loops_ccw;
  vector<shared_ptr<CurveLoop> > in_loops_cw;
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      // Check for degeneracy
      bool degen = BoundedUtils::loopIsDegenerate(loop_cvs[ki], tol_);
      if (degen)
	continue;
      
      shared_ptr<CurveLoop> curr_loop(new CurveLoop(loop_cvs[ki], tol_));
      bool inside = true;
      if (test_inside)
	{
	  shared_ptr<CurveBoundedDomain> trim_domain(new CurveBoundedDomain(curr_loop));
	  double upar, vpar;
	  trim_domain->getInternalPoint(upar, vpar);
	  inside = cell_domain->isInDomain(Vector2D(upar,vpar), tol_);
	  if (inside)
	    {
	      // Test also a boundary points
	      for (size_t kj=0; kj<loop_cvs[ki].size(); ++kj)
		{
		  double t1 = 0.5*(loop_cvs[ki][kj]->startparam() +
				   loop_cvs[ki][kj]->endparam());
		  Point bd_pt = loop_cvs[ki][kj]->point(t1);
		  bool bd_in = cell_domain->isInDomain(Vector2D(bd_pt[0], bd_pt[1]), tol_);
		  if (!bd_in)
		    {
		      inside = false;
		      break;
		    }
		}
	    }
	}
      if (inside)
	{
	  bool ccw = LoopUtils::loopIsCCW(loop_cvs[ki], tol_, int_tol);
	  if (ccw)
	    in_loops_ccw.push_back(curr_loop);
	  else
	    in_loops_cw.push_back(curr_loop);
	}
    }

  // Combine loops
  for (size_t ki=0; ki<in_loops_ccw.size(); ++ki)
    {
      vector<shared_ptr<CurveLoop> > outer_loop;
      outer_loop.push_back(in_loops_ccw[ki]);
      trim_loops.push_back(outer_loop);
    }

  // Distribute inner loops
  for (size_t kj=0; kj<in_loops_cw.size(); ++kj)
    {
      shared_ptr<CurveBoundedDomain> innerdom(new CurveBoundedDomain(in_loops_cw[kj]));
      double upar, vpar;
      innerdom->getInternalPoint(upar, vpar);
      size_t ki;
      for (ki=0; ki<trim_loops.size(); ++ki)
	{
	  shared_ptr<CurveBoundedDomain> outdom(new CurveBoundedDomain(trim_loops[ki][0]));
	  bool inside = outdom->isInDomain(Vector2D(upar, vpar), tol_);
	  if (inside)
	    {
	      trim_loops[ki].push_back(in_loops_cw[kj]);
	      break;
	    }
	}
    }
  int stop_break = 1;
}

//==============================================================================
void
CutCellQuad::quadraturePoints(vector<shared_ptr<CurveLoop> >& cell_loops,
			      vector<double>& quadraturepoints,
			      vector<double>& pointsweights,
			      vector<vector<shared_ptr<ParamCurve> > >& unresolved_cells,
			      vector<double>& curvequads,
			      vector<double>& curvenorms,
			      vector<double>& crvptweights,
			      vector<vector<shared_ptr<ParamCurve> > >& short_curves)
//==============================================================================
{
  #ifdef DEBUG
  std::ofstream of2("cell_curr.g2");
  for (size_t ki=0; ki<cell_loops.size(); ++ki)
    {
      int nmb = cell_loops[ki]->size();
      for (int ka=0; ka<nmb; ++ka)
	{
	  (*cell_loops[ki])[ka]->writeStandardHeader(of2);
	  (*cell_loops[ki])[ka]->write(of2);
	}
    }
#endif
  
  // Cell domain
  shared_ptr<CurveBoundedDomain> cvdom(new CurveBoundedDomain(cell_loops));
  RectDomain domain = cvdom->containingDomain();
  if (domain.umax() - domain.umin() < min_cell_size_ ||
      domain.vmax() - domain.vmin() < min_cell_size_)
    {
      vector<shared_ptr<ParamCurve> > unresolved_cvs;
      for (size_t ki=0; ki<cell_loops.size(); ++ki)
	{
	  int nmb = cell_loops[ki]->size();
	  for (int ka=0; ka<nmb; ++ka)
	    unresolved_cvs.push_back((*cell_loops[ki])[ka]);
	}
      unresolved_cells.push_back(unresolved_cvs);
      return;   // Cell not handled
    }

  // First check if a height direction can be found directly
  int basedir = (cell_loops.size() > 1) ? -1 : heightDirection(cell_loops[0], domain);;

  // Collect split points and select parameter direction
  vector<double> splitpar;
  //int basedir = -1;
  if (basedir <  0)
    splitPars(cell_loops, cvdom, domain, splitpar, basedir);

  if (splitpar.size() > 0)
    {
      // Split
      vector<shared_ptr<SplineCurve> > div_lines;
      for (size_t ki=0; ki<splitpar.size(); ++ki)
	{
	  Point pt1, pt2;
	  if (basedir == 1)
	    {
	      pt1 = Point(splitpar[ki], domain.vmin());
	      pt2 = Point(splitpar[ki], domain.vmax());
	    }
	  else
	    {
	      pt1 = Point(domain.umin(), splitpar[ki]);
	      pt2 = Point(domain.umax(), splitpar[ki]);
	    }
	  shared_ptr<SplineCurve> splitcv(new SplineCurve(pt1,pt2));
	  div_lines.push_back(splitcv);
	}

      vector<vector<shared_ptr<CurveLoop> > > trim_loops;
      createCutCell(cell_loops, div_lines, trim_loops, false);
      if (trim_loops.size() <= 1)
	{
	  vector<shared_ptr<ParamCurve> > unresolved_cvs;
	  for (size_t ki=0; ki<cell_loops.size(); ++ki)
	    {
	      int nmb = cell_loops[ki]->size();
	      for (int ka=0; ka<nmb; ++ka)
		unresolved_cvs.push_back((*cell_loops[ki])[ka]);
	    }
	  unresolved_cells.push_back(unresolved_cvs);
	  return;   // Cell not handled
	}
      for (size_t ki=0; ki<trim_loops.size(); ++ki)
	{
	  quadraturePoints(trim_loops[ki], quadraturepoints, pointsweights, unresolved_cells,
			   curvequads, curvenorms, crvptweights, short_curves);
	  int stop_break = 1;
	}
    }
  else
    {
      // // Find height direction, 1 = runs in u, 2 = runs in v
      // // There should be only the outer loop left now
      // if (basedir < 0)
      // 	basedir = heightDirection(cell_loops[0], domain);
 
      // Compute quadrature points
      vector<double> quadpts;
      vector<double> ptsweights;
      try {
	computeQuadraturePoints(cvdom, domain, basedir, quadpts, ptsweights);
      }
      catch (...)
	{
	  vector<shared_ptr<ParamCurve> > unresolved_cvs;
	  for (size_t ki=0; ki<cell_loops.size(); ++ki)
	    {
	      int nmb = cell_loops[ki]->size();
	      for (int ka=0; ka<nmb; ++ka)
		unresolved_cvs.push_back((*cell_loops[ki])[ka]);
	    }
	  unresolved_cells.push_back(unresolved_cvs);
	  return;   // Cell not handled
 	}
      quadraturepoints.insert(quadraturepoints.end(), quadpts.begin(), quadpts.end());
      pointsweights.insert(pointsweights.end(), ptsweights.begin(), ptsweights.end());

      // Quadrature points corresponding to boundary curves
      vector<vector<shared_ptr<ParamCurve> > > bd_segs;
      identifyBdCvs(loops_, cvdom, bd_segs);
      if (bd_segs.size() > 0)
	{
	  for (size_t ki=0; ki<bd_segs.size(); ++ki)
	    {
	      // Compute curve length
	      double len = 0.0;
	      for (size_t kj=0; kj<bd_segs[ki].size(); ++kj)
		len += bd_segs[ki][kj]->length(tol_);
	      if (len < min_cell_size_)
		short_curves.push_back(bd_segs[ki]);
	      else
		{
		  vector<double> crvquads;
		  vector<double> crvnorms;
		  vector<double> crvptwgts;
		  computeQuadraturePoints(bd_segs[ki], crvquads, crvnorms, crvptwgts);
		  curvequads.insert(curvequads.end(), crvquads.begin(), crvquads.end());
		  curvenorms.insert(curvenorms.end(), crvnorms.begin(), crvnorms.end());
		  crvptweights.insert(crvptweights.end(), crvptwgts.begin(), crvptwgts.end());
		}
	    }
	}
    }
}

//==============================================================================
void
CutCellQuad::computeQuadraturePoints(vector<shared_ptr<ParamCurve> >& bd_seg,
				     vector<double>& quadraturepoints,
				     vector<double>& quadraturenorms,
				     vector<double>& weights)
//==============================================================================
{
  // Reparameterize curves with respect to curve length and let parameterization
  // be sequential
  double tstart = bd_seg[0]->startparam();
  for (size_t ki=0; ki<bd_seg.size(); ++ki)
    {
      double cvlen = bd_seg[ki]->length(tol_);
      bd_seg[ki]->setParameterInterval(tstart, tstart+cvlen);
      tstart += cvlen;
    }

  // Make initial quadrature parameters
  double tend = tstart;
  tstart = bd_seg[0]->startparam();
  double del = tend - tstart;
  vector<double> quadpar(quadpar_.size());
  vector<double> wgt(quadpar_.size());
  for (size_t ki=0; ki<quadpar_.size(); ++ki)
    {
      quadpar[ki] = tstart + quadpar_[ki]*del;
      wgt[ki] = weights_[ki]*del;
    }

   // Compute quadrature points
  size_t ix = 0;
  vector<double> currquadpts;
  for (size_t ki=0; ki<quadpar.size(); ++ki)
    {
      while (bd_seg[ix]->endparam() < quadpar[ki])
	++ix;
      vector<Point> quadder(2);
      bd_seg[ix]->point(quadder, quadpar[ki], 1);
      quadraturepoints.push_back(quadder[0][0]);
      quadraturepoints.push_back(quadder[0][1]);
      Point quadnorm(quadder[1][1], -quadder[1][0]);
      (void)quadnorm.normalize_checked();
      quadraturenorms.push_back(quadnorm[0]);
      quadraturenorms.push_back(quadnorm[1]);
      weights.push_back(wgt[ki]*quadder[1].length());
     }
 }

//==============================================================================
int CutCellQuad::heightDirection(shared_ptr<CurveLoop> cell_loop,
				 const RectDomain& domain)
//==============================================================================
{
  int basedir = -1;  // Default is height direction not found

  // Distribute curves according to coincidence with the bounding box of the
  // cut cell domain
  vector<vector<shared_ptr<ParamCurve> > > sorted_curves;
  sorted_curves.resize(5);  // left, right, lower, upper, not along domain boundary
  int nmbcvs = cell_loop->size();
  double umaxlen = 0.0, vmaxlen = 0.0;
  for (int ka=0; ka<nmbcvs; ++ka)
    {
      shared_ptr<ParamCurve> cv = (*cell_loop)[ka];
      Point dir;
      if (cv->isLinear(dir, angtol_))
  	{
  	  Point pt1 = cv->point(cv->startparam());
  	  Point pt2 = cv->point(cv->endparam());
  	  int bd_ix = domain.whichBoundary(Vector2D(pt1[0],pt1[1]),
  					   Vector2D(pt2[0],pt2[1]), tol_);
  	  if (bd_ix == 0 || bd_ix == 1)
  	    {
	      vmaxlen = std::max(vmaxlen, fabs(pt1[1]-pt2[1]));
	      sorted_curves[bd_ix].push_back(cv);
  	    }
  	  else if (bd_ix == 2 || bd_ix == 3)
  	    {
	      umaxlen = std::max(umaxlen, fabs(pt1[0]-pt2[0]));
	      sorted_curves[bd_ix].push_back(cv);
  	    }
	  else
	    {
	      sorted_curves[4].push_back(cv);
	    }
  	}
      else
	{
	  sorted_curves[4].push_back(cv);
	}
    }

  if (sorted_curves[0].size() > 0 && sorted_curves[1].size() > 0 && sorted_curves[2].size() > 0 &&
      sorted_curves[3].size() > 0 && sorted_curves[4].size() == 0)
    {
      // Height direction defined from size of baseline
      return (umaxlen >= vmaxlen) ? 0 : 1;
    }

  if (sorted_curves[0].size() > 0 && sorted_curves[1].size() > 0)
    {
      // Potential for height direction along vertical edge
      // Collect remaining curves and sort into continuous sequences
      vector<shared_ptr<ParamCurve> > cvs;
      for (int ka=2; ka<5; ++ka)
	if (sorted_curves[ka].size() > 0)
	  cvs.insert(cvs.end(), sorted_curves[ka].begin(), sorted_curves[ka].end());
      if (checkSideCvs(cvs))
	basedir = 1;
      int stop_break1 = 1;
    }

  if (sorted_curves[2].size() > 0 && sorted_curves[3].size() > 0)
    {
      // Potential for height direction along horizontal edge
      vector<shared_ptr<ParamCurve> > cvs;
      for (int ka=0; ka<3; ++ka)
	{
	  if (ka == 2)
	    ka = 4;
	  if (sorted_curves[ka].size() > 0)
	    cvs.insert(cvs.end(), sorted_curves[ka].begin(), sorted_curves[ka].end());
	}
      if (checkSideCvs(cvs))
	basedir = 2;
      int stop_break2 = 1;
    }

  return basedir;
}

//==============================================================================
bool CutCellQuad::checkSideCvs(vector<shared_ptr<ParamCurve> >& cvs)
//==============================================================================
{
  // Sort curves. The curves are oriented correctly, but the sequence is not organized
  vector<vector<shared_ptr<ParamCurve> > > cv_seq;
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      Point pos1 = cvs[ki]->point(cvs[ki]->startparam());
      Point pos2 = cvs[ki]->point(cvs[ki]->endparam());
      size_t kj;
      for (kj=0; kj<cv_seq.size(); ++kj)
	{
	  Point p1 = cv_seq[kj][0]->point(cv_seq[kj][0]->startparam());
	  Point p2 = cv_seq[kj][cv_seq[kj].size()-1]->point(cv_seq[kj][cv_seq[kj].size()-1]->endparam());
	  double d1 = pos2.dist(p1);
	  double d2 = pos1.dist(p2);
	  if (d1 <= d2 && d1 < tol_)
	    {
	      cv_seq[kj].insert(cv_seq[kj].begin(), cvs[ki]);
	      break;
	    }
	  else if (d2 < tol_)
	    {
	      cv_seq[kj].push_back(cvs[ki]);
	      break;
	    }
	}
      if (kj == cv_seq.size())
	{
	  vector<shared_ptr<ParamCurve> > tmp;
	  tmp.push_back(cvs[ki]);
	  cv_seq.push_back(tmp);
	}
    }

  if (cv_seq.size() != 2)
    return false;  // Not a clear configuration

  // For each group of curves, check for corners 
  for (size_t ki=0; ki<cv_seq.size(); ++ki)
    {
      for (size_t kj=1; kj<cv_seq[ki].size(); ++kj)
	{
	  vector<Point> der1(2), der2(2);
	  cv_seq[ki][kj-1]->point(der1, cv_seq[ki][kj-1]->endparam(), 1);
	  cv_seq[ki][kj]->point(der2, cv_seq[ki][kj]->startparam(), 1);
	  double ang = der1[1].angle(der2[1]);
	  if (ang > angtol_)
	    return false;
	}
    }

  // Check for high curvature
  double anglim = M_PI/6.0;
   for (size_t ki=0; ki<cv_seq.size(); ++ki)
    {
      if (cv_seq[ki].size() == 0)
	continue;
      DirectionCone cone = cv_seq[ki][0]->directionCone();
      for (size_t kj=1; kj<cv_seq[ki].size(); ++kj)
	{
	  DirectionCone cone2 = cv_seq[ki][kj]->directionCone();
	  cone.addUnionWith(cone2);   // NB! Bug! Fix!
	}
      if (cone.angle() > anglim)
	return false;
    }

   
  return true;
}

#if 0
//==============================================================================
int CutCellQuad::heightDirection(shared_ptr<CurveLoop> cell_loop,
				 const RectDomain& domain)
//==============================================================================
{
  int basedir = -1;
  int nmbcvs = cell_loop->size();
  double umaxlen = 0.0, vmaxlen = 0.0;
  BoundingBox cvbox;
  for (int ka=0; ka<nmbcvs; ++ka)
    {
      shared_ptr<ParamCurve> cv = (*cell_loop)[ka];
      Point dir;
      if (cv->isLinear(dir, angtol_))
  	{
  	  Point pt1 = cv->point(cv->startparam());
  	  Point pt2 = cv->point(cv->endparam());
  	  int bd_ix = domain.whichBoundary(Vector2D(pt1[0],pt1[1]),
  					   Vector2D(pt2[0],pt2[1]), tol_);
  	  if (bd_ix == 0 || bd_ix == 1)
  	    {
	      vmaxlen = std::max(vmaxlen, fabs(pt1[1]-pt2[1]));
  	      // if (fabs(std::min(pt1[1],pt2[1])-domain.vmin()) < tol_ &&
  	      // 	  fabs(std::max(pt1[1],pt2[1])-domain.vmax()) < tol_)
  	      // 	{
  	      // 	  basedir = 2;
  	      // 	  break;
  	      // 	}
  	    }
  	  else if (bd_ix == 2 || bd_ix == 3)
  	    {
	      umaxlen = std::max(umaxlen, fabs(pt1[0]-pt2[0]));
  	      // if (fabs(std::min(pt1[0],pt2[0])-domain.umin()) < tol_ &&
  	      // 	  fabs(std::max(pt1[0],pt2[0])-domain.umax()) < tol_)
  	      // 	{
  	      // 	  basedir = 1;
  	      // 	  break;
  	      // 	}
  	    }
	  else
	    {
	      BoundingBox bb = cv->boundingBox();
	      if (cvbox.dimension() == 0)
		cvbox = bb;
	      else
		cvbox.addUnionWith(bb);
	    }
  	}
      else
	{
	  BoundingBox bb = cv->boundingBox();
	  if (cvbox.dimension() == 0)
	    cvbox = bb;
	  else
	    cvbox.addUnionWith(bb);
	}
    }
  
  if (cvbox.dimension() == 2)
    {
      double uboxlen = cvbox.high()[0] - cvbox.low()[0];
      double vboxlen = cvbox.high()[1] - cvbox.low()[1];
      double boxfraclim = 0.9;
      if (std::min(uboxlen,vboxlen)/std::max(uboxlen,vboxlen) < boxfraclim)
	basedir = (uboxlen < vboxlen) ? 2 : 1;
      else
	basedir = (umaxlen >= vmaxlen) ? 1 : 2;
    }
  else
    basedir = (umaxlen >= vmaxlen) ? 1 : 2;
  if (basedir < 0)
    {
      // For the time being
      THROW("Unresolved cell");
    }
  return basedir;
}
#endif

//==============================================================================
void
CutCellQuad::computeQuadraturePoints(shared_ptr<CurveBoundedDomain> cvdom,
				     const RectDomain& domain,
				     int dir,
				     vector<double>& quadraturepoints,
				     vector<double>& pointsweights)
//==============================================================================
{
  // Adapt quadrature values
  int dir2 = 3 - dir;
  double min1 = (dir == 1) ? domain.umin() : domain.vmin();
  double max1 = (dir == 1) ? domain.umax() : domain.vmax();
  double min2 = (dir == 1) ? domain.vmin() : domain.umin();
  double max2 = (dir == 1) ? domain.vmax() : domain.umax();
  double del1 = max1 - min1;
  vector<double> quadpar1(quadpar_.size());
  vector<double> wgt1(quadpar_.size());
  for (size_t ki=0; ki<quadpar_.size(); ++ki)
    {
      quadpar1[ki] = min1 + quadpar_[ki]*del1;
      wgt1[ki] = weights_[ki]*del1;
    }

  // Fetch constant parameter curves
  quadraturepoints.resize(2*quadpar_.size()*quadpar_.size());
  pointsweights.resize(quadpar_.size()*quadpar_.size());
  for (size_t ki=0; ki<quadpar1.size(); ++ki)
    {
      Point pt1, pt2;
      if (dir == 1)
	{
	  pt1 = Point(quadpar1[ki], min2);
	  pt2 = Point(quadpar1[ki], max2);
	}
      else
	{
	  pt1 = Point(min2, quadpar1[ki]);
	  pt2 = Point(max2, quadpar1[ki]);
	}
      shared_ptr<SplineCurve> crv(new SplineCurve(pt1, pt2));
      vector<double> inside;
      cvdom->findPcurveInsideSegments(*crv, tol_, inside, true, true);
      if (inside.size() != 2)
	{
	  THROW("Height function not properly defined");
	}

      // Compute quadrature
      Point pt3 = crv->ParamCurve::point(inside[0]);
      Point pt4 = crv->ParamCurve::point(inside[1]);
      int ix = (dir == 1) ? 1 : 0;
      double del2 = pt4[ix] - pt3[ix];
      for (size_t kj=0; kj<quadpar_.size(); ++kj)
	{
	  double quad2 = pt3[ix] + quadpar_[kj]*del2;
	  double wgt2 = weights_[kj]*del2;
	  if (dir == 1)
	    {
	      quadraturepoints[2*(kj*quadpar_.size()+ki)] = quadpar1[ki];
	      quadraturepoints[2*(kj*quadpar_.size()+ki)+1] = quad2;
	    }
	  else
	    {
	      quadraturepoints[2*(ki*quadpar_.size()+kj)] = quad2;
	      quadraturepoints[2*(ki*quadpar_.size()+kj)+1] = quadpar1[ki];
	    }
	  pointsweights[ki*quadpar_.size()+kj] = wgt1[ki]*wgt2;
	}
    }
}

//==============================================================================
void CutCellQuad::splitPars(vector<shared_ptr<CurveLoop> >& cell_loops,
			    shared_ptr<CurveBoundedDomain> cvdom,
			    const RectDomain& domain, vector<double>& splitpar,
			    int& dir)
//==============================================================================
{
  // Fetch loop corners including derivatives
  vector<Point> corner;
  int nmb_outer_corner;
  for (size_t ki=0; ki<cell_loops.size(); ++ki)
    {
      int nmb_corner = 0;
      int nmb = cell_loops[ki]->size();
      int kj, kk;
      for (kj=0, kk=1; kj<nmb; ++kj, kk=(kk+1)%nmb)
	{
	  vector<Point> der1(2), der2(2);
	  (*cell_loops[ki])[kj]->point(der1, (*cell_loops[ki])[kj]->endparam(), 1);
	  (*cell_loops[ki])[kk]->point(der2, (*cell_loops[ki])[kk]->startparam(), 1);
	  bool incorner = domain.isOnCorner(Vector2D(der1[0][0],der1[0][1]), tol_);
	  if (der1[1].angle(der2[1]) > angtol_)
	    {
	      if (!incorner)
		{
		  corner.insert(corner.end(), der1.begin(), der1.end());
		  corner.push_back(der2[1]);
		}
	      nmb_corner++;
	    }
	}
      if (ki == 0)
	nmb_outer_corner = nmb_corner;
    }
  
#ifdef DEBUG
  if (corner.size() > 0)
    {
      std::ofstream cp("corner.g2");
      cp << "400 1 0 4 255 0 0 255" << std::endl;
      cp << corner.size()/3 << std::endl;
      for (size_t ki=0; ki<corner.size(); ki+=3)
	cp << corner[ki] << "  0.0" << std::endl;
    }
#endif
  
  // Search for non-monotone parameter direction of trimming curves
  // Fetch also the bounding box containing the curves that are not parallel with
  // the parameter axis
  vector<vector<Point> > turnpts1, turnpts2;
  vector<vector<Point> > angpts;
  turnpts1.resize(cell_loops.size());
  turnpts2.resize(cell_loops.size());
  angpts.resize(cell_loops.size());
  BoundingBox cvbox;
  for (size_t ki=0; ki<cell_loops.size(); ++ki)
    {
      int nmb = cell_loops[ki]->size();
      for (int ka=0; ka<nmb; ++ka)
	{
	  shared_ptr<ParamCurve> cv = (*cell_loops[ki])[ka];
	  Point dir;
	  int bd_ix = -1;
	  if (cv->isLinear(dir, tol_))
	    {
	      Point pt1 = cv->point(cv->startparam());
	      Point pt2 = cv->point(cv->endparam());
	      bd_ix = domain.whichBoundary(Vector2D(pt1[0],pt1[1]),
					   Vector2D(pt2[0],pt2[1]), tol_);
	    }
	  else
	    {
	      fetchTurningPoints(cv, turnpts1[ki], turnpts2[ki], angpts[ki]);
	    }
	  if (bd_ix < 0)
	    {
	      BoundingBox bb = cv->boundingBox();
	      if (cvbox.dimension() == 0)
		cvbox = bb;
	      else
		cvbox.addUnionWith(bb);
	    }
	}
    }

  size_t nangpts = angpts[0].size();
  for (size_t ki=1; ki<angpts.size(); ++ki)
    nangpts += angpts[ki].size();
  
#ifdef DEBUG
  if (nangpts > 0)
    {
      std::ofstream cp("angpts.g2");
      cp << "400 1 0 4 100 155 0 255" << std::endl;
      cp << nangpts/3 << std::endl;
      for (size_t kj=0; kj<angpts.size(); ++kj)
	for (size_t ki=0; ki<angpts[kj].size(); ki+=3)
	  cp << angpts[kj][ki] << "  0.0" << std::endl;
    }
#endif
  
  size_t nturn1 = turnpts1[0].size();
  size_t nturn2 = turnpts2[0].size();
  for (size_t ki=1; ki<turnpts1.size(); ++ki)
    {
      nturn1 += turnpts1[ki].size();
      nturn2 += turnpts2[ki].size();
    }

#ifdef DEBUG
  if (nturn1+nturn2 > 0)
    {
      std::ofstream cp("turnpts.g2");
      cp << "400 1 0 4 0 155 100 255" << std::endl;
      cp << nturn1+nturn2 << std::endl;
      for (size_t kj=0; kj<turnpts1.size(); ++kj)
	{
	  for (size_t ki=0; ki<turnpts1[kj].size(); ++ki)
	    cp << turnpts1[kj][ki] << "  0.0" << std::endl;
	  for (size_t ki=0; ki<turnpts2[kj].size(); ++ki)
	    cp << turnpts2[kj][ki] << "  0.0" << std::endl;
	}
    }
#endif
  
  for (size_t ki=1; ki<angpts.size(); ++ki)
    corner.insert(corner.end(), angpts[ki].begin(), angpts[ki].end());

  vector<double> candpar;
  // int preferdir = (turnpts1.size() > 0 && turnpts2.size() == 0) ? 2 :
  //   ((turnpts2.size() > 0 && turnpts1.size() == 0) ? 1 : 0);
  int preferdir = (nturn1>0 && nturn2 == 0) ? 2 :
    ((nturn2 > 0 && nturn1 == 0) ? 1 : 0);
  defineSplits1(corner, cvdom, nmb_outer_corner, domain, cvbox,
		preferdir, splitpar, dir, candpar);

  if (splitpar.size() == 0)
    {
      // Check for closed inner loops
      vector<double> cell_split1, cell_split2;
      int cell_dir = -1;
      for (size_t ki=1; ki<cell_loops.size(); ++ki)
	{
	  // Check if any candidate parameter opens the loop
	  shared_ptr<CurveBoundedDomain> hole_domain(new CurveBoundedDomain(cell_loops[ki]));
	  RectDomain hole_box = hole_domain->containingDomain();
	  size_t kj;
	  for (kj=0; kj<candpar.size(); ++kj)
	    {
	      Point pt1, pt2;
	      if (dir == 1)
		{
		  pt1 = Point(hole_box.umin(), candpar[kj]);
		  pt2 = Point(hole_box.umax(), candpar[kj]);
		}
	      else
		{
		  pt1 = Point(candpar[kj], hole_box.vmin());
		  pt2 = Point(candpar[kj], hole_box.vmax());
		}
	      shared_ptr<SplineCurve> cv(new SplineCurve(pt1, pt2));
	      if (hole_domain->doIntersect(*cv, tol_))
		{
		  cell_split1.push_back(candpar[kj]);
		  cell_dir = 3 - dir;
		  break;
		}
	    }
	  if (kj == candpar.size())
	    {
	      double upar, vpar;
	      hole_domain->getInternalPoint(upar, vpar);
	      if (cell_dir == 1 || dir == 2)
		{
		  cell_split1.push_back(upar);
		  cell_split2.push_back(vpar);
		}
	      else
		{
		  cell_split1.push_back(vpar);
		  cell_split2.push_back(upar);
		  if (cell_dir < 0)
		    cell_dir = 1;
		}
	    }
	}
      if (cell_split1.size() > cell_split2.size())
	{
	  splitpar = cell_split1;
	  dir = cell_dir;
	}
      else if (cell_split2.size() > 0)
	{
	  splitpar = cell_split2;
	  if (dir <= 0)
	    dir = cell_dir;
	}
    }

  // if (splitpar.size() == 0 && turnpts1.size() + turnpts2.size() > 0)
  //   {
  //     defineSplits2(turnpts1, turnpts2, cvdom, domain, splitpar, dir, candpar);
  //   }

  if (dir == 1 && splitpar.size() == 0 && nturn1+nturn2+nangpts > 0)
    {
      vector<double> tmppar;
      for (size_t ki=0; ki<angpts.size(); ++ki)
	{
	  for (size_t kj=0; kj<angpts[ki].size(); kj+=3)
	    tmppar.push_back(angpts[ki][kj][0]);
	  for (size_t kj=0; kj<turnpts1[ki].size(); ++kj)
	    tmppar.push_back(turnpts1[ki][kj][0]);
	}
      double mid = 0.5*(domain.umin()+domain.umax());
      if (tmppar.size() == 0)
	splitpar.push_back(mid);
      else
	{
	  std::sort(tmppar.begin(), tmppar.end());
	  size_t ix = tmppar.size()/2; 
	  if (tmppar.size()%2 == 0 && fabs(tmppar[ix-1]-mid) < fabs(tmppar[ix]-mid))
	    --ix;
	  splitpar.push_back(tmppar[ix]);
	}
    }
  else if (dir == 2 && splitpar.size() == 0 && nturn1+nturn2+nangpts > 0)
    {
      vector<double> tmppar;
      for (size_t ki=0; ki<angpts.size(); ++ki)
	{
	  for (size_t kj=0; kj<angpts[ki].size(); kj+=3)
	    tmppar.push_back(angpts[ki][kj][1]);
	  for (size_t kj=0; kj<turnpts2[ki].size(); ++kj)
	    tmppar.push_back(turnpts2[ki][kj][1]);
	}
      double mid = 0.5*(domain.vmin()+domain.vmax());
      if (tmppar.size() == 0)
	splitpar.push_back(mid);
      else
	{
	  std::sort(tmppar.begin(), tmppar.end());
	  size_t ix = tmppar.size()/2; 
	  if (tmppar.size()%2 == 0 && fabs(tmppar[ix-1]-mid) < fabs(tmppar[ix]-mid))
	    --ix;
	  splitpar.push_back(tmppar[ix]);
	}
    }

  if (dir < 0 && splitpar.size() == 0)
    {
      // Subdivide in direction with highest curvature
      vector<shared_ptr<ParamCurve> > cvs;
        for (size_t ki=0; ki<cell_loops.size(); ++ki)
	  {
	    cvs.insert(cvs.end(), cell_loops[ki]->begin(), cell_loops[ki]->end());
	  }
      double curv1, curv2;
      computeDirCurvature(cvs, curv1, curv2);
      curv1 -= 2.0*(domain.umax() - domain.umin());
      curv2 -= 2.0*(domain.vmax() - domain.vmin());
      if (curv1 <= tol_ && curv2 <= tol_)
	{
	  // No subdivision
	  dir = (domain.umax() - domain.umin() >= domain.vmax() - domain.vmin()) ? 1 : 2;
	}
      else if (curv1 >= curv2)
	{
	  splitpar.push_back(0.5*(domain.vmin()+domain.vmax()));
	  dir = 2;
	}
      else
	{
	  splitpar.push_back(0.5*(domain.umin()+domain.umax()));
	  dir = 1;
	}
    }
  int stop_break = 1;
}

//==============================================================================
void CutCellQuad::defineSplits1(const vector<Point>& corner,
				shared_ptr<CurveBoundedDomain> cvdom,
				int nmb_outer_corner,
				const RectDomain& domain,
				const BoundingBox& cvbox,
				int preferdir,
				vector<double>& splitpar, int& dir,
				vector<double>& candpar)
//==============================================================================
{
  double pihalf = 0.5*M_PI;
  double umin = domain.umin();
  double umax = domain.umax();
  double vmin = domain.vmin();
  double vmax = domain.vmax();
  Point udir(1.0, 0.0);
  Point vdir(0.0, 1.0);
  vector<double> split_u1, split_u2, split_uin, split_uout;
  vector<double> split_v1, split_v2, split_vin, split_vout;

  // Sort corner information
  double eps = 1.0e-8;
  for (size_t ki=0; ki<corner.size(); ki+=3)
    {
      int on_v = (corner[ki][0]-umin < tol_) ? 1 : ((umax-corner[ki][0] < tol_) ? 2 : 0);
      int on_u = (corner[ki][1]-vmin < tol_) ? 1 : ((vmax-corner[ki][1] < tol_) ? 2 : 0);
      double uang1 = corner[ki+1].angle(udir);
      uang1 = std::min(uang1, M_PI-uang1);
      double uang2 = corner[ki+2].angle(udir);
      uang2 = std::min(uang2, M_PI-uang2);
      double vang1 = corner[ki+1].angle(vdir);
      vang1 = std::min(vang1, M_PI-vang1);
      double vang2 = corner[ki+2].angle(vdir);
      vang2 = std::min(vang2, M_PI-vang2);
      if (on_u == 1 && (uang1<angtol_ || uang2<angtol_))
	{
	  if (fabs(std::max(uang1,uang2)-pihalf) > angtol_)
	    split_u1.push_back(corner[ki][0]);
	}
      else if (on_u == 2 && (uang1<angtol_ || uang2<angtol_))
	{
	  if (fabs(std::max(uang1,uang2)-pihalf) > angtol_)
	    split_u2.push_back(corner[ki][0]);
	}
      else if (on_v == 1 && (vang1<angtol_ || vang2<angtol_))
	{
	  if (fabs(std::max(vang1,vang2)-pihalf) > angtol_)
	    split_v1.push_back(corner[ki][1]);
	}
      else if (on_v == 2 && (vang1<angtol_ || vang2<angtol_))
	{
	  if (fabs(std::max(vang1,vang2)-pihalf) > angtol_)
	    split_v2.push_back(corner[ki][1]);
	}
      else
	{
	  Point pt1 = corner[ki] + 2*tol_*udir;
	  Point pt2 = corner[ki] - 2*tol_*udir;
	  Point pt3 = corner[ki] + 2*tol_*vdir;
	  Point pt4 = corner[ki] - 2*tol_*vdir;
	  bool in1 = cvdom->isInDomain(Vector2D(pt1[0], pt1[1]), eps);
	  bool in2 = cvdom->isInDomain(Vector2D(pt2[0], pt2[1]), eps);
	  bool in3 = cvdom->isInDomain(Vector2D(pt3[0], pt3[1]), eps);
	  bool in4 = cvdom->isInDomain(Vector2D(pt4[0], pt4[1]), eps);
	  if (in1 || in2)
	    split_vin.push_back(corner[ki][1]);
	  else
	    split_vout.push_back(corner[ki][1]);
	  if (in3 || in4)
	    split_uin.push_back(corner[ki][0]);
	  else
	    split_uout.push_back(corner[ki][0]);
	}
    }


  size_t num_u = split_u1.size() + split_u2.size() + split_uin.size();
  size_t num_v = split_v1.size() + split_v2.size() + split_vin.size();
  int preferdir2 = preferdir;
  if (cvbox.dimension() == 2)
    {
      double uboxlen = cvbox.high()[0] - cvbox.low()[0];
      double vboxlen = cvbox.high()[1] - cvbox.low()[1];
      double boxfraclim = 0.9;
      if (preferdir == 0 && std::min(uboxlen,vboxlen)/std::max(uboxlen,vboxlen) < boxfraclim)
	preferdir2 = (uboxlen < vboxlen) ? 2 : 1;
    }
  if (num_u == 0 && num_v > 0)
    {
      bool samecv = false;
      if (split_v1.size() == 1 && split_v2.size() == 1)
	{
	  Vector2D pt1(umin, split_v1[0]);
	  Vector2D pt2(umax, split_v2[0]);
	  samecv = cvdom->onSmoothBdSeg(pt1, pt2, tol_, angtol_);
	}
      if (preferdir == 1)
	{
	  dir = 1;
	  candpar.insert(candpar.end(), split_v1.begin(), split_v1.end());
	  candpar.insert(candpar.end(), split_v2.begin(), split_v2.end());
	  std::sort(candpar.begin(), candpar.end());
	}
      else if (((split_v1.size() == 1 && split_v2.size() == 1 && samecv) ||
	   (split_v1.size()+split_v2.size() == 1 && nmb_outer_corner <= 4)) &&
	  split_vin.size() == 0 && split_vout.size() == 0)
	{
	  dir = 1;
	  candpar.insert(candpar.end(), split_v1.begin(), split_v1.end());
	  candpar.insert(candpar.end(), split_v2.begin(), split_v2.end());
	  std::sort(candpar.begin(), candpar.end());
	}
      else if (split_v1.size() > 0 && split_v2.size() > 0)
	{
	  dir = 2;
	  double lim = 0.2*(vmax-vmin);
	  double mindist = vmax - vmin;
	  for (size_t ki=0; ki<split_v1.size(); ++ki)
	    for (size_t kj=0; kj<split_v2.size(); ++kj)
	      {
		double dd = fabs(split_v1[ki]-split_v2[kj]);
		mindist = std::min(mindist, dd);
	      }
	  splitpar.insert(splitpar.end(), split_v1.begin(), split_v1.end());
	  splitpar.insert(splitpar.end(), split_v2.begin(), split_v2.end());
	  if (split_uin.size() > 0)
	    splitpar.insert(splitpar.end(), split_vin.begin(), split_vin.end());
	  std::sort(splitpar.begin(), splitpar.end());
	  if (splitpar.size() > 0)
	    checkSplits2(splitpar, vmin, vmax, lim);
	  if (splitpar.size() == 0 || mindist < lim)
	    {
	      splitpar.push_back(0.5*(umin+umax));
	      dir = 1;
	    }
	}
      else 
	{
	  dir = 2;
	  if (split_v1.size() > split_v2.size())
	    splitpar.insert(splitpar.end(), split_v1.begin(), split_v1.end());
	  else
	    splitpar.insert(splitpar.end(), split_v2.begin(), split_v2.end());
	  if (split_uin.size() > 0)
	    splitpar.insert(splitpar.end(), split_vin.begin(), split_vin.end());
	  std::sort(splitpar.begin(), splitpar.end());
	}
    }
  else if (num_v == 0 && num_u > 0)
    {
      bool samecv = false;
      if (split_u1.size() == 1 && split_u2.size() == 1)
	{
	  Vector2D pt1(split_u1[0], vmin);
	  Vector2D pt2(split_u2[0], vmax);
	  samecv = cvdom->isOnSameBdCrv(pt1, pt2, tol_);
	}
      if (preferdir == 2)
      	{
	  dir = 2;
	  candpar.insert(candpar.end(), split_u1.begin(), split_u1.end());
	  candpar.insert(candpar.end(), split_u2.begin(), split_u2.end());
	  std::sort(candpar.begin(), candpar.end());
	}
      else if (((split_u1.size() == 1 && split_u2.size() == 1 && samecv) ||
	   (split_u1.size()+split_u2.size() == 1 && nmb_outer_corner <= 4)) &&
	  split_uin.size() == 0 && split_vout.size() == 0)
	{
	  dir = 2;
	  candpar.insert(candpar.end(), split_u1.begin(), split_u1.end());
	  candpar.insert(candpar.end(), split_u2.begin(), split_u2.end());
	  std::sort(candpar.begin(), candpar.end());
	}
      else if (split_u1.size() > 0 && split_u2.size() > 0)
	{
	  dir = 1;
	  double lim = 0.2*(umax-umin);
	  double mindist = umax - umin;
	  for (size_t ki=0; ki<split_u1.size(); ++ki)
	    for (size_t kj=0; kj<split_u2.size(); ++kj)
	      {
		double dd = fabs(split_u1[ki]-split_u2[kj]);
		mindist = std::min(mindist, dd);
	      }
	  splitpar.insert(splitpar.end(), split_u1.begin(), split_u1.end());
	  splitpar.insert(splitpar.end(), split_u2.begin(), split_u2.end());
	  if (split_uin.size() > 0)
	    splitpar.insert(splitpar.end(), split_uin.begin(), split_uin.end());
	  std::sort(splitpar.begin(), splitpar.end());
	  
	  // Check for close split parameters
	  if (splitpar.size() > 1)
	    checkSplits2(splitpar, umin, umax, 0.2*(umax-umin));
	  if (splitpar.size() == 0 || mindist < lim)
	    {
	      splitpar.push_back(0.5*(vmin+vmax));
	      dir = 2;
	    }
	}
     else 
	{
	  dir = 1;
	  if (split_u1.size() > split_u2.size())
	    splitpar.insert(splitpar.end(), split_u1.begin(), split_u1.end());
	  else
	    splitpar.insert(splitpar.end(), split_u2.begin(), split_u2.end());
	  if (split_uin.size() > 0)
	    splitpar.insert(splitpar.end(), split_uin.begin(), split_uin.end());
	  std::sort(splitpar.begin(), splitpar.end());

	  // // Check for close split parameters
	  // bool OK = checkSplits(splitpar);
	  // if (!OK)
	  //   {
	  //     splitpar.clear();
	  //     splitpar.push_back(0.5*(vmin+vmax));
	  //     dir = 2;
	  //   }
	}
    }
  else if (num_u > 0 && num_v > 0)
    {
      if (split_uin.size() > 0 && !(split_vin.size() > 0 && preferdir2 == 2))
	{
	  dir = 1;
	  splitpar = split_uin;
	}
      else if (split_vin.size() > 0)
	{
	  dir = 2;
	  splitpar = split_vin;
	}
      else if (split_u1.size() == 0 && split_u2.size() > std::max(split_v1.size(), split_v2.size()))
	{
	  dir = 1;
	  splitpar = split_u2;
	}
      else if (split_u2.size() == 0 && split_u1.size() > std::max(split_v1.size(), split_v2.size()))
	{
	  dir = 1;
	  splitpar = split_u1;
	}
      else if (split_v1.size() == 0 && split_v2.size() > std::max(split_u1.size(), split_u2.size()))
	{
	  dir = 2;
	  splitpar = split_v2;
	}
      else if (split_v2.size() == 0 && split_v1.size() > std::max(split_u1.size(), split_u2.size()))
	{
	  dir = 2;
	  splitpar = split_v1;
	}
      else if (num_u == 1 && split_v1.size() > 0 && split_v2.size() > 0)
	{
	  dir = 1;
	  splitpar = (split_u1.size() > 0) ? split_u1 : split_u2;
	}
      else if (num_v == 1 && split_u1.size() > 0 && split_u2.size() > 0)
	{
	  dir = 2;
	  splitpar = (split_v1.size() > 0) ? split_v1 : split_v2;
	}
      else if (num_u == 1 && num_v > 1)
	{
	  vector<double> tmpsplit = split_v1;
	  tmpsplit.insert(tmpsplit.end(), split_v2.begin(), split_v2.end());
	  std::sort(tmpsplit.begin(), tmpsplit.end());
	  dir = 2;
	  if (split_u1.size() == 1)
	    splitpar.push_back(tmpsplit[0]);
	  else
	    splitpar.push_back(tmpsplit[tmpsplit.size()-1]);
	}
      else if (num_u == 0 && num_v == 1)
	{
	  vector<double> tmpsplit = split_u1;
	  tmpsplit.insert(tmpsplit.end(), split_u2.begin(), split_u2.end());
	  std::sort(tmpsplit.begin(), tmpsplit.end());
	  dir = 1;
	  if (split_v1.size() == 1)
	    splitpar.push_back(tmpsplit[0]);
	  else
	    splitpar.push_back(tmpsplit[tmpsplit.size()-1]);
	}
      else if (num_u == 1 && num_v == 1)
	{
	  if (preferdir == 1 || (preferdir == 0 && umax-umin > vmax-vmin))
	    {
	      dir = 1;
	      splitpar.insert(splitpar.begin(), split_u1.begin(), split_u1.end());
	      splitpar.insert(splitpar.begin(), split_u2.begin(), split_u2.end());
	    }
	    else
	    {
	      dir = 2;
	      splitpar.insert(splitpar.begin(), split_v1.begin(), split_v1.end());
	      splitpar.insert(splitpar.begin(), split_v2.begin(), split_v2.end());
	    }
	}
      else
	{
	  if (preferdir == 1 || (preferdir == 0 && umax-umin > vmax-vmin))
	    {
	      dir = 1;
	      vector<double> tmpsplit = split_u1;
	      tmpsplit.insert(tmpsplit.end(), split_u2.begin(), split_u2.end());
	      double mid = 0.5*(domain.umin()+domain.umax());
	      size_t ix = 0;
	      double len = fabs(tmpsplit[0]-mid);
	      for (size_t ki=1; ki<tmpsplit.size(); ++ki)
		{
		  double len2 = fabs(tmpsplit[ki]-mid);
		  if (len2 < len)
		    {
		      ix = ki;
		      len = len2;
		    }
		}
	      splitpar.push_back(tmpsplit[ix]);
	    }
	  else
	    {
	      dir = 2;
	      vector<double> tmpsplit = split_v1;
	      tmpsplit.insert(tmpsplit.end(), split_v2.begin(), split_v2.end());
	      double mid = 0.5*(domain.vmin()+domain.vmax());
	      size_t ix = 0;
	      double len = fabs(tmpsplit[0]-mid);
	      for (size_t ki=1; ki<tmpsplit.size(); ++ki)
		{
		  double len2 = fabs(tmpsplit[ki]-mid);
		  if (len2 < len)
		    {
		      ix = ki;
		      len = len2;
		    }
		}
	      splitpar.push_back(tmpsplit[ix]);
	    }
	}

      std::sort(splitpar.begin(), splitpar.end());
      if (dir == 1 && splitpar.size() > 1)
	checkSplits2(splitpar, umin, umax, 0.2*(umax-umin));
      else if (dir == 2 && splitpar.size() > 1)
	checkSplits2(splitpar, vmin, vmax, 0.2*(vmax-vmin));
      if (split_uin.size() == 0 && split_vin.size() == 0 && preferdir == 0)
	{
	  // Check split parameters
	  vector<double> tmpsplit(splitpar.begin(), splitpar.end());
	  double minp = (dir == 1)  ? domain.umin() : domain.vmin();
	  double maxp = (dir == 1)  ? domain.umax() : domain.vmax();
	  double lim = std::max(min_cell_size_, 0.1*(maxp - minp));
	  checkSplits2(tmpsplit, minp, maxp, lim);
	  if (tmpsplit.size() == 0)
	    {
	      // Try the other parameter direction
	      if (dir == 1)
		{
		  if (split_v1.size() > split_v2.size())
		    tmpsplit.insert(tmpsplit.begin(), split_v1.begin(), split_v1.end());
		  else
		    tmpsplit.insert(tmpsplit.begin(), split_v2.begin(), split_v2.end());
		}
	      else
		{
		  if (split_u1.size() > split_u2.size())
		    tmpsplit.insert(tmpsplit.begin(), split_u1.begin(), split_u1.end());
		  else
		    tmpsplit.insert(tmpsplit.begin(), split_u2.begin(), split_u2.end());
		}
	      checkSplits2(tmpsplit, minp, maxp, lim);
	      if (tmpsplit.size() > 0)
		{
		  splitpar = tmpsplit;
		  dir = 3 - dir;
		}
	    }
	}
 
      //  vector<double> tmpsplit1, tmpsplit2;
      // tmpsplit1.insert(tmpsplit1.end(), split_u1.begin(), split_u1.end());
      // tmpsplit1.insert(tmpsplit1.end(), split_u2.begin(), split_u2.end());
      // if (split_uin.size() > 0)
      // 	tmpsplit1.insert(tmpsplit1.end(), split_uin.begin(), split_uin.end());
      // std::sort(tmpsplit1.begin(), tmpsplit1.end());
      // tmpsplit2.insert(tmpsplit2.end(), split_v1.begin(), split_v1.end());
      // tmpsplit2.insert(tmpsplit2.end(), split_v2.begin(), split_v2.end());
      // if (split_vin.size() > 0)
      // 	tmpsplit2.insert(tmpsplit2.end(), split_vin.begin(), split_vin.end());
      // std::sort(tmpsplit2.begin(), tmpsplit2.end());

      // double u1 = umin;
      // double minu = umax - umin;
      // for (size_t ki=0; ki<tmpsplit1.size(); ++ki)
      // 	{
      // 	  minu = std::min(minu, tmpsplit1[ki]-u1);
      // 	  u1 = tmpsplit1[ki];
      // 	}
      // minu = std::min(minu, umax-u1);
			  
      // double v1 = vmin;
      // double minv = vmax - vmin;
      // for (size_t ki=0; ki<tmpsplit2.size(); ++ki)
      // 	{
      // 	  minv = std::min(minv, tmpsplit2[ki]-v1);
      // 	  v1 = tmpsplit2[ki];
      // 	}
      // minv = std::min(minv, vmax-v1);
      // double ufrac = minu/(umax - umin);
      // double vfrac = minv/(vmax - vmin);
      // if (fabs(ufrac-vfrac) < 0.1)
      // 	ufrac = (umax-umin > vmax-vmin) ? 2*vfrac : 0.5*vfrac;

      // double uboxlen = cvbox.high()[0] - cvbox.low()[0];
      // double vboxlen = cvbox.high()[1] - cvbox.low()[1];
      // double boxfraclim = 0.9;
      // if (preferdir == 0 && std::min(uboxlen,vboxlen)/std::max(uboxlen,vboxlen) < boxfraclim)
      // 	preferdir = (uboxlen < vboxlen) ? 2 : 1;
      
      // if (preferdir == 1 || split_uin.size() > 0 ||
      // 	  ((num_u > num_v || (num_u == num_v && ufrac > vfrac)) && split_vin.size() == 0))
      // 	{
      // 	  dir = 1;
      // 	  splitpar = tmpsplit1;

      // 	  // Check for close split parameters
      // 	  bool OK = checkSplits(splitpar);
      // 	  if (!OK)
      // 	    {
      // 	      // Try the other direction
      // 	      dir = 2;
      // 	      splitpar = tmpsplit2;

      // 	      // Check for close split parameters
      // 	      bool OK2 = checkSplits(splitpar);
      // 	      if (!OK2)
      // 		{
      // 		  std::cout << "Must select a good split parameter" << std::endl;
      // 		}
      // 	    }
      // 	}
      // else
      // 	{
      // 	  dir = 2;
      // 	  splitpar = tmpsplit2;

      // 	  // Check for close split parameters
      // 	  bool OK = checkSplits(splitpar);
      // 	  if (!OK)
      // 	    {
      // 	      // Try the other direction
      // 	      dir = 1;
      // 	      splitpar = tmpsplit1;

      // 	      // Check for close split parameters
      // 	      bool OK2 = checkSplits(splitpar);
      // 	      if (!OK2)
      // 		{
      // 		  std::cout << "Must select a good split parameter" << std::endl;
      // 		}
      // 	    }
      // 	}
    }
  else if (split_uout.size() > 0 || split_vout.size() > 0)
    {
      if (split_vout.size() > split_uout.size())
	{
	  dir = 2;
	  splitpar = split_uout;
	}
      else
	{
	  dir = 1;
	  splitpar = split_vout;
	}
    }

  if (splitpar.size() > 1)
    {
      std::sort(splitpar.begin(), splitpar.end());
      double minp = (dir == 1)  ? domain.umin() : domain.vmin();
      double maxp = (dir == 1)  ? domain.umax() : domain.vmax();
      double midp = 0.5*(minp + maxp);
      double lim = std::max(min_cell_size_, 0.1*(maxp - minp));
      if (splitpar[0] - minp < lim && maxp - splitpar[splitpar.size()-1] >= lim)
	splitpar.erase(splitpar.begin(), splitpar.begin()+1);
      else if (maxp - splitpar[splitpar.size()-1] < lim && splitpar[0] - minp >= lim)
	splitpar.pop_back();
      for (size_t ki=1; ki<splitpar.size(); )
	{
	  if (splitpar[ki] - splitpar[ki-1] < lim)
	    {
	      size_t ix = (fabs(midp - splitpar[ki-1]) < fabs(midp - splitpar[ki])) ? ki : ki-1;
	      splitpar.erase(splitpar.begin()+ix);
	    }
	  else
	    ++ki;
	}
    }
      
  
  // Check candidate parameters
  bool OKcand = checkSplits(candpar);
  if (!OKcand)
    candpar.clear();
}

//==============================================================================
void CutCellQuad::defineSplits2(const vector<Point>& turnpts1,
				const vector<Point>& turnpts2,
				shared_ptr<CurveBoundedDomain> cvdom,
				const RectDomain& domain,
				vector<double>& splitpar, int& dir,
				vector<double>& candpar)
//==============================================================================
{
  if (turnpts1.size() == 0)
    {
      dir = 1;
    }
  else if (turnpts2.size() == 0)
    {
      dir = 2;
    }
  else
    {
      if (dir == 1 && candpar.size() > 0)
	{
	  dir = 2;

	  // Find the candidate parameter closest to a turning point
	  size_t ix = 0;
	  double len = fabs(candpar[0] - turnpts1[0][0]);
	  for (size_t ki=0; ki<turnpts1.size(); ++ki)
	    for (size_t kj=0; kj<candpar.size(); ++kj)
	      {
		double len2 = fabs(candpar[kj] - turnpts1[ki][0]);
		if (len2 < len)
		  {
		    len = len2;
		    ix = kj;
		  }
	      }
	  splitpar.push_back(candpar[ix]);
	}
      else if (dir == 2 && candpar.size() > 0)
	{
	  dir = 1;

	  // Find the candidate parameter closest to a turning point
	  size_t ix = 0;
	  double len = fabs(candpar[0] - turnpts2[0][1]);
	  for (size_t ki=0; ki<turnpts2.size(); ++ki)
	    for (size_t kj=0; kj<candpar.size(); ++kj)
	      {
		double len2 = fabs(candpar[kj] - turnpts2[ki][1]);
		if (len2 < len)
		  {
		    len = len2;
		    ix = kj;
		  }
	      }
	  splitpar.push_back(candpar[ix]);
	}
      else
	{
	  double umid = 0.5*(domain.umin() + domain.umax());
	  double vmid = 0.5*(domain.vmin() + domain.vmax());
	  dir = 1;
	  size_t ix = 0;
	  double len = fabs(turnpts1[0][0] - umid);
	  for (size_t ki=1; ki<turnpts1.size(); ++ki)
	    {
	      double len2 = fabs(turnpts1[ki][0] - umid);
	      if (len2 < len)
		{
		  len = len2;
		  ix = ki;
		}
	    }
	  
	  for (size_t ki=0; ki<turnpts2.size(); ++ki)
	    {
	      double len2 = fabs(turnpts2[ki][1] - vmid);
	      if (len2 < len)
		{
		  dir = 2;
		  len = len2;
		  ix = ki;
		}
	    }

	  splitpar.push_back((dir == 1) ? turnpts1[ix][0] : turnpts2[ix][1]);
	}
    }
}
 
//==============================================================================
bool CutCellQuad::checkSplits(vector<double>& splitpar, double del)
//==============================================================================
{
  double eps = 1.0e-10;
  size_t kj;
  double del2 = (del > 0) ? del : min_cell_size_;
  for (kj=1; kj<splitpar.size();)
    {
      if (splitpar[kj] - splitpar[kj-1] <= eps)
	{
	  splitpar[kj-1] = 0.5*(splitpar[kj-1]+splitpar[kj]);
	  splitpar.erase(splitpar.begin()+kj);
	}
      else if (splitpar[kj] - splitpar[kj-1] < del2)
	break;
      else
	++kj;
    }
  return (kj < splitpar.size()) ? false : true;
}

//==============================================================================
void CutCellQuad::checkSplits2(vector<double>& splitpar, double minp, double maxp,
			       double lim)
//==============================================================================
{
  double midp = 0.5*(minp + maxp);
  double dd1 = splitpar[0] - minp;
  double dd2 = maxp - splitpar[splitpar.size()-1];
  if (dd1 < lim && (dd1 < dd2 || splitpar.size() > 2))
    splitpar.erase(splitpar.begin(), splitpar.begin()+1);
  if (dd2 < lim && splitpar.size() > 0)
    splitpar.pop_back();
  for (size_t ki=1; ki<splitpar.size(); )
    {
      if (splitpar[ki] - splitpar[ki-1] < lim)
	{
	  size_t ix = (fabs(midp - splitpar[ki-1]) < fabs(midp - splitpar[ki])) ? ki : ki-1;
	  splitpar.erase(splitpar.begin()+ix);
	}
      else
	++ki;
    }
}

//==============================================================================
void CutCellQuad::fetchTurningPoints(shared_ptr<ParamCurve> cv,
				     vector<Point>& turnpts1, vector<Point>& turnpts2,
				     vector<Point>& angpts)
//==============================================================================
{
  double eps = 1.0e-12;
  shared_ptr<SplineCurve> splcv(cv->geometryCurve());
  DirectionCone cone = splcv->directionCone();
  double anglim = M_PI/6.0;
  if (cone.angle() > anglim)
    {
      anglim = std::min(anglim, 0.5*cone.angle());
      vector<double>::iterator c1 = splcv->coefs_begin();
      vector<double>::iterator c2 = c1 + 2;
      int numcf = splcv->numCoefs();
      Point vec = Point(c2, c2+2) - Point(c1, c1+2);
      DirectionCone cone2;
      cone2.setFromArray(vec.begin(), vec.begin()+2, 2);
      c1 += 2;
      c2 += 2;
      for (int ki=2; ki<numcf; ++ki, c1+=2, c2+=2)
	{
	  Point vec = Point(c2, c2+2) - Point(c1, c1+2);
	  cone2.addUnionWith(vec);
	  if (cone2.angle() > anglim)
	    {
	      
	      double t1 = splcv->basis().grevilleParameter(ki-2);
	      double t2 = splcv->basis().grevilleParameter(ki);
	      vector<Point> der(2);
	      splcv->point(der, 0.5*(t1+t2), 1);
	      angpts.push_back(der[0]);
	      angpts.push_back(der[1]);
	      angpts.push_back(der[1]);
	      break;
	    }
	}
   }

  for (int pardir=0; pardir<2; ++pardir)
    {
      int prev_sgn = 0;
      int sgn;
      vector<double>::iterator c1 = splcv->coefs_begin();
      vector<double>::iterator c2 = c1 + 2;
      int numcf = splcv->numCoefs();
      for (int ki=1; ki<numcf; ++ki, c1+=2, c2+=2)
	{
	  double dd = *(c2+pardir) - *(c1+pardir);
	  sgn = (fabs(dd) < eps) ? 0 : ((dd < 0) ? -1 : 1);
	  if (sgn*prev_sgn < 0)
	    {
	      double t1 = splcv->basis().grevilleParameter(ki-1);
	      double t2 = splcv->basis().grevilleParameter(ki);
	      Point pos = splcv->ParamCurve::point(0.5*(t1+t2));
	      if (pardir == 0)
		turnpts1.push_back(pos);
	      else
		turnpts2.push_back(pos);
	    }
	  if (sgn != 0)
	    prev_sgn = sgn;
	}
    }
}

//==============================================================================
void CutCellQuad::computeDirCurvature(vector<shared_ptr<ParamCurve> > cvs,
				      double& curv1, double& curv2)
//==============================================================================
{
  curv1 = curv2 = 0.0;
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      shared_ptr<SplineCurve> splcv(cvs[ki]->geometryCurve());
      vector<double>::iterator c1 = splcv->coefs_begin();
      vector<double>::iterator c2 = c1 + 2;
      int numcf = splcv->numCoefs();
      for (int ki=1; ki<numcf; ++ki, c1+=2, c2+=2)
	{
	  double dd1 = *(c2) - *(c1);
	  double dd2 = *(c2+1) - *(c1+1);
	  curv1 += fabs(dd1);
	  curv2 += fabs(dd2);
	}
    }
}
