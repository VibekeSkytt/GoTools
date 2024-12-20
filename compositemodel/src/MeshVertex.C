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

#include "GoTools/compositemodel/MeshVertex.h"

using namespace Go;

//===========================================================================
MeshVertex::MeshVertex()
  :ftSamplePoint(), mean_edge_len_(0.0), area_(0.0), marked_(false)
//===========================================================================
{
  prev1_ = prev2_ = 0;
}

//===========================================================================
MeshVertex::MeshVertex(RevEngPoint* pt)
  :ftSamplePoint(pt->getPoint(), 0), mean_edge_len_(0.0), area_(0.0),
   marked_(false)
//===========================================================================
{
  setIndex(pt->getIndex());
  int num = pt->getNmbNeighbour();
  neighbour_idx_.resize(num);
  vector<ftSamplePoint*> next = pt->getNeighbours();
  for (int ka=0; ka<num; ++ka)
    neighbour_idx_[ka] = next[ka]->getIndex();

  pos_ = Point(xyz_[0], xyz_[1], xyz_[2]);
  Point norm1 = pt->getLocFuncNormal();
  Point norm2 = pt->getTriangNormal();
  if (norm1*norm2 < 0.0)
    norm2 *= -1;
  norm_ = 0.5*(norm1+norm2);
  norm_.normalize();
  dir_ = pt->minCurvatureVec();

  Point tmp_vec = norm_.cross(dir_);
  dir_ = tmp_vec.cross(norm_);
  dir_.normalize();
  tmp_vec = norm_.cross(dir_);
  tmp_vec.normalize();

  double eps = 1.0e-6;
  pos2_ = pos_ + eps*dir_ + eps*tmp_vec;
  dir2_ = Point(0.0, 0.0, 0.0);
  prev1_ = prev2_ = up_ = 0;
}

//===========================================================================
MeshVertex::MeshVertex(MeshVertex *v1, MeshVertex *v2)
  :ftSamplePoint(), mean_edge_len_(0.0), area_(0.0), marked_(false),
   prev1_(v1), prev2_(v2), up_(0)
//===========================================================================
{
  double eps = 1.0e-9;
  double area1 = v1->getArea();
  double area2 = v2->getArea();
  area_ = area1 + area2;
  xyz_ = (area_ > eps) ?
    (area1*v1->getPoint() + area2*v2->getPoint())/area_ :
    0.5*(v1->getPoint() + v2->getPoint());
  norm_ = area1*v1->getNormal() + area2*v2->getNormal();
  norm_.normalize_checked();
  Point dir1 = v1->getDirection();
  Point dir2 = v2->getDirection();
  if (dir1*dir2 < 0.0)
    dir2 *= -1;
  dir_ = area1*dir1 + area2*dir2;
  dir_.normalize_checked();

  pos_ = Point(xyz_[0], xyz_[1], xyz_[2]);
  Point tmp_vec = norm_.cross(dir_);
  dir_ = tmp_vec.cross(norm_);
  dir_.normalize();
  tmp_vec = norm_.cross(dir_);
  tmp_vec.normalize();

  pos2_ = pos_ + eps*dir_ + eps*tmp_vec;
  dir2_ = Point(0.0, 0.0, 0.0);

  // Transfer neighbourhood
  size_t ki, kj, kr, kh;
  vector<ftSamplePoint*> v1_next = v1->getNeighbours();
  vector<ftSamplePoint*> v2_next = v2->getNeighbours();
  /*  for (ki=0; ki<v1_next.size(); ++ki)
    {
      MeshVertex *pt = dynamic_cast<MeshVertex*>(v1_next[ki]);
      if (pt->up_)
	v1_next[ki] = pt->up_;
    }
  for (ki=0; ki<v1_next.size(); ++ki)
    {
      for (kj=ki+1; kj<v1_next.size(); )
	{
	  if (v1_next[ki] == v1_next[kj])
	    v1_next.erase(v1_next.begin()+kj);
	  else
	    ++kj;
	}
    }
    
  
  for (ki=0; ki<v2_next.size(); ++ki)
    {
      MeshVertex *pt = dynamic_cast<MeshVertex*>(v2_next[ki]);
      if (pt->up_)
	v2_next[ki] = pt->up_;
    }
  for (ki=0; ki<v2_next.size(); ++ki)
    {
      for (kj=ki+1; kj<v2_next.size(); )
	{
	  if (v2_next[ki] == v2_next[kj])
	    v2_next.erase(v2_next.begin()+kj);
	  else
	    ++kj;
	}
	}*/
    
  // Remove extra connections
  for (size_t ki=0; ki<v1_next.size(); )
    {
      size_t kj;
      for (kj=0; kj<v2_next.size(); ++kj)
	if (v1_next[ki] == v2_next[kj])
	  break;
      if (kj < v2_next.size())
	{
	  size_t kr = (ki==0) ? v1_next.size()-1 : ki-1;
	  size_t kh = (ki+1)%v1_next.size();
	  if (v1_next[kr] == v2 || v1_next[kh] == v2)
	    ++ki;
	  else
	    {
	      v1_next.erase(v1_next.begin()+ki);
	      v2_next.erase(v2_next.begin()+kj);
	    }
	}
      else
	++ki;
    }
  
  size_t v1_num = v1_next.size();
  size_t v2_num = v2_next.size();
  for (ki=0; ki<v1_num; ++ki)
    {
      if (next_.size() > 0 && next_[0] == v1_next[ki])
	break;
      if (v1_next[ki] == v2)
	continue;
      
      for (kj=0; kj<v2_num; ++kj)
	if (v1_next[ki] == v2_next[kj])
	  break;

      if (kj == v2_num || (ki>0 && v1_next[ki-1] == v2) ||
	  (ki==0 && v1_next[v1_next.size()-1] == v2))
	next_.push_back(v1_next[ki]);
      else
	{
	  next_.push_back(v2_next[kj]);
	  for (kr=(kj+1)%v2_num; kr!=kj; kr=(kr+1)%v2_num)
	    {
	      if (v2_next[kr] == v1)
		continue;
	      for (kh=(ki+1)%v1_num; kh!=ki; kh=(kh+1)%v1_num)
		if (v1_next[kh] == v2_next[kr])
		  break;
	      if (kh == ki)
		next_.push_back(v2_next[kr]);
	      else
		{
		  if (next_.size() > 0 && next_[0] == v1_next[kh])
		    ki = v1_num;
		  else
		    {
		      next_.push_back(v1_next[kh]);
		      ki = kh;
		    }
		  break;
		}
	    }
	}
    }

  weight_.resize(next_.size(), -1.0);
  setMeanEdgeLen();
  int stop_break = 1;
}

//===========================================================================
MeshVertex::MeshVertex(Vector3D xyz, int bnd)
  : ftSamplePoint(xyz, bnd), mean_edge_len_(0.0), area_(0.0),
    marked_(false)
//===========================================================================
{
  prev1_ = prev2_ = up_ = 0;
}

//===========================================================================
MeshVertex::~MeshVertex()
//===========================================================================
{
}

//===========================================================================
MeshVertex* MeshVertex::clone()
//===========================================================================
{
  MeshVertex* vx = new MeshVertex(xyz_, at_boundary_);
  vx->pos_ = pos_;
  vx->norm_ = norm_;
  vx->dir_ = dir_;
  vx->pos2_ = pos2_;
  vx->dir2_ = dir2_;
  vx->mean_edge_len_ = mean_edge_len_;
  vx->weight_.insert(vx->weight_.begin(), weight_.begin(), weight_.end());
  vx->area_ = area_;
  vx->prev1_ = prev1_;
  vx->prev2_ = prev2_;
  vx->up_ = up_;
  vx->next_.insert(vx->next_.end(), next_.begin(), next_.end());
  
  return vx;
}

//===========================================================================
double MeshVertex::setMeanEdgeLen()
//===========================================================================
{
  double fac = (next_.size() == 0) ? 0.0 : 1.0/(double)next_.size();
  for (size_t ki=0; ki<next_.size(); ++ki)
    mean_edge_len_ += pntDist(next_[ki]);

  return mean_edge_len_;
}

//===========================================================================
void MeshVertex::addNeighbourInfo(MeshVertex* next)
//===========================================================================
{
  size_t nn = next_.size();
  addNeighbour(next);
  if (next_.size() > nn)
    weight_.push_back(-1.0);  // Not set
}

//===========================================================================
void MeshVertex::setWeights()
//===========================================================================
{
  double eps = 1.0e-9;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      if (weight_[ki] >= 0.0)
	continue;  // Already set

      MeshVertex *pt = dynamic_cast<MeshVertex*>(next_[ki]);  // The next vertex along the edge

      // Identify the third vertex in triangles adjacent to the edge
      ftSamplePoint *pt1 = 0, *pt2 = 0;
      for (size_t kr=0; kr<next_.size(); ++kr)
	{
	  if (kr == ki)
	    continue;
	  for (size_t kj=0; kj<pt->next_.size(); ++kj)
	    {
	      if (pt->next_[kj] == next_[kr])
		{
		  if (!pt1)
		    pt1 = next_[kr];
		  else
		    {
		      pt2 = next_[kr];
		      break;
		    }
		}
	    }
	  if (pt1 && pt2)
	    break;
	}

      double cot_alpha = 0.0;
      if (pt1)
	{
	  Vector3D pt1_2 = pt1->getPoint();
	  Vector3D v1 = getPoint() - pt1_2;
	  Vector3D v2 = pt->getPoint() - pt1_2;
	  double tmp = (v1 % v2).length();
	  if (tmp > eps)
	    cot_alpha = v1*v2/tmp;
	}

      if (pt2)
	{
	  Vector3D pt2_2 = pt2->getPoint();
	  Vector3D v2 = getPoint() - pt2_2;
	  Vector3D v1 = pt->getPoint() - pt2_2;
	  double tmp = (v1 % v2).length();
	  if (tmp > eps)
	    cot_alpha += v1*v2/tmp;
	}

      cot_alpha *= 0.5;
      weight_[ki] = cot_alpha;

      for (size_t kj=0; kj<pt->next_.size(); ++kj)
	if (pt->next_[kj] == this)
	  pt->weight_[kj] = cot_alpha;
    }
}

      
//===========================================================================
void MeshVertex::computeArea()
//===========================================================================
{
  area_ = 0.0;
  if (next_.size() <= 2)
    return;  // Area zero

  Vector3D pos1 = next_[next_.size()-1]->getPoint();
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      Vector3D pos2 = next_[ki]->getPoint();
      Vector3D centr = (xyz_ + pos1 + pos2)/3.0;
      Vector3D mid1 = 0.5*(xyz_ + pos1);
      Vector3D mid2 = 0.5*(xyz_ + pos2);

      area_ += 0.5*(((xyz_-mid1)%(xyz_-centr)).length() +
		    ((xyz_-mid2)%(xyz_-centr)).length());

      pos1 = pos2;
    }
}

//===========================================================================
void MeshVertex::setAdjacencyUp()
//===========================================================================
{
  size_t ki, kj;
  for (ki=0; ki<next_.size(); ++ki)
    {
      MeshVertex *pt = dynamic_cast<MeshVertex*>(next_[ki]);
      if (pt->up_)
	next_[ki] = pt->up_;
    }
  
  for (ki=0; ki<next_.size(); ++ki)
    {
      for (kj=ki+1; kj<next_.size(); )
	{
	  if (next_[ki] == next_[kj])
	    next_.erase(next_.begin()+kj);
	  else
	    ++kj;
	}
    }
}
