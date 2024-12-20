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

#include "GoTools/compositemodel/ReMesh.h"
#include <fstream>
#include <iostream>

using namespace Go;
using std::vector;
using std::set;

//===========================================================================
ReMesh::ReMesh(shared_ptr<ftPointSet>& tri_sf)
  :tri_sf_(tri_sf), s0rot_(4)
//===========================================================================
{
  int num_levels = 6;
  initiateMeshVertices(num_levels);
}

//===========================================================================
ReMesh::~ReMesh()
//===========================================================================
{
}

//===========================================================================
void ReMesh::initiateMeshVertices(int num_levels)
//===========================================================================
{
  vertices_.resize(num_levels);
  int num = tri_sf_->size();
  vertices_[0].resize(num);
  for (int ka=0; ka<num; ++ka)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>((*tri_sf_)[ka]);
      vertices_[0][ka] = shared_ptr<MeshVertex>(new MeshVertex(pt));
      vertices_[0][ka]->setIndex(pt->getIndex());
    }

  // Transfer edge information
  mean_edge_len_.resize(1);
  mean_edge_len_[0] = 0.0;
  double fac = 1.0/(double)num;
  for (int ka=0; ka<num; ++ka)
    {
      int num_ix = vertices_[0][ka]->getNumNeighbourIdx();
      for (int kb=0; kb<num_ix; ++kb)
	{
	  int ix = vertices_[0][ka]->getNeighbourIdx(kb);
	  vertices_[0][ka]->addNeighbourInfo(vertices_[0][ix].get());
	  int stop_break0 = 1;
	}

      double len = vertices_[0][ka]->setMeanEdgeLen();
      mean_edge_len_[0] += fac*len;
    }

  // Compute weights associated each edge
  for (int ka=0; ka<num; ++ka)
    {
      vertices_[0][ka]->setWeights();
      vertices_[0][ka]->computeArea();
    }

  downSample();
  int stop_break = 1;
}

struct edge_info
{
  MeshVertex *v1_, *v2_;
  double wgt_;

  edge_info(MeshVertex* v1, MeshVertex *v2, double ratio)
  {
    v1_ = v1;
    v2_ = v2;
    wgt_ = ratio;
  }
};
int compare_edge_info(const edge_info& e1, const edge_info& e2)
{
  return e1.wgt_ < e2.wgt_;
}
  
//===========================================================================
void ReMesh::downSample()
//===========================================================================
{
  for (size_t ki=0; ki<vertices_[0].size(); ++ki)
    {
      vector<ftSamplePoint*> next = vertices_[0][ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  vector<ftSamplePoint*> next2 = next[kj]->getNeighbours();
	  size_t kr;
	  for (kr=0; kr<next2.size(); ++kr)
	    if (next2[kr] == vertices_[0][ki].get())
	      break;
	  if (kr == next2.size())
	    std::cout << "Incompatible next pointers, init, ki=" << ki << std::endl;
	}
    }
  
  for (size_t level=1; level<vertices_.size(); ++level)
    {
      // Collect information on edge "importance"
      vector<edge_info> cand_edgs;
      for (size_t ki=0; ki<vertices_[level-1].size(); ++ki)
	{
	  Point norm1 = vertices_[level-1][ki]->getNormal();
	  double area1 = vertices_[level-1][ki]->getArea();
	  vector<ftSamplePoint*> next = vertices_[level-1][ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      MeshVertex *pt = dynamic_cast<MeshVertex*>(next[kj]);
	      if (pt->isMarked())
		continue;

	      double dp = norm1*pt->getNormal();
	      double area2 = pt->getArea();
	      double ratio = (area1 > area2) ? area1/area2 : area2/area1;
	      cand_edgs.push_back(edge_info(vertices_[level-1][ki].get(),
					    pt, dp*ratio));
	    }
	  vertices_[level-1][ki]->setMarked(true);
	}
      for (size_t ki=0; ki<vertices_[level-1].size(); ++ki)
	vertices_[level-1][ki]->setMarked(false);

      std::sort(cand_edgs.begin(), cand_edgs.end(), compare_edge_info);

      size_t nc = 0;
      size_t num_marked = 0;
      for (size_t ki=0; ki<cand_edgs.size(); ++ki)
	{
	  if (cand_edgs[ki].v1_->isMarked() || cand_edgs[ki].v2_->isMarked())
	    continue;
	  cand_edgs[ki].v1_->setMarked(true);
	  cand_edgs[ki].v2_->setMarked(true);
	  cand_edgs[nc++] = cand_edgs[ki];
	  num_marked += 2;
	}

      size_t nvx = vertices_[level-1].size() - nc;
      vertices_[level].resize(nvx);
      for (size_t ki=0; ki<nc; ++ki)
	{
	  shared_ptr<MeshVertex> vx_merge(new MeshVertex(cand_edgs[ki].v1_,
							 cand_edgs[ki].v2_));
	  vx_merge->setIndex((int)ki);
	  vertices_[level][ki] = vx_merge;
	  cand_edgs[ki].v1_->setUp(vx_merge.get());
	  cand_edgs[ki].v2_->setUp(vx_merge.get());
	}

      for (size_t ki=0, kj=nc; ki<vertices_[level-1].size(); ++ki)
	{
	  if (vertices_[level-1][ki]->isMarked())
	    continue;
	  vertices_[level][kj] =
	    shared_ptr<MeshVertex>(vertices_[level-1][ki]->clone());
	  vertices_[level][kj]->setIndex((int)kj);
	  vertices_[level][kj]->setPrev(vertices_[level-1][ki].get(),
					  vertices_[level-1][ki].get());
	  vertices_[level-1][ki]->setUp(vertices_[level][kj].get());
	  ++kj;
	}

      for (size_t ki=0; ki<vertices_[level].size(); ++ki)
	vertices_[level][ki]->setAdjacencyUp();
	  
       // Check
      for (size_t ki=0; ki<vertices_[0].size(); ++ki)
	{
	  MeshVertex *up;
	  vertices_[0][ki]->getUp(up);
	  if (!up)
	    std::cout << "Missing up pointer, ki = " << ki << std::endl;
	}
      for (size_t kj=1; kj<=level; ++kj)
	{
	  for (size_t ki=0; ki<vertices_[kj].size(); ++ki)
	    {
	      MeshVertex *prev1, *prev2;
	      vertices_[kj][ki]->getPrev(prev1, prev2);
	      if (!(prev1 && prev2))
		{
		  std::cout << "Missing previous pointer, ki = " << ki;
		  std::cout << "kj = " << kj << std::endl;
		}
	    }
	}
  
      for (size_t ki=0; ki<vertices_[level].size(); ++ki)
	{
	  vector<ftSamplePoint*> next = vertices_[level][ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      vector<ftSamplePoint*> next2 = next[kj]->getNeighbours();
	      size_t kr;
	      for (kr=0; kr<next2.size(); ++kr)
		if (next2[kr] == vertices_[level][ki].get())
		  break;
	      if (kr == next2.size())
		{
		  //std::cout << "Incompatible next pointers, ki=" << ki << std::endl;
		  vertices_[level][ki]->removeNeighbour(next[kj]);
		}
	    }
	}
      
      for (size_t ki=0; ki<vertices_[level-1].size(); ++ki)
	vertices_[level-1][ki]->setMarked(false);
      
      for (size_t ki=0; ki<vertices_[level].size(); ++ki)
	vertices_[level][ki]->setWeights();  // Not the most effective solution

    }

}

//===========================================================================
void ReMesh::updateFields()
//===========================================================================
{
  for (int level=(int)vertices_.size()-1; level>=0; --level)
    updateLevelFields(level);
}

//===========================================================================
void ReMesh::updateLevelFields(int level)
//===========================================================================
{
  int num_iter = 6;
  std::ofstream of1_0("pos_field0.g2");
  of1_0 << "400 1 0 4 100 155 0 255" << std::endl;
  of1_0 << vertices_[level].size() << std::endl;
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    of1_0 << vertices_[level][ki]->getPos2() << std::endl;
      
  std::ofstream of2_0("dir_field0.g2");
  of2_0 << "410 1 0 4 0 155 100 255" << std::endl;
  of2_0 << vertices_[level].size() << std::endl;
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    {
      double len = 0.75*vertices_[level][ki]->getMeanEdgeLen();
      Point vec = len*vertices_[level][ki]->getDirection();
      of2_0 << vertices_[level][ki]->getPos2() << " ";
      of2_0 << vertices_[level][ki]->getPos2()+vec << std::endl;
    }
  
  for (int ka=0; ka<num_iter; ++ka)
    {
      updateOrientField(level);

      updatePositionField(level);

      int stop_break = 1;
    }

  std::ofstream of1("pos_field.g2");
  of1 << "400 1 0 4 100 155 0 255" << std::endl;
  of1 << vertices_[level].size() << std::endl;
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    of1 << vertices_[level][ki]->getPos2() << std::endl;
      
  std::ofstream of2("dir_field.g2");
  of2 << "410 1 0 4 0 155 100 255" << std::endl;
  of2 << vertices_[level].size() << std::endl;
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    {
      double len = 0.75*vertices_[level][ki]->getMeanEdgeLen();
      Point vec = len*vertices_[level][ki]->getDirection();
      of2 << vertices_[level][ki]->getPos2() << " ";
      of2 << vertices_[level][ki]->getPos2()+vec << std::endl;
    }
      
 transferToLevelDown(level);
}

//===========================================================================
void ReMesh::updateOrientField(int level)
//===========================================================================
{
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    {
      double wgt_sum = 0.0;
      Point dir = vertices_[level][ki]->getDirection();
      Point norm1 = vertices_[level][ki]->getNormal();
      vector<ftSamplePoint*> next = vertices_[level][ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  double wgt = vertices_[level][ki]->getWeight(kj);
	  if (wgt <= 0.0)
	    continue;
	  
	  MeshVertex *pt = dynamic_cast<MeshVertex*>(next[kj]);

	  Point firstdir, seconddir;
	  getDirections4(dir, norm1, pt->getDirection(), pt->getNormal(),
			 firstdir, seconddir);
	  dir = wgt_sum*firstdir + wgt*seconddir;
	  Point tmp = norm1*(norm1*dir);
	  dir -= tmp;
	  wgt_sum += wgt;
	  dir.normalize_checked();
	}
    
      if (wgt_sum > 0.0)
	vertices_[level][ki]->setDirection(dir);
    }
}

//===========================================================================
void ReMesh::updatePositionField(int level)
//===========================================================================
{
  double eps = 1.0e-9;
  double scale = 0.5, inv_scale = 1.0/scale;
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    {
      double wgt_sum = 0.0;
      Point pos1 = vertices_[level][ki]->getPos();
      Point pos1_2 = vertices_[level][ki]->getPos2();
      Point norm1 = vertices_[level][ki]->getNormal();
      Point dir1 = vertices_[level][ki]->getDirection();
      
      vector<ftSamplePoint*> next = vertices_[level][ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  double wgt = vertices_[level][ki]->getWeight(kj);
	  if (wgt <= 0.0)
	    continue;
	  
	  MeshVertex *pt = dynamic_cast<MeshVertex*>(next[kj]);
	  Point pos2 = pt->getPos();
	  Point pos2_2 = pt->getPos2();
	  Point norm2 = pt->getNormal();
	  Point dir2 = pt->getDirection();
	  
	  Point firstpos, secondpos;
	  getPositions4(pos1, norm1, dir1, pos1_2, pos2, norm2, dir2, pos2_2,
			scale, inv_scale, firstpos, secondpos);
	  pos1_2 = wgt_sum*firstpos + wgt*secondpos;
	  wgt_sum += wgt;
	  if (wgt_sum > eps)
	    pos1_2 /= wgt_sum;

	  Point tmp = (norm1*(pos1_2 - pos1))*norm1;
	  pos1_2 -= tmp;
	}

      if (wgt_sum > eps)
	{
	  pos1_2 = round4(pos1_2, dir1, norm1, pos1, scale, inv_scale);
	  vertices_[level][ki]->setPos2(pos1_2);
	}
    }
}

//===========================================================================
void ReMesh::transferToLevelDown(int level)
//===========================================================================
{
  if (level == 0)
    return; // No level down
  
  for (size_t ki=0; ki<vertices_[level].size(); ++ki)
    {
      MeshVertex *v1, *v2;
      vertices_[level][ki]->getPrev(v1, v2);
      Point norm1 = v1->getNormal();
      Point norm2 = v2->getNormal();
      Point dir = vertices_[level][ki]->getDirection();
      Point dir1 = dir - (norm1*dir)*norm1;
      Point dir2 = dir - (norm2*dir)*norm2;
      v1->setDirection(dir1);
      v2->setDirection(dir2);

      Point pos = vertices_[level][ki]->getPos2();
      Point pos1 = pos - (norm1*(pos - v1->getPos()))*norm1;
      Point pos2 = pos - (norm2*(pos - v2->getPos()))*norm2;
      v1->setPos2(pos1);
      v2->setPos2(pos2);
    }
}

int compare_wgtedge(const wgtEdge& e1, const wgtEdge& e2)
{
  return e1.second < e2.second;
}
  
int compare_snap_cand(const std::tuple<double, int, int, int>& c1,
		      const std::tuple<double, int, int, int>& c2)
{
  return (std::get<0>(c1) < std::get<0>(c2));
}
  

//===========================================================================
void ReMesh::extractGraph()
//===========================================================================
{
  double scale = 0.5, inv_scale = 1.0/scale;
  vector<vector<flagEdgeIx> > adj(vertices_[0].size());
  vector<wgtEdge> collapse;
  vector<int> nmb_collapse;
  vector<set<int> > d_set(vertices_[0].size());
  classifyEdges(adj, collapse, d_set, scale, inv_scale);
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      Point pt = vertices_[0][ki]->getPos2();
      std::ofstream l1("locvx1.g2");
      l1 << "400 1 0 4 0 0 0 255" << std::endl;
      l1 << "1" << std::endl;
      l1 << pt << std::endl;
      if (d_set[ki].size() > 0)
	{
	  l1 << "400 1 0 4 0 0 255 255" << std::endl;
	  l1 << d_set[ki].size() << std::endl;
	  for (auto it=d_set[ki].begin(); it!=d_set[ki].end(); ++it)
	    l1 << vertices_[0][*it]->getPos2() << std::endl;
	}
      if (adj[ki].size() > 0)
	{
	  l1 << "410 1 0 4 255 0 0 255" << std::endl;
	  l1 << adj[ki].size() << std::endl;
	  for (size_t kj=0; kj<adj[ki].size(); ++kj)
	    l1 << pt << " " << vertices_[0][adj[ki][kj].first]->getPos2() << std::endl;
	}
      int stop_loc1 = 1;
    }
	      
  std::sort(collapse.begin(), collapse.end(), compare_wgtedge);

  collapseEdges(adj, d_set, collapse, nmb_collapse);
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      Point pt = vertices_[0][ki]->getPos2();
      std::ofstream l1("locvx2.g2");
      l1 << "400 1 0 4 0 0 0 255" << std::endl;
      l1 << "1" << std::endl;
      l1 << pt << std::endl;
      if (d_set[ki].size() > 0)
	{
	  l1 << "400 1 0 4 0 0 255 255" << std::endl;
	  l1 << d_set[ki].size() << std::endl;
	  for (auto it=d_set[ki].begin(); it!=d_set[ki].end(); ++it)
	    l1 << vertices_[0][*it]->getPos2() << std::endl;
	}
      if (adj[ki].size() > 0)
	{
	  l1 << "410 1 0 4 255 0 0 255" << std::endl;
	  l1 << adj[ki].size() << std::endl;
	  for (size_t kj=0; kj<adj[ki].size(); ++kj)
	    l1 << pt << " " << vertices_[0][adj[ki][kj].first]->getPos2() << std::endl;
	}
      int stop_loc2 = 1;
    }

  size_t nvx = 0;
  double avg_collapse = 0.0;
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if (adj[ki].size() == 0)
	continue;
      avg_collapse += (double)nmb_collapse[ki];
      ++nvx;
    }
  avg_collapse /= (double)nvx;

  // Remove spurious vertices
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if (adj[ki].size() == 0)
	continue;
      if (nmb_collapse[ki] <= (int)(avg_collapse/10.0))
	{
	  adj[ki].clear();
	  --nvx;
	}
    }

  vector<Point> tmp_edgs;
  std::ofstream of1("cvx.g2");
  of1 << "400 1 0 4 0 0 0 255" << std::endl;
  of1 << nvx << std::endl;
  for (size_t ki=0; ki<adj.size(); ++ki)
    if (adj[ki].size() > 0)
      {
	Point pt = vertices_[0][ki]->getPos2();
	of1 << pt << std::endl;
	for (size_t kj=0; kj<adj[ki].size(); ++kj)
	  {
	    if (adj[ki][kj].first < (int)ki)
	      continue;
	    tmp_edgs.push_back(pt);
	    tmp_edgs.push_back(vertices_[0][adj[ki][kj].first]->getPos2());
	  }
      }
  std::ofstream of2("cedg.g2");
  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  of2 << tmp_edgs.size()/2 << std::endl;
  for (size_t ki=0; ki<tmp_edgs.size(); ki+=2)
    of2 << tmp_edgs[ki] << " " << tmp_edgs[ki+1] << std::endl;

    
  // Compute vertex position
  vector<Point> vx_pt(adj.size(), Point(0.0, 0.0, 0.0));
  vector<Point> vx_norm(adj.size(), Point(0.0, 0.0, 0.0));
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if (adj[ki].size() == 0)
	continue;

      double wgt_sum = 0.0;
      for (auto it=d_set[ki].begin(); it!=d_set[ki].end(); ++it)
	{
	  MeshVertex *vx = vertices_[0][*it].get();
	  Point pos = vx->getPos();
	  Point pos2 = vx->getPos2();
	  Point norm = vx->getNormal();
	  double wgt = std::exp(-pos2.dist2(pos)*inv_scale*inv_scale*9);
	  vx_pt[ki] += wgt*pos2;
	  vx_norm[ki] += wgt*norm;
	  wgt_sum += wgt;
	}
      vx_pt[ki] /= wgt_sum;
      vx_norm[ki].normalize_checked();
    }
  
  vector<Point> tmp_vx3;
  vector<Point> tmp_edgs3;
  for (size_t ki=0; ki<adj.size(); ++ki)
    if (adj[ki].size() > 0)
      {
	tmp_vx3.push_back(vx_pt[ki]);
	for (size_t kj=0; kj<adj[ki].size(); ++kj)
	  {
	    if (adj[ki][kj].first < (int)ki)
	      continue;
	    tmp_edgs3.push_back(vx_pt[ki]);
	    int id = vertices_[0][adj[ki][kj].first]->getIndex();
	    tmp_edgs3.push_back(vx_pt[id]);
	  }
      }

  std::ofstream of3("cvx2.g2");
  of3 << "400 1 0 4 0 0 255 255" << std::endl;
  of3 << tmp_vx3.size() << std::endl;
  for (size_t ki=0; ki<tmp_vx3.size(); ++ki)
    of3 << tmp_vx3[ki] << std::endl;

  std::ofstream of3_2("cedg2.g2");
  of3_2 << "410 1 0 4 255 0 0 255" << std::endl;
  of3_2 << tmp_edgs3.size()/2 << std::endl;
  for (size_t ki=0; ki<tmp_edgs3.size(); ki+=2)
    of3_2 << tmp_edgs3[ki] << " " << tmp_edgs3[ki+1] << std::endl;

  double thresh = 0.3*scale;
  while (true) 
    {
      vector<std::tuple<double, int, int, int> > cand;
      getSnapCandidates(adj, vx_pt, thresh, cand);
      std::sort(cand.begin(), cand.end(), compare_snap_cand);

      bool changed = false;
      for (size_t ki=0; ki<cand.size(); ++ki)
	{
	  int ix = std::get<1>(cand[ki]);
	  int jx = std::get<2>(cand[ki]);
	  int rx = std::get<3>(cand[ki]);
	  auto e1 = std::find_if(adj[ix].begin(), adj[ix].end(),
				 [jx](const flagEdgeIx& ed){return ed.first==jx;});
	  auto e2 = std::find_if(adj[jx].begin(), adj[jx].end(),
				 [rx](const flagEdgeIx& ed){return ed.first==rx;});
	  auto e3 = std::find_if(adj[rx].begin(), adj[rx].end(),
				 [ix](const flagEdgeIx& ed){return ed.first==ix;});
	  if (e1 == adj[ix].end() || e2 == adj[jx].end())
	    continue;

	  Point p1 = vx_pt[ix];
	  Point p2 = vx_pt[jx];
	  Point p3 = vx_pt[rx];
	  double a = p2.dist(p3);
	  double b = p1.dist(p2);
	  double c = p1.dist(p3);
	  double s = 0.5*(a + b + c);
	  double height = 2.0*sqrt(s*(s-a)*(s-b)*(s-c))/a;
	  if (height != std::get<0>(cand[ki]))
	    continue;

	  if (b < thresh || c < thresh)
	    {
	      // Merge ix with jx or rx
	      int mx = (b < thresh) ? jx : rx;
	      vx_pt[ix] = 0.5*(vx_pt[ix] + vx_pt[mx]);
	      vx_norm[ix] = 0.5*(vx_norm[ix] + vx_norm[mx]);
	      vx_norm[ix].normalize_checked();

	      // Transfer edges
	      std::set<int> adj2;
	      for (size_t kh=0; kh<adj[mx].size(); ++kh)
		{
		  int hx = adj[mx][kh].first;
		  if (hx != ix)
		    {
		      adj2.insert(hx);
		      for (size_t kv=0; kv<adj[hx].size(); ++kv)
			if (adj[hx][kv].first == mx)
			  adj[hx][kv].first = ix;
		    }
		}

	      for (size_t kh=0; kh<adj[ix].size(); ++kh)
		adj2.insert(adj[ix][kh].first);
	      adj2.erase(ix);
	      adj2.erase(mx);
	      adj[mx].clear();
	      adj[ix].clear();
	      for (auto it=adj2.begin(); it!=adj2.end(); ++it)
		adj[ix].push_back(std::make_pair(*it, 0));
	    }
	  else
	    {
	      // Merge jx and rx and replace ix
	      vx_pt[ix] = 0.5*(vx_pt[jx] + vx_pt[rx]);
	      vx_norm[ix] = 0.5*(vx_norm[jx] + vx_norm[rx]);
	      vx_norm[ix].normalize_checked();

	      if (e2 != adj[jx].end())
		adj[jx].erase(e2);
	      auto e4 = std::find_if(adj[rx].begin(), adj[rx].end(),
				     [jx](const flagEdgeIx& ed){return ed.first==jx;});
	      if (e4 != adj[rx].end())
		adj[rx].erase(e4);

	      if (e3 == adj[rx].end())
		{
		  adj[ix].push_back(std::make_pair(rx,0));
		  adj[rx].push_back(std::make_pair(ix,0));
		}
	    }

	  changed = true;
	}
      if (!changed)
	break;
    }
  
  vector<Point> tmp_vx2;
  vector<Point> tmp_edgs2;
  for (size_t ki=0; ki<adj.size(); ++ki)
    if (adj[ki].size() > 0)
      {
	Point pt = vertices_[0][ki]->getPos2();
	tmp_vx2.push_back(pt);
	for (size_t kj=0; kj<adj[ki].size(); ++kj)
	  {
	    if (adj[ki][kj].first < (int)ki)
	      continue;
	    tmp_edgs2.push_back(pt);
	    tmp_edgs2.push_back(vertices_[0][adj[ki][kj].first]->getPos2());
	  }
      }
  std::ofstream of4("cvx3.g2");
  of4 << "400 1 0 4 0 0 0 255" << std::endl;
  of4 << tmp_vx2.size() << std::endl;
  for (size_t ki=0; ki<tmp_vx2.size(); ++ki)
    of4 << tmp_vx2[ki] << std::endl;
  std::ofstream of5("cedg3.g2");
  of5 << "410 1 0 4 255 0 0 255" << std::endl;
  of5 << tmp_edgs2.size()/2 << std::endl;
  for (size_t ki=0; ki<tmp_edgs2.size(); ki+=2)
    of5 << tmp_edgs2[ki] << " " << tmp_edgs2[ki+1] << std::endl;
  int stop_break = 1;
}


//===========================================================================
void ReMesh::getSnapCandidates(vector<vector<flagEdgeIx> >& adj,
			       vector<Point>& vx, double thresh,
			       vector<std::tuple<double, int, int, int> >& cand)
//===========================================================================
{
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if (adj[ki].size() == 0)
	continue;
      Point p1 = vx[ki];
      for (size_t kj=0; kj<adj[ki].size(); ++kj)
	{
	  int jx = adj[ki][kj].first;
	  Point p2 = vx[jx];
	  for (size_t kr=0; kr<adj[jx].size(); ++kr)
	    {
	      int rx = adj[jx][kr].first;
	      if (rx == ki)
		continue;
	      
	      Point p3 = vx[rx];
	      double a = p2.dist(p3);
	      double b = p1.dist(p2);
	      double c = p1.dist(p3);
	      if (a > std::max(b, c))
		{
		  // Diagonal
		  double s = 0.5*(a + b + c);
		  double height = 2.0*sqrt(s*(s-a)*(s-b)*(s-c))/a;
		  if (height < thresh)
		    cand.push_back(std::make_tuple(height, (int)ki, jx, rx));
		}
	    }
	}
    }
}

//===========================================================================
void ReMesh::classifyEdges(vector<vector<flagEdgeIx> >& adj,
			   vector<wgtEdge>& collapse,
			   vector<set<int> >& d_set,
			   double scale, double inv_scale)
//===========================================================================
{
  for (size_t ki=0; ki<vertices_[0].size(); ++ki)
    {
      int idx1 = vertices_[0][ki]->getIndex();
      d_set[ki].insert(idx1);
      Point pos1 = vertices_[0][ki]->getPos();
      Point pos1_2 = vertices_[0][ki]->getPos2();
      Point dir1 = vertices_[0][ki]->getDirection();
      Point norm1 = vertices_[0][ki]->getNormal();
      vector<ftSamplePoint*> next = vertices_[0][ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  MeshVertex *pt = dynamic_cast<MeshVertex*>(next[kj]);
	  int idx2 = pt->getIndex();
	  if (idx2 < idx1)
	    continue;
	  
	  Point pos2 = pt->getPos();
	  Point pos2_2 = pt->getPos2();
	  Point norm2 = pt->getNormal();
	  Point dir2 = pt->getDirection();
	  
	  Point firstdir, seconddir;
	  getDirections4(dir1, norm1, dir2, norm2, firstdir, seconddir);

	  double err = 0.0;
	  Vector2i firstpos_ix, secondpos_ix;
	  getPositions4_idx(pos1, norm1, firstdir, pos1_2, pos2, norm2,
			    seconddir, pos2_2,
			    scale, inv_scale, firstpos_ix, secondpos_ix, err);
	  Vector2i diff = firstpos_ix - secondpos_ix;
	  diff[0] = abs(diff[0]);
	  diff[1] = abs(diff[1]);
	  if (std::max(diff[0],diff[1]) > 1 || diff == Vector2i(1, 1))
	    continue;

	  if (diff[0]+diff[1] == 0)
	    {
	      collapse.push_back(std::make_pair(std::make_pair(idx1, idx2), err));
	      d_set[idx1].insert(idx2);
	      d_set[idx2].insert(idx1);
	    }
	  else
	    adj[idx1].push_back(std::make_pair(idx2, 0));
	}
      int stop_break1 = 1;
    }
  int stop_break2 = 1;
}


//===========================================================================
void ReMesh::collapseEdges(vector<vector<flagEdgeIx> >& adj,
			   vector<set<int> >& d_set,
			   vector<wgtEdge>& collapse,
			   vector<int>& nmb_collapse)
//===========================================================================
{
  nmb_collapse.resize(vertices_[0].size(), 0);
  for (size_t ki=0; ki<collapse.size(); ++ki)
    {
      Edge edge = collapse[ki].first;  // Two indices
      if (edge.first > edge.second)
	std::swap(edge.first, edge.second);

      set<flagEdgeIx> merged(adj[edge.first].begin(), adj[edge.first].end());
      merged.insert(adj[edge.second].begin(), adj[edge.second].end());
      int ix = (d_set[edge.second].size() > d_set[edge.first].size()) ?
	edge.first : edge.second;
      adj[edge.first].clear();
      adj[edge.second].clear();
      adj[ix].insert(adj[ix].end(), merged.begin(), merged.end());
      int n_c = nmb_collapse[edge.first] + nmb_collapse[edge.second];
      nmb_collapse[ix] = n_c + 1;
      set<int> tmp_d(d_set[edge.first].begin(), d_set[edge.first].end());
      tmp_d.insert(d_set[edge.second].begin(), d_set[edge.second].end());
      d_set[ix] = tmp_d;
      adj[edge.first].shrink_to_fit();
      adj[edge.second].shrink_to_fit();
    }
  int stop_break = 1;
}

//===========================================================================
void ReMesh::getDirections4(Point& dir1, Point& norm1, Point& dir2,
			    Point& norm2, Point& res1, Point& res2)
//===========================================================================
{
  Point A[2] = {dir1, norm1.cross(dir1)};
  Point B[2] = {dir2, norm2.cross(dir2)};

  double max_val = std::numeric_limits<double>::lowest();
  int ix1 = -1, ix2 = -1;

  for (int ka=0; ka<2; ++ka)
    for (int kb=0; kb<2; ++kb)
      {
	double val = fabs(A[ka]*B[kb]);
	if (val > max_val)
	  {
	    ix1 = ka;
	    ix2 = kb;
	    max_val = val;
	  }
      }

  int sgn = (A[ix1]*B[ix2] >= 0.0) ? 1 : -1;
  res1 = A[ix1];
  res2 = sgn*B[ix2];
}

//===========================================================================
void ReMesh::getPositions4(Point& pos1, Point& norm1, Point& dir1, Point& o1,
			   Point& pos2, Point& norm2, Point& dir2, Point& o2,
			   double scale, double inv_scale, Point& res1, Point& res2)
//===========================================================================
{
  Point vec1 = norm1.cross(dir1);
  Point vec2 = norm2.cross(dir2);
  Point mid = middle(pos1, norm1, pos2, norm2);
  Point o1p = floor4(o1, dir1, norm1, mid, scale, inv_scale);
  Point o2p = floor4(o2, dir2, norm2, mid, scale, inv_scale);
  
  double min_val = std::numeric_limits<double>::max();
  int ix1 = -1, ix2 = -1;

  for (int ka=0; ka<4; ++ka)
    {
      Point o1t = o1p + (dir1*(ka&1) + vec1*((ka&2)>>1))*scale;
      for (int kb=0; kb<4; ++kb)
	{
	  Point o2t = o2p + (dir2*(kb&1) + vec2*((kb&2)>>1))*scale;
	  double val = (o1t - o2t).length2();
	  if (val < min_val)
	    {
	      ix1 = ka;
	      ix2 = kb;
	      min_val = val;
	    }
	}
    }

  res1 = o1p + (dir1*(ix1&1) + vec1*((ix1&2)>>1))*scale;
  res2 = o2p + (dir2*(ix2&1) + vec2*((ix2&2)>>1))*scale;
}


//===========================================================================
void ReMesh::getPositions4_idx(Point& pos1, Point& norm1, Point& dir1, Point& o1,
			       Point& pos2, Point& norm2, Point& dir2, Point& o2,
			       double scale, double inv_scale, Vector2i& res1,
			       Vector2i& res2, double& err)
//===========================================================================
{
  Point vec1 = norm1.cross(dir1);
  Point vec2 = norm2.cross(dir2);
  Point mid = middle(pos1, norm1, pos2, norm2);
  Vector2i o1p = floor4_idx(o1, dir1, norm1, mid, scale, inv_scale);
  Vector2i o2p = floor4_idx(o2, dir2, norm2, mid, scale, inv_scale);
  
  double min_val = std::numeric_limits<double>::max();
  int ix1 = -1, ix2 = -1;

  for (int ka=0; ka<4; ++ka)
    {
      Point o1t = o1 + (dir1*((ka&1) + o1p[0]) +
			vec1*(((ka&2)>>1) + o1p[1]))*scale;
      for (int kb=0; kb<4; ++kb)
	{
	  Point o2t = o2 + (dir2*((kb&1) + o2p[0]) +
			    vec2*(((kb&2)>>1) + o2p[1]))*scale;
	  double val = (o1t - o2t).length2();
	  if (val < min_val)
	    {
	      ix1 = ka;
	      ix2 = kb;
	      min_val = val;
	    }
	}
    }

  err = min_val;
  res1 = Vector2i((ix1&1) + o1p[0], ((ix1&2) >> 1) + o1p[1]);
  res2 = Vector2i((ix2&1) + o2p[0], ((ix2&2) >> 1) + o2p[1]);
}


//===========================================================================
Point ReMesh::middle(Point& pos1, Point& norm1, Point& pos2, Point& norm2)
//===========================================================================
{
  double eps = 1.0e-4;
  double n1p1 = norm1*pos1;
  double n1p2 = norm1*pos2;
  double n2p1 = norm2*pos1;
  double n2p2 = norm2*pos2;
  double n1n2 = norm1*norm2;

  double div = 1.0 - n1n2*n1n2 + eps;
  double lambda1 = 2.0*(n1p2 - n1p1 - n1n2*(n2p1 - n2p2))/div;
  double lambda2 = 2.0*(n2p1 - n2p2 - n1n2*(n1p2 - n1p1))/div;
  Point mid = 0.5*(pos1+pos2) - 0.25*(norm1*lambda1 + norm2*lambda2);

  return mid;
}

//===========================================================================
Point ReMesh::round4(Point& pos2, Point& dir, Point& norm, Point& pos,
		     double scale, double inv_scale)
//===========================================================================
{
  Point tmp = norm.cross(dir);
  Point dvec = pos - pos2;
  Point res = pos2 + scale*dir*std::round(inv_scale*(dir*dvec)) +
    scale*tmp*std::round(inv_scale*(tmp*dvec));

  return res;
}

//===========================================================================
Point ReMesh::floor4(Point& pos2, Point& dir, Point& norm, Point& pos,
		     double scale, double inv_scale)
//===========================================================================
{
  Point tmp = norm.cross(dir);
  Point dvec = pos - pos2;
  Point res = pos2 + scale*dir*std::floor(inv_scale*(dir*dvec)) +
    scale*tmp*std::floor(inv_scale*(tmp*dvec));

  return res;
}

//===========================================================================
Vector2i ReMesh::floor4_idx(Point& pos2, Point& dir, Point& norm, Point& pos,
			    double scale, double inv_scale)
//===========================================================================
{
  Point tmp = norm.cross(dir);
  Point dvec = pos - pos2;
  Vector2i res((int)std::floor((dir*dvec)*inv_scale),
	       (int)std::floor((tmp*dvec)*inv_scale));
  return res;
}
