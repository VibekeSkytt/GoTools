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

#ifndef _REMESH_H
#define _REMESH_H

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/MeshVertex.h"
#include "GoTools/utils/Point.h"
#include <set>

namespace Go
{
  typedef Array<int, 2> Vector2i;
  typedef std::pair<int, int> Edge;
  typedef std::pair<Edge, double> wgtEdge;
  typedef std::pair<int, int> flagEdgeIx;

  class ReMesh
  {
  public:
    ReMesh(shared_ptr<ftPointSet>& tri_sf);

    ~ReMesh();

    void updateLevelFields(int level);
    
    void updateFields();

    void extractGraph();
    
  private:
    shared_ptr<ftPointSet> tri_sf_;
    
    // To give room for a hierarchical mesh
    std::vector<std::vector<shared_ptr<MeshVertex> > > vertices_;
    std::vector<double> mean_edge_len_;
    int s0rot_;

    void initiateMeshVertices(int num_levels);

    void downSample();

    void updateOrientField(int level);

    void updatePositionField(int level);

    void transferToLevelDown(int level);

    void getSnapCandidates(std::vector<std::vector<flagEdgeIx> >& adj,
			   std::vector<Point>& vx, double thresh,
			   std::vector<std::tuple<double, int, int, int> >& cand);

    void classifyEdges(std::vector<std::vector<flagEdgeIx> >& adj,
		       std::vector<wgtEdge>& collapse,
		       std::vector<std::set<int> >& d_set,
		       double scale, double inv_scale);

    void collapseEdges(std::vector<std::vector<flagEdgeIx> >& adj,
		       std::vector<std::set<int> >& d_set,
		       std::vector<wgtEdge>& collapse,
		       std::vector<int>& nmb_collapse);

    void getDirections4(Point& dir1, Point& norm1, Point& dir2,
			Point& norm2, Point& res1, Point& res2);
    
    void getPositions4(Point& pos1, Point& norm1, Point& dir1, Point& o1,
		       Point& pos2, Point& norm2, Point& dir2, Point& o2,
		       double scale, double inv_scale, Point& res1, Point& res2);

    void getPositions4_idx(Point& pos1, Point& norm1, Point& dir1, Point& o1,
			   Point& pos2, Point& norm2, Point& dir2, Point& o2,
			   double scale, double inv_scale, Vector2i& res1,
			   Vector2i& res2, double& err);

    Point middle(Point& pos1, Point& norm1, Point& pos2, Point& norm2);

    Point round4(Point& pos2, Point& dir, Point& norm, Point& pos,
		 double scale, double inv_scale);
    
    Point floor4(Point& pos2, Point& dir, Point& norm, Point& pos,
		 double scale, double inv_scale);
    
    Vector2i floor4_idx(Point& pos2, Point& dir, Point& norm, Point& pos,
			double scale, double inv_scale);
  };
}


#endif // _REMESH_H
