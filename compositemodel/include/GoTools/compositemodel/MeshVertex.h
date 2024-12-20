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

#ifndef _MESHVERTEX_H
#define _MESHVERTEX_H

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  class MeshVertex : public ftSamplePoint
  {
  public:
    MeshVertex();

    MeshVertex(RevEngPoint* pt);
    
    MeshVertex(MeshVertex *v1, MeshVertex *v2);
    
    MeshVertex(Vector3D xyz, int bnd);
    
    ~MeshVertex();

    MeshVertex* clone();

    int getNumNeighbourIdx()
    {
      return (int)neighbour_idx_.size();
    }

    int getNeighbourIdx(int ki)
    {
      if (ki < 0 || ki >= (int)neighbour_idx_.size())
	return -1;
      return neighbour_idx_[ki];
    }

    double setMeanEdgeLen();

    double getMeanEdgeLen()
    {
      return mean_edge_len_;
    }

    void addNeighbourInfo(MeshVertex* next);

    void setWeights();  // Must be performed after all edges are transferred

    double getWeight(size_t ix)
    {
      return (ix < weight_.size()) ? weight_[ix] : -1.0;
    }

    void computeArea();    // Must be performed after all edges are transferred

    void setArea(double area)
    {
      area_ = area;
    }
    
    double getArea()
    {
      return area_;
    }
    
    Point getPos()
    {
      return Point(xyz_[0], xyz_[1], xyz_[2]);
    }
    
    Point& getPos2()
    {
      return pos2_;
    }
    
    void setPos2(Point pos2)
      {
	pos2_ = pos2;
      }
    Point& getDirection()
      {
	return dir_;
      }

    void setDirection(Point dir)
      {
	dir_ = dir;
      }

    Point& getNormal()
      {
	return norm_;
      }

    bool isMarked()
    {
      return marked_;
    }
    
    void setMarked(bool marked)
    {
      marked_ = marked;
    }

    void getPrev(MeshVertex*& prev1, MeshVertex*& prev2)
    {
      prev1 = prev1_;
      prev2 = prev2_;
    }

    void setPrev(MeshVertex* prev1, MeshVertex* prev2)
    {
      prev1_ = prev1;
      prev2_ = prev2;
    }

    void getUp(MeshVertex*& up)
    {
      up = up_;
    }

    void setUp(MeshVertex* up)
    {
      up_ = up;
    }

    void setAdjacencyUp();
    
  private:
    Point pos_, norm_, dir_;
    Point pos2_, dir2_;
    vector<int> neighbour_idx_;
    double mean_edge_len_;
    vector<double> weight_;
    double area_;
    MeshVertex *prev1_, *prev2_, *up_;
    bool marked_;
  };
    
}

#endif
