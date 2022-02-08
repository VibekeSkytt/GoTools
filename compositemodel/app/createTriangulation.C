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

#include "GoTools/utils/Point.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <GoTools/compositemodel/ttlTriang.h>
#include <GoTools/compositemodel/ttlPoint.h>
#include <GoTools/compositemodel/ftPointSet.h>
#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h> // For atof()

#define DEBUG

using namespace std;
using namespace Go;
using Go::Point;


int main( int argc, char* argv[] )
{
  if (argc != 5) {
    std::cout << "Input parameters : Input triangulation point cloud (.g2), output points(.g2), output edges(.g2), max length boundary"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream input(argv[1]);
  std::ofstream os1(argv[2]);
  std::ofstream os2(argv[3]);
  double maxlen = atoi(argv[4]);

    // Read points
  ObjectHeader header;
  PointCloud3D points;
  try {
    header.read(input);
    points.read(input);
  }
  catch (...)
    {
      std::cout << "ERROR: Not a valid point file" << std::endl;
      return 1;
    }
  int nmb_pts = points.numPoints();
  vector<double> pts(points.rawData(), points.rawData()+3*nmb_pts);

  shared_ptr<ftPointSet> pointset(new ftPointSet());
  // Define nodes

  for (int ki=0; ki<nmb_pts; ++ki)
    {
      Vector3D curr(pts[3*ki], pts[3*ki+1], pts[3*ki+2]);
      shared_ptr<ftSamplePoint> sample(new ftSamplePoint(curr, 0));
      (void)pointset->addEntry(sample);
    }
  vector<ttlPoint*> ttlpoints;
  for (int ki=0; ki<pointset->size(); ++ki)
    {
      ftSamplePoint *ftpt = (*pointset)[ki];
      Vector3D xyz = ftpt->getPoint();
      ttlpoints.push_back(new ttlPoint(ftpt, xyz[0], xyz[1], xyz[2]));
    }
  hetriang::Triangulation triang;
  triang.createDelaunay(ttlpoints.begin(), ttlpoints.end());
  triang.optimizeDelaunay();

  //pointset->write(os2);
  list<hetriang::Edge*>* edgs = triang.getEdges();
  vector<Point> node1;
  vector<Point> node2;
  for (auto it = (*edgs).begin(); it != (*edgs).end(); ++it)
    {
      shared_ptr<hetriang::Node> source = (*it)->getSourceNode();
      shared_ptr<hetriang::Node> target = (*it)->getTargetNode();
      Point pt1(source->x(), source->y(), source->z());
      Point pt2(target->x(), target->y(), target->z());
      double len = pt1.dist(pt2);
      hetriang::Edge* twin = (*it)->getTwinEdge();
      if (len > maxlen && (!twin))
	continue;
      node1.push_back(pt1);
      node2.push_back(pt2);
    }
  os2 << "410 1 0 0" << std::endl;
  os2 << node1.size() << std::endl;
  os1 << "400 1 0 4 255 0 0 255" << std::endl;
  os1 << node1.size() << std::endl;
  for (size_t kj=0; kj<node1.size(); ++kj)
    {
      os1 << node1[kj] << std::endl;
      os2 << node1[kj] << " " << node2[kj] << std::endl;
    }
//hed::Edge* startedg = (*edgs)->begin();
  int stop_break = 1;
  // streamsize prev1 = os1.precision(15);
  // streamsize prev2 = os3.precision(15);
  // streamsize prev3 = os3.precision(15);

  // os1 << "400 1 0 0" << std::endl;
  // os1 << nvertices << std::endl;
  // for (ki=0; ki<nvertices; ++ki)
    
  // for (ki=0; ki<nvertices; ++ki)
  //     os1 << vertices[3*ki] << " " << vertices[3*ki+1] << " " << vertices[3*ki+2] << std::endl;
  // os2 << "410 1 0 0" << std::endl;
  // os2 << 3*ntriangles << std::endl;
  // for (ki=0; ki<ntriangles; ++ki)
  //   {
  //     int kj = triangles[3*ki];
  //     int kr = triangles[3*ki+1];
  //     int kh = triangles[3*ki+2];
  //     os2 << vertices[3*kj] << " " << vertices[3*kj+1] << " " << vertices[3*kj+2] << " ";
  //     os2 << vertices[3*kr] << " " << vertices[3*kr+1] << " " << vertices[3*kr+2] << std::endl;
  //     os2 << vertices[3*kj] << " " << vertices[3*kj+1] << " " << vertices[3*kj+2] << " ";
  //     os2 << vertices[3*kh] << " " << vertices[3*kh+1] << " " << vertices[3*kh+2] << std::endl;
  //     os2 << vertices[3*kh] << " " << vertices[3*kh+1] << " " << vertices[3*kh+2] << " ";
  //     os2 << vertices[3*kr] << " " << vertices[3*kr+1] << " " << vertices[3*kr+2] << std::endl;
  //   }

  // os3 << "410 1 0 0" << std::endl;
  // os3 << 3*ntriangles << std::endl;
  // for (ki=0; ki<ntriangles; ++ki)
  //   {
  //     int kj = triangles[3*ki];
  //     int kr = triangles[3*ki+1];
  //     int kh = triangles[3*ki+2];
  //     Point pt1(vertices[3*kj],vertices[3*kj+1],vertices[3*kj+2]);
  //     Point pt2(vertices[3*kr],vertices[3*kr+1],vertices[3*kr+2]);
  //     Point pt3(vertices[3*kh],vertices[3*kh+1],vertices[3*kh+2]);
  //     Point dir1 = pt2 - pt1;
  //     Point dir2 = pt3 - pt1;
  //     Point dir3 = pt3 - pt2;
  //     Point norm1 = dir1.cross(dir2);
  //     os3 << pt1 << " " << pt1+norm1 << std::endl;
  //     Point norm2 = dir3.cross(dir1);
  //     os3 << pt2 << " " << pt2-norm2 << std::endl;
  //     Point norm3 = dir3.cross(dir2);
  //     os3 << pt3 << " " << pt3-norm3 << std::endl;
  //   }
   
}
