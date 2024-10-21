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

#ifndef _REVERSEENGINEERING_DOXYMAIN_H
#define _REVERSEENGINEERING_DOXYMAIN_H

/**
\page reverseengineering_doc Reverse engineering, from triangulated point cloud to CAD model.

The engine in the reverse engineering functionality is the class
\link Go::RevEng RevEng\endlink. The class is instanciated with a triangulated surface 
represented as an \link Go::ftPointSet ftPointSet\endlink, which represents a triangulation
in GoTools. Assume that the triangulation
is represented as a collection of vertices and triangles in the format
\verbatim
  vector<double> vertices;
  vector<int> triangles;
\endverbatim
The points are stored consequetive in vertices as (x,y,z) while three indices in triangles
represents one triangle. The triangles are expected to have a consistent orientation. The
vertices and triangles are transferred to ftPointSet with the following code snippet:
\verbatim
  size_t nvertices = vertices.size()/3;
  size_t ntriangles = triangles.size()/3;
  shared_ptr<ftPointSet> tri_sf = shared_ptr<ftPointSet>(new ftPointSet());
  for (size_t ki=0; ki<nvertices; ++ki)
    {
      Vector3D xyz(vertices[3*ki], vertices[3*ki+1], vertices[3*ki+2]);
      shared_ptr<RevEngPoint> currpt(new RevEngPoint(xyz, -1));
      currpt->setIndex(ki);
      tri_sf->addEntry(currpt);
    }

  for (size_t ki=0; ki<ntriangles; ++ki)
    {
      ftSamplePoint* pt1 = (*tri_sf)[triangles[3*ki]];
      ftSamplePoint* pt2 = (*tri_sf)[triangles[3*ki+1]];
      ftSamplePoint* pt3 = (*tri_sf)[triangles[3*ki+2]];
      pt1->addTriangle(pt2, pt3);
      pt2->addTriangle(pt3, pt1);
      pt3->addTriangle(pt1, pt2);
    }
\endverbatim
\link Go::RevEngPoint RevEngPoint\endlink represents a triangulation vertex and is inherited
from \link Go::ftSamplePoint ftSamplePoint\endlink. RevEngPoint is enhanced with information
such as estimated surface normal and curvature as well as associated functionality.

The reverse engineering process is organized as a sequence of operations that together 
consistute a work flow. The process is as follows:

 * <ol>
 * <li> Enhance points
 * <li> Classify points according to Gauss and mean curvature
 * <li> Segment point cloud into regions
 * <li> Surface creation
 * <li> Compute global properties such as main axes and update surfaces accordingly
 * <li> Edge creation
 * <li> Define blend surfaces
 * <li> Trim surfaces with respect to identified edges, blend surfaces and adjacent regions
 * <li> Create CAD model 

Point 4 to 6 are repeated three times, each time differently. The process
is automated, but organized as a sequence of commands to RevEng. This allows for storing
the state at a number of locations to resume the computation at a convenient time. Note that
storing and reading the state can be time consuming. In the following, we will describe
the process in some detail.

The first function to call is RevEng::enhancePoints. The points are approximated by a surface
in a local neighbourhood. Surface normal and principal curvature estimates are 
computed from this surface. Approximation errors are registered and used to set an
approximation tolerance for the proceeding computations. An additional surface normal is
computed from the triangulation. The two versions of the surface normal have different
pros and cons, and both are used in the computations.

Classification is performed in RevEng::classifyPoints. It is based on the size and
sign of estimated Gauss and mean curvature in the points. Very small curvature values are
deciphered as zero. A small curvature radius compared to the average distance between
triangle vertices indicates that the point is a part of an edge. As the expected typical 
measured objects has rounded edges is
further edge detection not a priority topic in the current version of the reverse engineering
functionality.

Next, the approximation tolerance is set by the call 
RevEng::setApproximationTolerance based on information from the preceeding computations. 
Alternatively, the application can use the
function RevEng::setApproxTol(double tol) if more control is preferred. As the given
point cloud is expected to be noisy, it is only required that a majority points associated 
to a surface will be fit by the surface within this tolerance in addition to requirements
on the average approximation error and normal direction.

RevEng::segmentIntoRegions collects connected groups of points with the same classification.
Identified edge points are excluded. Each group is stored in an instance of
\link Go::RevEngRegion RevEngRegion\endlink. 
*/
#endif // _REVERSEENGINEERING_DOXYMAIN_H
