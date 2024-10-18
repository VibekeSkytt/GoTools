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

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/BoundingBox.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  ObjectHeader head;
  head.read(file1);
  shared_ptr<SplineCurve> curve(new SplineCurve());
  curve->read(file1);

  Point pos0(3), axis0(3);
  std::cout << "Axis position: " << std::endl;
  for (int ka=0; ka<3; ++ka)
    std::cin >> pos0[ka];
  
  std::cout << "Axis direction: " << std::endl;
  for (int ka=0; ka<3; ++ka)
    std::cin >> axis0[ka];

  SweepSurfaceCreator sweep;
  double ang = 2*M_PI;
  shared_ptr<SplineSurface> surf(sweep.rotationalSweptSurface(*curve, ang, pos0, axis0));

  std::ofstream of0("rotational.g2");
  surf->writeStandardHeader(of0);
  surf->write(of0);
  
  // Sample
  int n1 = 40, n2 =  40;
  vector<shared_ptr<RevEngPoint> > points;

  BoundingBox bbox = surf->boundingBox();
  Point low = bbox.low();
  Point high = bbox.high();
  RectDomain dom = surf->parameterDomain();
  double u1 = dom.umin();
  double u2 = u1 + 0.4*(dom.umax()-u1);
  double udel = (u2 - u1)/(double)n1;  
  double v1 = dom.vmin();
  double v2 = dom.vmax();
  double vdel = (v2 - v1)/(double)n2;
  double upar, vpar;
  int ka, kb;
  std::ofstream osph("unitsphere.g2");
  osph << "400 1 0 4 100  0 155 255" << std::endl;
  osph << n1*n2 << std::endl;
  
  std::ofstream of("points.g2");
  of << "400 1 0 4 0 0 255 255" << std::endl;
  of << n1*n2 << std::endl;
  for (kb=0, vpar=v1; kb<n2; ++kb, vpar+=vdel)
    {
      for (ka=0, upar=u1; ka<n1; ++ka, upar+=udel)
	{
	  Point pos, norm;
	  surf->point(pos, upar, vpar);
	  surf->normal(norm, upar, vpar);

	  shared_ptr<RevEngPoint> rpt(new RevEngPoint(Vector3D(pos.begin()), -1));
	  rpt->addMongeInfo(norm, norm, 0.0, norm, 0.0, 0.0, 0.0);
	  points.push_back(rpt);

	  of << pos << std::endl;
	  osph << norm << std::endl;
	}
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(osph);
  sph.write(osph);
  
  vector<RevEngPoint*> points2(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    points2[ki] = points[ki].get();

  Point axis, Cx, Cy;

  vector<pair<vector<RevEngPoint*>::iterator,vector<RevEngPoint*>::iterator> > points3;
  points3.push_back(make_pair(points2.begin(), points2.end()));
  RevEngUtils::computeAxis(points3, axis, Cx, Cy);

  Point pos;
  double rad;
  RevEngUtils::computeCylPosRadius(points3, low, high,
				   axis, Cx, Cy, pos, rad);

  of0 << "400 1 0 4 255 0 0 255" << std::endl;
  of0 << "1" << std::endl;
  of0 << pos << std::endl;
  of0 << "410 1 0 4 255 0 0 255" << std::endl;
  of0 << "1" << std::endl;
  of0 << pos-axis << " " << pos+axis << std::endl;

  
  Point axis2, Cx2, Cy2, apex, apex2;
  double phi, phi2;
  RevEngUtils::coneAxis(points3, axis2, Cx2, Cy2);

  RevEngUtils::coneApex(points3, axis2, apex, phi);
  RevEngUtils::coneApex(points3, axis, apex2, phi2);
  
  int stop_break = 1;
}

