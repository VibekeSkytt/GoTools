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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/Utils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: point cloud in (.g2), 2D curves in (.g2), point cloud out(g2) " << std::endl;
    return -1;
  }

  std::ifstream pointsin(argv[1]);
  std::ifstream cvsin(argv[2]);
  std::ofstream fileout(argv[3]);

  ObjectHeader header;
  header.read(pointsin);
  PointCloud3D points;
  points.read(pointsin);

  BoundingBox box1 = points.boundingBox();
  printf("Point domain: [ %13.3f , %13.3f ] x [ %13.3f , %13.3f ] \n",
	 box1.low()[0], box1.high()[0], box1.low()[1] ,box1.high()[1]);

  vector<shared_ptr<ParamCurve> > cvs;
  BoundingBox box2;
  while (!cvsin.eof())
    {
      header.read(cvsin);
      shared_ptr<SplineCurve> tmpcv(new SplineCurve());
      tmpcv->read(cvsin);
      if (cvs.size() == 0)
	box2 = tmpcv->boundingBox();
      else
	{
	  BoundingBox box3 = tmpcv->boundingBox();
	  box2.addUnionWith(box3);
	}
      cvs.push_back(tmpcv);
      Utils::eatwhite(cvsin);
    }
  printf("Curve domain: [ %13.3f , %13.3f ] x [ %13.3f , %13.3f ] \n",
	 box2.low()[0], box2.high()[0], box2.low()[1] ,box2.high()[1]);

  double u1 = box2.low()[0];
  double u2 = box2.high()[0];
  double v1 = box2.low()[1];
  double v2 = box2.high()[1];

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  vector<double> data2;
  double eps = 1.0e-4;
  shared_ptr<CurveLoop> loop(new CurveLoop(cvs, eps));
  CurveBoundedDomain cvdom(loop);
  
  // Sort the points according to the u-parameter
  qsort(&data[0], nmb_pts, 3*sizeof(double), compare_u_par);

  // Traverse points
  int pp0, pp1;
  for (pp0=0; pp0<(int)data.size() && data[pp0]<u1; pp0+=3);
  for (pp1=pp0; pp1<(int)data.size() && data[pp1]<u2; pp1+=3);

  // Sort the current sub set of points according to the v-parameter
  qsort(&data[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_v_par);

  // Identify sub set of points
  int pp2, pp3;
  for (pp2=pp0; pp2<pp1 && data[pp2+1]<v1; pp2+=3);
  for (pp3=pp2; pp3<pp1 && data[pp3+1]<v2; pp3+=3);

  // Split sub set
  int num = 10;
  // Sort the points according to the u-parameter
  qsort(&data[0]+pp2, (pp3-pp2)/3, 3*sizeof(double), compare_u_par);
  int pp4 = pp2, pp5 = pp2;
  double udel = (u2 - u1)/(double)num;
  double ustop = u1 + udel;
  double vdel = (v2 - v1)/(double)num;
  for (int ki=0; ki<num; ++ki)
    {
      for (; pp5<pp3 && data[pp5]<ustop; pp5+=3);
      qsort(&data[0]+pp4, (pp5-pp4)/3, 3*sizeof(double), compare_v_par);
      int pp6=pp4, pp7=pp4;
      double vstop = v1 + vdel;
      for (int kj=0; kj<num; ++kj)
	{
	  for (; pp7<pp5 && data[pp7+1]<vstop; pp7+=3);

	  // Define parameter box as spline curves
	  vector<shared_ptr<SplineCurve> > pcvs(4);
	  pcvs[0] = shared_ptr<SplineCurve>(new SplineCurve(Point(ustop-udel,vstop-vdel), Point(ustop,vstop-vdel)));
	  pcvs[1] = shared_ptr<SplineCurve>(new SplineCurve(Point(ustop,vstop-vdel), Point(ustop,vstop)));
	  pcvs[2] = shared_ptr<SplineCurve>(new SplineCurve(Point(ustop-udel,vstop-vdel), Point(ustop-vdel,vstop)));
	  pcvs[3] = shared_ptr<SplineCurve>(new SplineCurve(Point(ustop-udel,vstop), Point(ustop,vstop)));

	  // Check overlap with curve domain
	  bool intersect = false;
	  for (int ka=0; ka<4; ++ka)
	    {
	      if (cvdom.doIntersect(*pcvs[ka], eps))
		{
		  intersect = true;
		  break;
		}
	    }
	  if (intersect)
	    {
	      // Check all points
	      for (int pp8=pp6; pp8<pp7; pp8+=3)
	      {
		if (cvdom.isInDomain(Vector2D(data[pp8], data[pp8+1]), eps))
		  data2.insert(data2.end(), &data[pp8], &data[pp8+3]);
	      }
	    }
	  else
	    {
	      // Check midpoint of domain
	      if (cvdom.isInDomain(Vector2D(ustop-0.5*udel,vstop-0.5*vdel),eps))
		data2.insert(data2.end(), &data[pp6], &data[pp7]);
	    }
	  pp6 = pp7;
	  vstop += vdel;
	}
      pp4 = pp5;
      ustop += udel;
    }
    

  // Collect output
  PointCloud3D points2(&data2[0], data2.size()/3);


  // Write to output
  points2.writeStandardHeader(fileout);
  points2.write(fileout);
}

