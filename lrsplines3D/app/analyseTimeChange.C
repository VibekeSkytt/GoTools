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

#include "GoTools/lrsplines3D/IdentifyChanges.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/lrsplines3D/LRSpline3DBezierCoefs.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace Go;
using namespace std;

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
      std::cout << "Usage: lr_spline_vol.g2, tri_variate_data.txt, derivative_limit, change_boxes.out" << std::endl;
      return 1;
  }

  std::ifstream ifvol(argv[1]);
  std::ifstream ifpts(argv[2]);
  double der_lim = atof(argv[3]);
  std::ofstream changebox(argv[4]);
  
  ObjectHeader oh;
  oh.read(ifvol);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifvol);
  if (!vol.get())
    {
      std::cout << "Missing volume" << std::endl;
      return 1;
    }

   // Read data pointsc
  int nmb_pts;
  ifpts >> nmb_pts;

  vector<double> pc4d;

  for (int ix=0; ix!=nmb_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifpts >> p0 >> p1 >> p2 >> q0;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
    }

  std::cout << "Finished reading" << std::endl;
  //vol->expandToFullTensorProduct();
  IdentifyChanges change(vol.get(), pc4d, true);

  vector<BoundingBox> bb;
  change.getChangeDom(der_lim, bb);

  vector<vector<double> > change_pts(bb.size());
  for (size_t kj=0; kj<pc4d.size(); kj+=4)
    {
      size_t ki;
      Point curr(pc4d[kj], pc4d[kj+1], pc4d[kj+2]);
      for (ki=0; ki<bb.size(); ++ki)
	if (bb[ki].containsPoint(curr))
	  break;
      if (ki < bb.size())
	{
	  change_pts[ki].push_back(pc4d[kj]);
	  change_pts[ki].push_back(pc4d[kj+1]);
	  change_pts[ki].push_back(pc4d[kj+3]);
	}
    }

  std::ofstream ofpt("change_pts.g2");
  for (size_t ki=0; ki<change_pts.size(); ++ki)
    {
      PointCloud3D cloud(&change_pts[ki][0], change_pts[ki].size()/3);
      cloud.writeStandardHeader(ofpt);
      cloud.write(ofpt);
    }

  changebox << bb.size() << std::endl;
  for (size_t ki=0; ki<bb.size(); ++ki)
    changebox << bb[ki].low() << " " << bb[ki].high() << std::endl;

  bool viz_sub = false;
  if (viz_sub)
    {
      for (size_t ki=0; ki<bb.size(); ++ki)
	{
	  double x1 = bb[ki].low()[0];
	  double y1 = bb[ki].low()[1];
	  double t1 = bb[ki].low()[2];
	  double x2 = bb[ki].high()[0];
	  double y2 = bb[ki].high()[1];
	  double t2 = bb[ki].high()[2];
	  double fuzzy = 1.0e-8;
	  shared_ptr<LRSplineVolume> subvol(vol->subVolume(x1, y1, t1,
							   x2, y2, t2, fuzzy));

	  std::string subb("sub_bezier.bb");
	  std::string subbder("sub_bezier_der.bb");
	  LRSpline3DBezierCoefs toBb(subvol.get());
	  toBb.getBezierCoefs(0.0, 1, 2, true, 0.0, false);
	  toBb.writeToFile(subb, 0);
	  toBb.writeToFile(subbder, 1);
	  int stop_break = 1;
	}
    }

  return 0;
}

  

