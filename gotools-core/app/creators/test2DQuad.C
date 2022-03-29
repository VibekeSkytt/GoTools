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

#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/CutCellQuad.h"
#include "GoTools/geometry/Utils.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;


int main(int argc, char** argv)
{
  if (argc != 7)
    {
      std::cout << "Input parameters: curve file, tolerance, ll[0], ll[1], ur[0], ur[1]" << std::endl;
      exit(1);
    }

  std::ifstream infile(argv[1]);
  ALWAYS_ERROR_IF(infile.bad(), "Bad or no input filename");
  double tol = atof(argv[2]);
  double ll0 = atof(argv[3]);
  double ll1 = atof(argv[4]);
  double ur0 = atof(argv[5]);
  double ur1 = atof(argv[6]);
  Point ll(ll0, ll1);
  Point ur(ur0, ur1);
  
  GoTools::init();

  vector<shared_ptr<ParamCurve> > curves;
  while (!infile.eof())
    {
      shared_ptr<ObjectHeader> header(new ObjectHeader());

      header->read(infile);
      shared_ptr<GeomObject> geom_obj(Factory::createObject(header->classType()));
      geom_obj->read(infile);
      shared_ptr<ParamCurve> cv =
	dynamic_pointer_cast<ParamCurve,GeomObject>(geom_obj);
      if (cv.get())
	curves.push_back(cv);

      Utils::eatwhite(infile);
    }

  CutCellQuad quad(curves, tol);
  int stat = quad.cellStat(ll, ur);
  std::cout << "Cell status: " << stat << std::endl;

  vector<double> quadval(4);
  quadval[0] = 0.069431844202973712388026755553595247452;
  quadval[1] = 0.33000947820757186759866712044837765640;
  quadval[2] = 0.66999052179242813240133287955162234360;
  quadval[3] = 0.93056815579702628761197324444640475255;
  vector<double> weights(4);
  weights[0] = 0.173927422568726928686531974610999703618;
  weights[1] = 0.326072577431273071313468025389000296382;
  weights[2] = 0.326072577431273071313468025389000296382;
  weights[3] = 0.173927422568726928686531974610999703618;
  double min_cell_size = 0.01;
  quad.setQuadratureInfo(quadval, weights, min_cell_size);

  vector<vector<shared_ptr<ParamCurve> > > unresolved_cells;
  vector<vector<shared_ptr<ParamCurve> > > short_cvs;
  vector<double> quadpt;
  vector<double> bdquad;
  vector<double> bdnorm;
  vector<double> ptweights;
  vector<double> bdweights;
  quad.quadrature(ll, ur, quadpt, ptweights, unresolved_cells,
		  bdquad, bdnorm, bdweights, short_cvs, stat);
  
  int stop_break = 1;
}

  
