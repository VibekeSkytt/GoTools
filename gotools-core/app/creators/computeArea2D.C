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
  if (argc != 11)
    {
      std::cout << "Input parameters: curve file, out file, tolerance, ll[0], ll[1], ur[0], ur[1], num1, num2, min_cell_size" << std::endl;
      exit(1);
    }

  std::ifstream infile(argv[1]);
  ALWAYS_ERROR_IF(infile.bad(), "Bad or no input filename");
  std::ofstream outfile(argv[2]);
  double tol = atof(argv[3]);
  double ll0 = atof(argv[4]);
  double ll1 = atof(argv[5]);
  double ur0 = atof(argv[6]);
  double ur1 = atof(argv[7]);
  int num1 = atoi(argv[8]);
  int num2 = atoi(argv[9]);
  double min_cell_size = atof(argv[10]);
   
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
  vector<vector<double> > quadpt;
  vector<vector<double> > bdquad;
  quad.setQuadratureInfo(quadval, weights, min_cell_size);

  Point pt1(0.0, 0.0);
  Point pt2(1.0, 0.0);
  shared_ptr<SplineCurve> lincv(new SplineCurve(pt1, 0.0, pt2, 2.0));
  double ta = lincv->startparam();
  double tb = lincv->endparam();
  double del = tb - ta;  
  double len0 = 0.0;
  for (int kb=0; kb<4; kb++)
    {
      vector<Point> der(2);
      lincv->point(der, ta + quadval[kb]*del, 1);
      double tmp = weights[kb]*del*der[1].length();
      len0 += tmp;
    }
  std::cout << "len0 : " << len0 << std::endl;

  vector<double> knots(6);
  knots[0] =  knots[1] = knots[2] = 0.0;
  knots[3] = knots[4] = knots[5] = 0.5*M_PI;
  vector<double> cfs(9);
  cfs[0] = 1;
  cfs[1] = 0;
  cfs[2] = 1;
  cfs[3] = 0.707106781186547;
  cfs[4] = 0.707106781186547;
  cfs[5] = 0.707106781186547;
  cfs[6] = 0;
  cfs[7] = 1;
  cfs[8] = 1;
  shared_ptr<SplineCurve> circseg(new SplineCurve(3, 3, &knots[0], &cfs[0], 2, true));
  double ta1 = circseg->startparam();
  double tb1 = circseg->endparam();
  del = tb1 - ta1;  
  double len1 = 0.0;
  for (int kb=0; kb<4; kb++)
    {
      vector<Point> der(2);
      circseg->point(der, ta1 + quadval[kb]*del, 1);
      double tmp = weights[kb]*del*der[1].length();
      len1 += tmp;
    }
  std::cout << "len1 : " << len1 << std::endl;

  
  int ki, kj;
  double del1 = (ur0 - ll0)/(double)num1;
  double del2 = (ur1 - ll1)/(double)num2;
  double u1, u2, v1, v2;

  std::ofstream mesh("mesh.g2");
  mesh << "410 1 0 4 255 0 0 255" << std::endl;
  mesh << num1 + num2 + 2 << std::endl;
  for (ki=0, u1=ll0; ki<=num1; ++ki, u1+=del1)
    {
      Point pt1(u1, ll1, 0.0);
      Point pt2(u1, ur1, 0.0);
      mesh << pt1 << " " << pt2 << std::endl;
    }
  for (ki=0, v1=ll1; ki<=num2; ++ki, v1+=del2)
    {
      Point pt1(ll0, v1, 0.0);
      Point pt2(ur0, v1, 0.0);
      mesh << pt1 << " " << pt2 << std::endl;
    }

  std::ofstream uncv("unresolved.g2");
  std::ofstream uncv2("shortcurves.g2");

  double area = 0.0;
  double ferenc = 0.0;
  for (kj=0, v1=ll1, v2=v1+del2; kj<num2; ++kj, v1=v2, v2+=del2)
    for (ki=0, u1=ll0, u2=u1+del1; ki<num1; ++ki, u1=u2, u2+=del1)
      {
	Point ll(u1,v1);
	Point ur(u2,v2);
	int stat = quad.cellStat(ll, ur);
	std::cout << "Cell nr " << kj*num1+ki+1 <<", status: " << stat << std::endl;


	if (stat> 0)
	  {
	    vector<vector<shared_ptr<ParamCurve> > > unresolved_cells;
	    vector<vector<shared_ptr<ParamCurve> > > short_cvs;
	    vector<double> quadpt;
	    vector<double> ptweights;
	    vector<double> bdquad;
	    vector<double> bdnorm;
	    vector<double> bdweights;
	    quad.quadrature(ll, ur, quadpt, ptweights, unresolved_cells,
			    bdquad, bdnorm, bdweights, short_cvs, stat);
  
	    outfile << "400 1 0 4 100 100 55 255" << std::endl;
	    outfile << ptweights.size() << std::endl;
	    for (size_t kr=0; kr<quadpt.size(); kr+=2)
	      {
		Point pt(quadpt[kr], quadpt[kr+1], 0.0);
		outfile << pt << std::endl;
		
		area += ptweights[kr/2];
	      }

	    for (size_t kk=0; kk<unresolved_cells.size(); ++kk)
	      for (size_t kr=0; kr<unresolved_cells[kk].size(); ++kr)
		{
		  unresolved_cells[kk][kr]->writeStandardHeader(uncv);
		  unresolved_cells[kk][kr]->write(uncv);
		}

	    if (bdquad.size() > 0)
	      {
		outfile << "400 1 0 4 155 0 100 255" << std::endl;
		outfile << bdweights.size() << std::endl;
		double tmp = 0.0;
		for (size_t kr=0; kr<bdquad.size(); kr+=2)
		  {
		    Point pt(bdquad[kr], bdquad[kr+1], 0.0);
		    outfile << pt << std::endl;
		    
		    tmp += bdweights[kr/2];
		  }
		ferenc += tmp;

		outfile << "410 1 0 4 155 0 100 255" << std::endl;
		outfile << bdweights.size() << std::endl;
		for (size_t kr=0; kr<bdquad.size(); kr+=2)
		  {
		    Point pt1(bdquad[kr], bdquad[kr+1], 0.0);
		    Point pt2(bdquad[kr]+bdnorm[kr], bdquad[kr+1]+bdnorm[kr+1], 0.0);
		    outfile << pt1 << " " << pt2 << std::endl;
		  }

		for (size_t kk=0; kk<short_cvs.size(); ++kk)
		  for (size_t kr=0; kr<short_cvs[kk].size(); ++kr)
		    {
		      short_cvs[kk][kr]->writeStandardHeader(uncv2);
		      short_cvs[kk][kr]->write(uncv2);
		    }
	    
	      }
	  }
	int stop_break = 1;
      }
  printf("Pi: %7.13f, 2Pi: %7.13f \n Area: %7.13f \n Circumference: %7.13f \n",M_PI, 2*M_PI, area, ferenc);
}

  
