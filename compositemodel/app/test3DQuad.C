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
#include "GoTools/compositemodel/CutCellQuad3D.h"
#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/Utils.h"
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;
using std::pair;
using std::string;

int compare(const char *str1, char str2[][8], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}


int main(int argc, char** argv)
{
  if (!(argc == 14 || argc == 15 || argc == 16))
    {
      std::cout << "Input parameters: model file, out file, tolerance, ll[0], ll[1], ll[2], ur[0], ur[1], ur[2], num1, num2, num3, min_cell_size <first cell> <last cell>" << std::endl;
      exit(1);
    }

  char* infile(argv[1]);
  std::ofstream outfile(argv[2]);
  double tol = atof(argv[3]);
  double ll0 = atof(argv[4]);
  double ll1 = atof(argv[5]);
  double ll2 = atof(argv[6]);
  double ur0 = atof(argv[7]);
  double ur1 = atof(argv[8]);
  double ur2 = atof(argv[9]);
  int num1 = atoi(argv[10]);
  int num2 = atoi(argv[11]);
  int num3 = atoi(argv[12]);
  double min_cell_size = atof(argv[13]);
  int firstcell = 0;
  int lastcell = -1;
  if (argc >= 15)
    firstcell = atoi(argv[14]);
  if (argc == 16)
    lastcell = atoi(argv[15]);
   
  // Check file type
  // Find file extension
  char *loc;
  char *last = 0;
  loc = strchr(infile, '.');
  while (loc != NULL)
    {
      last = loc;
      loc = strchr(loc+1, '.');
    }
  if (last == NULL)
    {
      std::cout << "Missing file extension of input file" << std::endl;
      exit(1);
    }
  char *input_type = last+1;

  // Check type of input file
  char keys[2][8] = {"g2", "g22"};
  int type_in;
  try {
    type_in = compare(input_type, keys, 2);
  }
  catch (...)
    {
      std::cout << "Invalide file extension of input file" << std::endl;
      exit(1);
    }
  if (type_in < 0)
    {
      std::cout << "Not a valid file type" << std::endl;
      exit(1);
    }

  double gap, neighbour, kink;
  shared_ptr<SurfaceModel> sfmodel;
  vector<shared_ptr<SurfaceModel> > voids;
  int material_id = -1;
  if (type_in == 1)
    {
      CompositeModelFileHandler fileread;
      shared_ptr<Body> body = fileread.readBody(infile);
      if (!body.get())
	{
	  sfmodel = fileread.readShell(infile);
	  if (!sfmodel.get())
	    exit(1);
	}
      else
	{
	  int nmb = body->nmbOfShells();
	  if (nmb == 1)
	    sfmodel = body->getOuterShell();
	  else
	    {
	      vector<shared_ptr<SurfaceModel> > all_shells = 
		body->getAllShells();		
	      sfmodel = all_shells[0];
	      voids.insert(voids.end(), all_shells.begin()+1,
			   all_shells.end());
	    }
	  material_id = body->getMaterial();
	}
      tpTolerances top = sfmodel->getTolerances();
      gap = top.gap;
      neighbour = top.neighbour;
      kink = top.kink;
    }
  else
    {
      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      double gap = 0.00001;
      double neighbour = 0.0001;
      double kink = 0.01;
      double approxtol = gap;

      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

      std::ifstream is(infile);
      CompositeModel *model = factory.createFromG2(is);

      sfmodel = 
	shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
      if (!sfmodel.get())
	{
	  std::cout << "No input model read" << std::endl;
	  exit(1);
	}
 
      if (sfmodel->nmbBoundaries() > 0)
	{
	  std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
	  exit(1);
	}
      
#ifdef DEBUG
      bool isOK = sfmodel->checkShellTopology();
      std::cout << "Shell topology: " << isOK << std::endl;
#endif
    }

  BoundingBox bb = sfmodel->boundingBox();
  Point low = bb.low();
  Point high = bb.high();
  std::cout << "Bounding box: " << low << " " << high << std::endl;

  CutCellQuad3D quad(sfmodel, tol);
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
  quad.setQuadratureInfo(quadval, weights, min_cell_size);

  
  int ki, kj, kh;
  double del1 = (ur0 - ll0)/(double)num1;
  double del2 = (ur1 - ll1)/(double)num2;
  double del3 = (ur2 - ll2)/(double)num3;
  double u1, u2, v1, v2, w1, w2;

  std::ofstream mesh("mesh.g2");
  mesh << "410 1 0 4 255 0 0 255" << std::endl;
  mesh << (num1+1)*(num2+num3+2) + (num2+1)*(num1+num3+2) + (num3+1)*(num1+num2+2) << std::endl;
  for (ki=0, u1=ll0; ki<=num1; ++ki, u1+=del1)
    {
      for (kj=0, v1=ll1; kj<=num2; ++kj, v1+=del2)
	{
	  Point pt1(u1, v1, ll2);
	  Point pt2(u1, v1, ur2);
	  mesh << pt1 << " " << pt2 << std::endl;
	}
      for (kj=0, w1=ll2; kj<=num3; ++kj, w1+=del3)
	{
	  Point pt1(u1, ll1, w1);
	  Point pt2(u1, ur1, w1);
	  mesh << pt1 << " " << pt2 << std::endl;
	}
    }
  for (ki=0, v1=ll1; ki<=num2; ++ki, v1+=del2)
    {
      for (kj=0, u1=ll0; kj<=num1; ++kj, u1+=del1)
	{
	  Point pt1(u1, v1, ll2);
	  Point pt2(u1, v1, ur2);
	  mesh << pt1 << " " << pt2 << std::endl;
	}
       for (kj=0, w1=ll2; kj<=num3; ++kj, w1+=del3)
	{
	  Point pt1(ll0, v1, w1);
	  Point pt2(ur0, v1, w1);
	  mesh << pt1 << " " << pt2 << std::endl;
	}
   }
  for (ki=0, w1=ll2; ki<=num3; ++ki, w1+=del3)
    {
      for (kj=0, u1=ll0; kj<=num1; ++kj, u1+=del1)
	{
	  Point pt1(u1, ll1, w1);
	  Point pt2(u1, ur1, w1);
	  mesh << pt1 << " " << pt2 << std::endl;
	}
      for (kj=0, v1=ll1; kj<=num2; ++kj, v1+=del2)
	{
	  Point pt1(ll0, v1, w1);
	  Point pt2(ur0, v1, w1);
	  mesh << pt1 << " " << pt2 << std::endl;
	}
    }

  std::ofstream uncv("unresolved.g2");
  std::ofstream uncv2("shortcurves.g2");

  if (lastcell < firstcell)
    lastcell = num1*num2*num3;
  int ncell = 0;
  for (kh=0, w1=ll2, w2=w1+del3; kh<num3; ++kh, w1=w2, w2+=del3)
    for (kj=0, v1=ll1, v2=v1+del2; kj<num2; ++kj, v1=v2, v2+=del2)
      for (ki=0, u1=ll0, u2=u1+del1; ki<num1; ++ki, u1=u2, u2+=del1, ++ncell)
	{
	  if (ncell < firstcell)
	    continue;
	  if (ncell >= lastcell)
	    break;
	      
	  Point ll(u1, v1, w1);
	  Point ur(u2, v2, w2);
	  int coinc = 0;
	  int stat = 0;
	  try {
	    stat = quad.cellStat(ll, ur, coinc);
	  }
	  catch (...)
	    {
	      std::cout << "Cell failed" << std::endl;
	    }
	std::cout << "Cell nr " << (kh*num2+kj)*num1+ki+1 <<", status: " << stat << std::endl;

	if (stat> 0)
	  {
	    vector<vector<shared_ptr<ParamSurface> > > unresolved_cells;
	    vector<vector<shared_ptr<ParamSurface> > > small_sfs;
	    vector<double> quadpt;
	    vector<double> ptweights;
	    vector<double> bdquad;
	    vector<double> bdnorm;
	    vector<double> bdweights;
	    try {
	      quad.quadrature(ll, ur, quadpt, ptweights, unresolved_cells,
			      bdquad, bdnorm, bdweights, small_sfs, stat, coinc);
  
	    }
	    catch (...)
	      {
		std::cout << "Cell failed" << std::endl;
	      }

	    outfile << "400 1 0 4 100 100 55 255" << std::endl;
	    outfile << ptweights.size() << std::endl;
	    for (size_t kr=0; kr<quadpt.size(); kr+=3)
	      {
		Point pt(quadpt[kr], quadpt[kr+1], quadpt[kr+2]);
		outfile << pt << std::endl;
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
		for (size_t kr=0; kr<bdquad.size(); kr+=3)
		  {
		    Point pt(bdquad[kr], bdquad[kr+1], bdquad[kr+2]);
		    outfile << pt << std::endl;
		  }

		outfile << "410 1 0 4 155 0 100 255" << std::endl;
		outfile << bdweights.size() << std::endl;
		for (size_t kr=0; kr<bdquad.size(); kr+=3)
		  {
		    Point pt1(bdquad[kr], bdquad[kr+1], bdquad[kr+2]);
		    Point pt2(bdquad[kr]+0.05*bdnorm[kr], bdquad[kr+1]+0.05*bdnorm[kr+1],
			      bdquad[kr+2]+0.05*bdnorm[kr+2]);
		    outfile << pt1 << " " << pt2 << std::endl;
		  }

		for (size_t kk=0; kk<small_sfs.size(); ++kk)
		  for (size_t kr=0; kr<small_sfs[kk].size(); ++kr)
		    {
		      small_sfs[kk][kr]->writeStandardHeader(uncv2);
		      small_sfs[kk][kr]->write(uncv2);
		    }
	    
	      }
	  }

	int stop_break = 1;
      }
}

  
