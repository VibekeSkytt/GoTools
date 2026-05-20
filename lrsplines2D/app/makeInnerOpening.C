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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 7 && argc != 8) {
    std::cout << "Input parameters : Input surface file(.g2), output surface (.g2), parameter direction(0, 1), displacement vector (x,y,z), (extra refinement (0/1) 1=default "  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream infile(argv[1]);
  std::ofstream outfile(argv[2]);
  int dir = atoi(argv[3]);
  Point vec(3);
  vec[0] = atof(argv[4]);
  vec[1] = atof(argv[5]);
  vec[2] = atof(argv[6]);
  bool refine = true;
  if (argc == 8)
    refine = atoi(argv[7]);

  // Read input surface
  ObjectHeader header;
  header.read(infile);

  shared_ptr<LRSplineSurface> lrsf(new LRSplineSurface());
  lrsf->read(infile);

  double eps = 1.0e-8;
  
  // Refine
  double umin = lrsf->startparam_u();
  double umax = lrsf->endparam_u();
  double vmin = lrsf->startparam_v();
  double vmax = lrsf->endparam_v();

  const Mesh2D mesh = lrsf->mesh();
  vector<double> knots1, knots2;
  for (auto it=mesh.knotsBegin(XFIXED); it!=mesh.knotsEnd(XFIXED); ++it)
    knots1.push_back(*it);
  for (auto it=mesh.knotsBegin(YFIXED); it!=mesh.knotsEnd(YFIXED); ++it)
    knots2.push_back(*it);
  if (refine)
    {
      vector<double> newknots1(knots1.size()-1);
      vector<double> newknots2(knots2.size()-1);
      for (size_t ki=1; ki<knots1.size(); ++ki)
	newknots1[ki-1] = 0.5*(knots1[ki-1]+knots1[ki]);
      for (size_t ki=1; ki<knots2.size(); ++ki)
	newknots2[ki-1] = 0.5*(knots2[ki-1]+knots2[ki]);

      vector<LRSplineSurface::Refinement2D> refs(newknots1.size()+newknots2.size());
      for (size_t ki=0; ki<newknots1.size(); ++ki)
	refs[ki].setVal(newknots1[ki], vmin, vmax, XFIXED, 1);
      for (size_t ki=0; ki<newknots2.size(); ++ki)
	refs[newknots1.size()+ki].setVal(newknots2[ki], umin, umax, YFIXED, 1);

      for (size_t ki=0; ki<refs.size(); ++ki)
	lrsf->refine(refs[ki]);
      knots1.insert(knots1.end(), newknots1.begin(), newknots1.end());
      std::sort(knots1.begin(), knots1.end());
      knots2.insert(knots2.end(), newknots2.begin(), newknots2.end());
      std::sort(knots2.begin(), knots2.end());
    }

   vector<LRBSpline2D*> bspl1, bspl2;
   double val;
   if (dir == 0)
     {
       // Multiple knot in the inner
       val = knots2[knots2.size()/2];
       double order = lrsf->degree(XFIXED) + 1;
       int ix1 = (int)knots1.size()/2 - order/2;
       int ix2 = ix1 + order;
       //std::cout << order << " " << ix1 << " " << ix2 << std::endl;

       LRSplineSurface::Refinement2D refs2;
       refs2.setVal(val, knots1[ix1], knots1[ix2], YFIXED, order-1);
       lrsf->refine(refs2);

       // Identify coefficients at opening
       lrsf->getBSplinesAtOpening(YFIXED, val, bspl1, bspl2);
     }
   else
     {
        // Multiple knot in the inner
       val = knots1[knots1.size()/2];
       double order = lrsf->degree(YFIXED) + 1;
       int ix1 = (int)knots2.size()/2 - order/2;
       int ix2 = ix1 + order;
       //std::cout << order << " " << ix1 << " " << ix2 << std::endl;

       LRSplineSurface::Refinement2D refs2;
       refs2.setVal(val, knots2[ix1], knots2[ix2], XFIXED, order-1);
       lrsf->refine(refs2);

       // Identify coefficients at opening
       lrsf->getBSplinesAtOpening(XFIXED, val, bspl1, bspl2);
     }
   std::ofstream of("opening_coefs.g2");
   if (bspl1.size() > 0)
     {
       of << "400 1 0 4 255 0 0 255" << std::endl;
       of << bspl1.size() << std::endl;
       for (size_t ki=0; ki<bspl1.size(); ++ki)
	 of << bspl1[ki]->Coef() << std::endl;
     }
   
   if (bspl2.size() > 0)
     {
       of << "400 1 0 4 0 255 0 255" << std::endl;
       of << bspl2.size() << std::endl;
       for (size_t ki=0; ki<bspl2.size(); ++ki)
	 of << bspl2[ki]->Coef() << std::endl;
     }

   vec *= 0.5;
   for (size_t ki=0; ki<bspl1.size(); ++ki)
     {
       Point coef = bspl1[ki]->Coef();
       coef -= vec;
       lrsf->setCoef(coef, bspl1[ki]);
     }
    for (size_t ki=0; ki<bspl2.size(); ++ki)
     {
       Point coef = bspl2[ki]->Coef();
       coef += vec;
       lrsf->setCoef(coef, bspl2[ki]);
     }
      
   std::ofstream of2("opened_coefs.g2");
   if (bspl1.size() > 0)
     {
       of2 << "400 1 0 4 255 0 0 255" << std::endl;
       of2 << bspl1.size() << std::endl;
       for (size_t ki=0; ki<bspl1.size(); ++ki)
	 of2 << bspl1[ki]->Coef() << std::endl;
     }
   
   if (bspl2.size() > 0)
     {
       of2 << "400 1 0 4 0 255 0 255" << std::endl;
       of2 << bspl2.size() << std::endl;
       for (size_t ki=0; ki<bspl2.size(); ++ki)
	 of2 << bspl2[ki]->Coef() << std::endl;
     }

   lrsf->writeStandardHeader(outfile);
   lrsf->write(outfile);

   if (dir == 0)
     {
       shared_ptr<LRSplineSurface> sub1(lrsf->subSurface(umin, vmin, umax, val, eps));
       shared_ptr<LRSplineSurface> sub2(lrsf->subSurface(umin, val, umax, vmax, eps));
   
       sub1->writeStandardHeader(outfile);
       sub1->write(outfile);
       sub2->writeStandardHeader(outfile);
       sub2->write(outfile);
     }
   else
     {
       shared_ptr<LRSplineSurface> sub1(lrsf->subSurface(umin, vmin, val, vmax, eps));
       shared_ptr<LRSplineSurface> sub2(lrsf->subSurface(val, vmin, umax, vmax, eps));
   
       sub1->writeStandardHeader(outfile);
       sub1->write(outfile);
       sub2->writeStandardHeader(outfile);
       sub2->write(outfile);
     }
}


