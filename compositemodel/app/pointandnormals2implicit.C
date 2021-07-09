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

#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "sisl.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
  if (argc != 6)
    {
    std::cout << "Input parameters : Input file, offset dist, mu, output file viz, degree"  << std::endl;
    exit(-1);
    }
  
  // Read the points and normals from file
  ifstream input(argv[1]);
  double offdist = atof(argv[2]);
  double mu = atof(argv[3]);
  ofstream outviz(argv[4]);
  int degree = atoi(argv[5]);

  if ( !input)
    {
      std::cerr << "Cannot open file for input\n";
      return 1;
    }

  // Read data
  double val;
  vector<double> points_extended;
  while (input >> val)
    {
      points_extended.push_back(val);
    }

  // Extract point data
  vector<double> points3D(points_extended.size()/5);
  vector<double> offpt3D(points_extended.size()/5);
  size_t kj = 0;
  for (size_t ki=0; ki<points_extended.size(); ki+=15)
    {
      Point norm(points_extended[ki+3], points_extended[ki+4], points_extended[ki+5]);
      norm.normalize();
      for (int ka=0; ka<3; ++ka, ++kj)
	{
	  points3D[kj] = points_extended[ki+ka];
	  offpt3D[kj] = points3D[kj] + offdist*norm[ka];
	}
    }  
  PointCloud3D cloud(&points3D[0], points3D.size()/3);
  PointCloud3D offcloud(&offpt3D[0], offpt3D.size()/3);

  std::ofstream of("IPA_cloud.g2");
  of << "400 1 0 4 0 0 255 255" << std::endl;
  cloud.write(of);
  of << "400 1 0 4 0 255 0 255" << std::endl;
  offcloud.write(of);
  
  // Implicitize
  ImplicitizePointCloudAlgo implicitize(cloud, degree);
  implicitize.perform();
  
  // Get result
  BernsteinTetrahedralPoly implicit;
  BaryCoordSystem3D bc;
  double sigma_min;
  implicitize.getResultData(implicit, bc, sigma_min);

  // Check accuracy
  double avdist = 0.0;
  double avdist2 = 0.0;
  double maxdist = 0.0;
  int dim = 3;
  double avdist_2 = 0.0;
  double avdist2_2 = 0.0;
  double maxdist_2 = 0.0;
  int numpt = cloud.numPoints();
  int nb = (degree+1)*(degree+2)*(degree+3)/6;
  for (int ki=0; ki<numpt; ++ki)
    {
      Vector3D curr = cloud.point(ki);
      Vector4D bary = bc.cartToBary(curr);
      double dist = implicit(bary);

      double dist2 = 0.0;
      int n = nb - 1;
      double bsum = 0.0;
      for (int i = 0; i <= degree; ++i)
	{
	  for (int j = 0; j <= degree-i; ++j) {
	    for (int k = 0; k <= degree-i-j; ++k, --n) {
	      {
		int l = degree-i-j-k;
		double bval = implicit.evalBasis(bary, i, j, k, l);
		bsum += bval;
		double coef = implicit[n];
		dist2 += bval*coef;
	      }
	    }
	  }
	}
      maxdist = std::max(maxdist, fabs(dist));
      avdist += dist;
      avdist2 += fabs(dist);
      maxdist_2 = std::max(maxdist, fabs(dist2));
      avdist_2 += dist2;
      avdist2_2 += fabs(dist2);
    }
  avdist /= (double)numpt;
  avdist2 /= (double)numpt;
  avdist_2 /= (double)numpt;
  avdist2_2 /= (double)numpt;
  std::cout << "Maximum distance: " << maxdist << std::endl;
  std::cout << "Average distance: " << avdist << std::endl;
  std::cout << "Average distance, absolute values: " << avdist2 << std::endl;
  std::cout << "Maximum distance: " << maxdist_2 << std::endl;
  std::cout << "Average distance: " << avdist_2 << std::endl;
  std::cout << "Average distance, absolute values: " << avdist2_2 << std::endl;

  // IPA
  vector<double> cfs(nb, 0.0);
  BernsteinTetrahedralPoly impl2(degree, cfs);
  //BernsteinTetrahedralPoly impl2 = implicit;

  int maxiter =  100000;
  vector<double> delta(2*numpt);
  vector<double> basisval(2*nb*numpt);

  for (int ki=0; ki<2*numpt; ++ki)
    {
      int n = nb - 1;
      Vector3D curr = (ki<numpt) ? cloud.point(ki) : offcloud.point(numpt-ki);
      Vector4D bary = bc.cartToBary(curr);
      for (int i = 0; i <= degree; ++i)
	{
	  for (int j = 0; j <= degree-i; ++j) {
	    for (int k = 0; k <= degree-i-j; ++k, --n) {
	      {
		int l = degree-i-j-k;
		double bval = impl2.evalBasis(bary, i, j, k, l);
		basisval[ki*nb+n] = bval;
	      }
	    }
	  }
	}
    }
  
  for (int iter=0; iter<maxiter; ++iter)
    {
      // Compute residuals
      for (int ki=0; ki<numpt; ++ki)
	{
	  Vector3D curr = cloud.point(ki);
	  Vector4D bary = bc.cartToBary(curr);
	  double dist = impl2(bary);
	  delta[ki] = -dist;
	}
      for (int ki=0; ki<numpt; ++ki)
	{
	  Vector3D curr = offcloud.point(ki);
	  Vector4D bary = bc.cartToBary(curr);
	  double dist = impl2(bary);
	  delta[numpt+ki] = 0.000001*(offdist - dist);
	}
      
      int n = nb - 1;
      vector<double> cfs2(nb);
      for (int i = 0; i <= degree; ++i)
	{
	  for (int j = 0; j <= degree-i; ++j) {
	    for (int k = 0; k <= degree-i-j; ++k, --n) {
	      {
		int l = degree-i-j-k;
		double del = 0.0;
		for (int ki=0; ki<2*numpt; ++ki)
		  {
		    double bval = basisval[ki*nb+n];
		    del += (bval*delta[ki]);
		  }
		double coef = impl2[n];
		cfs2[n] = coef+mu*del;
	      }
	    }
	  }
	}
      BernsteinTetrahedralPoly impl3(degree, cfs2);
      impl2 = impl3;
      int stop_break = 1;
    }
  
  // Check accuracy
   avdist = 0.0;
   maxdist = 0.0;
   avdist_2 = 0.0;
   for (int ki=0; ki<numpt; ++ki)
    {
      Vector3D curr = cloud.point(ki);
      Vector4D bary = bc.cartToBary(curr);
      double dist = impl2(bary);

      maxdist = std::max(maxdist, fabs(dist));
      avdist += dist;
      avdist2 += fabs(dist);
    }
   avdist /= (double)numpt;
   avdist2 /= (double)numpt;
   std::cout << "IPA maximum distance: " << maxdist << std::endl;
   std::cout << "IPA average distance: " << avdist << std::endl;
   std::cout << "IPA average distance, absolute values: " << avdist2 << std::endl;

   BoundingBox ptbox = cloud.boundingBox();
    BoundingBox ptbox2 = offcloud.boundingBox();
    ptbox.addUnionWith(ptbox2);
    Point low = ptbox.low();
    Point high = ptbox.high();
    Point bmid = 0.5*(low + high);

    Point dir(3);
    int nmb_sample;
    std::cout << "Box min: " << low << ", box max: " << high << std::endl;
    std::cout << "Give view direction: " << std::endl;
    std::cin >> dir;
    std::cout << "Give number of points: " << std::endl;
    std::cin >> nmb_sample;

    // Define bounding box as a surface model
    double gap = 1.0e-6;
    Point xdir(1.0, 0.0, 0.0);
    Point ydir(0.0, 1.0, 0.0);
    Point zdir(0.0, 0.0, 1.0);
    CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
    shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low, xdir, ydir, high[0]-low[0],
							  high[1]-low[1], high[2]-low[2]));
    
    // Find the coordinate direction with the largest angle with the view direction
    double a1 = xdir.angle(dir);
    double a2 = ydir.angle(dir);
    double a3 = zdir.angle(dir);
    Point dir2;
    if (a1 > std::min(a2, a3))
      dir2 = xdir;
    else if (a2 > a3)
      dir2 = ydir;
    else
      dir2 = zdir;
    Point dir3 = dir%dir2;
    dir2 = dir%dir3;
    if (dir2*(high-low) < 0.0)
      dir2 *= -1.0;
    if (dir3*(high-low) < 0.0)
      dir3 *= -1.0;
    dir2.normalize();
    dir3.normalize();
    double len = low.dist(high);
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-0.5*len*dir2, 0.0, bmid+0.5*len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-0.5*len*dir3, 0.0, bmid+0.5*len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ka, kb, kc;
    // Evaluate line

    int ik = degree + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (ka=0; ka<ik; ++ka)
      et[ik+ka] = 1.0;

    vector<double> points;
    vector<double> tmpline;
    for (kb=0, p2=0.0; kb<nmb_sample; ++kb, p2+=del)
      {
	for (ka=0, p1=0.0; ka<nmb_sample; ++ka, p1+=del)
	  {
	    // Compute barysentric coordinates of end points of line
	    // First cartesian
	    Point sfpos = ssf->ParamSurface::point(p1,p2);

	    Point cart1 = sfpos + len*dir;
	    Point cart2 = sfpos - len*dir;
	    tmpline.insert(tmpline.end(), cart1.begin(), cart1.end());
	    tmpline.insert(tmpline.end(), cart2.begin(), cart2.end());

	    Vector4D bary1 = bc.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
	    Vector4D bary2 = bc.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

	    Vector3D tp1 = bc.baryToCart(bary1);
	    Vector3D tp2 = bc.baryToCart(bary2);
	    
	    // Pick line
	    BernsteinPoly line = impl2.pickLine(bary1, bary2);

	    // Compute zeroes of bernstein polynomial
	    // First make sisl curve
	    vector<double> ecoef(line.coefsBegin(), line.coefsEnd());
	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
	    double zero = 0.0;

	    // Intersect
	    double eps = 1.0e-6;
	    int kstat = 0;
	    int kcrv=0, kpt=0;
	    double *epar = 0;
	    SISLIntcurve **intcv = 0;
	    if (qc)
	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);

	    // Compute cartesian points and curves associated with intersections
	    for (kc=0; kc<kpt; ++kc)
	      {
		Vector4D barypt = (1.0 - epar[kc])*bary1 + epar[kc]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		Vector3D pos = bc.baryToCart(barypt);
		points.insert(points.end(), pos.begin(), pos.end());
	      }
	  }
      }

  // Output
    if (points.size() > 0)
      {
	PointCloud3D ptcloud(&points[0], points.size()/3);
	outviz << "400 1 0 4 255 0 0 255" << std::endl;
	ptcloud.write(outviz);
      }
}
