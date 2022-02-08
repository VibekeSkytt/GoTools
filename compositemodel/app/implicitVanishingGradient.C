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
    std::cout << "Input parameters : Input file, output file, output file viz, degree, gradient factor"  << std::endl;
    exit(-1);
  }
  
    // Read the point cloud from file
    ifstream input(argv[1]);
    ofstream output(argv[2]);
    ofstream outviz(argv[3]);
    int degree = atoi(argv[4]);
    double gradfac = atof(argv[5]);
    ObjectHeader header;
    PointCloud3D cloud;
    input >> header >> cloud;

    BoundingBox ptbox = cloud.boundingBox();
    Point low = ptbox.low();
    Point high = ptbox.high();
    Point diag = high - low;
    double diaglen = diag.length();

    // Implicitize
    ImplicitizePointCloudAlgo implicitize(cloud, degree);
    implicitize.perform();

    // Get result
    BernsteinTetrahedralPoly implicit;
    BaryCoordSystem3D bc;
    double sigma_min;
    implicitize.getResultData(implicit, bc, sigma_min);

    // Differentiate
    Vector3D origo(0.0, 0.0, 0.0);
    Vector3D x1(1.0, 0.0, 0.0);
    Vector3D y1(0.0, 1.0, 0.0);
    Vector3D z1(0.0, 0.0, 1.0);
    Vector4D borigo = bc.cartToBary(origo);
    Vector4D bx1 = bc.cartToBary(x1);
    Vector4D by1 = bc.cartToBary(y1);
    Vector4D bz1 = bc.cartToBary(z1);
    Vector4D bvecx = bx1 - borigo;
    Vector4D bvecy = by1 - borigo;
    Vector4D bvecz = bz1 - borigo;
    BernsteinTetrahedralPoly bderx, bdery, bderz;
    implicit.deriv(1, bvecx, bderx);
    implicit.deriv(1, bvecy, bdery);
    implicit.deriv(1, bvecz, bderz);
    
    int ik = degree + 1;
    vector<double> et(2*ik, 0.0);  // Knot vector of line curve
    for (int ki=0; ki<ik; ++ki)
      et[ik+ki] = 1.0;

    // Check accuracy
    double avdist = 0.0;
    double maxdist = 0.0;
    double avdist0 = 0.0;
    double maxdist0 = 0.0;
    double avdist1 = 0.0;
    double maxdist1 = 0.0;
    int nmb0 = 0;
    int numpt = cloud.numPoints();
    vector<Vector3D> ptvecs;
    vector<double> projpts;
    vector<double> gradval1(numpt, -1.0);
    vector<double> gradval2;
    double lfac = 0.01;
    double minlen1 = std::numeric_limits<double>::max();
    double maxlen1 = 0.0;
    double minlen2 = std::numeric_limits<double>::max();
    double maxlen2 = 0.0;
    for (int ki=0; ki<numpt; ++ki)
      {
	Vector3D curr = cloud.point(ki);
	Vector4D bary = bc.cartToBary(curr);
	double dist = implicit(bary);
	double fdist, bdist;
	maxdist = std::max(maxdist, fabs(dist));
	avdist += fabs(dist);
	double dx1 = bderx(bary);
	double dy1 = bdery(bary);
	double dz1 = bderz(bary);
	Vector3D grad(dx1,dy1,dz1);
	double gradlen = grad.length();
	gradval1[ki] = gradlen;
	minlen1 = std::min(minlen1, gradlen);
	maxlen1 = std::max(maxlen1, gradlen);
	Vector3D curr2 = curr + grad;
	Vector4D bary2 = bc.cartToBary(curr2);
	Vector4D dv = bary2 - bary;
	bdist = implicit(bary2);
	Vector4D dir = (dv.length() < 1.0e-15) ? dv : dv/dv.length();
	int sgn = (dist > 0) ? -1 : 1;
	dir *= sgn*lfac*diaglen;
	Vector4D baryx = bary + dir;;
	double distx = implicit(baryx);

	BernsteinPoly line = implicit.pickLine(bary, baryx);
	vector<double> ecoef(line.coefsBegin(), line.coefsEnd());

	double tpar;
	int kstat = 0;
	double dist0 = 1.0e8;
	int found = 1;
	if (std::max(fabs(dist), fabs(distx)) > 1.0e-15)
	  {
	    // Scale up the coefficients to avoid too small number in the iteration
	    double sfac = std::max(100.0/std::max(dist, distx), 1.0);
	    for (size_t kr=0; kr<ecoef.size(); ++kr)
	      ecoef[kr] *= sfac;
	    SISLCurve *qc = newCurve(ik, ik, &et[0], &ecoef[0], 1, 1, 1);
	    double zero = 0.0;
	    SISLPoint *zpt = newPoint(&zero, 1, 0);

	    // Intersect
	    double eps = 1.0e-6;
	    int kcrv=0, kpt=0;
	    double *epar = 0;
	    SISLIntcurve **intcv = 0;
	    if (qc)
	      s1871(qc, &zero, 1, eps, &kpt, &epar, &kcrv, &intcv, &kstat);
	    if (kpt == 0)
	      found = 0;
	    if (kpt > 1)
	      std::sort(&epar[0], &epar[kpt]);
	    if (kpt > 0)
	      {
		Vector4D barypt = (1.0 - epar[0])*bary + epar[0]*baryx;
		Vector3D cartpt = bc.baryToCart(barypt);
		double ptdist = curr.dist(cartpt);
		tpar = epar[0];
	      }
	    else
	      {
		tpar = (fabs(dist) < fabs(distx)) ? et[ik-1] : et[ik];
		found = 0;
	      }
	    if (qc) freeCurve(qc);
	    if (intcv) freeIntcrvlist(intcv, kcrv);
	    if (epar) free(epar);
	  }
	else
	  {
	    tpar = (fabs(dist) < fabs(distx)) ? et[ik-1] : et[ik];
	  }
	
	Vector3D proj;
	if (found)
	  {
	    Vector4D barypt = (1.0 - tpar)*bary + tpar*baryx;
	    fdist = implicit(barypt);
	    Vector3D cartpt = bc.baryToCart(barypt);
	    dist0 = curr.dist(cartpt);
	    maxdist0 = std::max(maxdist0, dist0);
	    avdist0 += dist0;
	    ++nmb0;
	    projpts.insert(projpts.end(), &cartpt[0], &cartpt[3]);
	    double dx12 = bderx(barypt);
	    double dy12 = bdery(barypt);
	    double dz12 = bderz(barypt);
	    Vector3D grad2(dx12,dy12,dz12);
	    double gradlen2 = grad2.length();
	    gradval2.push_back(gradlen2);
	    minlen2 = std::min(minlen2, gradlen2);
	    maxlen2 = std::max(maxlen2, gradlen2);

	  }
	double dist2 = implicit(bary2);
	Vector3D currpt2 = bc.baryToCart(bary2);
	Vector3D norm = currpt2 - curr;
	double len = norm.length();
	if (len > 1.0e-10)
	  norm /= len;
	ptvecs.push_back(norm);
	double dist1 = fabs(dist)/gradlen;
	maxdist1 = std::max(maxdist1, dist1);
	avdist1 += dist1;
      }
    avdist /= (double)numpt;
    avdist0 /= (double)nmb0;
    avdist1 /= (double)numpt;
    std::cout << "Maximum field: " << maxdist << std::endl;
    std::cout << "Average field: " << avdist << std::endl;
    std::cout << "Maximum projected distance: " << maxdist0 << std::endl;
    std::cout << "Average projected distance: " << avdist0 << std::endl;
    std::cout << "numpt-nmb0: " << numpt - nmb0 << std::endl;
    std::cout << "Maximum estimated distance: " << maxdist1 << std::endl;
    std::cout << "Average esitmated distance: " << avdist1 << std::endl;
    std::cout << "Minimum gradient length point: " << minlen1 << std::endl;
    std::cout << "Maximum gradient length point: " << maxlen1 << std::endl;
    std::cout << "Minimum gradient length projection: " << minlen2 << std::endl;
    std::cout << "Maximum gradient length projection: " << maxlen2 << std::endl;

    // Write out implicit function
    std::cout << "Sigma_min: " << sigma_min << std::endl;
    output << implicit << endl
	   << bc << endl;
    cout << "Data written to output" << endl;

    vector<double> smallgrad1, smallgrad2;
    double gradlim1 = minlen1 + gradfac*(maxlen1 - minlen1);
    double gradlim2 = minlen2 + gradfac*(maxlen2 - minlen2);
    for (int ka=0; ka<numpt; ++ka)
      {
	if (gradval1[ka] >= 0.0 && gradval1[ka] <= gradlim1)
	  {
	    Vector3D pos = cloud.point(ka);
	    smallgrad1.insert(smallgrad1.end(), &pos[0], &pos[3]);
	  }
      }
    for (size_t kr=0; kr<projpts.size()/3; ++kr)
      {
	if (gradval2[kr] >= 0.0 && gradval2[kr] <= gradlim2)
	  {
	    smallgrad2.insert(smallgrad2.end(), &projpts[3*kr], &projpts[3*kr+3]);
	  }
      }
   
    // Fetch points on the implicit surface
    Point dir(3);
    int nmb_sample;
    std::cout << "Box min: " << low << ", box max: " << high << std::endl;
    std::cout << "Give view direction: " << std::endl;
    std::cin >> dir;
    std::cout << "Give number of points: " << std::endl;
    std::cin >> nmb_sample;

    // Find extension of view area
    dir.normalize();
    Point bmid = 0.5*(low + high);

    // Define bounding box as a surface model
    double gap = 1.0e-6;
    Point xdir(1.0, 0.0, 0.0);
    Point ydir(0.0, 1.0, 0.0);
    Point zdir(0.0, 0.0, 1.0);
    CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
    shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-0.5*diag, xdir, ydir, 2*diag[0],
							  2*diag[1], 2*diag[2]));
    
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
    shared_ptr<SplineCurve> cv1(new SplineCurve(bmid-len*dir2, 0.0, bmid+len*dir2, 1.0));
    shared_ptr<SplineCurve> cv2(new SplineCurve(bmid-len*dir3, 0.0, bmid+len*dir3, 1.0));
    SweepSurfaceCreator sweep;
    shared_ptr<SplineSurface> ssf(sweep.linearSweptSurface(*cv1, *cv2, bmid));
    double del = 1.0/(double)(nmb_sample-1);
    double p1, p2;
    int ki, kj, kr;

    vector<double> points;
    vector<double> vecs;
    vector<double> linesegs;
    vector<double> der;
    vector<double> der2;
    vector<double> lineder;
    // Evaluate line
    vector<double> tmpline;
    for (kj=0, p2=0.0; kj<nmb_sample; ++kj, p2+=del)
      {
	for (ki=0, p1=0.0; ki<nmb_sample; ++ki, p1+=del)
	  {
	    // Compute barysentric coordinates of end points of line
	    // First cartesian
	    Point sfpos = ssf->ParamSurface::point(p1,p2);
	    Point cart1 = sfpos + len*dir;
	    Point cart2 = sfpos - len*dir;

	    Vector4D bary1 = bc.cartToBary(Vector3D(cart1[0], cart1[1], cart1[2]));
	    Vector4D bary2 = bc.cartToBary(Vector3D(cart2[0], cart2[1], cart2[2]));

	    Vector3D tp1 = bc.baryToCart(bary1);
	    Vector3D tp2 = bc.baryToCart(bary2);
	    
	    // Pick line
	    BernsteinPoly line = implicit.pickLine(bary1, bary2);

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
	    for (kr=0; kr<kpt; ++kr)
	      {
		Vector4D barypt = (1.0 - epar[kr])*bary1 + epar[kr]*bary2;
		int kb;
		for (kb=0; kb<4; ++kb)
		  if (barypt[kb] < -0.001 || barypt[kb] > 1.001)
		    break;
		if (kb < 4)
		  continue;
		
		Vector3D pos = bc.baryToCart(barypt);
		points.insert(points.end(), pos.begin(), pos.end());

		double dx1 = bderx(barypt);
		double dy1 = bdery(barypt);
		double dz1 = bderz(barypt);
		Vector3D grad(dx1,dy1,dz1);
		vecs.insert(vecs.end(), grad.begin(), grad.end());
		
		// Check
		double dist2 = implicit(barypt);
		if (dist2 > 1.0e-4)
		  std::cout << "dist: " << dist2 << ", ki= " << ki << ", kj= " << kj << std::endl;
	      }
	    for (kr=0; kr<kcrv; ++kr)
	      {
		int ipt = intcv[kr]->ipoint;
		double par1 = intcv[kr]->epar1[0];
		double par2 = intcv[kr]->epar1[ipt-1];
		Point pp1 = (1.0-par1)*cart1 + par1*cart2;
		Point pp2 = (1.0-par2)*cart1 + par2*cart2;
		linesegs.insert(linesegs.end(), pp1.begin(), pp1.end());
		linesegs.insert(linesegs.end(), pp2.begin(), pp2.end());
		  
		// Check
		Point pmid = 0.5*pp1 + 0.5*pp2;
		Vector4D bary = bc.cartToBary(Vector3D(pmid[0], pmid[1], pmid[2]));
		double dist = implicit(bary);
		if (dist > 1.0e-4)
		  std::cout << "line dist: " << dist << ", ki= " << ki << ", kj= " << kj << std::endl;
	      }
	    
	    
	    if (qc) freeCurve(qc);
	    if (intcv) freeIntcrvlist(intcv, kcrv);
	    if (epar) free(epar);
	  }
      }

    // outviz << "410 1 0 4 155 0 100 255" << std::endl;
    // outviz << numpt << std::endl;
    // for (int ki=0; ki<numpt; ++ki)
    //   {
    // 	Vector3D curr = cloud.point(ki);
    // 	outviz << curr << " " << curr+0.2*ptvecs[ki] << std::endl;
    //   }
    
    // Output
    if (points.size() > 0)
      {
	PointCloud3D ptcloud(&points[0], points.size()/3);
	outviz << "400 1 0 4 255 0 0 255" << std::endl;
	ptcloud.write(outviz);
      }

    if (der.size() > 0)
      {
	PointCloud3D ptcloud(&der[0], der.size()/3);
	outviz << "400 1 0 4 0 255 0 255" << std::endl;
	ptcloud.write(outviz);
      }

    if (projpts.size() > 0)
      {
	PointCloud3D projcloud(&projpts[0], projpts.size()/3);
	outviz << "400 1 0 4 0 155 100 255" << std::endl;
	projcloud.write(outviz);
      }

    if (smallgrad1.size() > 0)
      {
	PointCloud3D grad1(&smallgrad1[0], smallgrad1.size()/3);
	outviz << "400 1 0 4 55 50 150 255" << std::endl;
	grad1.write(outviz);
      }

        if (smallgrad2.size() > 0)
      {
	PointCloud3D grad2(&smallgrad2[0], smallgrad2.size()/3);
	outviz << "400 1 0 4 155 100 0 255" << std::endl;
	grad2.write(outviz);
      }
return 0;
}
