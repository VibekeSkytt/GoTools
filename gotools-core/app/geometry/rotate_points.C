#include "GoTools/geometry/FileUtils.h"

#include <iostream>
#include <fstream>
#include <cmath>  /* This is to get access to sin, cos, ... from C++ */

void rotate_point_seq(double *pts1, int num_pts, double plane_normal[3],
		      double *pts2)
{
  /* It is assumed that the input points given in pts1 are represented as
     x1, y1, z1, x2, y2, z2, ...., and that the points lie in a plane.
     The number of points is num_pts. plane_normal specifies the normal
     of the plane in which the points are to be rotated. The resulting points
     are stored in pts2. This array must be pre-allocated. The lengths of
     pts1 and pts2 are 3*num_pts  */

  double midp[3];  /* Point in the plane where the input points lie */
  double norm[3];  /* Plane normal */
  double norm1[3]; /* Plane normal computed from two input points and
		      the found point in the plane. Cross product of the 
		      vectors from the two input points to the point in the
		      plane */
  double vec1[3], vec2[3]; /* Vectors between input points and
			      the found point in the plane */
  double vec3[3];  /* Help vector in computation of rotated points */
  int ki, kj, kr, kh;  /* Counters */
  double num_fac = 1.0/(double)num_pts;  /* Factor in computing average from
					    all input points */
  int num4 = num_pts/4;     /* Number of points between the two points
			       in the computation of one plane normal */
  int num10 = num_pts/10;   /* Every tenth point is used to compute plane
			       normal */
  int numd = 3*num_pts;                  /* Size of input data array */
  int nnorm = 0;                         /* Number of samples in plane normal
					    computation */
  double len1 = 0.0, len2 = 0.0;  /* To normalize vectors */
  double nom = 0.0;               /* Nominator in angle computation */
  double angle;                   /* Angle between the two plane normals */
  double eps = 1.0e-9;            /* Tolerance to check if the input and 
				     output plane are the same */
  double cosang, sinang;   /* Cosinus and sinus of angle between planes */
  double rotmat[9];        /* Rotational matrix */
  double tmp[3];           /* Help vector in rotation */
  double tmp2;             /* Used in the computation of rotated points */
  double scpr;             /* Scalar product of computed and output plane 
			      normal. If negative, the computed normal
			      shoult be negated to avoid turning the point
			      sequence */

  /* Compute point in the plane containing the given points */
  midp[0] = midp[1] = midp[2] = 0.0;
  for (ki=0; ki<numd; ki+=3)
    for (kj=0; kj<3; ++kj)
      midp[kj] += num_fac*pts1[ki+kj];

  /* Compute plane axis */
  /* A large number of points is assumed. If the number is less that
     ~100, num10 should be redused */
  norm[0] = norm[1] = norm[2] = 0.0;
  nnorm = 0;
  for (ki=0; ki<numd; ki+=(3*num10))
    {
      kj = (ki+3*num4)%numd;
      for (kr=0; kr<3; ++kr)
	{
	  vec1[kr] = pts1[ki+kr] - midp[kr];
	  vec2[kr] = pts1[kj+kr] - midp[kr];
	}
      norm1[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
      norm1[1] = -vec1[0]*vec2[2] + vec1[2]*vec2[0];
      norm1[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

      len1 = 0.0;
      for (kr=0; kr<3; ++kr)
	len1 += (norm1[kr]*norm1[kr]);
      len1 = sqrt(len1);

      if (len1 > eps)
	{
	  for (kr=0; kr<3; ++kr)
	    norm[kr] += (norm1[kr]/len1);
	  ++nnorm;
	}
    }

  if (nnorm == 0)
    return;  /* Input points lie at a line (or bad choice of sample points).
		Computing the covariance matrix and finding the nullspace is
		more robust, but more involved. */

  /* Normalize axis */
  for (kr=0; kr<3; ++kr)
    norm[kr] /= (double)nnorm;

  /* Ensure that the point cloud is rotated the shortest path to the 
     new plane */
  scpr = 0.0;
  for (kr=0; kr<3; ++kr)
    scpr += plane_normal[kr]*norm[kr];
  if (scpr < 0.0)
    {
      for (kr=0; kr<3; ++kr)
	norm[kr] *= -1.0;
    }
  
  /* The angle between the input and output plane */
  len1 = len2 = nom = 0.0;
  for (kr=0; kr<3; ++kr)
    {
      len1 += (norm[kr]*norm[kr]);
      len2 += (plane_normal[kr]*plane_normal[kr]);
      nom += (norm[kr]*plane_normal[kr]);
    }
  len1 = sqrt(len1);
  len2 = sqrt(len2);
  if (len1 < eps || len2 < eps)
    return;

  angle = acos(nom/(len1*len2));
  if (angle < eps)
    {
      /* Input and output plane is the same. Copy points */
      for (ki=0; ki<numd; ++ki)
	pts2[ki] = pts1[ki];
      return;
    }
	       
  /* Define rotation matrix */
  vec3[0] = norm[1]*plane_normal[2] - norm[2]*plane_normal[1];
  vec3[1] = -norm[0]*plane_normal[2] + norm[2]*plane_normal[0];
  vec3[2] = norm[0]*plane_normal[1] - norm[1]*plane_normal[0];
  len1 = 0.0;
  for (kr=0; kr<3; ++kr)
    len1 += (vec3[kr]*vec3[kr]);
  len1 = sqrt(len1);
  for (kr=0; kr<3; ++kr)
    vec3[kr] /= len1;
  
  cosang = cos(angle);
  sinang = sin(angle);
  rotmat[0] = vec3[0]*vec3[0] + cosang*(1.0-vec3[0]*vec3[0]);
  rotmat[1] = vec3[0]*vec3[1]*(1.0-cosang) - vec3[2]*sinang;
  rotmat[2] = vec3[0]*vec3[2]*(1.0-cosang) + vec3[1]*sinang;
  rotmat[3] = vec3[0]*vec3[1]*(1.0-cosang) + vec3[2]*sinang;
  rotmat[4] = vec3[1]*vec3[1] + cosang*(1.0-vec3[1]*vec3[1]);
  rotmat[5] = vec3[1]*vec3[2]*(1.0-cosang) - vec3[0]*sinang;
  rotmat[6] = vec3[0]*vec3[2]*(1.0-cosang) - vec3[1]*sinang;
  rotmat[7] = vec3[1]*vec3[2]*(1.0-cosang) + vec3[0]*sinang;
  rotmat[8] = vec3[2]*vec3[2] + cosang*(1.0-vec3[2]*vec3[2]);
	
  /* Rotate */
  for (ki=0; ki<numd; ki+=3)
    {
      for (kr=0; kr<3; ++kr)
	tmp[kr] = pts1[ki+kr] - midp[kr];
      for (kj=0; kj<3; ++kj)
	{
	  tmp2 = 0;
	  for (kr=0; kr<3; ++kr)
	    tmp2 += rotmat[3*kj+kr]*tmp[kr];
	  pts2[ki+kj] = tmp2 + midp[kj]; 
	}
    }
}

int main(int argc, char* argv[])
{
  if (argc != 6)
    exit(0);
  std::ifstream input(argv[1]);
  double norm1[3];
  norm1[0] = atof(argv[2]);
  norm1[1] = atof(argv[3]);
  norm1[2] = atof(argv[4]);
  std::ofstream output(argv[5]);

  std::vector<double> data;
  std::vector<double> extent(6);
  int nmb_pts;  /* Number of points, i.e. the length of the data array
		   divided by 3 (del) */
  int del = 3;

  /* Read from file. Uses functionality from SINTEF's GoTools library */
  FileUtils::readTxtPointFile(input, del, data, nmb_pts, extent);
  std::vector<double> data2(data.size());

  /* Call c-funtion to rotate the point sequence to the specified plane */
  rotate_point_seq(&data[0], nmb_pts, norm1, &data2[0]);

  /* Write to file, first the input points so the rotated points. The
     format is native to GoTools */
  output << "400 1 0 4 0 0 255 255" << std::endl;
  output << nmb_pts << std::endl;
  for (int ka=0; ka<nmb_pts; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	output << data[3*ka+kb] << " ";
      output << std::endl;
    }
  
  output << "400 1 0 4 255 0 0 255" << std::endl;
  output << nmb_pts << std::endl;
  for (int ka=0; ka<nmb_pts; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	output << data2[3*ka+kb] << " ";
      output << std::endl;
    }
  int stop_break = 1;
}

