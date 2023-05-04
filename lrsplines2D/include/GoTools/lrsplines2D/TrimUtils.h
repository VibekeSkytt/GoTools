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

#ifndef TRIMUTILS_H
#define TRIMUTILS_H

#include "GoTools/geometry/BoundedSurface.h"
#include <vector>
#include <math.h>

namespace Go
{
  /// Utility functionality related to trimming of LR B-spline surfaces with
  /// respect to a corresponding point cloud
  class TrimUtils
  {
  public:
    /// Constuctor given a parameterized point cloud. The points are given as (u,v,x,y,z)
    /// or (x,y,z) if the z-coordinate is parameterized by its x- and y-values.
    /// Note that the sequence of the points will be changed
    TrimUtils(double* points, int nmb_pts, int dim);

    /// Constuctor given a parameterized point cloud. The points are given as (u,v,x,y,z)
    /// or (x,y,z) if the z-coordinate is parameterized by its x- and y-values.
    /// Note that the sequence of the points will be changed.
    /// domain - parameter domain of surface to be trimmed. Must be equal to or larger
    /// than the bounding domain of the parameter points
    TrimUtils(double* points, int nmb_pts, int dim, double domain[]);

    /// Constuctor given a number of parameterized point cloud. The points are given as 
    /// (u,v,x,y,z) or (x,y,z) if the z-coordinate is parameterized by its x- and y-values.
    /// Note that the sequence of the points will be changed.
    /// domain - parameter domain of surface to be trimmed. Must be equal to or larger
    /// than the bounding domain of the parameter points
    /// The parameter all indicates if all point clouds need to have points if an area
    /// is regarded as being inside the trimming loop
    TrimUtils(std::vector<double*>& points, std::vector<int>& nmb_pts, int dim,
	      double domain[], bool all);

    ~TrimUtils();

    /// Using recusion, compute a polygon in the parameter domain surrounding the point 
    /// cloud and one polygon for each hole in the point cloud if outer_only = false.
    /// \param max_level Maximum recursion level (typical 1-3)
    /// \param nmb_div Number of boxes in which the domain will be divided into at
    /// each recursion level (typically a number between 8 and 20)
    /// \param seqs Sets of equences of parameter points defining one polygon for the 
    /// outer boundary and one for each hole. The outer boundary will be represented 
    /// by the first sequence
    void computeTrimSeqs(int max_level, int nmb_div,
			 std::vector<std::vector<double> >& seqs,
			 bool outer_only = true);

    /// Using recusion, compute one or more polygons in the parameter domain surrounding
    // isolated parts of the the point 
    /// cloud and one polygon for each hole in the point cloud.
    /// \param max_level Maximum recursion level (typical 1-3)
    /// \param nmb_div Number of boxes in which the domain will be divided into at
    /// each recursion level (typically a number between 8 and 20)
    /// \param seqs One set of sequences of parameter points for each isolated sub cloud.
    /// The sequences define one polygon for the outer boundary
    /// and one for each hole. For each set, the outer boundary will be represented by 
    /// the first sequence
    void computeAllTrimSeqs(int max_level, int nmb_div,
			    std::vector<std::vector<std::vector<double> > >& seqs);

    double getDomainDiag()
    {
      return sqrt(del_u_*del_u_ + del_v_*del_v_);
    }

    void getDomainLengths(double& del_u, double& del_v)
    {
      del_u = del_u_;
      del_v = del_v_;
    }

    /// Check for points left outside of the trimming loop.
    void extractOutsidePoint(shared_ptr<BoundedSurface>& surf,
			     std::vector<std::vector<double> >& outpoints,
			     int max_nmb);

  private:
    struct SubCloud
    {
      SubCloud()
      {
	processed_ = false;
	limit_[0] = limit_[1] = limit_[2] = limit_[3] = 0;
	num_with_pts_ = 0;
      }

      void setSize(size_t num)
      {
	nmb_pts_.resize(num);
	ix1_.resize(num);
	ix2_.resize(num);
	bb_.resize(num);
      }
      
      void setInfo(std::vector<int>& nmb_pts, 
		   std::vector<int>& ix1, std::vector<int>& ix2,
		   double dom[], std::vector<std::vector<double> >& bb)
      {
	nmb_pts_ = nmb_pts;
	ix1_ = ix1;
	ix2_ = ix2;
	bb_ = bb;
	for (int ki=0; ki<4; ++ki)
	  {
	    dom_[ki] = dom[ki];
	  }
	num_with_pts_ = 0;
	for (size_t ki=0; ki<nmb_pts_.size(); ++ki)
	  if (nmb_pts_[ki] > 0)
	    num_with_pts_++;
      }

      void setInfo(int idx, int nmb_pts, int ix1, int ix2,
		   double dom[], std::vector<double>& bb)
      {
	if (idx >= 0 && idx < (int)nmb_pts_.size())
	  {
	    if (nmb_pts_[idx] > 0)
	      num_with_pts_--;
	    if (nmb_pts > 0)
	      num_with_pts_++;
	    nmb_pts_[idx] = nmb_pts;
	    ix1_[idx] = ix1;
	    ix2_[idx] = ix2;
	    bb_[idx] = bb;
	    for (int ki=0; ki<4; ++ki)
	      {
		dom_[ki] = dom[ki];
	      }
	  }
      }

      void setDomain(double dom[])
      {
	for (int ki=0; ki<4; ++ki)
	  {
	    dom_[ki] = dom[ki];
	  }
      }

      void leftBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[0]);
	bd.push_back(dom_[3]);
	bd.push_back(dom_[0]);
	bd.push_back(dom_[2]);
      }
    
      void rightBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[1]);
	bd.push_back(dom_[2]);
	bd.push_back(dom_[1]);
	bd.push_back(dom_[3]);
      }

      void lowerBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[0]);
	bd.push_back(dom_[2]);
	bd.push_back(dom_[1]);
	bd.push_back(dom_[2]);
      }

      void upperBd(std::vector<double>& bd)
      {
	bd.push_back(dom_[1]);
	bd.push_back(dom_[3]);
	bd.push_back(dom_[0]);
	bd.push_back(dom_[3]);
      }

      bool limitedSupport(double frac)
      {
	double bb0 = bb_[0][0];
	double bb1 = bb_[0][1];
	double bb2 = bb_[0][2];
	double bb3 = bb_[0][3];
	for (size_t kj=1; kj<bb_.size(); ++kj)
	  {
	    bb0 = std::min(bb0, bb_[kj][0]);
	    bb1 = std::max(bb1, bb_[kj][1]);
	    bb2 = std::min(bb2, bb_[kj][2]);
	    bb3 = std::max(bb3, bb_[kj][3]);
	  }
	if ((bb1-bb0)/(dom_[1]-dom_[0]) < frac)
	  return true;
	if ((bb3-bb2)/(dom_[3]-dom_[2]) < frac)
	  return true;
	return false;
      }

      double dom_[4];
      std::vector<std::vector<double> > bb_;
      std::vector<int> nmb_pts_;
      std::vector<int> ix1_;
      std::vector<int> ix2_;
      bool processed_;
      short limit_[4];   // Sequence: left, right, lower, upper
      int num_with_pts_;
    };

    double eps2_;   // Equality tolerance
    std::vector<double*> points_;
    std::vector<int> nmb_pts_;
    int dim_;
    double domain_[4];
    double del_u_;
    double del_v_;
    bool all_;

    void computeTrimInfo(SubCloud& cloud,
			 int max_level, int nmb_div,
			 std::vector<std::vector<double> >& seqs);

    void distributePointCloud(std::vector<int>& ix1, std::vector<int>& ix2,
			      double domain[4],
			      int nmb_u, int nmb_v,
			      std::vector<SubCloud>& sub_clouds);

    void setSubSeq(SubCloud& cloud,
		   std::vector<double>& seq);
    
    void removeFalseTrimSeqs(std::vector<double>& limitseqs,
			     std::vector<std::vector<double> >& seqs);

    void mergeTrimSeqs(std::vector<std::vector<double> >& seqs,
		       std::vector<std::vector<double> >& seqs2);

    void cleanTrimResults(int limitsize,
			  std::vector<std::vector<double> >& seqs);

    void reOrganizeSeqs(std::vector<std::vector<double> >& seqs,
			bool outer_only);

    void reOrganizeSeqs(std::vector<std::vector<double> >& seqs,
			std::vector<int>& outer);

   void extractOutsidePoint(shared_ptr<BoundedSurface>& surf,
			    SubCloud& cloud,
			    std::vector<std::vector<double> >& outpoints,
			    int max_nmb);

   };
};

#endif
