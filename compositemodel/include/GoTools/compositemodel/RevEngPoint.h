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

#ifndef _REVENGPOINT_H
#define _REVENGPOINT_H

#include "GoTools/compositemodel/ftPointSet.h"
//#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/DirectionCone.h"

namespace Go
{
  class RevEngRegion;
  
  // Enumerations representing classification results
  enum
  {
   PCA_EDGE_UNDEF, PCA_NOT_EDGE, PCA_CLOSE_EDGE, PCA_EDGE
  };

  enum
  {
   C1_EDGE_UNDEF, C1_NOT_EDGE, C1_CLOSE_EDGE, C1_EDGE
  };

  enum
  {
   C2_EDGE_UNDEF, C2_NOT_EDGE, C2_CLOSE_EDGE, C2_EDGE
  };

  enum
  {
   RP_EDGE_UNDEF, RP_NOT_EDGE, RP_CLOSE_EDGE, RP_EDGE
  };

  enum
  {
   PCA_UNDEF, PCA_PLANAR, PCA_LINEAR, PCA_OTHER
  };

  enum
  {
   C1_UNDEF, C1_PEAK, C1_RIDGE, C1_SRIDGE, C1_NONE, C1_FLAT, C1_MINSURF, C1_PIT, C1_VALLEY, C1_SVALLEY 
  };
    
  enum
  {
   SI_UNDEF, SI_PLANE, SI_SCUP, SI_TRO, SI_RUT, SI_SRUT, SI_SAD, SI_SRID, SI_RID, SI_DOM, SI_SCAP
   // Plane, spherical cup, trough, rut, saddle rut, saddle, saddle ridge,
   // ridge, dome, spherical cap
  };

  enum
  {
   PS_UNDEF, PS_PLANE, PS_UC, PS_EC, PS_PC, PS_HC, PS_HS, PS_HX, PS_PX, PS_EX, PS_UX
  };

  /** Subclass of ftSamplePoint; Enhanced point for use in reverse engineering
   */
  class RevEngPoint : public ftSamplePoint
  {

  public:
    /// Constructor
    RevEngPoint();
    
   RevEngPoint(Vector3D xyz, int bnd);

    virtual ~RevEngPoint();

    void computeTriangNormal(double lim);

    void addCovarianceEigen(Point& eigen1, double lambda1, Point& eigen2,
			    double lambda2, Point& eigen3, double lambda3);

    void addMongeInfo(Point& norm, Point& mincvec, double minc, Point& maxcvec,
		      double maxc, double currdist, double avdist, double eps);

    Point getTriangNormal()
    {
      Point zero(0.0,0.0,0.0);
      return normalcone_.greaterThanPi() ? zero :
	      normalcone_.centre();
    }

    double getTriangAngle()
    {
      return normalcone_.greaterThanPi() ? 2.0*M_PI : normalcone_.angle();
    }

    double getMeanEdgLen();

    double getMeanEdgLen(double maxlen);

    const Point& getMongeNormal()
    {
      return Mongenormal_;
    }

    void turnMongeNorm()
      {
	Mongenormal_ *= -1;
      }

    void resetPointAssociation();
    
    const Point& getPCAEigen1()
    {
      return eigen1_;
    }

    const Point& getPCANormal()
    {
      return eigen3_;
    }

    void turnPCA()
    {
      eigen2_ *= -1;
      eigen3_ *= -1;
    }

    Point fetchClosePoints(double radius, int min_nmb, int max_nmb,
			   std::vector<Point>& nearpts);


    void fetchClosePoints2(double radius, int min_nmb, int max_nmb,
			   std::vector<RevEngPoint*>& nearpts,
			   RevEngRegion *region = 0);

    void fetchConnected(RevEngRegion *region, int max_nmb,
			std::vector<RevEngPoint*>& group);

    void setVisited()
    {
      visited_ = 1;
    }

    
    void unsetVisited()
    {
      visited_ = 0;
    }

    bool visited()
    {
      return (visited_ > 0);
    }
    
    void setMoved()
    {
      moved_ = 1;
    }

    
    void unsetMoved()
    {
      moved_ = 0;
    }

    bool moved()
    {
      return (moved_ > 0);
    }

    int nmbEigen()
    {
      return nmb_eigen_;
    }
    
    int nmbMonge()
    {
      return nmb_monge_;
    }
    
    double getSurfaceVariation()
    {
      return sfvariation_;
    }

    double getShapeIndex()
    {
      return shapeindex_;
    }

    double maxPrincipalCurvature()
    {
      return kmax_;
    }

    Point maxCurvatureVec()
    {
      return kvecmax_;
    }

    double minPrincipalCurvature()
    {
      return kmin_;
    }

    Point minCurvatureVec()
    {
      return kvecmin_;
    }

    double GaussCurvature()
    {
      return gausscurv_;
    }

    double meanCurvature()
      {
	return meancurv_;
      }

    double GaussCurvature0()
    {
      return gausscurv0_;
    }

    double meanCurvature0()
      {
	return meancurv0_;
      }

    void setMeanCurvature(double mean)
    {
      meancurv_ = mean;
    }

    void setGaussCurvature(double Gauss)
    {
      gausscurv_ = Gauss;
    }

    void updateCurvature()
    {
      meancurv0_ = meancurv_;
      gausscurv0_ = gausscurv_;
    }

    double getCurvedness()
    {
      return curvedness_;
    }

    void setRp(double rp[2])
    {
      for (int ka=0; ka<2; ++ka)
	rp_[ka] = rp[ka];
    }

    double getRp(int ix)
    {
      return (ix < 0 || ix > 2) ? 0.0 : rp_[ix];
    }

    double getfpa()
    {
      return fpa_;
    }

    double getspa()
    {
      return spa_;
    }
    
    void setEdgeClassification(int c1_edge, int c2_edge, int pca_edge, int rp_edge)
    {
      edge_[0] = pca_edge;
      edge_[1] = c1_edge;
      edge_[2] = c2_edge;
      edge_[3] = rp_edge;
    }

    void setClassification(int ctype, int si_type, int ps_type)
    {
      surf_[0] = PCA_UNDEF;
      surf_[1] = ctype;
      surf_[2] = si_type;
      surf_[3] = ps_type;
    }

    bool isEdge(int edge_class_type)
    {
      if (edge_class_type == 1)
	return (edge_[0] == PCA_EDGE);
      else if  (edge_class_type == 2)
	return (edge_[1] == C1_EDGE);
      else if  (edge_class_type == 3)
	return (edge_[2] == C2_EDGE);
      else if  (edge_class_type == 4)
	return (edge_[3] == RP_EDGE);
      else
	return false;
      // return (edge_[0] == PCA_EDGE || edge_[1] == C1_EDGE ||
      // 	      edge_[2] == C2_EDGE);
    }

    bool closeEdge(int edge_class_type)
    {
      if (edge_class_type == 1)
	return (edge_[0] >= PCA_CLOSE_EDGE);
      else if  (edge_class_type == 2)
	return (edge_[1] >= C1_CLOSE_EDGE);
      else if  (edge_class_type == 3)
	return (edge_[2] >= C2_CLOSE_EDGE);
      else if  (edge_class_type == 4)
	return (edge_[3] >= RP_CLOSE_EDGE);
      else
	return false;
      // return (edge_[0] >= PCA_CLOSE_EDGE || edge_[1] >= C1_CLOSE_EDGE ||
      // 	      edge_[2] >= C2_CLOSE_EDGE);
    }

    bool notEdge(int edge_class_type)
    {
      if (edge_class_type == 1)
	return (edge_[0] <= PCA_NOT_EDGE);
      else if  (edge_class_type == 2)
	return (edge_[1] <= C1_NOT_EDGE);
      else if  (edge_class_type == 3)
	return (edge_[2] <= C2_NOT_EDGE);
      else if  (edge_class_type == 4)
	return (edge_[3] <= RP_NOT_EDGE);
      else
	return true;
      // return (edge_[0] <= PCA_NOT_EDGE && edge_[1] <= C1_NOT_EDGE &&
      // 	      edge_[2] <= C2_NOT_EDGE);
    }

    bool isolatedEdge(int edge_class_type, int nmb, bool close);

    void setEdgeUndef()
    {
      edge_[0] = PCA_EDGE_UNDEF;
      edge_[1] = C1_EDGE_UNDEF;
      edge_[2] = C2_EDGE_UNDEF;
      edge_[3] = RP_EDGE_UNDEF;
    }

    double getPointDistance()
    {
      return ptdist_;
    }

    double getAveragePointDistance()
    {
      return avdist_;
    }

    int surfaceClassification(int classification_type) const;
    
    void adjustWithTriangNorm(double anglim);

    bool hasRegion()
    {
      return (region_ != 0);
    }

    RevEngRegion* region()
    {
      return region_;
    }

    void setRegion(RevEngRegion* region)
    {
      region_ = region;
    }

    void unsetRegion()
    {
      region_ = 0;
    }

    void unsetSurfInfo()
    {
      sfdist_ = -1.0;
      sfang_ = -1.0;
      uv_ = Vector2D(0.0, 0.0);
    }

    void adjacentRegions(std::vector<RevEngRegion*>& adj) const;

    int C1_surf()
    {
      return surf_[1];
    }

    int SI_surf()
    {
      return surf_[2];
    }

    int RP_surf()
    {
      return surf_[3];
    }

    void setGaussRad(double rad)
    {
      Grad_ = rad;
    }

    double getGaussRad()
    {
      return Grad_;
    }

    void setOutlier()
    {
      outlier_ = true;
    }

    void unsetOutlier()
    {
      outlier_ = false;
    }

    bool isOutlier()
    {
      return outlier_;
    }

    int nmbSameClassification(int classification_type) const;

    void setSurfaceDist(double dist, double ang)
    {
      sfdist_ = dist;
      sfang_ = ang;
    }

    double  getSurfaceDist()
    {
      return sfdist_;
    }
    
    void getSurfaceDist(double& dist, double& ang)
    {
      dist = sfdist_;
      ang = sfang_;
    }

    bool isNeighbour(RevEngRegion* reg) const;

    void addMove()
    {
      nmb_move_++;
    }

    int numMove()
    {
      return nmb_move_;
    }

    void getAdjInfo(double mean_edge_len, std::vector<RevEngRegion*>& adj_reg,
		    std::vector<RevEngPoint*>& adj_pt);

    std::vector<RevEngRegion*> getAdjacentRegions() const;
    int numAdjacentRegions() const;
    
    bool mergeWithAdjacent(double mean_edge_len);
    
    void store(std::ostream& os) const;
    void read(std::istream& is, double eps, vector<int>& next_ix);
    
  private:
    double avedglen_;
    
    // How to mark that the properties are not computed? Can check dimension
    // of eigenvectors and principal curvature vectors. Compute derived info
    // once the primary one is in place
    Point eigen1_, eigen2_, eigen3_;  // Eigenvectors of covariance matrix
    double lambda1_, lambda2_, lambda3_; // Eigenvalues of covariance matrix
    // lambda1_>= lambda2_ >= lambda3_
    Point Mongenormal_, kvecmin_, kvecmax_;  // Normal and principal curvature vectors
    double kmin_, kmax_;  // Principal curvatures from Monge's patch
    double ptdist_, avdist_;  // Distance to Monge's patch
    int nmb_eigen_, nmb_monge_;
    DirectionCone normalcone_;  // Span of surface normals computed from triangulation

    double meancurv0_, meancurv_;
    double gausscurv0_, gausscurv_;
    double sfvariation_;  // Surface variation computed from eigenvalues
    double curvedness_;   // Curvedness computed from principal curvatures
    double shapeindex_;   // Computed from principal curvatures
    double rp_[2];        // Smoothness indicator
    double fpa_;          // Flat point possibility association
    double spa_;          // Non-flat point possibility association

    // Parameters corresponding to classification results
    int edge_[4];   // Results of edge classification.
    // Sequence: Surface variation (PCA), curvature (C1), curveness (C2), smoothness indicator (RP)
    int surf_[4];   // Results of surface classification
    // Sequence: PCA, curvature (C1), shape index (SI), point possibility association (PS)
    double Grad_;
    
    // Group (segment) of classified points
    RevEngRegion* region_;
    bool outlier_;

    double sfdist_;  // Distance beteen point and surface representation
    double sfang_;   // Angle between estimated normal and surface normal
    int nmb_move_;

    mutable int visited_;
    mutable int moved_;

    void getNearby(Vector3D xyz, double radius, int max_nmb,
		   std::vector<RevEngPoint*>& near,
		   RevEngRegion *region = 0);
  };
}

#endif
