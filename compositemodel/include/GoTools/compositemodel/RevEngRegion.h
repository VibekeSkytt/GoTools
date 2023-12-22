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

#ifndef _REVENGREGION_H
#define _REVENGREGION_H

#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/utils/BoundingBox.h"
//#include "GoTools/compositemodel/HedgeSurface.h"
#include <set>


namespace Go
{
  class HedgeSurface;
  //class RevEngPoint;
  class Circle;
  class SplineCurve;
  class Plane;
  class Cylinder;
  class Sphere;
  class Cone;
  class Torus;
  class SplneSurface;
  
  enum
    {
     CLASSIFICATION_UNDEF, CLASSIFICATION_CURVATURE, CLASSIFICATION_SHAPEINDEX, CLASSIFICATION_POINTASSOCIATION
    };

  // Method used in edge classification
  enum
    {
     TRIANGULATION_EDGE, PCATYPE_EDGE, CURVATURE_EDGE, CNESS_EDGE, RPFAC_EDGE
    };

  // Preference for elementary. Level
  enum
    {
     ALWAYS_ELEM, PREFER_ELEM, BEST_ACCURACY
    };

  // Surface flag
  enum
    {
     ACCURACY_OK, ANGULAR_DEVIATION, ACCURACY_POOR, PROBABLE_HELIX, FEW_POINTS, NOT_SET
    };
  
  struct SweepData
  {
    int type_;  // Linear = 1, rotational = 2, cylinderlike = 3, conelike = 4
    shared_ptr<SplineCurve> profile_;
    Point location_;
    Point added_info_;
    double radius_;
    double angle_;
    double maxdist_;
    double avdist_;
    int num_in_;

    SweepData(int type, shared_ptr<SplineCurve> profile, Point location, 
	      Point info2, double maxdist, double avdist, int num_in,
	      double radius = -1.0, double angle = 0.0)
    {
      type_ = type;
      profile_ = profile;
      location_ = location;
      added_info_ = info2;
      maxdist_ = maxdist;
      avdist_ = avdist;
      num_in_ = num_in;
      radius_ = radius;
      angle_ = angle;
    }
  };

  
  class RevEngRegion
  {
  public:
    // Constructor
    RevEngRegion(int edge_class_type);

    RevEngRegion(int classification_type, int edge_class_type);

    RevEngRegion(int classification_type, int edge_class_type,
		 std::vector<RevEngPoint*>& points);

    ~RevEngRegion();
    
     int getClassificationType()
    {
      return classification_type_;
    }

    bool isCompatible(ClassType classtype, int sfcode);
    
    // Extend group
    void addPoint(RevEngPoint* point);

    // Remove. NB! Leaves overview information invalid.
    void removePoint(RevEngPoint* point);

    void addRegion(RevEngRegion* reg,
		   std::vector<std::pair<double, double> >& dist_ang,
		   double maxd=0.0, double avd=0.0, int num_inside=-1,
		   int num_inside2=-1);
    
    // Update overview information
    void updateInfo();

    // Extend region with adjacent points having the same classification
    void collect(RevEngPoint *pt, RevEngRegion* prev=0);

    int numPoints()
    {
      return (int)group_points_.size();
    }

    RevEngPoint* getPoint(int ix)
    {
      if (ix < 0 || ix >= (int)group_points_.size())
	return 0;
      else
	return group_points_[ix];
    }

    std::vector<RevEngPoint*>::iterator pointsBegin()
    {
      return group_points_.begin();
    }
    
    std::vector<RevEngPoint*>::iterator pointsEnd()
    {
      return group_points_.end();
    }
    
    std::vector<RevEngPoint*> getPoints()
    {
      return group_points_;
    }

    const BoundingBox& getBbox()
    {
      return bbox_;
    }

    BoundingBox getParameterBox();
    
    RevEngPoint* seedPointPlane(int min_next, double rfac, double angtol);
    
    void growLocalPlane(Point mainaxis[3], double tol,
			std::vector<RevEngPoint*>& plane_pts,
			shared_ptr<Plane>& plane_out);
    void growPlaneOrCyl(Point mainaxis[3], int min_pt_reg,
			double tol, double angtol,
			std::vector<RevEngRegion*>& grown_regions,
			std::vector<HedgeSurface*>& adj_surfs,
			std::vector<std::vector<RevEngPoint*> >& added_groups);
    // void growLocal(RevEngPoint* seed, double tol, double radius, int min_close,
    // 		   std::vector<RevEngPoint*>& out);

    void growWithSurf(Point mainaxis[3], int max_nmb, int min_pt_reg,
		      double tol, double angtol,
		      std::vector<RevEngRegion*>& grown_regions,
		      std::vector<HedgeSurface*>& adj_surfs,
		      bool use_base=false);

    bool mergePlanarReg(double zero_H, double zero_K, double tol,
			Point mainaxis[3],
			std::vector<RevEngRegion*>& grown_regions);
    
    void mergeAdjacentSimilar(double tol, double angtol,
			      std::vector<RevEngRegion*>& grown_regions,
			      std::vector<HedgeSurface*>& adj_surfs);

    void
    initPlaneCyl(int min_point, int min_pt_reg,
		 double tol, double angtol, Point mainaxis[3], 
		 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		 std::vector<vector<RevEngPoint*> >& out_groups,
		 std::vector<RevEngPoint*>& single_pts, bool& repeat);
    
    void
    splitComposedRegions(int classtype,
			 std::vector<shared_ptr<RevEngRegion> >& added_groups,
			 std::vector<RevEngPoint*>& single_pts);
    
    void
    splitWithShapeIndex(std::vector<shared_ptr<RevEngRegion> >& updated_regions);
    
    void
    splitFromSurfaceNormals(std::vector<RevEngPoint*>& smallrad,
			    std::vector<std::vector<RevEngPoint*> >& separate_group);
    void splitRegion(std::vector<std::vector<RevEngPoint*> >& separate_groups);

    void updateRegion(double approx_tol, double anglim,
		      std::vector<RevEngRegion*>& adapted_regions,
		      std::vector<shared_ptr<RevEngRegion> >& outdiv_regions);

     void joinRegions(Point mainaxis[3], double approx_tol, double anglim,
		     std::vector<RevEngRegion*>& adapted_regions);

    void extractOutPoints(std::vector<std::pair<double, double> >& dist_ang,
			  double tol,
			  std::vector<std::vector<RevEngPoint*> >& out_groups);

     void identifyAngPoints(std::vector<std::pair<double, double> >& dist_ang,
			    double tol, double disttol,
			    std::vector<RevEngPoint*>& ang_points);
    void identifyDistPoints(std::vector<std::pair<double, double> >& dist_ang,
			    double tol, double maxd, double avd,
			    std::vector<RevEngPoint*>& dist_points);
    
    bool isConnected();
   void connectedGroups(std::vector<RevEngPoint*>& move,
			std::vector<std::vector<RevEngPoint*> >& out_groups,
			bool outer=false);
   void extractSpesPoints(std::vector<RevEngPoint*>& move,
			   std::vector<std::vector<RevEngPoint*> >& out_groups,
			   bool outer=false);
    
    int getClassification()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return group_points_[0]->C1_surf();
      else if (classification_type_ == CLASSIFICATION_SHAPEINDEX)
	return group_points_[0]->SI_surf();
      else if (classification_type_ == CLASSIFICATION_POINTASSOCIATION)
	return group_points_[0]->RP_surf();
      else
	return -1;
    }
    
   bool cylindertype()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return (group_points_[0]->C1_surf() == C1_RIDGE ||
		group_points_[0]->C1_surf() == C1_VALLEY);
      else if (classification_type_ == CLASSIFICATION_SHAPEINDEX)
	return (group_points_[0]->SI_surf() == SI_RUT ||
		group_points_[0]->SI_surf() == SI_RID);
      else
	return false;
    }
    
   bool planartype()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return (group_points_[0]->C1_surf() == C1_FLAT);
      else if (classification_type_ == CLASSIFICATION_SHAPEINDEX)
	return (group_points_[0]->SI_surf() == SI_PLANE);
      else
	return false;
    }
    
    bool extractPlane(Point mainaxis[3],
		      double tol, int min_pt, int min_pt_reg, double angtol,
		      int prefer_elementary,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool extractCylinder(double tol, int min_pt, int min_pt_reg, 
			 double angtol, int prefer_elementary,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups,
			 bool& repeat);

    bool feasiblePlane(double zero_H, double zero_K) const;

    bool feasibleCylinder(double zero_H, double zero_K) const;

    bool extractSphere(double tol, int min_pt, int min_pt_reg, 
		       double angtol, int prefer_elementary,
		       std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		       std::vector<HedgeSurface*>& prevsfs,
		       std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool extractLinearSweep(double tol, int min_pt, int min_pt_reg, 
			    double angtol, int prefer_elementary,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs);

    bool extractCone(double tol, int min_pt, int min_pt_reg, 
		     double angtol, int prefer_elementary,
		     std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		     std::vector<HedgeSurface*>& prevsfs,
		     std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool extractTorus(double tol, int min_pt, int min_pt_reg, 
		      double angtol, int prefer_elementary,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool contextTorus(Point mainaxis[3],
		      double tol, int min_pt, int min_pt_reg, 
		      double angtol, int prefer_elementary,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool contextCylinder(Point mainaxis[3],
			 double tol, int min_pt, int min_pt_reg,
			 double angtol, int prefer_elementary,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool contextCylinder(Point mainaxis[3],
			 double tol, int min_pt, int min_pt_reg, 
			 double angtol, int prefer_elementary,
			 std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups);

    bool adjacentToCylinder(Point mainaxis[3],
			    double tol, int min_pt, int min_pt_reg,
			    double angtol, int prefer_elementary,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::vector<std::vector<RevEngPoint*> >& out_groups);

    RevEngPoint* closestPoint(const Point& pos, double& dist);

    std::vector<RevEngPoint*> extractNextToAdjacent(RevEngRegion* reg);

    std::vector<RevEngPoint*> extractBdPoints();

    std::vector<RevEngPoint*> extractBranchPoints();

    void growFromNeighbour(Point mainaxis[3], int min_pt_reg,
			   std::vector<RevEngPoint*>& seed, double tol,
			   double angtol, RevEngRegion *neighbour);
    

    bool tryOtherSurf(int prefer_elementary, bool replace);
    
    bool extractFreeform(double tol, int min_pt, int min_pt_reg,
			 double angtol, int prefer_elementary,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups);

    void implicitizeSplit();

    void segmentByPlaneGrow(Point mainaxis[3], double tol,
			    double angtol, int min_pt, 
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::vector<std::vector<RevEngPoint*> >& out_groups);
    
    bool segmentByPlaneAxis(Point mainaxis[3], int min_point_in,
			    int min_pt_reg,
			    double tol, double angtol, int prefer_elementary,
			    std::vector<RevEngRegion*>& adj_planar,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<shared_ptr<RevEngRegion> >& added_reg,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::vector<std::vector<RevEngPoint*> >& added_groups);
    
    void setHedge(HedgeSurface* surface)
    {
      associated_sf_.clear();
      associated_sf_.push_back(surface);
      computeDomain();
    }

    void addHedge(HedgeSurface* surface)
    {
      associated_sf_.push_back(surface);
      computeDomain();
    }

    bool hasSurface()
    {
      return (associated_sf_.size() > 0);
    }

    int numSurface()
    {
      return (int)associated_sf_.size();
    }

    HedgeSurface* getSurface(int ix)
    {
      return (ix < 0 || ix>=(int)associated_sf_.size()) ? 0 : associated_sf_[ix];
    }

    void clearSurface()
    {
      if (associated_sf_.size() > 0)
	associated_sf_.clear();
      maxdist_ = avdist_ = variance_ = 0.0;
      num_inside_ = num_inside2_ = 0;
    }

    int getSurfaceFlag()
    {
      return surfflag_;
    }

    void setSurfaceFlag(int surfflag)
    {
      surfflag_ = surfflag;
    }

    void getDomain(double dom[4])
    {
      for (int ka=0; ka<4; ++ka)
	dom[ka] = domain_[ka];
    }

    
    void setDomain(double dom[4])
    {
      for (int ka=0; ka<4; ++ka)
	domain_[ka] = dom[ka];
    }

    bool accuracyOK(int min_pt, double tol, int num_inside, double avdist)
    {
      return (num_inside > min_pt && num_inside > (int)group_points_.size()/2 &&
	      avdist <= tol);
    }

    bool accuracyOK(int num_points, int min_pt, double tol, int num_inside,
		    double avdist)
    {
      return (num_inside > min_pt && num_inside > num_points/2 &&
	      avdist <= tol);
    }

    int defineSfFlag(int min_point, double tol, int num_in,
		     int num_in2, double avd, bool type_cyl);
    
    int defineSfFlag(int num_points, int min_point, double tol, int num_in,
		     int num_in2, double avd, bool type_cyl);
    
    void computeDomain();
    
    const BoundingBox& boundingBox()
    {
      return bbox_;
    }

    const DirectionCone& getNormalCone()
    {
      return normalcone_;
    }

    void getPrincipalCurvatureInfo(double& mink1, double& maxk1, double& mink2, double& maxk2)
    {
      mink1 = mink1_;
      maxk1 = maxk1_;
      mink2 = mink2_;
      maxk2 = maxk2_;
    }

    void getAvCurvatureInfo(double& avH, double& avK, double& MAH, double& MAK)
    {
      avH = avH_;
      avK = avK_;
      MAH = MAH_;
      MAK = MAK_;
    }

    Point getMeanNormal()
    {
      return avnorm_;
    }
    
    void setAccuracy(double maxdist, double avdist, int num_inside,
		     int num_inside2);

    // Maximum distance, average absolute distance, number of points
    // satisfying distance and angular tolerance, number of points
    // satisfying distance tolerance
    void getAccuracy(double& maxdist, double& avdist, int& num_inside,
		     int& num_inside2)
    {
      maxdist = maxdist_;
      avdist = avdist_;
      num_inside = num_inside_;
      num_inside2 = num_inside2_;
    }

    double getAverageDist()
    {
      return avdist_;
    }

    int getNumInside()
    {
      return num_inside_;
    }
    
    int getNumInside2()
    {
      return num_inside2_;
    }
    
    double getVariance()
    {
      return variance_;
    }

    double getMaxSfDist()
    {
      return maxdist_;
    }

    void setVisited(bool visited)
    {
      visited_ = visited;
    }

    bool visited()
    {
      return visited_;
    }

    bool possiblePlane(double angtol, double inlim);
    bool possibleCylinder(double angtol, double inlim);
    bool possibleCone(double angtol, double inlim);
    bool possibleTorus(double angtol, double inlim);

    bool hasSweepInfo()
    {
      return (sweep_.get() != 0);
    }

    int sweepType()
    {
      return (sweep_.get() ? sweep_->type_ : 0);
    }

    bool hasDivideInfo()
    {
      return false;
    }

    void setRegionAdjacency();

    void updateRegionAdjacency();

    bool integrateInAdjacent(double mean_edge_len, int min_next,
			     int max_next, double tol, double angtol,
			     int max_nmb_outlier, RevEngRegion* taboo=0);

    bool adjustWithCylinder(Point mainaxis[3],
			    double tol, double angtol, int min_pt_reg,
			    std::vector<std::vector<RevEngPoint*> >& out_groups,
			    std::vector<RevEngRegion*>& grown_regions,
			    std::vector<HedgeSurface*>& adj_surfs);
    
    void getAdjInsideDist(shared_ptr<ParamSurface> surf, double dom[4],
			  double tol, RevEngRegion* reg,
			  double& avd, double& ava, int& nn,
			  std::vector<RevEngPoint*>& adjpts,
			  std::vector<double>& par_and_dist);

    void addAdjacentRegion(RevEngRegion* adj_reg)
    {
      //adj_reg->addAdjacentRegion(this);
      adjacent_regions_.insert(adj_reg);
    }
    
    void removeAdjacentRegion(RevEngRegion* adj_reg)
    {
      //adj_reg->removeAdjacentRegion(this);
      if (std::find(adjacent_regions_.begin(), adjacent_regions_.end(), adj_reg) != adjacent_regions_.end())
	adjacent_regions_.erase(adj_reg);
      // else
      // 	std::cout <<"Something wrong in adjacent regions" << std::endl;
    }
    
    void clearRegionAdjacency()
    {
      adjacent_regions_.clear();
    }
    
    bool hasAdjacentRegion(RevEngRegion* adj_reg)
    {
      return (adjacent_regions_.find(adj_reg) != adjacent_regions_.end());
    }

    int numAdjacentRegions()
    {
      return (int)adjacent_regions_.size();
    }

    void getAdjacentRegions(std::vector<RevEngRegion*>& adjacent)
    {
      adjacent.insert(adjacent.end(), adjacent_regions_.begin(), adjacent_regions_.end());
    }

    std::vector<RevEngRegion*> fetchAdjacentPlanar();

    std::vector<RevEngRegion*> fetchAdjacentCylindrical();

    bool segmentByAdjSfContext(Point mainaxis[3], int min_point_in, 
			       int min_pt_reg, double tol, double angtol,
			       std::vector<RevEngRegion*>& adj_planar,
			       std::vector<std::vector<RevEngPoint*> >& added_groups);

    void sortByAxis(vector<Point>& axis, double tol,
		    std::vector<std::vector<RevEngPoint*> >& groups1,
		    std::vector<std::vector<RevEngPoint*> >& groups2,
		    std::vector<RevEngPoint*>& remaining);

    bool extractCylByAxis(Point mainaxis[3], int min_point, int min_pt_reg,
			  double tol, double angtol, int prefer_elementary,
			  std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			  std::vector<shared_ptr<RevEngRegion> >& added_reg,
			  std::vector<std::vector<RevEngPoint*> >& out_groups,
			  std::vector<RevEngPoint*>& single_pts);

    Point directionFromAdjacent(double angtol);
    
    bool segmentByDirectionContext(int min_point_in, double tol,
				   const Point& direction, double angtol,
				   std::vector<std::vector<RevEngPoint*> >& added_groups);

     
    void getAdjacentElemInfo(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*>  >& adj_elem,
			     std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*>  >& adj_elem_base);

    void setPreviousReg(RevEngRegion *prev)
    {
      prev_region_ = prev;
    }
    
    void adjustBoundaries(double mean_edge_len, double tol, double angtol);
    
    std::vector<RevEngPoint*> sortPtsSeq(double mean_edge_len,
					 std::vector<RevEngPoint*>& seq_pts,
					 std::vector<RevEngPoint*>& sub_pts);

    vector<RevEngPoint*>  extractBdOutPoints(shared_ptr<SplineCurve>& crv,
					     std::vector<RevEngPoint*>& seq_pts,
					     double tol);

    void adjustWithSurf(Point mainaxis[3], int min_pt_reg, double tol, double angtol);

    double getFracNorm()
    {
      return frac_norm_in_;
    }
    
    bool hasBaseSf()
    {
      return basesf_.get();
    }

    void setBaseSf(shared_ptr<ParamSurface> base, double maxd, double avd,
		   int num_in, int num_in2)
    {
      basesf_ = base;
      maxdist_base_ = maxd;
      avdist_base_ = avd;
      num_in_base_ = num_in;
      num_in_base2_ = num_in2;
    }

    shared_ptr<ParamSurface> getBase()
    {
      return basesf_;
    }
    
    void getBase(shared_ptr<ParamSurface>& base, double& maxd, double& avd,
		 int& num_in, int& num_in2)
    {
      base = basesf_;
      maxd = maxdist_base_;
      avd = avdist_base_;
      num_in = num_in_base_;
      num_in2 = num_in_base2_;
    }
    
    void getBaseDist(double& maxd, double& avd, int& num_in, int& num_in2)
    {
      maxd = maxdist_base_;
      avd = avdist_base_;
      num_in = num_in_base_;
      num_in2 = num_in_base2_;
    }

    bool hasEdgeBetween(RevEngRegion* adj);

    void writeRegionInfo(std::ostream& of);
    void writeRegionPoints(std::ostream& of);
    void writeAdjacentPoints(std::ostream& of);
    void writeUnitSphereInfo(std::ostream& of);
    void writeSubTriangulation(std::ostream& of);
    void writeSurface(std::ostream& of);

    void store(std::ostream& os) const;
    void read(std::istream& is, shared_ptr<ftPointSet>& tri_sf,
	      std::vector<int>& associated_sf_id);
	      
  private:
    std::vector<RevEngPoint*> group_points_;   // Points belonging to classified segment
    int classification_type_;
    int edge_class_type_;
    std::vector<HedgeSurface*> associated_sf_;  // Can be two due to split along
    // seam of closed surface (should be fixed)
    int surfflag_;
    double domain_[4];
    shared_ptr<ImplicitApprox> impl_;
    double mink1_, maxk1_, mink2_, maxk2_;
    double avH_, avK_, MAH_, MAK_;
    BoundingBox bbox_;
    DirectionCone normalcone_;
    double frac_norm_in_;
    Point avnorm_;
    double maxdist_, avdist_, variance_;
    int num_inside_, num_inside2_;
    std::set<RevEngRegion*> adjacent_regions_;
    RevEngRegion* prev_region_;
    shared_ptr<ParamSurface> basesf_;
    double maxdist_base_, avdist_base_;
    int num_in_base_, num_in_base2_;

    shared_ptr<SweepData> sweep_;
    bool visited_;

    //Point& pluckerAxis();
    void extendWithGaussRad();
    void extendWithGaussRad2();
    void analyseNormals(double tol, Point& normal, Point& centre, double& radius);
    void analysePlaneProperties(Point avnorm, double angtol,
				std::vector<RevEngPoint*>& in,
				std::vector<RevEngPoint*> out);
    void analyseCylinderProperties(Point avvec, double angtol,
				   std::vector<RevEngPoint*>& in,
				   std::vector<RevEngPoint*> out);
    void  curveApprox(std::vector<Point>& points, double tol,
		      shared_ptr<Circle> circle,
		      std::vector<double>& param,
		      shared_ptr<SplineCurve>& curve, Point& xpos);
    void  curveApprox(std::vector<Point>& points,
		      shared_ptr<ParamCurve> cvin,
		      int ik, int in, 
		      shared_ptr<SplineCurve>& curve);
    void configSplit(std::vector<RevEngPoint*>& points,
		     std::vector<double>& param,
		     shared_ptr<Cylinder> cyl,
		     shared_ptr<SplineCurve> spl, double tol,
		     std::vector<std::vector<RevEngPoint*> >& configs);
    shared_ptr<Plane> computePlane(std::vector<RevEngPoint*>& points,
				   const Point& norm_dir, Point mainaxis[3]);
    shared_ptr<Cylinder>
    computeCylinder(std::vector<RevEngPoint*>& points, double tol);
    void analyseCylProject(shared_ptr<Cylinder> cyl, double tol,
			   std::vector<std::vector<RevEngPoint*> >& configs);
    void analyseCylRotate(shared_ptr<Cylinder> cyl, double tol,
			  double avdist, int num_in, double& avdist_lin,
			  int& num_in_lin, double& avdist_cub,
			  int& num_in_cub, shared_ptr<Cone>& cone);

    shared_ptr<Sphere> computeSphere(std::vector<RevEngPoint*>& points);
    shared_ptr<Cone> computeCone(std::vector<RevEngPoint*>& points, Point& apex);
    shared_ptr<Torus> computeTorus(std::vector<RevEngPoint*>& points,
				   double tol, shared_ptr<Torus>& torus2);
    shared_ptr<SplineSurface> computeLinearSwept(double tol, shared_ptr<SplineCurve>& profile,
						 Point& pt1, Point& pt2);
    shared_ptr<SplineSurface> computeFreeform(double tol);
    shared_ptr<SplineSurface> updateFreeform(double tol);
    void getPCA(double lambda[3], Point& eigen1, Point& eigen2, Point& eigen3);
    shared_ptr<SplineSurface> surfApprox(vector<RevEngPoint*>& points,
					 const BoundingBox& bbox);
    void splitCylinderRad(const Point& pos, const Point& axis,
			  const Point& Cx, const Point& Cy,
			  int nmb_split, std::vector<Point>& centr,
			  std::vector<double>& rad);
    void approximationAccuracy(std::vector<RevEngPoint*>& points,
			       shared_ptr<ParamSurface> surf,
			       double tol, double angtol,
			       double& maxd, double& avd,
			       std::vector<RevEngPoint*>& in,
			       std::vector<RevEngPoint*>& out);
    void checkReplaceSurf(Point mainaxis[3], int min_pt_reg, double tol,
			  double angtol, bool always=false);
    bool parameterizeOnSurf(shared_ptr<ParamSurface> surf,
			    std::vector<double>& data,
			    std::vector<double>& param,
			    int& inner1, int& inner2, bool& close1, bool& close2);
    bool reparameterize(std::vector<double>& param, std::vector<double>& param2,
			double& umin, double& umax, double& vmin, double& vmax);

    void getParExtent(double curr[2], int pdir, std::vector<std::vector<int> >& raster,
		      int& i1, int& i2);

    void defineRaster(std::vector<double>& param, int nmb_div,
		      std::vector<std::vector<int> >& raster, double& umin,
		      double& umax, double& vmin, double& vmax);

    void extendInCorner(std::vector<double>& data, std::vector<double>& param,
			double umin, double umax, double vmin, double vmax);

    void includeAdjacentRegion(RevEngRegion* reg, double maxd, double avd,
			       int num_inside, int num_inside2,
			       std::vector<double>& parvals,
			       std::vector<std::pair<double, double> >& dist_ang,
			       std::vector<RevEngRegion*>& added_adjacent);

    bool computeIntegrateInfo(std::vector<RevEngPoint*>& points, RevEngRegion *adj_reg,
			      double tol, double angtol, double radius, bool local_approx, 
			      int min_next, int max_next, int max_nmb_outlier, 
			      bool& outlier, int& nmb_pt_adj, double& maxdist, 
			      double& avdist, int& nmb_in, double& maxdist_adj, 
			      double& avdist_adj, int& nmb_in_adj,
			      int& nmb_in_adj2);

        bool
    analyseTorusContext(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
			double tol, double angtol, std::vector<size_t>& adjacent_ix,
			int& plane_ix, int& cyl_ix, Point& pos, Point& axis,
			Point& Cx, double& R1, double& R2, double cyl_dom[4],
			bool& outer);
    bool
    analyseCylinderContext(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
			   double tol, double angtol, Point mainaxis[3],
			   std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_planar,
			   Point& pos, Point& axis, Point& Cx, double& rad);

    bool planarComponent(Point vec, int min_point, int min_pt_reg, double tol,
			 double angtol, Point mainaxis[3],
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<vector<RevEngPoint*> >& out_groups,
			 std::vector<RevEngPoint*>& single_pts);

    bool defineHelicalInfo(shared_ptr<Cylinder> cyl,  double tol,
			   double angtol, int min_point, int min_pt_reg,
			   double avdist, int num_in1, int num_in2,
			   std::vector<std::pair<double,double> >& dist_ang,
			   std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			   std::vector<std::vector<RevEngPoint*> >& out_groups,
			   std::vector<RevEngPoint*>& single_pts);
    
    bool integratePlanarPoints(std::vector<Point>& dir,
			       std::vector<std::vector<RevEngPoint*> >& groups,
			       std::vector<std::pair<shared_ptr<ElementarySurface>,RevEngRegion*> >& adj_elem,
			       double tol, double angtol,
			       std::vector<RevEngPoint*>& remaining);
    
    bool defineCylindricalRegs(Point mainaxis[3],
			       std::vector<std::vector<RevEngPoint*> >& groups,
			       int min_point, int min_pt_reg,
			       double tol, double angtol,
			       std::vector<shared_ptr<RevEngRegion> >& added_reg,
			       std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			       std::vector<RevEngPoint*>& remaining);
    
    void axisFromAdjacent(double angtol, std::vector<Point>& axis);
  };
}

#endif
