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

#ifndef _REVENG_H
#define _REVENG_H

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"
#include <string.h>

namespace Go
{
  class RevEngPoint;
  //class RevEngRegion;
  class HedgeSurface;
  class RevEngEdge;

  // Elementary surface types to recognize (omitting currently
  // ellipsoid, elliptic cylinder, ...)
  enum
  {
   PLANE, CYLINDER, SPHERE, CONE, TORUS
  };

  // Characterization of model surface
  enum
  {
   SMOOTH=0, MEDIUM_ROUGH, ROUGH
  };
  
  struct SurfaceProperties
  {
    int sfix_;
    ClassType type_;
    ClassType prev_type_;
    Point dir_, loc_;
    double rad1_, rad2_;
    int num_points_;
    int surfflag_;

    SurfaceProperties(int ix, ClassType type, int num, int surfflag, Point& dir, 
		      Point& loc, ClassType prev_type=Class_Unknown, double rad1=-1,
		      double rad2=-1)
    {
      sfix_ = ix;
      type_ = type;
      num_points_ = num;
      surfflag_ = surfflag;
      dir_ = dir;
      loc_ = loc;
      prev_type_ = prev_type;
      rad1_ = rad1;
      rad2_ = rad2;
    }
  };

  struct AxisInfo
  {
    Point axis_;
    std::vector<std::pair<Point,int> > plane_loc_;
    std::vector<std::pair<Point,int> > rotational_loc_;

    AxisInfo(Point& axis)
    {
      axis_ = axis;
    }

    void addPlaneLocation(Point& loc, int num)
    {
      plane_loc_.push_back(std::make_pair(loc,num));
    }
    
    void addRotationalLocation(Point& loc, int num)
    {
      rotational_loc_.push_back(std::make_pair(loc,num));
    }
  };
    
  class RevEng
  {
  public:
    RevEng();
    
    RevEng(shared_ptr<ftPointSet> tri_sf, double mean_edge_len);

    ~RevEng();

    // Should this class have an option to run all operations in one
    // sequence without being started from outside?
      
    void enhancePoints();

    void edgeClassification();
    void classifyPoints();

    void segmentIntoRegions();

    void initialSurfaces();

    void growSurfaces();

    void updateAxesAndSurfaces();

    void firstEdges();

    void surfaceCreation(int pass);

    void manageBlends1();

    void manageBlends2();

    shared_ptr<SurfaceModel> createModel();

    void storeClassified(std::ostream& os) const;
    void readClassified(std::istream& is);
    void storeGrownRegions(std::ostream& os);
    void readGrownRegions(std::istream& is);
    void curvatureFilter();

    double getInitApproxTol();
    void setApproxTolerance();
    void setApproxTol(double eps)
    {
      approx_tol_= eps;
    }
    
    double getApproxTol()
    {
      return approx_tol_;
    }

    void setEdgeClassificationParams();
    int getEdgeClassificationType()
    {
      return edge_class_type_;
    }
    
    void setEdgeClassificationType(int edge_class_type)
    {
      edge_class_type_ = edge_class_type;
    }
    
    double getCfac()
    {
      return cfac_;
    }
    void setCfac(double cfac)
    {
      cfac_ = cfac;
    }

    void setClassificationParams();
    int getClassificationType()
    {
      return classification_type_;
    }
    void setClassificationType(int classification_type)
    {
      classification_type_ = classification_type;
    }

    double getMeanCurvatureZero()
    {
      return zero_H_;
    }

    void setMeanCurvatureZero(double zero_H)
    {
      zero_H_ = zero_H;
    }

    
    double getGaussCurvatureZero()
    {
      return zero_K_;
    }

    void setGaussCurvatureZero(double zero_K)
    {
      zero_K_ = zero_K;
    }

    int getElementaryPreferLevel()
    {
      return prefer_elementary_;
    }

    void setElementaryPreferLevel(int preferlevel)
    {
      prefer_elementary_ = preferlevel;
    }
    
    int getModelCharacterization()
    {
      return model_character_;
    }

    void setModelCharacterization(int character)
    {
      //model_character_ = std::min(ROUGH, std::max(SMOOTH, character));
      model_character_ = std::min(2, std::max(0, character));
    }

    void setMainAxis(Point mainaxis[3])
    {
      mainaxis_[0] = mainaxis[0];
      mainaxis_[1] = mainaxis[1];
      mainaxis_[2] = mainaxis[2];
    }

    void getMainAxis(Point mainaxis[3])
    {
      mainaxis[0] = mainaxis_[0];
      mainaxis[1] = mainaxis_[1];
      mainaxis[2] = mainaxis_[2];
    }

    // Prelimenary results
    int numRegions()
    {
      return (int)regions_.size();
    }

    shared_ptr<RevEngRegion> getRegion(int ix)
    {
      return regions_[ix];
    }

    int numSurfaces()
    {
      return (int)surfaces_.size();
    }

    shared_ptr<HedgeSurface> getSurface(int ix)
    {
      return surfaces_[ix];
    }

    // Could be used to elect if a group of points should be visualized
    double getMinPointRegion()
    {
      return min_point_region_;
    }

    void trimSurfaces();

    
    void updateRegionsAndSurfaces(int& ix, std::vector<RevEngRegion*>& grown_regions,
				  std::vector<HedgeSurface*>& adj_surfs);

   void smallRegionSurfaces();

    void growSmallRegionSurface(int& ix);

    void adaptToMainAxis();

   private:
    int model_character_;
    shared_ptr<ftPointSet> tri_sf_;
    double mean_edge_len_;
    std::vector<shared_ptr<RevEngRegion> > regions_;
    std::vector<RevEngPoint*> single_points_;
    std::vector<shared_ptr<HedgeSurface> > surfaces_;  // I think the 
    // surfaces must be collected here to have a stable storage
    // The surfaces can be freeform as well as primary. The collection
    // will be build gradually. The number of surfaces will increase and
    // decrease based on recognition, merging and splitting by trimming
    std::vector<shared_ptr<RevEngEdge> > edges_;  // Intersection curves
    // between surfaces with additional information
    shared_ptr<SurfaceModel> sfmodel_;
    BoundingBox bbox_;
    int min_next_;  // Minimum number of neighbouring points
    int max_next_;  // Estimate for maximum number of neighbouring points
    double rfac_;   // Factor for radius in which to search for neighbouring points
    int edge_class_type_ = CURVATURE_EDGE;
    int classification_type_ = CLASSIFICATION_CURVATURE;
    double cfac_;   // Edge points from curvature is given by
    // cfac_ times the average length of triangulation edges in a vertex
    double norm_ang_lim_; // Limit for when the cone angle corresponding
    // to triangle normals indicate an edge
    double norm_plane_lim_;  // Limit for when the cone angle corresponding
    // to triangle normals indicate a plane
    double zero_H_;  // When mean curvature is considered zero
    double zero_K_;  // When Gauss curvature is considered zero
    int min_point_region_;
    double approx_tol_;  // Approximation tolerance in region growing
    double int_tol_;  // Intersection tolerance
    double anglim_;
    int max_nmb_outlier_;

    int prefer_elementary_; // 0 = always, 1 = preferred, 2 = best accuracy

    Point mainaxis_[3];
    std::vector<AxisInfo> model_axis_;
    
    struct SmallSurface
    {
      int axis_ix_, pos_ix_, lev_ix_;
      vector<vector<RevEngPoint*> > assos_points_;
      vector<shared_ptr<ElementarySurface> > surfs_;
      vector<BoundingBox> bbox_;
      int type_;   // 1 = plane1, 2=plane2, 3=rotational

      SmallSurface(int ix1, int ix2, int ix3, int type,
		   vector<shared_ptr<ElementarySurface> >& surfs)
      {
	axis_ix_ = ix1;
	pos_ix_ = ix2;
	lev_ix_ = ix3;
	type_ = type;
	surfs_ = surfs;
      }

      void addPoints(vector<RevEngPoint*>& points, BoundingBox bb)
      {
	assos_points_.push_back(points);
	bbox_.push_back(bb);
      }
    };
  
    void initParameters();
    void updateParameters();
    bool recognizeOneSurface(int& ix, int min_point_in, double angtol,
			     int pass);
    void recognizeSurfaces(int min_point_in, int pass);
    void recognizeEdges(bool only_curve=false);
    bool createBlendSurface(int ix);
    void adjustPointRegions(int min_point_in);

    void computeAxisFromCylinder(Point initaxis[3], int min_num, double max_ang,
				 Point axis[3], int num_points[3]);
    
    void computeAxisFromPlane(Point initaxis[3], int min_num, double max_ang,
			      Point axis[3], int num_points[3]);
    
    void surfaceExtractOutput(int idx,
			      std::vector<std::vector<RevEngPoint*> > out_groups,
			      std::vector<HedgeSurface*> prev_surfs);

    bool segmentComposite(int& ix, int min_point_in, double angtol);
    bool segmentByPlaneGrow(int ix, int min_point_in, double angtol);
    bool segmentByAxis(int ix, int min_point_in);
    bool segmentByContext(int ix, int min_point_in, double angtol, bool first);
    void growSurface(int& ix, int pass = 1);
    void growBlendSurface(int& ix);
    void growMasterSurface(int& ix);

    void defineSmallRegionSurfaces();
    
    bool identifySmallRotational(std::vector<RevEngPoint*>& points,
				 Point midp, Point loc, Point axis, Point Cx,
				 double ppar1, double ppar2,
				 std::vector<shared_ptr<ElementarySurface> >& sfs);

    bool identifySmallPlanar(std::vector<RevEngRegion*>& groups,
			     Point loc, Point axis, Point Cx,
			     double ppar1, double ppar2, double delta,
			     std::vector<shared_ptr<ElementarySurface> >& sfs);

    bool identifySmallPlanar(std::vector<RevEngPoint*>& groups,
			     Point loc, Point axis, Point Cx,
			     double ppar1, double ppar2, double delta,
			     std::vector<shared_ptr<ElementarySurface> >& sfs);
    
    void planarAtPlane(shared_ptr<Plane> axis_plane,
		       std::vector<RevEngPoint*>& points,
		       std::vector<HedgeSurface*>& sfs,
		       std::vector<shared_ptr<RevEngRegion> >& plane_sf_reg,
		       std::vector<shared_ptr<HedgeSurface> >& plane_sf_hedge);
    
    void integrateInSmallSurfs(std::vector<shared_ptr<RevEngRegion> >& small_sf_reg,
			       std::vector<RevEngRegion*>& nosf_reg,
			       std::vector<RevEngRegion*>& include_reg);
    
    void extractSmallSurfs(SmallSurface& small_surf,
			   std::vector<shared_ptr<RevEngRegion> >& small_sf_reg,
			   std::vector<shared_ptr<HedgeSurface> >& small_sf_hedge,
			   std::vector<RevEngRegion*>& nosf_reg,
			   std::vector<RevEngPoint*>& non_assigned_pts);

    void doAdaptToAxis();

    bool axisUpdate(int ix, double max_ang, double angtol);
    
    void adjustWithMainAxis(std::vector<Point>& axes, std::vector<int>& num_pts);
    Point planarFit(std::vector<int>& sf_ix, Point axis);

    Point rotationalFit(std::vector<int>& sf_ix, Point axis, Point Cx,
			std::vector<RevEngEdge*>& nopar_edgs);

    void collectAxis(std::vector<SurfaceProperties>& sfprop);

    void computeMonge(RevEngPoint* pt, std::vector<RevEngPoint*>& points,
		      Point& vec1, Point& vec2, double radius);

    int setSmallRegionNumber();
    
    void storeParams(std::ostream& os) const;
    void readParams(std::istream& is);
    void setBoundingBox();
    
    void writeRegionStage(std::ostream& of, std::ostream& ofm, std::ostream& ofs) const;
    void writeRegionWithSurf(std::ostream& of) const;
    void checkConsistence(std::string text) const;

    std::vector<shared_ptr<RevEngEdge> >
    defineEdgesBetween(size_t ix1,shared_ptr<ElementarySurface>& surf1,
		       Point& dir1, size_t ix2, shared_ptr<ElementarySurface>& surf2,
		       Point& dir2, bool only_curve=false,
		       double lenlim=-1.0, bool check_common = true);
    
    shared_ptr<RevEngEdge> 
    defineOneEdge(size_t ix1, shared_ptr<ElementarySurface>& surf1,
		  Point& dir1, size_t ix2,
		  shared_ptr<ElementarySurface>& surf2, Point& dir2,
		  shared_ptr<CurveOnSurface>& int_cv1,
		  shared_ptr<CurveOnSurface>& int_cv2,
		  double width, std::vector<RevEngRegion*>& common_reg,
		  bool only_curve);
    
    RevEngPoint* getDistantPoint(shared_ptr<CurveOnSurface>& cv,
				 double tmin, double tmax, double dd,
				 double width,
				 std::vector<RevEngPoint*>& points);

    void extendBlendAssociation(size_t ix);
    
    bool setBlendEdge(size_t ix);
    
    double
    computeCylinderRadius(std::vector<std::vector<RevEngPoint*> > blend_pts,
			  double width, const Point& pos, const Point& axis,
			  const Point& dir1, const Point& dir2);
   shared_ptr<Cylinder>
    createCylinderBlend(std::vector<std::vector<RevEngPoint*> > blend_pts,
			double rad1, const Point& pos, const Point& axis,
			const Point& dir1, const Point& dir2, int sgn);
    
    double
    computeTorusRadius(std::vector<std::vector<RevEngPoint*> >& blend_pts,
		       shared_ptr<CurveOnSurface>& cv,
		       const Point& locp, const Point& normal,
		       shared_ptr<ElementarySurface> rotational,
		       double width, bool plane_out, bool rot_out);
    
    void getTorusParameters(shared_ptr<ElementarySurface> planar,
			    shared_ptr<ElementarySurface> rotational,
			    double radius, int sgn1, int sgn2, double& Rrad, 
			    Point& centre, Point& normal, Point& Cx);
    
    double
    computeTorusRadius(std::vector<std::vector<RevEngPoint*> >& blend_pts,
		       shared_ptr<CurveOnSurface>& cv,
		       shared_ptr<ElementarySurface> elem1,
		       shared_ptr<ElementarySurface> elem2,
		       double width, bool out1, bool out2, int sgn,
		       double& d2);
    
  bool
    getTorusParameters(shared_ptr<ElementarySurface> elem1,
		       shared_ptr<ElementarySurface> elem2, Point pos,
		       double radius, double d2, bool out1, bool out2, int sgn,
		       double& Rrad, Point& centre, Point& normal, Point& Cx);

  shared_ptr<Torus>
    torusBlend(std::vector<std::vector<RevEngPoint*> >& blend_pts,
	       std::vector<shared_ptr<CurveOnSurface> >& cvs,
	       const Point& locp, const Point& normal,
	       shared_ptr<ElementarySurface> rotational,
	       double width, bool plane_out, bool rot_out);
    
    shared_ptr<Torus>
    torusBlend(std::vector<std::vector<RevEngPoint*> >& blend_pts,
	       shared_ptr<CurveOnSurface>& cv,
	       shared_ptr<ElementarySurface> elem1,
	       shared_ptr<ElementarySurface> elem2,
	       double width, bool out1, bool out2, int sgn);
    
    void setBlendBoundaries(RevEngRegion *reg);

    void equalizeBlendRadii();

    void equalizeAdjacent(size_t ix1, size_t ix2);
    
    void updateBlendRadius(size_t ik, double radius);
    
    bool defineTorusCorner(size_t ix);

    void defineMissingCorner(std::vector<RevEngRegion*>& cand_adj);

    bool createTorusBlend(size_t ix);

    bool suitcaseCorner(std::vector<RevEngRegion*>& adj_blends,
			RevEngEdge *rev_edge);

    void extractOutPoints(std::vector<RevEngPoint*>& points, shared_ptr<ParamSurface> surf,
			  std::vector<int>& cv_ix,
			  double tol, double angtol,
			  std::vector<std::vector<RevEngPoint*> >& move2adj,
			  std::vector<RevEngPoint*>& remain);
  };

} // namespace Go

#endif // _REVENG_H
