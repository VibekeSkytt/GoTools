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

#ifndef _CHANGEDETECTION_H
#define _CHANGEDETECTION_H

#include "GoTools/utils/BoundingBox.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"
#include "GoTools/geometry/SurfaceTools.h"
#include <vector>

namespace Go
{

  class BoundedSurface;
  
  class ChangeDetection
  {
    struct SubSeqsInfo
    {
      std::vector<std::vector<double> > subseq_pts_;
      std::vector<int> pts_seq_ix_;
      std::vector<shared_ptr<LRSplineSurface> > subseq_sfs_;
      std::vector<double> maxdist_;
      std::vector<double> avdist_;
      std::vector<double> avdistout_;
      std::vector<int> nmb_out_;
      std::vector<shared_ptr<LRSplineSurface> > subseq_diffsfs_;
      std::vector<double> min_diff_;
      std::vector<double> max_diff_;
      std::vector<double> av_diff_;
    };
    
  public:
    ChangeDetection(std::vector<std::vector<double> >& point_clouds, int del);

    double suggestTimeStep();

    void defineTrivariate(double time_step);

    int getInitNumCoefs()
    {
      return init_num_coef_;
    }

    void setInitNumCoefs(int num_coefs)
    {
      init_num_coef_ = num_coefs;
    }

    int getDegree()
    {
      return degree_;
    }

    void setDegree(int degree)
    {
      degree_ = degree;
    }

    void volApprox(double tol, int max_iter);

    shared_ptr<LRSplineVolume> getApproxVol()
    {
      return vol_;
    }

    void identifyChanges(double limit);

    void extractChangeData();
    
    void surfApprox(double tol, int max_iter);

    void differenceSurfaces();

    void analyseDiffSurfaces(double threshold, double eps);

  private:
    std::vector<std::vector<double> > point_seqs_;   // Point clouds for each
    // acquisition
    int del_;   // Number of double for each point
    std::vector<BoundingBox> seq_box_;  // Bounding box for each acquisition
    BoundingBox av_box_;   // Average box for all time steps
    double time_delta_;  // Time interval between acquisitions in trivariate
    // data set
    std::vector<double> time_;  // Time parameter associated to point cloud

    std::vector<double> tri_pts_;  // Trivariate point cloud
    double domain_[6];             // Trivariate parameter domain
    int init_num_coef_;            // Initial number of coefficients in each parameter
    // direction (prior to distribution)
    int degree_;                   // Polynomial degree (same in all parameter directions)
    double mba_level_;             // Initial height in approximation
    shared_ptr<LRSplineVolume> vol_;  // Volume approximation
    double max_dist_, av_dist_out_, av_dist_all_;  // Approximation accuracy
    int num_out_;

    double der_lim_;
    std::vector<BoundingBox> change_box_;

    std::vector<SubSeqsInfo> change_pts_;

    void analyseOneDiffSurface(shared_ptr<BoundedSurface> bdsurf,
			       int nmb_area, int nmb_diff,
			       double delta, double eps);
    void
    materialMove(shared_ptr<BoundedSurface> bdsurf,
		 const std::vector<CurveVec>& curves,
		 std::vector<double>& isovals, int ix, int sgn,
		 std::vector<std::pair<std::vector<shared_ptr<ParamCurve> >, double> > mat);

    void
    getContourLoops(shared_ptr<BoundedSurface> bdsurf,
		    const CurveVec& curves, double isoval, int sgn, double eps,
		    std::vector<std::vector<shared_ptr<CurveOnSurface> > >& loops);

  };

} // end namespace Go

#endif  // _CHANGEDETECTION_H
