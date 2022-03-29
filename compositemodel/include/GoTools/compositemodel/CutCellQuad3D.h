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

#ifndef _CUTCELLQUAL3D_H_
#define  _CUTCELLQUAL3D_H_

#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/Body.h"

namespace Go
{
  class Body;

  class CutCellQuad3D
  {
  public:
    // Constructor
    CutCellQuad3D(shared_ptr<SurfaceModel> sf_model, double tol);

    // Define quadrature information
    void setQuadratureInfo(std::vector<double>& quadpar,
			   std::vector<double> weights,
			   double min_cell_size)
    {
      quadpar_ = quadpar;
      weights_ = weights;
      min_cell_size_ = min_cell_size;
    }
    
    // Check cell status
    int cellStat(const Point& ll, const Point& ur, int& coinc);

    // Compute quadrature points and weights
    void quadrature(const Point& ll, const Point& ur,
		    std::vector<double>& quadraturepoints,
		    std::vector<double>& pointsweights,
		    std::vector<std::vector<shared_ptr<ParamSurface> > >& unresolved_cells,
		    std::vector<double>& surfquads,
		    std::vector<double>& surfnorms,
		    std::vector<double>& sfptweights,
		    std::vector<std::vector<shared_ptr<ParamSurface> > >& small_sfs,
		    int stat = -1, int coinc = -1);

 private:
    double tol_;
    shared_ptr<Body> body_;
    std::vector<double> quadpar_;
    std::vector<double> weights_;
    double min_cell_size_;

    void createCutCell(const Point& ll, const Point& ur,
		       std::vector<shared_ptr<Body> >& cutcell);

    void createCellSfs(const Point& ll, const Point& ur,
		       std::vector<shared_ptr<SplineSurface> >& cell_sfs);
    void quadraturePoints(const Point& ll, const Point& ur,
			  shared_ptr<Body> body,
			  std::vector<double>& quadraturepoints,
			  std::vector<double>& pointsweights,
			  std::vector<std::vector<shared_ptr<ParamSurface> > >& unresolved_cells,
			  std::vector<double>& surfquads,
			  std::vector<double>& surfnorms,
			  std::vector<double>& sfptweights,
			  std::vector<std::vector<shared_ptr<ParamSurface> > >& small_sfs);


    int selectBaseDir(Point& ll, Point& ur,
		      shared_ptr<Body> body, int& splitdir, double& splitval,
		      shared_ptr<SplineCurve>& rulecv, Point& ruledir);

    
    void splitCell(shared_ptr<Body> body,
		   shared_ptr<ParamSurface> splitsf,
		   std::vector<shared_ptr<Body> >& subcell);

    void removeCoincFaces(shared_ptr<SurfaceModel>& mod1,
			  shared_ptr<SurfaceModel>& mod2,
			  shared_ptr<SurfaceModel>& mod3,
			  double tol);

    void fetchSharpEdges(shared_ptr<Body> body, std::vector<ftEdge*>& convex,
			 std::vector<ftEdge*>& concave);
   };
} // end namespace Go

#endif
