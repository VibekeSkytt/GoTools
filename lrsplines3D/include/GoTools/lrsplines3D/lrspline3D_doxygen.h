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

#ifndef _LRSPLINE3D_DOXYGEN_H
#define _LRSPLINE3D_DOXYGEN_H

/**
\page lrspline3D_doc GoTools LR Splines 3D Module

\section intro_sec Introduction
The lrsplines3D module in GoTools provides functionalities for working with LR (Locally Refined) B-spline volumes. These volumes are a powerful tool for representing complex 3D geometries with local refinement capabilities, offering flexibility in controlling the level of detail. The core of this module revolves around the LRSplineVolume class, which combines an underlying Mesh3D with a set of LRBSpline3D basis functions.

\section structure_sec Module Structure

The key components and their relationships within the lrsplines3D module are as follows:


\link Go::LRSplineVolume \endlink  represents a three-dimensional LR B-spline volume.

This class is the central component for defining, manipulating, and evaluating 3D LR spline volumes. An LRSplineVolume is essentially a sum of LRBSpline3D basis functions, each weighted by a control point, and defined over an adaptive Mesh3D.

\details
The LRSplineVolume manages the geometric and topological information of the spline volume. It holds references to the underlying Mesh3D structure, which defines the knot intervals and element connectivity. It also maintains a collection of LRBSpline3D basis functions, which are the fundamental building blocks of the spline volume. The class provides methods for evaluating the volume at specific parameter points, evaluating on a grid, and performing refinement operations.

\sa Go::Mesh3D
\sa Go::LRBSpline3D

\par Key Data Members:

    Go::Mesh3D& mesh_: The underlying 3D mesh that defines the parametric domain and local refinement structure.
    std::vector<LRBSpline3D*> b_splines_: A collection of 3D LR B-spline basis functions that form the volume.
    BSplineMap bsplines_: An internal map of individual B-spline basis functions.
    ElementMap emap_: A map holding information about the individual mesh elements and their associated basis functions.

\par Core Evaluation Methods:

    void point(Point& pt, double upar, double vpar, double wpar) const: \brief Evaluates the LRSplineVolume at a given (upar, vpar, wpar) parameter triplet. \param pt A Go::Point object where the evaluated 3D point will be stored. \param upar The u-parameter for evaluation. \param vpar The v-parameter for evaluation. \param wpar The w-parameter for evaluation.
    void point(Point& pt, double upar, double vpar, double wpar, Element3D* elem) const: \brief Evaluates the LRSplineVolume at (upar, vpar, wpar) using an optional Element3D hint for faster lookup. \param pt A Go::Point object where the evaluated 3D point will be stored. \param upar The u-parameter for evaluation. \param vpar The v-parameter for evaluation. \param wpar The w-parameter for evaluation. \param elem A pointer to an Element3D that might contain (upar, vpar, wpar), used to optimize the evaluation process.
    void point(std::vector<Point>& pts, double upar, double vpar, double wpar, int derivs, bool u_from_right = true, bool v_from_right = true, bool w_from_right = true, double resolution = 1.0e-12) const: \brief Evaluates the LRSplineVolume and its partial derivatives up to a specified order. \param pts A std::vector<Go::Point> to store the evaluated point and its derivatives. The order of results is typically (V, Vu, Vv, Vw, Vuu, Vuv, Vuw, Vvv, Vvw, Vww) up to the specified derivs order. \param upar The u-parameter for evaluation. \param vpar The v-parameter for evaluation. \param wpar The w-parameter for evaluation. \param derivs The maximum order of derivatives to compute (e.g., 0 for just the point, 1 for first derivatives, etc.). \param u_from_right Boolean indicating evaluation from the right in u-direction for boundary conditions. \param v_from_right Boolean indicating evaluation from the right in v-direction for boundary conditions. \param w_from_right Boolean indicating evaluation from the right in w-direction for boundary conditions. \param resolution A tolerance for numerical computations.
    void elementGridEvaluate(Element3D *element, std::vector<double>& upar, std::vector<double>& vpar, std::vector<double>& wpar, std::vector<double>& points) const: \brief Evaluates the volume on a grid of points within a specific Element3D. \param element A pointer to the Element3D within which to evaluate the grid. \param upar A std::vector<double> containing the u-parameter values for the grid. \param vpar A std::vector<double> containing the v-parameter values for the grid. \param wpar A std::vector<double> containing the w-parameter values for the grid. \param points A std::vector<double> where the flattened (x,y,z, x,y,z, ...) evaluated 3D points will be stored.
    void elementGridEvaluate(Element3D *element, double* upar, int usize, double* vpar, int vsize, double *wpar, int wsize, int deriv, std::vector<double>& points) const: \brief Evaluates the volume and its derivatives on a grid of points within an Element3D. \param element A pointer to the Element3D within which to evaluate the grid. \param upar A C-style array of u-parameter values. \param usize The size of the upar array. \param vpar A C-style array of v-parameter values. \param vsize The size of the vpar array. \param wpar A C-style array of w-parameter values. \param wsize The size of the wpar array. \param deriv The maximum order of derivatives to compute. \param points A std::vector<double> where the flattened evaluated points and their derivatives (F, Fu, Fv, Fw, etc.) will be stored.

\par Other Key Methods:

    LRSplineSurface* defineBivariate(int pardir, double parval, bool at_start) const: \brief Defines a constant parameter surface (a 2D LR spline surface) from the 3D volume. \param pardir The parameter direction (e.g., XDIR, YDIR, ZDIR) for the constant plane. \param parval The parameter value at which the surface is defined. \param at_start Boolean indicating if the surface is defined at the start of the domain in the given direction. \return A pointer to the newly created LRSplineSurface.

\link Go::LRBSpline3D \endlink represents a single 3D LR B-spline basis function.

These are the fundamental basis functions used to construct an LRSplineVolume. Each LRBSpline3D has its own local support defined by knot vectors in three parametric directions (u, v, w).

\details
An LRBSpline3D is characterized by its degree and knot vector indices for each parametric direction. They are combined linearly to form the LRSplineVolume.

\sa Go::LRSplineVolume

\section utility_classes Utility Classes and Functions

\subsection LRBSpline3DUtils
\brief Provides utility functions for operations specific to LRBSpline3D objects.
\details This header (LRBSpline3DUtils.h) contains functions like split_function which can split an LRBSpline3D into two new basis functions along a new knot, useful for local refinement.

\subsection LRSpline3DUtils
\brief Provides general utility functions for 3D LR splines.
\details This header (LRSpline3DUtils.h) includes helpers for common tasks related to LRSplineVolume objects, such as collecting basis functions that overlap a specified area. It also defines comparison functors for LRBSpline3D based on their support.

\subsection LRSplinePlotUtils3D
\brief Provides utilities for plotting and visualization of LRSplineVolume objects.
\details This header (LRSplinePlotUtils3D.h) includes functions such as writeElementLineCloud to output element grid lines and writePostscriptMesh to generate PostScript representations of the parametric mesh.

\subsection LRBenchmarkUtils3D
\brief Provides utilities for benchmarking performance, particularly related to refinement operations on LRSplineVolume.
\details This header (LRBenchmarkUtils3D.h) contains functions like benchmarkVolRefinement to measure the performance of volume refinement operations.

*/


#endif
