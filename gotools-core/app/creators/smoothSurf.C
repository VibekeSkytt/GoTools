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

#include "GoTools/creators/ModifySurf.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <fstream>


using namespace Go;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 3) {
	MESSAGE("Usage: inputfile outputfile.");
	return 0;
    }

    // Read input arguments
    std::ifstream filein(argv[1]);
    ALWAYS_ERROR_IF(filein.bad(), "Input file not found or file corrupt");

    std::ofstream fileout(argv[2]);

    // Input surface should be a Go::SplineSurface
    ObjectHeader header;
    header.read(filein);
    shared_ptr<SplineSurface> surf(new SplineSurface());
    surf->read(filein);

    DirectionCone cone1 = surf->normalCone();
    std::cout << "Input normal cone angle: " << cone1.angle() << std::endl;

    if (surf->numCoefs_u() == surf->order_u())
      {
	surf->insertKnot_u(0.75*surf->startparam_u()+0.25*surf->endparam_u());
	surf->insertKnot_u(0.5*(surf->startparam_u()+surf->endparam_u()));
	surf->insertKnot_u(0.25*surf->startparam_u()+0.75*surf->endparam_u());
      }

    if (surf->numCoefs_v() == surf->order_v())
      {
	surf->insertKnot_v(0.75*surf->startparam_v()+0.25*surf->endparam_v());
	surf->insertKnot_v(0.5*(surf->startparam_v()+surf->endparam_v()));
	surf->insertKnot_v(0.25*surf->startparam_v()+0.75*surf->endparam_v());
      }

    ModifySurf::smoothSurface(surf, 1);

    DirectionCone cone2 = surf->normalCone();
    std::cout << "Output normal cone angle: " << cone1.angle() << std::endl;

    surf->writeStandardHeader(fileout);
    surf->write(fileout);
}

