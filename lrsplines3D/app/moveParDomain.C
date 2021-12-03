#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#include <vector>

using std::vector;
using namespace Go;

int main (int argc, char *argv[]) {

  if (argc != 6) {
    std::cout << "usage: ./moveParDomain <input volume> <output volume> <lower left (x,y,z)> " << std::endl;
    return -1;
  }

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  double x0 = atof(argv[3]);
  double y0 = atof(argv[4]);
  double z0 = atof(argv[5]);

  ObjectHeader oh;
  oh.read(ifs);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs);

 
  // Update parameter domain
  double umin = vol->paramMin(XDIR);
  double umax = vol->paramMax(XDIR);
  double vmin = vol->paramMin(YDIR);
  double vmax = vol->paramMax(YDIR);
  double wmin = vol->paramMin(ZDIR);
  double wmax = vol->paramMax(ZDIR);
  
  vol->setParameterDomain(x0, x0+umax-umin, y0, y0+vmax-vmin,
			  z0, z0+wmax-wmin);

  vol->writeStandardHeader(ofs);
  vol->write(ofs);
}


