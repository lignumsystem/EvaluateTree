#ifndef SETSHOOTPOSITIONS_H
#define SETSHOOTPOSITIONS_H
//#include <cmath>
//#include <cstdio>
//#include <fstream>
//#include <utility>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
//#include <XMLTree.h>

extern int ran3_seed;

bool inEllipsoid(const Point &p, const Point& pc, const LGMdouble& rh, const LGMdouble& rv) {
  LGMdouble rh2 = pow(rh, 2.0);
  return (sqrt(pow(p.getX()-pc.getX(),2.0)/rh2+pow(p.getY()-pc.getY(),2.0)/rh2+
	       pow(p.getZ()-pc.getZ(),2.0)/pow(rv,2.0)) <= 1.0);
} 

class  SetShootPositions {
 public:
 SetShootPositions(const Point& pci, const LGMdouble& rhi, const LGMdouble& rvi, const double hc, const bool be,
		   const bool wp) :
  pc(pci), rh(rhi), rv(rvi), Hc(hc), both_ends(be), woody_parts(wp) {
    xmin = pc.getX() - rh;
    xmax = pc.getX() + rh;
    ymin = pc.getY() - rh;
    ymax = pc.getY() + rh;
    zmin = pc.getZ() - rv;
    zmax = pc.getZ() + rv;
  }

  TreeCompartment<ETCfSegment,ETCfBud>*
    operator()(TreeCompartment<ETCfSegment,ETCfBud>* tc)const
    {
      if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
	if(GetPoint(*ts).getZ() < Hc) {
	  return tc;
	}

	if(GetValue(*ts, LGAAf) > R_EPSILON || woody_parts) {
	  bool not_in_ellipsoid = true;
	  Point p_mid;    //Mid point of the segment
	  Point p_start;
	  Point p_end;
	  PositionVector dir;

	  while(not_in_ellipsoid) { 
	    p_mid.setX(2.0* rh * ran3(&ran3_seed) + xmin);
	    p_mid.setY(2.0* rh * ran3(&ran3_seed) + ymin);
	    p_mid.setZ(2.0* rv * ran3(&ran3_seed) + zmin);

	    LGMdouble theta = asin(ran3(&ran3_seed));      //inclination theta is in the range [0,pi/2] 
	    if(ran3(&ran3_seed) < 0.5)
	      theta = -theta;

	    LGMdouble fii = 2.0*PI_VALUE*ran3(&ran3_seed);
	    dir = PositionVector(cos(fii)*cos(theta),sin(fii)*cos(theta),sin(theta));
	    LGMdouble l = GetValue(*ts,LGAL);
	    p_start = p_mid - (l/2.0) * (Point)dir;

	    if(!both_ends) {
	      if(!inEllipsoid(p_mid, pc, rh, rv)) {
		continue;
	      }
	    } else {
	      if(!inEllipsoid(p_start, pc, rh, rv)) {
		continue;
	      }
	      p_end = p_mid + (l/2.0) * (Point)dir;
	      if(!inEllipsoid(p_end, pc, rh, rv)) {
		continue;
	      }
	    }
	    not_in_ellipsoid = false;  //is in the ellipsoid


	  } // while(not_in_ellipsoid) ...

	  SetPoint(*ts,p_start);
	  SetDirection(*ts,dir);
	} //if(GetValue(*ts, LGAAf) ...
      }
      return tc;
    }
 private:
  Point pc;
  LGMdouble rh, rv;
  LGMdouble xmin, xmax, ymin, ymax, zmin,zmax;
  double Hc;
  bool both_ends, woody_parts;
};
#endif
