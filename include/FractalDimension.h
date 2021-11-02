#ifndef FRACTALDIMENSION_H
#define FRACTALDIMENSION_H

#include <EvaluateTree.h>

bool in_box(const LGMdouble& x0, const LGMdouble& y0, const LGMdouble& z0,
           const LGMdouble& dd, const LGMdouble& x, const LGMdouble& y, const LGMdouble& z);

void determineFractalDimension(Tree<ETCfSegment,ETCfBud>& tree, int argc, char** argv, int no_points,
			       string tre_file);

// Sets Points on the surfaces of shoot cylinders (needle or woody).
// End disk of only last cylinders in the axis is covered with Points

// In order to find out if is the last segment in the Axis, must do with
// PropagateUp with the last segmment the Axis

class  PointSurface {
 public:
 PointSurface(list<Point>& s, const int np) : point_surface(s), n_points(np) {}
  ETCfSegment*
    operator()(ETCfSegment* last_seg, TreeCompartment<ETCfSegment,ETCfBud>* tc)const
  {
    if(Axis<ETCfSegment,ETCfBud>* ax = 
       dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(tc)) {
      last_seg = dynamic_cast<ETCfSegment*>(GetLastTreeSegment(*ax));
      return last_seg;
    }
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      LGMdouble r, l;
      if(GetValue(*ts, LGAAf) > R_EPSILON) {
        r = GetValue(*ts, LGARf);
      }
      else {
        r = GetValue(*ts, LGAR);
      }
      if(r < R_EPSILON)
        return last_seg;

      l = GetValue(*ts,LGAL);

      //Division of points on outer surface (= rectangle 2*PI*R x L)
      int nr = (int)sqrt((double)n_points*2.0*PI_VALUE*r/l);
      int nl = (int)((double)n_points/(double)nr);

      LGMdouble dl = l/(double)nl;
      LGMdouble dfii = 2.0*PI_VALUE/(double)nr;

      PositionVector dir = GetDirection(*ts);
      dir.normalize();

      PositionVector u;
      PositionVector up(0.0,0.0,1.0);
      if((dir == up) || ((PositionVector(0.0,0.0,0.0)-dir) == up))
        u = PositionVector(0.0,1.0,0.0);
      else
        u = up;

      PositionVector pointer = Cross(dir,u);
      pointer.normalize();

      Point base = GetPoint(*ts);

      PositionVector zero(0.0,0.0,0.0);
      LGMdouble rot = dfii/2.0;
      cout << endl;
      for(LGMdouble fii = dfii/2.0; fii < 2.0*PI_VALUE; fii += dfii) {
        pointer.rotate(zero,dir,rot);
        pointer.normalize();
        for(LGMdouble x = dl/2.0; x < l; x += dl) {
          Point xpo = Point(PositionVector(base) + x * dir);
          Point po = Point(PositionVector(xpo) + r*pointer);
          point_surface.push_back(po);
          //      cout << PositionVector(base) << " # " << po;
        }
        rot = dfii;
      }

      //Points to the end disk of last tree segment in an Axis
      if(ts == last_seg) {
        //in the proportion point on the outer surface
        int n_end_points = (int)(((PI_VALUE*r*r)/(2.0*PI_VALUE*r*l))*(double)n_points);
        //End points are generated to a square, correct number of points
        //with proportion of their areas
        n_end_points = (int)((4.0/PI_VALUE)*(double)n_end_points);
        //Finally, points to x and y direction
        int nx_y = sqrt((double)n_end_points) + 1;
	if(nx_y < 2)
	  nx_y = 2;

        //x-direction = direction of pointer
        PositionVector y_vec = Cross(pointer,dir);
        y_vec.normalize();
        LGMdouble dxy = 2.0*r/((double)nx_y-1.0);

        for(int i = 0; i < nx_y; i++) {
          for(int j = 0; j < nx_y; j++) {
            LGMdouble xi = (double)i*dxy - r;   //values are -r, -r+dxy, ..., r
            LGMdouble yj = (double)j*dxy - r;
            if(sqrt(xi*xi+yj*yj) > r)
              continue;
            Point p = Point(PositionVector(base)+l*dir + xi*pointer
                            + yj*y_vec);
            point_surface.push_back(p);
          }
        }
      }  //if(ts == last_seg) ...
    }
    return last_seg;
  }
 private:
  list<Point>& point_surface;
  int n_points;         //Desired number of points on surface of one segment
};


#endif
