#include <FractalDimension.h>
//#include <TMatrix3D.h>

//Checks if the point (x,y,z) is in box that has lower left corner
//(x0,y0,z0) and upper right corner (x0+dd,y0+dd,z0+dd)

bool in_box(const LGMdouble& x0, const LGMdouble& y0, const LGMdouble& z0,
           const LGMdouble& dd, const LGMdouble& x, const LGMdouble& y, const LGMdouble& z) {

  if(x < x0)
    return false;
  if(y < y0)
    return false;
  if(z < z0)
    return false;
  if(x > x0 + dd)
    return false;
  if(y > y0 + dd)
    return false;
  if(z > z0 + dd)
    return false;

  return true;
}

//=================================================================================
// Fractal Dimension
//
//=================================================================================

void determineFractalDimension(Tree<ETCfSegment,ETCfBud>& tree, int argc, char** argv, int no_points,
			       string tree_file) {

  // Store the points on surfaces of segments (needle cylinder or woody cylinder
  // if no needles) to a list os points

  list<Point> surfaces;                     //All points on the surfaces of shoot cylinders
  PointSurface ps(surfaces, no_points);

// Only top of an Axis covered with points --> must start with
// PropagateUp with the last segment of the main axis

  ETCfSegment* ts = 
    dynamic_cast<ETCfSegment*>(GetLastTreeSegment(GetAxis(tree)));
  PropagateUp(tree, ts, ps);

  cout << "The tree is covered by " << surfaces.size() << " points." << endl;
  cout << endl;

  //Box spanning the tree
  BoundingBox bb;
  FindCfBoundingBox<ETCfSegment,ETCfBud> fb;      //Bounding box including also
  //the segments that have no foliage
  bb = Accumulate(tree, bb, fb);
      
  Point ll = bb.getMin();
  Point ur = bb.getMax();
  LGMdouble xspan = ur.getX() - ll.getX();
  LGMdouble yspan = ur.getY() - ll.getY();
  LGMdouble zspan = ur.getZ() - ll.getZ();
  LGMdouble xmin = ll.getX();
  LGMdouble ymin = ll.getY();
  LGMdouble zmin = ll.getZ();
  LGMdouble xmax = ur.getX();
  LGMdouble ymax = ur.getY();
  LGMdouble zmax = ur.getZ();
     
  LGMdouble minbox = 150.0;

  string cla;
  if(ParseCommandLine(argc,argv,"-fractalDm",cla)) {
    minbox = atof(cla.c_str());
  }

  LGMdouble span = min(xspan,min(yspan, zspan));
  cout << "Original span is " << span << endl;

  LGMdouble d = 1.0, dprev = 1.0;
  cla.clear();
  if(ParseCommandLine(argc,argv,"-fractalD0",cla)) { //this starts "d iteration" from
    d = atof(cla.c_str());                           // a given value
    dprev = d/10.0;              // an ad hoc value
  }
  else {
    d = 1.0;
    dprev = d;
  }

string::size_type p = tree_file.find(".xml");
string::size_type slash = tree_file.rfind("/");
string f_name;
if(slash == string::npos) {
  f_name = tree_file.substr(0,p) + "-f.dat";
 }
 else {
   f_name = tree_file.substr(slash+1,p-slash-1) + "-f.dat";
 }

ofstream ff(f_name.c_str(), ofstream::trunc);

  cout << "d  1/d log(1/d) n_boxes n_occupied log(n_occupied) span npoints" << endl;
  ff << "d  1/d log(1/d) n_boxes n_occupied log(n_occupied) npoints span npoints" << endl;
  LGMdouble dd = 100.0;
  while(d < minbox && dd > 0.05) {   //box side must be longer than 0.05 m = 5 cm
    long int n_occupied = 0;
    long int ix = 0, iy = 0, iz = 0;
    dd = span/d;
    for(LGMdouble x = xmin; x <= xmax; x += dd) {
      ix++;
      for(LGMdouble y = ymin; y <= ymax; y += dd) {
	iy++;
	for(LGMdouble z = zmin; z <= zmax; z += dd) {
	  iz++;
	  //              bool already = false;
	  list<Point>::iterator I;
	  for(I = surfaces.begin(); I != surfaces.end(); I++) {
	    if(in_box(x,y,z,dd,(*I).getX(),(*I).getY(),(*I).getZ())) {
	      n_occupied++;
	      break;
	    }
	  }
	}
      }
    }
    cout << dd << " " << 1.0/dd << " " << log(1.0/dd) << " " << ix*iy*iz << " " << n_occupied 
	 << " " << log(n_occupied) <<  " " << span << " " << no_points <<  endl;
    ff << dd << " " << 1.0/dd << " " << log(1.0/dd) << " " << ix*iy*iz << " " << n_occupied 
       << " " << log(n_occupied) <<  " " << span << " " << no_points << endl;

    if(dd < 0.03)
      break;

    LGMdouble dnew = d + dprev;
    dprev = d;
    d = dnew;
  }
  ff.close();
  return;
}

