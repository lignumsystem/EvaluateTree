#ifndef STAR_MEAN_H
#define STAR_MEAN_H
#include <iostream>
#include <string>
#include <sstream>
#include <mathsym.h>
#include <vector>

using namespace std;
using namespace cxxadt;
#include <ParseCommandLine.h>
#include <Point.h>
#include <PositionVector.h>
#include <ParametricCurve.h>

#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1


/* int CylinderBeamShading(const Point& r0, const PositionVector& b,  */
/*                         const Point& rs, const PositionVector& a, */
/*                         double Rs, double Rw, double L,  */
/*                         double& distance ); */

int CylinderBeamShading(const Point& r0_1, const PositionVector& b, 
                        const Point& rs_1, const PositionVector& a,
                        double Rs, double L, double& distance );

int CylinderBeamShading(const Point& r0_1, const PositionVector& b, 
                        const Point& rs_1, const PositionVector& a,
                        double Rs, double L);

extern int ran3_seed;

class sm_main_program {
 public:
  sm_main_program() {
    int dummy = 0;
    K = ParametricCurve("0.000 0.24 0.262 0.214 0.524 0.195"
			" 0.785 0.178 1.047 0.170 1.309 0.171 1.571 0.190",
			dummy);
  }

  void usage()
  {
    cout << "Usage:  ./starmean runs [-Rs <value=0.02>] [-L <value=0.2>]" << endl;
    cout << "[-seed <value>] [-fol_den <Value=150>]" << endl;
 
  }
  bool read_command_line(int argc, char** argv);
  void create_segment(double theta);
  void evaluateWithRw(int& ran3_seed);
  void evaluate(int& ran3_seed);
  void create_randomsegment();
  void randomevaluate();

  //-----------------------------------
  int runs;
  double Rs, L, Rw;
  double fol_den, Af;
  ParametricCurve K;
  int hit, no_hit;
  vector<double> starm;    //starmean in 20 inclinations
  PositionVector seg_dir;
  Point seg_base;
  double tot_len;
  double star_sum;
  double proj;
};

bool sm_main_program::read_command_line(int argc, char** argv)
{
  bool success = false;

  if(argc < 2) {
    usage();
   return success;
 }

  runs = atoi(argv[1]);
 
 string clarg;
 
 Rs = 0.02;
 if (ParseCommandLine(argc,argv,"-Rs", clarg))
   if(clarg.length() > 0)
     Rs = atof(clarg.c_str());

 L = 0.2;
 if (ParseCommandLine(argc,argv,"-L", clarg))
   if(clarg.length() > 0)
     L = atof(clarg.c_str());

 if (ParseCommandLine(argc,argv,"-seed", clarg)){
   if (clarg.length() > 0){
     ran3_seed = atoi(clarg.c_str());
     ran3_seed = -abs(ran3_seed);
   }
 }
 ran3(&ran3_seed);
 
 fol_den = 50;
 if (ParseCommandLine(argc,argv,"-fol_den", clarg)){
   if (clarg.length() > 0){
     fol_den = atof(clarg.c_str());
   }
 }

 Af = PI_VALUE*Rs*Rs*L*fol_den;
 
 success = true;
 return success;

}  //END READ COMMAND LINE

void sm_main_program::create_segment(double theta) {

  seg_dir = PositionVector(cos(theta),0.0,sin(theta));
  seg_base = Point(Rs, Rs, 1.0);
 
}


void sm_main_program::evaluate(int& ran3_seed)
{
  double x = (L + 2.0*Rs) * ran3(&ran3_seed);
  double y = (2.0*Rs) * ran3(&ran3_seed);
  Point r0(x,y,0.0);
  PositionVector up(0.0,0.0,1.0);
  //  double Rw = 0.0;
  double distance = 0.0;
/*   int res = CylinderBeamShading(r0, up, seg_base, seg_dir, */
/* 				Rs, Rw, L, distance ); */
  int res = CylinderBeamShading(r0, up, seg_base, seg_dir,
				Rs, L, distance);

  if(res > 0){        //HIT_THE_FOLIAGE  == 1
       hit++;
       double a_dot_b = Dot(seg_dir,up);
       double incl = PI_DIV_2 - acos(fabs(a_dot_b));
       double k_ext = K(incl);
       star_sum += exp(-distance*fol_den* k_ext);
       tot_len += distance;
     }
  else
    no_hit++;
}

//This gives the same result (if fol_den > 0) as evaluate
void sm_main_program::evaluateWithRw(int& ran3_seed)
{
  double x = (L + 2.0*Rs) * ran3(&ran3_seed);
  double y = (2.0*Rs) * ran3(&ran3_seed);
  Point r0(x,y,0.0);
  PositionVector up(0.0,0.0,1.0);
  double distance = 0.0;
  int res = 0; //CylinderBeamShading(r0, up, seg_base, seg_dir,
               //	Rs, Rw, L, distance);
  if(fol_den == 0) {
    res = CylinderBeamShading(r0, up, seg_base, seg_dir,Rw, L);
    if(res < 0) {
      hit++;
      star_sum += 0.0;
	} else {
      no_hit++;
    }
  } 
  else {
    //If foliage, wood radius = 0
    res = CylinderBeamShading(r0, up, seg_base, seg_dir, Rs, L, distance);
    if(res > 0){
      hit++;
      double a_dot_b = Dot(seg_dir,up);
      double incl = PI_DIV_2 - acos(fabs(a_dot_b));
      double k_ext = K(incl);
      star_sum += exp(-distance*fol_den* k_ext);
      tot_len += distance;
    }
    else {
      no_hit++;
    }
  }
}


void sm_main_program::create_randomsegment() {

  double fii = 2.0*PI_VALUE*ran3(&ran3_seed);
  double theta = asin(ran3(&ran3_seed));
  if(ran3(&ran3_seed) < 0.5)
    theta = -theta;
  PositionVector pv(cos(fii)*cos(theta),sin(fii)*cos(theta),sin(theta));

  seg_dir =  PositionVector(cos(fii)*cos(theta),sin(fii)*cos(theta),sin(theta));
  seg_base = Point(L+Rs,L+Rs, 1.0);
 
  proj = 2.0*L*Rs*cos(theta)+PI_VALUE*Rs*Rs*sin(abs(theta));
}

void sm_main_program::randomevaluate()
{
  double x = 2.0*(L + Rs) * ran3(&ran3_seed);
  double y = 2.0*(L + Rs) * ran3(&ran3_seed);
  Point r0(x,y,0.0);
  PositionVector up(0.0,0.0,1.0);
  double Rw = 0.0;
  double distance = 0.0;
  //  int res = CylinderBeamShading(r0, up, seg_base, seg_dir,
  //			Rs, Rw, L, distance );
  int res = CylinderBeamShading(r0, up, seg_base, seg_dir,
				Rs, L, distance);
  if(res > 0){
       hit++;
       double a_dot_b = Dot(seg_dir,up);
       double incl = PI_DIV_2 - acos(fabs(a_dot_b));
       double k_ext = K(incl);
       star_sum += (1.0-exp(-distance*fol_den* k_ext))*proj/Af;
       tot_len += distance;
/*        cout << "dist fold k star" << distance << " " << fol_den << " " << k_ext */
/* 	    << " " << (1.0-exp(-distance*fol_den* k_ext)) << endl; */
     }
  else
    no_hit++;
 
}

double SAc(double phi, double r, double l)
{
  return 2 * l * cos(phi) * r + PI_VALUE * r * r * sin(phi);
}



//
//
//
double S(double phi, double Sf, double Wf, double r, double l)
{
  
    int dummy = 0;
    ParametricCurve  K("0.000 0.24 0.262 0.214 0.524 0.195"
			" 0.785 0.178 1.047 0.170 1.309 0.171 1.571 0.190",
			dummy);
	
	if (Sf * Wf == 0)
		return 0;
	if (SAc(phi, r, l) == 0)
		return 0;
	
	return SAc(phi, r, l)*(1 - exp(-K(phi)*Sf*Wf/SAc(phi, r, l)))/(Sf*Wf);
}


//Main program has been transformed to function that calculates STAR_mean
double star_mean(int runs, double Rs, double L, double fol_den, int& r_seed)
{

  sm_main_program sm;
  
  sm.runs = runs;
  sm.Rs = Rs;
  sm.L = L;
  sm.fol_den = fol_den;
  sm.Af =  PI_VALUE*Rs*Rs*L*fol_den;

  vector<pair<double, double> > star_and_cos(11);

  for(int k = 0; k < 11 ; k++) {
    double theta = PI_VALUE * (double)k /(2.0 * 10.0);
    sm.create_segment(theta);
    sm.star_sum = 0.0;
    sm.hit = 0;
    sm.no_hit = 0;
    sm.tot_len = 0.0;
    sm.proj = 2.0*sm.L*sm.Rs*cos(theta)+PI_VALUE*sm.Rs*sm.Rs*sin(theta);

    for(int i = 0; i < sm.runs; i++) {
      sm.evaluate(r_seed);
    }
    if(sm.hit > 0)
      star_and_cos[k].first = (1.0-sm.star_sum/(double)sm.hit)/(sm.Af/sm.proj);
    else
      star_and_cos[k].first = 0.0;
    star_and_cos[k].second = cos(theta);
  }
  double sum = 0, csum = 0;
  for(int k = 0; k < 11 ; k++) {
    sum += star_and_cos[k].first*star_and_cos[k].second;
    csum += star_and_cos[k].second;
  }

  double star;
  if(csum > 0.0)
    star = sum/csum;
  else star = 0.0;

  return star;
}

//Function that calculates STAR value 
//Modified from star_mean; the inclination of the shoot is given as input
//Also radius of woody part is considered

double star(int runs, double Rs, double Rw, double L, double fol_den, double incl, int& r_seed)
// Rs = radius of needle cylinder
// Rw = radius of woody part
// L = length of shoot
// fol_den = foliage density in the shoot (needle area / needle volume (m2/m3))
// fii = inclination angle of the shoot as defined by Oker-Blom and Smolander 1988
//       (fii = 0 means horizontal shoot, radiation is coming from zenith)  
{

  sm_main_program sm;
  
  sm.runs = runs;
  sm.Rs = Rs;
  sm.Rw = Rw;
  sm.L = L;
  sm.fol_den = fol_den;
  sm.Af =  PI_VALUE*(Rs*Rs-Rw*Rw)*L*fol_den;

    sm.create_segment(incl);
    sm.star_sum = 0.0;
    sm.hit = 0;
    sm.no_hit = 0;
    sm.tot_len = 0.0;
    sm.proj = 2.0*sm.L*sm.Rs*cos(incl)+PI_VALUE*sm.Rs*sm.Rs*sin(incl);

    for(int i = 0; i < sm.runs; i++) {
      sm.evaluateWithRw(r_seed);
    }
    double star;
    if(sm.hit > 0)
      star = (1.0-sm.star_sum/(double)sm.hit)/(sm.Af/sm.proj);
    else
      star = 0.0;

  return star;
}

#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD


#endif
