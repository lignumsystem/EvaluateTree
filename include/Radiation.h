#ifndef RADIATION_H
#define RADIATION_H
#include <star_mean.h>



//Calculates STAR and other values of segments and writes to a file. 
//Argument nruns in the constructor specifies the number of trials (beams of
//radiation shot through the segment randomly) for the accurate STAR calculation.
//If it is high, it takes a long time to calculate.

class CalculateSTAR {
 public:
  CalculateSTAR(const string& filename, int nruns): fname(filename), runs(nruns)
    {
      ofstream f(fname.c_str() , ofstream::trunc);
      f << "Height Dist_top Dist_stem age Rf L Af Wf Vf fol_den STAR STAR_eq" << endl;
      f.close();
    }

    TreeCompartment<ETCfSegment,ETCfBud>*
      operator()(TreeCompartment<ETCfSegment,ETCfBud>* tc) const
      {
        if (ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
          LGMdouble Wf = GetValue(*ts,LGAWf);
          if(Wf > 0.0) {
            LGMdouble Rf = GetValue(*ts,LGARf);
            LGMdouble L = GetValue(*ts,LGAL);
            LGMdouble Age = GetValue(*ts, LGAage);
            LGMdouble Af = GetValue(*ts, LGAAf);
            Tree<ETCfSegment,ETCfBud>& t = GetTree(*ts);
            LGMdouble H = GetValue(t, LGAH);
            Point pp = GetPoint(*ts);
            LGMdouble ver_dist = H - pp.getZ();
            Axis<ETCfSegment,ETCfBud>& ax = GetAxis(t);
            TreeCompartment<ETCfSegment,ETCfBud>* fc = GetFirstTreeCompartment(ax);
            LGMdouble hor_dist = sqrt(pow(pp.getX()-GetPoint(*fc).getX(),2.0) +
                                      pow(pp.getY()-GetPoint(*fc).getY(),2.0));
	    LGMdouble Vf = GetValue(*ts,LGAVf);
	    LGMdouble fol_den = Af/Vf;
	    LGMdouble star = star_mean(runs,Rf,L,fol_den, ran3_seed);
	    LGMdouble s_sum = 0.0, c_sum = 0.0;
	    LGMdouble Sf = GetValue(*ts,LGAsf);
	    LGMdouble Wf = GetValue(*ts, LGAWf);
	    for(int k = 0; k < 11 ; k++) {
	      double theta = PI_VALUE * (double)k /(2.0 * 10.0);
	      s_sum += cos(theta)*S(theta,Sf,Wf,Rf,L);
	      c_sum += cos(theta);
	    }
	    LGMdouble star_eq = s_sum/c_sum;

        
            ofstream f(fname.c_str() , ofstream::app);
            f << pp.getZ() << " " << ver_dist << " " << hor_dist << " " << Age <<  " " << 100.0*Rf
	      << " " << 100.0*L << " " << 10000.0*Af << " " << 2000.0*Wf << " " 
	      << 1000000.0*Vf << " " << fol_den << " " << star << " " << star_eq << endl;
            f.close();
          }
        }
        return tc;
      }

 private:
    string fname;
    int runs;
};

#endif
