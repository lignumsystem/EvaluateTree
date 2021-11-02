#ifndef SAPWOODFUNCTORS_H
#define SAPWOODFUNCTORS_H

// Usage: AccumulateDown(tree,PSAD_Info,PartSapwoodAreaDown(vector<double> frac),SetRh())


class PSAD_Info {
 public:
 PSAD_Info(const double sw) : order(0), sw_area(sw) {}
  int order;
  double sw_area;
};

class PartSapwoodAreaDown{
public:
 PartSapwoodAreaDown(const vector<double> frac) : fraction(frac) {}
 PartSapwoodAreaDown(const PartSapwoodAreaDown& swdown):fraction(swdown.fraction){}
 PSAD_Info& operator()(PSAD_Info& ps1, PSAD_Info& ps2)const
  {
    int o1 = ps1.order;
    int o2 = ps2.order;
    if (o1 < o2){
      double percent = fraction[o2-1];  //vector<double> indexing starts from 0
      double As = ps2.sw_area;
      double Asdown = percent*As;
      ps1.sw_area += Asdown;
    }
    else{
      ps1.sw_area += ps2.sw_area;
    }
    return ps1;
  }
private:
 const vector<double>fraction;
};


class SetRh{
public:
 SetRh(const bool cr) : correct_r(cr) {}
  PSAD_Info& operator()(PSAD_Info& data,TreeCompartment<ETCfSegment,ETCfBud>* tc)const
  {
      if (ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
	//sapwood demand of own foliage
	//sapwood area of new shoots in the tree of Nikinmaa et al 2014 (7d.xml)
	//is relation Wf(kg C) = 2.895752e+02 x (sapwood area (m2)), thus
	double own_swa = GetValue(*ts, LGAWf) / 2.895752e+02;
	double swa = own_swa + data.sw_area;            //+ demand from above

	data.sw_area = swa;                            //this goes downwards
	data.order = static_cast<int>(GetValue(*ts,LGAomega));

	double R = GetValue(*ts, LGAR);
	if(swa >= PI_VALUE*R*R) {
	  SetValue(*ts, LGARh, 0.0);
	  if(correct_r) {
	    SetValue(*ts,LGAR, sqrt(swa/PI_VALUE));
	  }
	} else {
	  SetValue(*ts, LGARh, sqrt((PI_VALUE*R*R-swa)/PI_VALUE));
	}
	
      }
      return data;
  }
 private:
  bool correct_r;
};

class PrintStemDiameters{
public:
  int& operator()(int& dummy,TreeCompartment<ETCfSegment,ETCfBud>* tc)const
  {
      if (ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
	if(static_cast<int>(GetValue(*ts,LGAomega)) == 1) {
	  double Z = GetMidPoint(*ts).getZ();
	  double R = GetValue(*ts, LGAR);
	  double Rh = GetValue(*ts,LGARh);
	  cout << Z << " " << R << " " << Rh << endl;
	}
      }
      return dummy;
  }
};


#endif

