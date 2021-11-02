#ifndef EVALUATEFUNCTORS_H
#define EVALUATEFUNCTORS_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
#include <XMLTree.h>

#include <EvaluateTree.h>

// ET = EvaluateTree

/* class ETCfBud; */
/* class ETHwBud; */


class CheckNeedlesLast {
 public:
  bool&  operator()
    (bool& last, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      if(last) {
	last = false;
	if(GetValue(*ts,LGAWf) < R_EPSILON) {
	  Point p = GetPoint(*ts);
	  LGMdouble r = GetValue(*ts, LGAR);
	  LGMdouble o = GetValue(*ts, LGAomega);
	  cout << p.getX() << " " << p.getY() << " " << p.getZ() << " "
	       << o << " " << r << endl;
	}
      }
    }
    return last;
  }
};


class MultWf10 {
 public:
  TreeCompartment<ETCfSegment,ETCfBud>* operator()
    (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      LGMdouble Wf = GetValue(*ts, LGAWf);
      if(Wf > R_EPSILON) {
	SetValue(*ts,LGAWf, 10.0*Wf);
      }
    }
    return tc;
  }
};


class ForwardGO {
 public:
  double&  operator()
    (double& last, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      SetValue(*ts, LGAomega, last);
    }
    if(Axis<ETCfSegment,ETCfBud>* ax = dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(tc)) {
      last += 1.0;
    }
    return last;
  }
};

class SetRf {
 public:
 SetRf(const double& setrf) : set_value(setrf) {}
  TreeCompartment<ETCfSegment,ETCfBud>*  operator()
    (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      if(GetValue(*ts, LGAWf) > R_EPSILON) {
	SetValue(*ts, LGARf, set_value);
      }
    }
    return tc;
  }
 private:
  double set_value;
};


class SetWf100 {
 public:
  TreeCompartment<ETCfSegment,ETCfBud>*  operator()
    (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      if(GetValue(*ts, LGAWf) > R_EPSILON) {
	SetValue(*ts, LGAWf, 1.3*GetValue(*ts, LGAWf));
      }
    }
    return tc;
  }
};


class AdjustWf {
 public:
 AdjustWf(const double& h) : treeh(h) {}
  TreeCompartment<ETCfSegment,ETCfBud>*  operator()
    (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      if(GetValue(*ts, LGAWf) > R_EPSILON) {
	SetValue(*ts, LGAWf, (1.0/(4.0-3.7*treeh/20.0))*GetValue(*ts, LGAWf));
      }
    }
    return tc;
  }
 private:
  double treeh;
};

class AdjustSf {
 public:
  TreeCompartment<ETCfSegment,ETCfBud>*  operator()
    (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      if(GetValue(*ts, LGAWf) > R_EPSILON) {
	SetValue(*ts, LGAsf, 1.15*GetValue(*ts, LGAsf));
      }
    }
    return tc;
  }
};


class SetSegmentNumber {
 public:
  int&  operator()
    (int& number, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      ts->setNumber(number);
      number++;
    }
    return number;
  }
};


class WriteCylinderDataStreit {
 public:
  int&  operator()
    (int& mother_number, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      LGMdouble r = GetValue(*ts, LGAR);
      LGMdouble l = GetValue(*ts, LGAL);
      Point p = GetPoint(*ts);
      PositionVector d = GetDirection(*ts);
      double order = GetValue(*ts, LGAomega);
      int my_number = ts->getNumber();
      cout << my_number << " " << r << " " << l << " "
	   << p.getX() << " " << p.getY() << " " << p.getZ() << " "
	   << d.getX() << " " << d.getY() << " " << d.getZ() << " "
	   << mother_number << " " <<  (int)order  << endl;
      mother_number = my_number;
    }
    return mother_number;
  }
};

//Adds a segment of length separation before the last segment in axis. This segment does not
//have foliage.
class SeparateTopShoots {
 public:
 SeparateTopShoots(const double sep, Tree<ETCfSegment,ETCfBud>& tt, const double s_sep = 0.0 ) : separation(sep), 
    my_tree(tt), second_separation(s_sep) {}
  TreeCompartment<ETCfSegment,ETCfBud>* operator() (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(Axis<ETCfSegment,ETCfBud>* ax = dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(tc)) {
      list<TreeCompartment<ETCfSegment,ETCfBud>*>& cmpls = GetTreeCompartmentList(*ax);
      int cmpls_size = static_cast<int>(cmpls.size());
      if(cmpls_size > 2) {  //at least one segment
	TreeSegment<ETCfSegment,ETCfBud>* last_seg = GetLastTreeSegment(*ax);
	if(last_seg == NULL) {  //maybe unnecessary
	  return tc;
	}
	
	int ordr =  static_cast<int>(GetValue(*last_seg,LGAomega));
	if(ordr < 2) {                    //not main stem
	  return tc;
	}
	if(ordr == 2 && cmpls_size < 6) { //Don't touch top branches
	  return tc;
	}

		//cout << "ORDER cmpls_size " << ordr << " " << cmpls_size << endl;

	// Add now new segment before the last one
	std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator I = cmpls.end();  //Points past last item
	I--; //Terminating bud
	I--; //Branching point
	I--; //Last segment
		//cout << "HEP" << endl;

	//	Tree<ETCfSegment,ETCfBud>my_tree =  GetTree(**I);
	Point base = GetPoint(**I);
		//cout << "HEP1" << endl;
		//cout << "X " << base;
	if(ETCfSegment* ee = dynamic_cast<ETCfSegment*>(*I)) {

	  	  //cout << "Isohep" << endl;
	}
	double len = GetValue(*(dynamic_cast<ETCfSegment*>(*I)),LGAL);
		//cout << "HEP2" << endl;

	PositionVector hdng = GetDirection(**I);
		//cout << "HEP3" << endl;

	Point top = base + (Point)(len*hdng);
	double rad = GetValue(*(dynamic_cast<ETCfSegment*>(*I)),LGAR);
		//cout << "HEP4" << endl;

	Point move = (Point)(separation*hdng);
	
	SetPoint(**I, base + move);       //Segment moved
	I++;
	SetPoint(**I, top + move);    //BranchingPoint moved
	I++;
	SetPoint(**I, top + move);    //Bud moved
	I--;   I--;
		//cout << "HEP" << endl;

	ETCfSegment* ts = new ETCfSegment(base, hdng, ordr, separation, rad, 0.0, &my_tree);
	SetValue(*ts, LGAWf, 0.0);
	SetValue(*ts, LGARf, rad);
	SetValue(*ts, LGAsf, 32.0);
	BranchingPoint<ETCfSegment,ETCfBud>* bp =
	  new BranchingPoint<ETCfSegment,ETCfBud>(base + move, hdng, ordr, &my_tree);
		//cout << "HEP" << endl;

	I = cmpls.insert(I, bp);
	I = cmpls.insert(I,ts);

/* 	std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator ss; */
/* 	int ii=1; */
/* 	for(ss=cmpls.begin(); ss != cmpls.end(); ss++) { */

	  //cout << ii << endl;
/* 	  if(BranchingPoint<ETCfSegment,ETCfBud>* bbb =dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*ss)) { */
/* 	    	    //cout << "On BP" << endl; */
/* 	  } */
/* 	  if(TreeSegment<ETCfSegment,ETCfBud>* ttt =dynamic_cast<TreeSegment<ETCfSegment,ETCfBud>*>(*ss)) { */
/* 	    	    //cout << "On TS" << endl; */
/* 	    	    //cout << GetValue(*ttt, LGAR) << endl; */
/* 	  } */
/* 	  if(Bud<ETCfSegment,ETCfBud>* uuu = dynamic_cast<Bud<ETCfSegment,ETCfBud>*>(*ss)) { */
/* 	    	    //cout << "On BUD" << endl; */
/* 	  } */
/* 	  if(Axis<ETCfSegment,ETCfBud>* aaa = dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(*ss)) { */
/* 	    	    //cout << "On BUD" << endl; */
/* 	  } */
 
	  //	  	  if(ss == I) //cout << "THIS" << endl;

	  //	  ii++;
	  //	}

	//	exit(0);
	//cout << "HEP" << endl;
/* 	  if(TreeSegment<ETCfSegment,ETCfBud>* ttt =dynamic_cast<TreeSegment<ETCfSegment,ETCfBud>*>(*I)) { */
/* 	    //cout << "On TS" << endl; */
/* 	    //cout << "L " << GetValue(*(dynamic_cast<TreeSegment<ETCfSegment,ETCfBud>*>(*I)), LGAL) << endl; */
/* 	  } */
	  //cout << (int)cmpls.size() << endl;
	  //if(I == cmpls.end())
	  	    //cout << "End" << endl;
	//		    cout << second_separation << " " << static_cast<int>(cmpls.size()) << endl;
	//		    cout << second_separation << " " << cmpls_size << endl;

		    //	  if(second_separation > 0.0 & static_cast<int>(cmpls.size()) > 6) {
	  if(second_separation > 0.0 &  cmpls_size > 6) {

	    //	    	cout << "HEP1 " <<   cmpls_size << endl;

	  I--;              //Points to BranchingPoint
	  //	  if(*I == NULL) {
	    //cout << "???" << endl;
	  //	  }
	  //	  Add second Segment only if there are sidebranches
	  	//cout << "HEP2 " <<   cmpls_size << endl;

	  //	  if(*I == NULL) {
	    //cout << "???" << endl;
	  //	  }
/* 	  if(BranchingPoint<ETCfSegment,ETCfBud>* bbb =dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*I)) { */
/* 	    //cout << "On BP" << endl; */
/* 	  } */
/* 	  if(TreeSegment<ETCfSegment,ETCfBud>* ttt =dynamic_cast<TreeSegment<ETCfSegment,ETCfBud>*>(*I)) { */
/* 	    //cout << "On TS" << endl; */
/* 	  } */
/* 	  if(Bud<ETCfSegment,ETCfBud>* uuu = dynamic_cast<Bud<ETCfSegment,ETCfBud>*>(*I)) { */
/* 	    //cout << "On BUD" << endl; */
/* 	  } */
/* 	  if(Axis<ETCfSegment,ETCfBud>* aaa = dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(*I)) { */
/* 	    //cout << "On BUD" << endl; */
/* 	  } */


	  std::list<Axis<ETCfSegment,ETCfBud>*>& ax_lst =
	    GetAxisList(*(dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*I)));
	  	  //cout << "Axlst  " << static_cast<int>(ax_lst.size()) << endl;
	  	//cout << "HEP3 " <<   cmpls_size << endl;

	  if(static_cast<int>(ax_lst.size()) > 0) {
	    	cout << "HEP4 " <<   static_cast<int>(ax_lst.size()) << endl;

	    I--;              //Points to TreeSegment
	    hdng = GetDirection(**I);
	    rad = GetValue(*(dynamic_cast<ETCfSegment*>(*I)),LGAR);
	    ts = new ETCfSegment(base, hdng, ordr, second_separation, rad, 0.0, &my_tree);
	    SetValue(*ts, LGAWf, 0.0);
	    SetValue(*ts, LGARf, rad);
	    SetValue(*ts, LGAsf, 32.0);
	    bp = new BranchingPoint<ETCfSegment,ETCfBud>(base, hdng, ordr, &my_tree);
	    I++;               //Points to BranchingPoint
	    I = cmpls.insert(I,ts);
	    I = cmpls.insert(I, bp);
	    cout << "0000000000000000  Added! " << endl;
	  }
	  }
      }
    }

    return tc;
  }
 private:
  double separation;
  Tree<ETCfSegment,ETCfBud>& my_tree;
  double second_separation;
};


class ChangeOrder {
 public:
 ChangeOrder() : d_order(0) {}
 ChangeOrder(const double d_o) : d_order(d_o) {}

  TreeCompartment<ETCfSegment,ETCfBud>*  operator()
    (TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)){
      double my_order = GetValue(*ts, LGAomega);
      SetValue(*ts, LGAomega, my_order + d_order);
      return tc;
    }
    if(BranchingPoint<ETCfSegment,ETCfBud>* bp =
       dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(tc)) {
      double my_order = GetValue(*bp, LGAomega);
      SetValue(*bp, LGAomega, my_order + d_order);
      return tc;
    }
    if(Bud<ETCfSegment,ETCfBud>* bud =
       dynamic_cast<Bud<ETCfSegment,ETCfBud>*>(tc)) {
      double my_order = GetValue(*bud, LGAomega);
      SetValue(*bud, LGAomega, my_order + d_order);
      return tc;
    }
    return tc;
  }
 private:
  double d_order;
};



class SetWoodyZero {
 public:
 SetWoodyZero() : Hc(R_HUGE) {}
 SetWoodyZero(const double hc) : Hc(hc) {}

  bool&  operator()
    (bool& set, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(BranchingPoint<ETCfSegment,ETCfBud>* bp =
       dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(tc)) {
      if(GetValue(*bp,LGAomega) < 2.0 && GetPoint(*bp).getZ() >= Hc ) {
	set = false;
      }
    }
    if( ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc) ) {
      if(set) {
	LGMdouble Wf = GetValue(*ts, LGAWf);
	if(Wf < R_EPSILON) {
	  SetValue(*ts,LGAR, 0.0);
	  SetValue(*ts,LGARh, 0.0);
	  SetValue(*ts,LGARf, 0.0);
	}
      }
    }
    return set;
  }
 private:
  double Hc;
};


class MaxCrownWidthInfo {
 public:
 MaxCrownWidthInfo(const double dist0, const Point point0): dist(dist0),
    p0(point0) {}
  double dist;
  Point p0;
};

class MaxCrownWidth {
 public:
 MaxCrownWidth() : foliage(false) {}
 MaxCrownWidth(const bool fi) : foliage(fi) {}

  MaxCrownWidthInfo&  operator()
    (MaxCrownWidthInfo& info, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(BranchingPoint<ETCfSegment,ETCfBud>* bp =
       dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(tc)) {
      if(GetValue(*bp,LGAomega) == 1.0) {
        info.dist = 0.0;    //distance starts from start of branch
	info.p0 = GetPoint(*bp); //start of branch
      }
    }
    if( ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc) ) {
      if(GetValue(*ts, LGAomega) > 1.0 ) {
	  if( !foliage || (GetValue(*ts, LGAWf) > R_EPSILON) ) {
	    Point pb = GetPoint(*ts);
	    Point pe = GetEndPoint(*ts);
	    double db = sqrt(pow(pb.getX()-(info.p0).getX(),2.0)+pow(pb.getY()-(info.p0).getY(),2.0));
	    double de = sqrt(pow(pe.getX()-(info.p0).getX(),2.0)+pow(pe.getY()-(info.p0).getY(),2.0));
	    if( db > info.dist ) info.dist = db;
	    if( de > info.dist ) info.dist = db;
	  }
	}
    }

    if( ETCfBud* b = dynamic_cast<ETCfBud*>(tc) ) {
     Tree<ETCfSegment,ETCfBud>& my_tree =  GetTree(*dynamic_cast<TreeCompartment<ETCfSegment,ETCfBud>*>(b));
     if(GetValue(my_tree, LGPyc) < info.dist ) {  //Using Tree parameter LGPyc here is a hack
       SetValue(my_tree, LGPyc, info.dist);
     }
    }
    return info;
  }
 private:
  bool foliage;
};


class PrintLastSegInfo {
 public:
  TreeCompartment<ETCfSegment,ETCfBud>*
    operator() (TreeCompartment<ETCfSegment,ETCfBud>* tc) const {
    if( Axis<ETCfSegment,ETCfBud> *ax = dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(tc) ) {
        TreeSegment<ETCfSegment,ETCfBud> *ls1 =  GetLastTreeSegment(*ax);
	if(!(ls1 == NULL)) {
	  ETCfSegment *ls = dynamic_cast<ETCfSegment*>(ls1);
	  Point p = GetPoint(*ls);
	  PositionVector d = GetDirection(*ls);
	  cout << p.getX() << " " << p.getY() << " " << p.getZ() << " " << d.getX() << " " << d.getY()
	       << " " << d.getZ() << " " << static_cast<int>(GetValue(*ls, LGAomega)) << " " 
	       << GetValue(*ls, LGAR) << " " << GetValue(*ls, LGARh) << " " << GetValue(*ls, LGARf)
	       << " " << GetValue(*ls, LGAL) << " " << GetValue(*ls, LGAWf)  << " "
	       << GetValue(*ls, LGAAf) << endl;
	}
    }
    return tc;
  }
}; 



#endif
