#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <vector>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
#include <XMLTree.h>
#include <EvaluateTree.h>

class a_info {
 public:
 a_info() : n_fork(0), orders(0), Qvs(0) {;}
  a_info& operator += (const a_info& add_this) {
    n_fork++;
    for(int i = 0; i < (int)(add_this.orders).size(); i++) {
      orders.push_back((add_this.orders)[i]);
    }
    for(int i = 0; i < (add_this.Qvs).size(); i++) {
      Qvs.push_back((add_this.Qvs)[i]);
    }

  }
  int n_fork;
  vector<double> orders;
  vector<double> Qvs;
};

class basi {
 public:
 basi(const double la_in) : la(la_in) {;}
  a_info& operator() (a_info& val, TreeCompartment<ETCfSegment,ETCfBud>* tc) const {
    vector<double>::iterator Io;
    vector<double>::iterator Iq;

    if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(tc)) {

      if((int)GetValue(*ts,LGAomega) == 1) {
	Iq = (val.Qvs).begin();
	*Iq += 1.0;
      }
      else {
	Iq = (val.Qvs).begin();
	*Iq += 0.25;
      }
      cout << "TS my order " << (int)GetValue(*ts,LGAomega) << endl;
      cout << "TS Orders" << endl;
      for(Io = (val.orders).begin(); Io != (val.orders).end(); Io++) {
	cout << (*Io) << " ";
      }
      cout << endl;
      cout << "TS Qvs" << endl;
      for(Iq = (val.Qvs).begin(); Iq != (val.Qvs).end(); Iq++) {
	cout << (*Iq) << " ";
      }
      cout << endl;

    }
    
    if(BranchingPoint<ETCfSegment,ETCfBud>* bp = dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(tc)) {
/*       cout << "my_order " << GetValue(*bp, LGAomega) << endl; */
/*       cout << "order Qv " << val.order << " " << val.Qv << endl; */
      cout << "BP my order " << (int)GetValue(*bp,LGAomega) << endl;
      cout << "BP Orders" << endl;
      for(Io = (val.orders).begin(); Io != (val.orders).end(); Io++) {
	cout << (*Io) << " ";
      }
      cout << endl;
      cout << "BP Qvs" << endl;
      for(Iq = (val.Qvs).begin(); Iq != (val.Qvs).end(); Iq++) {
	cout << (*Iq) << " ";
      }
      cout << endl;

      double sum2 = 0.0;
      double my_order = GetValue(*bp,LGAomega);
      for(int i = 0; i < (int)(val.orders).size(); i++) {
	if((val.orders)[i] > my_order) {
	  sum2 += (1.0 - la) * (val.Qvs)[i];
	    } else {
	  sum2 += la * (val.Qvs)[i];
	}
      }
      SetValue(*bp, LGAMaxD,sum2);
      cout << "BP I store " << sum2 << " " << GetValue(*bp,LGAMaxD) << endl;

      double sum = accumulate((val.Qvs).begin(),(val.Qvs).end(),0.0);
      (val.orders).clear();
      (val.Qvs).clear();
      (val.orders).push_back(GetValue(*bp,LGAomega));
      (val.Qvs).push_back(sum);

      cout << "BP I send; size, order, sum " << (int)(val.orders).size() << " " << GetValue(*bp,LGAomega)
	   <<  " " << sum << endl;
    }


    return val;
  }

 private:
  double la;
};


// bool streit for writing file in the format Katarina Streit sent tree data
// An example can be seen in /Users/matrps/Riston-D/E/Hankkeet/LIGNUM
// /Erillishankkeet/ERC-yhteistyo/2016/Henri-Riihimaki/Streit-trees/rauh/rauh0.txt
// Tiedostossa on yhden sylinterin (segmentin) tiedot yhdella rivilla, rivi:
// The following data on each line:
// 1.   number 
// 2.   radius
// 3.   length 
// 4.   start_x
// 5.   start_y
// 6.   start_z
// 7.   direction_x
// 8.   direction_y
// 9.   direction_z
// 10.  number_of_mother_internode
// 11   order

template<class TS, class BUD>
void writeCylinderData(int& line_number, int parent_no, int& my_branchno, 
		       Axis<TS,BUD>& axis, bool streit = false) {

  list<TreeCompartment<TS,BUD>*>& tc_ls =
    GetTreeCompartmentList(axis);
  if((int)tc_ls.size() < 3) {  //No TreeSegments - nothing to do
    return;
  }

  //1. Order of Axis
  TreeSegment<TS,BUD>* ts  = GetFirstTreeSegment(axis);
  int order = (int)GetValue(*ts,LGAomega);
  if(!streit) {
    order--;                        //Stem = order 0 except streit
  }

  // 2) Number of segments in the axis
  typename list<TreeCompartment<TS,BUD>*>::iterator Ic;
  int n_segs = 0;
  for(Ic = tc_ls.begin(); Ic != tc_ls.end(); Ic++) {
    if(TS* ts = dynamic_cast<TS*>(*Ic)){
      n_segs++;
    }
  }

  // 3) Go through the axis, write if TreeSegment

  int line_number_of_first_segment = line_number + 1;  //Store for later use
  int seg_no_in_axis = 0;
  for(Ic = tc_ls.begin(); Ic != tc_ls.end(); Ic++) {
    if(TS* ts = dynamic_cast<TS*>(*Ic)){
      seg_no_in_axis++;
      line_number++;
      int next = line_number + 1;
      if(seg_no_in_axis == n_segs) {
	next = 0;
      }
      int parent = 0;
      if(seg_no_in_axis == 1) {
	parent = parent_no;
      }
      else {
	parent = line_number - 1;
      }
      LGMdouble r = 1000.0*GetValue(*ts, LGAR);
      LGMdouble l = 1000.0*GetValue(*ts, LGAL);
      Point p = GetPoint(*ts);
      PositionVector d = GetDirection(*ts);
      LGMdouble rh = 1000.0*GetValue(*ts, LGARh);
      if(!streit) {
	cout << r << " " << l << " "
	     << p.getX() << " " << p.getY() << " " << p.getZ() << " "
	     << d.getX() << " " << d.getY() << " " << d.getZ() << " "
	     << parent << " " << next << " " << my_branchno << " "
	     << order << " " << seg_no_in_axis << " " << 0 << " " << rh << endl;
      } else {
	cout << line_number - 1 << " " << r/1000.0 << " " << l/1000.0 << " "
	     << p.getX() << " " << p.getY() << " " << p.getZ() << " "
	     << d.getX() << " " << d.getY() << " " << d.getZ() << " "
	     << parent -1 << " " <<  (int)order  << endl;
      }
    }
  }

  // 3) Process the BranchingPoints of this Axis
//      int new_branchno = my_branchno;

  seg_no_in_axis = 0;
  for(Ic = tc_ls.begin(); Ic != tc_ls.end(); Ic++) {
    if(TS* ts = dynamic_cast<TS*>(*Ic)){
      seg_no_in_axis++;
    }

    if(BranchingPoint<TS,BUD>* bp =
       dynamic_cast<BranchingPoint<TS,BUD>*>(*Ic)) {
      list<Axis<TS,BUD>*>& ax_lst = GetAxisList(*bp);
      int parent_no = line_number_of_first_segment + seg_no_in_axis - 1;
      typename list<Axis<TS,BUD>*>::iterator ai;
      for(ai = ax_lst.begin(); ai != ax_lst.end(); ai++) {
	list<TreeCompartment<TS,BUD>*>& tc_ls2 =
	  GetTreeCompartmentList(**ai);
	if((int)tc_ls2.size() >= 3) {
	  my_branchno++;
	  writeCylinderData(line_number, parent_no, my_branchno, **ai, streit);
	}
      }   //for(ai = ax_lst ...
    }     //if(BranchingPoint<TS,BUD>* bp = ...
    
  }
  return;
}


template<class TS, class BUD>
void writeBranchData(int parent_no, int& my_branchno, PositionVector parent_dir,
		     Axis<TS,BUD>& axis) {

  list<TreeCompartment<TS,BUD>*>& tc_ls =
    GetTreeCompartmentList(axis);
  if((int)tc_ls.size() < 3) {  //No TreeSegments - nothing to do
    return;
  }

  //1. Order of Axis
  TreeSegment<TS,BUD>* ts  = GetFirstTreeSegment(axis);
  int order = (int)GetValue(*ts,LGAomega);
  order--;                        //Stem = order 0


  // 3) Go through the axis

  LGMdouble vol = 0.0;
  LGMdouble length = 0.0;

  typename list<TreeCompartment<TS,BUD>*>::iterator Ic;
  for(Ic = tc_ls.begin(); Ic != tc_ls.end(); Ic++) {
    if(TS* ts = dynamic_cast<TS*>(*Ic)){
      vol += GetValue(*ts, LGAV);
      length += GetValue(*ts, LGAL);
    }
  }

  TreeSegment<TS,BUD>* first = GetFirstTreeSegment(axis);
  TreeSegment<TS,BUD>* last  = GetLastTreeSegment(axis);

  Point start_p = GetPoint(*first);
  Point end_p = GetEndPoint(*last);
  PositionVector dir = GetDirection(*first);
  LGMdouble angle = acos(Dot(dir,parent_dir))*180.0/PI_VALUE;
  LGMdouble r = GetValue(*first,LGAR);

  cout << order << " " << parent_no << " " << 1000.0*vol << " " << length << " "
       << 1000.0*r << " " << angle << " " 
       << start_p.getX() << " " << start_p.getY() << " " << start_p.getZ() << " "
       << end_p.getX() << " " << end_p.getY() << " " << end_p.getZ() << endl;

  int local_branchno = my_branchno;

  // 3) Process the BranchingPoints of this Axis

  for(Ic = tc_ls.begin(); Ic != tc_ls.end(); Ic++) {
    if(TS* ts = dynamic_cast<TS*>(*Ic)){
      parent_dir = GetDirection(*ts);
    }

    int new_parent_no = local_branchno;

    if(BranchingPoint<TS,BUD>* bp =
       dynamic_cast<BranchingPoint<TS,BUD>*>(*Ic)) {
      list<Axis<TS,BUD>*>& ax_lst = GetAxisList(*bp);
      typename list<Axis<TS,BUD>*>::iterator ai;
      for(ai = ax_lst.begin(); ai != ax_lst.end(); ai++) {
	list<TreeCompartment<TS,BUD>*>& tc_ls2 =
	  GetTreeCompartmentList(**ai);
	if((int)tc_ls2.size() >= 3) {
	  my_branchno++;
	  writeBranchData(new_parent_no, my_branchno, parent_dir, **ai);
	}
      }   //for(ai = ax_lst ...
    }     //if(BranchingPoint<TS,BUD>* bp = ...
    
  }

  return;
}

//template <class TS, class BUD>
  void simplify_tree(Axis<ETCfSegment,ETCfBud>& axis);

#endif
