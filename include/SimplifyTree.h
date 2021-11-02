#ifndef SIMPLIFYTREE_H
#define SIMPLIFYTREE_H

#include <EvaluateTree.h>


class SimplifyTree {
 public:
  bool&  operator() (bool& seq_removed, TreeCompartment<ETCfSegment,ETCfBud>* tc)const {
    if(seq_removed) {
      return seq_removed;
    }

    if(Axis<ETCfSegment,ETCfBud>* ax = dynamic_cast<Axis<ETCfSegment,ETCfBud>*>(tc)) {
      list<TreeCompartment<ETCfSegment,ETCfBud>*>& cmpls = GetTreeCompartmentList(*ax);

      if(static_cast<int>(cmpls.size()) < 5) {  //at least two segments (+ two BPs and 1 Bud) 
	seq_removed = false;
	return seq_removed;
      }

      seq_removed = false;
      std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator I_last;
      TreeSegment<ETCfSegment, ETCfBud>* last_seg = GetLastTreeSegment(*ax);
      std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator I_first = cmpls.begin();
	    
      while(I_first != cmpls.end() && *I_first != last_seg && !seq_removed) {
	int n_empty = 0;
	if(ETCfSegment* ts = dynamic_cast<ETCfSegment*>(*I_first)) {
	  I_last = I_first;
	  I_last++;           //As *I_first = TS, I_last cannot be past end
	  n_empty = 0;
	  bool find_end2 = false;

	  while((I_last != cmpls.end()) && (!find_end2)) {
	    if(BranchingPoint<ETCfSegment,ETCfBud>* bp =
	       dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*I_last)) {
	      std::list<Axis<ETCfSegment,ETCfBud>*>& axl = GetAxisList(*bp);
	      if(static_cast<int>(axl.size()) == 0) {
		n_empty++;
	      }
	      else { 
		find_end2 = true;
	      }
	    } // if(BranchingPoint<ETCfSegment,ETCfBud>* bp = ...

	    if(!find_end2) {
	      I_last++;
	    }
	  } // while((I_last != cmpls.end()) && (!find_end2)) ...

	  if(n_empty > 0) {
	    if(I_last == cmpls.end()) {
	      I_last--; I_last--; I_last--;   //while ended with I_last = cmpls.end()
	    } else {
	      I_last--;                       //while ended with find_end2 == true
	    }
	    seq_removed = true;
	    double L = GetPoint(**I_first) || GetEndPoint(*(dynamic_cast<ETCfSegment*>(*I_last)));
	    SetValue(*(dynamic_cast<ETCfSegment*>(*I_first)), LGAL, L);
	    SetValue(*(dynamic_cast<ETCfSegment*>(*I_first)), LGARTop,
		     GetValue(*(dynamic_cast<ETCfSegment*>(*I_last)),LGAR));
	    SetDirection(**I_first, PositionVector(GetEndPoint(*(dynamic_cast<ETCfSegment*>(*I_last)))
						   -GetPoint(**I_first)));
	    std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator I;
	    std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator Istart = ++I_first;
	    std::list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator Iend = I_last++;
	    for(I = Istart; I != Iend; I++) {
	      delete *I;
	    } 
	    I = cmpls.erase(Istart,I_last);
	    if(seq_removed)  break;        //break from the outer while loop and exit
	  }
	} //if(ETCfSegment* ts = ...
	I_first++;
      }  //  while(I_first != cmpls.end() && ...
    }
    return seq_removed;
  }

};

#endif
