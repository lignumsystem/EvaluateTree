#ifndef EVALUATETREE_H
#define EVALUATETREE_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
#include <XMLTree.h>

// ET = EvaluateTree

class ETCfBud;
class ETHwBud;

class ETCfSegment : public CfTreeSegment<ETCfSegment,ETCfBud>{
 public:
 ETCfSegment(const Point& p,const PositionVector& pv,
                     const LGMdouble go,const METER l, const METER r,
                     const METER rn, Tree<ETCfSegment,ETCfBud>* tree)
   :CfTreeSegment<ETCfSegment,ETCfBud>(p,pv,go,l,r,rn,tree)
    {}

  int getNumber() {return number;}
  void setNumber(const int& nu) {number = nu;}

 private:
  int number;
};

class ETHwSegment : public HwTreeSegment<ETHwSegment,ETHwBud,Triangle>{
 public:
 ETHwSegment(const Point& p,const PositionVector& pv,
                     const LGMdouble go,const METER l, const METER r,
                     const METER rn, Tree<ETHwSegment,ETHwBud>* tree)
   :HwTreeSegment<ETHwSegment,ETHwBud,Triangle>(p,pv,go,l,r,rn,tree)
    {}
};

class ETCfBud:public Bud<ETCfSegment,ETCfBud>{
 public:
 ETCfBud(const Point& p, const PositionVector& d, 
             const LGMdouble go, Tree<ETCfSegment,ETCfBud>* tree)
   :Bud<ETCfSegment,ETCfBud>(p,d,go,tree){}
};

class ETHwBud:public Bud<ETHwSegment,ETHwBud>{
 public:
 ETHwBud(const Point& p, const PositionVector& d, 
             const LGMdouble go, Tree<ETHwSegment,ETHwBud>* tree)
   :Bud<ETHwSegment,ETHwBud>(p,d,go,tree){}
};


#endif
