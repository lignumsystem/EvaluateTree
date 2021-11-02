#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <Lignum.h>
#include <Bisection.h>
#include <Shading.h>
#include <VoxelSpace.h>

#include <XMLTree.h>

#include <EvaluateTree.h>
#include <EvaluateFunctors.h>
#include <FractalDimension.h>
#include <Lacunarity.h>
#include <SetShootPositions.h>
#include <LGMVisualization.h>
#include <Functions.h>
#include <Radiation.h>
#include <SimplifyTree.h>
#include <SapwoodFunctors.h>

int ran3_seed;

void Usage()
{
  cout << endl;
  cout << "This program evaluates a number of things from a tree that is read in from the xml file." << endl;
  cout << "If input tree is a deciduous one, Triangle is the leaf shape. Analyzing trees with other leaf" << endl;
  cout << "shapes is not possible (without further programming) if foliage characteristcs are needed in the" << endl;
  cout << "analysis." << endl;
  cout << endl;
  cout << "Usage: ./evaluate  <tree.xml> [-deciduous] [-mantleArea <value> [-foliage] [-crown] [-sum] ]" << endl;
  cout << "       [-foliageAreaMass [-byOrder] [-noTitle]]" << endl;
  cout << "       [-heights] [-noCompartments]" << endl;
  cout << "       [-fractalDimension <no_points> [-fractalDm <minbox>] [-fractalD0 <value>] ]"<< endl;
  cout << "       [-lacunarity <box> [-no_points <number>]]" << endl;
  cout << "       [-ellipsoid] [-randomizeShoots -Zc <value> -rh <value> -rv <value> [-seed <value>]" << endl;
  cout << "       [-woodyParts] [-bothEnds]]"<< endl;
  cout << "       [-treeSegmentInformation] [-needlesLast] [-multWf10]" << endl;
  cout << "       [-foliageMassEqs [-noTitle]]  [-raumonen] [-raumonenB]" << endl;
  cout << "       [-sideBranches]" << endl;
  cout << "       [-STARvalues [-runs <value>] [-starFile <name>]]" << endl;
  cout << "       [-makeTree [-nShoot <n>] [-Rf <value>] [-len <value>] [-fol_den <value>]]" << endl;
  cout << "       [-moveTree [-X <value>] [-Y <value>] [-Z <value>]]" << endl;
  cout << "       [-segmentLengths [-byOrder] [-noTitle]] [-setRf <value>]" << endl;
  cout << "       [-monster] [-adjustWf] [-adjustSf] [-streit] [-countTreeSegments] [-simplify]" << endl;
  cout << "       [-separateTopShoots <value> [-secondSep <value>] ]" "[-displayStructure [-byAxis]]"<< endl;
  cout << "       [-verticalLA [-step <value>] [-min <value>] [-max <value>] ]" << endl;
  cout << "       [-changeOrder <value>] [-setOrder [-iniVal <value>]] [-setWoodyZero [-Hc <value>] ]" << endl;
  cout << "       [-maxCrownWidth [-foliage]] [-collectSegmentLengths [-order <value>]]" << endl;
  cout << "       [-noBranchWhorls] [-meanBranchAngle] [-setSapwoodArea [-correctR]]" << endl;
  cout << "       [-printLastSegInfo] [-printStemDiameters]" << endl;
  cout << endl;
  cout << "Input tree is assumed to be a conifer, -deciduous changes it to deciduous." << endl;
  cout << "-mantleArea    evaluates surface area of trunk and branches without foliage by gravelius orders" << endl;
  cout << "               and writes to console. With -foliage also parts with foliage are included." << endl;
  cout << "               If -crown is set, only mantle area of stemsegments inside crown are evaluated," << endl;
  cout << "               segments order == value are evaluated, if -sum is set order >= value segments evaluated." << endl;
  cout << "-foliageAreaMass   Evaluates foliage area (total in conifers, one-sided in deciduous) and mass" << endl;
  cout << "                   and writes to console. If flag -byOrder is set writes also by Gravelius orders" << endl;
  cout << "                   and total. If flag -noTitle only values are written" << endl;
  cout << "-heights           Tree height and height of crown base (only conifer at the moment)" << endl;
  cout << "-noCompartments    Number of TreeCompartment in different classes to console." << endl;
  cout << "-vizTree           LIGNUM visualization of tree. Leaves of deciduous tree shown as triangles." << endl;
  cout << "-fractalDimension  writes to console no. of occupied cubes as a function of cube size." << endl; 
  cout << "              That is; box counting. Shoot cylinders are covered with points that are used to" << endl;
  cout << "              test whether a cube is occupied. <no_points> specifies how many to cover a cylinder." << endl; 
  cout << "         Only for pines at the moment." << endl;
  cout << "-fractalDm <minbox>    minimum box side lentgth given as part of original box size." << endl;
  cout << "         Default = 1000 meaning minimun side length is 1/1000 of original." << endl;
  cout << "-fractalD0 <value>     Starts d iteration from a given value = span/<value>" << endl;
  cout << "         span is written out when this program starts." << endl;
  cout << "-lacunarity <box>" << endl;
  cout << "         Only for pines at the moment." << endl;
  cout << "-CACReflect  If part of gliding box window goes outside of grid if Centered Allain-Cloitre calculation," << endl;
  cout << "             it is reflected to other side of the grid (periodical boundary)." << endl;
  cout << "-ellipsoid   The crown volume (between highest and lowest TreeSegment with foliage) is evaluated in" << endl;
  cout << "             20 cm high circular disks. The radius of each disk is maximum distance of foliage TreeSegment" << endl;
  cout << "             (farthest of the ends). The radius of a circular ellpsoid that is between crown base" << endl;
  cout << "             and tree top is calculated and the result is written to console (tree_file height of" << endl;
  cout << "             ellipse center horizontal radius vertical radius)." << endl;
  cout << "-randomizeShoots   Moves segments with foliage above crown base into random positions in ellipsoid with center," << endl;
  cout << "             horizontal, and vertical radii given by -Zc, -rh, and  -rv. If -woodyParts is specified, also" << endl;
  cout << "             segments without foliage are moved. If -bothEnds is set both ends of segment are checked for" << endl;
  cout << "             being in ellipsoid. Normally only midpoint is checked. Writes the ellipsoid tree into file" << endl;
  cout << "             input_file-ellisoid.xml." << endl;
  cout << "-treeSegmentInformation      Writes (all) information about TreeSegments to file treesegmentinformation.dat" << endl;
  cout << "             using functor ../stl-lignum/include/PrintTreeSegmentInformationToFile<TS,BUD>" << endl;
  cout << "-needlesLast Writes out information about last segments of axes that do not have foliage." << endl;
  cout << "-multWf10    Multply needle mass by 10 and save tree in <filename>-10.xml" << endl;
  cout << "-foliageMassEqs              Evaluates foliage biomass of the tree and claculates fol. biomass with the aid of" << endl;
  cout << "             dbh and height using several equations: Repola, Marklund, Wirth and Ilvesniemi." << endl;
  cout << "             The equations are taken from ../DigitalTree/cf-main.cc" << endl;
  cout << "             Writes to console, if -noTitle writes only values." << endl;
  cout << "-sideBranches   Writes to console information along the main Axis (meant for fine root analysis)" << endl;
  cout << "             Dist. from base, # segments, # apexes, total segment length" << endl;
  cout << endl;
  cout << "-STARvalues  Writes STAR values and other information of segments to a file. Works only for conifers";
  cout << endl;
  cout << "             If filename is not given in -starFile <name> writes to starvalues.dat." << endl;
  cout << "             STAR is analyzed using MC, default repetitions = 1000, can be changed by -runs <value>";
  cout << endl;
  cout << "-makeTree    Creates a tree (axis & treesegments) consisting of given number (-nShoot) of identical shoots."
       << endl;
  cout << "             Default values are used if parameters are not specified." << endl;
  cout << "             -Rf = radius (m) of shoot cylinder, -len = length(m), -fol_den = foliage density (m2/m3) in"
       << endl;
  cout << "             shoot cylinder (all-sided needle area), woody radius = 0.0001 m" << endl;
  cout << endl;
  cout << "-moveTree    Moves a tree from its (base) current position to one given by -X <value> -Y <value>" << endl;
  cout << "             -Z <value> (if omitted position = (0, 0, 0,)). Writes to file inputfile_moved.xml" << endl;
  cout << "-segmentLengths    Works like -foliageAreaMass but for lengths of TreeSegments." << endl;
  cout << "-setRf <value>     Sets foliage radius to value in segments that have foliage" << endl;
  cout << "-monster     Sets Wf of shoots (that have foliage) = 100 kg C, their transparency is thus = 0," << endl;
  cout << "             monster trees. Writes tree to file *_monster.xml" << endl;
  cout << "-adjustWf    Adjusts Wf (at the moment Wf *= 1.0/(4.0-3.7*treeh/20.0) ), writes to *_adjwf.xml file." << endl;
  cout << "-adjustSf    Adjusts sf (at the moment sf *= 1.15, writes to *_adjsf.xml file." << endl;
  cout << "-simplify    Creates output tree in which adjacent segments have been lumped together if they had not" << endl;
  cout << "             later branches between. Output tree contains only architectural information (lengths, " << endl;
  cout << "             radii, orientation)" << endl; 
  cout << "-separateTopShoots <value>  Creates a shoot of length separation (cm) before the top shoot." << endl;
  cout << "             Writes to <treefile>_separateXXcm.xml where XX = separation" << endl;
  cout << "-displayStructure   Writes the structure of the tree on console (by DisplayStructure of" << endl;
  cout << "             TreeFunctorI.h). If -byAxis writes contents of each Axis on own line (by " << endl;
  cout << "             DisplayStructureAxis of TreeFunctorI.h)" << endl;
  cout << "-verticalLA  Evaluates vertical needle area (all-sided)/leaf area distribution in 20 cm slices," << endl;
  cout << "             -step <value> changes this default." << endl;
  cout << "-changeOrder Changes order of compartments by <value>. No questions asked." << endl;
  cout << "             Then writes the tree to file treefile-neworder.xml" << endl;
  cout << "-setOrder    Sets the Gravelius order of Segments, BranchingPoints and Buds starting from 1" << endl;
  cout << "             so that parts belonging to an axis forking off from the present one get order one higher." << endl;
  cout << "             If -iniVal <value> starts from that value." << endl;
  cout << "             Then writes the tree to file treefile-neworder.xml" << endl;
  cout << "-setWoodyZero       Sets radius of all segments without foliage = 0, if -Hc <value> is set, only parts" << endl;
  cout << "             below Hc (in branches starting below Hc) are affected." << endl;
  cout << "-maxCrownWidth      Prints out max horizontal distance from start of the branch of segments" << endl;
  cout << "                    in branches (all segments order >= 2). If -foliage is set only segments" << endl;
  cout << "                    with foliage are considered." << endl;
  cout << "-collectSegmentLengths Prints sum (m) of segment lengths of all orders. If -order <value> is set," << endl;
  cout << "                       only for that order. " << endl;
  cout << "-noBranchWhorls     Prints number of branch whorls (= branches leaving from exactly same spot," << endl;
  cout << "                    this should be developed to allow for small vertical distance between)" << endl;
  cout << "-meanBranchAngle    Prints mean branch angle and chord angle of branches. " << endl;
  cout << "-setSapwoodArea     Sets heartwood radius so that sapwood area of a segment corresponds to SA above." << endl;
  cout << "                    Branch to stem 0.75, higher order forking 0.89 only so much of SA matched - cf." << endl;
  cout << "                    Sievanen et al. 2008. Writes result to file sapwoodarea.xml" << endl;
  cout << "                    If -correctR follows, sets R of segment to match sapwood requirement in case it doesn't" 
       << endl;
  cout << "-printLastSegInfo   Prints out information about last segment in each axis" << endl;
  cout << "-printStemDiameters    Prints out taper curve" << endl;
  cout << endl;

}


int main(int argc, char** argv) {

  if(argc < 2) {
    Usage();
    exit(0);
  }

  Tree<ETCfSegment,ETCfBud> conifer(Point(0,0,0.0),PositionVector(0,0,1));
  Tree<ETHwSegment,ETHwBud> deciduous(Point(0,0,0.0),PositionVector(0,0,1));

  bool is_conifer = true;
  if(CheckCommandLine(argc,argv,"-deciduous"))
    is_conifer = false;

    string tree_file = argv[1];

  if(!CheckCommandLine(argc,argv,"-makeTree")) {
    if(is_conifer) {
      XMLDomTreeReader<ETCfSegment,ETCfBud> tree_reader;
      tree_reader.readXMLToTree(conifer, tree_file);
    }
    else {
      XMLDomTreeReader<ETHwSegment,ETHwBud,Triangle> tree_reader;
      tree_reader.readXMLToTree(deciduous, tree_file);
    }
  }


  //==========================================================================
  //
  // Surface area of woody parts
  //
  //=========================================================================
  string cla;
  if(ParseCommandLine(argc,argv,"-mantleArea",cla)) {
    int value = atoi(cla.c_str());
  
    LGMdouble cb = 0.0;
    bool crown = false;
    if(CheckCommandLine(argc,argv,"-crown")) {
      crown = true; 
      DCLData dcl;
      if(is_conifer) {
	AccumulateDown(conifer,dcl,AddBranchWf(),DiameterCrownBase<ETCfSegment,ETCfBud>());
      }
      else {
	AccumulateDown(deciduous,dcl,AddBranchWf(),DiameterCrownBase<ETHwSegment,ETHwBud>());
      }
      cb = dcl.HCrownBase();
    }

bool foliage = false;
    if(CheckCommandLine(argc,argv,"-foliage")) {
      foliage = true;
    }

    LGMdouble area = 0.0;
    LGMdouble ssum = 0.0;
    bool s_sum = false;
    if(CheckCommandLine(argc,argv,"-sum")) {
      s_sum = true;
    }

    if(!s_sum) {
      area = 0.0;
      if(is_conifer) {
	if(crown) {
	  CollectMantleArea<ETCfSegment,ETCfBud> cma(value,foliage,cb);
	  area = Accumulate(conifer,area,cma);
	} else {
	  CollectMantleArea<ETCfSegment,ETCfBud> cma(value,foliage);
	  area = Accumulate(conifer,area,cma);
	}
      } else {
	if(crown) {
	  CollectMantleArea<ETHwSegment,ETHwBud> cma(value,foliage,cb);
	  area = Accumulate(deciduous,area,cma);
	} else {
	  CollectMantleArea<ETHwSegment,ETHwBud> cma(value,foliage);
	  area = Accumulate(deciduous,area,cma);
	}
      }
      cout << tree_file << " " << value << " " << area << endl;
      exit(0);
    } else {
      for(int i = value; i < 8; i++) {
	area = 0.0;
	if(is_conifer) {
	  if(crown) {
	    CollectMantleArea<ETCfSegment,ETCfBud> cma(i,foliage,cb);
	    area = Accumulate(conifer,area,cma);
	  } else {
	    CollectMantleArea<ETCfSegment,ETCfBud> cma(i,foliage);
	    area = Accumulate(conifer,area,cma);
	  }
	} else {
	  if(crown) {
	    CollectMantleArea<ETHwSegment,ETHwBud> cma(i,foliage,cb);
	    area = Accumulate(deciduous,area,cma);
	  } else {
	    CollectMantleArea<ETHwSegment,ETHwBud> cma(i,foliage);
	    area = Accumulate(deciduous,area,cma);
	  }
	}
	ssum += area;
      }
      cout << tree_file << " " << ssum << endl;
      exit(0);
    }
  }


  //===============================================================================
  //
  // Tree height and height of crown base
  // (only conifer at the moment)
  //
  //==============================================================================

  if(CheckCommandLine(argc,argv,"-heights")) {
    if(is_conifer) {
      DCLData dcl;
      AccumulateDown(conifer,dcl,AddBranchWf(),DiameterCrownBase<ETCfSegment,ETCfBud>());
      cout << tree_file << " " << GetValue(conifer,LGAH) << " " << dcl.HCrownBase() << endl;
    }
    exit(0);
  }

  //===============================================================================
  //
  // Number of different classes of TreeCompartments
  //
  //==============================================================================

  if(CheckCommandLine(argc,argv,"-noCompartments")) {
    if(is_conifer) {
      int n = 0;
      n = Accumulate(conifer,n,CountBuds<ETCfSegment,ETCfBud>());
      cout << "Buds " << n << endl;
      n = 0;
      n = Accumulate(conifer,n,CountTreeSegments<ETCfSegment,ETCfBud>());
      cout << "Segments " << n << endl;
      n = 0;
      n = Accumulate(conifer,n,CountCfTreeSegmentsWithFoliage<ETCfSegment,ETCfBud>());
      cout << "Segments with foliage " << n << endl;
      n = 0;
      n = Accumulate(conifer,n,CountBranchingPoints<ETCfSegment,ETCfBud>());
      cout << "Branching points " << n << endl;
      n = 0;
      n = Accumulate(conifer,n,CountAxes<ETCfSegment,ETCfBud>());
      cout << "Axes " << n << endl;
    } else {
      int n = 0;
      n = Accumulate(deciduous,n,CountBuds<ETHwSegment,ETHwBud>());
      cout << "Buds " << n << endl;
      n = 0;
      n = Accumulate(deciduous,n,CountTreeSegments<ETHwSegment,ETHwBud>());
      cout << "Segments " << n << endl;
      n = 0;
      n = Accumulate(deciduous,n,CountHwTreeSegmentsWithFoliage<ETHwSegment,ETHwBud,Triangle>());
      cout << "Segments with foliage " << n << endl;
      n = 0;
      n = Accumulate(deciduous,n,CountBranchingPoints<ETHwSegment,ETHwBud>());
      cout << "Branching points " << n << endl;
      n = 0;
      n = Accumulate(deciduous,n,CountAxes<ETHwSegment,ETHwBud>());
      cout << "Axes " << n << endl;
    }
    exit(0);
  }

  //===============================================================================
  //
  // Vizualize tree with LIGNUM visualization -vizTree
  //
  //==============================================================================

  if(CheckCommandLine(argc,argv,"-vizTree")) {
    LGMVisualization viz;
    viz.InitVisualization(argc,argv);
    // textures 512x512
    if(is_conifer) {
      viz.AddCfTree<ETCfSegment,ETCfBud>(conifer, "conifer_stem.bmp", "needle.tga");
      viz.ResetCameraPosition(GetValue(conifer,LGAH));
    } else {
      viz.AddHwTree<ETHwSegment,ETHwBud,Triangle>(deciduous, "deciduous_stem.bmp", "leaf.tga");
      viz.ResetCameraPosition(GetValue(deciduous,LGAH));
    }
    viz.SetMode(SOLID);
    viz.StartVisualization();
    exit(0);
  }

  //===============================================================================
  //
  // Foliage area and mass by Gravelius orders and total
  //
  //==============================================================================

  if(CheckCommandLine(argc,argv,"-foliageAreaMass")) {
    bool by_order = false;
    if(CheckCommandLine(argc,argv,"-byOrder")) {
      by_order = true;
    }
    bool title = true;
    if(CheckCommandLine(argc,argv,"-noTitle")) {
      title = false;
    }
    if(title) {
      if(by_order) {
	cout << " Tree     order(stem = 1)    foliage area (m2)  mass (kg C)" << endl;
      }
      else {
	cout << " Tree  foliage area (m2)  mass (kg C)" << endl;
      }
    }

    LGMdouble area = 0.0, mass = 0.0;
    LGMdouble tarea = 0.0, tmass = 0.0;
       
    for(int i = 1; i < 8; i++) {
      area = 0.0;
      mass = 0.0;
      if(is_conifer) {
	CollectFoliageArea<ETCfSegment,ETCfBud> cfa((double)i);
	area = Accumulate(conifer,area,cfa);
	tarea += area;
	CollectFoliageMass<ETCfSegment,ETCfBud> cfm((double)i);
	mass = Accumulate(conifer,mass,cfm);
	tmass += mass;
      }
      else {
	CollectFoliageArea<ETHwSegment,ETHwBud> cfa((double)i);
	area = Accumulate(deciduous,area,cfa);
	tarea += area;
	CollectFoliageMass<ETHwSegment,ETHwBud> cfm((double)i);
	mass = Accumulate(deciduous,mass,cfa);
	tmass += mass;
      }
      if(by_order) {
	cout << tree_file << " " << i << " " << area << " " << mass << endl;
      }
    }
      cout << tree_file << " " << tarea << " " << tmass << endl;

    exit(0);
  }

  //===============================================================================
  //
  // Sum of segment lengths by Gravelius orders and total

  //==============================================================================

  if(CheckCommandLine(argc,argv,"-segmentLengths")) {

    double go = 0.0;
    ForwardGO gg;

    PropagateUp(conifer,go,gg);

    bool by_order = false;
    if(CheckCommandLine(argc,argv,"-byOrder")) {
      by_order = true;
    }
    bool title = true;
    if(CheckCommandLine(argc,argv,"-noTitle")) {
      title = false;
    }
    if(title) {
      if(by_order) {
	cout << " Tree     order(stem = 1)    total length (m)" << endl;
      }
      else {
	cout << " Tree  total length (m)" << endl;
      }
    }

    LGMdouble len = 0.0;
    LGMdouble tlen = 0.0;
       
    for(int i = 1; i < 8; i++) {
      len = 0.0;
      if(is_conifer) {
	CollectSegmentLengths<ETCfSegment,ETCfBud> csl((double)i);
	len = Accumulate(conifer,len,csl);
	tlen += len;
      }
      else {
	CollectSegmentLengths<ETHwSegment,ETHwBud> csl((double)i);
	len = Accumulate(deciduous,len,csl);
	tlen += len;
      }
      if(by_order) {
	cout << tree_file << " " << i << " " << len << endl;
      }
    }
      cout << tree_file << " " << tlen << endl;

    exit(0);
  }


  //=================================================================================
  // Fractal Dimension  -- only conifers at the moment!
  //=================================================================================
  cla.clear();
  if(ParseCommandLine(argc,argv,"-fractalDimension",cla)) {
    if(!is_conifer) {
      cout << "Fractal dimension possible only for conifers!" << endl;
      exit(0);
    }

    int no_points = atoi(cla.c_str());

    determineFractalDimension(conifer, argc, argv, no_points, tree_file);

    exit(0);
  }


  //=================================================================================
  // Lacunarity  -- only conifers at the moment!
  //=================================================================================
  cla.clear();
  if(ParseCommandLine(argc,argv,"-lacunarity",cla)) {
    if(!is_conifer) {
      cout << "Lacunarity possible only for conifers!" << endl;
      exit(0);
    }

    LGMdouble edge = atof(cla.c_str());

    lacunarity(conifer, argc, argv, edge, tree_file);

    exit(0);
  }


  //==============================================================================================
  // Dimensions of crown ellipsoid

  // The crown volume (between highest and lowest TreeSegment with foliage) is evaluated in
  // 20 cm high circular disks. The radius of each disk is maximum distance of foliage TreeSegment
  // (farthest of the ends). The radius of a circular ellpsoid that is between crown base
  // and tree top is calculated and the result is written to console.
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-ellipsoid")) {

    Point ex_point(0.0,0.0,R_HUGE);
    ex_point = Accumulate(conifer,ex_point,LowestCfSegmentWithFoliage<ETCfSegment,ETCfBud>());
    LGMdouble H_low = ex_point.getZ();
    ex_point = Point(0.0,0.0,-R_HUGE);
    ex_point = Accumulate(conifer,ex_point,HighestCfSegmentWithFoliage<ETCfSegment,ETCfBud>());
    LGMdouble H_high = ex_point.getZ();
    
    if(H_high <= H_low) {
      cout << "H_high <= H_low " << H_high << " " << H_low << endl;
      exit(0);
    }

    LGMdouble dist = H_high - H_low;
    LGMdouble step = 0.2;             //Thickness of disk = 0.2 m

    int layers = (int)(dist/step) + 1;

    Axis<ETCfSegment,ETCfBud>& ax = GetAxis(conifer);
    Point treeBase = GetPoint(*GetFirstTreeCompartment(ax));
    
    LGMdouble crown_volume = 0.0;
    for(int i = 0; i < layers; i++) {
      LGMdouble minH = (double)i * step + H_low;
      LGMdouble maxH = minH + step;
      LGMdouble angle = PI_VALUE / 2.0;
      LGMdouble maxR = 0.0;
      for(int j = 0; j < 4; j++) {
        double dir = (double)j * angle;
        FindRFunctor<ETCfSegment,ETCfBud>  findR(minH, maxH, dir, angle,
						 treeBase);
        double R = 0.0;
        R = Accumulate(conifer, R, findR);
	if(R > maxR)
	  maxR = R;
      }
      crown_volume += PI_VALUE * maxR * maxR * step;
    }

 
    // Base of ellipsoid = crown base
    // Top of ellipsoid = tree top

    LGMdouble H = GetValue(conifer, LGAH);
    DCLData     dcl;
    AccumulateDown(conifer,dcl,AddBranchWf(),DiameterCrownBase<ETCfSegment,ETCfBud>());
    LGMdouble Hc = dcl.HCrownBase();

    //Volume of ellipsoid = (4/3)*pi*r^2*(Htop-Hbase)/2
    LGMdouble r_ellipsoid = sqrt(2.0*crown_volume/((4.0/3.0)*PI_VALUE*(H-Hc)));

    // Also needle area

    double area = 0.0;
    CollectFoliageArea<ETCfSegment,ETCfBud> cfa;
    area = Accumulate(conifer,area,cfa);

    cout << tree_file << " " << (H + Hc)/2.0 << " " << r_ellipsoid << " " << (H - Hc)/2.0 << " " << area << endl;

    exit(0);
  }

  //==============================================================================================
  // Set positions and orientation of shoots (Segments with needles) randomly in a given 
  // ellipsoid.
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-randomizeShoots")) {

    //Read cl arguments that must be there
    string cla;
    LGMdouble zc = 0.0;
    if(ParseCommandLine(argc,argv,"-Zc",cla))
      zc = atof(cla.c_str());
    else {
      cout << "CLA -Zc is missing!" << endl;
      exit(0);
    }

    Point pc = GetPoint(conifer); // Position of tree = x,y position of center
                                  // of ellipsoid
    pc.setZ(zc);    

    LGMdouble rv = 0.0, rh = 0;
    if(ParseCommandLine(argc,argv,"-rh",cla))
      rh = atof(cla.c_str());
    else {
      cout << "CLA -rh is missing!" << endl;
      exit(0);
    }
    if(ParseCommandLine(argc,argv,"-rv",cla))
      rv = atof(cla.c_str());
    else {
      cout << "CLA -rv is missing!" << endl;
      exit(0);
    }

    ran3_seed = -3924678;
    if (ParseCommandLine(argc,argv,"-seed", cla)){
     ran3_seed = atoi(cla.c_str());
     ran3_seed = -abs(ran3_seed);
    }

    bool woody_parts = false;
    if( CheckCommandLine(argc,argv,"-woodyParts") ) {
      woody_parts = true;
    }

    bool both_ends = false;
    if( CheckCommandLine(argc,argv,"-bothEnds") ) {
      woody_parts = true;
    }

    DCLData     dcl;
    AccumulateDown(conifer,dcl,AddBranchWf(),DiameterCrownBase<ETCfSegment,ETCfBud>());
    LGMdouble Hc = dcl.HCrownBase();

    SetShootPositions ssp(pc, rh, rv, Hc, both_ends, woody_parts);
    ForEach(conifer,ssp);

    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "-ellipsoid.xml";

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);
    cout << "Tree " << o_file << " written." << endl;


    exit(0);
  }


  //==============================================================================================
  //TreeSegmentinformation
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-treeSegmentInformation")) {
    bool title = true;
    PrintTreeSegmentInformationToFile<ETCfSegment,ETCfBud> ptsif("treesegmentinformation.dat", title);
    ForEach(conifer, ptsif);
  }


  //==============================================================================================
  // CheckNeedlesLast - write out last segments of axes that do not have needles
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-needlesLast")) {

    bool initial = true;
    AccumulateDown(conifer, initial, CheckNeedlesLast());

    exit(0);
  }


  //==============================================================================================
  // -multWf10   Multply needle mass by 10 and save tree in <filename>-10.xml
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-multWf10")) {

    ForEach(conifer,  MultWf10());


    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;

    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "-10.xml";
    writer.writeTreeToXML(conifer, o_file);

    exit(0);
  }


  //==============================================================================================
  // -foliageMassEqs   Multply needle mass by 10 and save tree in <filename>-10.xmlFoliage mass 
  // of the tree and fol. mass by equations.
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-foliageMassEqs")) {

    CollectFoliageMass<ETCfSegment,ETCfBud> cfm(-1);
    LGMdouble mass = Accumulate(conifer,mass,cfm);                 //kg C


    LGMdouble dbh = GetValue(conifer,LGADbh);            //m
    LGMdouble h = GetValue(conifer,LGAH);                //m

    LGMdouble dk = 2+1.25*100.0*dbh;
    LGMdouble Wf_Repola =  exp(-5.007 + 15.289*dk/(dk+6.0) - 5.896*h/(h+1.0) + 0.5*(0.097 + 0.123));
    LGMdouble Wf_Marklund = exp(-3.4781+12.1095*100.0*dbh/(100.0*dbh+7.0)+0.0413*h-1.565*log(h));
    LGMdouble Wf_Wirth = 271.154*pow(dbh,2.0916);
    LGMdouble Wf_Ilvnm =0.1628+0.0233*PI_VALUE*(100.0*dbh/2.0)*(100.0*dbh/2.0);

    if(!CheckCommandLine(argc,argv,"-noTitle")) {
      cout << tree_file << " Dbh Wf Wf_Repola Wf_Marklund Wf_Wirth Wf_Ilvnm" << endl;
    }
    cout << tree_file << " " << 100.0*dbh << " " << 2.0*mass << " " << Wf_Repola << " " <<  Wf_Marklund << " "
	 <<  Wf_Wirth << " " << Wf_Ilvnm << endl;
  }

  //==============================================================================================
  // -raumonen
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-raumonen")) {

    int line_number = 0;
    int parent_no = 0;
    int my_branchno = 1;

    if(is_conifer) {
      Axis<ETCfSegment,ETCfBud>& my_axis = GetAxis(conifer);
      writeCylinderData(line_number, parent_no, my_branchno, my_axis);
     }
    else {
      Axis<ETHwSegment,ETHwBud>& my_axis = GetAxis(deciduous);
      writeCylinderData(line_number, parent_no, my_branchno, my_axis);
    }
  }



  //==============================================================================================
  // -raumonenB
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-raumonenB")) {

    int parent_no = 0;
    int my_branchno = 1;
    PositionVector vi(0.0,0.0,1.0);

    if(is_conifer) {

      Axis<ETCfSegment,ETCfBud>& my_axis = GetAxis(conifer);
      writeBranchData(parent_no, my_branchno, vi, my_axis);
    }
    else {
      Axis<ETHwSegment,ETHwBud>& my_axis = GetAxis(deciduous);
      writeBranchData(parent_no, my_branchno, vi, my_axis);
    }

  }

  //==============================================================================================
  //  -sideBranches
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-sideBranches")) {

    int e = tree_file.find(".xml");
    int a = tree_file.rfind("/") + 1;
    string o_file = tree_file.substr(a,e-a);

    if(is_conifer) {
      Axis<ETCfSegment,ETCfBud>& my_axis = GetAxis(conifer);
      list<TreeCompartment<ETCfSegment,ETCfBud>*> cl = GetTreeCompartmentList(my_axis);
      list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator aI;
      double base_dist = 0.0;
      for(aI = cl.begin(); aI != cl.end(); aI++) {
	if(ETCfSegment* seg = dynamic_cast<ETCfSegment*>(*aI) ) {
	  base_dist += GetValue(*seg, LGAL);
	}
      
	if(BranchingPoint<ETCfSegment,ETCfBud>* bp =
	   dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*aI)){
	  list<Axis<ETCfSegment,ETCfBud>*>& axis_ls = GetAxisList(*bp);
	  list<Axis<ETCfSegment,ETCfBud>*>::iterator bpI;
	    int apexes = 0, segs = 0;
	    double link_length = 0.0;

	  for(bpI = axis_ls.begin(); bpI != axis_ls.end(); bpI++) {
	    Tree<ETCfSegment,ETCfBud> this_tree(**bpI);
	    int a_ini = 0, seg_ini = 0;
	    double link_lini = 0.0;
	    apexes += Accumulate(this_tree, a_ini, CountAxesProper<ETCfSegment,ETCfBud>());
	    segs += Accumulate(this_tree, seg_ini, CountTreeSegments<ETCfSegment,ETCfBud>());
	    link_length += Accumulate(this_tree, link_lini, TotalSegmentLength<ETCfSegment,ETCfBud>());
	  }
	  if(apexes > 0 || segs > 0 || link_length > 0.0) {
	    cout << o_file << " " << base_dist << " " << apexes << " " << segs << " " << link_length << endl;
	  }
	}
      }

      //     else {
      //       Axis<ETHwSegment,ETHwBud>& my_axis = GetAxis(deciduous);
      //       writeBranchData(parent_no, my_branchno, vi, my_axis);
      //     }
  

    }
    exit(0);
  }

  //==============================================================================================
  //  STAR and other values of segments to file
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-STARvalues")) {

    if(is_conifer) {
      string starfile = "starfile.dat";
      int runs = 1000;
      if(ParseCommandLine(argc,argv,"-runs",cla)) {
	runs = atoi(cla.c_str());
      }
      if(ParseCommandLine(argc,argv,"-starFile",cla)) {
	starfile = cla;
      }    
      ForEach(conifer, CalculateSTAR(starfile, runs)); 
    }
    exit(0);
  }


  //==============================================================================================
  //  -makeTree
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-makeTree")) {
    LGMdouble Rf = 0.04;
    LGMdouble L = 0.15;
    LGMdouble fol_den = 130.0;  //Mean value from P. Stenberg and S. Palmroth and B. Bond and 
    //D. Sprugel and H. Smolander Tree Physiology  21  805-814  (2001)
    int n_s = 1000;
    if(ParseCommandLine(argc,argv,"-nShoot",cla)) {
      n_s = atoi(cla.c_str());
    }
    if(ParseCommandLine(argc,argv,"-Rf",cla)) {
      Rf = atof(cla.c_str());
    }
    if(ParseCommandLine(argc,argv,"-len",cla)) {
      L = atof(cla.c_str());
    }
    if(ParseCommandLine(argc,argv,"-fol_den",cla)) {
      fol_den = atof(cla.c_str());
    }

    LGMdouble Af = PI_VALUE * Rf * Rf * L * fol_den;
    LGMdouble sf = 28;     // m2 / kg C
    LGMdouble Wf = Af/sf;
    LGMdouble R = 0.0001;
    LGMdouble order = 1.0;
    LGMdouble Rh = 0.0;

    PositionVector up(0.0,0.0,1.0);
    Point zero(0.0,0.0,0.0);

    Axis<ETCfSegment,ETCfBud>& axis = GetAxis(conifer);
    ETCfBud* bud = new ETCfBud(zero, up, order, &conifer); 
    InsertTreeCompartment(axis, bud);

    for(int i = 0; i < n_s; i++) {
      ETCfSegment* ts = new ETCfSegment(zero, up, order, L, R, Rh, &conifer);
      //        SetValue(*ts,LGAR,R);
      //        SetValue(*ts,LGARh,Rh);
      //        SetValue(*ts,LGAL,L);
      SetValue(*ts,LGAsf,sf);
      SetValue(*ts,LGAWf,Wf);
      SetValue(*ts,LGARf,Rf);

      BranchingPoint<ETCfSegment,ETCfBud>* bp = new  BranchingPoint<ETCfSegment,
	ETCfBud>(zero,up,order,&conifer);
      
      InsertTreeCompartmentSecondLast(axis, ts);
      InsertTreeCompartmentSecondLast(axis, bp);
    }

    //After this the tree should be an upright sequence of shoot cylinders
    PropagateUp(conifer, zero, ConnectTree<ETCfSegment,ETCfBud>());

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, "simple.xml");


    exit(0);
  }


  //==============================================================================================
  //  -moveTree 
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-moveTree")) {
    Point p(0.0, 0.0, 0.0);
    if(ParseCommandLine(argc,argv,"-X",cla)) {
      p.setX(-atoi(cla.c_str()));
    }
    if(ParseCommandLine(argc,argv,"-Y",cla)) {
      p.setY(atoi(cla.c_str()));
    }
    if(ParseCommandLine(argc,argv,"-Z",cla)) {
      p.setZ(atoi(cla.c_str()));
    }

    if(is_conifer) {
      MoveTree<ETCfSegment,ETCfBud> move(p-GetPoint(conifer), conifer);
      ForEach(conifer, move);
    } else {
      MoveTree<ETHwSegment,ETHwBud> move(p-GetPoint(deciduous), deciduous);
      ForEach(deciduous, move);
    }
 
    // iterator to find the . position in the filename
    std::string::iterator beg = std::find(tree_file.begin(), tree_file.end(), '.'); 
    // iterator to find all the characters from . till the end of file
    std::string::iterator end = std::find(beg, tree_file.end(), ' ');
    string end_string = "_move.xml"; 
    tree_file.replace(beg, end, end_string);

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, tree_file);

    exit(0);
  }


  //==============================================================================================
  //  -setRf <value>
  //
  //===============================================================================================

  if(ParseCommandLine(argc,argv,"-setRf",cla)) {
    double set_rf = atof(cla.c_str());

    ForEach(conifer, SetRf(set_rf));


    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, "setrf.xml");


    exit(0);
  }

  //==============================================================================================
  //  -monster
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-monster")) {

    //When Wf has a large value, also Af is large
    ForEach(conifer, SetWf100());

    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "_monster.xml";
    cout << o_file << endl;
    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);
    exit(0);
  }

  //==============================================================================================
  //  -adjustWf
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-adjustWf")) {
  
    //at the moment Wf *= 1.0/(4.0-3.7*treeh/20.0)
    double h = GetValue(conifer, LGAH);

    ForEach(conifer, AdjustWf(h));

    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "_adjwf.xml";
    cout << o_file << endl;
    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);
    exit(0);
  }

  //==============================================================================================
  //  -adjustSf
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-adjustSf")) {
  
    //at the moment sf *= 1.15
    ForEach(conifer, AdjustSf());
  
    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "_adjsf.xml";
    cout << o_file << endl;
    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);
    exit(0);
  }


  //==============================================================================================
  // -streit
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-streit")) {

    if(is_conifer) {
      int initial_number = 0;
      Accumulate(conifer,initial_number, SetSegmentNumber());
      int base_number = -1;
      PropagateUp(conifer, base_number, WriteCylinderDataStreit());
     }
    else {
      cout << "Not yet implemented for deciduous trees." << endl;
    }

    exit(0);
  }


  //==============================================================================================
  // -countTreeSegments
  //
  //===============================================================================================



  if(CheckCommandLine(argc,argv,"-countTreeSegments")) {

    if(is_conifer) {
      unsigned int n_segments = 0;
      n_segments = Accumulate(conifer, n_segments, CountTreeSegments<ETCfSegment,ETCfBud>());
      cout << "There are " << n_segments << " segments." << endl;

    }
    else {
      //      n_segments = Accumulate(conifer, n_segments, CountTreeSegments<ETHwSegment,ETHwBud>());
      cout << "Not implemented for deciduous trees!" << endl;
    }

    exit(0);
  }


  //==============================================================================================
  //
  //  -simplify
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-simplify")) {
    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "_simple.xml";
    cout << o_file << endl;

    bool removed;
    int c = 0;
    do {
      removed = false;
      c++;
      removed = Accumulate(conifer, removed, SimplifyTree());
    } while(removed);
    cout << c - 1 << " removal rounds" << endl;

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);
  }


  //==============================================================================================
  //
  //  -separateTopShoots
  //
  //===============================================================================================

  cla.clear();
  if(ParseCommandLine(argc,argv,"-separateTopShoots",cla)) {
    double sep = atof(cla.c_str());
    sep *= 0.01;     // cm --> m
    
    double second_sep = 0.0;
    cla.clear();
    if(ParseCommandLine(argc,argv,"-secondSep",cla)) {
      second_sep = atof(cla.c_str());
      second_sep *= 0.01;    //cm --> m
    }

    SeparateTopShoots sts(sep, conifer, second_sep);
    //    bool oh = false;
    ForEach(conifer, sts);
    //Accumulate(conifer, oh, sts);

    Axis<ETCfSegment,ETCfBud>& axis = GetAxis(conifer);
    Point treeBase = GetPoint(*GetFirstTreeCompartment(axis));

    PropagateUp(conifer, treeBase, ConnectTree<ETCfSegment,ETCfBud>());

    stringstream ss;
    ss << static_cast<int>(100.0*sep);  //sep as part of file name
    string sep_str = ss.str();

    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "_separate" + sep_str + "cm.xml";
    cout << o_file << endl;

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);

    exit(0);
  }

  //==============================================================================================
  //
  //  -displayStructure
  //
  //===============================================================================================


  if(CheckCommandLine(argc,argv,"-displayStructure")) {

    if(is_conifer) {
      if(CheckCommandLine(argc,argv,"-byAxis")) {
	ForEach(conifer, DisplayStructureAxis<ETCfSegment,ETCfBud>());
      } else {
	DisplayStructure<ETCfSegment,ETCfBud>(conifer);
      }
    } else {
      if(CheckCommandLine(argc,argv,"-byAxis")) {
	ForEach(deciduous, DisplayStructureAxis<ETHwSegment,ETHwBud>());
      } else {
	DisplayStructure<ETHwSegment,ETHwBud>(deciduous);
      }
    }
    exit(0);
  }

  //==============================================================================================
  //
  //  -generateAxiom
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-generateAxiom")) {

    if(is_conifer) {
      if(CheckCommandLine(argc,argv,"-byAxis")) {
	ForEach(conifer, DisplayStructureAxis<ETCfSegment,ETCfBud>());
      } else {
	DisplayStructure<ETCfSegment,ETCfBud>(conifer);
      }
    } else {
      cout << "-generateAxiom not available for deciduous tree at the moment." << endl;
    }
    exit(0);
  }


  //==============================================================================================
  //
  //  -verticalLA
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-verticalLA")) {
    double step = 0.2;
    string cla;
    if(ParseCommandLine(argc,argv,"-step",cla)) {
      step = atof(cla.c_str());
    }
    double min_z = 0.0, max_z = 0.0;
  
    cla.clear();
    if(ParseCommandLine(argc,argv,"-min",cla)) {
      min_z = atof(cla.c_str());
    } else {
      if(is_conifer) {
	DCLData     dcl;
	AccumulateDown(conifer,dcl,AddBranchWf(),DiameterCrownBase<ETCfSegment,ETCfBud>());
	min_z = dcl.HCrownBase();
      } else {
	Axis<ETHwSegment,ETHwBud> ax =  GetAxis(deciduous);
	min_z  = (GetPoint(*GetFirstTreeCompartment(ax))).getZ();
      }
    }

    cla.clear();
    if(ParseCommandLine(argc,argv,"-max",cla)) {
      max_z = atof(cla.c_str());
    } else {
      if(is_conifer) {
	Axis<ETCfSegment,ETCfBud>& ax =  GetAxis(conifer);
	max_z  = (GetPoint(*GetTerminatingBud(ax))).getZ();
      } else {
	Axis<ETHwSegment,ETHwBud>& ax =  GetAxis(deciduous);
	max_z  = (GetPoint(*GetTerminatingBud(ax))).getZ();
      }
    }
    
    int dist_size = static_cast<int>((max_z - min_z)/step) + 1;
    vector<double> distn(dist_size, 0.0);
  
    if(is_conifer) {
      VerticalLeafAreaDistribution<ETCfSegment,ETCfBud> vlad(min_z, max_z, step);
      distn = Accumulate(conifer, distn, vlad);
    }  else {
      VerticalLeafAreaDistribution<ETHwSegment,ETHwBud> vlad(min_z, max_z, step);
      distn = Accumulate(deciduous, distn, vlad);
    }
    for(int i = 0; i < dist_size; i++) {
      cout << (static_cast<double>(i) + 0.5)*step << " " << distn[i] << endl;
    }
    exit(0);
  }

  //==============================================================================================
  //  -changeOrder <value>
  //
  //===============================================================================================

  cla.clear();
  if(ParseCommandLine(argc,argv,"-changeOrder",cla)) {
    if(argc < 3) {
      cout << "No value to change order!" << endl;
      exit(0);
    }
    double d_order = atof(cla.c_str());

    ForEach(conifer, ChangeOrder(d_order));

    int e = tree_file.find(".xml");
    string o_file = tree_file.substr(0,e) + "-neworder.xml";

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, o_file);
    cout << "Tree " << o_file << " written." << endl;
    exit(0);
  }


  //==============================================================================================
  //  -setOrder <value>
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-setOrder")) {
    int ini_val = 0.0; 
    string cla;
    if(ParseCommandLine(argc,argv,"-iniVal",cla)) {
      ini_val = atoi(cla.c_str()) - 1.0;
    }
    //Initial value to PropagateUp must be value of first axis - 1 !

    if(is_conifer) {
      PropagateUp(conifer, ini_val, SetGraveliusOrder<ETCfSegment,ETCfBud>());
      XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
      writer.writeTreeToXML(conifer, "gravelius-ordered.xml" );
    } else {
      PropagateUp(deciduous, ini_val, SetGraveliusOrder<ETHwSegment,ETHwBud>());
      XMLDomTreeWriter <ETHwSegment,ETHwBud>writer;
      writer.writeTreeToXML(deciduous, "gravelius-ordered.xml" );
    }

    cout << "Tree gravelius-ordered.xml written." << endl;
    exit(0);
  }


  //==============================================================================================
  //
  //  -maxCrownWidth 
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-maxCrownWidth")) {

    bool foliage = false;
    if(CheckCommandLine(argc,argv,"-foliage")) {
      foliage = true;
    }
    //Using Tree parameter LGPyc to store the distance is a hack
    if(is_conifer) {
      MaxCrownWidthInfo mcw_info(0.0, Point(0.0, 0.0, 0.0) ); 
      MaxCrownWidth mcw(foliage);
      SetValue(conifer, LGPyc, 0.0);
      PropagateUp(conifer, mcw_info, mcw);
      cout << GetValue(conifer, LGPyc) << endl;
     }
    else {
      cout << "Not yet implemented for deciduous trees." << endl;
    }

    exit(0);
  }

  //==============================================================================================
  //
  //  -collectSegmentLengths [-order <value>]
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-collectSegmentLengths")) {

    double order = -1.0;
    cla.clear();
    if(ParseCommandLine(argc,argv,"-order",cla)) {
      order = atof(cla.c_str());
    }

    double sum = 0.0;
    if(is_conifer) {
      sum = Accumulate(conifer, sum, CollectSegmentLengths<ETCfSegment,ETCfBud>(order));
     }
    else {
      sum = Accumulate(deciduous, sum, CollectSegmentLengths<ETHwSegment,ETHwBud>(order));
    }
    cout << sum << endl;

    exit(0);
  }



  //==============================================================================================
  //
  //  -noBranchWhorls
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-noBranchWhorls")) {
    int no_whorls = 0;
    if(is_conifer) {
      Axis<ETCfSegment,ETCfBud>& my_axis = GetAxis(conifer);
      list<TreeCompartment<ETCfSegment,ETCfBud>*> cl = GetTreeCompartmentList(my_axis);
      list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator aI;
      for(aI = cl.begin(); aI != cl.end(); aI++) {
      
	if( BranchingPoint<ETCfSegment,ETCfBud>* bp =
	   dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*aI) ) {
	  list<Axis<ETCfSegment,ETCfBud>*>& axis_ls = GetAxisList(*bp);
	  if( static_cast<int>(axis_ls.size()) > 1) {
	    int ind = 0;
	    list<Axis<ETCfSegment,ETCfBud>*>::iterator bpI;
	    for(bpI = axis_ls.begin(); bpI != axis_ls.end(); bpI++) {
	      list<TreeCompartment<ETCfSegment,ETCfBud>*> cmpl = GetTreeCompartmentList(**bpI);
	      if( static_cast<int>(cmpl.size()) > 2 ) {   //at least one segment in axis
		ind++;
	      }
	    }
	    if( ind > 1 ) {                //at least two branches in BranchingPoint
	      no_whorls++;
	    }
	  }
	}    //if( BranchingPoint<ETCfSegment ...
      }  //for(aI = cl.begin() ...
    } else {                      //Deciduous tree ... they seldom have branch whorls
      Axis<ETHwSegment,ETHwBud>& my_axis = GetAxis(deciduous);
      list<TreeCompartment<ETHwSegment,ETHwBud>*> cl = GetTreeCompartmentList(my_axis);
      list<TreeCompartment<ETHwSegment,ETHwBud>*>::iterator aI;
      for(aI = cl.begin(); aI != cl.end(); aI++) {
      
	if( BranchingPoint<ETHwSegment,ETHwBud>* bp =
	   dynamic_cast<BranchingPoint<ETHwSegment,ETHwBud>*>(*aI) ) {
	  list<Axis<ETHwSegment,ETHwBud>*>& axis_ls = GetAxisList(*bp);
	  if( static_cast<int>(axis_ls.size()) > 1) {
	    int ind = 0;
	    list<Axis<ETHwSegment,ETHwBud>*>::iterator bpI;
	    for(bpI = axis_ls.begin(); bpI != axis_ls.end(); bpI++) {
	      list<TreeCompartment<ETHwSegment,ETHwBud>*> cmpl = GetTreeCompartmentList(**bpI);
	      if( static_cast<int>(cmpl.size()) > 2 ) {   //at least one segment in axis
		ind++;
	      }
	    }
	    if( ind > 1 ) {                //at least two branches in BranchingPoint
	      no_whorls++;
	    }
	  }
	}    //if( BranchingPoint<ETHwSegment ...
      }  //for(aI = cl.begin() ...
    }
 
    cout << no_whorls << endl;
    exit(0);
  }

  //==============================================================================================
  //
  //  -meanBranchAngle
  //
  //===============================================================================================

  if(CheckCommandLine(argc,argv,"-meanBranchAngle")) {
    double mean_branch_angle = 0.0, mean_chord_angle = 0.0;
    PositionVector current_direction;
    int n_branch = 0;
    double chord_sum = 0.0;

    if(is_conifer) {
      Axis<ETCfSegment,ETCfBud>& my_axis = GetAxis(conifer);
      list<TreeCompartment<ETCfSegment,ETCfBud>*> cl = GetTreeCompartmentList(my_axis);
      list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator aI;

      for(aI = cl.begin(); aI != cl.end(); aI++) {
      
	if( ETCfSegment* seg = dynamic_cast<ETCfSegment*>(*aI) ) {
	  current_direction = GetDirection(*seg);
	}

	if( BranchingPoint<ETCfSegment,ETCfBud>* bp =
	   dynamic_cast<BranchingPoint<ETCfSegment,ETCfBud>*>(*aI) ) {
	  list<Axis<ETCfSegment,ETCfBud>*>& axis_ls = GetAxisList(*bp);
	    list<Axis<ETCfSegment,ETCfBud>*>::iterator bpI;
	    for(bpI = axis_ls.begin(); bpI != axis_ls.end(); bpI++) {
	      list<TreeCompartment<ETCfSegment,ETCfBud>*> cmpl = GetTreeCompartmentList(**bpI);
	      if( static_cast<int>(cmpl.size()) > 2 ) {   //at least one segment in axis
		TreeSegment<ETCfSegment,ETCfBud>* first_seg = GetFirstTreeSegment(**bpI);
		if(first_seg == NULL) {
		  cout << "First segment not present in Branch!" << endl;
		  exit(0);
		}

		//comparison with stem segment leaving from the same BranchingPoint -
		//if it exists, if not with the previous
		list<TreeCompartment<ETCfSegment,ETCfBud>*>::iterator cur = aI;
		cur++;
		if(cur != cl.end()) {
		  if( ETCfSegment* cur_seg = dynamic_cast<ETCfSegment*>(*cur) ) {
		    current_direction = GetDirection(*cur_seg);
		  }
		}

		PositionVector dir = GetDirection(*first_seg);
		double ac = Dot(current_direction, dir);
		mean_branch_angle += 180.0 * acos(ac) / PI_VALUE;
		n_branch++;
		Point branch_start = GetPoint(*first_seg);

		TreeSegment<ETCfSegment,ETCfBud>* last_seg = GetLastTreeSegment(**bpI);
		if(last_seg == NULL) {
		  cout << "Last segment not present in Branch!" << endl;
		  exit(0);
		}
		Point chord_end = GetEndPoint(*last_seg);
		PositionVector chord_dir(chord_end - branch_start);
		chord_dir.normalize();
		ac = Dot(current_direction, chord_dir);
		double chord_len = branch_start || chord_end;
		mean_chord_angle += chord_len * 180.0 * acos(ac) / PI_VALUE;
		chord_sum += chord_len;
	      }
	    }
	}    //if( BranchingPoint<ETCfSegment ...
      }  //for(aI = cl.begin() ...
    } else {                      //Deciduous tree ... they seldom have branch whorls
//       Axis<ETHwSegment,ETHwBud>& my_axis = GetAxis(deciduous);
//       list<TreeCompartment<ETHwSegment,ETHwBud>*> cl = GetTreeCompartmentList(my_axis);
//       list<TreeCompartment<ETHwSegment,ETHwBud>*>::iterator aI;
//       for(aI = cl.begin(); aI != cl.end(); aI++) {
      
// 	if( BranchingPoint<ETHwSegment,ETHwBud>* bp =
// 	   dynamic_cast<BranchingPoint<ETHwSegment,ETHwBud>*>(*aI) ) {
// 	  list<Axis<ETHwSegment,ETHwBud>*>& axis_ls = GetAxisList(*bp);
// 	  if( static_cast<int>(axis_ls.size()) > 1) {
// 	    int ind = 0;
// 	    list<Axis<ETHwSegment,ETHwBud>*>::iterator bpI;
// 	    for(bpI = axis_ls.begin(); bpI != axis_ls.end(); bpI++) {
// 	      list<TreeCompartment<ETHwSegment,ETHwBud>*> cmpl = GetTreeCompartmentList(**bpI);
// 	      if( static_cast<int>(cmpl.size()) > 2 ) {   //at least one segment in axis
// 		ind++;
// 	      }
// 	    }
// 	    if( ind > 1 ) {                //at least two branches in BranchingPoint
// 	      no_whorls++;
// 	    }
// 	  }
// 	}    //if( BranchingPoint<ETHwSegment ...
//       }  //for(aI = cl.begin() ...
    }
 
    if(n_branch > 0 && chord_sum > 0.0) {
      cout << mean_branch_angle / static_cast<double>(n_branch) << " "
	   << mean_chord_angle / chord_sum << endl;
	} else {
      cout << "0.0  0.0" << endl;
    }

    exit(0);
  }



  if(CheckCommandLine(argc,argv,"-printLastSegInfo")) {
    cout << "posx posy posz dirx diry dirz order R Rh Rf L Wf Af" << endl;  

    ForEach(conifer,PrintLastSegInfo());
    
    exit(0);
  }

  if(CheckCommandLine(argc,argv,"-setSapwoodArea")) {
    bool correct_r = false;
    if(CheckCommandLine(argc,argv,"-correctR")) {
      correct_r = true;
    }

   // Amount of sapwood going down from a lateral axis as afunction of its 
   // Gravelius order (stem == 1)
   double vv[] = {1.0, 0.75, 0.89, 0.89, 0.89, 0.89, 0.89};
   std::vector<double> fractions(begin(vv), end(vv));
   //   vector<double> fractions {1.0, 0.75, 0.89, 0.89, 0.89, 0.89, 0.89};

   PSAD_Info ini(0.0);     //Inital sapwood going down == 0
   AccumulateDown(conifer, ini, PartSapwoodAreaDown(fractions),SetRh(correct_r));

    XMLDomTreeWriter<ETCfSegment,ETCfBud> writer;
    writer.writeTreeToXML(conifer, "sapwoodarea.xml");
    cout << "Tree  to file sapwoodarea.xml" << endl;

    exit(0);
 }


 if(CheckCommandLine(argc,argv,"-printStemDiameters")) {

   int dummy = 0;
   cout << "Height  R  Rh" << endl;
   PropagateUp(conifer, dummy, PrintStemDiameters());

   exit(0);
 }

}
