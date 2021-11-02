# Useful analyses and modifications of a Lignum tree

This project must reside under lignum-core, that is, in the directory lignum-core/

Evaluatetree collects a number of short programs that either analyze of modify a Lignum tree. This collection has developed as a way to access/modify trees for various purposes. This is in a way a front end for the functors in lignum-core/stl-lignum/TreeFunctors.h (only a small part of them are used). <br />
The tree is read from a xml file, default is conifer. If input tree is a deciduous one, Triangle is the leaf shape. Analyzing deciduous trees with other leaf shapes is not possible at the moment. Only a limited number analyses/actions are possible for deciduous trees. <br />

Usage: ./evaluate  tree.xml [-deciduous] [-mantleArea value [-foliage] [-crown] [-sum] ] [-foliageAreaMass [-byOrder] [-noTitle]] <br />
       [-heights] [-noCompartments] [-fractalDimension no_points [-fractalDm minbox] [-fractalD0 value] ] <br />
       [-lacunarity box [-no_points number]] [-ellipsoid] [-randomizeShoots -Zc value -rh value -rv value] [-seed value]] <br />
       [-woodyParts] [-bothEnds]] [-treeSegmentInformation] [-needlesLast] [-multWf10] [-foliageMassEqs [-noTitle]]  [-raumonen] [-raumonenB] <br />
       [-sideBranches] [-STARvalues [-runs value] [-starFile name]] [-makeTree [-nShoot n] [-Rf value] [-len value] [-fol_den value]] <br />
       [-moveTree [-X value] [-Y value] [-Z value]] [-segmentLengths [-byOrder] [-noTitle]] [-setRf value] <br />
       [-monster] [-adjustWf] [-adjustSf] [-streit] [-countTreeSegments] [-simplify] <br />
       [-separateTopShoots value [-secondSep value] ] [-displayStructure [-byAxis]] [-verticalLA [-step value] [-min value] [-max value] ] <br />
       [-changeOrder value] [-setOrder [-iniVal value]] [-setWoodyZero [-Hc value] ] [-maxCrownWidth [-foliage]] [-collectSegmentLengths [-order value]] <br />
       [-noBranchWhorls] [-meanBranchAngle] [-setSapwoodArea [-correctR]] [-printLastSegInfo] [-printStemDiameters] <br />


Input tree is assumed to be a conifer, -deciduous changes it to deciduous. <br />
> -mantleArea    evaluates surface area of trunk and branches without foliage by gravelius orders and writes to console. With -foliage also parts with foliage are included. <br />
>               If -crown is set, only mantle area of stemsegments inside crown are evaluated, segments order == value are evaluated, if -sum is set order >= value segments evaluated. <br />
> -foliageAreaMass   Evaluates foliage area (total in conifers, one-sided in deciduous) and mass and writes to console. If flag -byOrder is set writes also by Gravelius orders and total. If flag -noTitle only values are written <br />
> -heights           Tree height and height of crown base (only conifer at the moment) <br />
> -noCompartments    Number of TreeCompartment in different classes to console. <br />
> -fractalDimension  writes to console no. of occupied cubes as a function of cube size. That is; box counting. Shoot cylinders are covered with points that are used to test whether a cube is occupied. no_points specifies how many to cover a cylinder. Only for pines at the moment. <br />
> -fractalDm <minbox>    minimum box side lentgth given as part of original box size. Default = 1000 meaning minimun side length is 1/1000 of original. <br />
> -fractalD0 value     Starts d iteration from a given value = span/value span is written out when this program starts. <br />
> -lacunarity box Only for pines at the moment. <br />
> -CACReflect  If part of gliding box window goes outside of grid if Centered Allain-Cloitre calculation, it is reflected to other side of the grid (periodical boundary). <br />
> -ellipsoid   The crown volume (between highest and lowest TreeSegment with foliage) is evaluated in 20 cm high circular disks. The radius of each disk is maximum distance of foliage TreeSegment (farthest of the ends). The radius of a circular ellpsoid that is between crown base and tree top is calculated and the result is written to console (tree_file height ofellipse center horizontal radius vertical radius). <br />
> -randomizeShoots   Moves segments with foliage above crown base into random positions in ellipsoid with center, horizontal, and vertical radii given by -Zc, -rh, and  -rv. If -woodyParts is specified, also segments without foliage are moved. If -bothEnds is set both ends of segment are checked for being in ellipsoid. Normally only midpoint is checked. -seed value sets seed value for the random number generator. Writes the ellipsoid tree into file input_file-ellisoid.xml. <br />
> -treeSegmentInformation      Writes (all) information about TreeSegments to file treesegmentinformation.dat using functor ../stl-lignum/include/PrintTreeSegmentInformationToFile<TS,BUD> <br />
> -needlesLast Writes out information about last segments of axes that do not have foliage. <br />
> -multWf10    Multply needle mass by 10 and save tree in <filename>-10.xml <br />
> -foliageMassEqs              Evaluates foliage biomass of the tree and claculates fol. biomass with the aid of dbh and height using several equations: Repola, Marklund, Wirth and Ilvesniemi. The equations are taken from ../DigitalTree/cf-main.cc. Writes to console, if -noTitle writes only values. <br />
> -sideBranches   Writes to console information along the main Axis (meant for fine root analysis) Dist. from base, # segments, # apexes, total segment length <br />
> -STARvalues  Writes STAR values and other information of segments to a file. Works only for conifers. If filename is not given in -starFile <name> writes to starvalues.dat. STAR is analyzed using MC, default repetitions = 1000, can be changed by -runs value <br />
> -makeTree    Creates a tree (axis & treesegments) consisting of given number (-nShoot) of identical shoots. Default values are used if parameters are not specified. -Rf = radius (m) of shoot cylinder, -len = length(m), -fol_den = foliage density (m2/m3) in shoot cylinder (all-sided needle area), woody radius = 0.0001 m <br />
> -moveTree    Moves a tree from its (base) current position to one given by -X value -Y value -Z value (if omitted position = (0, 0, 0,)). Writes to file inputfile_moved.xml <br />
> -segmentLengths    Works like -foliageAreaMass but for lengths of TreeSegments. <br />
> -setRf value     Sets foliage radius to value in segments that have foliage <br />
> -monster     Sets Wf of shoots (that have foliage) = 100 kg C, their transparency is thus = 0: monster trees. Writes tree to file *_monster.xml <br />
> -adjustWf    Adjusts Wf (at the moment Wf *= 1.0/(4.0-3.7*treeh/20.0) ), writes to *_adjwf.xml file. <br />
> -adjustSf    Adjusts sf (at the moment sf *= 1.15, writes to *_adjsf.xml file. <br />
> -simplify    Creates output tree in which adjacent segments have been lumped together if they had not lateral branches between. Output tree contains only architectural information (lengths, radii, orientation) <br /> 
> -separateTopShoots value  Creates a shoot of length separation (cm) before the top shoot. Writes to <treefile>_separateXXcm.xml where XX = separation <br />
> -displayStructure   Writes the structure of the tree on console (by DisplayStructure of TreeFunctorI.h). If -byAxis writes contents of each Axis on own line (by DisplayStructureAxis of TreeFunctorI.h) <br />
> -verticalLA  Evaluates vertical needle area (all-sided)/leaf area distribution in 20 cm slices, -step value changes this default. <br />
> -changeOrder Changes order of compartments by value. No questions asked. <br /> Then writes the tree to file treefile-neworder.xml <br />
> -setOrder    Sets the Gravelius order of Segments, BranchingPoints and Buds starting from 1 so that parts belonging to an axis forking off from the present one get order one higher. If -iniVal value starts from that value. Then writes the tree to file treefile-neworder.xml <br />
> -setWoodyZero       Sets radius of all segments without foliage = 0, if -Hc value is set, only parts below Hc (in branches starting below Hc) are affected. <br />
> -maxCrownWidth      Prints out max horizontal distance from start of the branch of segments in branches (all segments order >= 2). If -foliage is set only segments with foliage are considered. <br />
> -collectSegmentLengths Prints sum (m) of segment lengths of all orders. If -order value is set, only for that order.  <br />
> -noBranchWhorls     Prints number of branch whorls (= branches leaving from exactly same spot, this should be developed to allow for small vertical distance between) <br />
> -meanBranchAngle    Prints mean branch angle and chord angle of branches.  <br />
> -setSapwoodArea     Sets heartwood radius so that sapwood area of a segment corresponds to SA above. Branch to stem 0.75, higher order forking 0.89 only so much of SA matched - cf. Sievanen et al. 2008. Writes result to file sapwoodarea.xml. If -correctR follows, sets R of segment to match sapwood requirement in case it doesn't. <br />
> -printLastSegInfo   Prints out information about last segment in each axis <br />
> -printStemDiameters    Prints out taper curve <br />
