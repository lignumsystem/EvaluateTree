#include <FractalDimension.h>
#include <Lacunarity.h>


//=================================================================================
// Lacunarity
//
//=================================================================================

// AC = C. Allain and M. Cloitre 1991 Characterizing the lacunarity of random and deterministic fractal sets.
// Physical Review A 44: 3552-3558

//CAC = Centered Allain and Coitre

//These definitions are from D. Da Silva 2008 Characterizing multiscale nature of plants using fractal geometry
// descriptors, application on light-interception modeling. PhD thesis, Informations, Structures, Systemes,
//Université Montpellier 2 Sciences et Techniques du Languedoc



void lacunarity(Tree<ETCfSegment,ETCfBud>& tree, int argc, char** argv, LGMdouble edge_len, string tree_file) {

  // Store the points on surfaces of segments (needle cylinder or woody cylinder
  // if no needles) to a list of points
  int no_points = 1000;
  string cla;
  if(ParseCommandLine(argc,argv,"-no_points",cla)) {
    no_points = atoi(cla.c_str());
  }

  list<Point> surfaces;                     //All points on the surfaces of shoot cylinders
  PointSurface ps(surfaces, no_points);

// Only top of an Axis covered with points --> must start with
// PropagateUp with the last segment the main axis
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

  int xN = static_cast<int>(xspan/edge_len) + 3;
  int yN = static_cast<int>(yspan/edge_len) + 3;
  int zN = static_cast<int>(zspan/edge_len) + 3;
  
  long int gridN = static_cast<long int>(xN)*static_cast<long int>(yN)*static_cast<long int>(zN);
  cout << "xN " << xN << "  yN " << yN << "  zN " << zN << endl;
  cout << "Number of grid points (voxels) is " << gridN << endl;

  //  bool occupy[xN][yN][zN];      //3D array indicating if this cell is occupied by the tree
  TMatrix3D<bool> occupy(xN,yN,zN);
  for(int i = 0; i < xN; i++) {
    for(int j = 0; j < yN; j++) {
      for(int k = 0; k < zN; k++) {
	occupy[i][j][k] = false;
      }
    }
  }

  long int no_occupied = 0;
  list<Point>::iterator I;
  for(I = surfaces.begin(); I != surfaces.end(); I++) {
    int ix = static_cast<int>(((*I).getX()-xmin)/edge_len) + 1;
    int iy = static_cast<int>(((*I).getY()-ymin)/edge_len) + 1;
    int iz = static_cast<int>(((*I).getZ()-zmin)/edge_len) + 1;
    if(!occupy[ix][iy][iz]) {
      occupy[ix][iy][iz] = true;
      no_occupied = no_occupied + 1;
    }
  }

  for(int k = 0; k < zN; k++) {
    cout << "------------------ " << k << " --------------" << endl;
    cout << "j ";
    for(int i = 0; i < xN; i++) 
      cout << i << " ";
    cout << endl;

    for(int j = 0; j < yN; j++) {
      cout << j << " ";
      for(int i = 0; i < xN; i++) {
	cout << occupy[i][j][k] << " ";
      }
      cout << endl;
    }
  }

  cout << "Occupied cells " << no_occupied << endl;

  //Then find the set of boxes that defines the hull of the tree

  vector<list<pair<int,int> > > hull;
  findHull(occupy, hull, xN, yN, zN);

  cout << hull.size() << endl;
 
  for(int k = 0; k < zN; k++){
    cout << endl;
    cout << " ==============================="  << k << " =========================" << endl;
    list<pair<int,int> >::iterator I;
    for(I = (hull[k]).begin(); I != (hull[k]).end(); I++) {
      cout << (*I).first << " " << (*I).second << endl;
    }
  }

  exit(0);

  int gbm = min(xN, min(yN,zN))   - 1;      //Gliding box maximum span

  cout << "Gliding box maximum span " << gbm << endl;

  string::size_type p = tree_file.find(".xml");
  string::size_type slash = tree_file.rfind("/");
  string f_name;
  if(slash == string::npos) {
    f_name = tree_file.substr(0,p);
  }
  else {
    f_name = tree_file.substr(slash+1,p-slash-1);
  }
  string AC_name = f_name+"-AC.dat";
  string CAC_name = f_name+"-CAC.dat";
  string CACR_name = f_name+"-CACR.dat";

  ofstream AC(AC_name.c_str(), ofstream::trunc);
  ofstream CAC(CAC_name.c_str(), ofstream::trunc);
  ofstream CACR(CACR_name.c_str(), ofstream::trunc);

  AC << "xN " << xN << "  yN " << yN << "  zN " << zN << endl;
  AC << "gbox L box" << endl;
  CAC << "xN " << xN << "  yN " << yN << "  zN " << zN << endl;
  CAC << "gbox L box" << endl;
  CACR << "xN " << xN << "  yN " << yN << "  zN " << zN << endl;
  CACR << "gbox L box" << endl;

  int gb_inc = 2;
  if(gbm > 100)
    gb_inc = 4;
  if(gb_inc > 150)
    gb_inc = 6;

  int gb = 1;

  while(gb < gbm) {
    cout << "gb " << gb << endl;
    int gb_max = gb*gb*gb;
    vector<int> values(gb_max + 1);
    for(int i = 0; i <= gb_max; i++) {
      values[i] = i;
    }
    vector<int> AC_freq(gb_max + 1,0);
    vector<int> CAC_freq(gb_max + 1,0);
    vector<int> CACR_freq(gb_max + 1,0);

    if(gb == 1) {
      AC_freq[0] = gridN - no_occupied;
      AC_freq[1] = no_occupied;
      CAC_freq[0] = 0;
      CAC_freq[1] = no_occupied;
      CACR_freq[0] = 0;
      CACR_freq[1] = no_occupied;
  
    
    } else {
      calculateACFrequencies(occupy,xN,yN,zN,gb,AC_freq);
      calculateCACFrequencies(occupy,xN,yN,zN,gb,CAC_freq);
      //    calculateCACFrequenciesReflect(occupy,xN,yN,zN,gb,CACR_freq);
      calculateCACFrequenciesX(occupy,xN,yN,zN,gb,hull,CACR_freq);

    }
  
    LGMdouble AC_maxpoints = static_cast<LGMdouble>((xN-gb+1)*(yN-gb+1)*(zN-gb+1));
    LGMdouble CAC_maxpoints = static_cast<LGMdouble>(no_occupied);

    vector<LGMdouble> AC_prob(gb_max + 1);
    vector<LGMdouble> CAC_prob(gb_max + 1);
    vector<LGMdouble> CACR_prob(gb_max + 1);

    for(int i = 0; i <= gb_max; i++) {
      AC_prob[i] = static_cast<LGMdouble>(AC_freq[i])/AC_maxpoints;
      CAC_prob[i] = static_cast<LGMdouble>(CAC_freq[i])/CAC_maxpoints;
      CACR_prob[i] = static_cast<LGMdouble>(CACR_freq[i])/CAC_maxpoints;
    }

    LGMdouble sum_mean = 0.0, sum_var = 0.0;
    for(int i = 0; i <= gb_max; i++) {
      sum_mean += static_cast<LGMdouble>(values[i]) * AC_prob[i];
      sum_var += pow(static_cast<LGMdouble>(values[i]),2.0) * AC_prob[i];
    }
  
    cout << "Allain - Cloitre Lacunarity, gliding box = " << gb << " voxels" << endl;
    cout << "m1 m2 " << sum_mean << " " << sum_var << endl;
    cout << "ACL = " << sum_var / pow(sum_mean, 2.0) << endl;
    AC << gb << " " << sum_var / pow(sum_mean, 2.0) << " "  << edge_len << endl;

    sum_mean = 0.0; sum_var = 0.0;
    for(int i = 0; i <= gb_max; i++) {
      sum_mean += static_cast<LGMdouble>(values[i]) * CAC_prob[i];
      sum_var += pow(static_cast<LGMdouble>(values[i]),2.0) * CAC_prob[i];
    }

    cout << "Centered Allain - Cloitre Lacunarity, gliding box = " << gb << " voxels" << endl;
    cout << "m1 m2 " << sum_mean << " " << sum_var << endl;
    cout << "CACL = " << sum_var / pow(sum_mean, 2.0) << endl;
    CAC << gb << " "  << sum_var / pow(sum_mean, 2.0) << " "  << edge_len << endl;

  
    sum_mean = 0.0; sum_var = 0.0;
    for(int i = 0; i <= gb_max; i++) {
      sum_mean += static_cast<LGMdouble>(values[i]) * CACR_prob[i];
      sum_var += pow(static_cast<LGMdouble>(values[i]),2.0) * CACR_prob[i];
    }
  
    cout << "Centered Allain - Cloitre Lacunarity Reflection, gliding box = " << gb << " voxels" << endl;
    cout << "m1 m2 " << sum_mean << " " << sum_var << endl;
    cout << "CACRL = " << sum_var / pow(sum_mean, 2.0) << endl;
    CACR << gb << " " << sum_var / pow(sum_mean, 2.0) << " "  << edge_len << endl;

    gb += gb_inc;

    cout << "======== " << gb << " " << gb_inc << endl;

    //   if(gb >= 7) {
    //     exit(0);
    //   }
  
  }     // while(gb < gbm) { ...

  AC.close();
  CAC.close();
  CACR.close();

}


void calculateACFrequencies(const TMatrix3D<bool>& occupy, int xN, int yN, 
			    int zN, int gb, vector<int>& AC_freq) {

  //Gliding box size is an odd number of cells
  int edge = (gb - 1) / 2;
  for(int i = edge; i < xN-edge; i++) {
    for(int j = edge; j < yN-edge; j++) {
      for(int k = edge; k < zN-edge; k++) {
	int n_occ = 0;
	for(int ig = i - edge; ig < i + edge; ig++) {
	  for(int jg = j - edge; jg < j + edge; jg++) {
	    for(int kg = k - edge; kg < k + edge; kg++) {
	      if(occupy[ig][jg][kg]) {
		n_occ++;
	      }
	    }
	  }
	}
	(AC_freq[n_occ])++;
      }
    }
  }
  

  return;
}

void calculateCACFrequencies(const TMatrix3D<bool>& occupy, int xN, int yN, 
			     int zN, int gb, vector<int>& CAC_freq) {

  //Gliding box size is an odd number of cells
  int edge = (gb - 1) / 2;

  for(int i = 0; i < xN; i++) {
    for(int j = 0; j < yN; j++) {
      for(int k = 0; k < zN; k++) {
	if(occupy[i][j][k]) {
	  int n_occ = 0;
	  for(int ig = i - edge; ig <= i + edge; ig++) {
	    if(ig < 0 || ig >= xN) {
	      continue;
	    }
	    for(int jg = j - edge; jg <= j + edge; jg++) {
	      if(jg < 0 || jg >= yN) {
		continue;
	      }
	      for(int kg = k - edge; kg <= k + edge; kg++) {
		if(kg < 0 || kg >= zN) {
		  continue;
		}
		if(occupy[ig][jg][kg]) {
		  n_occ++;
		}
	      } //kg
	    } //jg
	  } //ig
	  (CAC_freq[n_occ])++;
	}  //if(occupy ... 
      }
    }
  }
  
  return;
}

void calculateCACFrequenciesReflect(const TMatrix3D<bool>& occupy, int xN, int yN, 
				    int zN, int gb, vector<int>& CAC_freq) {

  //Gliding box size is an odd number of cells
  int edge = (gb - 1) / 2;

  for(int i = 0; i < xN; i++) {
    for(int j = 0; j < yN; j++) {
      for(int k = 0; k < zN; k++) {
	if(occupy[i][j][k]) {
	  int n_occ = 0;
	  for(int ig = i - edge; ig <= i + edge; ig++) {
	    for(int jg = j - edge; jg <= j + edge; jg++) {
	      for(int kg = k - edge; kg <= k + edge; kg++) {
		int igr = ig;
		if(igr < 0)
		  igr = xN + igr;
		if(igr > xN - 1)
		  igr = igr - xN;
		int jgr = jg;
		if(jgr < 0)
		  jgr = yN + jgr;
		if(jgr > yN - 1)
		  jgr = jgr - yN;
		int kgr = kg;
		if(kgr < 0)
		  kgr = zN + kgr;
		if(kgr > zN - 1)
		  kgr = kgr - zN;
		if(occupy[igr][jgr][kgr]) {
		  n_occ++;
		}
	      } //kg
	    } //jg
	  } //ig
	  (CAC_freq[n_occ])++;
	}  //if(occupy ... 
      }
    }
  }
  
  return;
}


void calculateCACFrequenciesX(const TMatrix3D<bool>& occupy, int xN, int yN, 
			      int zN, int gb, vector<list<pair<int,int> > >&  hull,
			      vector<int>& CACR_freq)
{

  //Gliding box size is an odd number of cells
  //  int edge = (gb - 1) / 2;
  //  vector<list<pair<int,int> > > limits;

  for(int k = 0; k < zN; k++) {
    for(int j = 0; j < yN; j++) {
      for(int i = 0; i < xN; i++) {
	if(occupy[i][j][k]) {
	  cout << k << " " << i << " " << j << " - ";
	}
      }
    }
    cout << endl;
  }

  for(int k = 0; k < zN; k++){
    cout << endl;
    cout << " ==============================="  << k << " =========================" << endl;
    list<pair<int,int> >::iterator I;
    for(I = (hull[k]).begin(); I != (hull[k]).end(); I++) {
      cout << (*I).first << " " << (*I).second << endl;
    }
  }
  
  exit(0); 

//   for(int i = 0; i < xN; i++) {
//     for(int j = 0; j < yN; j++) {
//       for(int k = 0; k < zN; k++) {
// 	if(occupy[i][j][k]) {
// 	  int n_occ = 0;
// 	  for(int ig = i - edge; ig <= i + edge; ig++) {
// 	    for(int jg = j - edge; jg <= j + edge; jg++) {
// 	      for(int kg = k - edge; kg <= k + edge; kg++) {


// 		int igr = ig;
// 		if(igr < 0)
// 		  igr = xN + igr;
// 		if(igr > xN - 1)
// 		  igr = igr - xN;
// 		int jgr = jg;
// 		if(jgr < 0)
// 		  jgr = yN + jgr;
// 		if(jgr > yN - 1)
// 		  jgr = jgr - yN;
// 		int kgr = kg;
// 		if(kgr < 0)
// 		  kgr = zN + kgr;
// 		if(kgr > zN - 1)
// 		  kgr = kgr - zN;
// 		if(occupy[igr][jgr][kgr]) {
// 		  n_occ++;
// 		}
// 	      } kg
// 	    } jg
// 	  } ig
// 	  (CAC_freq[n_occ])++;
// 	}  if(occupy ... 
//       }
//     }
//   }
  
  return;
}


void findHull(const TMatrix3D<bool>& occupy, vector<list<pair<int,int> > >& hull,
	      int xN, int yN, int zN) {
  for(int k = 0;k < zN; k++) {
    list<pair<int,int> > level_k;
    bool found = false;
    int jf = 0;
    int i = 0, j = 0;
    for(i = 0; i < xN; i++) {
      found = false;
      for(j = 0; j < yN; j++) {
	if(occupy[i][j][k]){
	  found = true;
	  break;
	}
      }
      if(found) {
	if(j > 0) {
	  jf = j - 1;
	}
	else {
	  jf = 0;
	}
	level_k.push_back(pair<int,int>(i,jf));
      }
      // end of j left to right

      found = false;
      jf = yN - 1;
      for(j = yN-1; j > -1; j--) {
	if(occupy[i][j][k]){
	  found = true;
	  break;
	}
      }
      if(found) {
	if(j < yN-1) {
	  jf = j + 1;
	}
	else {
	  jf = yN - 1;
	}
	level_k.push_back(pair<int,int>(i,jf));
      }
      // end of j right to left

    }   // i iteration

    for(j = 0; j < yN; j++) {
      found = false;
      bool found_k = false;
      int i_f = 0;
      list<pair<int,int> >::iterator I;
      for(i = 0; i < xN; i++) {
	if(occupy[i][j][k]){
	  found = true;
	  break;
	}
      }
      if(found) {
	if(i > 0) {
	  i_f = i - 1;
	}
	else {
	  i_f = 0;
	}
	found_k = false;
	for(I = level_k.begin(); I != level_k.end(); I++) {
	  if((*I).first == i_f && (*I).second == j) {
	    found_k = true;
	    break;
	  }
	}
	if(!found_k) {
	  level_k.push_back(pair<int,int>(i_f,j));
	}
      }
      // end of i left to right

      found = false;
      i_f = xN - 1;
      for(i = xN-1; i > -1; i--) {
	if(occupy[i][j][k]){
	  found = true;
	  break;
	}
      }
      if(found) {
	if(i < xN - 1) {
	  i_f = i + 1;
	}
	else {
	  i_f = xN - 1;
	}
	found_k = false;
	for(I = level_k.begin(); I != level_k.end(); I++) {
	  if((*I).first == i_f && (*I).second == j) {
	    found_k = true;
	    break;
	  }
	}
	if(!found_k) {
	  level_k.push_back(pair<int,int>(i_f,j));
	}
      }
      // end of i right to left
    }   // j iteration

    hull.push_back(level_k);
  }   // End of k loop

//   for(int k = 0; k < zN; k++){
//     cout << endl;
//     cout << " ==============================="  << k << " =========================" << endl;
//     list<pair<int,int> >::iterator I;
//     for(I = (hull[k]).begin(); I != (hull[k]).end(); I++) {
//       cout << (*I).first << " " << (*I).second << endl;
//     }
//   }


}





