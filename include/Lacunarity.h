#ifndef LACUNARITY_H
#define LACUNARITY_H

#include <TMatrix3D.h>
#include <EvaluateTree.h>
#include <FractalDimension.h>



// AC = C. Allain and M. Cloitre 1991 Characterizing the lacunarity of random and deterministic fractal sets.
// Physical Review A 44: 3552-3558

//CAC = Centered Allain and Cloitre

//These definitions are from D. Da Silva 2008 Characterizing multiscale nature of plants using fractal geometry
// descriptors, application on light-interception modeling. PhD thesis, Informations, Structures, Systemes,
//Université Montpellier 2 Sciences et Techniques du Languedoc



void findHull(const TMatrix3D<bool>& occupy, vector<list<pair<int,int> > >& hull,
	      int xN, int yN, int zN);

void calculateACFrequencies(const TMatrix3D<bool>& occupy, int xN, int yN, 
			    int zN, int gb, vector<int>& AC_freq);

void calculateCACFrequencies(const TMatrix3D<bool>& occupy, int xN, int yN, 
			     int zN, int gb, vector<int>& CAC_freq);

void calculateCACFrequenciesReflect(const TMatrix3D<bool>& occupy, int xN, int yN, 
				    int zN, int gb, vector<int>& CAC_freq);

void calculateCACFrequenciesX(const TMatrix3D<bool>& occupy, int xN, int yN, 
			      int zN, int gb, vector<list<pair<int,int> > >&  hull,
			      vector<int>& CACR_freq);



void lacunarity(Tree<ETCfSegment,ETCfBud>& tree, int argc, char** argv, LGMdouble edge_len, string tree_file);



#endif
