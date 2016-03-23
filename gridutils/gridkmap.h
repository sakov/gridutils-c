/******************************************************************************
 *
 * File:           gridkmap.h
 *  
 * Created:        16 March 2016
 *  
 * Author:         Pavel Sakov
 *                 BoM
 *  
 * Purpose:        Calculates transformations between physical and index
 *                 space for a numerical grid using kd-tree
 *
 * Revisions:      
 *
 *****************************************************************************/

#if !defined(_GRIDKMAP_H)
#define _GRIDKMAP_H

struct gridkmap;
typedef struct gridkmap gridkmap;

gridkmap* gridkmap_build(int nce1, int nce2, double** gx, double** gy);
void gridkmap_destroy(gridkmap* gm);
int gridkmap_xy2ij(gridkmap* gm, double x, double y, int* i, int* j);
int gridkmap_getnce1(gridkmap* gm);
int gridkmap_getnce2(gridkmap* gm);
double** gridkmap_getxnodes(gridkmap* gm);
double** gridkmap_getynodes(gridkmap* gm);

#endif
