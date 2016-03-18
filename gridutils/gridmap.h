/******************************************************************************
 *
 * File:           gridmap.h
 *  
 * Created:        16 March 2016
 *  
 * Author:         Pavel Sakov
 *                 BoM
 *
 *                 A fair bit of this code is based on the original code by
 *                 Daniel Delbourgo with contributions from Stephen Walker and
 *                 Jason Waring (CSIRO Division of Oceanography, 1992-1998).
 *  
 * Purpose:        Calculates transformations between physical and index
 *                 space within a numerical grid. Mapping xy->ij can now
 *                 be conducted by one of two algorithms: via rendering grid
 *                 into a spatial binary tree and via kd-tree with grid nodes.
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_GRIDMAP_H)
#define _GRIDMAP_H

#define GRIDMAP_TYPE_BINARY 0
#define GRIDMAP_TYPE_KDTREE 1
#define GRIDMAP_TYPE_DEF GRIDMAP_TYPE_BINARY

extern int gridmaptype;

struct gridmap;
typedef struct gridmap gridmap;

gridmap* gridmap_build(int nce1, int nce2, double** gx, double** gy);
void gridmap_destroy(gridmap* gm);
int gridmap_fij2xy(gridmap* gm, double fi, double fj, double* x, double* y);
int gridmap_xy2ij(gridmap* gm, double x, double y, int* i, int* j);
int gridmap_xy2fij(gridmap* gm, double x, double y, double* fi, double* fj);
int gridmap_getnce1(gridmap* gm);
int gridmap_getnce2(gridmap* gm);

#endif
