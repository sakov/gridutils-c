/******************************************************************************
 *
 * File:           gridbmap.h
 *  
 * Created:        Thu Jan 23 14:00:00 EST 1997
 *  
 * Author:         Daniel Delbourgo/Stephen Walker/Jason Waring
 *                 CSIRO Marine Research
 *  
 * Purpose:        Calculates transformations between physical and index
 *                 space within a numerical grid
 *
 * Revisions:      110898 JRW
 *                 Added IJtoXY conversion and
 *                 fractional XYtoIJ and IJtoXY conversion
 *
 *                 2000 Pavel Sakov
 *                 Added branch calculation to handle both right- and
 *                 left-handed grids (calc_branch()).
 *
 *                 April 2002 Pavel Sakov
 *                 Major mods to handle topologically non-rectangular grids
 *
 *                 16032016 Pavel Sakov
 *                 -- Changed the base name from "gridmap" to "gridbmap".
 *                 -- Moved some stuff to gridmap.h.
 *
 *****************************************************************************/

#if !defined(_GRIDBMAP_H)
#define _GRIDBMAP_H

struct gridbmap;
typedef struct gridbmap gridbmap;

gridbmap* gridbmap_build(int nce1, int nce2, double** gx, double** gy);
void gridbmap_destroy(gridbmap* gm);
int gridbmap_xy2ij(gridbmap* gm, double x, double y, int* i, int* j);
int gridbmap_getnce1(gridbmap* gm);
int gridbmap_getnce2(gridbmap* gm);
double** gridbmap_getxnodes(gridbmap* gm);
double** gridbmap_getynodes(gridbmap* gm);
void gridbmap_getextent(gridbmap* gm, double* xmin, double* xmax, double* ymin, double* ymax);

#endif
