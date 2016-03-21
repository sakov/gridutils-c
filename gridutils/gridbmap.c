/******************************************************************************
 *
 * File:           gridmap.c
 *  
 * Created:        Fri Feb 19 11:36:24 EST 1992 (as xytoij.c)
 *  
 * Author:         Daniel Delbourgo/Stephen Walker
 *                 CSIRO Division of Oceanography
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
 *                 left-handed grids.
 *
 *                 April 2002 Pavel Sakov
 *                 Major mods to handle topologically non-rectangular grids
 *                 Renamed to gridmap.c
 *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "poly.h"
#include "gridbmap.h"
#include "gucommon.h"

#define EPS_COMPACT 1.0e-10

typedef struct subgrid {
    gridbmap* gmap;              /* gridf map this subgrid belongs to */
    poly* bound;                /* boundary polygon */
    int mini;                   /* minimal i index within the subgrid */
    int maxi;                   /* maximal i index within the subgrid */
    int minj;                   /* minimal j index within the subgrid */
    int maxj;                   /* maximal j index within the subgrid */
    struct subgrid* half1;      /* child 1 */
    struct subgrid* half2;      /* child 2 */
} subgrid;

struct gridbmap {
    poly* bound;                /* boundary polygon */
    subgrid* trunk;             /* binary tree trunk */
    int nleaves;                /* for debugging purposes */
    int nce1;                   /* number of cells in e1 direction */
    int nce2;                   /* number of cells in e2 direction */
    double** gx;                /* reference to array of X coords
                                 * [nce1+1][nce2+1] */
    double** gy;                /* reference to array of Y coords
                                 * [nce1+1][nce2+1] */
};

/** Creates a subgrid.
 * @param gm Grid map
 * @param pl Boundary polygon for the subgrid
 * @param mini Minimal i index within the subgrid
 * @param maxi Maximal i index within the subgrid
 * @param minj Minimal j index within the subgrid
 * @param maxj Maximal j index within the subgrid
 * @return Subgrid
 */
static subgrid* subgrid_create(gridbmap* gm, poly* pl, int i1, int i2, int j1, int j2)
{
    subgrid* l = malloc(sizeof(subgrid));

    double** gx = gm->gx;
    double** gy = gm->gy;
    int n = pl->n;
    double x, y;
    int i = i1;                 /* to eliminate warning */
    int j, ii;

    l->bound = pl;
    l->gmap = gm;
    l->mini = INT_MAX;
    l->maxi = INT_MIN;
    l->minj = INT_MAX;
    l->maxj = INT_MIN;
    l->half1 = NULL;
    l->half2 = NULL;

    if (n == 0)
        return l;

    x = pl->x[0];
    y = pl->y[0];

    for (j = j1; j <= j2; ++j)
        for (i = i1; i <= i2; ++i)
            if (x == gx[j][i] && y == gy[j][i])
                goto aftersearch;

  aftersearch:

    if (j > j2)
        gu_quit("subgrid_create(): boundary vertex not in the grid");

    l->mini = i;
    l->maxi = i;
    l->minj = j;
    l->maxj = j;

    for (ii = 1; ii < n; ++ii) {
        x = pl->x[ii];
        y = pl->y[ii];

        if (i > i1 && x == gx[j][i - 1] && y == gy[j][i - 1])
            i--;
        else if (i < i2 && x == gx[j][i + 1] && y == gy[j][i + 1])
            i++;
        else if (j > j1 && x == gx[j - 1][i] && y == gy[j - 1][i])
            j--;
        else if (j < j2 && x == gx[j + 1][i] && y == gy[j + 1][i])
            j++;
        else if (x == gx[j][i] && y == gy[j][i])
            continue;
        else
            gu_quit("subgrid_create(): boundary vertex not in the grid");

        if (l->mini > i)
            l->mini = i;
        if (l->maxi < i)
            l->maxi = i;
        if (l->minj > j)
            l->minj = j;
        if (l->maxj < j)
            l->maxj = j;
    }

    return l;
}

/** Destroys a subgrid.
 * @param l Subgrid
 */
static void subgrid_destroy(subgrid* l)
{
    if (l->half1 != NULL)
        subgrid_destroy(l->half1);
    if (l->half2 != NULL)
        subgrid_destroy(l->half2);
    poly_destroy(l->bound);
    free(l);
}

/** Cuts boundary polygon in two. 
 * The cut goes either horizontally ([fixed][changes]) or vertically 
 * ([changes][fixed]) in index space; the physical nodes are given by
 * input double arrays; first two intersections of the cutting polyline
 * with the polygon are used to form the new polygons.
 * @param pl Original polygon
 * @param gx Array of x cell corner coordinates
 * @param gy Array of y cell corner coordinates
 * @param horiz flag: 1 for horizontal cut; 0 otherwise
 * @param index Value of "fixed" index
 * @param start Start value of "variable" index
 * @param end End value of "variable" index
 * @param pl1 Output polygon 1
 * @param pl1 Output polygon 2
 */
static void cut_boundary(poly* pl, double** gx, double** gy, int horiz, int index, int start, int end, poly** pl1, poly** pl2)
{
    int n = pl->n;
    int i = -1;
    int i1 = -1;                /* array index of the first intersection */
    int i2 = -1;                /* array index of the second intersection */
    int ii1 = -1;               /* polygon index of the first intersection */
    int ii2 = -1;               /* polygon index of the second intersection */
    int tmp = -1;

    /*
     * if the polygon has been explicitely closed, ignore the last point 
     */
    if (poly_isclosed(pl, 1.0e-15))
        n--;

    if (horiz) {                /* horizontal cut */
        /*
         * find first intersection 
         */
        for (i = start; i < end; ++i) {
            ii1 = poly_findindex(pl, gx[index][i], gy[index][i]);
            if (ii1 < 0)
                /*
                 * this node does not belong to the boundary 
                 */
                continue;
            /*
             * this node belongs to the boundary polygon To accept it as a
             * valid intersection, we must ensure that the very next node is 
             * either inside the polygon or belongs to the boundary but is
             * not the next boundary node. 
             */
            tmp = poly_findindex(pl, gx[index][i + 1], gy[index][i + 1]);
            if ((tmp < 0 && poly_containspoint(pl, gx[index][i + 1], gy[index][i + 1])) || (tmp >= 0 && abs(tmp - ii1) > 1 && abs(tmp - ii1) < n - 1))
                break;
        }

        /*
         * how the for-cycle ended 
         */
        if (i < end)            /* ok */
            i1 = i;
        else                    /* no intersection found */
            return;

        /*
         * find second intersection start from the node next to the first
         * intersection 
         */
        for (i = i1 + 1; i <= end; ++i) {
            ii2 = poly_findindex(pl, gx[index][i], gy[index][i]);
            if (ii2 >= 0)
                /*
                 * this node must be inside the boundary polygon -- skip 
                 */
                break;
        }

        if (ii2 < 0)
            /*
             * no intersection found 
             */
            return;
        else
            /*
             * ok 
             */
            i2 = i;

        /*
         * we found all necessary details, now form the new polygons 
         */
        *pl1 = poly_create();
        *pl2 = poly_create();

        /*
         * add the portion of perimeter 
         */
        for (i = ii1; i != ii2; i = (i + 1) % n)
            poly_addpoint(*pl1, pl->x[i], pl->y[i]);
        /*
         * add the cutting section 
         */
        for (i = i2; i > i1; --i)
            poly_addpoint(*pl1, gx[index][i], gy[index][i]);

        /*
         * add the portion of perimeter 
         */
        for (i = ii2; i != ii1; i = (i + 1) % n)
            poly_addpoint(*pl2, pl->x[i], pl->y[i]);
        /*
         * add the cutting section 
         */
        for (i = i1; i < i2; ++i)
            poly_addpoint(*pl2, gx[index][i], gy[index][i]);

    } else {                    /* vertical cut */
        for (i = start; i < end; ++i) {
            ii1 = poly_findindex(pl, gx[i][index], gy[i][index]);
            if (ii1 < 0)
                continue;
            tmp = poly_findindex(pl, gx[i + 1][index], gy[i + 1][index]);
            if ((tmp < 0 && poly_containspoint(pl, gx[i + 1][index], gy[i + 1][index])) || (tmp >= 0 && abs(tmp - ii1) > 1 && abs(tmp - ii1) < n - 1))
                break;
        }

        if (i < end)
            i1 = i;
        else
            return;

        for (i = i1 + 1; i <= end; ++i) {
            ii2 = poly_findindex(pl, gx[i][index], gy[i][index]);
            if (ii2 >= 0)
                break;
        }

        if (ii2 < 0)
            return;
        else
            i2 = i;

        *pl1 = poly_create();
        *pl2 = poly_create();

        for (i = ii1; i != ii2; i = (i + 1) % n)
            poly_addpoint(*pl1, pl->x[i], pl->y[i]);
        for (i = i2; i > i1; --i)
            poly_addpoint(*pl1, gx[i][index], gy[i][index]);

        for (i = ii2; i != ii1; i = (i + 1) % n)
            poly_addpoint(*pl2, pl->x[i], pl->y[i]);
        for (i = i1; i < i2; ++i)
            poly_addpoint(*pl2, gx[i][index], gy[i][index]);
        /*
         * There used to be closure of the polylines here:
         *        poly_close(*pl1);
         *        poly_close(*pl2);
         * While this does not harm, it causes the boundary of a quadrilateral
         * to contain 5 points rather than 4, which caused allocation for 8
         * points... Anyway, considering the way search in poly_containspoint()
         * works, this is not necessary.
         */
    }
}

/* Divides a subgrid in two.
 * @param sg The subgrid to divide
 * @param subgrid1 Output subgrid 1
 * @param subgrid2 Output subgrid 2
 * @param gx Array of x cell corner coordinates
 * @param gy Array of y cell corner coordinates
 */
static void subgrid_divide(subgrid* sg, subgrid** sg1, subgrid** sg2)
{
    poly* pl1 = NULL;
    poly* pl2 = NULL;
    gridbmap* gm = sg->gmap;
    int index;

    if ((sg->maxi <= sg->mini + 1) && (sg->maxj <= sg->minj + 1)) {
        *sg1 = *sg2 = NULL;
        return;
    }

    if (sg->maxi - sg->mini > sg->maxj - sg->minj) {
        /*
         * divide "vertically" 
         */
        index = (sg->mini + sg->maxi) / 2;
        cut_boundary(sg->bound, gm->gx, gm->gy, 0, index, sg->minj, sg->maxj, &pl1, &pl2);
    } else {
        /*
         * divide "horizontally" 
         */
        index = (sg->minj + sg->maxj) / 2;
        cut_boundary(sg->bound, gm->gx, gm->gy, 1, index, sg->mini, sg->maxi, &pl1, &pl2);
    }

    if (pl1 == NULL || pl2 == NULL)
        gu_quit("dividesubgrid(): could not cut the boundary");

    *sg1 = subgrid_create(gm, pl1, sg->mini, sg->maxi, sg->minj, sg->maxj);
    *sg2 = subgrid_create(gm, pl2, sg->mini, sg->maxi, sg->minj, sg->maxj);
}

/**
 */
static void gridbmap_subdivide(gridbmap* gm, subgrid* sg)
{
    subgrid* sg1 = NULL;
    subgrid* sg2 = NULL;

    subgrid_divide(sg, &sg1, &sg2);

    if (sg1 != NULL) {
        sg->half1 = sg1;
        ++(gm->nleaves);
        gridbmap_subdivide(gm, sg1);
    }
    if (sg2 != NULL) {
        gridbmap_subdivide(gm, sg2);
        sg->half2 = sg2;
        ++(gm->nleaves);
    }
    poly_compact(sg->bound, EPS_COMPACT);
}

/** Builds a grid map structure to facilitate conversion from coordinate
 * to index space.
 *
 * @param gx array of X coordinates (of size (nce1+1)*(nce2+1))
 * @param gy array of Y coordinates (of size (nce1+1)*(nce2+1))
 * @param nce1 number of cells in e1 direction
 * @param nce2 number of cells in e2 direction
 * @return a map tree to be used by xy2ij
 */
gridbmap* gridbmap_build(int nce1, int nce2, double** gx, double** gy)
{
    gridbmap* gm = malloc(sizeof(gridbmap));
    poly* bound;
    subgrid* trunk;

    gm->nce1 = nce1;
    gm->nce2 = nce2;
    gm->gx = gx;
    gm->gy = gy;

    bound = poly_formbound(nce1, nce2, gx, gy);
    trunk = subgrid_create(gm, bound, 0, nce1, 0, nce2);

    gm->bound = bound;
    gm->trunk = trunk;
    gm->nleaves = 1;

    gridbmap_subdivide(gm, trunk);       /* recursive */

    return gm;
}

/** Destroys a grid map.
 * @param gm Grid map
 */
void gridbmap_destroy(gridbmap* gm)
{
    subgrid_destroy(gm->trunk);
    free(gm);
}

/** Calculates indices (i,j) of a grid cell containing point (x,y).
 *
 * @param gm Grid map
 * @param x X coordinate
 * @param y Y coordinate
 * @param i pointer to returned I indice value
 * @param j pointer to returned J indice value
 * @return 1 if successful, 0 otherwhile
 */
int gridbmap_xy2ij(gridbmap* gm, double x, double y, int* i, int* j)
{
    subgrid* sg = gm->trunk;

    /*
     * check if point is in grid outline 
     */
    sg = gm->trunk;
    if (!poly_containspoint(sg->bound, x, y))
        return 0;

    /*
     * do the full search 
     */
    while (sg->half1 != NULL) {
        /*
         * Test on the point being within the boundary polyline is the most
         * expensive part of the mapping; therefore, perform it in a branch
         * that contains a smaller (number of points-wise) polyline.
         */
        if (sg->half1->bound->n <= sg->half2->bound->n)
            sg = (poly_containspoint(sg->half1->bound, x, y)) ? sg->half1 : sg->half2;
        else
            sg = (poly_containspoint(sg->half2->bound, x, y)) ? sg->half2 : sg->half1;
    }

    *i = sg->mini;
    *j = sg->minj;

    return 1;
}

/**
 */
int gridbmap_getnce1(gridbmap* gm)
{
    return gm->nce1;
}

/**
 */
int gridbmap_getnce2(gridbmap* gm)
{
    return gm->nce2;
}

/**
 */
void gridbmap_getextent(gridbmap* gm, double* xmin, double* xmax, double* ymin, double* ymax)
{
    extent* e = &gm->trunk->bound->e;

    *xmin = e->xmin;
    *xmax = e->xmax;
    *ymin = e->ymin;
    *ymax = e->ymax;
}

/**
 */
double** gridbmap_getxnodes(gridbmap* gm)
{
    return gm->gx;
}

/**
 */
double** gridbmap_getynodes(gridbmap* gm)
{
    return gm->gy;
}
