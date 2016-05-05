/******************************************************************************
 *
 * File:           gridbmap.c
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
 *                 16 March 2016 Pavel Sakov
 *                 Moved most of the functionality related with spatial binary
 *                 trees to gridbmap.c. The code of this file is basically
 *                 common for "gridkmap" and "gridbmap" objects.
 *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "poly.h"
#include "gridnodes.h"
#include "gridmap.h"
#include "gridbmap.h"
#include "gridkmap.h"
#include "gucommon.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5

struct gridmap {
    void* map;
    int type;
    int sign;
};

/**
 */
gridmap* gridmap_build(int nce1, int nce2, double** gx, double** gy, int type)
{
    gridmap* gm = malloc(sizeof(gridmap));

    gm->type = type;
    if (gm->type == GRIDMAP_TYPE_BINARY)
        gm->map = gridbmap_build(nce1, nce2, gx, gy);
    else if (gm->type == GRIDMAP_TYPE_KDTREE)
        gm->map = gridkmap_build(nce1, nce2, gx, gy);
    else
        gu_quit("grid map type = %d: unknown type", type);
    gm->sign = 0;

    return gm;
}

/**
 */
gridmap* gridmap_build2(gridnodes* gn)
{
    gridmap* gm = malloc(sizeof(gridmap));
    int nce1 = gridnodes_getnce1(gn);
    int nce2 = gridnodes_getnce2(gn);
    double** gx = gridnodes_getx(gn);
    double** gy = gridnodes_gety(gn);
    int type = gridnodes_getmaptype(gn);

    gm->type = type;
    if (gm->type == GRIDMAP_TYPE_BINARY)
        gm->map = gridbmap_build(nce1, nce2, gx, gy);
    else if (gm->type == GRIDMAP_TYPE_KDTREE)
        gm->map = gridkmap_build(nce1, nce2, gx, gy);
    else
        gu_quit("grid map type = %d: unknown type", type);
    gm->sign = 0;

    return gm;
}

/**
 */
void gridmap_destroy(gridmap* gm)
{
    if (gm->type == GRIDMAP_TYPE_BINARY)
        gridbmap_destroy(gm->map);
    else if (gm->type == GRIDMAP_TYPE_KDTREE)
        gridkmap_destroy(gm->map);

    free(gm);
}

/**
 */
int gridmap_xy2ij(gridmap* gm, double x, double y, int* i, int* j)
{
    int success = 0;

    *i = -1;
    *j = -1;
    if (!isfinite(x + y))
        return success;

    if (gm->type == GRIDMAP_TYPE_BINARY)
        success = gridbmap_xy2ij(gm->map, x, y, i, j);
    else if (gm->type == GRIDMAP_TYPE_KDTREE)
        success = gridkmap_xy2ij(gm->map, x, y, i, j);

    return success;
}

/** Calculates (x,y) coordinates for a point within the grid with
 * specified fractional indices (i,j).
 *
 * The transformation used to compute the coords is a forward
 * tetragonal bilinear texture mapping.
 *
 * @param gm a tree structure returned from xytoij_init
 * @param fi I indice value
 * @param fj J indice value
 * @param x Pointer to returned X coordinate
 * @param y Pointer to returned Y coordinate
 * @return non-zero if successful
 */
int gridmap_fij2xy(gridmap* gm, double fi, double fj, double* x, double* y)
{
    int status = 1;
    double** gx = NULL;
    double** gy = NULL;
    int nce1 = -1, nce2 = -1;
    int i, j;
    double u, v;
    double a, b, c, d, e, f, g, h;

    if (gm->type == GRIDMAP_TYPE_BINARY) {
        gx = gridbmap_getxnodes(gm->map);
        gy = gridbmap_getynodes(gm->map);
        nce1 = gridbmap_getnce1(gm->map);
        nce2 = gridbmap_getnce2(gm->map);
    } else if (gm->type == GRIDMAP_TYPE_KDTREE) {
        gx = gridkmap_getxnodes(gm->map);
        gy = gridkmap_getynodes(gm->map);
        nce1 = gridkmap_getnce1(gm->map);
        nce2 = gridkmap_getnce2(gm->map);
    }

    /*
     * Trim I to range 0 to nce1 
     */
    if (fi < 0) {
        fi = 0.0;
        status = 0;
    }

    if (fi > nce1) {
        fi = nce1 - EPS;
        status = 0;
    }

    /*
     * Trim J to range 0 to nce2 
     */
    if (fj < 0) {
        fj = 0.0;
        status = 0;
    }

    if (fj > nce2) {
        fj = nce2 - EPS;
        status = 0;
    }

    i = (int) fi;
    j = (int) fj;

    u = fi - i;
    v = fj - j;

    if (u == 0.0 && v == 0.0) {
        *x = gx[j][i];
        *y = gy[j][i];
    } else if (u == 0.0) {
        *x = gx[j + 1][i] * v + gx[j][i] * (1.0 - v);
        *y = gy[j + 1][i] * v + gy[j][i] * (1.0 - v);
    } else if (v == 0.0) {
        *x = gx[j][i + 1] * u + gx[j][i] * (1.0 - u);
        *y = gy[j][i + 1] * u + gy[j][i] * (1.0 - u);
    } else {
        a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        b = gx[j][i + 1] - gx[j][i];
        c = gx[j + 1][i] - gx[j][i];
        d = gx[j][i];
        e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        f = gy[j][i + 1] - gy[j][i];
        g = gy[j + 1][i] - gy[j][i];
        h = gy[j][i];

        *x = a * u * v + b * u + c * v + d;
        *y = e * u * v + f * u + g * v + h;
    }

    return status;
}

/** Calculates the branch of sqrt() to be taken in gridbmap_xy2fij(). Has to be
 * called only once for a grid.
 * 
 * @param gm Grid map
 * @param x X coordinate
 * @param y Y coordinate
 * @return 1 or -1 if successful; 0 otherwhile
 */
static int calc_branch(gridmap* gm, double x, double y)
{
    int i, j;
    double** gx = NULL;
    double** gy = NULL;
    int sign = 1;
    double error[2];

    /*
     * normally one tries xytoij() before calling calc_branch() 
     */
    if (gridmap_xy2ij(gm, x, y, &i, &j) == 0)
        return 0;               /* failed */

    if (gm->type == GRIDMAP_TYPE_BINARY) {
        gx = gridbmap_getxnodes(gm->map);
        gy = gridbmap_getynodes(gm->map);
    } else if (gm->type == GRIDMAP_TYPE_KDTREE) {
        gx = gridkmap_getxnodes(gm->map);
        gy = gridkmap_getynodes(gm->map);
    }

    {
        double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        double b = gx[j][i + 1] - gx[j][i];
        double c = gx[j + 1][i] - gx[j][i];
        double d = gx[j][i];
        double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        double f = gy[j][i + 1] - gy[j][i];
        double g = gy[j + 1][i] - gy[j][i];
        double h = gy[j][i];

        double A = a * f - b * e;

        double B, C;
        int k;

        /*
         * normally one checks A before calling calc_branch() 
         */
        if (fabs(A) < EPS_ZERO)
            return 0;           /* failed */

        B = e * x - a * y + a * h - d * e + c * f - b * g;
        C = g * x - c * y + c * h - d * g;

        for (k = 0; k < 2; ++k) {
            double u = (-B + sign * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
            double v_denom = a * u + c;
            double v = (fabs(v_denom) < EPS_ZERO) ? (y - f * u - h) / (e * u + g) : (x - b * u - d) / v_denom;

            error[k] = 0.0;

            if (u < 0.0)
                error[k] -= u;
            else if (u > 1.0)
                error[k] += (u - 1.0);
            if (v < 0.0)
                error[k] -= v;
            else if (v > 1.0)
                error[k] += (v - 1.0);

            sign = -1;
        }
    }

    if (error[0] < error[1])
        return 1;
    return -1;
}

/** Calculates (x,y) coordinates for a point within a numerical grid specified
 * by fractional indices (i,j).
 *
 * The transformation used to compute the indices is an inverse
 * tetragonal bilinear texture mapping.
 *
 * At the moment we assume that either there is only one grid in the model or
 * all the grids use a uniform branch in xy2fij() convertion.
 *
 * @param gm Grid map
 * @param x X coordinate
 * @param y Y coordinate
 * @param x Pointer to returned fractional I index
 * @param y Pointer to returned fractional J index
 * @return 1 if successful, 0 otherwise
 */
int gridmap_xy2fij(gridmap* gm, double x, double y, double* fi, double* fj)
{
    int i, j;
    double** gx = NULL;
    double** gy = NULL;

    if (gridmap_xy2ij(gm, x, y, &i, &j) == 0)
        return 0;               /* failed */

    if (gm->type == GRIDMAP_TYPE_BINARY) {
        gx = gridbmap_getxnodes(gm->map);
        gy = gridbmap_getynodes(gm->map);
    } else if (gm->type == GRIDMAP_TYPE_KDTREE) {
        gx = gridkmap_getxnodes(gm->map);
        gy = gridkmap_getynodes(gm->map);
    }

    {
        double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        double b = gx[j][i + 1] - gx[j][i];
        double c = gx[j + 1][i] - gx[j][i];
        double d = gx[j][i];
        double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        double f = gy[j][i + 1] - gy[j][i];
        double g = gy[j + 1][i] - gy[j][i];
        double h = gy[j][i];

        double A = a * f - b * e;
        double B = e * x - a * y + a * h - d * e + c * f - b * g;
        double C = g * x - c * y + c * h - d * g;

        double u, v, d1, d2;

        if (fabs(A) < EPS_ZERO)
            u = -C / B * (1.0 + A * C / B / B);
        else {
            if (gm->sign == 0) {
                gm->sign = calc_branch(gm, x, y);
                if (gm->sign == 0)
                    return 0;   /* failed */
            }
            u = (-B + gm->sign * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
        }
        d1 = a * u + c;
        d2 = e * u + g;
        v = (fabs(d2) > fabs(d1)) ? (y - f * u - h) / d2 : (x - b * u - d) / d1;

        if (u < 0.0)
            u = 0.0;
        else if (u >= 1.0)
            u = 1.0 - EPS;
        if (v < 0.0)
            v = 0.0;
        else if (v >= 1.0)
            v = 1.0 - EPS;

        *fi = i + u;
        *fj = j + v;
    }

    return 1;
}

/**
 */
int gridmap_getnce1(gridmap* gm)
{
    int nce1 = -1;

    if (gm->type == GRIDMAP_TYPE_BINARY)
        nce1 = gridbmap_getnce1(gm->map);
    else if (gm->type == GRIDMAP_TYPE_KDTREE)
        nce1 = gridkmap_getnce1(gm->map);

    return nce1;
}

/**
 */
int gridmap_getnce2(gridmap* gm)
{
    int nce2 = -1;

    if (gm->type == GRIDMAP_TYPE_BINARY)
        nce2 = gridbmap_getnce2(gm->map);
    else if (gm->type == GRIDMAP_TYPE_KDTREE)
        nce2 = gridkmap_getnce2(gm->map);

    return nce2;
}
