#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "kdtree.h"
#include "poly.h"
#include "gridkmap.h"
#include "gucommon.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5

struct gridkmap {
    int nce1;                   /* number of cells in e1 direction */
    int nce2;                   /* number of cells in e2 direction */
    double** gx;                /* reference to array of X coords
                                 * [nce2+1][nce1+1] */
    double** gy;                /* reference to array of Y coords
                                 * [nce2+1][nce1+1] */
    kdtree* kd;                 /* kd tree with grid nodes */
};

/** Builds a grid map structure to facilitate conversion from coordinate
 * to index space.
 *
 * @param gx array of X coordinates (of size (nce1+1)*(nce2+1))
 * @param gy array of Y coordinates (of size (nce1+1)*(nce2+1))
 * @param nce1 number of cells in e1 direction
 * @param nce2 number of cells in e2 direction
 * @return a map tree to be used by xy2ij
 */
gridkmap* gridkmap_build(int nce1, int nce2, double** gx, double** gy)
{
    gridkmap* gm = malloc(sizeof(gridkmap));
    int* ids = NULL;
    int n, ii;

    gm->nce1 = nce1;
    gm->nce2 = nce2;
    gm->gx = gx;
    gm->gy = gy;

    nce1++;
    nce2++;
    n = nce1 * nce2;
    ids = malloc(n * sizeof(int));
    for (ii = 0; ii < n; ++ii)
        ids[ii] = ii;
    shuffle(n, ids);

    gm->kd = kd_create(2);
    for (ii = 0; ii < n; ++ii) {
        int id = ids[ii];
        double pos[2];
        int i, j;

        i = id % nce1;
        j = id / nce1;
        pos[0] = gx[j][i];
        pos[1] = gy[j][i];

        if (!isnan(pos[0]) && !isnan(pos[1]))
            kd_insert(gm->kd, pos, (void*) (long int) id);
    }

    free(ids);

    return gm;
}

/**
 */
void gridkmap_destroy(gridkmap* gm)
{
    kd_free(gm->kd);
    free(gm);
}

/**
 */
int gridkmap_xy2ij(gridkmap* gm, double x, double y, int* iout, int* jout)
{
    double pos[2] = { x, y };
    kdnode* nearest = kd_nearest(gm->kd, pos);
    size_t id = (size_t) kdnode_data(nearest);
    poly* p = poly_create();
    int success = 0;
    int i, j, i1, i2, j1, j2;

    *iout = -1;
    *jout = -1;

    j = id / (gm->nce1 + 1);
    i = id % (gm->nce1 + 1);

    i1 = (i > 0) ? i - 1 : i;
    i2 = (i < gm->nce1) ? i + 1 : i;
    j1 = (j > 0) ? j - 1 : j;
    j2 = (j < gm->nce2) ? j + 1 : j;

    /*
     * TODO?: looks a bit heavyweight
     */
    for (j = j1; j <= j2 - 1; ++j)
        for (i = i1; i <= i2 - 1; ++i) {
            poly_addpoint(p, gm->gx[j][i], gm->gy[j][i]);
            poly_addpoint(p, gm->gx[j][i + 1], gm->gy[j][i + 1]);
            poly_addpoint(p, gm->gx[j + 1][i + 1], gm->gy[j + 1][i + 1]);
            poly_addpoint(p, gm->gx[j + 1][i], gm->gy[j + 1][i]);
            poly_addpoint(p, gm->gx[j][i], gm->gy[j][i]);

            if (poly_containspoint(p, x, y)) {
                success = 1;
                *iout = i;
                *jout = j;
                goto finish;
            }
            poly_clear(p);
        }

  finish:

    poly_destroy(p);

    return success;
}

/**
 */
int gridkmap_getnce1(gridkmap* gm)
{
    return gm->nce1;
}

/**
 */
int gridkmap_getnce2(gridkmap* gm)
{
    return gm->nce2;
}

/**
 */
double** gridkmap_getxnodes(gridkmap* gm)
{
    return gm->gx;
}

/**
 */
double** gridkmap_getynodes(gridkmap* gm)
{
    return gm->gy;
}
