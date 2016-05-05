/******************************************************************************
 *
 *  File:           xy2ij.c
 *  
 *  Created         22/01/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Given a numerical grid, converts point coordinates from
 *                  (X,Y) to (I,J) space
 *
 *  Revisions:      none.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "guquit.h"
#include "gridnodes.h"
#include "gridmap.h"
#include "gucommon.h"

#define BUFSIZE 10240

static int reverse = 0;
static int force = 0;
static NODETYPE nt = NT_DD;
static int gridmaptype = GRIDMAP_TYPE_DEF;

typedef int (*mapfn) (void*, double, double, double*, double*);

/**
 */
static void version()
{
    printf("  xy2ij/libgu version %s\n", gu_version);
    exit(0);
}

/**
 */
static void usage()
{
    printf("  Usage: xy2ij [-i {DD|CO}] [-f] [-k] [-r] [-v] -g <grid file> -o <point file>\n");
    printf("  Run \"xy2ij -h\" for more information.\n");

    exit(0);
}

/**
 */
static void info()
{
    printf("  Usage: xy2ij [options] -g <grid file> -o <point file>\n");
    printf("  Where:\n");
    printf("    <grid file> -- text file with node coordinates (see remarks below)\n");
    printf("      (use \"stdin\" or \"-\" for standard input)\n");
    printf("    <point file> -- text file with coordinates to be converted (first two\n");
    printf("      columns used as point coordinates) (use \"stdin\" or \"-\" for standard input)\n");
    printf("  Options:\n");
    printf("    -f -- do not exit with error for points outside grid\n");
    printf("    -i <node type> -- input node type\n");
    printf("    -k -- use kd-tree for mapping\n");
    printf("    -r -- make convertion from index to physical space\n");
    printf("    -v -- verbose / version\n");
    printf("  Node types:\n");
    printf("    DD -- double density nodes (default) \n");
    printf("    CO -- cell corner nodes\n");
    printf("  Description:\n");
    printf("    `xy2ij' reads grid nodes from a file. After that, it reads points from\n");
    printf("     standard input, converts them from (X,Y) to (I,J) space or vice versa,\n");
    printf("     and writes results to the standard output.\n");
    printf("  Remarks:\n");
    printf("    1. The input file must contain header describing the node array dimension:\n");
    printf("       ## <nx> x <ny>\n");
    printf("       where for double density nodes nx = nce1 * 2 + 1, ny = nce2 * 2 + 1;\n");
    printf("       for corner nodes  nx = nce1 + 1, ny = nce2 + 1; and for center nodes\n");
    printf("       nx = nce1, ny = nce2.\n");
    printf("    2. After the header, the grid file must contain (nx * ny) lines with X and\n");
    printf("       Y node coordinates.\n");
    printf("    3. An empty or commented line in the input grid file as well as NaNs for\n");
    printf("       node coordinates indicate an invalid node.\n");
    printf("    4. A grid cell is valid if all corner nodes are valid (not NaNs). Only\n");
    printf("       points in valid cells may be converted between physical and index\n");
    printf("       space.\n");
    printf("    5. The grid (union of all valid grid cells) must be simply connected both in\n");
    printf("       physical and index space.\n");
    printf("       NOTE: from v.1.01.0 this condition can sometimes be bypassed by using\n");
    printf("       \"-k\" option.\n");
    printf("    6. Although `xy2ij' can not use grids of CE (cell centres) type, such grids\n");
    printf("         can be converted to CO or DD type with `getnodes'.\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** gfname, char** ofname)
{
    int i;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-')
            usage();
        else {
            switch (argv[i][1]) {
            case 'i':
                i++;
                if (i == argc)
                    gu_quit("no node type found after \"-i\"");
                if (strcasecmp("dd", argv[i]) == 0)
                    nt = NT_DD;
                else if (strcasecmp("ce", argv[i]) == 0)
                    gu_quit("cell centre node type is not supported by xy2ij");
                else if (strcasecmp("co", argv[i]) == 0)
                    nt = NT_COR;
                else
                    gu_quit("input node type \"%s\" not recognised", argv[i]);
                i++;
                break;
            case 'f':
                i++;
                force = 1;
                break;
            case 'g':
                i++;
                *gfname = argv[i];
                i++;
                break;
            case 'h':
                info();
                break;
            case 'k':
                gridmaptype = GRIDMAP_TYPE_KDTREE;
                i++;
                break;
            case 'o':
                i++;
                *ofname = argv[i];
                i++;
                break;
            case 'r':
                i++;
                reverse = 1;
                break;
            case 'v':
                i++;
                gu_verbose = 1;
                break;
            default:
                usage();
                break;
            }
        }
    }

    if (gu_verbose && argc == 2)
        version();

    if (*gfname == NULL || *ofname == NULL)
        usage();
}

/**
 */
int main(int argc, char* argv[])
{
    char* gfname = NULL;
    char* ofname = NULL;
    FILE* of = NULL;
    gridnodes* gn = NULL;
    void* map = NULL;
    mapfn fn = NULL;
    char buf[BUFSIZE];
    int count, count_success;

    parse_commandline(argc, argv, &gfname, &ofname);
    if (nt == NT_DD) {
        gridnodes* gndd = gridnodes_read(gfname, NT_DD);

        gridnodes_validate(gndd);
        gn = gridnodes_transform(gndd, NT_COR);
        gridnodes_destroy(gndd);
    } else {
        gn = gridnodes_read(gfname, NT_COR);
        gridnodes_validate(gn);
    }

    /*
     * build grid map 
     */
    if (gu_verbose)
        fprintf(stderr, "## parsing the grid into %s...", (gridmaptype == GRIDMAP_TYPE_KDTREE) ? "kd-tree" : "binary tree");
    map = gridmap_build(gridnodes_getnce1(gn), gridnodes_getnce2(gn), gridnodes_getx(gn), gridnodes_gety(gn), gridmaptype);
    if (gu_verbose)
        fprintf(stderr, "done\n");

    if (strcmp(ofname, "stdin") == 0 || strcmp(ofname, "-") == 0)
        of = stdin;
    else
        of = gu_fopen(ofname, "r");

    fn = (reverse) ? (mapfn) gridmap_fij2xy : (mapfn) gridmap_xy2fij;

    /*
     * read points to be mapped, do the mapping and write results to stdout 
     */
    if (gu_verbose)
        fprintf(stderr, "## mapping the points: ");
    count = 0;
    count_success = 0;
    while (fgets(buf, BUFSIZE, of) != NULL) {
        char rem[BUFSIZE] = "";
        double xc, yc, ic, jc;

        if (sscanf(buf, "%lf %lf %[^\n]", &xc, &yc, rem) >= 2) {
            if (fn(map, xc, yc, &ic, &jc)) {
                if (!isnan(ic)) {
                    count_success++;
                    printf("%.15g %.15g %s\n", ic, jc, rem);
                } else
                    printf("NaN NaN %s\n", rem);
            } else {
                if (!force)
                    gu_quit("could not convert (%.15g, %.15g) from %s to %s space", xc, yc, (reverse) ? "index" : "physical", (reverse) ? "physical" : "index");
                else
                    printf("NaN NaN %s\n", rem);
            }
            count++;
            if (gu_verbose && count % 1000 == 0)
                fprintf(stderr, ".");
        } else
            printf("%s", buf);
    }
    if (gu_verbose) {
        fprintf(stderr, "\n");
        fprintf(stderr, "## total mappings: %d\n", count);
        fprintf(stderr, "##   successful: %d\n", count_success);
        fprintf(stderr, "##   unsuccessful: %d\n", count - count_success);
    }

    if (of != stdin)
        fclose(of);
    gridmap_destroy(map);
    gridnodes_destroy(gn);

    return 0;
}
