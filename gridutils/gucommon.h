/******************************************************************************
 *
 *  File:           gucommon.h
 *  
 *  Created         22/01/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Header file with some common stuff for grid utilities
 *  Revisions:      none
 *
 *****************************************************************************/

#if !defined(_GUCOMMON_H)
#define _GUCOMMON_H

#include "guquit.h"

extern int gu_verbose;          /* set verbosity from your application */
extern char* gu_version;

FILE* gu_fopen(const char* path, const char* mode);
void* gu_alloc2d(size_t nj, size_t ni, size_t unitsize);
void gu_free2d(void* dummy);
int** gu_readmask(char* fname, int nx, int ny);

#endif
