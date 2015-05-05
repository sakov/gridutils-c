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

typedef void (*gu_quitfn) (char* format, ...);

extern int gu_verbose;          /* set verbosity from your application */
extern char* gu_version;
extern gu_quitfn gu_quit;

FILE* gu_fopen(const char* path, const char* mode);
void* gu_alloc2d(int n1, int n2, size_t size);
void gu_free2d(void* dummy);
int** gu_readmask(char* fname, int nx, int ny);
void gu_setquitfn(gu_quitfn quitfn);

#endif
