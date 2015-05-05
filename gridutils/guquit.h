/******************************************************************************
 *
 *  File:           guquit.h
 *  
 *  Created         05/05/2015
 *  
 *  Author:         Pavel Sakov
 *                  BoM
 *  
 *  Purpose:        Header file for quit function
 *  Revisions:      none
 *
 *****************************************************************************/

#if !defined(_GUQUIT_H)
#define _GUQUIT_H

typedef void (*gu_quitfn) (char* format, ...);
extern gu_quitfn gu_quit;
void gu_setquitfn(gu_quitfn quitfn);

#endif
