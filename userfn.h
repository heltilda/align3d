/*
 *  userfn.h
 *
 *  Cicada
 *  Copyright (C) 2006 Brian C. Ross
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or any later version.

 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef UserFunctions_h
#define UserFunctions_h

#include "lnklst.h"
#include <stdlib.h>



// Comment out the following line if Cicada is to be embedded in a larger C/C++ application
// (i.e. if Cicada's own main() function is not needed).

#define CicadaMainProgram



// The required type of user-defined C or C++ functions

typedef struct {
    const char *functionName;
    int(*functionPtr)(int, char **);
} userFunction;


// arg_info:  a data type for the final argv parameter to each user function call.
// This gives a report on the sizes, etc. of each argument.

typedef struct {
    ccInt argType;
    ccInt argIndices;
} arg_info;



// Prototypes & misc definitions

#define byValue(a) NULL,a
#define fromArg(n) NULL,NULL,n
#define endArgs NULL,NULL,-1


extern userFunction UserFunctions[];
extern const ccInt userFunctionsNum;

extern int pass2nums(int, char **);

extern void getArgs(int, char **, ...);

#endif

