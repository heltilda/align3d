/*
 *  userfn.c(pp) - for the user's own C/C++ functions
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

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "intrpt.h"
#include "ccmain.h"
#include "userfn.h"
#include "MonteCarlo.h"
#include "Interpolation.h"
#include "LinMath.h"
#include "align3d.h"


// **********************************************
// #include any user-defined header files here








// **********************************************
// UserFunctions defines the { names, function addresses } of user's C routines.
// Each routine must be of the form:  int RoutineName(int argc, char **argv)


userFunction UserFunctions[] = { { "cicada", &runCicada },
                
                    // Wormulator C/C++ routines
                
                { "SampleDistribution", &call_SampleDistribution }, { "Interpolate", &call_Interpolate },
                { "Measure_R_2N", &call_Measure_R_2N }, { "Measure_R_dot_u0", &call_Measure_R_dot_u0 },
                { "PropagateChain", &call_PropagateChain }, { "MakeChains", &call_MakeChains },
                
                    // align3d C/C++ routines
                
                { "GetNeighbors", &call_GetNeighbors },
                { "IterateProbs", &call_IterateProbs },
                { "GaussianChain", &call_GaussianChain }, { "InitGaussRand", &call_InitGaussRand }, { "GaussRand", &call_GaussRand },
                { "Entropy", &call_Entropy },
                
                    // generally useful C/C++ routines
                
                { "StartTimer", &call_StartTimer }, { "StopTimer", &call_StopTimer },
                { "SetRND", &call_SetRND }, { "GetRND", &call_GetRND }          };


const ccInt userFunctionsNum = (ccInt) (sizeof(UserFunctions)/sizeof(userFunction));      // for Cicada's own records





// **********************************************
// C routines may be inserted here, or included in separate source files.




// Example of a user-written C routine.
// Try it out with:  call("pass2nums", 5, pi)

int pass2nums(int argc, char **argv)
{
    ccInt *a;
    ccFloat b;
    
    getArgs(argc, argv, &a, byValue(&b));
    printf("passed %i by reference and %g by value\n", *a, b);
    
    return passed;
}





// Misc -- useful for loading arguments into the user's C routines.
// 
// example:
//
// int myFunction(int argc, char **argv)
// {
//     int *var1;
//     ccBool var2, *var3;
//     
//     getArgs(argc, argv, &var1, byValue(&var2), &var3)
//     
//     ...     // assumes call() passed (int, bool, bool)
//     
// }

void getArgs(int argc, char **argv, ...)
{
    int loopArg;
    va_list theArgs;
    char **nextarg;
    arg_info *myArgInfo = (arg_info *) argv[argc];
    
    va_start(theArgs, argv);
    for (loopArg = 0; loopArg < argc; loopArg++)  {
        nextarg = va_arg(theArgs, char **);
        if (nextarg != NULL)
            *nextarg = argv[loopArg];
        else  {
            nextarg = va_arg(theArgs, char **);
            if (nextarg == NULL)  {
                loopArg = va_arg(theArgs, int) - 1;
                if (loopArg < -1)  break;       }
            else  {
                size_t numBytes = (size_t) myArgInfo[loopArg].argIndices*typeSizes[myArgInfo[loopArg].argType];
                if (numBytes > 0)  memcpy((void *) nextarg, (void *) argv[loopArg], numBytes);
    }   }   }
    va_end(theArgs);
    
    return;
}
