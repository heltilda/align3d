/*
 *  Interpolation.c
 *  Yazoo
 *
 *  Created by Brian Ross on 1/4/06.
 *  Copyright 2006 Brian Ross. All rights reserved.
 *
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "intrpt.h"
#include "Interpolation.h"


const double pi = 3.141592653589793238462643383;



// Checks the arguments of functions run by call()

ccBool CheckArgInfo(arg_info *ArgInfo, const int *targetTypes, const int *targetIndices, int argc, int targetArgc, const char *functionName)
{
    ccInt c1, prevIndices;
    int prevNegNum = 0;
    
	if (argc != targetArgc)  {
        printf(" *** call(\"%s\") error:  wrong number of args (%i vs expected %i) ***\n", functionName, argc, targetArgc);
        return 1;       }
    
    for (c1 = 0; c1 < argc; c1++)  {
        if ( (targetTypes[c1] >= 0) && (ArgInfo[c1].argType != targetTypes[c1]) )  {
            printf(" *** call(\"%s\") error (wrong arg type:  %i vs expected %i) on argument %i ***\n", functionName, ArgInfo[c1].argType, targetTypes[c1], c1);
            return 1;       }
        if ( targetIndices[c1] >= 0 )  {
            if ( ArgInfo[c1].argIndices != targetIndices[c1] )  {
                printf(" *** call(\"%s\") error (wrong arg # of indices:  %i vs expected %i) on argument %i ***\n", functionName, ArgInfo[c1].argIndices, targetIndices[c1], c1);
                return 2;
        }   }
        else    {
            if (prevNegNum == targetIndices[c1])  {
                if (ArgInfo[c1].argIndices != prevIndices)  {
                    printf(" *** call(\"%s\") error (wrong arg # of indices:  %i vs expected %i) on argument %i ***\n", functionName, ArgInfo[c1].argIndices, prevIndices, c1);
                    return 3;
            }   }
            else  {
                prevNegNum = targetIndices[c1];
                prevIndices = ArgInfo[c1].argIndices;
        }   }
    }

    return passed;
}



// Next 2 routines:  draw random samples from a discretized distribution.
// Uses both a random()->sample mapping (CDF_Inverse) complemented by the CDF for cases in which CDF_I has a steep slope, hence interpolation is unreliable.
// For the latter case, obtain the sample from the CDF by the same bracketing procedure as in MonteCarlo.zoo for generating CDF_Inverse.

int call_SampleDistribution(int argc, char **argv)
{
    double *sampleVector, max_ID, *workspace, **vector_ptrs;
    int *derivatives;
    ccInt numSamples, sampleVectorSize, loopCDF, c1, *tableIDs;
    InterpolationTables *invCDF;
    ccBool ifInverse;
    arg_info *ArgInfo;
    
    ArgInfo = (arg_info *) *(argv+argc);
    
    if ( (argc < 3) || ((ArgInfo+argc-4)->argType != int_type) || ((ArgInfo+argc-3)->argType != double_type)
                || ((ArgInfo+argc-2)->argType != int_type) )      {
        printf("usage: call(\"SampleDistribution\", (InterpolationTables) CDF_Inverse.export, (ccInt) tableIDs[], (double) sample vector[], ");
        printf("(bool) if_inverse_tables\n");
        return 1;       }
    
    sampleVectorSize = (argc-4)/4;
    
    tableIDs = (ccInt *) *(argv+argc-4);
    sampleVector = (double *) *(argv+argc-3);
    ifInverse = *(int *) *(argv+argc-2);
    
    numSamples = (ArgInfo+argc-4)->argIndices;
    max_ID = 0;
    for (c1 = 0; c1 < numSamples; c1++)    {
    if (*(tableIDs+c1) > max_ID)     {
        max_ID = *(tableIDs+c1);
    }}
    
    if (4*sampleVectorSize != argc - 4)  {
        printf("SampleDistribution() error:  some problems with the interpolation tables that were passed\n");
        return 1;       }
    
    if (sampleVectorSize * numSamples != (ArgInfo+argc-4)->argIndices)  {
        printf("SampleDistribution() error:  num_sample_vector_elements must be dimensions * num_TableID_elements\n");
        return 1;       }
    
    invCDF = (InterpolationTables *) malloc(sampleVectorSize*sizeof(InterpolationTables));
    workspace = (double *) malloc((size_t) sampleVectorSize*sizeof(double));
    vector_ptrs = (double **) malloc((size_t) sampleVectorSize*sizeof(double *));
    derivatives = (int *) calloc((size_t) sampleVectorSize, (size_t) 1);
    
    if ((invCDF == NULL) || (workspace == NULL) || (vector_ptrs == NULL) || (derivatives == NULL))       {
        printf("SampleDistribution() error:  out of memory\n");
        if (invCDF != NULL)  free((void *) invCDF);
        if (workspace != NULL)  free((void *) workspace);
        if (vector_ptrs != NULL)  free((void *) vector_ptrs);
        if (derivatives != NULL)  free((void *) derivatives);
        return 1;       }
    
    for (loopCDF = 0; loopCDF < sampleVectorSize; loopCDF++)     {
        *(vector_ptrs + loopCDF) = sampleVector + loopCDF * numSamples;
        if (loadInterpolationTable(argv + 4*loopCDF, ArgInfo + 4*loopCDF, invCDF + loopCDF) != passed)  break;
        if ((invCDF + loopCDF)->TablesNum <= max_ID)      {
            printf("SampleDistribution() error:  too few interpolation tables passed (InterpolationTables not initialized?)\n");
    }   }
    
    if (loopCDF == sampleVectorSize)     // i.e. no break
        SampleDistribution(invCDF, tableIDs, vector_ptrs, sampleVectorSize, numSamples, derivatives, workspace, ifInverse);
    
    free((void *) invCDF);
    free((void *) workspace);
    free((void *) vector_ptrs);
    free((void *) derivatives);
    
    return 0;
}


// Samples an N-dimensional distribution given a series of 1-N-dimensional CDF-inverse interpolation tables.

void SampleDistribution(InterpolationTables *invCDF, ccInt *tableIDs, double **samples, ccInt sampleVectorSize, ccInt numSamples, int *derivatives,
                        double *workspace, ccBool ifInverse)
{
    ccInt loopSample, loopVector;
    
    for (loopSample = 0; loopSample < numSamples; loopSample++)       {
        
        for (loopVector = 0; loopVector < sampleVectorSize; loopVector++)
            workspace[loopVector] = ((double) rand())/RAND_MAX;
        
        if (ifInverse)  {
            
            for (loopVector = 0; loopVector < sampleVectorSize; loopVector++)
                *(workspace + loopVector) = Interpolate(invCDF+loopVector, *(tableIDs+loopSample), workspace, derivatives, NULL);
            
            for (loopVector = 0; loopVector < sampleVectorSize; loopVector++)
                samples[loopVector][loopSample] = workspace[loopVector];       }
        
        else  {
        for (loopVector = 0; loopVector < sampleVectorSize; loopVector++)  {
            samples[loopVector][loopSample] = InterpolateInverse(invCDF, *(tableIDs+loopSample), workspace, loopVector+1);
    }   }}
}



// Next 3 functions:  evaluate a function using a linear interpolation table.

int call_Interpolate(int argc, char **argv)
{
    InterpolationTables tempIT;
    ccInt table_num, dimensions, cD;
    int ifInverse, *derivatives;
    signed char *bounds_flag, **bounds_flag_pass;
    double *argument, *answer;
    arg_info *ArgInfo = (arg_info *) *(argv+argc);
    ccBool no_derivs;

    const int ArgTypes[] = { double_type, double_type, int_type, double_type, int_type, int_type, int_type, double_type,    
                                                            int_type, int_type, double_type, int_type };
    const int ArgIndices[] = { -1, -1, -2, -3, 1, 1, 1, -4, -4, -4, 1, 1 };
    
    if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "Interpolate") != passed)  return 1;
    
    if (loadInterpolationTable(argv, ArgInfo, &tempIT) != passed)  return 1;
    dimensions = *(ccInt *) *(argv+4);
    table_num = *(ccInt *) *(argv+5);
    ifInverse = *(int *) *(argv+6);
    argument = (double *) *(argv+7);
    derivatives = (int *) *(argv+8);
    bounds_flag = (signed char *) *(argv+9);
    answer = (double *) *(argv+10);
    
    if (table_num >= tempIT.TablesNum)  {  printf("Interpolate_ID() error:  table ID is greater than the number of tables\n");  return 1;  }
    
    no_derivs = ccTrue;
    for (cD = 0; cD < dimensions; cD++)  {
        if (*(derivatives+cD) > 0)  no_derivs = ccFalse;
        if (*(derivatives+cD) > 2)  {
            printf("Interpolate() error:  derivative must be 0 (no derivative), 1 or 2 (first/second derivative)\n");
            return 2;
    }   }
    
    bounds_flag_pass = (signed char **) malloc(dimensions*sizeof(signed char *));
    for (cD = 0; cD < dimensions; cD++)  *(bounds_flag_pass+cD) = bounds_flag+cD;
    
    if (tempIT.EntriesNum == 0)  *answer = 0;
    else if (!ifInverse)  {
        *answer = Interpolate(&tempIT, table_num, argument, derivatives, bounds_flag_pass);       }
    
    else if (no_derivs)        {
        
        double low, middle, high, low_eval, middle_eval, high_eval, Target, *top_dim;
        const double Tolerance = 1e-8;
        
        top_dim = argument+dimensions-1;
        Target = *top_dim;
        
        *top_dim = low = *(tempIT.StartVal + dimensions - 1);
        low_eval = Interpolate(&tempIT, table_num, argument, derivatives, bounds_flag_pass);
        *top_dim = high = low + (*(tempIT.StepVal + dimensions - 1)) * (*(tempIT.RowSize + dimensions - 1) - 1);
        high_eval = Interpolate(&tempIT, table_num, argument, derivatives, bounds_flag_pass);
        
        while (high-low >= Tolerance)      {
            *top_dim = middle = (low + high)/2;
            middle_eval = Interpolate(&tempIT, table_num, argument, derivatives, bounds_flag_pass);

            if (middle_eval < Target)       {
                if (low == middle)  break;
                low = middle;
                low_eval = middle_eval;     }
            else        {
                if (high == middle)  break;
                high = middle;
                high_eval = middle_eval;    }
        }
        
        if (low_eval == high_eval)  *answer = low;
        else  *answer = low + (high-low)*(Target-low_eval)/(high_eval-low_eval);
    }
    
    else  {  printf("Interpolate() error:  all derivatives must be 0 when inverting a CDF");  return 2;  }
    
    free(bounds_flag_pass);

    return passed;
}



// Calculates f(x), f'(x), or f''(x) using an interpolation table.  The order of the derivative is the last argument.
// If drawing from CDF, then CDF^-1 is the first argument; CDF (2nd arg) is used instead if CDF^-1 is likely to be inaccurate (large slope).

double Interpolate(InterpolationTables *theTable, ccInt TableID, double *x, int *derivatives, signed char **bounds_flag)
{
    signed char *null_bf = NULL;
    if (bounds_flag == NULL)  bounds_flag = &null_bf;
    
    return Interpolate_1_dim(theTable->StartVal + TableID*theTable->dimensions, theTable->StepVal + TableID*theTable->dimensions, theTable->RowSize,
                    theTable->DataPtr + TableID*theTable->EntriesNum, x, theTable->EntriesNum / (*(theTable->RowSize)), derivatives, bounds_flag);
}


// Interpolates over one dimension of the input vector

double Interpolate_1_dim(double *StartVal, double *StepVal, ccInt *RowSize, double *TableEntries, double *x,
                            ccInt EntriesPerIndex, int *deriv, signed char **bounds_flag)
{
    int idx, ci, total_evals;
    double real_index, remainder, val_1, val_2, f_evals[4];
    signed char **next_bf, *null_bf = NULL;
    
    total_evals = *deriv + 2;
    
    if (*RowSize <= *deriv)  return 0;
    else if (*RowSize == *deriv + 1)  total_evals--;

    real_index = ((*x) - (*StartVal)) / (*StepVal) - 0.5*(*deriv);
    idx = floor(real_index);
    
    if (real_index < 0.)     {
        real_index = 0.;
        idx = 0;
        if (*bounds_flag != NULL)  {
            if (*x < *StartVal)  **bounds_flag = -1;
            else if (*RowSize > 3)  **bounds_flag = -2;
    }   }
    else if (idx + total_evals > *RowSize)  {
        idx = *RowSize - 2 - (*deriv);
        real_index = idx + 1;
        if (*bounds_flag != NULL)  **bounds_flag = 1;       }

    if (*RowSize == *deriv + 1)     {
        real_index = 0.;
        idx = 0;        }

    remainder = real_index - idx;
    
    if (*bounds_flag == NULL)  next_bf = &null_bf;
    else  next_bf = bounds_flag+1;
    
    for (ci = 0; ci < total_evals; ci++)  {
        if (EntriesPerIndex == 1)   f_evals[ci] = *(TableEntries + idx + ci);
        else  f_evals[ci] = Interpolate_1_dim(StartVal+1, StepVal+1, RowSize+1, TableEntries+EntriesPerIndex*((ccInt) idx + ci),
                                    x+1, EntriesPerIndex / (*(RowSize+1)), deriv+1, next_bf);      }
    
    if (*deriv == 0)     {
        val_1 = f_evals[0];
        val_2 = f_evals[1];         }
    
    else if (*deriv == 1)  {
        val_1 = ( f_evals[1] - f_evals[0] ) / (*StepVal);
        val_2 = ( f_evals[2] - f_evals[1] ) / (*StepVal);            }
    
    else if (*deriv == 2)        {
        val_1 = ( f_evals[0] + f_evals[2] - 2 * f_evals[1] ) / ((*StepVal) * (*StepVal));
        val_2 = ( f_evals[1] + f_evals[3] - 2 * f_evals[2] ) / ((*StepVal) * (*StepVal));      }
    
    if (remainder == 0)  return val_1;                  // if *RowSize == *deriv + 1, then val2 may be NaN or inf
    else  return val_1*(1. - remainder) + val_2*remainder;
}



// Calculates f^-1(x) using numerical interpolation tables.

double InterpolateInverse(InterpolationTables *theTable, ccInt TableID, double *x, ccInt dimensions)
{
    return InterpolateInverse_1_dim(theTable, TableID, 0, x, dimensions);
}


// Interpolates over one dimension of the input vector

double InterpolateInverse_1_dim(InterpolationTables *theTable, ccInt TableID, ccInt OldEntryOffset, double *x, ccInt remaining_dims)
{
    ccInt StartOffset = (TableID+1)*theTable->dimensions - 1;
    ccInt OneRowSize = *(theTable->RowSize+theTable->dimensions-1);
    ccInt bracket[2] = { 0, OneRowSize - 1 }, middle, EntryOffset = OldEntryOffset*OneRowSize;
    double *TableEntries = theTable->DataPtr + TableID*theTable->EntriesNum + EntryOffset;
    double remainder, f_evals[2];
    int ci;
    
    while (bracket[1] - bracket[0] > 1)  {
        middle = (int) floor((bracket[0]+bracket[1])/2);
        if (*(TableEntries + middle) < *x)  bracket[0] = middle;
        else  bracket[1] = middle;        }
    
    if (bracket[0] == bracket[1])  remainder = 0;
    else  remainder = ((*x) - (*(TableEntries + bracket[0]))) / ((*(TableEntries + bracket[1])) - (*(TableEntries + bracket[0])));
    
    for (ci = 0; ci < 2; ci++)  {
        if (remaining_dims == 1)   f_evals[ci] = (*(theTable->StartVal + StartOffset)) + bracket[ci]*(*(theTable->StepVal + StartOffset));
        else  f_evals[ci] = InterpolateInverse_1_dim(theTable+1, TableID, EntryOffset + bracket[ci], x+1, remaining_dims-1);      }
    
    return f_evals[0]*(1. - remainder) + f_evals[1]*remainder;
}



// Loads an interpolation table from a call() routine:  4 consecutive arguments

int loadInterpolationTable(char **args_start, arg_info *arg_info_start, InterpolationTables *theTable)
{
    ccInt c1, check_entries;
    
    theTable->dimensions = (arg_info_start+2)->argIndices;
    if (theTable->dimensions == 0)  theTable->TablesNum = 0;
    else  theTable->TablesNum = arg_info_start->argIndices / theTable->dimensions;
    if (theTable->TablesNum == 0)  theTable->EntriesNum = 0;
    else  theTable->EntriesNum = ((arg_info_start+3)->argIndices) / theTable->TablesNum;
    if (args_start != NULL)  {
        theTable->StartVal = (double *) *(args_start+0);
        theTable->StepVal = (double *) *(args_start+1);
        theTable->RowSize = (ccInt *) *(args_start+2);
        theTable->DataPtr = (double *) *(args_start+3);     }
    
    if ((theTable->TablesNum) * (theTable->EntriesNum) != (arg_info_start+3)->argIndices)      {
        printf("Error:  each interpolation table needs to be of the same length\n");
        return 1;           }

    if ((theTable->dimensions) * (theTable->TablesNum) != (arg_info_start)->argIndices)      {
        printf("Error:  each start/step vector needs to be of length dimensions*tables_num\n");
        return 1;           }
    
    check_entries = 1;
    for (c1 = 0; c1 < theTable->dimensions; c1++)  check_entries *= *(theTable->RowSize + c1);
    if (theTable->dimensions == 0)  check_entries = 0;

    if (check_entries != theTable->EntriesNum)      {
        printf("Error:  the table is the wrong size for the given RowSize vector\n");
        return 1;           }    
    
    return passed;
}
