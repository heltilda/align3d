/*
 *  Interpolation.h
 *  Yazoo
 *
 *  Created by Brian Ross on 1/4/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef NumericsHeader
#define NumericsHeader

#include "userfn.h"

extern const double pi;


typedef struct {
    ccInt dimensions;       // the dimensionality of the table
    ccInt TablesNum;        // how many consecutive tables are stored
    ccInt EntriesNum;       // the number of elements in each table; equals the product of RowSizes
    double *StartVal;       // the list of starting values for the rows of each table (one number for each table)
    double *StepVal;        // the list of discretization steps for the rows of each table (one number for each table)
    ccInt *RowSize;         // the number of points in each dimension
    double *DataPtr;        // the EntriesNum*TablesNum values stored in the tables
} InterpolationTables;


extern ccBool CheckArgInfo(arg_info *, const int *, const int *, int, int, const char *);
extern int call_SampleDistribution(int argc, char **argv);
extern void SampleDistribution(InterpolationTables *, ccInt *, double **, ccInt, ccInt, int *, double *, ccBool);
extern int call_Interpolate(int, char **);
extern double Interpolate(InterpolationTables *, ccInt, double *, int *, signed char **);
extern double Interpolate_1_dim(double *, double *, ccInt *, double *, double *, ccInt, int *, signed char **);
extern double InterpolateInverse(InterpolationTables *, ccInt, double *, ccInt);
extern double InterpolateInverse_1_dim(InterpolationTables *, ccInt, ccInt, double *, ccInt);
extern int loadInterpolationTable(char **, arg_info *, InterpolationTables *);

#endif
