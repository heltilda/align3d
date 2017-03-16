/*
 *  MonteCarlo.h
 *  Yazoo
 *
 *  Created by Brian Ross on 1/4/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MonteCarloHeader
#define MonteCarloHeader

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include "LinMath.h"
#include "Interpolation.h"


extern gsl_rng *gsl_rand_gen;

extern int call_MakeChains(int, char **);
extern int call_PropagateChain(int, char **);
extern void PropagateChain(vector_3 *, vector_3 *, vector_3 *, vector_3 *, double **, ccInt, ccBool, ccBool);
extern void Displace(vector_3 *, double **, ccInt, vector_3 *, vector_3 *, vector_3 *);
extern void Renorm(vector_3 *);
extern void RenormAxes(vector_3 *, vector_3 *, vector_3 *);
extern void v3Rotate(vector_3 *, vector_3 *, double, vector_3 *);
extern double MakeChain(ccInt, void **);
extern void GetStiffness(double *, gsl_vector *, gsl_matrix *);
extern void GetConstraintDeriv(gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *, vector_3 *, vector_3 *, vector_3 *, vector_3 *, char);
//extern void StoreCostDeriv(gsl_matrix *, double, vector_3 *, size_t, ccInt, double);
extern void StoreCostDeriv(gsl_matrix *, gsl_matrix *, gsl_vector *, double, vector_3 *, size_t, size_t, ccInt, double, vector_3 **, char);
extern void FixEnergies(double *, gsl_vector *, gsl_matrix *);
extern void FixDerivs(gsl_matrix *, gsl_vector *, gsl_matrix *, gsl_vector *);
extern void test_grad(gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_vector *, gsl_vector *, vector_3 *, vector_3 *, vector_3 *, vector_3 *, char);
extern void sample_grads(gsl_vector *, gsl_matrix *, vector_3 *, vector_3 *, vector_3 *, vector_3 *,
                  double *, ccInt, double, ccBool);
extern double sg_angle(double *, double);
extern double MeasureConcentration(ccInt, void **);
extern int call_Measure_R_2N(int, char **);
extern double Measure_R_2N(ccInt, void **);
extern int call_Measure_R_dot_u0(int, char **);
extern double Measure_R_dot_u0(ccInt, void **);
extern void GetMeanAndError(ccInt *, ccInt, ccInt, ccBool, ccBool,
                    double(*SampleFunction)(ccInt, void **), void **, double *, double *);
extern int call_CountChainCrossings(int, char **);
extern int call_StartTimer(int, char **);
extern int call_StopTimer(int, char **);
extern int call_SetRND(int, char **);
extern int call_GetRND(int, char **);
extern int SetGetRND(int, char **, ccInt);

#endif
