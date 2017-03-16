/*
 *  align3d.h
 *  align3d
 *
 *  Created by Brian Ross on 2/15/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <gsl/gsl_multimin.h>
#include "cicada.h"
#include "lnklst.h"

extern int call_GetNeighbors(int, char **);
extern ccInt getNeighbors(ccInt, double, linkedlist *);
extern void forEachNeighbor(ccInt, double, ccInt, void(*)(ccInt));

extern int call_IterateProbs(int, char **);
extern int GetC(double *, ccInt, double, int);
extern double GetZ(const gsl_vector *, void *);
extern void GetGradZ(const gsl_vector *, void *, gsl_vector *);
extern void GetZAndGradZ(const gsl_vector *, void *, double *, gsl_vector *);
extern void load_f(const gsl_vector *);
extern void save_grad_f(const gsl_vector *);
extern double IterateProbs(int);
extern void PropArray(linkedlist *, double *, void(*)(ccInt), int, int);
extern void FillArray(int);
extern void RenormZ(double *, ccInt, double *, double);
extern void Z_prop(ccInt);
extern void Z_bridge(ccInt);
extern void CountNeighbors(ccInt);
extern void InitW();
extern void GetAllChains(ccInt, ccInt, double, ccBool);

extern int call_GaussianChain(int, char **);
extern double GaussProb(double, ccInt, ccInt, char);
extern double EffectiveL2(double);
extern double SqDotDistance(ccInt, ccInt, double, double, double);
extern double AddLog(double, double, double);
extern double max(double, double);
extern int call_InitGaussRand(int, char **);
extern int call_GaussRand(int, char **);

extern int call_Entropy(int, char **);
extern ccBool AreSpots(ccInt);
