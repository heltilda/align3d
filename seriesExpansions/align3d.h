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
#include "gmp.h"

enum allPropsArg { doFillZ, doBridgeZ, doBridgeS, doPropS, doBridge_dCdZ, doGetGradW };
enum allElsArg { doClearP, doRenormZ, doCalcP, doDivZ, do_dCdZ, doSetS1x, doAddS1x, doSetSxN, doAddSxN };

extern int call_GetNeighbors(int, char **);
extern void getNeighbors(ccInt, ccFloat);
extern void forEachNeighbor(ccInt, ccInt, ccFloat, ccInt, void(*)(ccInt));

extern int call_IterateProbs(int, char **);
extern ccBool setMasks();
extern ccBool getOptState(ccFloat *);
extern double getC(const gsl_vector *, void *);
extern void getGradC(const gsl_vector *, void *, gsl_vector *);
extern void getCAndGradC(const gsl_vector *, void *, double *, gsl_vector *);
extern void load_fw(const gsl_vector *);
extern void save_grad_fw(const gsl_vector *);
extern ccFloat IterateProbs(ccBool);
extern void forEachPropagator(linkedlist *, ccFloat *, mpf_t **, void(*)(ccInt), int, allPropsArg);
extern void forEachElement(allElsArg);
extern void RenormZ(ccFloat *, ccInt, ccFloat *, ccFloat);
extern void Z_prop(ccInt);
extern void Z_bridge(ccInt);
extern void CountNeighbors(ccInt);
extern void GetAllChains(ccInt, ccInt, ccFloat, ccBool, ccBool);

extern int call_GaussianChain(int, char **);
extern ccFloat GaussProb(ccFloat, ccInt, ccInt, char);
extern ccFloat sqSpotDistance(ccInt, ccInt, ccFloat, ccFloat, ccFloat);

extern int call_Entropy(int, char **);
extern ccBool mappingIsAllowed(ccInt, ccInt, ccInt, ccBool);
extern ccBool hasSpots(ccInt);

extern int call_clock(int, char **);
