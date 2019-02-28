/*
 *  align3d.cpp
 *  align3d
 *
 *  Created by Brian Ross on 2/15/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include "gmp.h"

#include "align3d.h"
#include "lnklst.h"
#include "intrpt.h"
#include "userfn.h"


ccInt numLoci, numSpots, numColors, numZterms, numNeighborLengths, loopZterm, propLocus, bridgingLocus, doublePrecision, n_skip_max, alpha_0, global_alpha_base, boundary, numMinParams;
ccInt *ZtermFirstLocus, *fixedLoci, *fixedSpots, *colorFirstSpot, *locusColor, *alphas, *calcs, *numSkippedColors, *neighbors, *neighborIdx, *numLociPerColor, numModelRegimes;
ccFloat max_pfn_mismatch, max_spot_overnorm, exaggeration, Zpropagated, l_step, multiplier;
ccFloat *model, *l, *p_fn, *N_difference, *extraOcc, lastC, propDl, *spotX, *spotY, *spotZ, *spotDx, *spotDy, *spotDz, *f, *w, *grad_fw, *nlog_fw, *fw;
ccFloat *Z_1x_norm, *Z_xN_norm, *sensitivity_norm;
ccFloat *ZtermWeights, *Z, *Z_norm, *dC_dZ, dC_dZ_factor, *locusMappingProbs, *spotMappingProbs;
ccBool *locusMask, *spotMask, *spotWasUsed, *avoidConstrainedSpot, *breakInContour, prop_hadContourBreak, setUnboundPenalty, overBoundary, propMode;
linkedlist *p, *Z_1x, *Z_xN, *ZtoProp, *sensitivity;
linkedlist **dp_arrays[] = { &p, &Z_1x, &Z_xN, &sensitivity };
ccFloat **dp_norms[] = { &Z_norm, &Z_1x_norm, &Z_xN_norm, &sensitivity_norm };
mpf_t **p_mp, **Z_1x_mp, **Z_xN_mp, **sensitivity_mp, *Z_mp, *dC_dZ_mp, *grad_fw_mp, **ZtoProp_mp, dC_dp_mp, dC_dZ_factor_mp, w_factor_mp, multiplier_mp, p_value_mp, tmp_mp, tmp2_mp;
mpf_t ***mp_arrays[] = { &p_mp, &Z_1x_mp, &Z_xN_mp, &sensitivity_mp };
mpf_t *mp_consts[] = { &dC_dp_mp, &dC_dZ_factor_mp, &w_factor_mp, &multiplier_mp, &p_value_mp, &tmp_mp, &tmp2_mp };

gsl_multimin_fdfminimizer *opt_struct;
const gsl_rng_type *gsl_rand_gen_type;

const ccFloat pi = 3.141592653589793238462643383;





// *********** GetNeigbors() ***********


// Finds the neighbors of each spot and stores them in the neighbors array

int call_GetNeighbors(int argc, char **argv)
{
    ccInt c1;
    ccFloat pCutoff;
    arg_info *argInfo = (arg_info *) *(argv+argc);
    
    getArgs(argc, argv, &spotX, &spotY, &spotZ, &spotDx, &spotDy, &spotDz, &colorFirstSpot, &neighbors, &neighborIdx,
                    byValue(&l_step), byValue(&pCutoff), &model);
    
    numSpots = argInfo[0].argIndices;
    numColors = argInfo[6].argIndices - 1;
    numModelRegimes = argInfo[11].argIndices / 3;
    numNeighborLengths = argInfo[8].argIndices / (numSpots*numColors) - 1;
    exaggeration = 1.;
    
    for (c1 = 0; c1 < numSpots; c1++)  getNeighbors(c1, pCutoff);
    
    return passed;
}



// getNeighbors() finds the neighbors of a single given image spot, storing them in neighborLL

void getNeighbors(ccInt fromSpot, ccFloat pCutoff)
{
    ccInt neighborCounter, loopImageSpot, loopLinkerLength, loopColor, dummy_calcs;
    ccBool doAddSpot;
    
    calcs = &dummy_calcs;
    
    neighborCounter = fromSpot*(numSpots-1);
    for (loopColor = 0; loopColor < numColors; loopColor++)  {
        
        ccInt *oneColorIdxList = neighborIdx + (fromSpot*numColors + loopColor)*(numNeighborLengths+1);
        
        for (loopLinkerLength = 0; loopLinkerLength < numNeighborLengths; loopLinkerLength++)  {
            oneColorIdxList[loopLinkerLength] = neighborCounter;
            for (loopImageSpot = colorFirstSpot[loopColor]; loopImageSpot < colorFirstSpot[loopColor+1]; loopImageSpot++)  {
            if (loopImageSpot != fromSpot)  {
                
                if (loopLinkerLength == numNeighborLengths-1)  doAddSpot = ccTrue;
                else  doAddSpot = (GaussProb((loopLinkerLength+1) * l_step, fromSpot, loopImageSpot, 1) >= pCutoff);
                if (loopLinkerLength > 0)  {
                    if (GaussProb(loopLinkerLength * l_step, fromSpot, loopImageSpot, 1) >= pCutoff)  doAddSpot = ccFalse;      }
                
                if (doAddSpot)  {
                    neighbors[neighborCounter] = loopImageSpot;
                    neighborCounter++;
        }   }}  }
        
        oneColorIdxList[numNeighborLengths] = neighborCounter;
    }
}




// forEachNeighbor() performs a given operation (toDo()) for each neighbor of imageSpot, where the neighbor range is determined by 'dl'

void forEachNeighbor(ccInt neighborLocus, ccInt fromSpot, ccFloat dl, ccInt neighborColor, void(*toDo)(ccInt))
{
    ccInt maxNeighborList, cN;
    
    if (overBoundary)  {
    for (cN = colorFirstSpot[neighborColor]; cN < colorFirstSpot[neighborColor + 1]; cN++)  {
        toDo(cN);
    }}
    
    else  {
        ccInt *oneColorIdxList = neighborIdx + (fromSpot*numColors + neighborColor)*(numNeighborLengths+1);
        
        maxNeighborList = floor(dl / l_step);
        if (maxNeighborList >= numNeighborLengths)  maxNeighborList = numNeighborLengths-1;
        for (cN = oneColorIdxList[0]; cN < oneColorIdxList[maxNeighborList+1]; cN++)   {
        if (mappingIsAllowed(neighborLocus, neighbors[cN], -1, ccFalse))  {
            toDo(neighbors[cN]);
    }   }}
}







// *********** IterateZ() ***********


#define mpDo(a,b) if (doublePrecision==0) {a} else {b}


// This is run by call("Iterate", ...).  It runs the overlap-potential optimizer for a given number of iterations,
// or until convergence.

int call_IterateProbs(int argc, char **argv)
{
    ccInt loopColor, i, alpha, c1, *iterations, maxIterations, logsize, optMethod;
    ccFloat *nlog_f, *nlog_w, *C, *LogSource, *LogOutput, *FE, init_step_size, tMax, *tElapsed;
    gsl_multimin_function_fdf opt_fun_grad;
    gsl_vector *start_f;
    int iteration_rtrn;
    clock_t t_start, t_end;
    arg_info *argInfo = (arg_info *) argv[argc];
    ccBool byEnumeration, calcExact, convergenceTest;
    
    const gsl_multimin_fdfminimizer_type *opt_alg_types[] = { gsl_multimin_fdfminimizer_conjugate_fr,
                                    gsl_multimin_fdfminimizer_conjugate_pr, gsl_multimin_fdfminimizer_vector_bfgs2, gsl_multimin_fdfminimizer_steepest_descent };
    
    getArgs(argc, argv, &spotX, &spotY, &spotZ, &spotDx, &spotDy, &spotDz, &spotMappingProbs, &nlog_f, &nlog_w, &grad_fw,
                    &colorFirstSpot, &l, &locusColor, &ZtermFirstLocus, &ZtermWeights, &fixedLoci, &fixedSpots, &breakInContour,
                    &p, &Z_1x, &Z_xN, &Z_1x_norm, &Z_xN_norm, &Z_norm,
                    &dC_dZ, &sensitivity, &sensitivity_norm, &locusMappingProbs, &neighbors, &neighborIdx,
                    byValue(&l_step), &model, byValue(&n_skip_max), &p_fn,
                    byValue(&doublePrecision), byValue(&exaggeration), byValue(&max_pfn_mismatch), byValue(&max_spot_overnorm), byValue(&init_step_size),
                    &C, &LogSource, &LogOutput, &calcs, &tElapsed, &FE, &iterations, byValue(&maxIterations), byValue(&tMax),
                    byValue(&optMethod), byValue(&setUnboundPenalty), byValue(&byEnumeration), byValue(&calcExact), &avoidConstrainedSpot);
    
    numSpots = argInfo[0].argIndices;
    numColors = argInfo[10].argIndices - 1;
    numLoci = argInfo[11].argIndices;
    numZterms = argInfo[13].argIndices - 1;
    numNeighborLengths = argInfo[29].argIndices / (numSpots*numColors) - 1;
    numModelRegimes = argInfo[31].argIndices / 3;
    logsize = argInfo[40].argIndices;
    
    if (numLoci*numSpots == 0)  return passed;
    
    f = (ccFloat *) malloc(numSpots * sizeof(ccFloat));
    w = (ccFloat *) malloc(numColors * sizeof(ccFloat));
    Z = (ccFloat *) malloc(numLoci * sizeof(ccFloat));
    alphas = (ccInt *) malloc(numLoci * sizeof(ccInt));
    locusMask = (ccBool *) malloc(numLoci * sizeof(ccBool));
    spotMask = (ccBool *) malloc(numSpots * sizeof(ccBool));
    spotWasUsed = (ccBool *) malloc(numSpots * sizeof(ccBool));
    N_difference = (ccFloat *) malloc(numColors * sizeof(ccFloat));
    extraOcc = (ccFloat *) malloc(numColors * sizeof(ccFloat));
    numSkippedColors = (ccInt *) malloc(numColors * sizeof(ccInt));
    numLociPerColor = (ccInt *) malloc(numColors * sizeof(ccInt));
    if ((f == NULL) || (w == NULL) || (Z == NULL) || (alphas == NULL) || (locusMask == NULL) || (spotMask == NULL)
            || (spotWasUsed == NULL) || (N_difference == NULL) || (extraOcc == NULL) || (numSkippedColors == NULL) || (numLociPerColor == NULL))  return 1;
    
    for (loopColor = 0; loopColor < numColors; loopColor++)  {  N_difference[loopColor] = 0.;  extraOcc[loopColor] = 0.;  }
    max_pfn_mismatch *= numColors;
    max_spot_overnorm *= numSpots;
    
    if (!setUnboundPenalty)  {  numMinParams = numSpots;  nlog_fw = nlog_f;  fw = f;  }
    else  {  numMinParams = numColors;  nlog_fw = nlog_w;  fw = w;  }
    
    for (c1 = 0; c1 < numSpots; c1++)  f[c1] = exp(-nlog_f[c1]);
    for (c1 = 0; c1 < numColors; c1++)  {  w[c1] = exp(-nlog_w[c1]);  numLociPerColor[c1] = 0;  }
    for (i = 0; i < numLoci; i++)  numLociPerColor[locusColor[i]]++;
    
    if ((doublePrecision != 0) && (!byEnumeration))  {
        
        mpf_set_default_prec(doublePrecision);
        
        for (c1 = 0; c1 < 7; c1++)  mpf_init(*mp_consts[c1]);
        
        Z_mp = (mpf_t *) malloc(numLoci*sizeof(mpf_t));
        dC_dZ_mp = (mpf_t *) malloc(numLoci*sizeof(mpf_t));
        grad_fw_mp = (mpf_t *) malloc(numMinParams*sizeof(mpf_t));
        if ((Z_mp == NULL) || (dC_dZ_mp == NULL) || (grad_fw_mp == NULL))  return 1;
        for (i = 0; i < numLoci; i++)  {  mpf_init(Z_mp[i]);  mpf_init(dC_dZ_mp[i]);  }
        for (alpha = 0; alpha < numMinParams; alpha++)  mpf_init(grad_fw_mp[alpha]);
        
        for (c1 = 0; c1 < 4; c1++)  {
            *mp_arrays[c1] = (mpf_t **) malloc(numLoci*sizeof(mpf_t *));
            if (*mp_arrays[c1] == NULL)  return 1;
            for (i = 0; i < numLoci; i++)  {
                ccInt alpha, one_color = locusColor[i];
                ccInt alpha_base = colorFirstSpot[one_color];
                
                (*mp_arrays[c1])[i] = (mpf_t *) malloc((colorFirstSpot[one_color+1] - alpha_base)*sizeof(mpf_t));
                if ((*mp_arrays[c1])[i] == NULL)  return 1;
                for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)    {
                    mpf_init((*mp_arrays[c1])[i][alpha-alpha_base]);
    }   }   }   }
    
    
    t_start = clock();
    
    if (byEnumeration)  {
        
        ccInt alpha;
        
        
        for (i = 0; i < numLoci; i++)  {
            
            ccInt one_color = locusColor[i];
            ccInt alpha_base = colorFirstSpot[one_color];
            
            for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)    {
                LL_Double(p + i, 1)[alpha - alpha_base] = 0.;
        }   }
        
        
        *Z = 0.;
        for (loopZterm = 0; loopZterm < numZterms; loopZterm++)  {
        if (setMasks())  {
            for (alpha = 0; alpha < numSpots; alpha++)  spotWasUsed[alpha] = ccFalse;
            
            GetAllChains(-1, 0, ZtermWeights[loopZterm], calcExact, ccTrue);
        }}
        
        
        if (*Z != 0.)  {
        for (i = 0; i < numLoci; i++)  {
            
            ccInt one_color = locusColor[i];
            ccInt alpha_base = colorFirstSpot[one_color];
            
            for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)    {
                LL_Double(p + i, 1)[alpha - alpha_base] /= *Z;
        }}  }
        
        *Z_norm = 0;
    }
    
    else if (maxIterations == 0)  {
        
        C[0] = IterateProbs(ccTrue);
        
        for (c1 = 0; c1 < numMinParams; c1++)  {
            mpDo( grad_fw[c1] = -grad_fw[c1] * fw[c1]; , grad_fw[c1] = -mpf_get_d(grad_fw_mp[c1]) * fw[c1]; ) }
        
        C[1] = 0.;
        for (c1 = 0; c1 < numMinParams; c1++)  mpDo( C[1] += (ccFloat) grad_fw[c1] * grad_fw[c1]; , C[1] += (ccFloat) mpf_get_d(grad_fw_mp[c1]) * mpf_get_d(grad_fw_mp[c1]); )
        C[1] = sqrt(C[1]);          }
    
    else if (optMethod == 0)  {
        
        C[2*(*iterations)+0] = IterateProbs(ccFalse);
        C[2*(*iterations)+1] = 0.;              }
    
    else  {
        
        opt_fun_grad.n = numMinParams;
        opt_fun_grad.f = &getC;
        opt_fun_grad.df = &getGradC;
        opt_fun_grad.fdf = &getCAndGradC;
        opt_fun_grad.params = NULL;
        
        start_f = gsl_vector_calloc(numMinParams);
        for (c1 = 0; c1 < numMinParams; c1++)
            gsl_vector_set(start_f, c1, (double) nlog_fw[c1]);
        
        opt_struct = gsl_multimin_fdfminimizer_alloc(opt_alg_types[optMethod - 1], numMinParams);
        if (opt_struct == 0)  return 10;
        gsl_multimin_fdfminimizer_set(opt_struct, &opt_fun_grad, start_f, (double) init_step_size, 0.1);
        
        convergenceTest = getOptState(C + 2*(*iterations));
        
        if (!convergenceTest)  {
        while (*iterations <= maxIterations)     {
            
            iteration_rtrn = gsl_multimin_fdfminimizer_iterate(opt_struct);
            
            convergenceTest = getOptState(C + 2*(*iterations));
            
            if ((iteration_rtrn == GSL_ENOPROG) && (!convergenceTest))  {       // the sigmoid shape of the cost function (flattening at high -log f/w) can confuse the minimizer
            if (lastC < C[2*(*iterations)])  {
                C[2*(*iterations)] = lastC;
                for (c1 = 0; c1 < numMinParams; c1++)  {
                    nlog_fw[c1] = -log(fw[c1]);
                    gsl_vector_set(start_f, c1, (double) nlog_fw[c1]);  }
                gsl_multimin_fdfminimizer_set(opt_struct, &opt_fun_grad, start_f, (double) init_step_size, 0.1);
                iteration_rtrn = GSL_SUCCESS;
            }}
            
            t_end = clock();
            *tElapsed = ((ccFloat) t_end) / CLOCKS_PER_SEC;
            
            if (logsize > 0)  {
            for (c1 = 0; c1 < logsize; c1++)  {
                LogOutput[logsize*(*iterations) + c1] = LogSource[c1];      }}
            
            if ((convergenceTest) || ((tMax > 0.) && (((ccFloat) (t_end - t_start)) / CLOCKS_PER_SEC > tMax)) || (iteration_rtrn == GSL_ENOPROG))  break;
            
            (*iterations)++;
        }}
        
        load_fw(gsl_multimin_fdfminimizer_x(opt_struct));         // make sure 'p' is the the most up-to-date mapping
        IterateProbs(ccTrue);
        
        gsl_multimin_fdfminimizer_free(opt_struct);
        gsl_vector_free(start_f);
    }
    
    
    if ((doublePrecision == 0) || (byEnumeration))  {
    for (i = 0; i < numLoci; i++)  {
        Z_norm[i] += log(fabs((double) Z[i]));
    }}
    
    
    free((void *) alphas);
    free((void *) locusMask);
    free((void *) spotMask);
    free((void *) spotWasUsed);
    free((void *) N_difference);
    free((void *) extraOcc);
    free((void *) numSkippedColors);
    free((void *) numLociPerColor);
    free((void *) f);
    free((void *) w);
    free((void *) Z);
    
    if ((doublePrecision != 0) && (!byEnumeration))  {
        
        for (c1 = 0; c1 < 4; c1++)  {
            for (i = 0; i < numLoci; i++)  {
                ccInt alpha, one_color = locusColor[i];
                ccInt alpha_base = colorFirstSpot[one_color];
                signed long int log2exp;
                ccFloat mantissa;
                
                if (dp_arrays[c1] == &p)  mpf_set(tmp_mp, Z_mp[i]);
                else  {
                    mpf_set_d(tmp_mp, 0.);
                    for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)
                        mpf_add(tmp_mp, tmp_mp, (*mp_arrays[c1])[i][alpha-alpha_base]);
                    mpf_abs(tmp_mp, tmp_mp);
                    
                    for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)  {
                        if (mpf_sgn(tmp_mp) != 0)  mpf_div((*mp_arrays[c1])[i][alpha-alpha_base], (*mp_arrays[c1])[i][alpha-alpha_base], tmp_mp);
                        LL_Double((*dp_arrays[c1]) + i, 1)[alpha-alpha_base] = mpf_get_d((*mp_arrays[c1])[i][alpha-alpha_base]);
                        mpf_clear((*mp_arrays[c1])[i][alpha-alpha_base]);
                }   }
                
                mantissa = mpf_get_d_2exp(&log2exp, tmp_mp);
                (*dp_norms[c1])[i] = log((double) mantissa) + log2exp*log(2);
                
                free((void *) (*mp_arrays[c1])[i]);       }
            
            free((void *) *mp_arrays[c1]);  }
        
        for (c1 = 0; c1 < 7; c1++)  mpf_clear(*mp_consts[c1]);
        
        for (i = 0; i < numLoci; i++)  {  mpf_clear(Z_mp[i]); mpf_clear(dC_dZ_mp[i]);  }
        for (alpha = 0; alpha < numMinParams; alpha++)  mpf_clear(grad_fw_mp[alpha]);
        free((void *) Z_mp);  free((void *) dC_dZ_mp);  free((void *) grad_fw_mp);
    }
    
    return passed;
}


// setMasks() sets the locus/image spot masks for a given Z term, and also writes the relevant locus->spot mappings
// (in case we're doing an enumerative calculation with no fixed overlaps).

ccBool setMasks()
{
    ccInt i, alpha, loopFixedLocus;
    
    for (i = 0; i < numLoci; i++)  locusMask[i] = ccTrue;
    for (alpha = 0; alpha < numSpots; alpha++)  spotMask[alpha] = ccTrue;
    
    for (loopFixedLocus = ZtermFirstLocus[loopZterm]; loopFixedLocus < ZtermFirstLocus[loopZterm+1]; loopFixedLocus++)  {
        ccInt oneFixedLocus = fixedLoci[loopFixedLocus];
        ccInt oneFixedSpot = fixedSpots[loopFixedLocus];
        
        if (!locusMask[oneFixedLocus])  return ccFalse;
        
        locusMask[oneFixedLocus] = ccFalse;
        if (oneFixedSpot >= 0)  spotMask[fixedSpots[loopFixedLocus]] = ccFalse;
    }
    
    return ccTrue;
}


// GetOptState() stores the state of the minimization (f (best guess), C), and returns true iff converged

ccBool getOptState(ccFloat *oneC)
{
    ccInt c1;
    gsl_vector *one_fw, *one_grad;
    
    oneC[0] = (ccFloat) gsl_multimin_fdfminimizer_minimum(opt_struct);
    one_grad = gsl_multimin_fdfminimizer_gradient(opt_struct);
    oneC[1] = (ccFloat) gsl_blas_dnrm2(one_grad);
    
    one_fw = gsl_multimin_fdfminimizer_x(opt_struct);
    for (c1 = 0; c1 < numMinParams; c1++)  {
        nlog_fw[c1] = (ccFloat) gsl_vector_get(one_fw, c1);
        mpDo( grad_fw[c1] = (ccFloat) gsl_vector_get(one_grad, c1); , mpf_set_d(grad_fw_mp[c1], gsl_vector_get(one_grad, c1)); )      }
    
    return (oneC[0] < 1.);
}


// Next six routines:  return C and/or grad-C (used by various optimization routines)

double getC(const gsl_vector *current_fw, void *dummy)
{
    load_fw(current_fw);
    
    lastC = IterateProbs(ccFalse);
    return (double) lastC;
}


void getGradC(const gsl_vector *current_fw, void *dummy, gsl_vector *grad_current_fw)
{
    load_fw(current_fw);
    
    IterateProbs(ccTrue);
    
    save_grad_fw(grad_current_fw);
}


void getCAndGradC(const gsl_vector *current_fw, void *dummy, double *C, gsl_vector *grad_current_fw)
{
    load_fw(current_fw);
    
    *C = lastC = (double) IterateProbs(ccTrue);
    
    save_grad_fw(grad_current_fw);
}


// load_f() converts (-log f) to f

void load_fw(const gsl_vector *current_fw)
{
    ccInt c1;
    
    for (c1 = 0; c1 < numMinParams; c1++)  {
        fw[c1] = exp(-gsl_vector_get(current_fw, c1));    }
}


// save_grad_f() converts dC/df to dC/d(-log f)

void save_grad_fw(const gsl_vector *grad_OF)
{
    ccInt c1;
    ccFloat one_grad_f;
    
    for (c1 = 0; c1 < numMinParams; c1++)  {
        mpDo( one_grad_f = -grad_fw[c1] * fw[c1]; , one_grad_f = -mpf_get_d(grad_fw_mp[c1]) * fw[c1]; )
        if (grad_OF != NULL)  gsl_vector_set((gsl_vector *) grad_OF, c1, (double) one_grad_f);      }
}



// Forward-propagates the Z matrices, and back-propagates the sensitivities.  Three modes:
// 0: do nothing -- just sum up the locus-wise and spot-wise normalizations and calculate the cost function
// 1: propagate Zs to calculate p, then do (0)
// 2: do (1 & 0), then calculate sensitivities

ccFloat IterateProbs(ccBool doGradient)
{
    ccInt i, alpha, alpha_base, loopColor;
    ccFloat *first_p, C;
    
    
            // Calculate Z and p
    
    forEachElement(doClearP);
    for (i = 0; i < numLoci; i++)  mpDo( Z[i] = 0.; , mpf_set_d(Z_mp[i], (double) 0.); )
    
    for (loopZterm = 0; loopZterm < numZterms; loopZterm++)  {
    if (setMasks())  {
        
        forEachPropagator(Z_1x, Z_1x_norm, Z_1x_mp, &Z_prop, 1, doFillZ);
        forEachPropagator(Z_xN, Z_xN_norm, Z_xN_mp, &Z_prop, -1, doFillZ);
        forEachElement(doRenormZ);
        forEachPropagator(Z_xN, Z_xN_norm, Z_xN_mp, &Z_bridge, -1, doBridgeZ);
        forEachElement(doCalcP);                        // p0 = Z_1x * Z_xN
    }}
    
    forEachElement(doDivZ);         // p = p0 / Z
    
    
            // Sum probs over the contour and image
    
    for (i = 0; i < numLoci; i++)  locusMappingProbs[i] = 0.;
    for (alpha = 0; alpha < numSpots; alpha++)  spotMappingProbs[alpha] = 0.;
    C = 0.;
    
    for (loopColor = 0; loopColor < numColors; loopColor++)  {
        
        ccFloat p_av, p_tot = 0.;
        alpha_base = colorFirstSpot[loopColor];
        
        for (i = 0; i < numLoci; i++)  {
        if (locusColor[i] == loopColor)  {
            
            if (hasSpots(i))  {
                first_p = LL_Double(p + i, 1);
                
                for (alpha = 0; alpha < colorFirstSpot[loopColor+1] - alpha_base; alpha++)  {
                    mpDo( , first_p[alpha] = mpf_get_d(p_mp[i][alpha]); )
                    
                    p_tot += first_p[alpha];
                    locusMappingProbs[i] += first_p[alpha];
                    spotMappingProbs[alpha+alpha_base] += first_p[alpha];
        }}  }   }
        p_av = p_tot / numLociPerColor[loopColor];
        
        N_difference[loopColor] = p_av - (1. - p_fn[loopColor]);
        C += 0.5 * N_difference[loopColor] * N_difference[loopColor] / max_pfn_mismatch / max_pfn_mismatch;
    }
    
    if (!setUnboundPenalty)  {
    for (alpha = 0; alpha < numSpots; alpha++)  {
    if (spotMappingProbs[alpha] > 1.)  {
        C += 0.5 * (spotMappingProbs[alpha] - 1.) * (spotMappingProbs[alpha] - 1.) / max_spot_overnorm / max_spot_overnorm;
    }}}
    
    
    if (!doGradient)  return C;
    
    
            // Compute the gradient dC/df and dC/dw
    
    for (i = 0; i < numLoci; i++)  mpDo( dC_dZ[i] = 0.; , mpf_set_d(dC_dZ_mp[i], 0.); )
    forEachElement(do_dCdZ);
    
    for (alpha = 0; alpha < numMinParams; alpha++)  mpDo( grad_fw[alpha] = 0.; , mpf_set_d(grad_fw_mp[alpha], 0.); )
    
    for (loopZterm = 0; loopZterm < numZterms; loopZterm++)  {
    if (setMasks())  {
        
        if (numZterms > 1)  {
            forEachPropagator(Z_1x, Z_1x_norm, Z_1x_mp, &Z_prop, 1, doFillZ);
            forEachPropagator(Z_xN, Z_xN_norm, Z_xN_mp, &Z_prop, -1, doFillZ);       }
        
        if (!setUnboundPenalty)  forEachPropagator(Z_xN, Z_xN_norm, Z_xN_mp, &Z_bridge, -1, doBridge_dCdZ);
        
        forEachElement(doSetS1x);
        forEachPropagator(Z_xN, Z_xN_norm, Z_xN_mp, &Z_bridge, -1, doBridgeS);
        forEachPropagator(sensitivity, sensitivity_norm, sensitivity_mp, &Z_prop, -1, doPropS);
        if (!setUnboundPenalty)  forEachElement(doAddS1x);
        if (setUnboundPenalty)  forEachPropagator(Z_1x, Z_1x_norm, Z_1x_mp, &Z_prop, 1, doGetGradW);
        
        forEachElement(doSetSxN);
        forEachPropagator(Z_1x, Z_1x_norm, Z_1x_mp, &Z_bridge, 1, doBridgeS);
        forEachPropagator(sensitivity, sensitivity_norm, sensitivity_mp, &Z_prop, 1, doPropS);
        if (!setUnboundPenalty)  forEachElement(doAddSxN);
        if (setUnboundPenalty)  forEachPropagator(Z_xN, Z_xN_norm, Z_xN_mp, &Z_prop, -1, doGetGradW);
    }}
    
    return C;
}


// forEachPropagator() fills the elements of arrays like Z_1x, etc. whose fill rule is like a layered neural network.

void forEachPropagator(linkedlist *prop_list, ccFloat *prop_list_norm, mpf_t **prop_list_mp, void(*NeighborFunction)(ccInt), int direction, allPropsArg mode)
{
    ccInt i, n_skip, one_color, loopColor, nSkip, boundary_overshoot, alpha, other_boundary, alpha_base, alpha_top, BICoffset;
    ccFloat *one_layer, *other_layer, *other_norm_list, *end_layer, *one_sense, normOffset, newNorm, otherNorm;
    ccBool setGradWSource;
    linkedlist *other_list;
    mpf_t **other_list_mp, *one_layer_mp, *other_layer_mp;
    
    if (direction == 1)  {  boundary = 0;  other_boundary = numLoci - 1; BICoffset = 0;  }
    else  {  boundary = numLoci - 1;  other_boundary = 0; BICoffset = -1;  }
    
    mpDo( ZtoProp = prop_list; , ZtoProp_mp = prop_list_mp; )
    propMode = mode;
    
    setGradWSource = ((mode == doBridgeS) && (direction < 0) && (setUnboundPenalty));
    
    boundary_overshoot = 0;
    if ((mode == doBridgeZ) || (mode == doBridge_dCdZ) || (mode == doBridgeS) || (mode == doGetGradW))
          {  if (mode != doGetGradW)  boundary_overshoot = 1;  }
    
    if ((prop_list == Z_1x) || ((propMode == doPropS) && (direction == 1)))
          {  other_list = Z_xN;  other_norm_list = Z_xN_norm;  other_list_mp = Z_xN_mp;  }
    else  {  other_list = Z_1x;  other_norm_list = Z_1x_norm;  other_list_mp = Z_1x_mp;  }
    
    for (i = boundary; (i+boundary_overshoot >= 0) && (i-boundary_overshoot < (ccInt) numLoci); i += direction)  {
    if (hasSpots(i))  {
        
        if ((i < 0) || (i >= numLoci))  {
            overBoundary = ccTrue;
            otherNorm = 0;
            alpha_base = 0;
            alpha_top = 1;     }
        
        else  {
            
            overBoundary = ccFalse;
            
            otherNorm = other_norm_list[i];
            
            mpDo( one_layer = LL_Double(prop_list + i, 1); , one_layer_mp = prop_list_mp[i]; )
            mpDo( other_layer = LL_Double(other_list + i, 1); , other_layer_mp = other_list_mp[i]; )
            
            one_sense = LL_Double(sensitivity + i, 1);
            one_color = locusColor[i];
            
            alpha_base = colorFirstSpot[one_color];
            alpha_top = colorFirstSpot[one_color+1];            }
        
        
        for (alpha = alpha_base; alpha < alpha_top; alpha++)  {
            
            ccInt alphaIdx = alpha - alpha_base;
            ccFloat w_factor = 0.;
            
            if (mode == doFillZ)  { mpDo( one_layer[alphaIdx] = 0.; , mpf_set_d(one_layer_mp[alphaIdx], 0.); ) }
            
            if (mappingIsAllowed(i, alpha, -1, ccFalse))  {
                alpha_0 = alpha;
                
                for (loopColor = 0; loopColor < numColors; loopColor++)  numSkippedColors[loopColor] = 0;
                
                prop_hadContourBreak = ccFalse;
                for (n_skip = 0; n_skip <= n_skip_max; n_skip++)  {
                    
                    if (abs(i - boundary) < n_skip)  break;
                    
                    propLocus = i - (n_skip + 1) * direction;
                    bridgingLocus = i;
                    
                    prop_hadContourBreak = prop_hadContourBreak || breakInContour[propLocus+BICoffset];
                    
                    if (hasSpots(propLocus))  {
                        
                        mpDo( , mpf_set_d(w_factor_mp, exp(w_factor)); )
                        
                        if (abs(i - boundary) == n_skip)  {
                            mpDo( Zpropagated = 1.; , mpf_set_d(p_value_mp, 1.); )
                            newNorm = 0.;       }
                        
                        else  {
                            mpDo( Zpropagated = 0.; , mpf_set_d(p_value_mp, 0.); )
                            propDl = fabs( (double) (l[propLocus] - l[i]) );
                            global_alpha_base = colorFirstSpot[locusColor[propLocus]];
                            
                            newNorm = prop_list_norm[propLocus];        }
                        
                        if (NeighborFunction == &Z_bridge)  {
                            
                            if (overBoundary)  mpDo( multiplier = 1.; , mpf_set_d(multiplier_mp, 1.); )
                            else if (mode == doBridgeS)  mpDo( multiplier = sqrt(f[alpha]); , mpf_set_d(multiplier_mp, sqrt(f[alpha])); )
                            else  mpDo( multiplier = sqrt(f[alpha]) * (LL_Double(Z_1x + i, 1)[alphaIdx]); ,
                                mpf_set_d(multiplier_mp, sqrt(f[alpha]));
                                mpf_mul(multiplier_mp, multiplier_mp, Z_1x_mp[i][alphaIdx]);  )
                            if (abs(i - boundary) == n_skip)  { mpDo( Zpropagated = multiplier; , mpf_set(p_value_mp, multiplier_mp); ) }
                            
                            if (mode == doBridge_dCdZ)  {
                                mpDo( dC_dZ_factor = 0.; , mpf_set_d(dC_dZ_factor_mp, 0.); )
                                for (nSkip = 1; nSkip <= n_skip; nSkip++)   {
                                    mpDo( dC_dZ_factor += 0.5 * ZtermWeights[loopZterm] * dC_dZ[i - nSkip*direction] * exp( w_factor + otherNorm + newNorm - Z_norm[i - nSkip*direction] ); ,
                                        mpf_set_d(tmp_mp, 0.5 * ZtermWeights[loopZterm]);
                                        mpf_mul(tmp_mp, tmp_mp, dC_dZ_mp[i - nSkip*direction]);
                                        mpf_mul(tmp_mp, tmp_mp, w_factor_mp);
                                        mpf_add(dC_dZ_factor_mp, dC_dZ_factor_mp, tmp_mp); )
                        }   }   }
                        
                        if (abs(i - boundary) != n_skip)
                            forEachNeighbor(propLocus, alpha, propDl, locusColor[propLocus], NeighborFunction);
                        
                        if (mode == doFillZ)  {
                            ccFloat nextNorm = 0.;
                            if (i != boundary)  nextNorm = prop_list_norm[i - direction];
                            mpDo( one_layer[alphaIdx] += Zpropagated * sqrt(f[alpha]) * exp(w_factor + newNorm - nextNorm); ,
                                mpf_set_d(tmp_mp, sqrt(f[alpha]));
                                mpf_mul(tmp_mp, tmp_mp, p_value_mp);
                                mpf_mul(tmp_mp, tmp_mp, w_factor_mp);
                                mpf_add(one_layer_mp[alphaIdx], one_layer_mp[alphaIdx], tmp_mp); )
                        }
                        
                        else if (mode == doGetGradW)  {
                            for (loopColor = 0; loopColor < numColors; loopColor++)   {
                                mpDo( grad_fw[loopColor] += Zpropagated * numSkippedColors[loopColor] * sqrt(f[alpha])
                                        * one_sense[alphaIdx] * exp(w_factor + newNorm + sensitivity_norm[i]) / w[loopColor]; ,
                                    mpf_set_d(tmp_mp, numSkippedColors[loopColor] * sqrt(f[alpha]) / w[loopColor]);
                                    mpf_mul(tmp_mp, tmp_mp, p_value_mp);
                                    mpf_mul(tmp_mp, tmp_mp, sensitivity_mp[i][alphaIdx]);
                                    mpf_mul(tmp_mp, tmp_mp, w_factor_mp);
                                    mpf_add(grad_fw_mp[loopColor], grad_fw_mp[loopColor], tmp_mp); )
                        }  }
                        
                        else if (mode == doBridgeZ)  {
                            for (nSkip = 1; nSkip <= n_skip; nSkip++)  {
                                mpDo( Z[i - nSkip*direction] += ZtermWeights[loopZterm] * Zpropagated * exp( w_factor + otherNorm + newNorm - Z_norm[i - nSkip*direction] ); ,
                                    mpf_set_d(tmp_mp, ZtermWeights[loopZterm]);
                                    mpf_mul(tmp_mp, tmp_mp, p_value_mp);
                                    mpf_mul(tmp_mp, tmp_mp, w_factor_mp);
                                    mpf_add(Z_mp[i - nSkip*direction], Z_mp[i - nSkip*direction], tmp_mp); )
                        }   }
                        
                        else if ((mode == doBridge_dCdZ) && (!overBoundary))  {
                            mpDo( grad_fw[alpha] += Zpropagated * dC_dZ_factor / f[alpha]; ,
                                mpf_set_d(tmp_mp, f[alpha]);
                                mpf_div(tmp_mp, p_value_mp, tmp_mp);
                                mpf_mul(tmp_mp, tmp_mp, dC_dZ_factor_mp);
                                mpf_add(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp); )
                        }
                        
                        else if (mode == doBridgeS)  {
                            ccFloat to_add;
                            
                            for (nSkip = 1; nSkip <= n_skip; nSkip++)  {
                                
                                if (!overBoundary)  {
                                    mpDo(
                                        end_layer = LL_Double(sensitivity + i, 1);
                                        
                                        to_add = ZtermWeights[loopZterm] * dC_dZ[i - nSkip*direction] * Zpropagated
                                                    * exp( w_factor + newNorm - Z_norm[i - nSkip*direction] - sensitivity_norm[i] );
                                        
                                        end_layer[alphaIdx] += to_add;     // next line:  -(1/2) Z^i (dC/dZ^i) / f_i  (or with Z_i)
                                        if (!setUnboundPenalty)  grad_fw[alpha] -= 0.5 * other_layer[alphaIdx] * to_add
                                                         * exp( otherNorm + sensitivity_norm[i] ) / f[alpha];  ,
                                        
                                        mpf_set_d(tmp_mp, ZtermWeights[loopZterm]);
                                        mpf_mul(tmp_mp, tmp_mp, dC_dZ_mp[i - nSkip*direction]);
                                        mpf_mul(tmp_mp, tmp_mp, p_value_mp);
                                        mpf_mul(tmp_mp, tmp_mp, w_factor_mp);
                                        
                                        mpf_add(sensitivity_mp[i][alphaIdx], sensitivity_mp[i][alphaIdx], tmp_mp);
                                        if (!setUnboundPenalty)  {
                                            mpf_set_d(tmp2_mp, 0.5 / f[alpha]);
                                            mpf_mul(tmp_mp, tmp_mp, tmp2_mp);
                                            mpf_mul(tmp_mp, tmp_mp, other_layer_mp[alphaIdx]);
                                            mpf_sub(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp);       }
                                )   }
                                
                                if (setGradWSource)  {
                                    mpDo( to_add = ZtermWeights[loopZterm] * dC_dZ[i - nSkip*direction]
                                            * Zpropagated * exp( w_factor + newNorm + otherNorm - Z_norm[i - nSkip*direction] ); ,
                                        mpf_set_d(tmp_mp, ZtermWeights[loopZterm]);
                                        mpf_mul(tmp_mp, tmp_mp, dC_dZ_mp[i - nSkip*direction]);
                                        mpf_mul(tmp_mp, tmp_mp, p_value_mp);
                                        mpf_mul(tmp_mp, tmp_mp, w_factor_mp); )
                                    if (!overBoundary)  { mpDo( to_add *= LL_Double(Z_1x + i, 1)[alphaIdx]; , mpf_mul(tmp_mp, tmp_mp, Z_1x_mp[i][alphaIdx]); ) }
                                    for (loopColor = 0; loopColor < numColors; loopColor++)   {
                                        mpDo( grad_fw[loopColor] += numSkippedColors[loopColor] * to_add / w[loopColor]; ,
                                            mpf_set_d(tmp2_mp, ((double) numSkippedColors[loopColor]) / w[loopColor]);
                                            mpf_mul(tmp2_mp, tmp_mp, tmp2_mp);
                                            mpf_add(grad_fw_mp[loopColor], grad_fw_mp[loopColor], tmp2_mp); )
                        }   }   }   }
                        
                        else if ((mode == doPropS) && (abs(i - boundary) != n_skip))  {
                            mpDo( one_layer[alphaIdx] += Zpropagated * sqrt(f[alpha]) * exp(w_factor + prop_list_norm[propLocus] - prop_list_norm[i]); ,
                                mpf_set_d(tmp_mp, sqrt(f[alpha]));
                                mpf_mul(tmp_mp, tmp_mp, p_value_mp);
                                mpf_mul(tmp_mp, tmp_mp, w_factor_mp);
                                mpf_add(one_layer_mp[alphaIdx], one_layer_mp[alphaIdx], tmp_mp); )
                    }   }
                    
                    if (!mappingIsAllowed(propLocus, -1, -1, ccFalse))  break;
                    
                    if (abs(i - boundary) != n_skip)  {
                        numSkippedColors[locusColor[propLocus]]++;
                        w_factor += log(w[locusColor[propLocus]]);
                }   }
        }   }
        
            // renormalize the Z/s arrays
        
        mpDo(
            if (mode == doFillZ)  {
                if (i == boundary)  normOffset = 0.;
                else  normOffset = prop_list_norm[i - direction];
                RenormZ(one_layer, alpha_top - alpha_base, prop_list_norm + i, normOffset);            }
            
            else if (mode == doPropS)  {
                normOffset = prop_list_norm[i];
                RenormZ(one_layer, alpha_top - alpha_base, prop_list_norm + i, normOffset);            } , )
    }}
}



// forEachElement() iterates some local operation over each element of some array (i.e. doesn't propagate information between loci)

void forEachElement(allElsArg mode)
{
    ccInt i, alpha, alpha_base, one_color;
    ccFloat *first_p, *one_Z_1x, *one_Z_xN, *one_sense, one_Z_term, dC_dp, pScaleFactor;
    mpf_t *first_p_mp, *one_Z_1x_mp, *one_Z_xN_mp, *one_sense_mp;
    
    for (i = 0; i < numLoci; i++)  {
    if (hasSpots(i))  {
        
        ccFloat normConversionFactor = exp(Z_1x_norm[i] + Z_xN_norm[i] - Z_norm[i]);
        
        mpDo( first_p = LL_Double(p + i, 1); , first_p_mp = p_mp[i]; )
        mpDo( one_Z_1x = LL_Double(Z_1x + i, 1); , one_Z_1x_mp = Z_1x_mp[i]; )
        mpDo( one_Z_xN = LL_Double(Z_xN + i, 1); , one_Z_xN_mp = Z_xN_mp[i]; )
        mpDo( one_sense = LL_Double(sensitivity + i, 1); , one_sense_mp = sensitivity_mp[i]; )
        one_color = locusColor[i];
        
        if (mode == doRenormZ)  {
            mpDo( ccFloat thisPNorm = Z_1x_norm[i] + Z_xN_norm[i];
            
            if (loopZterm == 0)  {  pScaleFactor = 0.;  Z_norm[i] = thisPNorm;  }
            else if (thisPNorm <= Z_norm[i])  pScaleFactor = 1.;
            else  {  pScaleFactor = exp(Z_norm[i] - thisPNorm);  Z_norm[i] = thisPNorm;  }
            
            Z[i] *= pScaleFactor; , )
        }
        
        alpha_base = colorFirstSpot[one_color];
        
        for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)  {
            
            ccInt alphaIdx = alpha-alpha_base;
            ccFloat dC_dp_p_fn_rate = 1.;
            
            if (extraOcc[one_color] != 0.)  dC_dp_p_fn_rate -= 1. / (extraOcc[one_color]*(extraOcc[one_color]-1.));
            dC_dp_p_fn_rate *= N_difference[one_color] / numLociPerColor[one_color] / max_pfn_mismatch / max_pfn_mismatch;
            
            mpDo( dC_dp = dC_dp_p_fn_rate; , mpf_set_d(dC_dp_mp, (double) dC_dp_p_fn_rate); )
            if ((!setUnboundPenalty) && (spotMappingProbs[alpha] > 1.))  {
                mpDo( dC_dp += (spotMappingProbs[alpha] - 1.) / max_spot_overnorm / max_spot_overnorm; ,
                    mpf_set_d(tmp_mp, (double) ((spotMappingProbs[alpha] - 1.) / max_spot_overnorm / max_spot_overnorm)); mpf_add(dC_dp_mp, dC_dp_mp, tmp_mp); )       }
            
            if (mode == doClearP)  {
                mpDo( first_p[alphaIdx] = 0.; , mpf_set_d(first_p_mp[alphaIdx], (double) 0.); )        }
            
            if (mode == doRenormZ)  {
                mpDo(first_p[alphaIdx] *= pScaleFactor; , )        }
            
            else if (mode == doCalcP)  {
                mpDo( one_Z_term = ZtermWeights[loopZterm] * one_Z_1x[alphaIdx] * one_Z_xN[alphaIdx] * normConversionFactor; ,
                    mpf_set_d(tmp_mp, (double) ZtermWeights[loopZterm]);
                    mpf_mul(tmp_mp, tmp_mp, one_Z_1x_mp[alphaIdx]);
                    mpf_mul(tmp_mp, tmp_mp, one_Z_xN_mp[alphaIdx]); )
                mpDo( Z[i] += one_Z_term; , mpf_add(Z_mp[i], Z_mp[i], tmp_mp); )
                mpDo( first_p[alphaIdx] += one_Z_term; , mpf_add(first_p_mp[alphaIdx], first_p_mp[alphaIdx], tmp_mp); )       }
            
            else if (mode == doDivZ)  {
                mpDo( if (Z[i] != 0.)  first_p[alphaIdx] /= Z[i]; , if (mpf_sgn(Z_mp[i]) != 0)  mpf_div(first_p_mp[alphaIdx], first_p_mp[alphaIdx], Z_mp[i] ); )     }
            
            else if (mode == do_dCdZ)  {
                mpDo( if (Z[i] != 0.)  dC_dZ[i] -= dC_dp * first_p[alphaIdx] / Z[i]; ,
                    if (mpf_sgn(Z_mp[i]) != 0)  {  mpf_mul(tmp_mp, dC_dp_mp, first_p_mp[alphaIdx]);  mpf_div(tmp_mp, tmp_mp, Z_mp[i]);  mpf_sub(dC_dZ_mp[i], dC_dZ_mp[i], tmp_mp);  } )     }
            
            else if (mode == doSetS1x)  {
                mpDo(
                    ccFloat to_add = dC_dZ[i];
                    if (Z[i] != 0.)
                        to_add += dC_dp / Z[i];
                    to_add *= ZtermWeights[loopZterm] * normConversionFactor * one_Z_xN[alphaIdx];
                    
                    one_sense[alphaIdx] = to_add;
                    if (!setUnboundPenalty)  grad_fw[alpha] -= 0.5 * one_Z_1x[alphaIdx] * to_add / f[alpha];  ,
                    
                    mpf_set(one_sense_mp[alphaIdx], dC_dZ_mp[i]);
                    if (mpf_sgn(Z_mp[i]) != 0)  {
                        mpf_div(tmp_mp, dC_dp_mp, Z_mp[i]);
                        mpf_add(one_sense_mp[alphaIdx], one_sense_mp[alphaIdx], tmp_mp);   }
                    
                    mpf_set_d(tmp_mp, (double) ZtermWeights[loopZterm]);
                    mpf_mul(tmp_mp, tmp_mp, one_Z_xN_mp[alphaIdx]);
                    mpf_mul(one_sense_mp[alphaIdx], one_sense_mp[alphaIdx], tmp_mp);
                    
                    if (!setUnboundPenalty)  {
                        mpf_mul(tmp_mp, one_sense_mp[alphaIdx], one_Z_1x_mp[alphaIdx]);
                        mpf_set_d(tmp2_mp, 0.5);
                        mpf_mul(tmp_mp, tmp_mp, tmp2_mp);
                        mpf_set_d(tmp2_mp, f[alpha]);
                        mpf_div(tmp_mp, tmp_mp, tmp2_mp);
                        mpf_sub(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp);      }
            )   }
            
            else if (mode == doAddS1x)  {
                mpDo( grad_fw[alpha] += one_Z_1x[alphaIdx] * one_sense[alphaIdx] * exp(Z_1x_norm[i] + sensitivity_norm[i]) / f[alpha]; ,
                    mpf_set_d(tmp_mp, (double) f[alpha]);
                    mpf_div(tmp_mp, one_Z_1x_mp[alphaIdx], tmp_mp);
                    mpf_mul(tmp_mp, tmp_mp, one_sense_mp[alphaIdx]);
                    mpf_add(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp); )        }
            
            else if (mode == doSetSxN)  {
                mpDo(
                    ccFloat to_add = dC_dZ[i];
                    if (Z[i] != 0.)
                        to_add += dC_dp / Z[i];
                    to_add *= ZtermWeights[loopZterm] * normConversionFactor * one_Z_1x[alphaIdx];
                    
                    one_sense[alphaIdx] = to_add;
                    if (!setUnboundPenalty)  grad_fw[alpha] -= 0.5 * one_Z_xN[alphaIdx] * to_add / f[alpha]; ,
                    
                    mpf_set(one_sense_mp[alphaIdx], dC_dZ_mp[i]);
                    if (mpf_sgn(Z_mp[i]) != 0)  {
                        mpf_div(tmp_mp, dC_dp_mp, Z_mp[i]);
                        mpf_add(one_sense_mp[alphaIdx], one_sense_mp[alphaIdx], tmp_mp);   }
                    
                    mpf_set_d(tmp_mp, (double) ZtermWeights[loopZterm]);
                    mpf_mul(tmp_mp, tmp_mp, one_Z_1x_mp[alphaIdx]);
                    mpf_mul(one_sense_mp[alphaIdx], one_sense_mp[alphaIdx], tmp_mp);
                    
                    if (!setUnboundPenalty)  {
                        mpf_mul(tmp_mp, one_sense_mp[alphaIdx], one_Z_xN_mp[alphaIdx]);
                        mpf_set_d(tmp2_mp, 0.5);
                        mpf_mul(tmp_mp, tmp_mp, tmp2_mp);
                        mpf_set_d(tmp2_mp, f[alpha]);
                        mpf_div(tmp_mp, tmp_mp, tmp2_mp);
                        mpf_sub(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp);      }
            )   }
            
            else if (mode == doAddSxN)  {
                mpDo( grad_fw[alpha] += one_Z_xN[alphaIdx] * one_sense[alphaIdx] * exp(Z_xN_norm[i] + sensitivity_norm[i]) / f[alpha]; ,
                    mpf_set_d(tmp_mp, (double) f[alpha]);
                    mpf_div(tmp_mp, one_Z_xN_mp[alphaIdx], tmp_mp);
                    mpf_mul(tmp_mp, tmp_mp, one_sense_mp[alphaIdx]);
                    mpf_add(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp); )        }
        }
    }}
    
        // set the initial normalizations of the Z/s arrays
    
    if (mode == doSetS1x)  {
    for (i = 0; i < numLoci; i++)  {
        sensitivity_norm[i] = -Z_1x_norm[i];     }}
    else if (mode == doSetSxN)  {
    for (i = 0; i < numLoci; i++)  {
        sensitivity_norm[i] = -Z_xN_norm[i];     }}
}



// RenormZ() rescales one column of a Z/s array along with its normalization factor, stored in a separate array.

void RenormZ(ccFloat *firstZ, ccInt numZs, ccFloat *norm, ccFloat normOffset)
{
    ccInt c1;
    ccFloat tot;
    
    tot = 0.;
    for (c1 = 0; c1 < numZs; c1++)  tot += fabs((double) firstZ[c1]);
    
    if (tot == 0.)  {  *norm = normOffset;  return;  }
    
    for (c1 = 0; c1 < numZs; c1++)  firstZ[c1] /= tot;
    *norm = log(tot) + normOffset;
}



// Z_prop() -- the basic routine propagating terms in the Z/s arrays

void Z_prop(ccInt alpha)
{
    ccFloat GP_rootf = sqrt(f[alpha]);
    if (!prop_hadContourBreak)  GP_rootf *= GaussProb(propDl, alpha_0, alpha, 0);
    
    mpDo( Zpropagated += LL_Double(ZtoProp + propLocus, 1)[alpha - global_alpha_base] * GP_rootf;  ,
        
        mpf_set_d(tmp_mp, (double) GP_rootf);
        mpf_mul(tmp_mp, tmp_mp, ZtoProp_mp[propLocus][alpha - global_alpha_base]);
        mpf_add(p_value_mp, p_value_mp, tmp_mp);    )
}


// Z_bridge() introduces terms coming from false negatives

void Z_bridge(ccInt alpha)
{
    if (abs(bridgingLocus - propLocus) > 1)  {
        
        mpDo(
            ccFloat to_add = LL_Double(ZtoProp + propLocus, 1)[alpha - global_alpha_base] * sqrt(f[alpha]) * multiplier;
            if ((!overBoundary) && (!prop_hadContourBreak))  to_add *= GaussProb(propDl, alpha_0, alpha, 0);
            
            if (propMode == doBridge_dCdZ)  {
                grad_fw[alpha] += to_add * dC_dZ_factor / f[alpha];    }
            Zpropagated += to_add;  ,
            
            mpf_set_d(tmp_mp, (double) sqrt(f[alpha]));
            mpf_mul(tmp_mp, tmp_mp, multiplier_mp);
            mpf_mul(tmp_mp, tmp_mp, ZtoProp_mp[propLocus][alpha - global_alpha_base]);
            if ((!overBoundary) && (!prop_hadContourBreak))  {
                mpf_set_d(tmp2_mp, (double) GaussProb(propDl, alpha_0, alpha, 0));
                mpf_mul(tmp_mp, tmp_mp, tmp2_mp);       }
            
            mpf_add(p_value_mp, p_value_mp, tmp_mp);
            if (propMode == doBridge_dCdZ)  {
                mpf_mul(tmp_mp, tmp_mp, dC_dZ_factor_mp);
                mpf_set_d(tmp2_mp, f[alpha]);
                mpf_div(tmp_mp, tmp_mp, tmp2_mp);
                mpf_add(grad_fw_mp[alpha], grad_fw_mp[alpha], tmp_mp);    }
    )   }
}


void CountNeighbors(ccInt alpha)
{
    mpDo( Zpropagated += 1.; , mpf_set_d(tmp_mp, 1.); mpf_add(p_value_mp, p_value_mp, tmp_mp); )
}

#undef mpDo




// GetAllChains() performs an exact calculation of the partition function (or simulates the align3d approximate calculation if noOverlapsAtAll == false)
// by explicitly enumerating each conformation.  Uses the mask to keep track of which spots have been used in the current chain.

void GetAllChains(ccInt prevLocus, ccInt i, ccFloat pMapping, ccBool noOverlapsAtAll, ccBool hadContourBreak)
{
    if (i == numLoci)  {
        
        ccInt j;
        
        *Z += pMapping;
        
        for (j = 0; j < numLoci; j++)  {
            
            ccInt j_color = locusColor[j];
            ccInt j_alpha_base = (ccInt) colorFirstSpot[j_color];
            ccInt j_alpha = alphas[j];
            
            if (j_alpha >= j_alpha_base)  {
//printf(" %i->%i", j+1, j_alpha+1);
                LL_Double(p + j, 1)[j_alpha - j_alpha_base] += pMapping;
        }   }
//printf(" (%g)\n", pMapping);
    }
    
    else  {
        
        ccInt one_color = locusColor[i];
        ccInt alpha_base = colorFirstSpot[one_color], alpha;
        
        alphas[i] = alpha_base-1;
        if (i - prevLocus <= n_skip_max)  {
        if (mappingIsAllowed(i, -1, -1, ccFalse))  {
            GetAllChains(prevLocus, i+1, pMapping * w[one_color], noOverlapsAtAll, hadContourBreak || breakInContour[i]);
        }}
        
        for (alpha = alpha_base; alpha < (ccInt) colorFirstSpot[one_color+1]; alpha++)  {
        if (mappingIsAllowed(i, alpha, prevLocus, noOverlapsAtAll))  {
            
            ccFloat extra_p = 1.;
            if (!hadContourBreak)  {
                ccFloat dl = fabs( (double) (l[i] - l[prevLocus]) );
                extra_p = GaussProb(dl, alpha, alphas[prevLocus], 0);       }
            
            alphas[i] = alpha;
            spotWasUsed[alpha] = ccTrue;
            GetAllChains(i, i+1, pMapping * extra_p, noOverlapsAtAll, breakInContour[i]);
            spotWasUsed[alpha] = ccFalse;
    }   }}
}





// *********** Misc ***********



int call_GaussianChain(int argc, char **argv)
{
    ccInt FD1, FD2, mode, dummy_calcs;
    ccFloat L, *result;
    arg_info *argInfo = (arg_info *) *(argv+argc);
    
    getArgs(argc, argv, byValue(&mode), &spotX, &spotY, &spotZ, &spotDx, &spotDy, &spotDz,
                byValue(&L), &model, byValue(&FD1), byValue(&FD2), &result);
    
    if ((FD1 == 0) || (FD2 == 0) || (FD1 > numSpots) || (FD2 > numSpots))  {
        printf("GaussianChain():  image spot 1 or 2 out of range\n");
        return 2;               }
    
    exaggeration = 1.;
    calcs = &dummy_calcs;
    numModelRegimes = argInfo[8].argIndices / 3;
    
    *result = GaussProb(L, FD1-1, FD2-1, mode);
    
    return passed;
}



ccFloat GaussProb(ccFloat dl, ccInt imageSpot1, ccInt imageSpot2, char GN_mode)
{
    ccInt loopModelRegime;
    ccFloat invDecayConstant2, ax, ay, az, ans, expectedDistance, *modelCopy = model;
    const ccFloat pi3 = pi*pi*pi;
    
    (*calcs)++;
    
    for (loopModelRegime = 1; loopModelRegime < numModelRegimes; loopModelRegime++)  {
        if (dl < modelCopy[3])  break;
        modelCopy += 3;     }
    
    expectedDistance = modelCopy[1] * pow(dl, modelCopy[2]);
    invDecayConstant2 = (2*expectedDistance*expectedDistance)/3;
    
    ax = 1./(invDecayConstant2 + 2*spotDx[imageSpot1]*spotDx[imageSpot1] + 2*spotDx[imageSpot2]*spotDx[imageSpot2]);
    ay = 1./(invDecayConstant2 + 2*spotDy[imageSpot1]*spotDy[imageSpot1] + 2*spotDy[imageSpot2]*spotDy[imageSpot2]);
    az = 1./(invDecayConstant2 + 2*spotDz[imageSpot1]*spotDz[imageSpot1] + 2*spotDz[imageSpot2]*spotDz[imageSpot2]);
    
    if (GN_mode == 0)
        ans = sqrt(ax*ay*az/pi3)*exp(-sqSpotDistance(imageSpot1, imageSpot2, ax, ay, az));
    else if (GN_mode == 1)
        ans = exp(-sqSpotDistance(imageSpot1, imageSpot2, ax, ay, az));
    else
        ans = sqrt(ax*ay*az/pi3);
    
    if (exaggeration != 1.)
        ans = pow(ans, exaggeration);
    
    return ans;
}



ccFloat sqSpotDistance(ccInt spot1, ccInt spot2, ccFloat wx, ccFloat wy, ccFloat wz)
{
    ccFloat dx, dy, dz;
    
    dx = spotX[spot1] - spotX[spot2];
    dy = spotY[spot1] - spotY[spot2];
    dz = spotZ[spot1] - spotZ[spot2];
    
    return wx*dx*dx + wy*dy*dy + wz*dz*dz;
}







// *********** Diagnostics ***********


int call_Entropy(int argc, char **argv)
{
    ccInt c1, alpha_base, alpha, one_color, *C2F;
    arg_info *argInfo = (arg_info *) *(argv+argc);
    ccFloat *info, prob_sum, *first_p, one_p, x_res, y_res, z_res, alpha_x, alpha_y, alpha_z, one_tot_p;
    ccBool IfAvg, IfZeroNorm, countFalseNegatives;
    
    getArgs(argc, argv, &spotX, &spotY, &spotZ, &p, &locusColor, &colorFirstSpot, &C2F,
                byValue(&x_res), byValue(&y_res), byValue(&z_res), byValue(&countFalseNegatives), byValue(&IfAvg), &info);
    
    numLoci = argInfo[3].argIndices;
    
    IfZeroNorm = ccFalse;
    if (x_res*y_res*z_res == 0)  IfZeroNorm = ccTrue;
    else  {
        alpha_x = 1. / (2 * x_res * x_res);
        alpha_y = 1. / (2 * y_res * y_res);
        alpha_z = 1. / (2 * z_res * z_res);         }
    
    *info = 0;
    for (c1 = 0; c1 < numLoci; c1++)  {
    if (hasSpots(c1))  {
        
        first_p = LL_Double(p + c1, 1);
        one_color = locusColor[c1];
        
        prob_sum = one_tot_p = 0.;
        alpha_base = colorFirstSpot[one_color];
        for (alpha = alpha_base; alpha < colorFirstSpot[one_color+1]; alpha++)  {
            one_p = first_p[alpha - alpha_base];
            prob_sum += one_p;
            
            if (IfAvg)  {
                if (one_p != 0.)  *info -= one_p*log(one_p);    }
            
            else  {
                if ((!IfZeroNorm) && (C2F[c1] >= 1))
                    one_tot_p += one_p * exp(-sqSpotDistance(alpha, C2F[c1] - 1, alpha_x, alpha_y, alpha_z));
                else if ((IfZeroNorm) && (alpha == C2F[c1] - 1))
                    one_tot_p += one_p;
        }   }
        
        if ((countFalseNegatives) && (prob_sum >= 0.))  {
            if (IfAvg)  {
                if (prob_sum < 1.)  *info -= (1. - prob_sum)*log(1. - prob_sum);         }
            else  {
                if (C2F[c1] == 0)  *info -= log(1. - prob_sum);
                else if (one_tot_p < 1.)  *info -= log(one_tot_p);        // screen out the > 1 case
        }   }
    }}
    
    return passed;
}


// mappingIsAllowed() returns true if either the locus/spot are both free (not fixed), or else the locus is fixed to map to the spot.
// If the spot is -1, then returns true if the locus is either free or fixed to a false negative. 

ccBool mappingIsAllowed(ccInt i, ccInt alpha, ccInt prevLocus, ccBool noOverlapsAtAll)
{
    ccInt loopFixedLocus;
    
    if ((i < 0) || (i >= numLoci))  return ccTrue;
    
    if (noOverlapsAtAll)  {
        if (spotWasUsed[alpha])  return ccFalse;        }
    else if (prevLocus >= 0)  {
        if (alpha == alphas[prevLocus])  return ccFalse;       }
    
    if (alpha >= 0)  {
        if ((!avoidConstrainedSpot[alpha]) && (locusMask[i]))  return ccTrue;
        if (locusMask[i] != spotMask[alpha])  return ccFalse;    }
    if (locusMask[i])  return ccTrue;
    
    for (loopFixedLocus = ZtermFirstLocus[loopZterm]; loopFixedLocus < ZtermFirstLocus[loopZterm+1]; loopFixedLocus++)  {
        if (fixedLoci[loopFixedLocus] == i)  {
            if (alpha < 0)  return (fixedSpots[loopFixedLocus] < 0);
            else if (fixedSpots[loopFixedLocus] == alpha)  return ccTrue;
    }   }
    
    return ccFalse;
}



// hasSpots() returns false only if a locus has a color that no spots in the image have

ccBool hasSpots(ccInt i)
{
    if ((i < 0) || (i >= numLoci))  return ccTrue;
    
    return (colorFirstSpot[locusColor[i]+1] - colorFirstSpot[locusColor[i]] > 0);
}


// clock() returns the system time in seconds since startup

int call_clock(int argc, char **argv)
{
    ccFloat *currentTime;
    
    getArgs(argc, argv, &currentTime);
    
    *currentTime = ((ccFloat) clock()) / CLOCKS_PER_SEC;
    
    return passed;
}
