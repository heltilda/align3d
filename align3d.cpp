/*
 *  align3d.c
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
#include "align3d.h"
#include "LinMath.h"
#include "Interpolation.h"
#include "MonteCarlo.h"
#include "lnklst.h"
#include "intrpt.h"
#include "userfn.h"


double tiny = 1.e-15, tinylog = -100;

typedef struct {
    double *x, *y, *z, *dx, *dy, *dz;
    ccInt numSpots, *color_bottoms, *color_tops;
} imageType;

typedef struct {
    double *l;
    int *color, *mask;
    ccInt numSpots;
} contour_type;

typedef struct {
    linkedlist *p, *Z_1x, *Z_xN, *prop_u_or_d, *Neighbors, *sensitivity;
    ccInt colors_num, neighbor_lengths_num, alpha_0, alpha_base, n_skip_max, alpha_0_base, f_or_w;
    double *CSumProbs, *FSumProbs, *Z_1x_norm, *Z_xN_norm, *Z, *dC_dZ, *sensitivity_norm, *calcs, *prop_u_or_d_norm;
    double *f, *grad_f, multiplier, dC_dZ_factor, lp, l_step, p_value, l, p_fn, fd_weight, other_norm;
    double NOF_length, w, K1, N_difference, *other_layer, w_norm, w_factor, exaggeration;
    ccInt boundary, direction, cSpot_0, cSpot, *alphas;
    ccBool over_boundary, prop_mode, pa_SetGradWSource, *alpha_mask;
double *pLR;
} a3d_params_type;

imageType image;
contour_type contour;
a3d_params_type a3d_params;

gsl_multimin_fminimizer *opt_struct_no_grad;
gsl_multimin_fdfminimizer *opt_struct_grad;
const gsl_rng_type *gsl_rand_gen_type;

const char pa_FillZ = 0, fa_CalcP = 0;
const char pa_BridgeZ = 1, fa_DivZ = 1;
const char pa_BridgeS = 2, fa_dC_dZ = 2;
const char pa_PropS = 3;
const char pa_Bridge_dC_dZ = 4;
const char pa_GetGradW = 5;
const char fa_SetS1x = 4;
const char fa_AddS1x = 5;
const char fa_SetSxN = 6;
const char fa_AddSxN = 7;

ccBool dbg = ccFalse;
ccInt dbg1, dbg2;

extern void do_pLR(ccInt);





// *********** GetNeigbors() ***********

// Finds the neighbors of each spot and stores them in the Neighbors linked lists

int call_GetNeighbors(int argc, char **argv)
{
    ccInt c1, rtrn;
    double pCutoff;
	arg_info *ArgInfo = (arg_info *) *(argv+argc);
    
	const int ArgTypes[] = { double_type, double_type, double_type, double_type, double_type, double_type,
                                                int_type, int_type, string_type, double_type, double_type, double_type };
	const int ArgIndices[] = { -1, -1, -1, -1, -1, -1, -2, -2, -3, 1, 1, 1 };
	
	if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "GetNeighbors") != passed)  return 1;
	
    getArgs(argc, argv, &image.x, &image.y, &image.z, &image.dx, &image.dy, &image.dz, &image.color_bottoms, &image.color_tops, &a3d_params.Neighbors,
                    byValue(&a3d_params.l_step), byValue(&pCutoff), byValue(&a3d_params.lp));
    
    a3d_params.exaggeration = 1.;
    
	image.numSpots = ArgInfo[0].argIndices;
    if (image.numSpots == 0)  {  printf("GetNeighbors error:  number of image spots must be greater than zero\n");  return 1;  }
    a3d_params.colors_num = ArgInfo[6].argIndices;
    a3d_params.neighbor_lengths_num = ArgInfo[8].argIndices / (image.numSpots*a3d_params.colors_num);
    
    if ( ArgInfo[8].argIndices < a3d_params.neighbor_lengths_num*image.numSpots*a3d_params.colors_num )     {
        printf("GetNeighbors() error: top(Neighbors[]) must be a multiple of image-spots-num * colors_num\n");
        return 3;           }
    
    for (c1 = 0; c1 < ArgInfo[8].argIndices; c1++)     {
        rtrn = passed;
        if ( a3d_params.Neighbors[c1].memory == NULL )
            rtrn = newLinkedList(a3d_params.Neighbors+c1, 0, 1, 200, ccTrue);
        else if (a3d_params.Neighbors[c1].elementNum > 0)
            rtrn = deleteElements(a3d_params.Neighbors+c1, 1, a3d_params.Neighbors[c1].elementNum);
        if (rtrn != passed)
            {  printf("GetNeighbors(): out of memory\n");  return 6;  }
    }
	
    for (c1 = 0; c1 < image.numSpots; c1++)     {
        rtrn = getNeighbors(c1, pCutoff, a3d_params.Neighbors + (c1*a3d_params.colors_num*a3d_params.neighbor_lengths_num));
        if (rtrn != passed)    
            {  printf("GetNeighbors():  out of memory\n");  return 8;  }           }
    
    for (c1 = 0; c1 < ArgInfo[8].argIndices; c1++)  {
    if (defragmentLinkedList(a3d_params.Neighbors+c1) != passed)  {
        printf("GetNeighbors():  out of memory\n");
        return 8;       }}
    
	return passed;
}


// getNeighbors() finds the neighbors of a single given image spot, storing them in neighborLL

ccInt getNeighbors(ccInt firstImageSpot, double pCutoff, linkedlist *neighborLL)
{
    ccInt loopImageSpot, loopLinkerLength, loopColor, longestList, *newImageSpot, rtrn;
    double L, calcs;
    linkedlist *theLL;
    
    a3d_params.calcs = &calcs;
    
    for (loopColor = 0; loopColor < a3d_params.colors_num; loopColor++)        {
    for (loopImageSpot = image.color_bottoms[loopColor] - 1; loopImageSpot < image.color_tops[loopColor]; loopImageSpot++)      {
    if (loopImageSpot != firstImageSpot)      {

        longestList = a3d_params.neighbor_lengths_num;
        for (loopLinkerLength = 0; loopLinkerLength < a3d_params.neighbor_lengths_num; loopLinkerLength++)       {
            L = (loopLinkerLength+1) * a3d_params.l_step;
            if (GaussProb(L, firstImageSpot, loopImageSpot, 1) >= pCutoff)  {
                longestList = loopLinkerLength;
                loopLinkerLength = a3d_params.neighbor_lengths_num;
        }   }

        if (longestList != a3d_params.neighbor_lengths_num)       {
            
            theLL = neighborLL + a3d_params.neighbor_lengths_num*loopColor + longestList;
            rtrn = addElements(theLL, sizeof(ccInt), ccFalse);
            if (rtrn != passed)  return rtrn;
            
            newImageSpot = LL_int(theLL, theLL->elementNum-sizeof(ccInt)+1);
            if (newImageSpot != (ccInt *) (((char *) element(theLL, theLL->elementNum)) - sizeof(ccInt) + 1))     {
                rtrn = defragmentLinkedList(theLL);
                if (rtrn != passed)  return 8;
                newImageSpot = LL_int(theLL, theLL->elementNum-sizeof(ccInt)+1);        }
            
            *newImageSpot = loopImageSpot;            }
    }}}
    
    return passed;
}




// forEachNeighbor() performs a given operation (toDo()) for each neighbor of imageSpot, where the neighbor range is determined by 'l'

void forEachNeighbor(ccInt imageSpot, double l, ccInt neighborColor, void(*toDo)(ccInt))
{
    ccInt MaxNeighborList, cNL, *loopNeighbor, ListTop, cN; 
    linkedlist *loopLL;
    
    MaxNeighborList = floor(l / a3d_params.l_step);
    
    if (a3d_params.over_boundary)  {
    for (cN = *(image.color_bottoms + neighborColor) - 1; cN < *(image.color_tops + neighborColor); cN++)  {
        toDo(cN);
    }}
    
    else if (MaxNeighborList < a3d_params.neighbor_lengths_num)  {
    for (cNL = 0; cNL <= MaxNeighborList; cNL++)   {
        loopLL = a3d_params.Neighbors + a3d_params.neighbor_lengths_num*(a3d_params.colors_num * imageSpot + neighborColor) + cNL;
        if (loopLL->elementNum > 0)      {
            loopNeighbor = LL_int(loopLL, 1);
            ListTop = loopLL->elementNum/sizeof(ccInt);
            for (cN = 0; cN < ListTop; cN++)   {
                toDo(*loopNeighbor);
                loopNeighbor++;
    }}  }   }
    
    else    {
    for (cN = *(image.color_bottoms + neighborColor) - 1; cN < *(image.color_tops + neighborColor); cN++)  {
        if (cN != imageSpot)  toDo(cN);
    }}
}







// *********** IterateZ() ***********



// This is run by call("Iterate", ...).  It runs the overlap-potential optimizer for a given number of iterations,
// or until convergence.

int call_IterateProbs(int argc, char **argv)
{
    ccInt loop_iteration, c1, colored_numSpots, *iterations, logsize, linkedLretrn;
    double *calc_time, *C, *LogSource, *LogOutput, *FE, init_step_size;
    double line_min_param, convergence_limit;
    int conv_test;
    int highest_color, opt_method;
    gsl_multimin_function opt_fun_no_grad;
    gsl_multimin_function_fdf opt_fun_grad;
    gsl_vector *start_f, *one_step_size;
    clock_t t_start, t_end;
	arg_info *ArgInfo = (arg_info *) *(argv+argc);
	ccBool calcExact;
    
    const gsl_multimin_fdfminimizer_type *opt_alg_types[] = { gsl_multimin_fdfminimizer_conjugate_fr,
                                    gsl_multimin_fdfminimizer_conjugate_pr, gsl_multimin_fdfminimizer_vector_bfgs2, gsl_multimin_fdfminimizer_steepest_descent };
    
	const int ArgTypes[] = { double_type, double_type, double_type, double_type, double_type, double_type, double_type, double_type, double_type,
                                    int_type, int_type, double_type, int_type, int_type,
                                    string_type, string_type, string_type, double_type, double_type, double_type, double_type,
                                    string_type, double_type, double_type, string_type, double_type,
                                    double_type, int_type, double_type, double_type, double_type, double_type, double_type, double_type,
                                    double_type, double_type, double_type, double_type, double_type, int_type, double_type, int_type, int_type,
                                    bool_type, double_type };
	const int ArgIndices[] = { -1, -1, -1, -1, -1, -1, -1, -2, -2, -3, -3, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -5,
                                            1, 1, 1, 1, 1, 1, 1, 1, -6, -7, -8, 1, 1, 1, 1, 1, 1, 1, 1, -9 };
	
	if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "IterateProbs") != passed)  return 1;
	
    getArgs(argc, argv, &image.x, &image.y, &image.z, &image.dx, &image.dy, &image.dz, &a3d_params.FSumProbs, &a3d_params.f, &a3d_params.grad_f,
                    &image.color_bottoms, &image.color_tops, &contour.l, &contour.color, &contour.mask,
                    &a3d_params.p, &a3d_params.Z_1x, &a3d_params.Z_xN, &a3d_params.Z_1x_norm, &a3d_params.Z_xN_norm, &a3d_params.Z,
                    &a3d_params.dC_dZ, &a3d_params.sensitivity, &a3d_params.sensitivity_norm, &a3d_params.CSumProbs, &a3d_params.Neighbors,
                    byValue(&a3d_params.l_step), byValue(&a3d_params.lp), byValue(&a3d_params.n_skip_max), byValue(&a3d_params.p_fn),
                    byValue(&a3d_params.exaggeration), byValue(&a3d_params.K1), byValue(&init_step_size), byValue(&line_min_param),
                    &C, &LogSource, &LogOutput, &a3d_params.calcs, &calc_time, &FE, &iterations,
                    byValue(&convergence_limit), byValue(&opt_method), byValue(&a3d_params.f_or_w), byValue(&calcExact), &a3d_params.pLR);
    
	image.numSpots = ArgInfo[0].argIndices;
	contour.numSpots = ArgInfo[11].argIndices;
    
    if (contour.numSpots*image.numSpots == 0)  return passed;
    
    if (ArgInfo[7].argIndices != image.numSpots+1)       {
        printf("IterateProbs() error:  top(f) and top(grad_f) must be numSpots+1\n");
        return 2;           }
	
    logsize = ArgInfo[34].argIndices;
    if (*iterations > 0)  {
    if ( (ArgInfo[33].argIndices != ((*iterations)+1)*2) || (ArgInfo[35].argIndices != logsize * ((*iterations)+1)) )  {
        printf("IterateProbs() error:  must have iterations+1 = top(Cs) and top(LogSource)*(iterations+1) = top(LogOutput)\n");
        return 3;           }}
    
    if (opt_method > 4)  {  printf("IterateProbs() error:  opt_method must be from 0-4\n");  return 3;  }
    
    highest_color = 0;
    for (c1 = 0; c1 < contour.numSpots; c1++)  {
    if (*(contour.color + c1) > highest_color)  {
        highest_color = *(contour.color + c1);      }}
    
    if (ArgInfo[9].argIndices < highest_color+1)  {
        printf("IterateProbs() error: lowest_color[], highest_color[] need to be of the same length > max(contour.color, image.color)\n");
        return 3;           }
    
    a3d_params.colors_num = ArgInfo[9].argIndices;
    a3d_params.neighbor_lengths_num = ArgInfo[24].argIndices / (image.numSpots*a3d_params.colors_num);
    
    if ( ArgInfo[24].argIndices < a3d_params.neighbor_lengths_num*image.numSpots*a3d_params.colors_num )  {
        printf("IterateProbs() error: top(Neighbors[]) must be a multiple of image-spots-num * colors_num\n");
        return 4;           }

    for (c1 = 0; c1 < contour.numSpots; c1++)  {
    if (AreSpots(c1))  {
        colored_numSpots = (image.color_tops[contour.color[c1]] + 1) - image.color_bottoms[contour.color[c1]];
        
        if (colored_numSpots > 0)  {
            if ( (a3d_params.p[c1].memory == 0)
                     || (a3d_params.Z_1x[c1].memory == 0) || (a3d_params.Z_xN[c1].memory == 0)
                     || (a3d_params.p[c1].elementNum != colored_numSpots*sizeof(double))
                     || (a3d_params.Z_1x[c1].elementNum != colored_numSpots*sizeof(double))
                     || (a3d_params.Z_xN[c1].elementNum != colored_numSpots*sizeof(double))  )  {
                printf("IterateProbs(): ODotProbs[%i], Z_1x/xN[%i] need to be initialized to the same size (%i)\n", 
                            c1, c1, colored_numSpots);
                return 5;           }

            linkedLretrn = defragmentLinkedList(a3d_params.p+c1);
            if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(a3d_params.Z_1x+c1);
            if (linkedLretrn == passed)  linkedLretrn = defragmentLinkedList(a3d_params.Z_xN+c1);
            if (linkedLretrn != passed)     {
                printf("IterateProbs(): out of memory\n");
                return 6;
    }}   }  }
    
    for (c1 = 0; c1 < ArgInfo[24].argIndices; c1++)     {
        linkedLretrn = defragmentLinkedList(a3d_params.Neighbors+c1);
        if (linkedLretrn != passed)     {
            printf("IterateProbs(): out of memory\n");
            return 6;
    }   }
	
    t_start = clock();
    
    if (calcExact)  {
        
        ccInt i, alpha;
        
        a3d_params.alphas = (ccInt *) malloc(contour.numSpots * sizeof(ccInt));
        a3d_params.alpha_mask = (ccBool *) malloc(contour.numSpots * sizeof(ccBool));
        
        a3d_params.f[image.numSpots] = exp(-(*(a3d_params.f + image.numSpots)));
        InitW();
        a3d_params.f[image.numSpots] = -log(*(a3d_params.f + image.numSpots));
        
        for (i = 0; i < contour.numSpots; i++)  {
            
            ccInt one_color = contour.color[i];
            ccInt alpha_base = image.color_bottoms[one_color] - 1;
            
            for (alpha = alpha_base; alpha < image.color_tops[one_color]; alpha++)    {
                *(LL_Double(a3d_params.p + i, 1) + alpha - alpha_base) = 0.;
        }   }
        
        for (alpha = 0; alpha < image.numSpots; alpha++)
            a3d_params.alpha_mask[alpha] = ccFalse;
        
        *(a3d_params.Z) = 0.;
        GetAllChains(0, -1, 1., (*iterations > 0));
        
        for (i = 0; i < contour.numSpots; i++)  {
            
            ccInt one_color = contour.color[i];
            ccInt alpha_base = image.color_bottoms[one_color] - 1;
            
            for (alpha = alpha_base; alpha < image.color_tops[one_color]; alpha++)    {
                LL_Double(a3d_params.p + i, 1)[alpha - alpha_base] /= *(a3d_params.Z);
        }   }
        
        free(a3d_params.alphas);
        free(a3d_params.alpha_mask);
    }
    
    else if (*iterations == 0)  {

        load_f(NULL);
        
        C[0] = IterateProbs(2);
        
        save_grad_f(NULL);
        for (c1 = 0; c1 < image.numSpots+1; c1++)  a3d_params.f[c1] = -log(a3d_params.f[c1]);
        
        C[1] = 0;
        for (c1 = 0; c1 < image.numSpots+1; c1++)  C[1] += a3d_params.grad_f[c1] * a3d_params.grad_f[c1];
        C[1] = sqrt(C[1]);          }
    
    else  {
        
        opt_fun_no_grad.n = opt_fun_grad.n = image.numSpots+1;
        opt_fun_no_grad.params = opt_fun_grad.params = 0;
        opt_fun_no_grad.f = opt_fun_grad.f = &GetZ;
        opt_fun_grad.df = &GetGradZ;
        opt_fun_grad.fdf = &GetZAndGradZ;
        
        start_f = gsl_vector_alloc(image.numSpots+1);
        for (c1 = 0; c1 < image.numSpots+1; c1++)
            gsl_vector_set(start_f, c1, a3d_params.f[c1]);
        
        if (opt_method == 0)  {
            opt_struct_no_grad = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, image.numSpots+1);
            if (opt_struct_no_grad == 0)  return 10;
            
            one_step_size = gsl_vector_alloc(image.numSpots+1);
            for (c1 = 0; c1 < image.numSpots+1; c1++)
                gsl_vector_set(one_step_size, c1, init_step_size);
            gsl_multimin_fminimizer_set(opt_struct_no_grad, &opt_fun_no_grad, start_f, one_step_size);         }
        else  {
            opt_struct_grad = gsl_multimin_fdfminimizer_alloc(opt_alg_types[opt_method - 1], image.numSpots+1);
            if (opt_struct_grad == 0)  return 10;
            gsl_multimin_fdfminimizer_set(opt_struct_grad, &opt_fun_grad, start_f, init_step_size, line_min_param);        }
        
        conv_test = GetC(C, 0, convergence_limit, opt_method);
        
        for (loop_iteration = 0; loop_iteration < *iterations; loop_iteration++)     {
            
            if (opt_method == 0)  {
                gsl_multimin_fminimizer_iterate(opt_struct_no_grad);        }
            else  {
                gsl_multimin_fdfminimizer_iterate(opt_struct_grad);         }
            
            conv_test = GetC(C, loop_iteration+1, convergence_limit, opt_method);
            
            if (logsize > 0)  {
            for (c1 = 0; c1 < logsize; c1++)  {
                LogOutput[logsize*(loop_iteration + 1) + c1] = LogSource[c1];      }}
            
            if (conv_test == GSL_SUCCESS)  {  *iterations = loop_iteration + 1;  break;  }
    }   }
    
    t_end = clock();
    *calc_time += ((double) t_end - t_start) / CLOCKS_PER_SEC;
    
    if ((!calcExact) && (*iterations > 0))        {
        
        if (opt_method == 0)        {
            gsl_multimin_fminimizer_free(opt_struct_no_grad);
            gsl_vector_free(one_step_size);         }
        else
            gsl_multimin_fdfminimizer_free(opt_struct_grad);
        
        gsl_vector_free(start_f);           }
    
	return passed;
}


// GetC() stores the state of the minimization (f (best guess), C), and returns true iff converged

int GetC(double *C, ccInt C_offset, double convergence_limit, int opt_method)
{
    ccInt c1;
    double conv_param, alg0_size;
    int conv_test;
    gsl_vector *one_f, *one_grad;
    
    if (opt_method == 0)     {
        one_f = gsl_multimin_fminimizer_x(opt_struct_no_grad);
        C[2*C_offset] = gsl_multimin_fminimizer_minimum(opt_struct_no_grad);
        conv_param = alg0_size = gsl_multimin_fminimizer_size(opt_struct_no_grad);
        conv_test = gsl_multimin_test_size(alg0_size, convergence_limit);               }
    else        {
        one_f = gsl_multimin_fdfminimizer_x(opt_struct_grad);
        C[2*C_offset] = gsl_multimin_fdfminimizer_minimum(opt_struct_grad);
        one_grad = gsl_multimin_fdfminimizer_gradient(opt_struct_grad);
        conv_param = gsl_blas_dnrm2(one_grad);
        conv_test = gsl_multimin_test_gradient(one_grad, convergence_limit);
        
        for (c1 = 0; c1 < image.numSpots+1; c1++)
            a3d_params.grad_f[c1] = gsl_vector_get(one_grad, c1);           }
    
    for (c1 = 0; c1 < image.numSpots+1; c1++)
        a3d_params.f[c1] = gsl_vector_get(one_f, c1);
    
    C[2*C_offset + 1] = conv_param;
    
    return conv_test;
}


// Next three routines:  return C and/or grad-C (used by various optimization routines)

double GetZ(const gsl_vector *OF, void *dummy)
{
    load_f(OF);
    
    return IterateProbs(1);
}


void GetGradZ(const gsl_vector *OF, void *dummy, gsl_vector *grad_OF)
{
    load_f(OF);
    
    IterateProbs(2);
    
    save_grad_f(grad_OF);
}


void GetZAndGradZ(const gsl_vector *OF, void *dummy, double *C, gsl_vector *grad_OF)
{
    load_f(OF);
    
    *C = IterateProbs(2);
    
    save_grad_f(grad_OF);
}


// load_f() converts (-log f) to f

void load_f(const gsl_vector *OF)
{
    ccInt c1;
    double one_f;
    
    for (c1 = 0; c1 < image.numSpots+1; c1++)  {
        if (OF == NULL)  one_f = a3d_params.f[c1];
        else  one_f = gsl_vector_get(OF, c1);
        
        a3d_params.f[c1] = exp(-one_f);            }
}


// save_grad_f() converts dC/df to dC/d(-log f)

void save_grad_f(const gsl_vector *grad_OF)
{
    ccInt c1;
    double one_grad_f;
    
    for (c1 = 0; c1 < image.numSpots+1; c1++)  {
        one_grad_f = -a3d_params.grad_f[c1] * a3d_params.f[c1];
        
        if (grad_OF == NULL)  a3d_params.grad_f[c1] = one_grad_f;
        else  gsl_vector_set((gsl_vector *) grad_OF, c1, one_grad_f);           }
}



// Forward-propagates the Z matrices, and back-propagates the sensitivities.  Three modes:
// 0: do nothing
// 1: prop. Zs to calculate p
// 2: (1) plus: sensitivities

double IterateProbs(int toDo)
{
    ccInt alpha, alpha_base, one_color;
    int i;
    double olf, *first_p, C, p_tot;

if (a3d_params.f_or_w == 1) * (a3d_params.f + image.numSpots) = 1;
else  for (alpha = 0; alpha < image.numSpots; alpha++)  *(a3d_params.f + alpha) = 1;
    
    InitW();
    
    if (toDo >= 1)     {
        
        PropArray(a3d_params.Z_1x, a3d_params.Z_1x_norm, &Z_prop, 1, pa_FillZ);
        PropArray(a3d_params.Z_xN, a3d_params.Z_xN_norm, &Z_prop, -1, pa_FillZ);
        FillArray(fa_CalcP);        // p0 = Z_1x * Z_xN
        PropArray(a3d_params.Z_xN, a3d_params.Z_xN_norm, &Z_bridge, -1, pa_BridgeZ);
        FillArray(fa_DivZ);         // p = p0 / Z
    }
//PropArray(a3d_params.Z_1x, a3d_params.Z_1x_norm, &do_pLR, 1, 100);
    
    
            // Sum probs along the contour
    
    p_tot = 0;
    for (i = 0; i < contour.numSpots; i++)  {
        
        a3d_params.CSumProbs[i] = 0;
        
        if (AreSpots(i))     {
            first_p = LL_Double(a3d_params.p + i, 1);
            one_color = contour.color[i];
            
            for (alpha = 0; alpha < image.color_tops[one_color] - image.color_bottoms[one_color] + 1; alpha++)     {
                a3d_params.CSumProbs[i] += first_p[alpha];
                p_tot += first_p[alpha];
    }   }   }
    
    
            // Sum probs over the image
    
    for (one_color = 0; one_color < a3d_params.colors_num; one_color++)  {
        alpha_base = image.color_bottoms[one_color] - 1;
        for (alpha = alpha_base; alpha < image.color_tops[one_color]; alpha++)  {
            
            olf = 0;
            for (i = 0; i < contour.numSpots; i++)  {
            if (contour.color[i] == one_color)  {
                olf += LL_Double(a3d_params.p + i, 1)[alpha - alpha_base];
            }}
            a3d_params.FSumProbs[alpha] = olf;
    }   }
    
    a3d_params.N_difference = p_tot - (1. - a3d_params.p_fn) * contour.numSpots;
    C = 0.5 * a3d_params.K1 * a3d_params.N_difference * a3d_params.N_difference;
    for (alpha = 0; alpha < image.numSpots; alpha++)  {
    if (a3d_params.FSumProbs[alpha] > 1)  {
        C += 0.5 * (a3d_params.FSumProbs[alpha] - 1.) * (a3d_params.FSumProbs[alpha] - 1.);
    }}
//    C /= contour.numSpots;
    
    
    if (toDo <= 1)  return C;
    
    
        // compute the gradient dC/df and dC/dw
    
    for (alpha = 0; alpha < image.numSpots; alpha++)
        a3d_params.grad_f[alpha] = 0;
    a3d_params.grad_f[image.numSpots] = 0;
    
    FillArray(fa_dC_dZ);
    PropArray(a3d_params.Z_xN, a3d_params.Z_xN_norm, &Z_bridge, -1, pa_Bridge_dC_dZ);
    
    FillArray(fa_SetS1x);
    PropArray(a3d_params.Z_xN, a3d_params.Z_xN_norm, &Z_bridge, -1, pa_BridgeS);
    PropArray(a3d_params.sensitivity, a3d_params.sensitivity_norm, &Z_prop, -1, pa_PropS);
    FillArray(fa_AddS1x);
    if (a3d_params.f_or_w == 2)  PropArray(a3d_params.Z_1x, a3d_params.Z_1x_norm, &Z_prop, 1, pa_GetGradW);/**/
    
    FillArray(fa_SetSxN);
    PropArray(a3d_params.Z_1x, a3d_params.Z_1x_norm, &Z_bridge, 1, pa_BridgeS);
    PropArray(a3d_params.sensitivity, a3d_params.sensitivity_norm, &Z_prop, 1, pa_PropS);
    FillArray(fa_AddSxN);
    if (a3d_params.f_or_w == 2)  PropArray(a3d_params.Z_xN, a3d_params.Z_xN_norm, &Z_prop, -1, pa_GetGradW);/**/
    
if (a3d_params.f_or_w == 1)  *(a3d_params.grad_f + image.numSpots) = 0;
else  for (alpha = 0; alpha < image.numSpots; alpha++)  *(a3d_params.grad_f + alpha) = 0;
    return C;
}



// PropArray() fills the elements of arrays like Z_1x, etc. whose fill rule is like a layered neural network.

void PropArray(linkedlist *prop_list, double *prop_list_norm, void(*NeighborFunction)(ccInt), int direction, int mode)
{
    int i, n_skip, one_color, skip_counter, boundary_overshoot;
    ccInt alpha, other_boundary, alpha_base, alpha_top;
    double *one_layer, *other_norm_list, *end_layer, *one_sense, norm_offset;
    double new_norm;
    linkedlist *other_list;
    
    if (direction == 1)  {  a3d_params.boundary = 0;  other_boundary = contour.numSpots - 1;  }
    else  {  a3d_params.boundary = contour.numSpots - 1;  other_boundary = 0;  }
    a3d_params.direction = direction;
    
    a3d_params.prop_u_or_d = prop_list;
    a3d_params.prop_u_or_d_norm = prop_list_norm;
    a3d_params.prop_mode = mode;
    
    if ((mode == pa_BridgeS) && (direction < 0) && (a3d_params.f_or_w == /*1*/2))  a3d_params.pa_SetGradWSource = ccTrue;
    else  a3d_params.pa_SetGradWSource = ccFalse;
    
    boundary_overshoot = 0;
    if ((mode == pa_BridgeZ) || (mode == pa_Bridge_dC_dZ) || (mode == pa_BridgeS) || (mode == pa_GetGradW))
          {  if (mode != pa_GetGradW)  boundary_overshoot = 1;  }
    
    if ((prop_list == a3d_params.Z_1x) || ((a3d_params.prop_mode == pa_PropS) && (direction == 1)))
          {  other_list = a3d_params.Z_xN;  other_norm_list = a3d_params.Z_xN_norm;  }
    else  {  other_list = a3d_params.Z_1x;  other_norm_list = a3d_params.Z_1x_norm;  }
    
    for (i = a3d_params.boundary; (i+boundary_overshoot >= 0) && (i-boundary_overshoot < (ccInt) contour.numSpots); i += direction)  {
    if (AreSpots(i))     {
        
        if ((i < 0) || (i >= contour.numSpots))  {
            a3d_params.over_boundary = ccTrue;
            a3d_params.other_norm = 0;
            alpha_base = 0;
            alpha_top = 1;     }
        
        else  {
            
            a3d_params.over_boundary = ccFalse;
            
            a3d_params.other_norm = *(other_norm_list + i);
            
            one_layer = LL_Double(prop_list + i, 1);
            a3d_params.other_layer = LL_Double(other_list + i, 1);
            
            one_sense = LL_Double(a3d_params.sensitivity + i, 1);
            one_color = contour.color[i];
            
            alpha_base = image.color_bottoms[one_color] - 1;
            alpha_top = image.color_tops[one_color];            }
        
        
        for (alpha = alpha_base; alpha < alpha_top; alpha++)     {
            
            if (mode == pa_FillZ)  one_layer[alpha - alpha_base] = 0.;
            a3d_params.alpha_0 = alpha;
            
            for (n_skip = 0; n_skip <= a3d_params.n_skip_max; n_skip++)        {
                
                if (abs(i - a3d_params.boundary) < n_skip)  break;
                
                a3d_params.cSpot = i - (n_skip + 1) * direction;
                a3d_params.cSpot_0 = i;
                
                if (AreSpots(a3d_params.cSpot))     {
                    
                    a3d_params.w_factor = n_skip*a3d_params.w;
                    
                    if (abs(i - a3d_params.boundary) == n_skip)  {
                        a3d_params.p_value = 0.;
                        for (one_color = 0; one_color < a3d_params.colors_num; one_color++)
                            forEachNeighbor(alpha, (contour.l[contour.numSpots - 1]*(abs(i-a3d_params.boundary)+1.))/contour.numSpots, one_color, &CountNeighbors);
a3d_params.p_value = 1.;
                        if (a3d_params.p_value == 0.)  a3d_params.p_value = 1.;
                        new_norm = 0.;       }
                    
                    else  {
                        a3d_params.p_value = 0;
                        a3d_params.l = fabs( contour.l[a3d_params.cSpot] - contour.l[i] );
                        a3d_params.alpha_base = image.color_bottoms[contour.color[a3d_params.cSpot]] - 1;
                        a3d_params.alpha_0_base = alpha_base;
                        
                        new_norm = prop_list_norm[a3d_params.cSpot];        }
                    
                    if (NeighborFunction == &Z_bridge)  {
                        
                        if (a3d_params.over_boundary)  a3d_params.multiplier = 1.;
                        else if (mode == pa_BridgeS)  a3d_params.multiplier = sqrt(a3d_params.f[alpha]);
                        else  a3d_params.multiplier = sqrt(a3d_params.f[alpha]) * (*(LL_Double(a3d_params.Z_1x + i, 1) + (alpha - alpha_base)));
                        if (abs(i - a3d_params.boundary) == n_skip)  a3d_params.p_value = a3d_params.multiplier;
                        
                        if (mode == pa_Bridge_dC_dZ)  {
                            a3d_params.dC_dZ_factor = 0.;
                            for (skip_counter = 1; skip_counter <= n_skip; skip_counter++)   {
                                a3d_params.dC_dZ_factor += 0.5 * a3d_params.dC_dZ[i - skip_counter*direction] * exp( a3d_params.w_factor + a3d_params.other_norm
                                        + new_norm - a3d_params.Z_1x_norm[ + i - skip_counter*direction] - a3d_params.Z_xN_norm[i - skip_counter*direction] );
                    }   }   }
                    
                    if (abs(i - a3d_params.boundary) != n_skip)
                        forEachNeighbor(alpha, a3d_params.l, contour.color[a3d_params.cSpot], NeighborFunction);
                    
                    if (mode == pa_FillZ)  {
                        double next_norm = 0., cSpot_norm = 0.;
                        if (i != a3d_params.boundary)  next_norm = prop_list_norm[i - direction];
                        if (abs(i - a3d_params.boundary) != n_skip)  cSpot_norm = prop_list_norm[a3d_params.cSpot];
                        one_layer[alpha - alpha_base] += a3d_params.p_value * sqrt(a3d_params.f[alpha])
                                                * exp(a3d_params.w_factor + new_norm - next_norm);        }
                    
                    else if (mode == pa_GetGradW)  {
                        a3d_params.grad_f[image.numSpots] += a3d_params.w_norm * n_skip * a3d_params.p_value * sqrt(a3d_params.f[alpha])
                                * one_sense[alpha - alpha_base] * exp(a3d_params.w_factor - a3d_params.w + new_norm + a3d_params.sensitivity_norm[i]);    }
                    
                    else if (mode == pa_BridgeZ)     {
                        for (skip_counter = 1; skip_counter <= n_skip; skip_counter++)  {
                            a3d_params.Z[i - skip_counter*direction] += a3d_params.p_value * exp( a3d_params.w_factor + a3d_params.other_norm + new_norm
                                    - a3d_params.Z_1x_norm[i - skip_counter*direction] - a3d_params.Z_xN_norm[i - skip_counter*direction] );
                    }   }
                    
                    else if ((mode == pa_Bridge_dC_dZ) && (!a3d_params.over_boundary))  {
                        a3d_params.grad_f[alpha] += a3d_params.p_value * a3d_params.dC_dZ_factor / a3d_params.f[alpha];          }
                    
                    else if (mode == pa_BridgeS)  {
                        double to_add;
                        
                        for (skip_counter = 1; skip_counter <= n_skip; skip_counter++)      {
                            
                            if (!a3d_params.over_boundary)  {
                                end_layer = LL_Double(a3d_params.sensitivity + i, 1);
                                
                                to_add = a3d_params.dC_dZ[i - skip_counter*direction] * a3d_params.p_value
                                            * exp( a3d_params.w_factor + new_norm - a3d_params.Z_1x_norm[i - skip_counter*direction]
                                                    - a3d_params.Z_xN_norm[i - skip_counter*direction] - a3d_params.sensitivity_norm[i] );
                                
                                end_layer[alpha - alpha_base] += to_add;     // next line:  -(1/2) Z^i (dC/dZ^i) / f_i  (or with Z_i)
                                a3d_params.grad_f[alpha] -= 0.5 * a3d_params.other_layer[alpha - alpha_base] * to_add
                                                 * exp( a3d_params.other_norm + a3d_params.sensitivity_norm[i] ) / a3d_params.f[alpha];   }
                            
                            if (a3d_params.pa_SetGradWSource)  {
                                to_add = a3d_params.w_norm * n_skip * a3d_params.dC_dZ[i - skip_counter*direction]
                                        * a3d_params.p_value * exp( a3d_params.w_factor - a3d_params.w + new_norm + a3d_params.other_norm
                                                - a3d_params.Z_1x_norm[i - skip_counter*direction] - a3d_params.Z_xN_norm[i - skip_counter*direction] );
                                if (!a3d_params.over_boundary)  to_add *= LL_Double(a3d_params.Z_1x + i, 1)[alpha - alpha_base];
                                a3d_params.grad_f[image.numSpots] += to_add;
                    }   }   }
                    
                    else if ((mode == pa_PropS) && (abs(i - a3d_params.boundary) != n_skip))  {
                        double to_add = a3d_params.p_value * sqrt(a3d_params.f[alpha])
                                            * exp(a3d_params.w_factor + prop_list_norm[a3d_params.cSpot] - prop_list_norm[i]);
                        one_layer[alpha - alpha_base] += to_add;                }
                }
            }
        }
        
            // renormalize the Z/s arrays
        
        if (mode == pa_FillZ)      {
            if (i == a3d_params.boundary)  norm_offset = 0;
            else  norm_offset = prop_list_norm[i - direction];
            RenormZ(one_layer, alpha_top - alpha_base, prop_list_norm + i, norm_offset);            }
        
        else if (mode == pa_PropS)      {
            norm_offset = prop_list_norm[i];
            RenormZ(one_layer, alpha_top - alpha_base, prop_list_norm + i, norm_offset);            }
    }}
}



// FillArray() iterates -- not propagates -- over each element of some array.

void FillArray(int mode)
{
    ccInt i, alpha, alpha_base, one_color;
    double *first_p, *one_Z_1x, *one_Z_xN, *one_sense, one_p, dC_dp;
    
    for (i = 0; i < contour.numSpots; i++)  {
    if (contour.mask[i] != 0)  {
    if (AreSpots(i))  {
        
        first_p = LL_Double(a3d_params.p + i, 1);
        one_Z_1x = LL_Double(a3d_params.Z_1x + i, 1);
        one_Z_xN = LL_Double(a3d_params.Z_xN + i, 1);
        one_sense = LL_Double(a3d_params.sensitivity + i, 1);
        one_color = contour.color[i];
        
        alpha_base = image.color_bottoms[one_color] - 1;
        
        if (mode == fa_CalcP)  a3d_params.Z[i] = 0;
        else if (mode == fa_dC_dZ)  a3d_params.dC_dZ[i] = 0;
        
        for (alpha = alpha_base; alpha < image.color_tops[one_color]; alpha++)  {
            
            dC_dp = a3d_params.K1 * a3d_params.N_difference;
            if (a3d_params.FSumProbs[alpha] > 1.)  dC_dp += a3d_params.FSumProbs[alpha] - 1.;
            
            if (mode == fa_CalcP)  {
                one_p = one_Z_1x[alpha - alpha_base] * one_Z_xN[alpha - alpha_base];
                first_p[alpha - alpha_base] = one_p;
                a3d_params.Z[i] += one_p;            }
            
            else if ((mode == fa_DivZ) && (*(a3d_params.Z + i) > 0.))  {
                first_p[alpha - alpha_base] /= a3d_params.Z[i];     }
            
            else if (mode == fa_dC_dZ)  {
                a3d_params.dC_dZ[i] -= dC_dp * first_p[alpha - alpha_base] / a3d_params.Z[i];     }
            
            else if (mode == fa_SetS1x)  {
                double to_add = a3d_params.dC_dZ[i] * one_Z_xN[alpha - alpha_base];
                if (one_Z_1x[alpha - alpha_base] > 0)
                    to_add += dC_dp * one_Z_xN[alpha - alpha_base] / a3d_params.Z[i];
                
                one_sense[alpha - alpha_base] = to_add;
                a3d_params.grad_f[alpha] -= 0.5 * one_Z_1x[alpha - alpha_base] * to_add / a3d_params.f[alpha];      }
            
            else if (mode == fa_AddS1x)  {
                a3d_params.grad_f[alpha] += one_Z_1x[alpha - alpha_base] * one_sense[alpha - alpha_base]
                                            * exp(a3d_params.Z_1x_norm[i] + a3d_params.sensitivity_norm[i]) / a3d_params.f[alpha];      }
            
            else if (mode == fa_SetSxN)  {
                double to_add = a3d_params.dC_dZ[i] * one_Z_1x[alpha - alpha_base];
                if (one_Z_xN[alpha - alpha_base] > 0)
                    to_add += dC_dp * one_Z_1x[alpha - alpha_base] / a3d_params.Z[i];
                
                one_sense[alpha - alpha_base] = to_add;
                a3d_params.grad_f[alpha] -= 0.5 * one_Z_xN[alpha - alpha_base] * to_add / a3d_params.f[alpha];      }
            
            else if (mode == fa_AddSxN)  {
                a3d_params.grad_f[alpha] += one_Z_xN[alpha - alpha_base] * one_sense[alpha - alpha_base]
                                            * exp(a3d_params.Z_xN_norm[i] + a3d_params.sensitivity_norm[i]) / a3d_params.f[alpha];  }
        }
    }}}
    
        // set the initial normalizations of the Z/s arrays
    
    if (mode == fa_SetS1x)  {
    for (i = 0; i < contour.numSpots; i++)  {
        a3d_params.sensitivity_norm[i] = -a3d_params.Z_1x_norm[i];     }}
    else if (mode == fa_SetSxN)  {
    for (i = 0; i < contour.numSpots; i++)  {
        a3d_params.sensitivity_norm[i] = -a3d_params.Z_xN_norm[i];     }}
}



// RenormZ() rescales one column of a Z/s array along with its normalization factor, stored in a separate array.

void RenormZ(double *first_z, ccInt zs_num, double *Norm, double norm_offset)
{
    ccInt c1;
    double tot, norm;
    
    tot = 0;
    for (c1 = 0; c1 < zs_num; c1++)  tot += fabs(first_z[c1]);
    
    norm = tot;
    if (norm == 0)  { *Norm = norm_offset;  return;  }
    
    for (c1 = 0; c1 < zs_num; c1++)  first_z[c1] /= norm;
    *Norm = log(norm) + norm_offset;
}



// Z_prop() -- the basic routine propagating terms in the Z/s arrays

void Z_prop(ccInt alpha)
{
    a3d_params.p_value += LL_Double(a3d_params.prop_u_or_d + a3d_params.cSpot, 1)[alpha - a3d_params.alpha_base]
                    * GaussProb(a3d_params.l, a3d_params.alpha_0, alpha, 0) * sqrt(a3d_params.f[alpha]);
}


// Z_bridge() introduces terms coming from false negatives

void Z_bridge(ccInt alpha)
{
    if (abs(a3d_params.cSpot_0 - a3d_params.cSpot) > 1)  {
        
        double to_add = LL_Double(a3d_params.prop_u_or_d + a3d_params.cSpot, 1)[alpha - a3d_params.alpha_base]
                        * sqrt(a3d_params.f[alpha]) * a3d_params.multiplier;
        if (!a3d_params.over_boundary)  to_add *= GaussProb(a3d_params.l, a3d_params.alpha_0, alpha, 0);
        
        if (a3d_params.prop_mode == pa_Bridge_dC_dZ)  {
            a3d_params.grad_f[alpha] += to_add * a3d_params.dC_dZ_factor / a3d_params.f[alpha];    }
        a3d_params.p_value += to_add;               }
}


void CountNeighbors(ccInt alpha)
{
    a3d_params.p_value += 1.;
}

void do_pLR(ccInt alpha)
{
    double L = a3d_params.l;
    double R = sqrt(SqDotDistance(alpha, a3d_params.alpha_0, 1., 1., 1.));
    ccInt x_bin = L / 2;
    ccInt y_bin = R / 2;
    
    if ((x_bin >= 0) && (y_bin >= 0) && (x_bin < 100) && (y_bin < 100))
        *(a3d_params.pLR + y_bin*100 + x_bin) += GaussProb(a3d_params.l, a3d_params.alpha_0, alpha, 0)*exp(a3d_params.w*(a3d_params.cSpot_0-a3d_params.cSpot-1))
                * (*(LL_Double(a3d_params.Z_1x + a3d_params.cSpot, 1) + (alpha - a3d_params.alpha_base)))
                * (*(LL_Double(a3d_params.Z_xN + a3d_params.cSpot_0, 1) + (a3d_params.alpha_0 - a3d_params.alpha_0_base)))
                * exp((*(a3d_params.Z_1x_norm + a3d_params.cSpot)) - (*(a3d_params.Z_1x_norm + a3d_params.cSpot_0)))
                / (*(a3d_params.Z + a3d_params.cSpot_0));
}


// InitW() computes the w-factor based on a 'w_norm' that compares sensibly with the mean terms connecting pairs of spots.

void InitW()
{
    double L_avg = contour.l[contour.numSpots - 1] / (contour.numSpots * (1. - a3d_params.p_fn));
    a3d_params.w_norm = GaussProb(L_avg, 0, 0, 2) * GaussProb(L_avg, 0, 0, 2) / GaussProb(L_avg*(2. - a3d_params.p_fn), 0, 0, 2);
    a3d_params.w = log( a3d_params.f[image.numSpots] * a3d_params.w_norm );
}


void GetAllChains(ccInt spotNum, ccInt lastSpot, double p, ccBool noOverlapsAtAll)
{
    if (spotNum == contour.numSpots)  {
        
        ccInt j;
        
        *(a3d_params.Z) += p;
        
        for (j = 0; j < contour.numSpots; j++)  {
            
            ccInt j_color = contour.color[j];
            ccInt j_alpha_base = (ccInt) image.color_bottoms[j_color] - 1;
            ccInt j_alpha = a3d_params.alphas[j];
            
            if (j_alpha >= j_alpha_base)  {
                LL_Double(a3d_params.p + j, 1)[j_alpha - j_alpha_base] += p;
        }   }
    }
    
    else  {
        
        ccInt one_color = contour.color[spotNum];
        ccInt alpha_base = image.color_bottoms[one_color] - 1, alpha;
        
        for (alpha = alpha_base-1; alpha < (ccInt) image.color_tops[one_color]; alpha++)  {
            
            ccBool if_masked = ccFalse;
            
            if (noOverlapsAtAll)  {
            if (alpha > alpha_base-1)  {
                if_masked = a3d_params.alpha_mask[alpha];        }}
            else  {
            if (spotNum > 0)  {
                if_masked = ((alpha == a3d_params.alphas[spotNum - 1]) && (alpha >= 0));       }}
            
            if (!if_masked)  {
                
                a3d_params.alphas[spotNum] = alpha;
                if (alpha > alpha_base-1)  a3d_params.alpha_mask[alpha] = ccTrue;
                
                if (alpha == alpha_base - 1)  {
                if (spotNum - lastSpot < a3d_params.n_skip_max)  {
                    GetAllChains(spotNum+1, lastSpot, p * exp(a3d_params.w) / a3d_params.w_norm, noOverlapsAtAll);
                }}
                
                else  {
                    
                    double extra_p = 1.;
                    if (lastSpot != -1)  {
                        double l = fabs( contour.l[spotNum] - contour.l[lastSpot] );
                        extra_p = GaussProb(l, alpha, a3d_params.alphas[lastSpot], 0) / a3d_params.w_norm;     }
                    
                    GetAllChains(spotNum+1, spotNum, p*extra_p, noOverlapsAtAll);        }
                
                if (alpha > alpha_base-1)  a3d_params.alpha_mask[alpha] = ccFalse;
    }   }   }
}





// *********** Misc ***********



int call_GaussianChain(int argc, char **argv)
{
	arg_info *ArgInfo = (arg_info *) *(argv+argc);
    ccInt FD1, FD2;
    double L, dummy_calcs, *result;

	const int ArgTypes[] = { double_type, double_type, double_type, double_type, double_type, double_type,
                                     double_type, double_type, int_type, int_type, double_type };
	const int ArgIndices[] = { -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1 };
	
	if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "GaussianChain") != passed)  return 1;
	
    getArgs(argc, argv, &(image.x), &(image.y), &(image.z), &(image.dx), &(image.dy), &(image.dz),
                byValue(&L), byValue(&(a3d_params.lp)), byValue(&FD1), byValue(&FD2), &result);
	
    if ((FD1 == 0) || (FD2 == 0) || (FD1 > image.numSpots) || (FD2 > image.numSpots))       {
        printf("GaussianChain():  image spot 1 or 2 out of range\n");
        return 2;               }
    
    a3d_params.exaggeration = 1.;
    a3d_params.calcs = &dummy_calcs;
    
    *result = GaussProb(L, FD1-1, FD2-1, 0);
    
    return passed;
}



double GaussProb(double l, ccInt imageSpot1, ccInt imageSpot2, char GN_mode)
{
    double alpha, ax, ay, az, ans;
    const double pi3 = pi*pi*pi;
    
    *(a3d_params.calcs) += 1;
    
    alpha = 3/(2*EffectiveL2(l));
    
    ax = 1./(1./alpha + 2*image.dx[imageSpot1]*image.dx[imageSpot1] + 2*image.dx[imageSpot2]*image.dx[imageSpot2]);
    ay = 1./(1./alpha + 2*image.dy[imageSpot1]*image.dy[imageSpot1] + 2*image.dy[imageSpot2]*image.dy[imageSpot2]);
    az = 1./(1./alpha + 2*image.dz[imageSpot1]*image.dz[imageSpot1] + 2*image.dz[imageSpot2]*image.dz[imageSpot2]);
    
    if (GN_mode == 0)
        ans = sqrt(ax*ay*az/pi3)*exp(-SqDotDistance(imageSpot1, imageSpot2, ax, ay, az));
    else if (GN_mode == 1)
        ans = exp(-SqDotDistance(imageSpot1, imageSpot2, ax, ay, az));
    else
        ans = sqrt(ax*ay*az/pi3);
    
    if (a3d_params.exaggeration != 1.)
        ans = pow(ans, a3d_params.exaggeration);
//    if (a3d_params.exaggeration != 1.)  ans = pow(ans, a3d_params.exaggeration);
    return ans;
}



double EffectiveL2(double l)
{
    if (l > 2*a3d_params.lp)  return 2*l*a3d_params.lp;
    else if (l*l > 0)  return l*l;
    else  return 1.;            // 1 bp
}


double SqDotDistance(ccInt spot1, ccInt spot2, double wx, double wy, double wz)
{
	double dx, dy, dz;
    
    dx = image.x[spot1] - image.x[spot2];
    dy = image.y[spot1] - image.y[spot2];
    dz = image.z[spot1] - image.z[spot2];
    
    return wx*dx*dx + wy*dy*dy + wz*dz*dz;
}

double AddLog(double x1, double x2, double sgn)
{
    if (x1 > x2)  return x1 + log(1 + sgn*exp(x2 - x1));
    else  return x2 + log(1 + sgn*exp(x1 - x2));
}

double max(double d1, double d2)
{
    if (d1 > d2)  return d1;
    else  return d2;
}



int call_InitGaussRand(int argc, char **argv)
{
	arg_info *ArgInfo = (arg_info *) argv[argc];
    clock_t rand_seed;
    
	const int ArgTypes[] = {  };
	const int ArgIndices[] = {  };
	
	if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "InitGaussRand") != passed)  return 1;
    
    rand_seed = clock();
    
    gsl_rng_env_setup();
    gsl_rand_gen_type = gsl_rng_default;
    gsl_rand_gen = gsl_rng_alloc(gsl_rand_gen_type);
    if ((void *) gsl_rand_gen == NULL)  return 1;
    gsl_rng_set(gsl_rand_gen, (ccInt) rand_seed);
    
    return passed;
}



int call_GaussRand(int argc, char **argv)
{
    double *result, sigma;
	arg_info *ArgInfo = (arg_info *) argv[argc];
    
	const int ArgTypes[] = { double_type, double_type };
	const int ArgIndices[] = { 1, 1 };
	
	if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "GaussRand") != passed)  return 1;
    
    getArgs(argc, argv, byValue(&sigma), &result);
    
    *result = gsl_ran_gaussian(gsl_rand_gen, sigma);
    
    return 0;
}







// *********** Diagnostics ***********


int call_Entropy(int argc, char **argv)
{
    ccInt c1, alpha_base, alpha, one_color, *C2F;
	arg_info *ArgInfo = (arg_info *) *(argv+argc);
    double *info, prob_sum, *first_p, one_p, x_res, y_res, z_res, alpha_x, alpha_y, alpha_z, one_tot_p;
	ccBool IfAvg, IfZeroNorm, countFalseNegatives;
    
	const int ArgTypes[] = { double_type, double_type, double_type, string_type, int_type,
                                    int_type, int_type, int_type, double_type, double_type, double_type, bool_type, double_type };
	const int ArgIndices[] = { -1, -1, -1, -2, -2, -3, -3, -4, 1, 1, 1, 1, 1 };
	
	if (CheckArgInfo(ArgInfo, ArgTypes, ArgIndices, argc, sizeof(ArgTypes)/sizeof(int), "Entropy") != passed)  return 1;
	
    getArgs(argc, argv, &(image.x), &(image.y), &(image.z), &a3d_params.p, &contour.color, &image.color_bottoms, &image.color_tops, &C2F,
                byValue(&x_res), byValue(&y_res), byValue(&z_res), byValue(&countFalseNegatives), &info);
    
    contour.numSpots = ArgInfo[3].argIndices;
    
    IfZeroNorm = ccFalse;
    if (x_res*y_res*z_res == 0)  IfZeroNorm = ccTrue;
    else  {
        alpha_x = 1. / (2 * x_res * x_res);
        alpha_y = 1. / (2 * y_res * y_res);
        alpha_z = 1. / (2 * z_res * z_res);         }
    
    IfAvg = ccFalse;
    if (ArgInfo[7].argIndices == 0)  IfAvg = ccTrue;
    else if (ArgInfo[7].argIndices != ArgInfo[3].argIndices)  {
        printf("Entropy() error:  C2F and p[] must be the same length");
        return 1;           }
    
    *info = 0;
    for (c1 = 0; c1 < contour.numSpots; c1++)  {
    if (AreSpots(c1))  {
    if (a3d_params.p[c1].elementNum > 0)  {
        if (a3d_params.p[c1].memory == 0)  {
            printf("Entropy() error:  p[] is not initialized\n");
            return 2;           }
        
        first_p = LL_Double(a3d_params.p + c1, 1);
        one_color = contour.color[c1];
        
        prob_sum = one_tot_p = 0;
        alpha_base = image.color_bottoms[one_color] - 1;
        for (alpha = alpha_base; alpha < image.color_tops[one_color]; alpha++)     {
            one_p = first_p[alpha - alpha_base];
            prob_sum += one_p;
            
            if (IfAvg)  {
                if (one_p > 0)  *info -= one_p*log(one_p);    }
            
            else  {
                if ((!IfZeroNorm) && (C2F[c1] >= 1))
                    one_tot_p += one_p * exp(-SqDotDistance(alpha, C2F[c1] - 1, alpha_x, alpha_y, alpha_z));
                else if ((IfZeroNorm) && (alpha == C2F[c1] - 1))
                    one_tot_p += one_p;
        }   }
        
        if ((countFalseNegatives) && (prob_sum >= 0))      {
            if (IfAvg)  {
                if (prob_sum < 1.)  *info -= (1. - prob_sum)*log(1. - prob_sum);         }
            else  {
                if (C2F[c1] == 0)  *info -= log(1. - prob_sum);
                else if (one_tot_p < 1.)  *info -= log(one_tot_p);        // screen out the > 1 case
        }   }
    }}}
    
	return passed;
}


ccBool AreSpots(ccInt i)
{
    if ((i < 0) || (i >= contour.numSpots))  return ccTrue;
    
    return (image.color_tops[contour.color[i]] - image.color_bottoms[contour.color[i]] + 1 > 0);
}
