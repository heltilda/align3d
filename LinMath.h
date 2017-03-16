/*
 *  LinMath.h
 *  Yazoo
 *
 *  Created by Brian Ross on 1/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LinMathHeader
#define LinMathHeader

#include "cicada.h"

typedef struct {
    double real;
    double imag;
} cdouble;

typedef struct {
    double x;
    double y;
    double z;
} vector_3;

typedef struct {
    cdouble x;
    cdouble y;
    cdouble z;
} cvector_3;

extern void vAdd(double *, double *, ccInt, double *);
extern void v3Add(vector_3 *, vector_3 *, vector_3 *);
extern void v3Sub(vector_3 *, vector_3 *, vector_3 *);
extern void cvAdd(cdouble *, cdouble *, ccInt, cdouble *);
extern void vSub(double *, double *, ccInt, double *);
extern void cvSub(cdouble *, cdouble *, ccInt, cdouble *);
extern void vScalarMult(double *, ccInt, double, double *);
extern void v3ScalarMult(vector_3 *, double, vector_3 *);
extern void cvScalarMult(cdouble *, ccInt, cdouble, cdouble *);
extern void vScalarDiv(double *, ccInt, double, double *);
extern void cvScalarvDiv(cdouble *, ccInt, cdouble, cdouble *);
extern double vDot(double *, double *, ccInt);
extern double dot(vector_3 *, vector_3 *);
extern vector_3 cross(vector_3 *, vector_3 *);
extern cdouble cDot(cvector_3 *, cvector_3 *);
extern cvector_3 cCross(cvector_3 *, cvector_3 *);
extern void GramSchmidt(double *, double *, ccInt);
extern void v3GramSchmidt(vector_3 *, vector_3 *);
extern cdouble cAdd(cdouble, cdouble);
extern cdouble cSub(cdouble, cdouble);
extern cdouble cMult(cdouble, cdouble);
extern cdouble cDiv(cdouble, cdouble);
extern cdouble crMult(cdouble, double);
extern cdouble cConj(cdouble);
extern void v3Get(vector_3 *, double *);
extern void v3Put(vector_3 *, double *);

extern const cdouble c_zero;
extern const cdouble c_one;
extern const cdouble c_i;
extern const cdouble c_minus_one;
extern const cdouble c_minus_i;

#endif
