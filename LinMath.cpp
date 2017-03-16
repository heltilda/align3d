/*
 *  LinMath.c
 *  Yazoo
 *
 *  Created by Brian Ross on 1/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "LinMath.h"

void vAdd(double *v1, double *v2, int v_length, double *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)
        *(result+c1) = (*(v1+c1))+(*(v2+c1));
}

void v3Add(vector_3 *v1, vector_3 *v2, vector_3 *result)
{
    result->x = (v1->x)+(v2->x);
    result->y = (v1->y)+(v2->y);
    result->z = (v1->z)+(v2->z);
}

void v3Sub(vector_3 *v1, vector_3 *v2, vector_3 *result)
{
    result->x = (v1->x)-(v2->x);
    result->y = (v1->y)-(v2->y);
    result->z = (v1->z)-(v2->z);
}

void cvAdd(cdouble *v1, cdouble *v2, int v_length, cdouble *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)
        *(result+c1) = cAdd(*(v1+c1), *(v2+c1));
}

void vSub(double *v1, double *v2, int v_length, double *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)
        *(result+c1) = (*(v1+c1))-(*(v2+c1));
}

void cvSub(cdouble *v1, cdouble *v2, int v_length, cdouble *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)
        *(result+c1) = cSub(*(v1+c1), *(v2+c1));
}

void vScalarMult(double *vector, int v_length, double scalar, double *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)   {
        *(result+c1) = (*(vector+c1))*scalar;   }
}

void v3ScalarMult(vector_3 *vector, double scalar, vector_3 *result)
{
    result->x = (vector->x)*scalar;
    result->y = (vector->y)*scalar;
    result->z = (vector->z)*scalar;
}

void cvScalarMult(cdouble *vector, int v_length, cdouble scalar, cdouble *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)   {
        *(result+c1) = cMult(*(vector+c1), scalar);     }
}

void vScalarDiv(double *vector, int v_length, double scalar, double *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)   {
        *(result+c1) = (*(vector+c1))/scalar;   }
}

void cvScalarDiv(cdouble *vector, int v_length, cdouble scalar, cdouble *result)
{
    int c1;
    
    for (c1 = 0; c1 < v_length; c1++)   {
        *(result+c1) = cDiv(*(vector+c1), scalar);      }
}

double vDot(double *v1, double *v2, int length)
{
    int c1;
    double dp;

    dp = 0;
    for (c1 = 0; c1 < length; c1 ++)
        dp += (*(v1+c1)) * (*(v2+c1));

    return dp;
}

double dot(vector_3 *v1, vector_3 *v2)
{
    return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

vector_3 cross(vector_3 *v1, vector_3 *v2)      // don't pass return value as a 3rd argument
{                                               // that would overwrite the args while we're using them
    vector_3 result;
    
    result.x = v1->y*v2->z - v1->z*v2->y;
    result.y = v1->z*v2->x - v1->x*v2->z;
    result.z = v1->x*v2->y - v1->y*v2->x;
    
    return result;
}

cdouble cDot(cvector_3 *v1, cvector_3 *v2)
{
    return cAdd(cAdd(cMult(v1->x, v2->x), cMult(v1->y, v2->y)), cMult(v1->z, v2->z));
}

cvector_3 cCross(cvector_3 *v1, cvector_3 *v2)
{
    cvector_3 result;
    
    result.x = cSub(cMult(v1->y, v2->z), cMult(v1->z, v2->y));
    result.y = cSub(cMult(v1->z, v2->x), cMult(v1->x, v2->z));
    result.z = cSub(cMult(v1->x, v2->y), cMult(v1->y, v2->x));
    
    return result;
}

void GramSchmidt(double *v1, double *v2, int length)
{
    int c1;
    double norm;

    norm = vDot(v2, v2, length);
    if (norm == 0.)  return;

    norm = vDot(v1, v2, length) / norm;

    for (c1 = 0; c1 < length; c1++)
        *(v1+c1) -= norm * (*(v2+c1));
}

void v3GramSchmidt(vector_3 *v1, vector_3 *v2)
{
    double norm;
    
    norm = dot(v2, v2);
    if (norm == 0.)  return;
    
    norm = dot(v1, v2) / norm;
    
    v1->x -= norm * v2->x;
    v1->y -= norm * v2->y;
    v1->z -= norm * v2->z;
}

cdouble cAdd(cdouble cd1, cdouble cd2)
{
    cdouble sum;
    
    sum.real = cd1.real + cd2.real;
    sum.imag = cd1.imag + cd2.imag;
    return sum;
}

cdouble cSub(cdouble cd1, cdouble cd2)
{
    cdouble diff;
    
    diff.real = cd1.real - cd2.real;
    diff.imag = cd1.imag - cd2.imag;
    return diff;
}

cdouble cMult(cdouble cd1, cdouble cd2)
{
    cdouble prod;
    
    prod.real = cd1.real*cd2.real - cd1.imag*cd2.imag;
    prod.imag = cd1.real*cd2.imag + cd1.imag*cd2.real;
    return prod;
}

cdouble cDiv(cdouble cd1, cdouble cd2)
{
    double norm = 1./(cd2.real*cd2.real + cd2.imag*cd2.imag);
    cdouble quot;
    
    quot.real = norm*(cd1.real*cd2.real + cd1.imag*cd2.imag);
    quot.imag = norm*(cd1.imag*cd2.real - cd1.real*cd2.imag);
    return quot;
}

cdouble crMult(cdouble cd1, double d1)
{
    cdouble prod;
    
    prod.real = cd1.real*d1;
    prod.imag = cd1.imag*d1;
    return prod;
}

cdouble cConj(cdouble cd1)
{
    cdouble conjugate;
    
    conjugate.real = cd1.real;
    conjugate.imag = -cd1.imag;
    return conjugate;
}

void v3Get(vector_3 *v3, double *vN)
{
    v3->x = *(vN);  
    v3->y = *(vN + 1);  
    v3->z = *(vN + 2);
}

void v3Put(vector_3 *v3, double *vN)
{
    *(vN) = v3->x;  
    *(vN + 1) = v3->y;  
    *(vN + 2) = v3->z;
}


const cdouble c_zero = { 0, 0 };
const cdouble c_one = { 1, 0 };
const cdouble c_i = { 0, 1 };
const cdouble c_minus_one = { 1, 0 };
const cdouble c_minus_i = { 0, 1 };
