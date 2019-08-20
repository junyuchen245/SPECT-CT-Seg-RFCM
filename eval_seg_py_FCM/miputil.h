/*
	miputil.h
	
	$Id: miputil.h 101 2009-03-03 15:00:38Z frey $
*/

#ifndef MIPUTIL_H
#define MIPUTIL_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#if MALLOC_H_NEEDED
#include <malloc.h>
#endif

#define NR_END 1
#define FREE_ARG char*

/* Prototypes */

#ifdef __cplusplus
extern "C" {
#endif

float dot(float *vec1, float *vec2, int length);
void set_float(float *array, int number, float constant);
void set_int(int *array, int number, int constant);
void set_double(double *array, int number, double constant);
double sum_double(double *array, int number);
float sum_float(float *array, int number);
void scale_float(float *array, int number, float fac);
int sum_int(int *array, int number);
float sum_short(short *array, int number);
int *ivector(int);
float *vector(int);
double *dvector(int);
void free_ivector(int *);
void free_vector(float *);
void free_dvector(double *);
float *pfDupVec(float *pfVec, int iVecLen);
void vZeroNeg(float *array, int number);

#ifdef __cplusplus
}
#endif

#endif /* MIPUTIL_H */
