/*
	miputil.c
	
	$Id: miputil.c 101 2009-03-03 15:00:38Z frey $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef MALLOC_H_NEEDED
#include <malloc.h>
#endif
#include "miputil.h"

float dot(float *vec1, float *vec2, int length)
{
  int     i;
  float dotProduct = 0.0;
  
  for (i=0; i<length; i++)
    dotProduct += vec1[i]*vec2[i];
  
  return(dotProduct);
}

void set_float(float *array, int number, float constant)
{
  int   i;
  for (i=0; i<number; i++)
    array[i] = constant;
}

float sum_float(float *array, int number)
{
  int   i;
  double sum = 0.0;
  for (i=0; i<number; i++)
    sum += array[i];
  return((float)sum);
}

float sum_short(short *array, int number)
{
  int   i;
  double sum = 0.0;
  for (i=0; i<number; i++)
    sum += (float)array[i];
  return((float)sum);
}

void scale_float(float *array, int number, float fac)
{
	int   i;
	for (i=0; i<number; i++)
		array[i] *= fac;
}

float *pfDupVec(float *pfVec, int iVecLen)
{
	float *pfOutVec;

	pfOutVec=vector(iVecLen);
	if (pfOutVec == NULL) return NULL;
	memcpy(pfOutVec, pfVec, iVecLen*sizeof(float));
	return pfOutVec;
}
	
void set_int(int *array, int number, int constant)
{  
  int   i;
  for (i=0; i<number; i++)
    array[i] = constant;
}

int sum_int(int *array, int number)
{
  int	i;
  int	sum = 0;
  for (i=0; i<number; i++)
    sum += array[i];
  return(sum);
}

void vZeroNeg(float *array, int number)
{
	int i;
	for(i=0; i<number; i++)
		if (array[i] < 0)
			array[i] = 0;
}

void set_double(double *array, int number, double constant)
{
  int   i;
  for (i=0; i<number; i++)
    array[i] = constant;
}

double sum_double(double *array, int number)
{
  int   i;
  double sum = 0.0;
  for (i=0; i<number; i++)
    sum += array[i];
  return(sum);
}

double *dvector(int iLen)
     /* allocates a double precision vector w/length iLen */
{
  return((double *)malloc(sizeof(double)*iLen));
}


void free_dvector(double *pdVec)
{
  if (pdVec != NULL) free(pdVec);
}

float *vector(int iLen)
     /* allocates a single precision vector w/length iLen */
{
  return((float *)malloc(sizeof(float)*iLen));
}

void free_vector(float *pfVec)
{
  if (pfVec != NULL) free(pfVec);
}

int *ivector(int iLen)
     /* allocates an integer vector w/length iLen */
{
  return((int *)malloc(sizeof(int)*iLen));
}

void free_ivector(int *piVec)
{
  if (piVec != NULL) free(piVec);
}

