#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "math.h"

#define MAX_MAT_SIZE 1000

typedef struct mat {

	int mat[MAX_MAT_SIZE];
	int n, m;
}mat;

typedef unsigned int uint;

void initTestMat(mat* A);
void printMat(mat* A);
void initMask(mat* A, uint k);
float distance(int a, int b, int c, int d);
void padMat(mat* A, uint xpad, uint ypad);
void printMatAddress(mat* A);
void unpadMat(mat* A, uint xpad, uint ypad);
int eleWiseMult(int* A, uint n_a, uint m_a, int* B, uint n_b, uint m_b, uint xMult, uint yMult);
mat* corr2D(mat* A, mat* B);
void occupiedSlots(mat* corr, mat* A);
uint findMax(mat* corr);

void main()
{
	mat A, mask, *corr;
	uint k,max;
	k = 1;
	printf("\nmat A:\n");
	initTestMat(&A);
	printMat(&A);
	initMask(&mask,k);
	printf("\nMask:\n");
	printMat(&mask);

	corr = corr2D(&A, &mask);
	occupiedSlots(corr, &A);
	printf("\nCorrelation:\n");
	printMat(corr);
	max = findMax(corr);
	printf("max value of correlation = %d\n", max);

}
void initTestMat(mat* A)
{
	(A->mat[0]) = 1;
	(A->mat[1]) = 0;
	(A->mat[2]) = 0;
	(A->mat[3]) = 0;

	(A->mat[4]) = 0;
	(A->mat[5]) = 1;
	(A->mat[6]) = 1;
	(A->mat[7]) = 0;

	(A->mat[8]) = 0;
	(A->mat[9]) = 0;
	(A->mat[10]) = 0;
	(A->mat[11]) = 1;

	A->n = 3;
	A->m = 4;
}
void printMat(mat* A)
{
	uint i, j;
	for (i = 0; i < A->n; i++)
	{
		for (j = 0; j < A->m; j++)
		{
			printf("| %d |", *(A->mat + A->m*i + j));
		}
		printf("\n");
	}
}
void printMatAddress(mat* A)
{
	uint i, j;
	for (i = 0; i < A->n; i++)
	{
		for (j = 0; j < A->m; j++)
		{
			printf("| %p |", (A->mat + A->m*i + j));
		}
		printf("\n");
	}
}
void initMask(mat* A, uint k)
{
	uint i, j;
	A->n = (k << 1) + 1;
	A->m = (k << 1) + 1;

	for (i = 0; i < A->n; i++)
	{
		for (j = 0; j < A->m; j++)
		{
			if (k >= distance(i, j, k, k))
			{
				*(A->mat + A->m*i + j) = 1;
			}
			else
			{
				*(A->mat + A->m * i + j) = 0;
			}
		}
	}
}
float distance(int a, int b, int c, int d)
{
	return (sqrt(pow((c - a), 2) + pow((d - b), 2)));
}
void padMat(mat* A, uint xpad, uint ypad)
{
	uint i,j;
	mat tmp;
	tmp.m = A->m;
	tmp.n = A->n;
	memcpy(tmp.mat, A->mat, sizeof(int) * ((A->m) * (A->n)));

	A->m += (xpad << 1);
	A->n += (ypad << 1);

	for (i = 0; i < A->n; i++)
	{
		for (j = 0; j < A->m; j++)
		{
			if ((i < ypad) || (j < xpad) || (i > (tmp.n - 1) + ypad) || (j > (tmp.m - 1) + xpad))
			{
				*(A->mat + (A->m) * i + j) = 0;
			}
			else
			{
				*(A->mat + (A->m) * i + j) = *(tmp.mat + tmp.m * (i - ypad) + (j - xpad));
			}
		}
	}
}
void unpadMat(mat* A, uint xpad, uint ypad)
{
	uint i, j, new_n, new_m;

	new_n = A->n - (ypad << 1);
	new_m = A->m - (xpad << 1);

	for (i = 0; i < new_n; i++)
	{
		for (j = 0; j < new_m; j++)
		{
			*(A->mat + new_m * i + j) = *(A->mat + (A->m) * (i + ypad) + (j + xpad));
		}
	}

	A->m = new_m;
	A->n = new_n;

}
int eleWiseMult(int* A,uint n_a, uint m_a, int* B,uint n_b, uint m_b, uint xMult, uint yMult)
{
	uint i, j,sum;

	sum = 0;
	for (i = 0; i < yMult; i++)
	{
		for (j = 0; j < xMult; j++)
		{
			sum += ((*(A + m_a * i + j)) * (*(B + m_b * i + j)));
		}
	}
	return sum;
}
mat* corr2D(mat* A, mat* B)
{
	uint i, j, xpad, ypad;
	mat corr;

	corr.n = A->n;
	corr.m = A->m;

	xpad = (B->m - 1)>>1;
	ypad = (B->n - 1)>>1;

	padMat(A, xpad, ypad);

	for (i = 0; i < corr.n; i++)
	{
		for (j = 0; j < corr.m; j++)
		{
			*(corr.mat + corr.m * i + j) = eleWiseMult((A->mat + A->m * i + j), A->n, A->m, B->mat, B->n, B->m, B->m, B->n);
		}
	}
	unpadMat(A, xpad, ypad);
	return &corr;
}
void occupiedSlots(mat* corr, mat* A)
{
	uint i, j;

	for (i = 0; i < corr->n; i++)
	{
		for (j = 0; j < corr->m; j++)
		{
			if (*(A->mat + A->m * i + j) == 1)
			{
				*(corr->mat + corr->m * i + j) = 0;
			}
		}
	}
}
uint findMax(mat* corr)
{
	uint i, j,max;

	max = *(corr->mat);
	for (i = 0; i < corr->n; i++)
	{
		for (j = 0; j < corr->m; j++)
		{
			if (*(corr->mat + corr->m * i + j) > max)
			{
				max = *(corr->mat + corr->m * i + j);
			}
		}
	}
	return max;
}