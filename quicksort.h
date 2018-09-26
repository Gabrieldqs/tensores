#ifndef QUICKSORT_H
#define	QUICKSORT_H

#include <stdio.h>
#include <stdlib.h>
#include "tensorglyph.h"
#include "tensorgrid.h"


extern int ordSize;

#define M 10
typedef gcgTENSORGLYPH* Item;
//#define key(A) ( (O[12] * (A)->k1) + (O[0] * (A)->k2) + (O[1] * (A)->k3) + (O[2] * (A)->cp) + (O[3] * (A)->cl) + (O[4] * (A)->cs) + (O[5] * (A)->m1) + (O[6] * (A)->m2) + (O[7] * (A)->A3) + (O[8] * (A)->J4) + (O[9] * (A)->FA) + (O[10] * (A)->RA) + (O[11] * (A)->dist2Obs))
//#define key(A) ((A)->A3)
//#define key(A) ((A)->k1 + (A)->k2 + (A)->k3 + (A)->A3) //+ (A)->k3
#define key(A) ((A)->key) //+ (A)->k3
//#define key(A) (ordSize - (A)->numberOfVisits)
//#define key(A) ((A)->numberOfVisits)
#define less(A, B) ( (key(A)) < (key(B)) )
#define exch(A, B) { Item t = A; A = B; B = t; }
#define compexch(A, B) if (less(B, A)) exch(A, B)

void bubbleSort(gcgTENSORGLYPH** arr, int n, bool verbose);
void nqs(gcgTENSORGLYPH** arr, int low, int high, bool verbose);
void mergesort(gcgTENSORGLYPH** a, int low, int high);
void NquickSort(gcgTENSORGLYPH** arr, int left, int right, float* O);

void quicksortM32(Item a[], int l, int r, float* O /* = 0 */);


void partitionQK(gcgTENSORGLYPH** A, int left, int right, int* i, int* j);

void sortQK(gcgTENSORGLYPH** A, int left, int right) ;

int partition(gcgTENSORGLYPH **v, int inicio, int fim) ;

void quicksort(gcgTENSORGLYPH **v, int inicio, int fim) ;
void quicksortM3(gcgTENSORGLYPH **v, int inicio, int fim) ;

void quickSort(gcgTENSORGLYPH **A, int n);

void remakeHeap(int left, int right, gcgTENSORGLYPH** A, char type);

void createHeap(gcgTENSORGLYPH** A, int n, char type);


void swap3(gcgTENSORGLYPH **a, gcgTENSORGLYPH **b) ;

int partition3(gcgTENSORGLYPH **vec, int left, int right, float obX, float obY, float obZ) ;

void quickSort3(gcgTENSORGLYPH **vec, int left, int right, float obX, float obY, float obZ);

void bubble(gcgTENSORGLYPH **v, int qtd);
void swapbubble( gcgTENSORGLYPH **v, int i);
void insertionSort(gcgTENSORGLYPH **v, int n);
#endif	/* QUICKSORT_H */

