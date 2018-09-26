/*
 * File:   gcgTensor.h
 * Author: joseluiz
 *
 * Created on October 11, 2012, 10:48 AM
 */

#ifndef GCGTENSOR_H
#define	GCGTENSOR_H


#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "gcg.h"
#include "tensorglyph.h"
#include "tensorgrid.h"
#include "optimizedGrid.h"
#include "tensorgenerate.h"
#include "tensorHeap.h"
#include "quicksort.h"



#define CEREBRO 0
#define HELICE 1
#define PONTAS 2
#define FIELD2D 3
#define CEREBRO2 4
#define CEREBRO1 5
#define HELICE1 6
#define ADDRESSVOLUME(a, b, c, fw, fh)  (((b) * (fw) ) + (a) + ((c) * (fw) * (fh)))
#define ADDRESSGLYPHS(a, b, c, fw, fh)  (((b) * (fw)  ) + (a) + ((c) * (fw)  * (fh)))


class gcgTensor {
public:

    int fwidth, fheight, fdepth;
    int totalSize, nParticles;
    MATRIX3 * field;
    VECTOR3 *velField;
    float * norms;
    VECTOR3  * normDerivatives;
    VECTOR3 maxDerivatives, minDerivatives;
    float minDerivN, maxDerivN;
    void changeDrawing();
    void storeField(const char* fileName, MATRIX3 * field);
    int createTensorlineDup(VECTOR3* vField, gcgTENSORGLYPH** glyphs, int *nglyphs, unsigned int width, unsigned int height, unsigned int depth, int seedI, int seedJ, int seedK, float bigger, int i /* = 0 */);
    float * confidence;
    float bigger;
    int noriginals;
    int ordSize;
    float calculateNorm(MATRIX3 t);
    float maxNorm, minNorm;
    void drawField();
    void drawNorms();
    gcgTensor();
    void derivateNorms();
    void derivateNorms3D();
    void derivateNorms3D3();
    void derivateNorms3D2();
    void derivateNorms3D4();
    void drawFewNorms();
    void applyGaussian();
    gcgTensor(int width, int height, int depth, int nPart);
    virtual ~gcgTensor();
    void createTensorGlyphs(gcgTENSORGLYPH** glyphs, int *nglyphs, MATRIX3* tensorField, float* confidence, unsigned int width, unsigned int height, unsigned int depth, float *bigger, bool isOriginal);

    gcgTENSORGLYPH** glyphs_original;

    gcgTENSORGLYPH** tensor_ord;

//private:



    gcgOptimizedGridSPH * geometry;
    void readField();
    void readTXTField(int tamanho, const char* fileName, MATRIX3 *field);
    void readRawField(int tamanho, const char* fileName, MATRIX3 *field, float *confidence);
    void readTXTField3(int a, const char* fileName, MATRIX3 * field, int tamanho);
    void readNewField(int tamanho, const char* fileName, MATRIX3 *field);

    void setGlobalVars();


};

#endif	/* GCGTENSOR_H */

