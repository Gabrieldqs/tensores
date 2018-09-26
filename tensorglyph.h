/*************************************************************************************
    GCG - GROUP FOR COMPUTER GRAPHICS
        Universidade Federal de Juiz de Fora
        Instituto de Ci�ncias Exatas
        Departamento de Ci�ncia da Computa��o

gcgTENSORGLYPH
   Marcelo Bernardes Vieira
   T�ssio Knop de Castro
**************************************************************************************/

#ifndef _gcgTENSORGLYPH_H_
#define _gcgTENSORGLYPH_H_


//#define MIN_LSCALE 0.1  //min length scale
//#define MIN_ALPHA  0.1
//#define MIN_BETA   0.1

#include "gcg.h"
#include "squadric.h"

class gcgTENSORGLYPH{
    public:


    //properties of the tensor matrix
    unsigned int hasParticle;
    unsigned int capacity;
        
    MATRIX3 eigenVectors;
    VECTOR3 eigenValues;
    
    float alpha, beta, gamma;
    float key;
    int vi, vj, vk;
    //anisotropy components
    float cl, cp, cs;                   //linear, planar and spherical

    //tensor matrix
    MATRIX3 tensor;
    //position of its center
    VECTOR3 pos;

    //k1 = 1 - |e1 . obs|
    float k1;

    //k2 = 1 - |e2 . obs|
    float k2;


    //k3 = |e3 . obs|
    float k3;

    //  1st central moment of eigenvalues
    float m1;

    //  2nd central moment of eigenvalues
    float m2;

    //skewness of eigenvalues
    float A3;

    //invariant of tensor D
    float J4;

    //fractional anisotropy
    float FA;

    //relative anisotropy
    float RA;


    //number of the tensorline
    int tnumber;

    //curvatura values
    float curv;

    //After tensor_ord keep the orignal
    gcgTENSORGLYPH *origin;

    //value computed from the eigenvalues to determine the color of the glyph: its 1/(sum(eigenvalues))
    float weight;

    float confidence;

    //Distance to the observer
    float dist2Obs;

    //Distance to the outermous
    float dist2Otm;

    //Its an outermoust?
    float otm;

    //Result of the equation
    float equ;
    
    int createdHere;
    
    bool unableSite;

    //constructor
    //obtains alpha and beta from a tensor, extracting its eigenvalues and computing
    //its anisotropy components (linear, planar, spherical)
    //if linear >= planar, we use the first parametrization; otherwise we use the second
//    gcgTENSORGLYPH(MATRIX3 tensor, VECTOR3 pos, float confidence);
    
    gcgTENSORGLYPH(MATRIX3 tensor, VECTOR3 pos, float confidence, bool isOriginal, int numSubGlyphs, gcgTENSORGLYPH*** subGlyphs, int voi, int voj, int vok, int fwidth, int fheight, int fdepth);
    void setTensor(MATRIX3 tensor);
    gcgTENSORGLYPH();
    ~gcgTENSORGLYPH();

    void draw(float scale = 1);
    
    
    /******SQUADRIC******/
//    bool parametrization;
//    gcgSUPERQUADRIC * sq;
//    void create(unsigned int slicesphi,unsigned int slicestheta, bool drawNormals, bool drawMesh);
//    void draw(bool drawMesh,bool drawNormals, float scale = 1);
    
};

//void gcgHeatColor(float normheat, VECTOR3 color);
#endif


