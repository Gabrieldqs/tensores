/*************************************************************************************
    GCG - GROUP FOR COMPUTER GRAPHICS
        Universidade Federal de Juiz de Fora
        Instituto de Ci�ncias Exatas
        Departamento de Ci�ncia da Computa��o

gcgTENSORGLYPH
   Marcelo Bernardes Vieira
   T�ssio Knop de Castro
**************************************************************************************/

#include "tensorglyph.h"
#include <math.h>
#include <iostream>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

bool changeA3 = false;
float globalGamma = 1.;
#define ADDRESSGLYPHS2(a, b, c, fw, fh)  (((b) * (fw)  ) + (a) + ((c) * (fw)  * (fh)))


gcgTENSORGLYPH::gcgTENSORGLYPH(){
}

//obtains alpha and beta from a tensor, extracting its eigenvalues and computing
//its anisotropy components (linear, planar, spherical)
//if linear >= planar, we use the first parametrization; otherwise we use the second
gcgTENSORGLYPH::gcgTENSORGLYPH(MATRIX3 tensor, VECTOR3 pos, float confidence, bool isOriginal, int numSubGlyphs, gcgTENSORGLYPH*** subGlyphs, int voi, int voj, int vok, int fwidth, int fheight, int fdepth){
    gcgCOPYMATRIX3(this->tensor,tensor);
    gcgEigenSymmetricMatrix3(tensor,eigenVectors,eigenValues);
    //eigenvalues sorted in decreasing order, but this function returns them in increasing order
    //obtaining eigenvectors
//    if(eigenValues[0] > 0) printf("%f - %f - %f\n", eigenValues[0], eigenValues[1], eigenValues[2]);



//    float aux = eigenValues[0];
//    eigenValues[0] = eigenValues[2];
//    eigenValues[2] = aux;
//    VECTOR3 tmp;
//    gcgCOPYVECTOR3(tmp, &eigenVectors[6]);
//    gcgCOPYVECTOR3(&eigenVectors[6], &eigenVectors[0]);
//    gcgCOPYVECTOR3(&eigenVectors[0], tmp);
    hasParticle = 0;
    unableSite = false;

    
    FLOAT aux = eigenValues[0];
    eigenValues[0] = eigenValues[2];
    eigenValues[2] = aux;
    VECTOR3 vaux = {eigenVectors[0],eigenVectors[1],eigenVectors[2]};
    eigenVectors[0] = eigenVectors[6];
    eigenVectors[1] = eigenVectors[7];
    eigenVectors[2] = eigenVectors[8];
    gcgCROSSVECTOR3(&eigenVectors[6], eigenVectors, &eigenVectors[3]);

    
    
    this->pos[0]   = pos[0];
    this->pos[1]   = pos[1];
    this->pos[2]   = pos[2];
    this->confidence = confidence;
    this->dist2Obs = 0.0;
    this->k1 = 0.0;
    this->k2 = 0.0;
    this->k3 = 0.0;
    this->curv = 0.0;
    this->tnumber = 0.0;
    this->otm = 0.0;
    this->dist2Otm = 0.0;
    this->vi = voi;
    this->vj = voj;
    this->vk = vok;
    this->createdHere = 0;
    this->alpha = 0;
    this->beta = 0;
    this->gamma = globalGamma;
    this->key = 0;
    this->capacity = 0;
    
    
    //anisotropy components
//    if(eigenValues[0] + eigenValues[1] + eigenValues[2] > 1.0) { // Para manter entre 0 e 1
    if(!FEQUAL((eigenValues[0] + eigenValues[1] + eigenValues[2]), 0.0)) { // P/ Zerar se der divisão por 0
      cl = (eigenValues[0] - eigenValues[1])/ (eigenValues[0] + eigenValues[1] + eigenValues[2]);
      cp = ((eigenValues[1] - eigenValues[2]) + (eigenValues[1] - eigenValues[2]))/ (eigenValues[0] + eigenValues[1] + eigenValues[2]);
      cs = (eigenValues[2] + eigenValues[2] + eigenValues[2])/ (eigenValues[0] + eigenValues[1] + eigenValues[2]);

      m1 = (eigenValues[0] + eigenValues[1] + eigenValues[2])/3;
      m2 = (((eigenValues[0] - m1)*(eigenValues[0] - m1)) + ((eigenValues[1] - m1)*(eigenValues[1] - m1)) + ((eigenValues[2] - m1)*(eigenValues[2] - m1)))/3;

      if(fabs(m2) < 0.002) A3 = 0;
      else A3 =   (((eigenValues[0] - m1)*(eigenValues[0] - m1)*(eigenValues[0] - m1))
            + ((eigenValues[1] - m1)*(eigenValues[1] - m1)*(eigenValues[1] - m1))
            + ((eigenValues[2] - m1)*(eigenValues[2] - m1)*(eigenValues[2] - m1)))
            / (3 * m2 * sqrt(m2));

      J4 =  (eigenValues[0]*eigenValues[0]) + (eigenValues[1]*eigenValues[1]) + (eigenValues[2]*eigenValues[2]);
      FA = (3/sqrt(2)) * sqrt(m2 / J4);
      RA = sqrt(m2/(2*m1));
    } else {
        J4 = FA = RA = A3 = cl = cp = 0; cs = 1;
    
    }

    if(cl > 1. || cs > 1. || cp > 1.) {J4 = FA = RA = A3 = cl = cp = 0; cs = 1; }

    //set the pointers to NULL
    this->origin = this;


    //making the eigenvalues vary from 0 to 1
    gcgNORMALIZEVECTOR3(eigenValues,eigenValues);
    
    //computing the tensor's weight
    weight =(1/ gcgLENGTHVECTOR3(eigenValues));
//     weight =gcgLENGTHVECTOR3(eigenValues);
    
    /*****SQUADRIC*****/
    
//    this->gamma = gamma;
//    
//    FLOAT alpha, beta;
//    if (cl >= cp){
//        parametrization = FIRST_FORM;
//        alpha = 1 - cp;
//        alpha = pow(alpha,gamma); //sharpening alpha
//        beta = 1 - cl;
//        beta = pow(beta,gamma);   //sharpening beta
//    }
//    else{
//        parametrization = SECOND_FORM;
//        alpha = 1 - cl;
//        alpha = pow(alpha,gamma); //sharpening alpha
//        beta = 1 - cp;
//        beta = pow(beta,gamma);   //sharpening beta
//    }
//    if (alpha < MIN_ALPHA) alpha = MIN_ALPHA;
//    if (beta < MIN_BETA) beta = MIN_BETA;
//    sq = new gcgSUPERQUADRIC(alpha,beta,pos,parametrization,gamma);
//    sq->create(7, 7, false, false);
    
}


void gcgTENSORGLYPH::setTensor(MATRIX3 tensor){
    gcgCOPYMATRIX3(this->tensor,tensor);
    gcgEigenSymmetricMatrix3(tensor,eigenVectors,eigenValues);
    //eigenvalues sorted in decreasing order, but this function returns them in increasing order
    //obtaining eigenvectors
//    if(eigenValues[0] > 0) printf("%f - %f - %f\n", eigenValues[0], eigenValues[1], eigenValues[2]);



//    float aux = eigenValues[0];
//    eigenValues[0] = eigenValues[2];
//    eigenValues[2] = aux;
//    VECTOR3 tmp;
//    gcgCOPYVECTOR3(tmp, &eigenVectors[6]);
//    gcgCOPYVECTOR3(&eigenVectors[6], &eigenVectors[0]);
//    gcgCOPYVECTOR3(&eigenVectors[0], tmp);
    
    FLOAT aux = eigenValues[0];
    eigenValues[0] = eigenValues[2];
    eigenValues[2] = aux;
    VECTOR3 vaux = {eigenVectors[0],eigenVectors[1],eigenVectors[2]};
    eigenVectors[0] = eigenVectors[6];
    eigenVectors[1] = eigenVectors[7];
    eigenVectors[2] = eigenVectors[8];
    gcgCROSSVECTOR3(&eigenVectors[6], eigenVectors, &eigenVectors[3]);
    
    
    hasParticle = 0;
    unableSite = false;

    this->pos[0]   = pos[0];
    this->pos[1]   = pos[1];
    this->pos[2]   = pos[2];
    this->confidence = confidence;
    this->dist2Obs = 0.0;
    this->k1 = 0.0;
    this->k2 = 0.0;
    this->k3 = 0.0;
    this->curv = 0.0;
    this->tnumber = 0.0;
    this->otm = 0.0;
    this->dist2Otm = 0.0;
    this->vi = 0;
    this->vj = 0;
    this->vk = 0;
    this->createdHere = 0;
    this->alpha = 0;
    this->beta = 0;
    this->gamma = globalGamma;
    this->key = 0;
    this->capacity = 0;
    

    //anisotropy components
//    if(eigenValues[0] + eigenValues[1] + eigenValues[2] > 1.0) { // Para manter entre 0 e 1
    if(!FEQUAL((eigenValues[0] + eigenValues[1] + eigenValues[2]), 0.0)) { 
      cl = (eigenValues[0] - eigenValues[1])/ (eigenValues[0] + eigenValues[1] + eigenValues[2]);
      cp = ((eigenValues[1] - eigenValues[2]) + (eigenValues[1] - eigenValues[2]))/ (eigenValues[0] + eigenValues[1] + eigenValues[2]);
      cs = (eigenValues[2] + eigenValues[2] + eigenValues[2])/ (eigenValues[0] + eigenValues[1] + eigenValues[2]);

      m1 = (eigenValues[0] + eigenValues[1] + eigenValues[2])/3;
      m2 = (((eigenValues[0] - m1)*(eigenValues[0] - m1)) + ((eigenValues[1] - m1)*(eigenValues[1] - m1)) + ((eigenValues[2] - m1)*(eigenValues[2] - m1)))/3;

      if(fabs(m2) < 0.002) A3 = 0;
      else A3 =   (((eigenValues[0] - m1)*(eigenValues[0] - m1)*(eigenValues[0] - m1))
            + ((eigenValues[1] - m1)*(eigenValues[1] - m1)*(eigenValues[1] - m1))
            + ((eigenValues[2] - m1)*(eigenValues[2] - m1)*(eigenValues[2] - m1)))
            / (3 * m2 * sqrt(m2));

      J4 =  (eigenValues[0]*eigenValues[0]) + (eigenValues[1]*eigenValues[1]) + (eigenValues[2]*eigenValues[2]);
      FA = (3/sqrt(2)) * sqrt(m2 / J4);
      RA = sqrt(m2/(2*m1));
    } else {J4 = FA = RA = A3 = cl = cp = 0; cs = 1; }

    if(cl > 1. || cs > 1. || cp > 1.) {J4 = FA = RA = A3 = cl = cp = 0; cs = 1; }
        
    
    //set the pointers to NULL
    this->origin = this;
    
    
    //making the eigenvalues vary from 0 to 1
    gcgNORMALIZEVECTOR3(eigenValues,eigenValues);
    
    //computing the tensor's weight
    weight =(1/ gcgLENGTHVECTOR3(eigenValues));
    
    /*****SQUADRIC*****/
    
//    this->gamma = gamma;
    
//    FLOAT alpha, beta;
//    if (cl >= cp){
//        parametrization = FIRST_FORM;
//        alpha = 1 - cp;
//        alpha = pow(alpha,gamma); //sharpening alpha
//        beta = 1 - cl;
//        beta = pow(beta,gamma);   //sharpening beta
//    }
//    else{
//        parametrization = SECOND_FORM;
//        alpha = 1 - cl;
//        alpha = pow(alpha,gamma); //sharpening alpha
//        beta = 1 - cp;
//        beta = pow(beta,gamma);   //sharpening beta
//    }
//    if (alpha < MIN_ALPHA) alpha = MIN_ALPHA;
//    if (beta < MIN_BETA) beta = MIN_BETA;
//    sq = new gcgSUPERQUADRIC(alpha,beta,pos,parametrization,gamma);
//    sq->create(7, 7, false, false);
}

//void gcgTENSORGLYPH::create(unsigned int slicesphi,unsigned int slicestheta, bool drawNormals, bool drawMesh){
//    //in diffusion tensors, the length scales are defined by the eigenvalues
//    FLOAT xScale = (eigenValues[0] > MIN_LSCALE)? eigenValues[0]: MIN_LSCALE,
//      yScale = (eigenValues[1] > MIN_LSCALE)? eigenValues[1]: MIN_LSCALE,
//      zScale = (eigenValues[2] > MIN_LSCALE)? eigenValues[2]: MIN_LSCALE;
//
//    sq->create(slicesphi,slicestheta, xScale*0.49,yScale*0.49,zScale*0.49); //deixando 0.01 de espaço entre um glifo e seu vizinho, no caso onde eles tenham comprimento maximo
//}

gcgTENSORGLYPH::~gcgTENSORGLYPH(){

}

/*
void gcgHeatColor(float normheat, VECTOR3 color) {
  #define GCG_HEAT_COLORS 6
  const static float colors[GCG_HEAT_COLORS][3] = { {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 0.0},
                                                    {1.0, 1.0, 0.0}, {1.0, 0.5, 0.0}, {1.0, 0.0, 0.0} };
  float colorpos = normheat * (GCG_HEAT_COLORS - 1);
  int   index = (int) colorpos;
  if(index < 0) {gcgCOPYVECTOR3(color, colors[0]); }
  else if(index > GCG_HEAT_COLORS - 2) { gcgCOPYVECTOR3(color, colors[GCG_HEAT_COLORS - 1]); }
       else {
         float  s1 = colorpos - (float) index, s0 = 1.0 - s1;
         gcgSETVECTOR3(color, colors[index][0] * s0 + colors[index + 1][0] * s1,
                              colors[index][1] * s0 + colors[index + 1][1] * s1,
                              colors[index][2] * s0 + colors[index + 1][2] * s1 );
       }
}
*/

//void gcgTENSORGLYPH::draw(float scale){
//    VECTOR3 color;
//
//    gcgHeatColor(weight,color);
//    glColor4f(color[0],color[1],color[2], confidence);
//    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//    glPushMatrix();
//    glTranslatef(pos[0],pos[1],pos[2]);
//    glScalef(scale,scale,scale);
//    MATRIX4 orientation = {eigenVectors[0], eigenVectors[1],eigenVectors[2], 0,
//                       eigenVectors[3],eigenVectors[4],eigenVectors[5], 0,
//                       eigenVectors[6],eigenVectors[7],eigenVectors[8], 0,
//                       0, 0, 0, 1};
//    glMultMatrixf(orientation);
//
//    glBegin(GL_LINE);
//    glutSolidSphere (0.05, 16, 16);
//    glEnd( );
//
//    glPopMatrix();
//}

void gcgTENSORGLYPH::draw(float scale){
    VECTOR3 color;
    /*    if (weight < 0.25){
        color[0] = 1;
        color[1] = 0;
        color[2] = 0;
    }
    else if (weight < 0.50){
        color[0] = 0;
        color[1] = 1;
        color[2] = 0;
    }
        else if (weight < 0.75){
           color[0] = 0;
           color[1] = 0;
           color[2] = 1;
        }
            else if (weight < 1){
                color[0] = 1;
                color[1] = 1;
                color[2] = 0;
            }~*/
    gcgHeatColor(weight,color);
    glColor4f(color[0],color[1],color[2], confidence);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    glTranslatef(pos[0],pos[1],pos[2]);
    glScalef(scale,scale,scale);
    MATRIX4 orientation = {eigenVectors[0], eigenVectors[1],eigenVectors[2], 0,
                       eigenVectors[3],eigenVectors[4],eigenVectors[5], 0,
                       eigenVectors[6],eigenVectors[7],eigenVectors[8], 0,
                       0, 0, 0, 1};
    glMultMatrixf(orientation);
//    sq->draw();
    glPopMatrix();
}
