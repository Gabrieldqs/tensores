/////*********************************************************************************************
/***************** SMOOTHED PARTICLE HYDRODYNAMICS ********************************************
*** T�ssio Knop de Castro
*** Gildo de Almeida Leonel
***********************************************************************************************
particle.h: defines the properties and procedures of each particle
***********************************************************************************************/

/*
 * Massa específica
 * Termo de Pressão
 * Força de Pressão

 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "gcg.h"
#include "geometryInterface.h"
#include "tensorglyph.h"
#include <math.h>
#include <cmath>

using namespace std;

class Octree;
void changeMi(bool add);
void changeK(bool add);
void changeDensity(bool add);
numberType calculateMass(numberType volume, numberType density, numberType nParticles);
void defineConstants();
char* printInverseTensor();
void changeTensorApplication();
float printViscosity();
float printStiffness();
float printRestDensity();
char * printEditMode();
void changeEditMode();
float frobeniusNorm(MATRIX3 t);

class gcgParticleSPH: public gcgPARTICLE{
  public:
  numberType radius;
  numberType mass;
  numberType massDensity;
  numberType lastMassDensity;
  numberType pressureTerm;
  numberType influenceRadius;
  numberType smoothWidth;
  numberType densityDerivate;
  VECTOR3 density;
  numberType densityVal;
  VECTOR3 gravity;
  VECTOR3 lastDensity;
  numberType densityDerivative;
  VECTOR3 pressure;
  numberType pressureVal;
  VECTOR3 lastPressure;
  VECTOR3 viscosity, lastViscosity;
  numberType direc;
  numberType simulate;
  numberType forces;
  VECTOR3 position, lastPosition;
  VECTOR3 dir;
  VECTOR3 externalForces;
//  bool isDead;
//  gcgTENSORGLYPH* ordRelated;
  VECTOR3 velocity; ///used in Euler integration.In Leapfrog, it must be computed from nextVel and prevVel
//  VECTOR3 nextVel, prevVel; ///nextVel = v(t + 1/2* timestep), prevVel = v(t - 1/2*timestep). Used in Leapfrog integration
  VECTOR3 correction;
  VECTOR3 acceleration;
  VECTOR3 color;
  VECTOR3 lastVelocity;
  VECTOR3 nextVelocity;
  VECTOR3 lastAcceleration;
  VECTOR3 normal;
  VECTOR3 surfaceTension;
  VECTOR3 buoyancy;
  VECTOR3 tensorForce;
  gcgTEXT text;
  int gridAdreess;

  void resetIteration();
  void checkCollision();

  /*************/
  VECTOR3 debugV;
  VECTOR3 debugDif;
  /*************/

//  VECTOR3 * getLennarJones();
  vector<gcgParticleSPH*> neighbours;
  bool neighbourCheck;
  gcgParticleSPH* next;
  int samePosCount;
  unsigned int simulationStep;
  VECTOR3 gradientVal;

  Octree * parentOctree;

  gcgParticleSPH(numberType mass, numberType radius, VECTOR3 pos,numberType direc);
  gcgParticleSPH(numberType mass, numberType radius, VECTOR3 pos, VECTOR3 initialVel);
  gcgParticleSPH();
  void addNeighbour(gcgParticleSPH* n);

  gcgTENSORGLYPH * tensor;

//private:
public:
  void defineParticle(numberType mass, numberType radius, VECTOR3 pos, VECTOR3 initialVel);


    void calculateViscosity();
    void calculateMassDensity(bool debug);
    void calculatePressure();
    void calculatePressureTerm();
    void calculateExternalForces(VECTOR3 v);
    void calculateGravity();
    void calculateAcceleration(gcgParticleSPH * neighbour);
    void computeAcceleration();
    void calculateNormal();
    void calculateSurfaceTension();
    void calculateBuoyancy(VECTOR3 v);
//    void calculateTensorForce(VECTOR3 direction, numberType sign, numberType cl);
    void calculateTensorForce();
    /******/
    
        /******SQUADRIC******/
    bool parametrization;
    gcgSUPERQUADRIC * sq;
    void defineSQuad(gcgTENSORGLYPH * glyph);
    void createSQuad(unsigned int slicesphi,unsigned int slicestheta, bool drawNormals, bool drawMesh, gcgTENSORGLYPH * glyph);
};
#endif
