/* 
 * File:   sph.h
 * Author: joseluiz
 *
 * Created on July 11, 2012, 2:02 PM
 */

#ifndef SPH_H
#define	SPH_H

#include "gcg.h"
#include "octree.h"
#include "particle.h"
#include <time.h>
#include <math.h>
#include <cmath>
#include "gcgTensor.h"
#include "newGrid.h"
#include "squadric.h"
#include "superquadric.h"


float printParticlesRadius();
void printDebugVals();
char* printColour();
numberType getSRadius();
void toggleExternals();
char* printUseExternals();
void toggleDistortion();
char* printDistortion();
void changeRadius(bool add);
float printRadius();
void changeGradientScale(bool add);
float printGradientScale();
void changeEditMode();
void changeEditValue(bool add);
void importState(char * fName);



class gcgSph {
public:
//    gcgOptimizedGridSPH* geometry;                     ///manages the neighbours detection and its data structure
    unsigned int nParticles;
    unsigned int width, height, depth;                  ///dimensions of the bounding box
    unsigned int boxWidth, boxHeight, boxDepth;                  ///dimensions of the bounding box
    unsigned int widthFluid, heightFluid, depthFluid;   ///dimensions of the bounding box containing the fluid at the initialization
    numberType soundVelocity;                                ///it is also 10*sqrt(2 * g * H), as defined in eq. 4.6 #define SOUNDVELOCITY 171
    bool isblack;
    bool usedVolume;
    numberType alpha;
    numberType initialDensity;
//    numberType timeStep;
    VECTOR3 gravity;
    numberType dynamicViscosity;                             ///As defined in p.72. it's a constant
    numberType smoothWidth;   ///smoothWidth = h: directly linked with the numerical precision. If it is too small, there aren't enough neighbours.
                        ///if it is too big, some local properties will be globally smoothed. usually h = 1.3 * delta(x)
    numberType influenceFactor;                              ///influenceFactor = k. usually k = 2
    numberType influenceRadius;                              /// k*h radius
    numberType pRadius;

//    vector<GRIDCELL> sample;                            ///grid used for the initialization
//    FLOAT initialDistance;                              ///the initial distance between the particles
    
    numberType direction;
    
    numberType maxvel;     ///the modulus of the maximum velocity at the current iteration, used to refresh the timestep
    numberType mindensity; ///the modulus of the minimum density at the current iteration, used to refresh the timestep
    numberType mindt;      ///the modulus of the minimum density at the current iteration, used to refresh the timestep
    numberType particlesMass;


//    int n_loops; /* n_loops equation-solving for 1 neighbour search performed */

    ///precomputed constants to increase performance
    numberType RDERIVATIVE; ///constant used in the kernel (pre-computed to increase performance)
    numberType B;
    numberType INVDENSITY;
    void toggleBox();
    void drawBox();
    bool randomPosition(VECTOR3 * position);
    void resetParticles();
    gcgParticleSPH ** particles;
    double simulationTime;
    
    newGrid * grid;
    
    void changeColourType();
    
  private:
    Octree * octree;
    
    void checkBoundaryCondition(gcgParticleSPH * p, VECTOR3 at, VECTOR3 vt);

  /****** METHODS *****/
  public:
//    gcgSph(float mass, int nParticles, float density, unsigned int fluidWidth, unsigned int fluidHeight, unsigned int fluidDepth);
    gcgSph(int nParticles, unsigned int fluidWidth, unsigned int fluidHeight, unsigned int fluidDepth, unsigned int width, unsigned int height, unsigned int depth);
    gcgSph(int pX, int pY, int pZ, unsigned int fluidWidth, unsigned int fluidHeight, unsigned int fluidDepth, unsigned int width, unsigned int height, unsigned int depth) ;
    void iterate();
    void drawParticles();
    void drawParticlesSuperQuadric();
    void drawParticlesSuperQuadricTassio();
    void dandleCollisions();
    void saveState();
    void insertParticles();
    void exportState();
    virtual ~gcgSph();
        void interpolateTensor(gcgParticleSPH * p);
    void interpolateTensor3D(gcgParticleSPH * p, bool debug, int i);
//        gcgBasicSPH(gcgOptimizedGridSPH* geometry, float initialDensity, float dynamicViscosity, unsigned int nParticles, float gravity,
//                unsigned int widthFluid, unsigned int heightFluid, unsigned int depthFluid,
//                unsigned int widthBox, unsigned int heightBox, unsigned int depthBox, float smoothWidth = 1.3, float influenceFactor = 2.0);
//    ~gcgBasicSPH();
//
//    float kernel(float *dist);
//    float kernelDerivative(float *dist);

    ///LOOP SPH
//    void iterate();
//
//    void computeGhostPotential(FLOAT* g, gcgParticleSPH* p);
//
//    inline void computePressure();
//
//
//    void computeDensityDerivative(gcgPARTICLE *p);
//    void computeDensityDerivative(gcgParticleSPH *p);
//
//    void computeAcceleration(gcgPARTICLE* p);
//    void computeAcceleration(gcgParticleSPH* p);
//
//    ///integrates velocity, density and position
//    void euler(gcgPARTICLE *p);
//    void euler(gcgParticleSPH *p);
//
//    void leapfrog(gcgPARTICLE *p);
//    void leapfrog(gcgParticleSPH *p);
//
//    ///XSPH
//    void correctVelocity(gcgParticleSPH* p);
//
//    void applyBoundaryConditions(gcgPARTICLE *p);
//    void applyBoundaryConditions(gcgParticleSPH *p);
//
//    void numericalCorrections(gcgPARTICLE* p);
//    void numericalCorrections(gcgParticleSPH* p);
//
//    void refreshTimeCFL();
};


gcgSph * loadState(int val);

#endif	/* SPH_H */

