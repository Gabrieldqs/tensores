/*********************************************************************************************
***************** SMOOTHED PARTICLE HYDRODYNAMICS ********************************************
*** T�ssio Knop de Castro
*** Gildo de Almeida Leonel
***********************************************************************************************
particle.cpp: defines the properties and procedures of each particle
***********************************************************************************************/

#include "particle.h"
#include "sph.h"
#include "squadric.h"


#define GAS 0
#define LIQUID 1

#define MIN_LSCALE 0.1  //min length scale
#define MIN_ALPHA  0.1
#define MIN_BETA   0.1

extern int fluidType;

VECTOR3 gravitationalConstant = {0.,-9.82, 0.};
extern float simulationStept;
extern bool supportDistortion;
numberType restDensity;
numberType buoyancyDiffusion;
numberType mi; /// Viscosity
numberType sphAlpha; ///Surface tension
int ESTIMATED_NEIGHBOURS;
numberType kConst; //Gas stiffness
bool inverseTensorApplication = false;
char useInverseTensor[6];
extern bool simula;

char* printInverseTensor(){
    return useInverseTensor;
}

float printViscosity(){
    return mi;
}

float printStiffness(){
    return kConst;
}

float printRestDensity(){
    return restDensity;
}

float frobeniusNorm(MATRIX3 t){
    
    return sqrtf((t[0]*t[0]) + (t[1]*t[1]) + (t[2]*t[2]) + (t[3]*t[3]) + (t[4]*t[4]) + (t[5]*t[5]) + (t[6]*t[6]) + (t[7]*t[7]) + (t[8]*t[8]));
    
}

void changeTensorApplication(){
    inverseTensorApplication = !inverseTensorApplication;
    if(inverseTensorApplication){
        sprintf(useInverseTensor, "YES");
    }
    else{
        sprintf(useInverseTensor, "NO");
    }
}

void printMatrix(MATRIX3 m){
    printf("%f %f %f %f %f %f %f %f %f\n", m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
}

void changeK(bool add){
    if(add){ kConst += 0.1;}
    else{
        kConst -= 0.1;
        if(kConst < 0.1) kConst = 0.1;
    }
    
    printf("Stiffness changed to: %f\n", kConst);
}

void changeMi(bool add){
    if(add){ mi += 0.1;}
    else{
        if(mi < 0.1) mi = 0.1;
        mi -= 0.1;
    }
    
    printf("Viscosity coef. changed to: %f\n", mi);
}

void changeDensity(bool add){
    if(add) restDensity += 0.01;
    else restDensity -= 0.01;
    if(restDensity < 0.01) restDensity = 0.1;
    printf("Rest density changed to: %f\n", restDensity);
}

void defineConstants(){
    if(inverseTensorApplication){
        sprintf(useInverseTensor, "true");
    }
    else{
        sprintf(useInverseTensor, "false");
    }
    /***********WATER***********/
if (fluidType == LIQUID){
    restDensity = 998.29;
    buoyancyDiffusion = 5;
    mi = 3.5; /// Viscosity
    sphAlpha = 0.0728; ///Surface tension
    ESTIMATED_NEIGHBOURS = 20; //3d
//    ESTIMATED_NEIGHBOURS = 7;
    kConst = 3;
}
else{
/***********STEAM***********/
    restDensity = 0.59;
    buoyancyDiffusion = 5;
//    mi = 0.01; /// Viscosity
//    mi = 0.6; /// Viscosity
    mi = 0.3; /// Viscosity
    sphAlpha = 0.0; ///Surface tension
    ESTIMATED_NEIGHBOURS = 12; //3d
//    ESTIMATED_NEIGHBOURS = 20; //3d
//    ESTIMATED_NEIGHBOURS = 9;
    kConst = 4; 
}

}

//#define Cr 0.0

//numberType supportRadius = 0.0457;

numberType supportRadius;
//numberType supportRadius = 0.84194515;
//numberType supportRadius = 0.9;


numberType influenceFactor = 1;
numberType influenceRadius;
//numberType influenceRadius = influenceFactor * supportRadius;
//numberType supportRadius = 0.0457;
//numberType supportRadius3 = supportRadius*supportRadius*supportRadius;
//numberType supportRadius6 = supportRadius3*supportRadius*supportRadius*supportRadius;
//numberType supportRadius9 = supportRadius6*supportRadius*supportRadius*supportRadius;
numberType supportRadius3, supportRadius6, supportRadius9;



extern numberType TIME_STEP;


//#define mi 0.000103 /// Viscosity

#define sphL sqrtf(restDensity/ESTIMATED_NEIGHBOURS)
#define Cr 0.0

////////////////////////////
//#include "tensorglyph.h"
float distortDistance(gcgParticleSPH * p, gcgParticleSPH * neighbour, VECTOR3 * v){
//    VECTOR3 r1, r2, tmpV, finalVec;
//    MATRIX3 tensor;
//    gcgADDMATRIX3(tensor, p->tensor->tensor, neighbour->tensor->tensor);
//    gcgSCALEMATRIX3(tensor, tensor, 0.5);
//    gcgSUBVECTOR3(tmpV, p->position, neighbour->position);
//    gcgAPPLYMATRIX3VECTOR3(finalVec, tensor, tmpV);
//    gcgCOPYVECTOR3(*v, finalVec);
////        
//    return gcgLENGTHVECTOR3(finalVec);
    
    
    
    VECTOR3 r1, r2, tmpV, finalVec;
    
//        
    gcgSUBVECTOR3(tmpV, p->position, neighbour->position);
    if(inverseTensorApplication){
        MATRIX3 t1, t2, mID;
        gcgIDENTITYMATRIX3(mID);
        gcgSUBMATRIX3(t1, mID, p->tensor->tensor);
        gcgSUBMATRIX3(t2, mID, neighbour->tensor->tensor);
//        gcgInverseMatrix3(t1, p->tensor->tensor);
//        gcgInverseMatrix3(t2, neighbour->tensor->tensor);
//        printf("Norma das inversas: %f <-> %f\n", frobeniusNorm(t1), frobeniusNorm(t2));
        gcgAPPLYMATRIX3VECTOR3(r1, t1, tmpV);
        gcgAPPLYMATRIX3VECTOR3(r2, t2, tmpV);
    }
    else{
        gcgAPPLYMATRIX3VECTOR3(r1, p->tensor->tensor, tmpV);
        gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);
    }
    
    gcgSETVECTOR3(finalVec, (r1[0]+r2[0])/2.0, (r1[1]+r2[1])/2.0, (r1[2]+r2[2])/2.0);
    gcgCOPYVECTOR3(*v, finalVec);
    return gcgLENGTHVECTOR3(finalVec);
    
}

numberType calculateMass(numberType volume, numberType density, int nParticles){
    return (density*(volume/(numberType)nParticles));
}

numberType kernelDefault(numberType dist){
    if(dist > supportRadius) return 0.0f;
    numberType value = (supportRadius*supportRadius)-(dist*dist);
    value = value*value*value;
    value *= (315/(64*M_PI*supportRadius9));
    return value;
}


void kernelDefaultDerivative(numberType dist, VECTOR3 * vec){
    if(dist > supportRadius){ gcgSETVECTOR3(&vec, 0,0,0);}
    else{
        numberType value = (supportRadius*supportRadius)-(dist*dist);
        value = value*value;
        value *= (-945/(32*M_PI*supportRadius9));
        gcgSCALEVECTOR3(*vec, *vec, value);
    }
}

numberType kernelDefaultLaplacian(numberType dist){
    if(dist > supportRadius) return 0.0f;
    else{
        numberType value = (supportRadius*supportRadius)-(dist*dist);
        value *= (-945/(32*M_PI*supportRadius9)) * ((3*supportRadius*supportRadius)-(dist*dist*7));
        return value;
    }
}

///Ver o support radius -> influcente radiuu

///Não usando o vetor - trocar
void kernelPressureDerivative(numberType dist, VECTOR3 * vec){
    if(dist > supportRadius){ gcgSETVECTOR3(&vec, 0,0,0);}
    else{
        numberType value = supportRadius-dist;
        value = value*value;
        value *= (-45.0/(M_PI*supportRadius6));
        gcgSCALEVECTOR3(*vec, *vec, value);
    }
}

numberType kernelViscosityLaplacian(numberType dist){
    if(dist > supportRadius) return 0.0f;
    numberType value = supportRadius-dist;
    value *= (45/(M_PI*supportRadius6));
    return value;
}

void gcgParticleSPH::defineSQuad(gcgTENSORGLYPH* glyph){
    
    
    float alpha, beta;
    if (glyph->cl >= glyph->cp){
        parametrization = FIRST_FORM;
        alpha = 1 - glyph->cp;
        alpha = pow(alpha,glyph->gamma); //sharpening alpha
        beta = 1 - glyph->cl;
        beta = pow(beta,glyph->gamma);   //sharpening beta
    }
    else{
        parametrization = SECOND_FORM;
        alpha = 1 - glyph->cl;
        alpha = pow(alpha,glyph->gamma); //sharpening alpha
        beta = 1 - glyph->cp;
        beta = pow(beta,glyph->gamma);   //sharpening beta
    }
    if (alpha < MIN_ALPHA) alpha = MIN_ALPHA;
    if (beta < MIN_BETA) beta = MIN_BETA;
    sq = new gcgSUPERQUADRIC(alpha,beta,glyph->pos,parametrization,glyph->gamma);
//    createSQuad(1, 1, false, false, glyph);
//    createSQuad(30, 30, false, false, glyph);
//    createSQuad(12, 12, false, false, glyph);
//    if(simula){
//        createSQuad(1, 1, false, false, glyph);
//    }
//    else{
//        createSQuad(20, 20, false, false, glyph);
//    }
    createSQuad(7, 7, false, false, glyph);
//    createSQuad(1, 1, false, false, glyph);
//    sq->create(7, 7, false, false);
    
}

void gcgParticleSPH::createSQuad(unsigned int slicesphi,unsigned int slicestheta, bool drawNormals, bool drawMesh, gcgTENSORGLYPH * glyph){
    //in diffusion tensors, the length scales are defined by the eigenvalues
    FLOAT xScale = (glyph->eigenValues[0] > MIN_LSCALE)? glyph->eigenValues[0]: MIN_LSCALE,
      yScale = (glyph->eigenValues[1] > MIN_LSCALE)? glyph->eigenValues[1]: MIN_LSCALE,
      zScale = (glyph->eigenValues[2] > MIN_LSCALE)? glyph->eigenValues[2]: MIN_LSCALE;

    sq->create(slicesphi,slicestheta, xScale*0.49,yScale*0.49,zScale*0.49); //deixando 0.01 de espaço entre um glifo e seu vizinho, no caso onde eles tenham comprimento maximo
}

void gcgParticleSPH::calculateTensorForce(){
    
    
    gcgSETVECTOR3(tensorForce, - tensor->tensor[0] + tensor->tensor[1] - tensor->tensor[3] + tensor->tensor[4], - tensor->tensor[0] - tensor->tensor[1] + tensor->tensor[3] + tensor->tensor[4], 0);
    
//    gcgSETVECTOR3(tensorForce, tensor->eigenVectors[0] * tensor->cl, tensor->eigenVectors[1] * tensor->cl, tensor->eigenVectors[2] * tensor->cl);
    
    
    
    
    
    
//    VECTOR3 vIn, vOut, vProp, v1, vAux, nextPos;
//    MATRIX2 scaleD;
//    int posX, posY, posZ;
//    numberType wPunct = 0.9;
//    numberType dt = 0.05;
//
//    gcgSETVECTOR3(vIn, 0.0, 0.0, 0.0);
////    gcgCOPYVECTOR3(vIn, velocity);
//    gcgSCALEMATRIX3(scaleD, glyphs[addressG]->tensor, (2.0f / bigger));
//    gcgAPPLYMATRIX3VECTOR3(vOut, scaleD, vIn);
//    
////    printf("Bigger: %f\n", bigger);
//
//    gcgSETVECTOR3(v1, glyphs[addressG]->eigenVectors[0], glyphs[addressG]->eigenVectors[1], glyphs[addressG]->eigenVectors[2])
//
//
//    gcgSCALEVECTOR3(v1, v1, glyphs[addressG]->cl);
//    gcgSCALEVECTOR3(vIn, vIn, 1.0 - wPunct);
//    gcgSCALEVECTOR3(vOut, vOut, wPunct);
//    gcgADDVECTOR3(vAux, vIn, vOut);
//    gcgSCALEVECTOR3(vAux, vAux, 1.0 - glyphs[addressG]->cl);
//    gcgADDVECTOR3(vProp, vAux, v1);
//    
//    gcgSCALEVECTOR3(vProp, vProp, sign * dt);
////    gcgSCALEVECTOR3(vProp, vProp,(-sign) * dt);
////    
//////     gcgSCALEVECTOR3(scaleVel, velField[pos], (p->direc * dt)); //dt = 0.2
//////    gcgADDVECTOR3(p->pos, p->pos, scaleVel);
//////    geometry->particles[i] = p;
//////    gcgCOPYVECTOR3(geometry->particles[i]->vel, scaleVel);
////    
//
//
//    
////    VECTOR3 vProp;
////    numberType dt = 0.2;
//    
////    gcgSCALEVECTOR3(vProp, direction,  dt);
//
//    
//    
//    gcgCOPYVECTOR3(tensorForce, vProp);
//    
////    printf("Force: %f %f %f\n", tensorForce[0],tensorForce[1], tensorForce[2]);
//    
//    
////    gcgSCALEVECTOR3(tensorForce, direction, sign * cl);
    
}

///Talvez rever e recalcular vizinhanca antes
void gcgParticleSPH::checkCollision(){
    
    VECTOR3 differenceVector;
    numberType dist;
    VECTOR3 neighbourNormal, thisNormal, x, cp, aux;
    numberType collisionDepth;
    
    VECTOR3 curPosition, curVelocity;
    
    gcgSETVECTOR3(curVelocity, 0.0, 0.0, 0.0);
    int count  = 0;
    
    gcgCOPYVECTOR3(curPosition, this->position);
    
    for(std::vector<gcgParticleSPH*>::iterator it = this->neighbours.begin(); it != this->neighbours.end(); ++it) {
    
        gcgParticleSPH * neighbour = *it;
        
//        gcgSETVECTOR3(differenceVector, (position[0]-neighbour->position[0]), (position[1]-neighbour->position[1]), (position[2]-neighbour->position[2]));
        gcgSUBVECTOR3(differenceVector, position, neighbour->position);
        
        dist = gcgLENGTHVECTOR3(differenceVector);
//        dist = sqrtf((differenceVector[0]*differenceVector[0])+(differenceVector[1]*differenceVector[1])+(differenceVector[2]*differenceVector[2]));
        
        if(dist < (this->radius + neighbour->radius)){
            
            
            gcgSCALEVECTOR3(neighbourNormal, differenceVector, 1./dist);
            gcgSCALEVECTOR3(thisNormal, neighbourNormal, -1.);
            
            
            gcgSCALEVECTOR3(cp, neighbourNormal, neighbour->radius);
            gcgADDVECTOR3(cp, cp, neighbour->position);
            
            gcgSCALEVECTOR3(x, thisNormal, this->radius);
            gcgADDVECTOR3(x, x, this->position);
            
            
            
            ///Calclar profundidade da colisao
            gcgSUBVECTOR3(aux, neighbour->position, x); /// (c - x) ---- para o vizinho
            collisionDepth = fabs(gcgLENGTHVECTOR3(aux)-neighbour->radius);
            
            ///Calcular nova velocidade
            if(FEQUAL(0.0, gcgLENGTHVECTOR3(this->nextVelocity))){
                gcgSETVECTOR3(aux, 0.0, 0.0, 0.0);
            }
            else{
                gcgSCALEVECTOR3(aux, neighbourNormal, (gcgDOTVECTOR3(this->nextVelocity, neighbourNormal))*(1.0 + (collisionDepth*Cr/(TIME_STEP*gcgLENGTHVECTOR3(this->nextVelocity)))));
            }
            
//            gcgSUBVECTOR3(this->nextVelocity, this->nextVelocity, aux);
            gcgSUBVECTOR3(curVelocity, curVelocity, aux);
            count++;
            
            gcgSCALEVECTOR3(curPosition, neighbourNormal, this->radius); /////Conferir a projeção
            
            gcgADDVECTOR3(curPosition, curPosition, cp);
            
//            break;
//            gcgSETVECTOR3(curPosition, )
            
            
        }
    
    }
//    if(count == 0) return;
    
//    gcgSCALEVECTOR3(curVelocity, curVelocity, 1./((numberType)count));
    
    gcgSUBVECTOR3(this->nextVelocity, this->nextVelocity, curVelocity);
//    gcgCOPYVECTOR3(this->position, curPosition);
    
}

void gcgParticleSPH::resetIteration(){
    neighbourCheck = false;
    neighbours.clear();
    lastMassDensity = massDensity;
    massDensity = 0.;
    gcgCOPYVECTOR3(lastViscosity, viscosity);
    gcgCOPYVECTOR3(lastPressure, pressure);
    gcgSETVECTOR3(pressure, 0., 0., 0.);
    gcgSETVECTOR3(viscosity, 0., 0., 0.);
    gcgSETVECTOR3(externalForces, 0.,0.,0.);
    gcgSETVECTOR3(tensorForce, 0.,0.,0.);
    
}


gcgParticleSPH::gcgParticleSPH(numberType mass, numberType radius, VECTOR3 pos, numberType direc /* = 0.0 */){
  static int counter = 0;
    this->index = counter++;
  ///this->pressure = 1.0; ///QUAL � A PRESS�O INICIAL??
  this->mass = mass;
  this->radius = radius;
  this->direc = direc;
  gcgSETVECTOR3(velocity,0.,0.,0.);
  gcgSETVECTOR3(lastVelocity,0.,0.,0.);
  gcgSETVECTOR3(nextVelocity,0.,0.,0.);
  gcgSETVECTOR3(acceleration,0,0,0);
  gcgSETVECTOR3(dir, 0.0, 0.0, 0.0);
  simulationStep = 0;
  gcgCOPYVECTOR3(this->position, pos);
  gcgCOPYVECTOR3(this->startPos, pos);
  samePosCount = 0;
  simulate = true;
  static gcgRANDOM random;
  gcgHeatColor(random.random(),color);
  parentOctree = NULL;
  next = NULL;
  gcgSETVECTOR3(viscosity, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(pressure, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(density, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(externalForces, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(gravity, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(surfaceTension, 0.0, 0.0, 0.0);
  massDensity = 0;
  lastMassDensity = 0;
}

gcgParticleSPH::gcgParticleSPH(numberType mass, numberType radius, VECTOR3 pos, VECTOR3 initialVel){
  static int counter = 0;
  this->index = counter++;
  ///this->pressure = 1.0; ///QUAL � A PRESS�O INICIAL??
  this->mass = mass;
//  this->mass = 25;
  this->radius = radius;
  
  this->tensor = new gcgTENSORGLYPH();
//  this->radius = pow(0.75*mass/(restDensity * M_PI), 1.0/3.0);
//  this->smoothWidth = 1.3 * 3.0 *radius;

  influenceRadius = 1 * influenceFactor;
//  printf("Influence: %f\n", influenceRadius);
  neighbourCheck = false;
  gcgCOPYVECTOR3(velocity,initialVel);
  gcgSETVECTOR3(lastVelocity,0.,0.,0.);
  gcgSETVECTOR3(nextVelocity,0.,0.,0.);
  gcgSETVECTOR3(acceleration,0,0,0);
  gcgCOPYVECTOR3(this->position, pos);
  gcgSETVECTOR3(color,0.0,0.0,1.0);
  next = NULL;
  gcgSETVECTOR3(viscosity, 0.0, 0., 0.);
  gcgSETVECTOR3(pressure, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(density, 0.0, 0., 0.);
  gcgSETVECTOR3(externalForces, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(gravity, 0.0, 0.0, 0.0);
  gcgSETVECTOR3(surfaceTension, 0.0, 0.0, 0.0);
  parentOctree = NULL;
  massDensity = 0;
  lastMassDensity = 0;
  neighbours.clear();
//  printf("Particula criada: %f %f %f\n", position[0], position[1], position[2]);
}


void gcgParticleSPH::defineParticle(numberType mass, numberType radius, VECTOR3 pos, VECTOR3 initialVel){
  printf("Vai criar\n");
  static int counter = 0;
  this->index = counter++;
  ///this->pressure = 1.0; ///QUAL � A PRESS�O INICIAL??
  this->mass = mass;
//  this->mass = 25;
  this->radius = radius;
//  this->radius = pow(0.75*mass/(restDensity * M_PI), 1.0/3.0);
//  this->smoothWidth = 1.3 * 3.0 *radius;

  influenceRadius = 1 * influenceFactor;
//  printf("Influence: %f\n", influenceRadius);
  neighbourCheck = false;
  gcgCOPYVECTOR3(velocity,initialVel);
  gcgSETVECTOR3(lastVelocity,0.,0.,0.);
  gcgSETVECTOR3(nextVelocity,0.,0.,0.);
  gcgSETVECTOR3(acceleration,0,0,0);
  gcgCOPYVECTOR3(this->position, pos);
  gcgSETVECTOR3(color,0.0,0.0,1.0);
  next = NULL;
  gcgSETVECTOR3(viscosity, 0.0, 0., 0.);
  gcgSETVECTOR3(pressure, 0.0, 0., 0.);
  gcgSETVECTOR3(density, 0.0, 0., 0.);
  parentOctree = NULL;
  printf("Criada\n");
  neighbours.clear();
//  printf("Particula criada: %f %f %f\n", position[0], position[1], position[2]);
}

gcgParticleSPH::gcgParticleSPH(){
 static int counter = 0;
  this->index = counter++;
  ///this->pressure = 1.0; ///QUAL � A PRESS�O INICIAL??
  this->mass = 0.;
  this->radius = 0.;
  gcgSETVECTOR3(velocity,0,0,0);
  gcgSETVECTOR3(acceleration,0,0,0);
  gcgSETVECTOR3(this->position, 0,0,0);
  gcgSETVECTOR3(color,0.0,0.0,1.0);
  next = NULL;
  gcgSETVECTOR3(viscosity, 0.0, 0., 0.);
  gcgSETVECTOR3(pressure, 0.0, 0., 0.);
  gcgSETVECTOR3(density, 0.0, 0., 0.);
  parentOctree = NULL;
}
//VECTOR3 * gcgParticleSPH::getLennarJones(){
//    
//}

///correct
void gcgParticleSPH::addNeighbour(gcgParticleSPH* n){
  this->neighbours.push_back(n);
}



void gcgParticleSPH::calculateMassDensity(bool debug){
//    lastDensity = density;
    VECTOR3 differenceVector;
    numberType dist;
    numberType sum = 0.0;
    for(std::vector<gcgParticleSPH*>::iterator it = this->neighbours.begin(); it != this->neighbours.end(); ++it) {
    
        gcgParticleSPH * neighbour = *it;
        
        gcgSETVECTOR3(differenceVector, (position[0]-neighbour->position[0]), (position[1]-neighbour->position[1]), (position[2]-neighbour->position[2]));
        dist = gcgLENGTHVECTOR3(differenceVector);
        
        /**********************************/
        if(supportDistortion){
            
        if(debug){printf("Antes: %f %f %f <> %f --- ", differenceVector[0], differenceVector[1], differenceVector[2], dist);}
        VECTOR3 r1, r2, tmpV, finalVec;
//        
        gcgSUBVECTOR3(tmpV, position, neighbour->position);
        gcgAPPLYMATRIX3VECTOR3(r1, tensor->tensor, tmpV);
////        printf("%f %f %f - ", r1[0], r1[1], r1[2]);
        gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);
////        printf("%f %f %f ----- ", r2[0], r2[1], r2[2]);
////        gcgSCALEVECTOR3(tmpV, tmpV, -1.);
////        gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);
//        
//        MATRIX3 tempTensor;
//        gcgADDMATRIX3(tempTensor, tensor->tensor, neighbour->tensor->tensor);
//        gcgSCALEMATRIX3(tempTensor, tempTensor, 0.5);
        gcgSETVECTOR3(finalVec, (r1[0]+r2[0])/2.0, (r1[1]+r2[1])/2.0, (r1[2]+r2[2])/2.0);
//        
        
        dist = gcgLENGTHVECTOR3(finalVec);
        numberType d1 = dist;
        
        
        dist = distortDistance(this, neighbour, &differenceVector);
//        
//        if(FEQUAL(dist, d1)){
//        printf("%f %f %f - %f\n", finalVec[0], finalVec[1], finalVec[2], d1);
//        printf("%f %f %f - %f\n", differenceVector[0], differenceVector[1], differenceVector[2], dist);
//        }
        
        if(debug){
            printf("P1: %f %f %f - P2: %f %f %f\n", position[0], position[1], position[2], neighbour->position[0], neighbour->position[1], neighbour->position[2]);
            printf("T1: "); printMatrix(tensor->tensor);
            printf("T2: "); printMatrix(neighbour->tensor->tensor);
//            printf("Depois: %f %f %f <> %f\n", finalVec[0], finalVec[1], finalVec[2], dist);
        }
        }
//        sum += neighbour->mass * kernelDefault(gcgLENGTHVECTOR3(finalVec));
        /*************************************************/
        
//        gcgAPPLYMATRIX3VECTOR3(r1, tensor->tensor, tmpV);
//        printf("%f %f %f - ", r1[0], r1[1], r1[2]);
//        printf("%f %f %f\n", r2[0], r2[1], r2[2]);
        
//      
        numberType val = neighbour->mass * kernelDefault(dist);
//        printf("MassDensity: %f\n", val);
        
//        sum += neighbour->mass * kernelDefault(dist);
        sum += val;
    
    }
    
    massDensity = sum;
//    if(!FEQUAL(massDensity, 0.0)){ printf("Mass-density %f\n", massDensity);}
//    if(massDensity > 0.0){ printf("Mass-density %f - Dist %f\n", massDensity, dist);}
}

void gcgParticleSPH::calculatePressureTerm(){
    pressureTerm = kConst*(massDensity-restDensity);
//    pressureTerm = kConst*massDensity;
//    if(fluidType == LIQUID) pressureTerm = kConst*(massDensity-restDensity);
//    else pressureTerm = kConst*massDensity;
}

numberType sign(numberType a){
    if(a > 0) return 1.;
    else return -1.;
}

void gcgParticleSPH::calculatePressure(){
    
    VECTOR3 sum = {0,0,0};
    
    numberType PiDiv = pressureTerm/(massDensity*massDensity);
    
    for(std::vector<gcgParticleSPH*>::iterator it = this->neighbours.begin(); it != this->neighbours.end(); ++it) {
        
        gcgParticleSPH * neighbour = *it;
        
        if(FEQUAL(neighbour->massDensity, 0.0)){continue;} ///evitando escape de vizinho errado
        
        VECTOR3 differenceVector;
        gcgSETVECTOR3(differenceVector, (position[0]-neighbour->position[0]), (position[1]-neighbour->position[1]), (position[2]-neighbour->position[2]));
//        gcgSCALEVECTOR3(differenceVector, differenceVector, -1.);
        
        /**********/
        //Jeito 1
        MATRIX3 tensors;
        gcgADDMATRIX3(tensors, tensor->tensor, neighbour->tensor->tensor);
        
        //Jeito 2
        VECTOR3 v1, v2;
        gcgAPPLYMATRIX3VECTOR3(v1, tensor->tensor, differenceVector);
        gcgAPPLYMATRIX3VECTOR3(v2, neighbour->tensor->tensor, neighbour->position);
        gcgCOPYVECTOR3(debugV, v1);
        gcgCOPYVECTOR3(debugDif, differenceVector);
        gcgSCALEVECTOR3(debugDif, debugDif, -1.);
        gcgCOPYVECTOR3(neighbour->debugV, v2);
        
//        printf("%f %f %f\n",debugV[0], differenceVector[0], sign(differenceVector[0]));
//        printf("%f %f %f\n",debugV[1], differenceVector[1], sign(differenceVector[1]));
//        printf("%f %f %f\n",debugV[2], differenceVector[2], sign(differenceVector[2]));
        
//        debugV[0] *= sign(differenceVector[0]);
//        debugV[1] *= sign(differenceVector[1]);
//        debugV[2] *= sign(differenceVector[2]);
        
        
        
        /**********/
    
        
        /**********/
        //Jeito 1
        VECTOR3 tmpVec;
        gcgAPPLYMATRIX3VECTOR3(tmpVec, tensors, differenceVector);
//        gcgCOPYVECTOR3(differenceVector, tmpVec);
        
        //Jeito 2
//        gcgSUBVECTOR3(differenceVector, v1, v2);
        
        /**********/
        
        ///trocar pelos length
        numberType dist = sqrt((differenceVector[0]*differenceVector[0])+(differenceVector[1]*differenceVector[1])+(differenceVector[2]*differenceVector[2]));
        
        ///Tentando arrumar a criacao aleatoria
        if(dist > supportRadius || FEQUAL(dist, 0)){
//        if(dist > supportRadius){
            continue;
    //        gcgSETVECTOR3(pressure, 0.0, 0.0, 0.0);
        }
        else{
            numberType calc = PiDiv  + (neighbour->pressureTerm/(neighbour->massDensity*neighbour->massDensity));
            calc = calc * neighbour->mass;
            
//            printf("Antes: %f %f %f <> %f --- ", differenceVector[0], differenceVector[1], differenceVector[2], dist);
            
            gcgSCALEVECTOR3(differenceVector, differenceVector, 1./dist);
            
//            printf("%f %f %f\n", differenceVector[0], differenceVector[1], differenceVector[2]);
            
            /***************/
            
            if(supportDistortion){
            VECTOR3 r1, r2, tmpV, finalVec;
//            
            gcgSUBVECTOR3(tmpV, position, neighbour->position);
            gcgAPPLYMATRIX3VECTOR3(r1, tensor->tensor, tmpV);
//    //        printf("%f %f %f - ", r1[0], r1[1], r1[2]);
            gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);
//    //        printf("%f %f %f ----- ", r2[0], r2[1], r2[2]);
//    //        gcgSCALEVECTOR3(tmpV, tmpV, -1.);
//    //        gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);
//
            gcgSETVECTOR3(finalVec, (r1[0]+r2[0])/2.0, (r1[1]+r2[1])/2.0, (r1[2]+r2[2])/2.0);
            dist = gcgLENGTHVECTOR3(finalVec);
            gcgCOPYVECTOR3(differenceVector, finalVec);
//            printf("Depois: %f %f %f <> %f --- ", differenceVector[0], differenceVector[1], differenceVector[2], dist);
            
//            printf("%f %f %f\n", differenceVector[0], differenceVector[1], differenceVector[2]);
            dist = distortDistance(this, neighbour, &differenceVector);
            gcgSCALEVECTOR3(differenceVector, differenceVector, 1./dist);
            if(FEQUAL(dist, 0.)){continue;}
            }
            
            /*************/
            
            
            kernelPressureDerivative(dist, &differenceVector);
                    
            gcgSCALEVECTOR3(differenceVector, differenceVector, calc);
            
//            printf("Pressure: %f %f %f <> %f\n", differenceVector[0], differenceVector[1], differenceVector[2], gcgLENGTHVECTOR3(differenceVector));
            gcgADDVECTOR3(sum, sum, differenceVector);
        }
    }
    
    gcgSCALEVECTOR3(pressure, sum, (-1.0*massDensity));
    ///Contornando erro numerico na integracao
//    if(FEQUAL(gcgLENGTHVECTOR3(pressure), 0.0)){
//    if(FEQUAL(pressure[0], 0.0)){
//        pressure[0] = 0;
//    }
//    if(FEQUAL(pressure[1], 0.0)){
//        pressure[1] = 0;
//    }
//    if(FEQUAL(pressure[2], 0.0)){
//        pressure[2] = 0;
//    }
//    printf("fim\n");
    
}

///Tem diferenca entre formulas
void gcgParticleSPH::calculateViscosity(){
    
    VECTOR3 sum = {0,0,0};
    
    for(std::vector<gcgParticleSPH*>::iterator it = this->neighbours.begin(); it != this->neighbours.end(); ++it){
        
        gcgParticleSPH * neighbour = * it;
        
        if(FEQUAL(neighbour->massDensity, 0.0)){continue;} ///evitando escape de vizinho errado
        
        VECTOR3 differenceVector, velocityDif;
        
        gcgSETVECTOR3(differenceVector, (position[0]-neighbour->position[0]), (position[1]-neighbour->position[1]), (position[2]-neighbour->position[2]));

        numberType dist = sqrt((differenceVector[0]*differenceVector[0])+(differenceVector[1]*differenceVector[1])+(differenceVector[2]*differenceVector[2]));
        
        gcgSUBVECTOR3(velocityDif, neighbour->velocity, this->velocity);
        
        
        /***********************/
        if(supportDistortion){
        VECTOR3 r1, r2, tmpV, finalVec;

        gcgSUBVECTOR3(tmpV, position, neighbour->position);
        gcgAPPLYMATRIX3VECTOR3(r1, tensor->tensor, tmpV);
//        printf("%f %f %f - ", r1[0], r1[1], r1[2]);
        gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);
//        printf("%f %f %f ----- ", r2[0], r2[1], r2[2]);
//        gcgSCALEVECTOR3(tmpV, tmpV, -1.);
//        gcgAPPLYMATRIX3VECTOR3(r2, neighbour->tensor->tensor, tmpV);

        gcgSETVECTOR3(finalVec, (r1[0]+r2[0])/2.0, (r1[1]+r2[1])/2.0, (r1[2]+r2[2])/2.0);
        
//        printf("%f %f %f - ", finalVec[0], finalVec[1], finalVec[2]);
        
        
//        MATRIX3 tc;
//        gcgADDMATRIX3(tc, tensor->tensor, neighbour->tensor->tensor);
//        gcgSCALEMATRIX3(tc, tc, 1.0/2.0);
        
        
        dist = gcgLENGTHVECTOR3(finalVec);
            
          dist = distortDistance(this, neighbour, &differenceVector);
        
//        gcgCOPYVECTOR3(differenceVector, finalVec);
//        gcgSCALEVECTOR3(differenceVector, differenceVector, 1./dist);/
        
//        gcgAPPLYMATRIX3VECTOR3(finalVec, tc, tmpV);
//        printf("%f %f %f\n", finalVec[0], finalVec[1], finalVec[2]);
        }
        
        /**************/
        VECTOR3 v1,v2;
        ///Jeito 1
        
        
//        gcgSETVECTOR3(v1, tensor->eigenVectors[0], tensor->eigenVectors[1], tensor->eigenVectors[2]);
//        gcgSETVECTOR3(v2, neighbour->tensor->eigenVectors[0], neighbour->tensor->eigenVectors[1], neighbour->tensor->eigenVectors[2]);
//        
//        gcgSCALEVECTOR3(v1, v1,tensor->cl);
//        gcgSCALEVECTOR3(v2, v2, neighbour->tensor->cl);
//        
//        if(tensor->cl >1. || tensor->cl < 0.0) printf("eitaa\n");
//        
////        gcgADDVECTOR3(v1, velocity, v1);
////        gcgADDVECTOR3(v2, neighbour->velocity, v2);
//        
//        gcgSUBVECTOR3(v1, velocity, v1);
//        gcgSUBVECTOR3(v2, neighbour->velocity, v2);
        
//        gcgSUBVECTOR3(velocityDif, v2, v1);
        
        ///Jeito 2
        
//        gcgSETVECTOR3(v1, tensor->eigenVectors[0], tensor->eigenVectors[1], tensor->eigenVectors[2]);
//        gcgSETVECTOR3(v2, neighbour->tensor->eigenVectors[0], neighbour->tensor->eigenVectors[1], neighbour->tensor->eigenVectors[2]);
//        
//        gcgADDVECTOR3(v1,v2, v1);
//        gcgSCALEVECTOR3(v1, v1, 0.5);
//        
//        gcgSUBVECTOR3(velocityDif, neighbour->velocity, this->velocity);
//        gcgADDVECTOR3(velocityDif, velocityDif, v1);
        
        /**************/
        
//                    printf("Aceleracao: %f %f %f\n", this->velocity[0], this->velocity[1], this->velocity[2]);
//                    printf("Aceleracao: %f %f %f\n", neighbour->velocity[0], neighbour->velocity[1], neighbour->velocity[2]);
        
//        if(dist > supportRadius){
////            printf("Aqui!\n");
//            continue;
////            return;
//    //        gcgSETVECTOR3(viscosity, 0.0, 0.0, 0.0);
//        }
//        
//        else{
//            printf("v1: %f %f %f\nv2: %f %f %f\n", v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
    //        gcgSCALEVECTOR3(viscosity, differenceVector, ((mi/massDensity)*(neighbour->mass*kernelViscosityFunction(dist))));
            gcgSCALEVECTOR3(velocityDif, velocityDif, ((neighbour->mass/neighbour->massDensity)*kernelViscosityLaplacian(dist))); ///4.17 se todas tem msma massa especifica
//            gcgSCALEVECTOR3(velocityDif, velocityDif, (neighbour->mass*kernelViscosityLaplacian(dist))); ///4.19 se todas tem msma massa especifica
//            printf("Viscosity: %f %f %f <> %f\n", velocityDif[0], velocityDif[1], velocityDif[2], gcgLENGTHVECTOR3(velocityDif));
            gcgADDVECTOR3(sum, sum, velocityDif);
            
//        }
        
    }
    
//    gcgSCALEVECTOR3(viscosity, sum, (mi/massDensity)); //4.19
    gcgSCALEVECTOR3(viscosity, sum, mi); //4.17
}

void gcgParticleSPH::calculateNormal(){
  
    VECTOR3 sum = {0,0,0};
            
    for(std::vector<gcgParticleSPH*>::iterator it = this->neighbours.begin(); it != this->neighbours.end(); ++it){
        
        gcgParticleSPH * neighbour = * it;
        
        VECTOR3 differenceVector;
        
        gcgSETVECTOR3(differenceVector, (position[0]-neighbour->position[0]), (position[1]-neighbour->position[1]), (position[2]-neighbour->position[2]));
        
        numberType dist = sqrt((differenceVector[0]*differenceVector[0])+(differenceVector[1]*differenceVector[1])+(differenceVector[2]*differenceVector[2]));
        
        kernelDefaultDerivative(dist, &differenceVector);
        
        gcgSCALEVECTOR3(differenceVector, differenceVector, (neighbour->mass*neighbour->massDensity));
        
        gcgADDVECTOR3(sum, sum, differenceVector);
        
    }
    
    gcgCOPYVECTOR3(normal, sum);
    
}

void gcgParticleSPH::calculateGravity(){
//    gcgSCALEVECTOR3(gravity, gravitationalConstant, 1);
    gcgSCALEVECTOR3(gravity, gravitationalConstant, massDensity);
//    printf("Gravidade: %f %f %f\n", gravity[0], gravity[1], gravity[2]);
}

void gcgParticleSPH::calculateBuoyancy(VECTOR3 v){
    
//    gcgSCALEVECTOR3(buoyancy, gravitationalConstant, (buoyancyDiffusion*(massDensity-restDensity)));
    gcgSCALEVECTOR3(buoyancy, v, (buoyancyDiffusion*(massDensity-restDensity)));
}

void gcgParticleSPH::calculateExternalForces(VECTOR3 v){
    
    
    if (fluidType == GAS){
//        printf("GAS");
        if(!FEQUAL(massDensity, 0.));
        this->calculateBuoyancy(v);
        gcgCOPYVECTOR3(externalForces, buoyancy);
    }
    else{
//        printf("LIQUID");
        this->calculateGravity();
        this->calculateSurfaceTension();
//        gcgCOPYVECTOR3(externalForces, gravity);
        gcgADDVECTOR3(externalForces, gravity, surfaceTension);
    }
}
 

void gcgParticleSPH::calculateSurfaceTension(){

    this->calculateNormal();
    if(gcgLENGTHVECTOR3(this->normal) < sphL){
        gcgSETVECTOR3(this->surfaceTension, 0, 0, 0);
    }
    else{
        numberType sum = 0.0;
            
        for(std::vector<gcgParticleSPH*>::iterator it = this->neighbours.begin(); it != this->neighbours.end(); ++it){

            gcgParticleSPH * neighbour = * it;

            VECTOR3 differenceVector;

            gcgSETVECTOR3(differenceVector, (position[0]-neighbour->position[0]), (position[1]-neighbour->position[1]), (position[2]-neighbour->position[2]));

            numberType dist = sqrt((differenceVector[0]*differenceVector[0])+(differenceVector[1]*differenceVector[1])+(differenceVector[2]*differenceVector[2]));

            sum += kernelDefaultLaplacian(dist)*(neighbour->mass/neighbour->massDensity);

        }
        
        gcgSCALEVECTOR3(this->surfaceTension, this->normal, (1.0/gcgLENGTHVECTOR3(this->normal)));
        gcgSCALEVECTOR3(this->surfaceTension, this->surfaceTension, sum);
        
        
    }
    
}
//

/*
 Antiga pressao
 
            numberType calc = -((neighbour->mass*(pressureTerm+neighbour->pressureTerm)/(2*neighbour->massDensity)));
            kernelPressureFunction(dist, &differenceVector);
            gcgSCALEVECTOR3(differenceVector, differenceVector, calc);
            VECTOR3 aux;
    //        gcgCOPYVECTOR3(pressure, differenceVector);
    //        gcgSCALEVECTOR3(pressure, pressure, -1.0); ///Negativando pra navier stokes
            gcgCOPYVECTOR3(aux, differenceVector);
            gcgSCALEVECTOR3(aux, aux, -1.0); ///Negativando pra navier stokes
            gcgADDVECTOR3(pressure, pressure, aux);
 
 */