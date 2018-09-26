/*
 * File:   gcgSph.cpp
 * Author: joseluiz
 *
 * Created on July 11, 2012, 2:02 PM
 */

#include <string.h>

#include "sph.h"

static gcgRANDOMGAUSSIAN gcgRandom;
float TIME_STEP = 0.05;
float TIME_STEP2 = TIME_STEP*TIME_STEP;
//bool supportDistortion = false;
unsigned int simulationStep = 0;
bool supportDistortion = false;
bool useExternals = true;
//bool supportDistortion = false;
#define radiusScale 1
#define VERLET 0
#define LEAP_FROG 1
#define ENERGY 0.7 // energia conservada
//#define ENERGY 1.0 // energia conservada
//#define INTEGRATOR VERLET
#define INTEGRATOR LEAP_FROG
#define GAS 0
#define LIQUID 1
//#define ADDRESS(a, b, c, dimX, dimY) (((c) * dimX * dimY + (b) * dimX + (a)) * 9)
#define ADDRESSVOLUME(a, b, c, fw, fh)  (((b) * (fw) ) + (a) + ((c) * (fw) * (fh)))
#define ADDRESSGLYPHS(a, b, c, fw, fh)  (((b) * (fw)  ) + (a) + ((c) * (fw)  * (fh)))
int state = 0;
int xport = 0;
bool filterParticles = false;
float filterT = 0.35;

//int fluidType = LIQUID;
int fluidType = GAS;
//#define INTEGRATOR VERLET
///Rever o leap nao funcionando

extern numberType supportRadius;
extern numberType influenceRadius;
extern numberType influenceFactor;
extern numberType supportRadius3;
extern numberType supportRadius6;
extern numberType supportRadius9;
extern int ESTIMATED_NEIGHBOURS;
extern numberType restDensity;
extern int tipoCampo;
extern numberType mi; /// Viscosity
extern numberType kConst;
extern VECTOR3 cutBounds1;
extern VECTOR3 cutBounds2;
extern bool cutVis;

gcgTensor * teste;

int translateW, translateH, translateD;

numberType minDebugVal = 9999999;
numberType maxDebugVal = -9999999;
VECTOR3 maxVec = {0,0,0}, minVec={0,0,0};
VECTOR3 debugTensor = {0,0,0};

unsigned int informationType = 3;
char colourId[15];
char areExternalsOn[4];
char isDistortionOn[4];
char editMode[11];
int editModeVal = 3;
float gradientScale = 6;
bool activateBoost = true;
float particlesRadius;
//bool slapX = true;
bool slapX = false;
//bool slapY = true;
bool slapY = false;
//bool slapZ = true;
bool slapZ = false;
bool slapEigen = false;
//bool slapEigen = false;
int slapMax = 20;

float printParticlesRadius(){
    return particlesRadius;
}

float printGradientScale(){
    return gradientScale;
}

void changeParticleRadius(bool add){
    if(add){
        particlesRadius += 0.003;
    }
    else{
        particlesRadius -= 0.003;
        if(particlesRadius < 0.001) particlesRadius = 0.001;
    }
}

void changeEditValue(bool add){
    switch(editModeVal){
        case 0:
            changeRadius(add);
            break;
        case 1:
            changeK(add);
            break;
        case 2:
            changeMi(add);
            break;
        case 3:
            changeGradientScale(add);
            break;
        case 4:
            changeParticleRadius(add);
            break;
    }
}

void changeGradientScale(bool add){
    if(add){
        gradientScale += 0.1;
    }
    else{
        gradientScale -= 0.1;
        if(gradientScale < 0.0) gradientScale = 0.0;
    }
}

void toggleExternals(){
    useExternals = !useExternals;
    if(useExternals){
        sprintf(areExternalsOn, "ON");
    }
    else{
        sprintf(areExternalsOn, "OFF");
    }
}

void toggleDistortion(){
    supportDistortion = !supportDistortion;
    if(supportDistortion){
        sprintf(isDistortionOn, "ON");
    }
    else{
        sprintf(isDistortionOn, "OFF");
    }
}

char* printColour(){
    return colourId;
}

void changeRadius(bool add){
    if(add){
        supportRadius += 0.1;
    }
    else{
        supportRadius -= 0.1;
        if(supportRadius < 0.1) supportRadius = 0.1;
    }
    supportRadius3 = supportRadius * supportRadius * supportRadius;
    supportRadius6 = supportRadius3 * supportRadius3;
    supportRadius9 = supportRadius6 * supportRadius3;
}

void setSupRadius(numberType r){
    
    supportRadius = r;
    supportRadius3 = supportRadius * supportRadius * supportRadius;
    supportRadius6 = supportRadius3 * supportRadius3;
    supportRadius9 = supportRadius6 * supportRadius3;
}

float printRadius(){
    return supportRadius;
}

char* printUseExternals(){
    return areExternalsOn;
}

char* printDistortion(){
    return isDistortionOn;
}

char* printEditMode(){
    return editMode;
}



numberType getSRadius(){
    return supportRadius;
}

void changeEditMode(){
    editModeVal++;
    editModeVal = editModeVal % 5;

    switch(editModeVal){
        case 0:
            sprintf(editMode, "Sup. Radius");
            break;
        case 1:
            sprintf(editMode, "Stiffness");
            break;
        case 2:
            sprintf(editMode, "Viscosity");
            break;
        case 3:
            sprintf(editMode, "Grad Scale");
            break;
        case 4:
            sprintf(editMode, "P. Radius");
            break;
    }

}

void gcgSph::insertParticles(){

    int nAdd = 400;
//    nParticles += nAdd;
    particles = (gcgParticleSPH**)realloc(particles, (nParticles+nAdd)*sizeof(gcgParticleSPH*));
    VECTOR3 pos, vel;
    gcgZEROVECTOR3(vel);
    VECTOR3 nPos;
    for(int i = 0; i < nAdd; i++){
        int v = nParticles  + i;
        bool create = false;
         while(!create){
             while(!randomPosition(&pos)){}
             if(grid->hasNeighbours(pos)){create = true;}

         }
//        gcgSETVECTOR3(pos, 15+(int)(i/10), 15+(int)(i%10), 0);
        if(pos[0] >= width -1){pos[0] = width-1.1;}
        if(pos[1] >= height -1){pos[1] = height-1.1;}
        if(tipoCampo != FIELD2D){
                if(pos[2] >= depth -1){pos[2] = depth-1.1;}
        }
        particles[v] = new gcgParticleSPH(particlesMass, pRadius, pos, vel);
        grid->addParticle(particles[v]);
        
//        int addr = ADDRESSGLYPHS(pos[0], pos[1], pos[2], width, height);
//                particles[count]->tensor->setTensor(teste->glyphs_original[v]->tensor);
        if(tipoCampo == FIELD2D)  interpolateTensor(particles[v]);
        else{interpolateTensor3D(particles[v], false, 0);}
        
        /**********SQUAD**********/
                
                
        particles[v]->defineSQuad(particles[v]->tensor);

        /********************/
    }

    nParticles += nAdd;    
    
}

void changeInformationType(){
    informationType++;
    informationType = informationType % 4;

    switch(informationType){
        case 0:
            sprintf(colourId, "Mass-Density");
            break;
        case 1:
            sprintf(colourId, "Acceleration");
            break;
        case 2:
            sprintf(colourId, "Viscosity");
            break;
        case 3:
            sprintf(colourId, "Pressure");
            break;
    }

}

void gcgSph::changeColourType(){
    changeInformationType();
    maxDebugVal = -99999999;
    minDebugVal = 99999999;

    for(int i = 0; i < nParticles; i++){
        numberType tmpVal;
        gcgParticleSPH * p = particles[i];
        switch(informationType){
            case 0:{
                tmpVal = p->lastMassDensity;
                break;
            }
            case 1:{
                tmpVal = gcgLENGTHVECTOR3(p->lastAcceleration);
                break;
            }
            case 2:{
                tmpVal = gcgLENGTHVECTOR3(p->lastViscosity);
                break;
            }
            case 3:{
                tmpVal = gcgLENGTHVECTOR3(p->lastPressure);
                break;
            }

        }

        if(tmpVal >= maxDebugVal){
            maxDebugVal = tmpVal;
        }
        if(tmpVal <= minDebugVal){
            minDebugVal = tmpVal;
        }
    }

    glutPostRedisplay();
}


//class RandomQueue{
//public:
//    unsigned int count = 0;
//    VECTOR3 ** points;
//    RandomQueue(int totalPoints){
//        count = totalPoints;
//        points = (VECTOR3**) malloc(count*sizeof(VECTOR3*));
//        for(int i = 0; i < count; i++){
//            points[i] = (VECTOR3*) malloc(sizeof(VECTOR3));
//            gcgSETVECTOR3(points[i], 0.,0.,0.);
//        }
//    }
//    ~RandomQueue(){
//        for(int i = 0; i < count; i++){
//            delete(points[i]);
//        }
//        delete(points);
//    }
//    void push(VECTOR3 point){
//        gcgCOPYVECTOR3(points[count], point);
//        count++;
//    }
//    void pop(VECTOR3 point){
//        int countMinus1 = (count - 1);
//        int index = rand() % (countMinus1);
//        gcgCOPYVECTOR3(point, points[index]);
//        if(index != (countMinus1)){
//            gcgCOPYVECTOR3(points[index], points[countMinus1]);
//        }
//        gcgSETVECTOR3(points[countMinus1], 0.,0.,0.);
//        count--;
//    }
//};


//void generate_poisson(width, height, depth, min_dist, new_points_count)
//{
//  //Create the grid
//  cellSize = min_dist/sqrt(2);
//
//  grid = Grid3D(Point(
//    (ceil(width/cell_size),         //grid width
//     ceil(height/cell_size),
//     ceil(depth/cell_size))));      //grid height
//
//  //RandomQueue works like a queue, except that it
//  //pops a random element from the queue instead of
//  //the element at the head of the queue
//  RandomQueue * processList = new RandomQueue(new_points_count);
//
//  samplePoints = List();
//
//  //generate the first point randomly
//  //and updates
//
//  firstPoint = Point(rand(width), rand(height));
//
//  //update containers
//  processList.push(firstPoint);
//  samplePoints.push(firstPoint);
//  grid[imageToGrid(firstPoint, cellSize)] = firstPoint;
//
//  //generate other points from points in queue.
//  while (not processList.empty())
//  {
//    point = processList.pop();
//    for (i = 0; i < new_points_count; i++)
//    {
//      newPoint = generateRandomPointAround(point, min_dist);
//      //check that the point is in the image region
//      //and no points exists in the point's neighbourhood
//      if (inRectangle(newPoint) and
//        not inNeighbourhood(grid, newPoint, min_dist,
//          cellSize))
//      {
//        //update containers
//        processList.push(newPoint);
//        samplePoints.push(newPoint);
//        grid[imageToGrid(newPoint, cellSize)] =  newPoint;
//      }
//    }
//  }
//  return samplePoints;
//}

/***********    DEBUG     *************/

void printDebugVals(){
//    printf("Maior aceleração: %f | Vetor: %f %f %f\n", maxAcc, maxVec[0], maxVec[1], maxVec[2]);
//    printf("Menor aceleração: %f | Vetor: %f %f %f\n", minAcc, minVec[0], minVec[1], minVec[2]);
}


/***********    DEBUG     *************/

bool gcgSph::randomPosition(VECTOR3 * position){

    float r = fabs(gcgRandom.random());
    while (r > 1.0 || r < 0.0) r = fabs(gcgRandom.random());
//    int   n = (int) ((float) teste->ordSize * (r / 2.0));///2.0 no lugar do alpha
    int   n = (int) ((float) teste->ordSize * (r));///2.0 no lugar do alpha

    VECTOR3 pos;
    int posX, posY, posZ;

//        n = i;


    pos[0] = teste->tensor_ord[n]->pos[0];
    pos[1] = teste->tensor_ord[n]->pos[1];
    pos[2] = teste->tensor_ord[n]->pos[2];
//        tensor_ord[n]->numberOfVisits++;
    posX = (int) pos[0];
    posY = (int) pos[1];
    posZ = (int) pos[2];
    int p = ADDRESSVOLUME(posX, posY, posZ, width, height);


    //Restrições de criação
    //        if ((p < totalSize) && (p >= 0) && (tensor_ord[p]->otm == true) && (tensor_ord[p]->cs <= 0.01)) {
    //        if ((p < totalSize) && (p >= 0) && (tensor_ord[p]->cs <= 0.01)) {
    if ((p < (width*height*depth)) && (p >= 0)) {
        gcgCOPYVECTOR3(*position, pos);
        return true;
    }
    else return false;

}

int signal(float x){
    if(x < 0) return -1;
    if(x == 0) return 0;
    return 1;
}

float getDistance(VECTOR3 p1, VECTOR3 p2){
    VECTOR3 differenceVector;
    gcgSETVECTOR3(differenceVector, (p1[0]-p2[0]), (p1[1]-p2[1]), (p1[2]-p2[2]));
    return (sqrtf((differenceVector[0]*differenceVector[0])+(differenceVector[1]*differenceVector[1])+(differenceVector[2]*differenceVector[2])));
}

void gcgSph::toggleBox(){
    usedVolume = !usedVolume;
    if(usedVolume){
        this->boxWidth = width;
        this->boxHeight = height;
        this->boxDepth = depth;
    }
    else{
        this->boxWidth = widthFluid;
        this->boxHeight = heightFluid;
        this->boxDepth = depthFluid;
    }
}

void gcgSph::resetParticles(){

    for(int i = 0; i < nParticles; i++){
         while(!randomPosition(&particles[i]->position)){

         }
         particles[i]->resetIteration();
         gcgZEROVECTOR3(particles[i]->velocity);
         gcgZEROVECTOR3(particles[i]->acceleration);
        #if (INTEGRATOR == LEAP_FROG)
            gcgSCALEVECTOR3(particles[i]->lastVelocity, particles[i]->acceleration, (-0.5*TIME_STEP));
            gcgADDVECTOR3(particles[i]->lastVelocity, particles[i]->lastVelocity, particles[i]->velocity);
            grid->moveParticle(particles[i]);
        #endif
    }

}

void gcgSph::interpolateTensor(gcgParticleSPH * p){


    float x1, x2, y1, y2, w, h, dx1, dx2, dy1, dy2;
//    float norm11, norm12, norm21, norm22;
    VECTOR3 norm11, norm12, norm21, norm22;

    x1 = floor(p->position[0]);

    if(x1 >= width) x1 -= 1;

    x2 = x1+1.;

    if(x2 >= width) x2 -= 1;

    y1 = floor(p->position[1]);

    if(y1 >= height) y1 -= 1;

    y2 = y1+1.;

    if(y2 >= height) y2 -= 1;

    w = x2 - x1;
    h = y2 - y1;

    dx1 = (float)(x2 - p->position[0]) / w;
    dx2 = 1.0 - dx1;
//    dx2 = (float) (p->position[0] - x1) / w;

    MATRIX3 R1, R2, P, temp1, temp2;

    //    R1 = ((x2 – x)/(x2 – x1))*Q11 + ((x – x1)/(x2 – x1))*Q21
    int glyphPos = ADDRESSGLYPHS(x1, y1, 0, width, height);
    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, dx1);
    gcgCOPYVECTOR3(norm11, teste->normDerivatives[glyphPos]);

    glyphPos = ADDRESSGLYPHS(x2, y1, 0, width, height);
    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, dx2);

    gcgADDMATRIX3(R1, temp1, temp2);
    gcgCOPYVECTOR3(norm21, teste->normDerivatives[glyphPos]);


    //    R2 = ((x2 – x)/(x2 – x1))*Q12 + ((x – x1)/(x2 – x1))*Q22
    glyphPos = ADDRESSGLYPHS(x1, y2, 0, width, height);
    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, dx1);
    gcgCOPYVECTOR3(norm12, teste->normDerivatives[glyphPos]);

    glyphPos = ADDRESSGLYPHS(x2, y2, 0, width, height);
    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, dx2);
    gcgADDMATRIX3(R2, temp1, temp2);

//    norm22 = teste->norms[glyphPos];
    gcgCOPYVECTOR3(norm22, teste->normDerivatives[glyphPos]);

    //    P = ((y2 – y)/(y2 – y1))*R1 + ((y – y1)/(y2 – y1))*R2
    dy1 = (y2 - p->position[1])/h;
    dy2 = (p->position[1] - y1)/h;

    gcgSCALEMATRIX3(temp1, R1, dy1);
    gcgSCALEMATRIX3(temp2, R2, dy2);

    gcgADDMATRIX3(P, temp1, temp2);

    p->tensor->setTensor(P);

    p->gradientVal[0] = (dy1*((norm11[0]*dx1)+(norm21[0]*dx2)))+(dy2*((norm12[0]*dx1)+(norm22[0]*dx2)));
    p->gradientVal[1] = (dy1*((norm11[1]*dx1)+(norm21[1]*dx2)))+(dy2*((norm12[1]*dx1)+(norm22[1]*dx2)));
    p->gradientVal[2] = 0;

    gcgSCALEVECTOR3(p->gradientVal, p->gradientVal, gradientScale);

    if(activateBoost){
        if(p->tensor->cs < 0.8){
            gcgSCALEVECTOR3(p->gradientVal, p->gradientVal, 1.0 + (1.0 - p->tensor->cs));
        }
    }

//    printf("%p - %f %f %f\n", p, p->gradientVal[0], p->gradientVal[1], p->gradientVal[2]);



}

void gcgSph::interpolateTensor3D(gcgParticleSPH * p, bool debug, int i){


    float xd, xd1, yd, yd1, zd, zd1, x1, x2, y1, y2, z1, z2, w, h, d, dx1, dx2, dy1, dy2;

//    float norm11, norm12, norm21, norm22;

    x1 = floor(p->position[0]);
    if(x1 >= width) x1 = width - 2;
//    if(x1 >= width) x1 -= 1;

    x2 = x1+1.;
    if(x2 >= width) x2 -= 1;

    y1 = floor(p->position[1]);
    if(y1 >= height) y1 = height - 2;
//    if(y1 >= height) y1 -= 1;

    y2 = y1+1.;
    if(y2 >= height) y2 -= 1;

    z1 = floor(p->position[2]);
    if(z1 >= depth) z1 = depth -2;
//    if(z1 >= depth) z1 -= 1;

    z2 = z1+1.;
    if(z2 >= depth) z2 -= 1;


    w = x2 - x1;
    h = y2 - y1;
    d = z2 - z1;

    xd = (p->position[0] - x1)/w;
    xd1 = 1-xd;
    yd = (p->position[1] - y1)/h;
    yd1 = 1-yd;
    zd = (p->position[2] - z1)/d;
    zd1 = 1-zd;

    ///Interpolate x
    MATRIX3 P, temp1, temp2;
    int glyphPos;
    MATRIX3 c00, c01, c11, c10;
    VECTOR3 vc00, vc01, vc10, vc11, vP, vtemp1, vtemp2;

    glyphPos = ADDRESSGLYPHS(x1, y1, z1, width, height);
    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, xd1);
    gcgSCALEVECTOR3(vtemp1, teste->normDerivatives[glyphPos], xd1);


    glyphPos = ADDRESSGLYPHS(x2, y1, z1, width, height);
    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, xd);
    gcgSCALEVECTOR3(vtemp2, teste->normDerivatives[glyphPos], xd);

    gcgADDMATRIX3(c00, temp1, temp2);
    gcgADDVECTOR3(vc00, vtemp1, vtemp2);


    glyphPos = ADDRESSGLYPHS(x1, y2, z1, width, height);
    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, xd1);
    gcgSCALEVECTOR3(vtemp1, teste->normDerivatives[glyphPos], xd1);

    glyphPos = ADDRESSGLYPHS(x2, y2, z1, width, height);
    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, xd);
    gcgSCALEVECTOR3(vtemp2, teste->normDerivatives[glyphPos], xd);

    gcgADDMATRIX3(c10, temp1, temp2);
    gcgADDVECTOR3(vc10, vtemp1, vtemp2);


    glyphPos = ADDRESSGLYPHS(x1, y1, z2, width, height);
    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, xd1);
    gcgSCALEVECTOR3(vtemp1, teste->normDerivatives[glyphPos], xd1);

    glyphPos = ADDRESSGLYPHS(x2, y1, z2, width, height);
    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, xd);
    gcgSCALEVECTOR3(vtemp2, teste->normDerivatives[glyphPos], xd);

    gcgADDMATRIX3(c01, temp1, temp2);
    gcgADDVECTOR3(vc01, vtemp1, vtemp2);


    glyphPos = ADDRESSGLYPHS(x1, y2, z2, width, height);
    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, xd1);
    gcgSCALEVECTOR3(vtemp1, teste->normDerivatives[glyphPos], xd1);

    glyphPos = ADDRESSGLYPHS(x2, y2, z2, width, height);
    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, xd);
    gcgSCALEVECTOR3(vtemp2, teste->normDerivatives[glyphPos], xd);

    gcgADDMATRIX3(c11, temp1, temp2);
    gcgADDVECTOR3(vc11, vtemp1, vtemp2);


    ///Interpolate y
    MATRIX3 c0, c1;
    VECTOR3 vc0, vc1;

    gcgSCALEMATRIX3(temp1, c00, yd1);
    gcgSCALEMATRIX3(temp2, c10, yd);
    gcgADDMATRIX3(c0, temp1, temp2);

    gcgSCALEMATRIX3(temp1, c01, yd1);
    gcgSCALEMATRIX3(temp2, c11, yd);
    gcgADDMATRIX3(c1, temp1, temp2);

    gcgSCALEVECTOR3(vtemp1, vc00, yd1);
    gcgSCALEVECTOR3(vtemp2, vc10, yd);
    gcgADDVECTOR3(vc0, vtemp1, vtemp2);

    gcgSCALEVECTOR3(vtemp1, vc01, yd1);
    gcgSCALEVECTOR3(vtemp2, vc11, yd);
    gcgADDVECTOR3(vc1, vtemp1, vtemp2);

    ///Interpolate z
    gcgSCALEMATRIX3(temp1, c0, zd1);
    gcgSCALEMATRIX3(temp2, c1, zd);
    gcgADDMATRIX3(P, temp1, temp2);

    gcgSCALEVECTOR3(vtemp1, vc0, zd1);
    gcgSCALEVECTOR3(vtemp2, vc1, zd);
    gcgADDVECTOR3(vP, vtemp1, vtemp2);


    p->tensor->setTensor(P);
    if(debug){
        printf("Interpolado: %f %f %f %f %f %f %f %f %f\n", P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
        if(i == 499){
            printf("Bounds: %f %f %f %f %f %f\n", x1, x2, y1, y2, z1, z2);
        }
    }
    gcgCOPYVECTOR3(p->gradientVal, vP);
    gcgSCALEVECTOR3(p->gradientVal, p->gradientVal, gradientScale);
//    printf("%p - %f %f %f\n", p, p->gradientVal[0], p->gradientVal[1], p->gradientVal[2]);



}

//void gcgSph::interpolateTensor(gcgParticleSPH * p){
//
//
//    float x1, x2, y1, y2, w, h, dx1, dx2;
//
//    x1 = floor(p->position[0]);
//
//    if(x1 >= width) x1 -= 1;
//
//    x2 = x1+1.;
//
//    if(x2 >= width) x2 -= 1;
//
//    y1 = floor(p->position[1]);
//
//    if(y1 >= height) y1 -= 1;
//
//    y2 = y1+1.;
//
//    if(y2 >= height) y2 -= 1;
//
//    w = x2 - x1;
//    h = y2 - y1;
//
//    dx1 = (float)(x2 - p->position[0]) / w;
//    dx2 = 1.0 - dx1;
////    dx2 = (float) (p->position[0] - x1) / w;
//
//    MATRIX3 R1, R2, P, temp1, temp2;
//
//    //    R1 = ((x2 – x)/(x2 – x1))*Q11 + ((x – x1)/(x2 – x1))*Q21
//    int glyphPos = ADDRESSGLYPHS(x1, y1, 0, width, height);
//    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, dx1);
//
//    glyphPos = ADDRESSGLYPHS(x2, y1, 0, width, height);
//    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, dx2);
//
//    gcgADDMATRIX3(R1, temp1, temp2);
//
//    //    R2 = ((x2 – x)/(x2 – x1))*Q12 + ((x – x1)/(x2 – x1))*Q22
//    glyphPos = ADDRESSGLYPHS(x1, y2, 0, width, height);
//    gcgSCALEMATRIX3(temp1, teste->glyphs_original[glyphPos]->tensor, dx1);
//
//    glyphPos = ADDRESSGLYPHS(x2, y2, 0, width, height);
//    gcgSCALEMATRIX3(temp2, teste->glyphs_original[glyphPos]->tensor, dx2);
//
//    gcgADDMATRIX3(R2, temp1, temp2);
//
//    //    P = ((y2 – y)/(y2 – y1))*R1 + ((y – y1)/(y2 – y1))*R2
//
//    gcgSCALEMATRIX3(temp1, R1, (y2 - p->position[1])/h);
//    gcgSCALEMATRIX3(temp2, R2, (p->position[1] - y1)/h);
//
//    gcgADDMATRIX3(P, temp1, temp2);
//
//    p->tensor->setTensor(P);
//
//
//}


gcgSph::gcgSph(int nParticles, unsigned int fluidWidth, unsigned int fluidHeight, unsigned int fluidDepth, unsigned int width, unsigned int height, unsigned int depth) {


    defineConstants();

    srand ( time(NULL) );

    int minX = 10, maxX = 20;
    srand ( time(NULL) );
    usedVolume = 0;
    simulationTime = 0.0;
  /* generate secret number: */
    this->width = width;
    this->height = height;
    this->depth = depth;
    this->boxWidth = fluidWidth;
    this->boxHeight = fluidHeight;
    this->boxDepth = fluidDepth;
    this->widthFluid = fluidWidth;
    this->heightFluid = fluidHeight;
    this->depthFluid = fluidDepth;

    translateW = -(int)(width / 2);
    translateH = -(int)(height / 2);
    translateD = -(int)(depth / 2);



    this->nParticles = nParticles;
    this->initialDensity = restDensity;
    VECTOR3 velocity, pos;
    gcgSETVECTOR3(velocity, 0.,0.,0.);
    VECTOR3 c1, c2;
    gcgSETVECTOR3(c1,0,0,0);
    gcgSETVECTOR3(c2,width,height,depth);
//    octree = new Octree(c1, c2, 0, NULL);
    unsigned int volume = fluidDepth * fluidHeight * fluidWidth;

    direction = 1.;

//    supportRadius = 0.85*pow((3.*(numberType)volume*(numberType)ESTIMATED_NEIGHBOURS/(4.*M_PI*nParticles)), (1./3.));
//    supportRadius = pow((3.*(numberType)volume*(numberType)ESTIMATED_NEIGHBOURS/(4.*M_PI*nParticles)), (1./3.));
    supportRadius = radiusScale*2*pow((3.*(numberType)volume*(numberType)ESTIMATED_NEIGHBOURS/(4.*M_PI*nParticles)), (1./3.));
//    supportRadius = 1.893664;
    supportRadius3 = supportRadius * supportRadius * supportRadius;
    supportRadius6 = supportRadius3 * supportRadius3;
    supportRadius9 = supportRadius6 * supportRadius3;
    influenceRadius = influenceFactor * supportRadius;

    if(tipoCampo == FIELD2D){
    grid = new newGrid(width, height, depth, supportRadius, true);
    }
    else{grid = new newGrid(width, height, depth, supportRadius, false);}
//    this->nParticles = (int)((60 * volume)/(4*M_PI*supportRadius3));

    printf("Particulas: %d\n", this->nParticles);

//    this->width = width;
//    this->height = height;
//    this->depth = depth;

    numberType mass = (numberType)((numberType)restDensity * (numberType)volume) / ((numberType)nParticles);
    this->particlesMass = mass;
    printf("Densidade: %f | Volume: %f | Massa: %f\n",restDensity, (numberType)volume, mass);

    numberType cres = 0.1;
//    numberType cres = 1;
    numberType x = 5;
    int count = 0;
    particles = (gcgParticleSPH**)malloc(nParticles*sizeof(gcgParticleSPH*));



    pRadius = radiusScale * pow(0.75*mass/(restDensity * M_PI), 1.0/3.0);
//    pRadius = 0.4136;
    particlesRadius = pRadius;
//    numberType pRadius = 0.5;
    printf("Raio: %f\n", pRadius);

//    int inix = 20, iniy = 15, iniz = 13;
//    int inix = 17, iniy = 15, iniz = 16; //ate27k
//    int inix = 20, iniy = 15, iniz = 13;
    int inix = 1, iniy = 1, iniz = 1;
//    int inix = 1, iniy = 1, iniz = 1;
//    int inix = width/2 - (fluidWidth/2), iniy = height/2 - (fluidHeight/2), iniz = depth/2  - (fluidDepth/2);
//    int inix = width/2 , iniy = height/2 - 5, iniz = depth/2 - 5;
    if(tipoCampo == FIELD2D) iniz = 0;
//    else iniz = 20.;
//    int inix = 31, iniy = 31;
    
//    cres = nParticles / ((float)volume);
    
//    while(count < nParticles){
//        
//        float i = count % fluidWidth;
//        float j = ((int)(floor(count / (fluidWidth)))) % fluidHeight;
//        float k = ((int)(floor(count / (fluidWidth*fluidHeight)))) % fluidDepth;
//        gcgSETVECTOR3(pos, inix+(i*cres), iniy+(j*cres), iniz+(k*cres));
//        
//        particles[count] = new gcgParticleSPH(mass, pRadius, pos, velocity);
//        printf("%f %f %f - %d\n", pos[0], pos[1], pos[2], count);
//                grid->addParticle(particles[count]);
//
////                particles[count]->defineParticle(mass, pRadius, pos, velocity);
////                particles[count].kernelDerivative(NULL);
////                octree->add(particles[count]);
//                gcgSETVECTOR3(particles[count]->color, 0.,0.8,0.0);
//
//                int v = ADDRESSGLYPHS(pos[0], pos[1], pos[2], width, height);
////                particles[count]->tensor->setTensor(teste->glyphs_original[v]->tensor);
//                if(tipoCampo == FIELD2D)  interpolateTensor(particles[count]);
//                else{interpolateTensor3D(particles[count], false, 0);}
////                interpolateTensor3D(particles[count]);
//
//
//
//
//
//                #if (INTEGRATOR == LEAP_FROG)
////                    gcgSCALEVECTOR3(particles[count]->nextVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
////                    gcgADDVECTOR3(particles[count]->nextVelocity, particles[count]->nextVelocity, particles[count]->velocity);
//                    gcgSCALEVECTOR3(particles[count]->lastVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
//                    gcgADDVECTOR3(particles[count]->lastVelocity, particles[count]->lastVelocity, particles[count]->velocity);
//                #endif
//
//                count++;
//    }

//    FILE * dbg = fopen("dbg.txt", "w");
    for(int i = 0; i < fluidWidth; i++){
//    for(int i = 0; i < fluidWidth; i++){
        if(count  >= nParticles) break;
        for(int j = 0; j < fluidHeight; j++){
            if(count  >= nParticles) break;
            for(int k = 0; k < nParticles; k++){
//            for(int i = 0; i < nParticles; i++){
                if(count  >= nParticles) break;
//                printf("count: %d\n", i);
                gcgSETVECTOR3(pos, inix+(i*cres), iniy+(j*cres), iniz+(k*cres));
                gcgSETVECTOR3(pos, 1,1,1);



                particles[count] = new gcgParticleSPH(mass, pRadius, pos, velocity);
//                fprintf(dbg, "%f %f %f\n", pos[0], pos[1], pos[2]);
                grid->addParticle(particles[count]);

//                particles[count]->defineParticle(mass, pRadius, pos, velocity);
//                particles[count].kernelDerivative(NULL);
//                octree->add(particles[count]);
                gcgSETVECTOR3(particles[count]->color, 0.,0.8,0.0);

                int v = ADDRESSGLYPHS(pos[0], pos[1], pos[2], width, height);
//                particles[count]->tensor->setTensor(teste->glyphs_original[v]->tensor);
                if(tipoCampo == FIELD2D)  interpolateTensor(particles[count]);
                else{interpolateTensor3D(particles[count], false, 0);}
                
//                particles[count]->tensor->setTensor(teste->glyphs_original[v]->tensor);
                /**********SQUAD**********/
                
                
                particles[count]->defineSQuad(particles[count]->tensor);
                
                /********************/
//                interpolateTensor3D(particles[count]);





                #if (INTEGRATOR == LEAP_FROG)
//                    gcgSCALEVECTOR3(particles[count]->nextVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
//                    gcgADDVECTOR3(particles[count]->nextVelocity, particles[count]->nextVelocity, particles[count]->velocity);
                    gcgSCALEVECTOR3(particles[count]->lastVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
                    gcgADDVECTOR3(particles[count]->lastVelocity, particles[count]->lastVelocity, particles[count]->velocity);
                #endif

                count++;

            }
        }
    }
//    fclose(dbg);
    
    printf("Alocou sph - Count = %d\n", count);
    if(useExternals){
        sprintf(areExternalsOn, "ON");
    }
    else{
        sprintf(areExternalsOn, "OFF");
    }
    if(supportDistortion){
        sprintf(isDistortionOn, "ON");
    }
    else{
        sprintf(isDistortionOn, "OFF");
    }
    changeEditMode();

}


gcgSph::gcgSph(int pX, int pY, int pZ, unsigned int fluidWidth, unsigned int fluidHeight, unsigned int fluidDepth, unsigned int width, unsigned int height, unsigned int depth) {

    int nParticles = pX*pY*pZ;
    defineConstants();

    srand ( time(NULL) );

    int minX = 10, maxX = 20;
    srand ( time(NULL) );
    usedVolume = 0;
    simulationTime = 0.0;
  /* generate secret number: */
    this->width = width;
    this->height = height;
    this->depth = depth;
    this->boxWidth = fluidWidth;
    this->boxHeight = fluidHeight;
    this->boxDepth = fluidDepth;
    this->widthFluid = fluidWidth;
    this->heightFluid = fluidHeight;
    this->depthFluid = fluidDepth;

    translateW = -(int)(width / 2);
    translateH = -(int)(height / 2);
    translateD = -(int)(depth / 2);



    this->nParticles = nParticles;
    this->initialDensity = restDensity;
    VECTOR3 velocity, pos;
    gcgSETVECTOR3(velocity, 0.,0.,0.);
    VECTOR3 c1, c2;
    gcgSETVECTOR3(c1,0,0,0);
    gcgSETVECTOR3(c2,width,height,depth);
//    octree = new Octree(c1, c2, 0, NULL);
    unsigned int volume = fluidDepth * fluidHeight * fluidWidth;

    direction = 1.;

//    supportRadius = 0.85*pow((3.*(numberType)volume*(numberType)ESTIMATED_NEIGHBOURS/(4.*M_PI*nParticles)), (1./3.));
//    supportRadius = pow((3.*(numberType)volume*(numberType)ESTIMATED_NEIGHBOURS/(4.*M_PI*nParticles)), (1./3.));
    supportRadius = radiusScale*2*pow((3.*(numberType)volume*(numberType)ESTIMATED_NEIGHBOURS/(4.*M_PI*nParticles)), (1./3.));
    supportRadius3 = supportRadius * supportRadius * supportRadius;
    supportRadius6 = supportRadius3 * supportRadius3;
    supportRadius9 = supportRadius6 * supportRadius3;
    influenceRadius = influenceFactor * supportRadius;

    if(tipoCampo == FIELD2D){
    grid = new newGrid(width, height, depth, supportRadius, true);
    }
    else{grid = new newGrid(width, height, depth, supportRadius, false);}
//    this->nParticles = (int)((60 * volume)/(4*M_PI*supportRadius3));

    printf("Particulas: %d\n", this->nParticles);

//    this->width = width;
//    this->height = height;
//    this->depth = depth;

    numberType mass = (numberType)((numberType)restDensity * (numberType)volume) / ((numberType)nParticles);
    this->particlesMass = mass;
    printf("Densidade: %f | Volume: %f | Massa: %f\n",restDensity, (numberType)volume, mass);

    numberType cres = 0.5;
//    numberType cres = 1;
    numberType x = 5;
    int count = 0;
    particles = (gcgParticleSPH**)malloc(nParticles*sizeof(gcgParticleSPH*));



    pRadius = radiusScale * pow(0.75*mass/(restDensity * M_PI), 1.0/3.0);
    particlesRadius = pRadius;
//    numberType pRadius = 0.5;
    printf("Raio: %f\n", pRadius);

    int inix = 1, iniy = 1, iniz = 1; //2d
//    int inix = 25, iniy = 15, iniz = 16; //ate5k
//    int inix = 25, iniy = 15, iniz = 16; //ate9k
//    int inix = 17, iniy = 15, iniz = 16; //ate27k
//    int inix = 10, iniy = 15, iniz = 13; //ate27k
//    int inix = 8, iniy = 15, iniz = 13; //60k
//    int inix = 40, iniy = 55, iniz = 45; //brin
//    int inix = 32, iniy = 17, iniz = 0;
//    int inix = 1, iniy = 1, iniz = 1;
//    int inix = 0, iniy = 0, iniz = 0;
//    int inix = width/2 - (fluidWidth/2), iniy = height/2 - (fluidHeight/2), iniz = depth/2  - (fluidDepth/2);
//    int inix = width/2 - (fluidWidth/2), iniy = height/2 - (fluidHeight/2) -11, iniz = depth/2  - (fluidDepth/2);
//    float inix = width/2 - (fluidWidth/2) + 0.5, iniy = height/2 - (fluidHeight/2) - 8.5, iniz = depth/2  - (fluidDepth/2);
//    int inix = width/2 , iniy = height/2 - 5, iniz = depth/2 - 5;
    if(tipoCampo == FIELD2D) iniz = 0;
//    else iniz = 20.;
//    int inix = 31, iniy = 31;
    
    float xCres = ((float)fluidWidth)/((float)pX);
    float yCres = ((float)fluidHeight)/((float)pY);
    float zCres = ((float)fluidDepth)/((float)pZ);
    
    while(count < nParticles){
        
        float i = count % fluidWidth;
        float j = ((int)(floor(count / (fluidWidth)))) % fluidHeight;
        float k = ((int)(floor(count / (fluidWidth*fluidHeight)))) % fluidDepth;
        gcgSETVECTOR3(pos, inix+(i*xCres), iniy+(j*yCres), iniz+(k*zCres));
        
        particles[count] = new gcgParticleSPH(mass, pRadius, pos, velocity);
//        printf("%f %f %f - %d\n", pos[0], pos[1], pos[2], count);
                grid->addParticle(particles[count]);

//                particles[count]->defineParticle(mass, pRadius, pos, velocity);
//                particles[count].kernelDerivative(NULL);
//                octree->add(particles[count]);
                gcgSETVECTOR3(particles[count]->color, 0.,0.8,0.0);

                int v = ADDRESSGLYPHS(pos[0], pos[1], pos[2], width, height);
//                particles[count]->tensor->setTensor(teste->glyphs_original[v]->tensor);
                if(tipoCampo == FIELD2D)  interpolateTensor(particles[count]);
                else{interpolateTensor3D(particles[count], false, 0);}
//                interpolateTensor3D(particles[count]);


                /**********SQUAD**********/
                
                
                particles[count]->defineSQuad(particles[count]->tensor);
                
                /********************/


                #if (INTEGRATOR == LEAP_FROG)
//                    gcgSCALEVECTOR3(particles[count]->nextVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
//                    gcgADDVECTOR3(particles[count]->nextVelocity, particles[count]->nextVelocity, particles[count]->velocity);
                    gcgSCALEVECTOR3(particles[count]->lastVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
                    gcgADDVECTOR3(particles[count]->lastVelocity, particles[count]->lastVelocity, particles[count]->velocity);
                #endif

                count++;
    }

//    FILE * dbg = fopen("dbg.txt", "w");
//    for(int i = 0; i < fluidWidth; i++){
//        if(count  >= nParticles) break;
//        for(int j = 0; j < fluidHeight; j++){
//            if(count  >= nParticles) break;
//            for(int k = 0; k < fluidDepth; k++){
////            for(int i = 0; i < nParticles; i++){
//                if(count  >= nParticles) break;
////                printf("count: %d\n", i);
//                gcgSETVECTOR3(pos, inix+(i*cres), iniy+(j*cres), iniz+(k*cres));
//
//
//
//                particles[count] = new gcgParticleSPH(mass, pRadius, pos, velocity);
//                fprintf(dbg, "%f %f %f\n", pos[0], pos[1], pos[2]);
//                grid->addParticle(particles[count]);
//
////                particles[count]->defineParticle(mass, pRadius, pos, velocity);
////                particles[count].kernelDerivative(NULL);
////                octree->add(particles[count]);
//                gcgSETVECTOR3(particles[count]->color, 0.,0.8,0.0);
//
//                int v = ADDRESSGLYPHS(pos[0], pos[1], pos[2], width, height);
////                particles[count]->tensor->setTensor(teste->glyphs_original[v]->tensor);
//                if(tipoCampo == FIELD2D)  interpolateTensor(particles[count]);
//                else{interpolateTensor3D(particles[count], false, 0);}
////                interpolateTensor3D(particles[count]);
//
//
//
//
//
//                #if (INTEGRATOR == LEAP_FROG)
////                    gcgSCALEVECTOR3(particles[count]->nextVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
////                    gcgADDVECTOR3(particles[count]->nextVelocity, particles[count]->nextVelocity, particles[count]->velocity);
//                    gcgSCALEVECTOR3(particles[count]->lastVelocity, particles[count]->acceleration, (-0.5*TIME_STEP));
//                    gcgADDVECTOR3(particles[count]->lastVelocity, particles[count]->lastVelocity, particles[count]->velocity);
//                #endif
//
//                count++;
//
//            }
//        }
//    }
//    fclose(dbg);
    
    printf("Alocou sph - Count = %d\n", count);
    if(useExternals){
        sprintf(areExternalsOn, "ON");
    }
    else{
        sprintf(areExternalsOn, "OFF");
    }
    if(supportDistortion){
        sprintf(isDistortionOn, "ON");
    }
    else{
        sprintf(isDistortionOn, "OFF");
    }
    changeEditMode();

}


void gcgSph::checkBoundaryCondition(gcgParticleSPH * p, VECTOR3 at, VECTOR3 vt){
    ///Colocar restricao


    VECTOR3 predpos;

    VECTOR3 aux;

    /******************/

//    gcgSCALEVECTOR3(p->tensorForce, p->tensorForce, 1./gcgLENGTHVECTOR3(p->tensorForce));
//    gcgSCALEVECTOR3(p->tensorForce, p->tensorForce, TIME_STEP);
//    gcgADDVECTOR3(p->position, p->position, p->tensorForce);
//    return;
    /*******************/

  #if (INTEGRATOR == VERLET)
    gcgADDVECTOR3(aux, vt, at);
  #else
        gcgSCALEVECTOR3(aux, p->nextVelocity, TIME_STEP);
//    gcgADDVECTOR3(aux, p->nextVelocity, p->lastVelocity);
//    gcgSCALEVECTOR3(aux, aux, 0.5);
  #endif

  gcgADDVECTOR3(predpos, p->position, aux);


  ///if the particle is outside of the bounding box, we invert its velocity
    if (predpos[0] <= 0){
       #if (INTEGRATOR == VERLET)
        p->velocity[0] = ENERGY*fabs(p->velocity[0]);
       #else
        p->nextVelocity[0] = ENERGY*fabs(p->nextVelocity[0]);
       #endif
     }

     if (predpos[0] >= boxWidth){
       #if (INTEGRATOR == VERLET)
        p->velocity[0] = -ENERGY*fabs(p->velocity[0]);
       #else
        p->nextVelocity[0] = -ENERGY*fabs(p->nextVelocity[0]);
       #endif
     }

     if ((predpos[1] <= 0)){// || (predpos[1] > height)){
        ///we should check this better in the future, because perhaps this can lead to strange results if the particle goes too far from the box
       #if (INTEGRATOR == VERLET)
        p->velocity[1] = ENERGY*fabs(p->velocity[1]);
       #else
        p->nextVelocity[1] = ENERGY*fabs(p->nextVelocity[1]);
       #endif
     }
     if ((predpos[1] >= boxHeight)){// || (predpos[1] > height)){
        ///we should check this better in the future, because perhaps this can lead to strange results if the particle goes too far from the box
       #if (INTEGRATOR == VERLET)
        p->velocity[1] = -ENERGY*fabs(p->velocity[1]);
       #else
        p->nextVelocity[1] = -ENERGY*fabs(p->nextVelocity[1]);
       #endif
     }


     if (predpos[2] <= 0){
       #if (INTEGRATOR == VERLET)
        p->velocity[2] = ENERGY*fabs(p->velocity[2]);
       #else
        p->nextVelocity[2] = ENERGY*fabs(p->nextVelocity[2]);
       #endif
     }

     if (predpos[2] >= boxDepth){
       #if (INTEGRATOR == VERLET)
        p->velocity[2] = -ENERGY*fabs(p->velocity[2]);
       #else
        p->nextVelocity[2] = -ENERGY*fabs(p->nextVelocity[2]);
       #endif
     }
  ///p->pos = p->pos + timeStep * p->vel; it's the new vel, as in the book. is it right?



  #if (INTEGRATOR == VERLET)
    gcgADDVECTOR3(aux, vt, at);
  #else
    gcgSCALEVECTOR3(aux, p->nextVelocity, TIME_STEP);
//        gcgADDVECTOR3(aux, p->nextVelocity, p->lastVelocity);
//        gcgSCALEVECTOR3(aux, aux, 0.5);
  #endif


//  gcgCOPYVECTOR3(predpos, p->position);

  gcgADDVECTOR3(p->position, p->position, aux);

  if(p->position[0] < 0) {p->position[0] = 0;}
  else if(p->position[0] >= width - 1){ p->position[0] = width -1.1;}

  if(p->position[1] < 0) {p->position[1] = 0;}
  else if(p->position[1] >= height - 1){ p->position[1] = height -1.1;}

  if(p->position[2] < 0) {p->position[2] = 0;}
  else{
      if(depth > 1 && p->position[2] >= depth - 1){ p->position[2] = depth -1.11;}
  }




//  octree->particleMoved(p, predpos);

//  p->checkCollision();

}


//void gcgSph::handleCollisions(){
//    for(int i = 0; i < nParticles; i++){
//        gcgParticleSPH * p = particles[i];
//        p->checkCollision();
//}
//Zerar as listas de neightbouts


void removeFromVector(gcgParticleSPH * source, gcgParticleSPH * element){
    int vectorIndex = 0;
    if(source->neighbours.size() > 0){
        for(std::vector<gcgParticleSPH*>::iterator it = source->neighbours.begin(); it != source->neighbours.end(); ++it){
            if(*it == element){
//                printf("Removendo %p de %p\n", element, source);
                source->neighbours.erase(source->neighbours.begin() + vectorIndex);
                break;
            }
            vectorIndex++;
        }
    }
}


void gcgSph::iterate() {
//    if(((simulationStep+1) % 300) == 0) saveState();
    #if (INTEGRATOR == VERLET)
        for(int i = 0; i < nParticles; i++) {
            gcgParticleSPH* p = particles[i];
            VECTOR3 vt, at;::itera::itera
            gcgSCALEVECTOR3(vt, p->velocity, TIME_STEP);
            gcgSCALEVECTOR3(at, p->acceleration, TIME_STEP2*0.5);
            checkBoundaryCondition(p, at, vt);
            ///Fazendo no boundary
//            gcgADDVECTOR3(p->position, p->position, vt);
//            gcgADDVECTOR3(p->position, p->position, at);
            gcgCOPYVECTOR3(p->lastAcceleration, p->acceleration);
            gcgSETVECTOR3(p->acceleration, 0.0, 0.0, 0.0);
            gcgCOPYVECTOR3(p->lastVelocity, p->velocity);
        }
    #endif

    /*************************************/
    ///Finding neighbours
    ////////////////////////////////////////
//    for(int i = 0; i < nParticles; i++){
//        gcgParticleSPH * p = particles[i];
//        p->neighbourCheck = true;
//        for(int j = i; j < nParticles; j++){ //OLHAR A VERIFICACAO DE VIZINHOS
//            gcgParticleSPH * p2 = particles[j];
//            if(!p2->neighbourCheck){
////                if(getDistance(p->position, p2->position) <= p->influenceRadius){
//                float dist = getDistance(p->position, p2->position);
//                if(dist <= supportRadius && !FEQUAL(dist, 0.0)){
////                    printf("Adicionando P: %p - %d | P2: %p - %d\n", p, p->neighbours.size(), p2, p2->neighbours.size());
//                    p->neighbours.push_back(p2);
//                    p2->neighbours.push_back(p);
//                }
//            }
//        }
//    }
     ///Olhar novamente o resetIteration
     maxDebugVal = -99999999;
     minDebugVal = 99999999;

    for(int i = 0; i < nParticles; i++){
        gcgParticleSPH * p = particles[i];
        grid->findNeighbours(p);
    }

    ///////////////////////////////////////////


    /////////////SPH/////////////////////////////
    ///Termo de pressao e massa especifica

    for(int i = 0; i < nParticles; i++){
        gcgParticleSPH * p = particles[i];
        if(p->neighbours.size() <= 0) {
            p->massDensity = 0;
        }
        else{
            bool dbg = false;
//            if(simulationTime >= 6.8 && i==399){dbg = true;}
            p->calculateMassDensity(dbg);
        }

        p->calculatePressureTerm();

    //        printf("I = %d | Termo de pressao: %f\n", i, particles[i]->pressureTerm);
    //        printf("Termo de pressao: %f %f %f\n", particles[i]->pressureTerm[0], particles[i]->pressureTerm[1], particles[i]->pressureTerm[2]);
    }

    ///Forcas
    for(int i = 0; i < nParticles; i++){

        ///Limite do suporte compacto
        gcgParticleSPH * p = particles[i];


        if(FEQUAL(p->massDensity*p->massDensity, 0.0)){
//            printf("Size %p: %d\n", p, p->neighbours.size());
//            printf("Size %p: %d | md: %f\n", p, p->neighbours.size(), p->massDensity);

//            if(p->neighbours.size() > 0){
//                printf("Removendo!\n");
                for(std::vector<gcgParticleSPH*>::iterator it = p->neighbours.begin(); it != p->neighbours.end(); ++it){
                    gcgParticleSPH * n1 = *it;
//                    printf("Removendo %p de %p\n", n1, p);
                    removeFromVector(n1, p);


                }
                p->neighbours.clear();
//            }
            
            
        }

        if(p->neighbours.size() > 0){

//            if(i == 108){
//                printf("dbg\n");
//            }
            p->calculatePressure();
            p->calculateViscosity();
        }



//        p->calculateGravity();

//        p->calculateTensorForce();
//        VECTOR3 eita;
//        p->calculateExternalForces(eita);

        /***************    TENSOR FORCE       ***************/
//        if(p->position[0] < 0 || p->position[0] >= width || p->position[1] < 0 || p->position[1] >= height || p->position[2] < 0 || p->position[2] >= depth){
//            gcgZEROVECTOR3(p->tensorForce);
//        }
//        else{
//
//            int posX = (unsigned int) p->position[0];// * (width  - 1)) / (unsigned int) width;
//            int posY = (unsigned int) p->position[1];// * (height - 1)) / (unsigned int) height;
//            int posZ = (unsigned int) p->position[2];// * (depth  - 1)) / (unsigned int) depth;
//            int v = ADDRESSVOLUME(posX, posY, posZ, width, height);
//
//            p->calculateExternalForces(teste->velField[v]);
//
//            VECTOR3 tt;
//            gcgSETVECTOR3(tt, teste->glyphs_original[v]->eigenVectors[0], teste->glyphs_original[v]->eigenVectors[1], teste->glyphs_original[v]->eigenVectors[2]);
//
//
////            printf("Pos: %d - Force: %f %f %f\n", v, tt[0],tt[1], tt[2]);
//
//            p->calculateTensorForce(teste->glyphs_original, teste->bigger, posX, posY, posZ, teste->velField[v], 1., teste->glyphs_original[v]->cl, v);
////            gcgADDVECTOR3(p->position, p->position, p->tensorForce);
////            return;
////            gcgCOPYVECTOR3(debugTensor, tt);///desenhar tensor
////            printf("Tensor: %f -- %f %f %f\n", gcgLENGTHVECTOR3(p->tensorForce), p->tensorForce[0], p->tensorForce[1], p->tensorForce[2]);
////            gcgADDVECTOR3(p->externalForces, p->externalForces, p->tensorForce);
//
//            gcgSETVECTOR3(p->tensorForce, 0., 0., 0.);
//
//        }
        /***************************************************/
    }


//    for(int i = 0; i < nParticles; i++){
//        gcgParticleSPH * p = particles[i];
//        if(p->neighbours.size() <= 0) {
//            p->massDensity = 0;
//        }
//        else{
//            for(std::vector<gcgParticleSPH*>::iterator it = p->neighbours.begin(); it != p->neighbours.end(); ++it) {
//    //        for(int j = 0; j < nParticles; j++){
//                p->calculateMassDensity(*it);
//    //                printf("I = %d | J = %d | Massa especifica: %f\n", i, j, particles[i]->massDensity);
//            }
//            p->calculatePressureTerm();
//    //        printf("I = %d | Termo de pressao: %f\n", i, particles[i]->pressureTerm);
//    //        printf("Termo de pressao: %f %f %f\n", particles[i]->pressureTerm[0], particles[i]->pressureTerm[1], particles[i]->pressureTerm[2]);
//        }
//    }
//
//    ///Forcas
//    for(int i = 0; i < nParticles; i++){
//        gcgParticleSPH * p = particles[i];
//        for(std::vector<gcgParticleSPH*>::iterator it = p->neighbours.begin(); it != p->neighbours.end(); ++it) {
//                p->calculatePressure(*it);
//                p->calculateViscosity(*it);
//        }
//        p->calculateGravity();
//    }

    /////////////////////////SPH////////////////////////
    ///Nova aceleracao
    #if (INTEGRATOR == VERLET)
        for(int i = 0; i < nParticles; i++) {
            gcgParticleSPH* p = particles[i];
            /////////Navier Stokes
            gcgADDVECTOR3(p->acceleration, p->pressure, p->viscosity);
            gcgADDVECTOR3(p->acceleration, p->acceleration, p->externalForces);
            gcgSCALEVECTOR3(p->acceleration, p->acceleration, 1./p->massDensity);
            printf("Aceleracao: %f %f %f\n", p->acceleration[0], p->acceleration[1], p->acceleration[2]);
    //        gcgCOPYVECTOR3(p->acceleration, p->gravity);
            ////////////////////////////
            VECTOR3 at;
            gcgADDVECTOR3(at, p->lastAcceleration, p->acceleration);
            gcgSCALEVECTOR3(at, at, 0.5*TIME_STEP);
            gcgADDVECTOR3(p->velocity, p->lastVelocity, at);

            ///Recalcula o raio
//            p->radius = radiusScale * (pow((3.*p->mass/(float)(4.*M_PI*p->massDensity)),(1./3.)));
        }

    #elif (INTEGRATOR == LEAP_FROG)

        for(int i = 0; i < nParticles; i++) {
            gcgParticleSPH* p = particles[i];
//            if(simulationStept >= 6.8 && i == 399){printf("Partícula %d - Massa específica: %f - Termo de pressão: %f\nPos: %f %f %f\nPressao: %f %f %f\nViscosidade: %f %f %f\nGrad: %f %f %f\n", i, p->massDensity, p->pressureTerm, p->position[0], p->position[1], p->position[2], p->pressure[0], p->pressure[1],p->pressure[2], p->viscosity[0], p->viscosity[1],p->viscosity[2], p->gradientVal[0], p->gradientVal[1], p->gradientVal[2]);}
            if(p->neighbours.size() > 0){
            gcgADDVECTOR3(p->acceleration, p->pressure, p->viscosity);
            }
            /************Desabilitando externas************/
//            gcgADDVECTOR3(p->acceleration, p->acceleration, p->gravity);
//            gcgADDVECTOR3(p->acceleration, p->acceleration, p->tensorForce);

//            gcgADDVECTOR3(p->acceleration, p->acceleration, p->externalForces);


            /********DEBUG***********/
//            float curVal = gcgLENGTHVECTOR3(p->acceleration);
//            if(curVal >= maxAcc){
//                gcgCOPYVECTOR3(maxVec, p->acceleration);
//                maxAcc = curVal;
//            }
//            if(curVal <= minAcc){
//                gcgCOPYVECTOR3(minVec, p->acceleration);
//                minAcc = curVal;
//            }
            /********DEBUG***********/

            if(!FEQUAL(p->massDensity, 0)){
            gcgSCALEVECTOR3(p->acceleration, p->acceleration, 1./p->massDensity);
            }

            /*******************************/
            if(useExternals){
                VECTOR3 gVal = {0,0,0};
                gcgSCALEVECTOR3(gVal, p->gradientVal, 1./p->mass);
//                if(!FEQUAL(p->massDensity, 0)){
//                        gcgSCALEVECTOR3(gVal, p->gradientVal, 1./p->massDensity);
//                }
                gcgADDVECTOR3(p->acceleration, p->acceleration, gVal);

            }
            VECTOR3 slapV = {0,0,0};
            if(((int)(simulationStep+1) % slapMax) == 0){
                if(slapX) slapV[0] = (float) (gcgRandom.random());
                if(slapY) slapV[1] = (float) (gcgRandom.random());
                if(tipoCampo != FIELD2D){ if(slapZ) slapV[2] = (float) (gcgRandom.random());}
                if(slapEigen) gcgSETVECTOR3(slapV, p->tensor->eigenVectors[0], p->tensor->eigenVectors[1], p->tensor->eigenVectors[2]);
                gcgSCALEVECTOR3(slapV, slapV, 20.0);
//                printf("%f %f %f %f %f %f\n", (gcgRandom.random()), (gcgRandom.random()), (gcgRandom.random()), (gcgRandom.random()), (gcgRandom.random()), (gcgRandom.random()));
                gcgADDVECTOR3(p->acceleration, p->acceleration, slapV);

            }
//            gcgCOPYVECTOR3(p->acceleration, p->tensorForce); ///tirar!!!
//            gcgCOPYVECTOR3(p->acceleration, p->gradientVal); ///tirar!!!
//            gcgCOPYVECTOR3(p->acceleration, p->buoyancy); ///tirar!!!
            /**********************************/

//            if(p->acceleration[0] > 50) {
//                printf("Aceleração: %f\n", p->acceleration[0]);
//            }

//            if(gcgLENGTHVECTOR3(p->acceleration) > 100){
//                printf("High acc\n");
//            }
            gcgSCALEVECTOR3(p->nextVelocity, p->acceleration, TIME_STEP);
            gcgADDVECTOR3(p->nextVelocity, p->nextVelocity, p->lastVelocity);


//            if(p->acceleration[0] != 0){
//                printf("Maior\n");
//            }


            checkBoundaryCondition(p, NULL, NULL);

            gcgADDVECTOR3(p->velocity, p->lastVelocity, p->nextVelocity);

            gcgSCALEVECTOR3(p->velocity, p->velocity, 0.5);

            gcgCOPYVECTOR3(p->lastVelocity, p->nextVelocity);



//            VECTOR3 pos;
//            gcgCOPYVECTOR3(pos, p->position);
//
//            gcgSCALEVECTOR3(p->position, p->nextVelocity, TIME_STEP);
//            gcgADDVECTOR3(p->position, p->position, pos);

            gcgCOPYVECTOR3(p->lastAcceleration, p->acceleration);
//            printf("Aceleracao: %f %f %f\n", p->acceleration[0], p->acceleration[1], p->acceleration[2]);
            gcgSETVECTOR3(p->acceleration, 0.0, 0.0, 0.0);


//            p->radius = radiusScale * (pow((0.75*p->mass/(float)(3*M_PI*p->massDensity)),(1./3.))); /// Fixo

            grid->moveParticle(p);

            p->resetIteration();
            if(tipoCampo == FIELD2D) interpolateTensor(p);
            else {
                bool dbg = false;
//                if(simulationTime >= 6.8) dbg = true;
                if(dbg) printf("i: %d --", i);
                interpolateTensor3D(p, dbg, i);
            }
            
            /**********SQUAD**********/
                
                
            p->defineSQuad(p->tensor);

            /********************/
//            interpolateTensor3D(p);
//            interpolateTensor(p);

            numberType tmpVal;
            switch(informationType){
                case 0:{
                    tmpVal = p->lastMassDensity;
                    break;
                }
                case 1:{
                    tmpVal = gcgLENGTHVECTOR3(p->lastAcceleration);
                    break;
                }
                case 2:{
                    tmpVal = gcgLENGTHVECTOR3(p->lastViscosity);
                    break;
                }
                case 3:{
                    tmpVal = gcgLENGTHVECTOR3(p->lastPressure);
                    break;
                }

            }

            if(tmpVal >= maxDebugVal){
                maxDebugVal = tmpVal;
            }
            if(tmpVal <= minDebugVal){
                minDebugVal = tmpVal;
            }



        }



//
//
//        for(int i = 0; i < nParticles; i++) {
//            gcgParticleSPH* p = particles[i];
//
//            if(p->neighbours.size() > 0){
//            gcgADDVECTOR3(p->acceleration, p->pressure, p->viscosity);
//            }
//
////            gcgADDVECTOR3(p->acceleration, p->acceleration, p->gravity);
//            gcgADDVECTOR3(p->acceleration, p->acceleration, p->externalForces);
//
//            if(p->neighbours.size() > 0){
//            gcgSCALEVECTOR3(p->acceleration, p->acceleration, 1./p->massDensity);
//            }
//
//            gcgSCALEVECTOR3(p->nextVelocity, p->acceleration, TIME_STEP);
//            gcgADDVECTOR3(p->nextVelocity, p->nextVelocity, p->lastVelocity);
//
//            checkBoundaryCondition(p, NULL, NULL);
//
//            gcgCOPYVECTOR3(p->lastVelocity, p->nextVelocity);
//
////            VECTOR3 pos;
////            gcgCOPYVECTOR3(pos, p->position);
////
////            gcgSCALEVECTOR3(p->position, p->nextVelocity, TIME_STEP);
////            gcgADDVECTOR3(p->position, p->position, pos);
//
//            gcgCOPYVECTOR3(p->lastAcceleration, p->acceleration);
////            printf("Aceleracao: %f %f %f\n", p->acceleration[0], p->acceleration[1], p->acceleration[2]);
////            printf("Velocidade: %f %f %f\n", p->acceleration[0], p->acceleration[1], p->acceleration[2]);
//            gcgSETVECTOR3(p->acceleration, 0.0, 0.0, 0.0);
//
////            p->radius = radiusScale * (pow((0.75*p->mass/(float)(3*M_PI*p->massDensity)),(1./3.))); /// Fixo
//            p->resetIteration();
//        }
    #endif

//     MATRIX3 eM, tmpM;
//     gcgParticleSPH * p1 = particles[0];
//     gcgParticleSPH * p2 = particles[1];
////    gcgSETMATRIX3(eM,
////               p1->tensor->eigenVectors[0] * p1->tensor->eigenValues[0], p1->tensor->eigenVectors[1] * p1->tensor->eigenValues[0], p1->tensor->eigenVectors[2] * p1->tensor->eigenValues[0],
////               p1->tensor->eigenVectors[3] * p1->tensor->eigenValues[1], p1->tensor->eigenVectors[4] * p1->tensor->eigenValues[1], p1->tensor->eigenVectors[5] * p1->tensor->eigenValues[1],
////               p1->tensor->eigenVectors[6] * p1->tensor->eigenValues[2], p1->tensor->eigenVectors[7] * p1->tensor->eigenValues[2], p1->tensor->eigenVectors[8] * p1->tensor->eigenValues[2]);
//    gcgSETMATRIX3(eM,
//               p1->tensor->eigenVectors[3] * p1->tensor->eigenValues[1], p1->tensor->eigenVectors[4] * p1->tensor->eigenValues[1], p1->tensor->eigenVectors[5] * p1->tensor->eigenValues[1],
//               p1->tensor->eigenVectors[0] * p1->tensor->eigenValues[0], p1->tensor->eigenVectors[1] * p1->tensor->eigenValues[0], p1->tensor->eigenVectors[2] * p1->tensor->eigenValues[0],
//               p1->tensor->eigenVectors[6] * p1->tensor->eigenValues[2], p1->tensor->eigenVectors[7] * p1->tensor->eigenValues[2], p1->tensor->eigenVectors[8] * p1->tensor->eigenValues[2]);
////    //                  gcgTRANSPOSEMATRIX3(tmpM2, tmpM);
////    //                  gcgCOPYMATRIX3(tmpM, tmpM2);
////    gcgAPPLYMATRIX3VECTOR3(pf, tmpM, p);
//     VECTOR3 v, pf, ptmp;
//     float lv = gcgLENGTHVECTOR3(p1->position);
//     gcgSUBVECTOR3(v, p2->position, p1->position);
//     gcgCOPYVECTOR3(p1->debugDif, v);
////     gcgSCALEVECTOR3(ptmp, v, supportRadius/lv);
////     gcgAPPLYMATRIX3VECTOR3(pf, p1->tensor->tensor, ptmp);
//     gcgAPPLYMATRIX3VECTOR3(pf, eM, v);
//
////     printf("%f < %f\n", lv, gcgLENGTHVECTOR3(pf));
//
//     gcgAPPLYMATRIX3VECTOR3(pf, p1->tensor->tensor, p2->position);
//
//     lv = ((((pf[0]-p1->position[0])*(pf[0]-p1->position[0]))/(p1->tensor->eigenValues[1]*p1->tensor->eigenValues[1])) +((pf[1]-p1->position[1])*(pf[1]-p1->position[1]))/(p1->tensor->eigenValues[0]*p1->tensor->eigenValues[0]));
////     lv = ((((pf[0])*(pf[0]))/(p1->tensor->eigenValues[0]*p1->tensor->eigenValues[0])) +((pf[1])*(pf[1]))/(p1->tensor->eigenValues[1]*p1->tensor->eigenValues[1]));
//     printf("%f < 1.\n", lv);

//    gcgSETMATRIX3(eM,
//     p1->tensor->eigenVectors[0], p1->tensor->eigenVectors[1], p1->tensor->eigenVectors[2],
//     p1->tensor->eigenVectors[3], p1->tensor->eigenVectors[4], p1->tensor->eigenVectors[5],
//     p1->tensor->eigenVectors[6], p1->tensor->eigenVectors[7], p1->tensor->eigenVectors[8]);

//    gcgSETMATRIX3(eM,
//     p1->tensor->eigenVectors[3], p1->tensor->eigenVectors[4], p1->tensor->eigenVectors[5],
//     p1->tensor->eigenVectors[0], p1->tensor->eigenVectors[1], p1->tensor->eigenVectors[2],
//     p1->tensor->eigenVectors[6], p1->tensor->eigenVectors[7], p1->tensor->eigenVectors[8]);

//    gcgAPPLYMATRIX3VECTOR3(pf, eM, v);
//    gcgCOPYVECTOR3(v, pf);
//    gcgAPPLYMATRIX3VECTOR3(v, p1->tensor->tensor, pf);


//    gcgAPPLYMATRIX3VECTOR3(pf, p1->tensor->tensor, v);
//    gcgCOPYVECTOR3(v, pf);
//    gcgAPPLYMATRIX3VECTOR3(v, eM, pf);



    //         if((gcgLENGTHVECTOR3(p)/supportRadius) < 0.) printf("Dentro!\n");
//    printf("%f = %f\n", gcgLENGTHVECTOR3(v), supportRadius);
//if((gcgLENGTHVECTOR3(p)/supportRadius) < 1.) printf("Dentro!\n");






     simulationTime += TIME_STEP;
     simulationStep++;

//     printf("T = %f\n", simulationTime);
}

void gcgSph::drawBox() {

    glPushMatrix();
    glTranslatef(translateW, translateH, translateD);

    glColor3f(0., 0., 0.8);
    glBegin(GL_LINE_LOOP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(boxWidth , 0.0, 0.0);
    glVertex3f(boxWidth , boxHeight , 0.0);
    glVertex3f(0, boxHeight , 0.0);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3f(0.0, 0.0, boxDepth );
    glVertex3f(boxWidth , 0.0, boxDepth );
    glVertex3f(boxWidth , boxHeight , boxDepth );
    glVertex3f(0, boxHeight , boxDepth );
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, boxDepth );

    glVertex3f(boxWidth , 0.0, 0.0);
    glVertex3f(boxWidth , 0.0, boxDepth );


    glVertex3f(0.0, boxHeight , 0.0);
    glVertex3f(0.0, boxHeight , boxDepth );

    glVertex3f(boxWidth , boxHeight , 0.0);
    glVertex3f(boxWidth , boxHeight , boxDepth );

    glEnd();

    glPopMatrix();
}

void DrawEllipse(numberType a1, numberType a2, VECTOR3 e1, MATRIX3 tensor, MATRIX3 evectors, VECTOR3 evals, gcgSph * sph, gcgParticleSPH * part)
{

	numberType x,y,z;
	int t;
        glPointSize(2.5);


        for(t = 0; t < 360; t +=10)
	{
          glPointSize(1.5);
          glBegin(GL_POINTS);
            glColor4f(0.5,0.5,0.5,1.);
          numberType dg = ((numberType)t) * M_PI / 180.0;
          x = supportRadius*sin(dg);
	  y = supportRadius*cos(dg);
	  z = 0;

          VECTOR3 p = {x, y, z};
          VECTOR3 pf;
          gcgAPPLYMATRIX3VECTOR3(pf, tensor, p);
          gcgCOPYVECTOR3(pf, p);

//	 glVertex3f(x,y,z);
//         if(t == 0){
//                glEnd();
//                glPointSize(4.);
//                glBegin(GL_POINTS);
//                glVertex3f(pf[0],pf[1],pf[2]);
//                printf("SPH - Sup: %f - Tensor: %p - Ang: %f Sin %f Cos %f| %f %f %f\n", supportRadius, tensor, dg, sin(dg), cos(dg), pf[0], pf[1], pf[2]);
//                glEnd();
//                glPointSize(2.);
//                glBegin(GL_POINTS);
//                continue;
//            }

	 glVertex3f(pf[0],pf[1],pf[2]);

         MATRIX3 tmpM, tmpM2;
                  gcgSETMATRIX3(tmpM,
                          part->tensor->eigenVectors[3] * part->tensor->eigenValues[1], part->tensor->eigenVectors[4] * part->tensor->eigenValues[1], part->tensor->eigenVectors[5] * part->tensor->eigenValues[1],
                          part->tensor->eigenVectors[0] * part->tensor->eigenValues[0], part->tensor->eigenVectors[1] * part->tensor->eigenValues[0], part->tensor->eigenVectors[2] * part->tensor->eigenValues[0],
                          part->tensor->eigenVectors[6] * part->tensor->eigenValues[2], part->tensor->eigenVectors[7] * part->tensor->eigenValues[2], part->tensor->eigenVectors[8] * part->tensor->eigenValues[2]);
////                  gcgTRANSPOSEMATRIX3(tmpM2, tmpM);
////                  gcgCOPYMATRIX3(tmpM, tmpM2);
//         gcgAPPLYMATRIX3VECTOR3(pf, tmpM, p);
         glEnd();
//         gcgSCALEMATRIX3(tmpM, tmpM, supportRadius);
         glPointSize(3.5);

         if(part == sph->particles[0]){
          glBegin(GL_POINTS);
         glColor4f(0.,1.,0.,1.);
//         gcgAPPLYMATRIX3VECTOR3(pf, sph->particles[0]->tensor->eigenVectors, sph->particles[0]->tensor->eigenValues);

//         gcgAPPLYMATRIX3VECTOR3(pf, tmpM, p);
//         gcgSCALEVECTOR3(pf, pf, supportRadius);
//         gcgSETMATRIX3(tmpM2, part->tensor->tensor[3], part->tensor->tensor[4], part->tensor->tensor[5], part->tensor->tensor[0], part->tensor->tensor[1], part->tensor->tensor[2], part->tensor->tensor[6], part->tensor->tensor[7], part->tensor->tensor[8]);
//         gcgAPPLYMATRIX3VECTOR3(p, part->tensor->tensor, pf);
//         gcgCOPYVECTOR3(pf, p);

//         gcgCOPYVECTOR3(p, sph->particles[1]->position);
         gcgSUBVECTOR3(p, sph->particles[1]->position, sph->particles[0]->position);
//         gcgAPPLYMATRIX3VECTOR3(pf, part->tensor->tensor, p);
         gcgAPPLYMATRIX3VECTOR3(pf, tmpM, p);


         glVertex3f(pf[0],pf[1],pf[2]);


//        gcgSETMATRIX3(tmpM,
//                   sph->particles[0]->tensor->eigenVectors[0], sph->particles[0]->tensor->eigenVectors[1], sph->particles[0]->tensor->eigenVectors[2],
//                   sph->particles[0]->tensor->eigentmpMVectors[3], sph->particles[0]->tensor->eigenVectors[4], sph->particles[0]->tensor->eigenVectors[5],
//                   sph->particles[0]->tensor->eigenVectors[6], sph->particles[0]->tensor->eigenVectors[7], sph->particles[0]->tensor->eigenVectors[8]);
//         gcgSUBVECTOR3(p, sph->particles[1]->position, sph->particles[0]->position);
//         glEnd();
//         gcgDrawVectorPyramid(0., 0., 0., p, 1.0);
//         gcgAPPLYMATRIX3VECTOR3(pf, tmpM, p);
//         glColor4f(0.,0.,1.,1.);
//         gcgAPPLYMATRIX3VECTOR3(p, sph->particles[0]->tensor->tensor, pf);
//         glBegin(GL_POINTS);
//         glVertex3f(p[0],p[1],p[2]);
////         if((gcgLENGTHVECTOR3(p)/supportRadius) < 0.) printf("Dentro!\n");
//         float aaa2 = gcgLENGTHVECTOR3(p)/supportRadius;
//         printf("%f = %f\n", gcgLENGTHVECTOR3(p), supportRadius);
//         if((gcgLENGTHVECTOR3(p)/supportRadius) < 1.) printf("Dentro!\n");
         glEnd();
         }
       }
       glEnd();
       glPointSize(1.0);
}

void gcgSph::drawParticles() {

    glPushMatrix();
//    glTranslatef(3, 0., 0.);
    glTranslatef(translateW, translateH, translateD);
    glEnable(GL_LIGHTING);
//    glDisable(GL_LIGHTING);
    numberType subDensity = maxDebugVal - minDebugVal;
    numberType min = 999999999, max = -99999999999;
    int midP = floor((float)nParticles / 2.0) + 1;
    for(int i = 0; i<nParticles; i++){


        gcgParticleSPH* p = particles[i];

        glPushMatrix();
        glTranslatef(p->position[0],p->position[1], p->position[2]);
//        printf("Posicao: %f %f %f\n", p->position[0],p->position[1], p->position[2]);
        VECTOR3 color;
        numberType heat;
        numberType tmpVal;

        switch(informationType){
            case 0:{
                tmpVal = p->lastMassDensity;
                break;
            }
            case 1:{
                tmpVal = gcgLENGTHVECTOR3(p->lastAcceleration);
                break;
            }
            case 2:{
                tmpVal = gcgLENGTHVECTOR3(p->lastViscosity);
                break;
            }
            case 3:{
                tmpVal = gcgLENGTHVECTOR3(p->lastPressure);
                break;
            }

        }

        heat = numberType((tmpVal - minDebugVal)/subDensity);
        if(heat > max){ max = heat; }
        if(heat < min){ min = heat; }
        gcgHeatColor(heat, color);
        glColor3f(color[0], color[1], color[2]);
//        glColor3f(1.0, 0., 0.);
//        glColor3f(1., 0., 0.);
//        glutSolidSphere(p->radius, 20, 20);
        
        if(cutVis){
            if(!((p->position[0] > cutBounds1[0]) && (p->position[0] < cutBounds2[0]) && (p->position[1] > cutBounds1[1]) && (p->position[1] < cutBounds2[1]) && (p->position[2] > cutBounds1[2]) && (p->position[2] < cutBounds2[2]))) {glPopMatrix(); continue;}
        }
        
        if(filterParticles && (heat < filterT)){glPopMatrix(); continue;}
        
        glutSolidSphere(particlesRadius, 20, 20);


//        glColor4f(0,0,1., 1.);
        VECTOR3 dir;
        gcgSETVECTOR3(dir, p->tensor->eigenVectors[0], p->tensor->eigenVectors[1], p->tensor->eigenVectors[2]);
        
//                gcgDrawVectorPyramid(0., 0., 0., dir, 1.0f);
//                gcgDrawVectorPyramid(0., 0., 0., p->gradientVal, 1.0f);

        VECTOR3 e1;
        gcgSETVECTOR3(e1, p->tensor->eigenVectors[0], p->tensor->eigenVectors[1], p->tensor->eigenVectors[2]);
//        DrawEllipse(p->tensor->eigenVectors[0], p->tensor->eigenVectors[1], e1, p->tensor->tensor, p->tensor->eigenVectors, p->tensor->eigenValues, this, p);
//        glPopMatrix();


        if(i == midP){
//            printf("Pos: %f %f %f\n", p->position[0], p->position[1], p->position[2]);
        }


        glPopMatrix();

    }
    glPopMatrix();

//    printf("Min %f - Max %f\n", min, max);

}


void gcgSph::drawParticlesSuperQuadricTassio() {
   
    glDisable(GL_LIGHTING);
//    glEnable(GL_LIGHTING);
    glPushMatrix();
    glTranslatef(translateW, translateH, translateD);
    
    vector<gcgParticleSPH*>::iterator it;
    float cs;
    int pos;
    VECTOR3 posAVector;
    MATRIX4 matrix;
    numberType heat;
    numberType tmpVal;
    MATRIX3 color;
    numberType subDensity = maxDebugVal - minDebugVal;

    
    glMatrixMode(GL_MODELVIEW);
    
    for (int i = 0; i < nParticles; i++) {
        gcgParticleSPH* p = particles[i];

        
        switch(informationType){
            case 0:{
                tmpVal = p->lastMassDensity;
                break;
            }
            case 1:{
                tmpVal = gcgLENGTHVECTOR3(p->lastAcceleration);
                break;
            }
            case 2:{
                tmpVal = gcgLENGTHVECTOR3(p->lastViscosity);
                break;
            }
            case 3:{
                tmpVal = gcgLENGTHVECTOR3(p->lastPressure);
                break;
            }
        }
        
        heat = numberType((tmpVal - minDebugVal)/subDensity);
         if(cutVis){
            if(!((p->position[0] > cutBounds1[0]) && (p->position[0] < cutBounds2[0]) && (p->position[1] > cutBounds1[1]) && (p->position[1] < cutBounds2[1]) && (p->position[2] > cutBounds1[2]) && (p->position[2] < cutBounds2[2]))) {continue;}
        }
        if((filterParticles && (heat < filterT))|| p->tensor->cs >0.7){continue;}
        gcgHeatColor(heat, color);
        
//        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        
        glPushMatrix();
        glTranslatef(p->position[0],p->position[1], p->position[2]);
//        pos = ADDRESSVOLUME((int) p->position[0], (int) p->position[1], (int) p->position[2], width, height);
        cs = 1.0;

//        glTranslatef(p->position[0], p->position[1], p->position[2]);
//            //VERIFICAR!!!
//            gcgSCALEVECTOR3(posAVector, velField[pos], p->direc);
//
//            //  Compute rotation matrix
        gcgSETVECTOR3(posAVector, p->tensor->eigenVectors[0], p->tensor->eigenVectors[1], p->tensor->eigenVectors[2]);
            
        
        MATRIX4 orientation = {p->tensor->eigenVectors[0], p->tensor->eigenVectors[1],p->tensor->eigenVectors[2], 0,
                       p->tensor->eigenVectors[3],p->tensor->eigenVectors[4],p->tensor->eigenVectors[5], 0,
                       p->tensor->eigenVectors[6],p->tensor->eigenVectors[7],p->tensor->eigenVectors[8], 0,
                       0, 0, 0, 1};
        glMultMatrixf(orientation);
        
        
//        gcgComputeAlignMatrix(matrix, posAVector);
//            glMultMatrixf(matrix);
//
////            if (original) gcgHeatColor(p->tensor->equ, p->color);
//            gcgSETVECTOR3(color)
//
//            glColor4f(p->color[0], p->color[1], p->color[2], 1.0);
      
//        gcgHeatColor(p->tensor->cl, color);
            glColor4f(color[0], color[1], color[2], 1.0);
       
//        glColor4f(0.6, 0.6, 0.6, 1.0);
        

//        if(p->position[1] <= 25)  p->sq->draw(true, false, true);
//        else  p->sq->draw(true, false, false);
        p->sq->draw(true, false, false);


        glPopMatrix();
    }
    
    glPopMatrix();

}

gcgSph::~gcgSph() {
}

void gcgSph::saveState(){
    
    char fileName[9];
    FILE * f;
    sprintf(fileName, "sphSave%d", state);
    state++;
    
    f = fopen(fileName, "w");
    
    
    fwrite(&width, 1, sizeof(unsigned int), f);
    fwrite(&height, 1, sizeof(unsigned int), f);
    fwrite(&depth, 1, sizeof(unsigned int), f);
    fwrite(&widthFluid, 1, sizeof(unsigned int), f);
    fwrite(&heightFluid, 1, sizeof(unsigned int), f);
    fwrite(&depthFluid, 1, sizeof(unsigned int), f);
    fwrite(&nParticles, 1, sizeof(unsigned int), f);
    fwrite(&supportRadius, 1, sizeof(numberType), f);
//    fwrite(&restDensity, 1, sizeof(numberType), f);
    fwrite(&mi, 1, sizeof(numberType), f);
    fwrite(&kConst, 1, sizeof(numberType), f);
    
    for(int i = 0; i < nParticles; i++){
        fwrite(particles[i]->position, 1, sizeof(VECTOR3), f);
        fwrite(particles[i]->lastPosition, 1, sizeof(VECTOR3), f);
        fwrite(&particles[i]->mass, 1, sizeof(numberType), f);
        fwrite(particles[i]->lastAcceleration, 1, sizeof(VECTOR3), f);
        fwrite(particles[i]->velocity, 1, sizeof(VECTOR3), f);
        fwrite(particles[i]->lastVelocity, 1, sizeof(VECTOR3), f);
        fwrite(particles[i]->gradientVal, 1, sizeof(VECTOR3), f);
        fwrite(particles[i]->tensor->tensor, 1, sizeof(MATRIX3), f);
    }
    
    fclose(f);
}


void gcgSph::exportState(){
    
    char fileName[12];
    FILE * f;
    sprintf(fileName, "sphExport%d", xport);
    xport++;
    
    f = fopen(fileName, "w");
    
    
    fwrite(&width, 1, sizeof(unsigned int), f);
    fwrite(&height, 1, sizeof(unsigned int), f);
    fwrite(&depth, 1, sizeof(unsigned int), f);
    fwrite(&nParticles, 1, sizeof(unsigned int), f);
    fwrite(&particlesMass, 1, sizeof(numberType), f);
    for(int i = 0; i < nParticles; i++){
        fwrite(particles[i]->position, 1, sizeof(VECTOR3), f);
        fwrite(&particles[i]->massDensity, 1, sizeof(numberType), f);
    }
    
    fclose(f);
}


void importState(char * fName){
    
    
    FILE * f = fopen(fName, "r");
    
    int width, height, depth;
    int nParticles;
    float mass, mDensity;
    VECTOR3 pos;
    
    fread(&width, 1, sizeof(unsigned int), f);
    fread(&height, 1, sizeof(unsigned int), f);
    fread(&depth, 1, sizeof(unsigned int), f);
    fread(&nParticles, 1, sizeof(unsigned int), f);
    fread(&mass, 1, sizeof(numberType), f);
    for(int i = 0; i < nParticles; i++){
        fread(pos, 1, sizeof(VECTOR3), f);
        fread(&mDensity, 1, sizeof(numberType), f);
    }
    
    fclose(f);
}




gcgSph * loadState(int val){
    
    char fileName[9];
    FILE * f;
    sprintf(fileName, "sphSave%d", val);
    
    f = fopen(fileName, "r");
    
    unsigned int w, h, d, wf, hf, df, nP;
    
    fread(&w, 1, sizeof(unsigned int), f);
    fread(&h, 1, sizeof(unsigned int), f);
    fread(&d, 1, sizeof(unsigned int), f);
    fread(&wf, 1, sizeof(unsigned int), f);
    fread(&hf, 1, sizeof(unsigned int), f);
    fread(&df, 1, sizeof(unsigned int), f);
    fread(&nP, 1, sizeof(unsigned int), f);
    
    gcgSph * sph = new gcgSph(nP, wf, hf, df, w, h, d);
            
    numberType tVal;
    fread(&tVal, 1, sizeof(numberType), f); ///Support Radius
    setSupRadius(tVal);
    
//    fread(&restDensity, 1, sizeof(numberType), f);
    fread(&mi, 1, sizeof(numberType), f);
    fread(&kConst, 1, sizeof(numberType), f);
    
    for(int i = 0; i < nP; i++){
        sph->particles[i]->resetIteration();
        fread(sph->particles[i]->position, 1, sizeof(VECTOR3), f);
        fread(sph->particles[i]->lastPosition, 1, sizeof(VECTOR3), f);
        fread(&sph->particles[i]->mass, 1, sizeof(numberType), f);
        fread(sph->particles[i]->lastAcceleration, 1, sizeof(VECTOR3), f);
        fread(sph->particles[i]->velocity, 1, sizeof(VECTOR3), f);
        fread(sph->particles[i]->lastVelocity, 1, sizeof(VECTOR3), f);
        fread(sph->particles[i]->gradientVal, 1, sizeof(VECTOR3), f);
        
        MATRIX3 tempT;
        fread(tempT, 1, sizeof(MATRIX3), f);
        
        sph->particles[i]->tensor->setTensor(tempT);
        sph->grid->moveParticle(sph->particles[i]);
    }
    
    fclose(f);
    
    return sph;
}