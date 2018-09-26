/*
 * File:   octree.cpp
 * Author: joseluiz
 *
 * Created on April 25, 2012, 1:43 PM
 */

#include <ios>

#include "octree.h"


bool drawOctreeBoxes = true;
bool print = false;

//int dbgcounter = 0;

//#include "optimizedGrid.h"

//octree::octree() {
//}
//gcgParticleSPH::gcgParticleSPH(){
//    gcgSETVECTOR3(pos, 0., 0., 0.);
//    radius = 1.;
//}


void Octree::buildTree(){
    if(depth < MAX_OCTREE_DEPTH){
        hasChildren = true;
        for(int x = 0; x < 2; x++){
            float minX;
            float maxX;
            if (x == 0) {
                minX = corner1[0];
                maxX = center[0];
            }
            else {
                minX = center[0];
                maxX = corner2[0];
            }
            for(int y = 0; y < 2; y++){
                float minY;
                float maxY;
                if (y == 0) {
                    minY = corner1[1];
                    maxY = center[1];
                }
                else {
                    minY = center[1];
                    maxY = corner2[1];
                }
                for(int z = 0; z < 2; z++){
                    float minZ;
                    float maxZ;
                    if (z == 0) {
                        minZ = corner1[2];
                        maxZ = center[2];
                    }
                    else {
                        minZ = center[2];
                        maxZ = corner2[2];
                    }
                    VECTOR3 c1, c2;
                    gcgSETVECTOR3(c1, minX, minY, minZ);
                    gcgSETVECTOR3(c2, maxX, maxY, maxZ);
                    children[x][y][z] = new Octree(c1, c2, depth+1, this);
                    children[x][y][z]->buildTree();
                }
            }
        }
    }
    else{
        hasChildren = false;
        return;
    }
}

void Octree::drawOctree(){
    
    if(drawOctreeBoxes){
        glLineWidth(2.);
        glColor3f(0.4, 0.4, 0.4);
        glBegin(GL_LINE_LOOP);
        glColor3f(1., 0., 0.);
        glVertex3f(corner1[0], corner1[1], corner1[2]);
        glVertex3f(corner2[0], corner1[1], corner1[2]);
        glColor3f(0.4, 0.4, 0.4);
        glVertex3f(corner2[0], corner1[1], corner2[2]);
        glVertex3f(corner1[0], corner1[1], corner2[2]);
        glEnd();

        glBegin(GL_LINE_STRIP);
        glVertex3f(corner1[0], corner2[1], corner1[2]);
        glVertex3f(corner2[0], corner2[1], corner1[2]);
        glVertex3f(corner2[0], corner2[1], corner2[2]);
        glColor3f(0., 0., 1.);
        glVertex3f(corner1[0], corner2[1], corner2[2]);
        glVertex3f(corner1[0], corner2[1], corner1[2]);
        glEnd();

        glBegin(GL_LINES);
        glColor3f(0., 1., 0.);
        glVertex3f(corner1[0], corner1[1], corner1[2]);
        glVertex3f(corner1[0], corner2[1], corner1[2]);
        glColor3f(0.4, 0.4, 0.4);
        glVertex3f(corner2[0], corner1[1], corner1[2]);
        glVertex3f(corner2[0], corner2[1], corner1[2]);
        glVertex3f(corner2[0], corner1[1], corner2[2]);
        glVertex3f(corner2[0], corner2[1], corner2[2]);
        glVertex3f(corner1[0], corner1[1], corner2[2]);
        glVertex3f(corner1[0], corner2[1], corner2[2]);
        glEnd();
    }

    if(!hasChildren){
        for(set<gcgParticleSPH*>::iterator it = particles.begin(); it != particles.end(); it++) {
            glPushMatrix();
//                glLoadIdentity();
            gcgParticleSPH* p = *it;
            glTranslatef(p->position[0],p->position[1], p->position[2]);
            glColor3f(p->color[0], p->color[1], p->color[2]);
            glutSolidSphere(p->radius, 10, 10);
            gcgDrawVectorPyramid(0., 0., 0., p->acceleration, 1.0f);
            p->text.textPosition(0.0, 0.0);
            if(print) p->text.gcgprintf("%d %d %d", p->acceleration[0], p->acceleration[1], p->acceleration[2]);
            print = !print;
            
            glPopMatrix();
//            if((p->position[0] < corner1[0]) || (p->position[0] > corner2[0]) || (p->position[1] < corner1[1]) || (p->position[1] > corner2[1]) || (p->position[2] < corner1[2]) || (p->position[2] > corner2[2])){
////                p->simulate = false;
//            }
        }

    }
    else if(hasChildren && (depth < MAX_OCTREE_DEPTH)){
        for(int x = 0; x < 2; x++){
            for(int y = 0; y < 2; y++){
                for(int z = 0; z < 2; z++){
                    children[x][y][z]->drawOctree();
                }
            }
        }
    }
}

//Adds a ball to or removes a from one to the children of this
void Octree::fileParticle(gcgParticleSPH * particle, VECTOR3 pos, bool addParticle) {
    //Figure out in which child(ren) the ball belongs
    for(int x = 0; x < 2; x++) {
        if (x == 0) {
            if (particle->position[0] - particle->radius > center[0]) {
                continue;
            }
        }
        else if (particle->position[0] + particle->radius < center[0]) {
            continue;
        }

        for(int y = 0; y < 2; y++) {
            if (y == 0) {
                if (particle->position[1] - particle->radius > center[1]) {
                    continue;
                }
            }
            else if (particle->position[1] + particle->radius < center[1]) {
                continue;
            }

            for(int z = 0; z < 2; z++) {
                if (z == 0) {
                    if (particle->position[2] - particle->radius > center[2]) {
                        continue;
                    }
                }
                else if (particle->position[2] + particle->radius < center[2]) {
                    continue;
                }

                //Add or remove the ball
                if (addParticle) {
                    children[x][y][z]->add(particle);
                }
                else {
                    children[x][y][z]->remove(particle, pos);
                }
            }
        }
    }
}

void Octree::haveChildren() {
    for(int x = 0; x < 2; x++) {
        float minX;
        float maxX;
        if (x == 0) {
            minX = corner1[0];
            maxX = center[0];
        }
        else {
            minX = center[0];
            maxX = corner2[0];
        }

        for(int y = 0; y < 2; y++) {
            float minY;
            float maxY;
            if (y == 0) {
                minY = corner1[1];
                maxY = center[1];
            }
            else {
                minY = center[1];
                maxY = corner2[1];
            }

            for(int z = 0; z < 2; z++) {
                float minZ;
                float maxZ;
                if (z == 0) {
                    minZ = corner1[2];
                    maxZ = center[2];
                }
                else {
                    minZ = center[2];
                    maxZ = corner2[2];
                }

                VECTOR3 c1, c2;
                gcgSETVECTOR3(c1, minX, minY, minZ);
                gcgSETVECTOR3(c2, maxX, maxY, maxZ);
                children[x][y][z] = new Octree(c1, c2, depth+1, this);
            }
        }
    }

    //Remove all balls from "balls" and add them to the new children
    for(set<gcgParticleSPH*>::iterator it = particles.begin(); it != particles.end();
            it++) {
        gcgParticleSPH* particle = *it;
        fileParticle(particle, particle->position, true);
    }
    particles.clear();

    hasChildren = true;
}

void Octree::destroyChildren() {
    for(int x = 0; x < 2; x++) {
        for(int y = 0; y < 2; y++) {
            for(int z = 0; z < 2; z++) {
                delete children[x][y][z];
            }
        }
    }

    hasChildren = false;
}

void Octree::collectParticles(set<gcgParticleSPH*> &bs) {
    if (hasChildren) {
            for(int x = 0; x < 2; x++) {
                    for(int y = 0; y < 2; y++) {
                            for(int z = 0; z < 2; z++) {
                                    children[x][y][z]->collectParticles(bs);
                            }
                    }
            }
    }
    else {
            for(set<gcgParticleSPH*>::iterator it = particles.begin(); it != particles.end(); it++) {
                    gcgParticleSPH* particle = *it;
                    bs.insert(particle);
            }
    }
}

void Octree::remove(gcgParticleSPH* particle, VECTOR3 pos) {
    this->numParticles--;

    if (hasChildren && numParticles < MIN_PARTICLES_PER_OCTREE) {
        destroyChildren();
    }

    if (hasChildren) {
        fileParticle(particle, pos, false);
    }
    else {
        particles.erase(particle);
    }
}

void Octree::checkCollision(){
    if(hasChildren){
        for(int x = 0; x < 2; x++)
            for(int y = 0; y < 2; y++)
                for(int z = 0; z < 2; z++)
                    children[x][y][z]->checkCollision();
    }
    else{
        set<gcgParticleSPH*>::iterator itEnd = particles.end();
        if(particles.size() > 1) itEnd--;
         for(set<gcgParticleSPH*>::iterator it = particles.begin(); it != itEnd; it++){
            gcgParticleSPH* particle = *it;
            if (particle->simulationStep == 0) continue;
            set<gcgParticleSPH*>::iterator nextElement = it;
            nextElement++;
            for(; nextElement != particles.end(); nextElement++){
                gcgParticleSPH* particle2 = *nextElement;
                if (particle2->simulationStep == 0) continue;
                float dist = sqrtf(pow(particle->position[0] - particle2->position[0],2) + pow(particle->position[1] - particle2->position[1],2) + pow(particle->position[2] - particle2->position[2],2));
                if(dist < (particle->radius + particle2->radius)) {particle->simulate = false; particle2->simulate = false;};//printf("Colidiu");
            }
         }
    }
}

double calcLennardJonesForce(double distance){
    double force;
    force = LJConst * (pow((potentialSigma/(distance)),14) - (0.5*pow((potentialSigma/(distance)),8)));
    
//    gcgSCALEVECTOR3(result[0], result[0], force);
    return force;
}

///laplaciano da gaussiana
double calcLaplaceGauss(double distance){
    double force;
    double dist2 = distance*distance;
    double dist4 = dist2*dist2;
//    derivate -(x^2-1)*e^(-((x^2)/2))
    force = -20*((dist2-1)*exp(-dist2/2));
//    force = -(dist2-1)*exp(-((dist2)/2));
//    force = exp(-((dist2)/2))*distance*(dist2-3);
//    force = 16.0/3.0*exp(-dist4/6.0)*distance*(dist4-dist2-3);
//    gcgSCALEVECTOR3(result[0], result[0], force);
    return force;
}
void Octree::checkCollision(gcgParticleSPH * particle, VECTOR3 * forceReturn){

    VECTOR3 LJForce;    
    gcgSETVECTOR3(LJForce, 0., 0., 0.);
    
    set<gcgParticleSPH*>::iterator itEnd = particles.end();
    
    
    for(set<gcgParticleSPH*>::iterator it = particles.begin(); it != itEnd; it++){
    gcgParticleSPH* particle2 = *it;

      if(particle == particle2) continue;
        
        VECTOR3 direction;//,force
        gcgSUBVECTOR3(direction, particle->position, particle2->position);
        float dist = sqrtf((direction[0]*direction[0]) + (direction[1]*direction[1]) * (direction[2]*direction[2]));
        if(dist <= 0.2) continue;
//        if((dist - particle2->radius) <= potentialRadius) {
//            vizNumb++;
            
            gcgNORMALIZEVECTOR3(direction, direction);
            
            gcgSCALEVECTOR3(direction, direction, calcLennardJonesForce(dist));
            gcgADDVECTOR3(LJForce, LJForce, direction);
    }
//    printf("Lennard Jones: %f %f %f\n", LJForce[0], LJForce[1], LJForce[2]);
    gcgCOPYVECTOR3((*forceReturn), LJForce);
    
}

void Octree::checkCollisionBruteForce(){

    
    
    set<gcgParticleSPH*>::iterator itEnd = particles.end();
    if(particles.size() > 1) itEnd--;

    for(set<gcgParticleSPH*>::iterator it = particles.begin(); it != itEnd; it++){
        gcgParticleSPH* particle = *it;
        set<gcgParticleSPH*>::iterator nextElement = it;
        nextElement++;
        for(; nextElement != particles.end(); nextElement++){
            gcgParticleSPH* particle2 = *nextElement;
            if(particle == particle2) continue;
            VECTOR3 direction;//,force
            gcgSUBVECTOR3(direction, particle->position, particle2->position);
            float dist = sqrtf((direction[0]*direction[0]) + (direction[1]*direction[1]) + (direction[2]*direction[2]));
//            if(dist <= 0.0) continue;
            gcgNORMALIZEVECTOR3(direction, direction);
//            gcgSCALEVECTOR3(direction, direction, calcLennardJonesForce(dist));
//            VECTOR3 diffVel;
//            gcgSUBVECTOR3(diffVel, particle2->velocity, particle->velocity);
//            gcgSCALEVECTOR3(direction, direction, calcLaplaceGauss(dist));
            gcgSCALEVECTOR3(direction, direction, calcLennardJonesForce(dist));
//            gcgADDVECTOR3(direction, direction, diffVel);
            gcgADDVECTOR3(particle->acceleration, particle->acceleration, direction);
//            gcgSCALEVECTOR3(direction, direction, ((1./2.0)*TIME_STEP2));
            gcgSCALEVECTOR3(direction, direction, -1.);
            gcgADDVECTOR3(particle2->acceleration, particle2->acceleration, direction);            

            
//            printf("Pos: %.16f %.16f %.16f\n", particle->position[0], particle->position[1], particle->position[2]);
//            printf("Pos2: %.16f %.16f %.16f\n", particle2->position[0], particle2->position[1], particle2->position[2]);
            printf("Dist: %f - Lennard Jones: %.16f %.16f %.16f\n", dist, direction[0], direction[1], direction[2]);
        }
    }
    

    
}

void Octree::particleMoved(gcgParticleSPH* particle, VECTOR3 oldPos) {
//    if(dbgcounter < 5){
//        printf("--------------------------------------\nDebug %d:\n", dbgcounter++);
//        printf("%f %f %f - ParticleMoved\n", oldPos[0], oldPos[1], oldPos[2]);
//    }
    //Verificar se sai do octante
    
    if((particle->position[0] < corner1[0]) || (particle->position[0] > corner2[0]) || (particle->position[1] < corner1[1]) || (particle->position[1] > corner2[1]) || (particle->position[2] < corner1[2]) || (particle->position[2] > corner2[2]))
    {
        VECTOR3 ppos;
        gcgCOPYVECTOR3(ppos, particle->position);
        gcgCOPYVECTOR3(particle->position, oldPos);
        remove(particle, particle->position);
        gcgCOPYVECTOR3(particle->position, ppos);
        add(particle);
    }
}

void Octree::add(gcgParticleSPH* particle) {
    numParticles++;
    if (!hasChildren && depth < MAX_OCTREE_DEPTH && numParticles > MAX_PARTICLES_PER_OCTREE) {
        haveChildren();
    }

    if (hasChildren) {
        fileParticle(particle, particle->position, true);
    }
    else {
        particles.insert(particle);
        particle->parentOctree = this;
    }
}

void Octree::remove(gcgParticleSPH * particle) {
        remove(particle, particle->position);
    }

Octree::Octree(VECTOR3 c1, VECTOR3 c2, int d, Octree * _parent) {
    gcgCOPYVECTOR3(corner1,c1);
    gcgCOPYVECTOR3(corner2,c2);
    VECTOR3 sum;
    gcgADDVECTOR3(sum, c1, c2);
    gcgSCALEVECTOR3(center, sum, (float)1/2);
    depth = d;
    numParticles = 0;
    hasChildren = false;
    this->parent = _parent;
    this->drawBox = drawOctreeBoxes;
}


Octree::~Octree()
{
    if (hasChildren) {
        destroyChildren();
    }
}

void toggleBoxingDraw(){
    drawOctreeBoxes = !drawOctreeBoxes;
}