/*
 * File:   octree.h
 * Author: joseluiz
 *
 * Created on April 25, 2012, 1:43 PM
 */

#ifndef OCTREE_H
#define	OCTREE_H

#include "gcg.h"
//#ifndef PARTICLE_H
#include "particle.h"
//#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <iostream>
#include <set>
#include <math.h>
using namespace std;

const int MAX_OCTREE_DEPTH = 7;
const int MIN_PARTICLES_PER_OCTREE = 2000;
const int MAX_PARTICLES_PER_OCTREE = 4500;
const float potentialRadius = 100.;
//#define potentialEpsilon 0.01
//#define potentialSigma 0.06
//#define potentialEpsilon 0.5
//#define potentialSigma 1.5
#define potentialEpsilon 0.5
#define potentialSigma 1.
#define LJConst ((potentialEpsilon*48.0)/(potentialSigma*potentialSigma))
//void calcLennardJonesForce(float distance, VECTOR3 * result);
/********
 Colocar ponteiro pra octree nas particuals
 *********/

//class gcgParticleSPH{
//public:
// VECTOR3 pos;
// float radius;
// VECTOR3 vel;
// VECTOR3 color;
// gcgParticleSPH();
//};

void toggleBoxingDraw();

class Octree {
public:
    Octree(VECTOR3 c1, VECTOR3 c2, int d, Octree * _parent);
    void add(gcgParticleSPH* particle);
    void particleMovedRealTime(gcgParticleSPH* particle, VECTOR3 oldPos);
    void particleMoved(gcgParticleSPH* particle, VECTOR3 oldPos);
    void buildTree();
    void drawOctree();
    void remove(gcgParticleSPH* particle, VECTOR3 pos);
    void checkCollision();
    void checkCollision(gcgParticleSPH * particle, VECTOR3 * forceReturn);
    void checkCollisionBruteForce();
    bool drawBox;
    virtual ~Octree();
private:
        void fileParticle(gcgParticleSPH * particle, VECTOR3 pos, bool addParticle);
        void haveChildren();
        void destroyChildren();
        void remove(gcgParticleSPH* particle);


        void collectParticles(set<gcgParticleSPH*> &bs);
        VECTOR3 corner1; //(minX, minY, minZ)
        VECTOR3 corner2; //(maxX, maxY, maxZ)
        VECTOR3 center;//((minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2)
        float radius;
        Octree * parent;
        Octree *children[2][2][2];
         //Whether this has children
        bool hasChildren;
        
        //The balls in this, if this doesn't have any children
//        set<Ball*> balls;
        //The depth of this in the tree
        int depth;
        //The number of balls in this, including those stored in its children
        int numParticles;

        //The balls in this, if this doesn't have any children
        set<gcgParticleSPH*> particles;

};


#endif	/* OCTREE_H */

