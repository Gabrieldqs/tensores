/* 
 * File:   newGrid.h
 * Author: joseluiz
 *
 * Created on November 21, 2012, 1:30 PM
 */

#ifndef NEWGRID_H
#define	NEWGRID_H

#include "gcg.h"
#include "particle.h"
#include <math.h>
#include <cmath>
#include <stdlib.h>

class gridCell{
public:
    vector<gcgParticleSPH*> particles;
    gridCell();
    virtual ~gridCell();
};

class newGrid {
private:
    bool is2D;
    float width, height, depth;
    float cellRadius;
//    gridCell *** cell;
    gridCell ** cell;
    int totalCells;
public:
    void addParticle(gcgParticleSPH * p);
    void moveParticle(gcgParticleSPH * p);
    void findNeighbours(gcgParticleSPH * p);
    int hasNeighbours(VECTOR3 pos);

    newGrid();
    newGrid(float w, float h, float d, float radius, bool m2D);
    virtual ~newGrid();
private:

};

#endif	/* NEWGRID_H */

