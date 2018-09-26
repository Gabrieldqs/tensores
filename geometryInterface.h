/*********************************************************************************************
***************** SMOOTHED PARTICLE HYDRODYNAMICS ********************************************
*** T�ssio Knop de Castro
*** Gildo de Almeida Leonel
***********************************************************************************************/

#ifndef GEOMETRYINTERFACE_H
#define GEOMETRYINTERFACE_H

#include "gcg.h"

class gcgPARTICLE {
public:
    VECTOR3 startPos;
    int index; //for use of gcgGEOMETRYINTERFACE
};

class gcgGEOMETRYINTERFACE {
public:

    virtual void beginIteration() = 0;

    virtual gcgPARTICLE *getFirstParticle() = 0;

    virtual int getNumberOfNeighbors(gcgPARTICLE *particle) = 0;

    virtual gcgPARTICLE *getNeighbor(gcgPARTICLE *particle, int index) = 0;

    virtual gcgPARTICLE *getNextParticle(gcgPARTICLE *particle) = 0;

    virtual gcgPARTICLE *notifyUpdate(gcgPARTICLE *particle) = 0;

    virtual void endIteration() = 0;

};



#endif
