/*********************************************************************************************
***************** SMOOTHED PARTICLE HYDRODYNAMICS ********************************************
*** Tássio Knop de Castro
*** Gildo de Almeida Leonel
***********************************************************************************************/

#ifndef NEIGHBORSPH_H
#define NEIGHBORSPH_H
#include <vector>
#include "geometryInterface.h"
#include "particle.h"


using namespace std;

typedef int particle_index;

typedef struct
{
	void* prealloc;
	int used_size;
	int capacity;
} sph_mem_pool;

typedef struct
{
	int index;
	float distsq;
	float dist;
	VECTOR3 vdist;
} sph_neighbour;

typedef struct
{
	int* sizes;
	sph_neighbour** p; /** **/
	void* pool;
	int used_by_rigid; // Neighbour list memory used by rigid body grid
	int n_poolused; // Used size of memory pool
} sph_neighbour_list;


typedef struct
{
    float x;
    float y;
    float z;

} vector3;



class gcgOptimizedGridSPH: public gcgGEOMETRYINTERFACE {
  public:
    //	int** grid; /** grid indexed by grid **/
      particle_index** grid; /** particles indexed by grid **/
      gcgParticleSPH** particles; /** array of particles **/
      int* sizes;      /** number of particles grid has **/
    //	particle_index* sizes;      /** number of grid grid has **/
      int* caps;       /** capacity of each grid **/
      int capacity;
      int width;
      int height;
      int depth;
      float grid_len;  /** length of grid edge */

      float minx;
      float miny;
      float minz;

      bool taged;

      int nloops; /** amount of loops using the same neighbours list**/

      sph_mem_pool mempool;
      sph_neighbour_list n_list;
      ///vector3* pos;      // Position of grid (procurar armazenar as posicoes de uma forma sequencial porque deve melhorar)
      int nparticle;
      float search_radius;
	//int* particle_to_grid; /** particle index to grid index map **/

    void beginIteration();

    gcgPARTICLE *getFirstParticle(){};

    int getNumberOfNeighbors(gcgPARTICLE *particle){};

    gcgPARTICLE *getNeighbor(gcgPARTICLE *particle, int index){};

    gcgPARTICLE *getNextParticle(gcgPARTICLE *particle){};

    gcgPARTICLE *notifyUpdate(gcgPARTICLE *particle){};

    void endIteration(){};

    void sph_neighbour_list_create(int n_grid);
    void grid_create(int n_grid, float grid_len);
    void sph_grid_clear();
    void sph_grid_alloc();
    void sph_grid_get_neighbours();
    ///void init(int width, int height, int depth, int ngrid, float smoothlen, float search_radius, int nloops);
    void init(int width, int height, int depth, int ngrid, float smoothlen, float search_radius, int nloops);
    void addParticle(gcgParticleSPH* p);
    void recomputeDistances();
    void removeParticle(int p);
    void resetParticles();

};


#endif
