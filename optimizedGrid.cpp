/*********************************************************************************************
***************** SMOOTHED PARTICLE HYDRODYNAMICS ********************************************
*** T�ssio Knop de Castro
*** Gildo de Almeida Leonel
***********************************************************************************************/

#include "optimizedGrid.h"
#include "particle.h"
#include "gcg.h"
#include <cstdlib>
#include <stdlib.h>
#include <string.h>

#define EMPTYCELL -1
#define GRIDX(pos) ((int)((pos)[0]/influenceRadius))
#define GRIDY(pos) ((int)((pos)[1]/influenceRadius))
#define GRIDZ(pos) ((int)((pos)[2]/influenceRadius))

#define INDEX(x,y,z) ((z) * gheight * gwidth + (y) * gwidth + (x))
#define GRIDINDEX(pos) (INDEX(GRIDX(pos),GRIDY(pos), GRIDZ(pos)))

#define INITIAL_GRID_CAPACITY 60


///se der o erro sph_mem_pool_alloc: Out of memory pool size, basta aumentar o fator MEMPOOL_SCALE para resolver
#define MEMPOOL_SCALE 1

void sph_mem_pool_create(sph_mem_pool* pool, int size)
{
	pool->capacity = size;
	pool->prealloc = malloc(size);
	pool->used_size = 0;
}

void sph_mem_pool_free(sph_mem_pool* pool)
{
	pool->used_size = 0;
}


void* sph_mem_pool_alloc(sph_mem_pool* pool, int size)
{
	void* ret;
	if (size + pool->used_size > pool->capacity)
	{
		printf("sph_mem_pool_alloc: Out of memory pool size %d used %d\n", size, pool->used_size);
	}

	ret = &(((char*)pool->prealloc)[pool->used_size]);
	pool->used_size += size;

	return ret;
}

void gcgOptimizedGridSPH::grid_create(int n_grid, float grid_len)
{
	int i;
	int g_size;

	this->capacity = 40*40*40*MEMPOOL_SCALE;
	g_size = this->capacity;

//	g->sizes = (int*)malloc(g_size*sizeof(int));
	this->sizes = (int*)malloc(g_size*sizeof(particle_index));
	this->caps = (int*)malloc(g_size*sizeof(int));

	this->grid = (int**)malloc(g_size*sizeof(int*));

	//g->grid[0] = (int*)malloc(INITIAL_GRID_CAPACITY*g_size*sizeof(int));

	for (i = 1; i < g_size; i++)
	{
		//g->grid[i] = &(g->grid[0][i*INITIAL_GRID_CAPACITY]);
		this->caps[i] = 0;//INITIAL_GRID_CAPACITY;
	}

	this->grid_len = grid_len;

	sph_mem_pool_create(&this->mempool, sizeof(particle_index)*80*80*INITIAL_GRID_CAPACITY*MEMPOOL_SCALE);

//#ifdef EXP_SORT
//	indices = (int*)malloc(n_grid*sizeof(int));
//#endif
}

void gcgOptimizedGridSPH::sph_neighbour_list_create(int n_grid)
{
	n_list.sizes = (int*)malloc(n_grid*sizeof(int));
	///each particle cannot have more than 250 neighbours, but it's perfectly feasible.
	n_list.pool = (void*)malloc(250*n_grid*sizeof(sph_neighbour));
	//sph_mem_pool_create(&l->pool, 120*n_grid*sizeof(sph_neighbour));
	n_list.p = (sph_neighbour**)malloc(n_grid*sizeof(sph_neighbour*));
}



/**
 * Clear grids hash and resize grids from the distribution of grid.
 */
void gcgOptimizedGridSPH::sph_grid_clear()
{
	float fmin_x;
	float fmax_x;
	float fmin_y;
	float fmax_y;
	float fmin_z;
	float fmax_z;
	int i;

	fmin_x = fmax_x =particles[0]->position[0];
	fmin_y = fmax_y =particles[0]->position[1];
	fmin_z = fmax_z =particles[0]->position[2];

	for (i = 0; i < nparticle; i++)
	{
		float *p = particles[i]->position;
		if (fmin_x > p[0])
			fmin_x = p[0];
		if (fmax_x < p[0])
			fmax_x = p[0];
		if (fmin_y > p[1])
			fmin_y = p[1];
		if (fmax_y < p[1])
			fmax_y = p[1];
		if (fmin_z > p[2])
			fmin_z = p[2];
		if (fmax_z < p[2])
			fmax_z = p[2];
	}

	this->width  = (int)((fmax_x - fmin_x + this->grid_len)/this->grid_len);
	this->height = (int)((fmax_y - fmin_y + this->grid_len)/this->grid_len);
	this->depth  = (int)((fmax_z - fmin_z + this->grid_len)/this->grid_len);

	this->minx = fmin_x;
	this->miny = fmin_y;
	this->minz = fmin_z;

	if (this->width*this->height*this->depth > this->capacity)
	{
		printf("Out of Capacity\n");
	}

	memset(this->sizes, 0, this->width*this->height*this->depth*sizeof(particle_index));
	memset(this->caps, 0, this->width*this->height*this->depth*sizeof(int));
	sph_mem_pool_free(&this->mempool);
}


void gcgOptimizedGridSPH::sph_grid_get_neighbours()
{
	int i;
	int j;
	int neighbour_grid;
	int gx;
	int gy;
	int gz;
	int gindex; /* grid index */
	int pindex; /* particle index */
	float search_radius2;

	float inv_glen;

	//int n_poolused; /** Used size of memory pool of neighbour list **/

	n_list.n_poolused = 0;//n_list.used_by_rigid;
	search_radius2 = search_radius*search_radius;

	inv_glen = 1.0f/this->grid_len;

	for (i = 0; i < nparticle; i++)
	{
		n_list.p[i] = &((sph_neighbour*)n_list.pool)[n_list.n_poolused];
		n_list.p[i][0].index = i;
		n_list.p[i][0].distsq = 0;
		n_list.p[i][0].dist = 0;
		n_list.sizes[i] = 1;

		float *pos = particles[i]->position;
		gx = (int)((pos[0] - this->minx)*inv_glen);
		gy = (int)((pos[1] - this->miny)*inv_glen);
		gz = (int)((pos[2] - this->minz)*inv_glen);

		gindex = gx + gy*this->width + gz*this->width*this->height;

		for (gz = -1; gz <= 1; gz++)
			for (gy = -1; gy <= 1; gy++)
				for (gx = -1; gx <= 1; gx++)
				{
					neighbour_grid = gindex + this->width*this->height*gz + this->width*gy + gx;
					if ((neighbour_grid < 0) || (neighbour_grid >= this->width*this->depth*this->height))
						continue;

					for (j = 0; j < this->sizes[neighbour_grid]; j++)
					{


						pindex = this->grid[neighbour_grid][j];
						//////distsq = vec3_distsq(&pos[i], &pos[pindex]);

            VECTOR3 dist;
            gcgSUBVECTOR3(dist, particles[i]->position, particles[pindex]->position);
            float distsq = gcgDOTVECTOR3(dist,dist);;

						if (distsq < search_radius2)
						{
							int nindex = n_list.sizes[i];
							n_list.p[i][nindex].index = pindex;
							gcgCOPYVECTOR3(n_list.p[i][nindex].vdist,dist);
							n_list.p[i][nindex].distsq = distsq;
							n_list.p[i][nindex].dist = sqrt(distsq);
							n_list.sizes[i]++;
						}
					}
				}

		if (this->sizes[gindex] == this->caps[gindex])
		{
//			this->grid[gindex] = (particle_index*)sph_mem_pool_alloc(&this->mempool, 50*sizeof(particle_index));
			this->grid[gindex] = (particle_index*)sph_mem_pool_alloc(&this->mempool, 150*sizeof(particle_index));
			this->caps[gindex] = 150;
			//gput_printf("sph_grid_get_neighbours: Out of Grid Capacity\n");
			//exit(1);
		}

		n_list.n_poolused += n_list.sizes[i];
		this->grid[gindex][this->sizes[gindex]] = i;
		this->sizes[gindex]++;
	}

	//gput_printf("Num. of interactions = %d\n", n_list.n_poolused);
}

void gcgOptimizedGridSPH::init(int width, int height, int depth, int nparticle, float smoothlen, float search_radius, int nloops){
  this->nloops = nloops;
  this->nparticle = nparticle;
  this->search_radius = search_radius;
  this->width = width;
  this->height = height;
  this->depth = depth;
  this->particles = (gcgParticleSPH**) malloc(sizeof(gcgParticleSPH*) * nparticle);
  this->taged = false;

//  sph_neighbour_list_create(nparticle);
//  grid_create(nparticle, smoothlen);

  ///sph_grid_clear();
  ///sph_grid_get_neighbours();

}


/**
 * Allocate grid to grids
 */
void gcgOptimizedGridSPH::sph_grid_alloc()
{
	int i;

	float inv_glen = 1.0f/grid_len;

	for (i = 0; i < nparticle; i++)
	{
		int gx;
		int gy;
		int gz;
		int gindex;
/*
		n_list->p[i] = &((sph_neighbour*)n_list->pool)[n_poolused];
		n_list->p[i][0].index = i;
		n_list->p[i][0].distsq = 0;
		n_list->sizes[i] = 1;
		*/
		float *pos = particles[i]->position;
		gx = (int)((pos[0] - minx) *inv_glen);
		gy = (int)((pos[1] - miny) *inv_glen);
		gz = (int)((pos[2] - minz) *inv_glen);

		gindex = gx + gy*width + gz*width*height;

		if (sizes[gindex] == caps[gindex])
		{
			grid[gindex] = (particle_index*)sph_mem_pool_alloc(&mempool, 50*sizeof(particle_index));
			caps[gindex] = 50;
			//gput_printf("sph_grid_get_neighbours: Out of Grid Capacity\n");
			//exit(1);
		}

//		gput_printf("%d %d %d\n", gx, gy, gz);

		grid[gindex][sizes[gindex]] = i;
		sizes[gindex]++;
	}
}



void gcgOptimizedGridSPH::recomputeDistances(){
  for (int i = 0; i < nparticle; i++){
    for (int j = 0; j < n_list.sizes[i]; j++){
      sph_neighbour* n = &n_list.p[i][j];
      gcgSUBVECTOR3(n->vdist, particles[i]->position, particles[n->index]->position);
      n->distsq = gcgDOTVECTOR3(n->vdist,n->vdist);
      n->dist = sqrt(n->distsq);
    }
  }
}


void gcgOptimizedGridSPH::beginIteration(){
  static int loop = 0;
  if (!loop){ ///if loop == 0
//    sph_mem_pool_free(&this->mempool);
    sph_grid_clear();
    sph_grid_alloc();
    sph_grid_get_neighbours();
  }
  else{
    recomputeDistances();
  }
  ++loop;
  if (loop == 2) loop = 0;

  };

void gcgOptimizedGridSPH::addParticle(gcgParticleSPH* p){
  static int counter = 0;
  particles[counter] = p;
  ++counter;
  this->nparticle = counter;
};


void gcgOptimizedGridSPH::resetParticles(){
 for(unsigned int i = 0; i < nparticle; i++){
   gcgCOPYVECTOR3(particles[i]->position, particles[i]->startPos);
//   gcgSETVECTOR3(particles[i]->color,0.0,0.0,1.0);
 }
};

void gcgOptimizedGridSPH::removeParticle(int p){
  for (int i = p; i < this->nparticle; i++ )
    particles[i] = NULL;
  
  this->nparticle = this->nparticle - 1;
};
