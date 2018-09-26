/* 
 * File:   newGrid.cpp
 * Author: joseluiz
 * 
 * Created on November 21, 2012, 1:30 PM
 */

#include "newGrid.h"

int cellsX, cellsY, cellsZ;

#define ADDRESSGRID(x, y, z)  (x) + (y*cellsX) + (z*cellsX*cellsY)

gridCell::gridCell(){
    
}

gridCell::~gridCell(){
    
}

newGrid::newGrid() {
    
}



newGrid::newGrid(float w, float h, float d, float radius, bool m2D) {
    
    width = w;
    height = h;
    depth = d;
    cellRadius = radius;
    cellsX = (int)(ceil(w/radius));
    cellsY = (int)(ceil(h/radius));
    is2D = m2D;
    if(is2D) {
        cellsZ = 0;
    }
    else {
        cellsZ = (int)(ceil(d/radius));
    }
    
//   cell = (gridCell ***) malloc (cellsX * sizeof(gridCell**));
    if(is2D) totalCells = cellsX * cellsY;
    else totalCells = cellsX * cellsY * cellsZ;
    
   cell = (gridCell **) malloc (totalCells * sizeof(gridCell*));
 
   for(int i = 0; i < totalCells; i++){
       cell[i] = (gridCell*)malloc(sizeof(gridCell));
   }
   
   
   
//   cell = (gridCell ***) malloc (cellsX * sizeof(gridCell**));
//   
//   for(int i = 0; i < cellsX; i++){
//       cell[i] = (gridCell **) malloc
//   }
   
    
}

void newGrid::addParticle(gcgParticleSPH* p){

    int v = ADDRESSGRID(floor(p->position[0]/cellRadius), floor(p->position[1]/cellRadius), floor(p->position[2]/cellRadius));
    p->gridAdreess = v;
    cell[v]->particles.push_back(p);
    
}

void newGrid::moveParticle(gcgParticleSPH* p){
    int v = ADDRESSGRID(floor(p->position[0]/cellRadius), floor(p->position[1]/cellRadius), floor(p->position[2]/cellRadius));
    
    if(v == p->gridAdreess) return;
    
//    printf("Procurando!\n");
    int vectorIndex = 0, oldV = p->gridAdreess;
    
    for(std::vector<gcgParticleSPH*>::iterator it = cell[oldV]->particles.begin(); it != cell[oldV]->particles.end(); ++it){
        if(*it == p){
//            printf("Encontrado!\n");
            cell[oldV]->particles.erase(cell[oldV]->particles.begin() + vectorIndex);
            break;
        }
        vectorIndex++;
    }
    
    p->gridAdreess = v;
    cell[v]->particles.push_back(p);
}


void newGrid::findNeighbours(gcgParticleSPH* p){
    
    p->neighbourCheck = true;
    
    int v;
    float zx, yx, xx;
    for(int z = - 1; z <= +1; z++){
        
        zx = (((float)z) * cellRadius) + p->position[2];
        if(zx < 0 || zx > depth) continue;
        
        
        for(int y = - 1; y <= +1; y++){
            
            yx = (((float)y) * cellRadius) + p->position[1];
            if(yx < 0 || yx > height) continue;
            
            for(int x = - 1; x <= +1; x++){
                
                xx = (((float)x) * cellRadius) + p->position[0];
                if(xx < 0 || xx > width) continue;
                
                v = ADDRESSGRID(floor(xx/cellRadius), floor(yx/cellRadius), floor(zx/cellRadius));
                
                for(std::vector<gcgParticleSPH*>::iterator it = cell[v]->particles.begin(); it != cell[v]->particles.end(); ++it){
                    gcgParticleSPH * p2 = *it;
                    if(!p2->neighbourCheck){
                        VECTOR3 diff;
                        gcgSUBVECTOR3(diff, p->position, p2->position);
                        float dist = gcgLENGTHVECTOR3(diff);
                        
                        if(dist <= cellRadius && !FEQUAL(dist, 0.0)){
//                    printf("Adicionando P: %p - %d | P2: %p - %d\n", p, p->neighbours.size(), p2, p2->neighbours.size());
                            p->neighbours.push_back(p2);
                            p2->neighbours.push_back(p);
                        }
                    }
                }
            }
        }
    }
}


int newGrid::hasNeighbours(VECTOR3 pos){
    
    
    int v;
    float zx, yx, xx;
    int s = 0;
    for(int z = - 1; z <= +1; z++){
        
        zx = (((float)z) * cellRadius) + pos[2];
        if(zx < 0 || zx > depth) continue;
        
        
        for(int y = - 1; y <= +1; y++){
            
            yx = (((float)y) * cellRadius) + pos[1];
            if(yx < 0 || yx > height) continue;
            
            for(int x = - 1; x <= +1; x++){
                
                xx = (((float)x) * cellRadius) + pos[0];
                if(xx < 0 || xx > width) continue;
                
                v = ADDRESSGRID(floor(xx/cellRadius), floor(yx/cellRadius), floor(zx/cellRadius));
                
                for(std::vector<gcgParticleSPH*>::iterator it = cell[v]->particles.begin(); it != cell[v]->particles.end(); ++it){
                    gcgParticleSPH * p2 = *it;
                        VECTOR3 diff;
                        gcgSUBVECTOR3(diff, pos, p2->position);
                        float dist = gcgLENGTHVECTOR3(diff);
                        
                        if(dist <= cellRadius && !FEQUAL(dist, 0.0)){
                            s++;
                        }
                }
            }
        }
    }
    if(s > 10) return 1;
    return 0;
}

newGrid::~newGrid() {
}

