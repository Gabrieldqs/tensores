#include "tensorgenerate.h"

gcgGCGTENSORGENERATE::gcgGCGTENSORGENERATE(int width, int height, int depth ){

    grid = new gcgTENSORGRID();
    grid->createGrid(width, height, depth);
    VECTOR3 v, v1, v2;
    MATRIX3 t;
    v1[0] = 2.0;
    v1[1] = 2.0;
    v1[2] = 2.0;

   for (unsigned int z = 0; z < depth; z++)
     for (unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++){
            gcgSETVECTOR3(v,x,y,z);
            gcgSCALEVECTOR3(v2,v,sqrt((x*x)+(y*y)+(z*z)));
            gcgNORMALIZEVECTOR3(v1,v2);
            gcgVECTOR3TOTENSOR(t,v1,v1);
           gcgCOPYVECTOR3(grid->data[z* width * height + y*width+ x],t) ;
    }
    

//   for (unsigned int z = 0; z < depth; z++)
//     for (unsigned int y = 0; y < height; y++)
//        for (unsigned int x = 0; x < width; x++){
//            gcgSETVECTOR3(v,x,y,z);
//            gcgSUBVECTOR3(v2,v1,v);
//            gcgSUBVECTOR3(v2,v2,v);
//            gcgSCALEVECTOR3(v1,v1, 1 / sqrt(((v1[0]-v[0])*(v1[0]-v[0]))+((v1[1]-v[1])*(v1[1]-v[1]))+((v1[2]-v[2])*(v1[2]-v[2]))));
//            gcgSCALEVECTOR3(v2,v2, 1 / sqrt(((v2[0]-v[0])*(v2[0]-v[0]))+((v2[1]-v[1])*(v2[1]-v[1]))+((v2[2]-v[2])*(v2[2]-v[2]))));
//            gcgADDVECTOR3(v,v1,v2);
//            gcgNORMALIZEVECTOR3(v,v);
//            gcgVECTOR3TOTENSOR(t,v,v);
//            gcgCOPYVECTOR3(grid->data[z* width * height + y*width+ x],t) ;
//    }
}

gcgGCGTENSORGENERATE::~gcgGCGTENSORGENERATE() {
//    grid->destroyGrid;
}
