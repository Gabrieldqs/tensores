/*
 * File:   gcgTensor.cpp
 * Author: joseluiz
 *
 * Created on October 11, 2012, 10:48 AM
 */

#include "gcgTensor.h"

#define THRESHOLD  0.2
int tipoCampo = HELICE;
//int tipoCampo = HELICE1;
//int tipoCampo = FIELD2D;
//int tipoCampo = CEREBRO;
//int tipoCampo = CEREBRO2;
//int tipoCampo = CEREBRO1;
int normShow = 5;



/** Retirar na limpeza **/
float initialDensity, totalSpace, mass, radius, initialDistance, smoothWidth, gridInfluenceRadius;
int totalSize;



float gcgTensor::calculateNorm(MATRIX3 t){

    return sqrtf((t[0]*t[0]) + (t[1]*t[1]) + (t[2]*t[2]) + (t[3]*t[3]) + (t[4]*t[4]) + (t[5]*t[5]) + (t[6]*t[6]) + (t[7]*t[7]) + (t[8]*t[8]));

}

int gcgTensor::createTensorlineDup(VECTOR3* vField, gcgTENSORGLYPH** glyphs, int *nglyphs, unsigned int width, unsigned int height, unsigned int depth, int seedI, int seedJ, int seedK, float bigger, int i /* = 0 */) {

    VECTOR3 vIn, vOut, vProp, v1, vAux, nextPos;
    MATRIX2 scaleD;
    int adressG, posX, posY, posZ, adressStart;
    float wPunct = 0.9;

    adressStart = ADDRESSGLYPHS(seedI, seedJ, seedK, width, height);
    adressG = adressStart;
    float dt = 0.1;
    if (glyphs[adressG]->tnumber == 0) {


        //FOWARD
        adressG = adressStart;
        gcgSETVECTOR3(vIn, 0.0, 0.0, 0.0);
        gcgSCALEMATRIX3(scaleD, glyphs[adressG]->tensor, (2.0f / bigger));
        gcgAPPLYMATRIX3VECTOR3(vOut, scaleD, vIn);


        gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[0], glyphs[adressG]->eigenVectors[1], glyphs[adressG]->eigenVectors[2])


        gcgSCALEVECTOR3(v1, v1, glyphs[adressG]->cl);
        gcgSCALEVECTOR3(vIn, vIn, 1.0 - wPunct);
        gcgSCALEVECTOR3(vOut, vOut, wPunct);
        gcgADDVECTOR3(vAux, vIn, vOut);
        gcgSCALEVECTOR3(vAux, vAux, 1.0 - glyphs[adressG]->cl);
        gcgADDVECTOR3(vProp, vAux, v1);


        if (gcgLENGTHVECTOR3(vField[adressG]) < gcgLENGTHVECTOR3(vProp))
            gcgCOPYVECTOR3(vField[adressG], vProp);



        //assumindo dt = 1
        gcgSCALEVECTOR3(vProp, vProp, dt);
        gcgADDVECTOR3(nextPos, glyphs[adressG]->pos, vProp);


        posX = (unsigned int) nextPos[0];// * (width - 1)) / (unsigned int) fwidth;
        posY = (unsigned int) nextPos[1];// * (height - 1)) / (unsigned int) fheight;
        posZ = (unsigned int) nextPos[2];// * (depth - 1)) / (unsigned int) fdepth;


        int aux = 1;

        while ((posX < (unsigned int) fwidth) && (posY < (unsigned int) fheight) && (posZ < (unsigned int) fdepth) && (posX > 0) && (posY > 0) && (posZ > 0) && (adressG < *nglyphs) && (adressG != ADDRESSGLYPHS(posX, posY, posZ, width, height)) && (gcgLENGTHVECTOR3( vField[ADDRESSGLYPHS(posX, posY, posZ, width, height)]) > 0.0) ){

            adressG = ADDRESSGLYPHS(posX, posY, posZ, width, height);
            gcgCOPYVECTOR3(vIn, vProp);
            if (gcgDOTVECTOR3(vIn, v1) < 0) gcgSCALEVECTOR3(v1, v1, -1);

            gcgAPPLYMATRIX3VECTOR3(vOut, glyphs[adressG]->tensor, vIn);

            if ((glyphs[adressG]->eigenValues[0] > glyphs[adressG]->eigenValues[1]) && (glyphs[adressG]->eigenValues[0] > glyphs[adressG]->eigenValues[2]))
                gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[0], glyphs[adressG]->eigenVectors[1], glyphs[adressG]->eigenVectors[2])
            else if ((glyphs[adressG]->eigenValues[1] > glyphs[adressG]->eigenValues[0]) && (glyphs[adressG]->eigenValues[1] > glyphs[adressG]->eigenValues[2]))
                gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[3], glyphs[adressG]->eigenVectors[4], glyphs[adressG]->eigenVectors[5])
            else if ((glyphs[adressG]->eigenValues[2] > glyphs[adressG]->eigenValues[0]) && (glyphs[adressG]->eigenValues[2] > glyphs[adressG]->eigenValues[1]))
                gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[6], glyphs[adressG]->eigenVectors[7], glyphs[adressG]->eigenVectors[8]);

            gcgSCALEVECTOR3(v1, v1, glyphs[adressG]->cl);
            gcgSCALEVECTOR3(vIn, vIn, 1.0 - wPunct);
            gcgSCALEVECTOR3(vOut, vOut, wPunct);
            gcgADDVECTOR3(vAux, vIn, vOut);
            gcgSCALEVECTOR3(vAux, vAux, 1.0 - glyphs[adressG]->cl);
            gcgADDVECTOR3(vProp, vAux, v1);


            //            if (gcgLENGTHVECTOR3(vField[adressG]) < gcgLENGTHVECTOR3(vProp))
            gcgCOPYVECTOR3(vField[adressG], vProp);

            //assumindo dt = 1
            gcgSCALEVECTOR3(vProp, vProp, dt);
            gcgADDVECTOR3(nextPos, glyphs[adressG]->pos, vProp);
            glyphs[adressG]->tnumber = adressStart;
            posX = ((unsigned int) nextPos[0] * (width - 1)) / (unsigned int) fwidth;
            posY = ((unsigned int) nextPos[1] * (height - 1)) / (unsigned int) fheight;
            posZ = ((unsigned int) nextPos[2] * (depth - 1)) / (unsigned int) fdepth;

            aux++;

        }

//        return aux;



        //BACKWARD

        adressG = adressStart;
        gcgSETVECTOR3(vIn, 0.0, 0.0, 0.0);
        gcgSCALEMATRIX3(scaleD, glyphs[adressG]->tensor, (2.0f / bigger));
        gcgAPPLYMATRIX3VECTOR3(vOut, scaleD, vIn);


        gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[0], glyphs[adressG]->eigenVectors[1], glyphs[adressG]->eigenVectors[2])

        gcgSCALEVECTOR3(v1, v1, glyphs[adressG]->cl);
        gcgSCALEVECTOR3(vIn, vIn, 1.0 - wPunct);
        gcgSCALEVECTOR3(vOut, vOut, wPunct);
        gcgADDVECTOR3(vAux, vIn, vOut);
        gcgSCALEVECTOR3(vAux, vAux, 1.0 - glyphs[adressG]->cl);
        gcgADDVECTOR3(vProp, vAux, v1);
        gcgSCALEVECTOR3(vProp, vProp, -1.0);

        //assumindo dt = 1
            gcgSCALEVECTOR3(vProp, vProp, dt);
        gcgADDVECTOR3(nextPos, glyphs[adressG]->pos, vProp);


        posX = (unsigned int) nextPos[0];// * (width - 1)) / (unsigned int) fwidth;
        posY = (unsigned int) nextPos[1];// * (height - 1)) / (unsigned int) fheight;
        posZ = (unsigned int) nextPos[2];// * (depth - 1)) / (unsigned int) fdepth;

        while ((posX < (unsigned int) fwidth) && (posY < (unsigned int) fheight) && (posZ < (unsigned int) fdepth) && (posX > 0) && (posY > 0) && (posZ > 0) && (adressG < *nglyphs) && (adressG != ADDRESSGLYPHS(posX, posY, posZ, width, height))&& (gcgLENGTHVECTOR3( vField[ADDRESSGLYPHS(posX, posY, posZ, width, height)]) > 0.0) ){


            adressG = ADDRESSGLYPHS(posX, posY, posZ, width, height);
            gcgCOPYVECTOR3(vIn, vProp);
            if (gcgDOTVECTOR3(vIn, v1) >= 0) gcgSCALEVECTOR3(v1, v1, -1);

            gcgAPPLYMATRIX3VECTOR3(vOut, glyphs[adressG]->tensor, vIn);

            if ((glyphs[adressG]->eigenValues[0] > glyphs[adressG]->eigenValues[1]) && (glyphs[adressG]->eigenValues[0] > glyphs[adressG]->eigenValues[2]))
                gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[0], glyphs[adressG]->eigenVectors[1], glyphs[adressG]->eigenVectors[2])
            else if ((glyphs[adressG]->eigenValues[1] > glyphs[adressG]->eigenValues[0]) && (glyphs[adressG]->eigenValues[1] > glyphs[adressG]->eigenValues[2]))
                gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[3], glyphs[adressG]->eigenVectors[4], glyphs[adressG]->eigenVectors[5])
            else if ((glyphs[adressG]->eigenValues[2] > glyphs[adressG]->eigenValues[0]) && (glyphs[adressG]->eigenValues[2] > glyphs[adressG]->eigenValues[1]))
                gcgSETVECTOR3(v1, glyphs[adressG]->eigenVectors[6], glyphs[adressG]->eigenVectors[7], glyphs[adressG]->eigenVectors[8]);

            gcgSCALEVECTOR3(v1, v1, glyphs[adressG]->cl);
            gcgSCALEVECTOR3(vIn, vIn, 1.0 - wPunct);
            gcgSCALEVECTOR3(vOut, vOut, wPunct);
            gcgADDVECTOR3(vAux, vIn, vOut);
            gcgSCALEVECTOR3(vAux, vAux, 1.0 - glyphs[adressG]->cl);
            gcgADDVECTOR3(vProp, vAux, v1);


            //            if (gcgLENGTHVECTOR3(vField[adressG]) < gcgLENGTHVECTOR3(vProp))
            gcgCOPYVECTOR3(vField[adressG], vProp);

            //assumindo dt = 1
            gcgSCALEVECTOR3(vProp, vProp, dt);
            gcgADDVECTOR3(nextPos, glyphs[adressG]->pos, vProp);
            glyphs[adressG]->tnumber = adressStart;

            posX = (unsigned int) nextPos[0];// * (width - 1)) / (unsigned int) fwidth;
            posY = (unsigned int) nextPos[1];// * (height - 1)) / (unsigned int) fheight;
            posZ = (unsigned int) nextPos[2];// * (depth - 1)) / (unsigned int) fdepth;

            aux++;

        }



        return aux;
    }

    return 0;
}

void gcgTensor::setGlobalVars(){

    initialDensity = 10;
    totalSpace = (float) fwidth * (float) fheight * (float) fdepth;
    totalSize = fwidth * fheight * fdepth;
    mass = initialDensity * totalSpace / (nParticles * 10.5);
    radius = pow(0.75 * mass / (initialDensity * M_PI), 1.0 / 3.0);

    initialDistance = 3.0 * radius;
    smoothWidth = 1.3 * initialDistance;
    gridInfluenceRadius = smoothWidth * 2; ///Influence factor = 2

}
/**************************/

inline void endian_swap(unsigned int *x) {
    *x = (*x >> 24) |
            ((*x << 8) & 0x00FF0000) |
            ((*x >> 8) & 0x0000FF00) |
            (*x << 24);
}

void gcgTensor::readNewField(int tamanho, const char* fileName, MATRIX3 *field) {

    FILE* file;
    puts("teste1");
    if (!(file = fopen(fileName, "rb+"))) puts("Erro ao importar campo tensorial!");

//    fscanf(file, "%d %d %d\n\n", &fwidth, &fheight, &fdepth);
    //stores the tensor field
    for (int i = 0; i < tamanho; ++i) {
        fread(&field[i], 1,  9 * sizeof(float), file);
//        fscanf(file, "%f %f %f\n%f %f %f\n%f %f %f\n", &field[i][0], &field[i][1], &field[i][2], &field[i][3], &field[i][4], &field[i][5], &field[i][6], &field[i][7], &field[i][8]);
        if (feof(file)) break;
    }
    fclose(file);
}

void gcgTensor::readTXTField(int tamanho, const char* fileName, MATRIX3 *field) {

    FILE* file;
     puts("teste2");
    if (!(file = fopen(fileName, "rb+"))) puts("Erro ao importar campo tensorial!");

    fscanf(file, "%d %d %d\n\n", &fwidth, &fheight, &fdepth);
    //stores the tensor field
    for (int i = 0; i < tamanho; ++i) {
        fscanf(file, "%f %f %f\n%f %f %f\n%f %f %f\n", &field[i][0], &field[i][1], &field[i][2], &field[i][3], &field[i][4], &field[i][5], &field[i][6], &field[i][7], &field[i][8]);
        if (feof(file)) break;
    }
    fclose(file);
}


void gcgTensor::readTXTField3(int a, const char* fileName, MATRIX3 * field, int tamanho) {

    FILE* file;
     puts("teste3");
    if (!(file = fopen(fileName, "rb"))) puts("Erro ao importar campo tensorial! txtf3");
    int w, h, d;
    fread(&w, sizeof(int), 1, file);
    fread(&h, sizeof(int), 1, file);
    fread(&d, sizeof(int), 1, file);

    
//    printf("%d %d %d", w, h, d);
//    exit(0);

//    printf("Tamanho %s: %d\n", fileName, tamanho);
    for (int i = 0; i < tamanho; ++i) {
        for(int k = 0; k < 9; k++){
           fread(&field[i][k], sizeof(float), 1, file) ;
        }

        if (feof(file)) break;
    }
    fclose(file);
}

void gcgTensor::storeField(const char* fileName, MATRIX3 * field) {

    FILE* file;
    int tSize = fwidth*fheight*fdepth;
    if (!(file = fopen(fileName, "wb"))) puts("Erro ao criar arquivo!");
    fwrite(&fwidth, sizeof(int), 1, file);
    fwrite(&fheight, sizeof(int), 1, file);
    fwrite(&fdepth, sizeof(int), 1, file);

//    printf("%d %d %d", w, h, d);
//    exit(0);

//    printf("Tamanho %s: %d\n", fileName, tamanho);
    for (int i = 0; i < tSize; ++i) {
        for(int k = 0; k < 9; k++){
           fwrite(&field[i][k], sizeof(float), 1, file) ;
        }

        if (feof(file)) break;
    }
    fclose(file);
}

void gcgTensor::readRawField(int tamanho, const char* fileName, MATRIX3 *field, float *confidence) {

    FILE* file;
    FILE* fdb = fopen("dbg2.txt", "w");
    if (!(file = fopen(fileName, "rb+"))) puts("Erro ao importar campo tensorial!");
    //  if(!(file = fopen("gk2-rcc-mask.raw", "rb+"))) puts("Erro ao importar campo tensorial!");

    // stores the tensor field
    for (int i = 0; i < tamanho; ++i) {
        fread(&confidence[i], sizeof (float), 1, file);
        fread(&field[i][0], sizeof (float) + sizeof (float) + sizeof (float), 1, file);
        fread(&field[i][4], sizeof (float) + sizeof (float), 1, file);
        fread(&field[i][8], sizeof (float), 1, file);
        field[i][3] = field[i][1];
        field[i][6] = field[i][2];
        field[i][7] = field[i][5];

        unsigned int *dummy = (unsigned int*) &confidence[i];
        endian_swap(&dummy[0]);
        // printf("%f \n",confidence[i]);

        dummy = (unsigned int*) &field[i][0];

        
        
        endian_swap(&dummy[0]);
        endian_swap(&dummy[1]);
        endian_swap(&dummy[2]);
        endian_swap(&dummy[3]);
        endian_swap(&dummy[4]);
        endian_swap(&dummy[5]);
        endian_swap(&dummy[6]);
        endian_swap(&dummy[7]);
        endian_swap(&dummy[8]);
        
        fprintf(fdb, "\n%f %f %f\n%f %f %f\n%f %f %f\n", field[i][0], field[i][1], field[i][2], field[i][3], field[i][4], field[i][5], field[i][6], field[i][7], field[i][8]);
        
//        fprintf(fdb, "\n%f %f %f\n%f %f %f\n%f %f %f\n", field[i][1], field[i][2], field[i][3], field[i][4], field[i][5], field[i][6], field[i][7], field[i][8], field[i][9]);

    }
    fclose(fdb);
    fclose(file);
}

gcgTensor::gcgTensor() {
}

gcgTensor::gcgTensor(int width, int height, int depth, int nPart) {

    setGlobalVars();

    //inicio da etapa de codificação
     //inicio da etapa de codificação
    fwidth = width;
    fheight = height;
    fdepth = depth;
    nParticles = nPart;
    totalSize = (fwidth) * (fheight) * (fdepth);

    velField = (VECTOR3*) malloc(sizeof (VECTOR3) * (int) (totalSize));

    glyphs_original = (gcgTENSORGLYPH**) malloc(sizeof (gcgTENSORGLYPH*) * (int) totalSize);

    tensor_ord = (gcgTENSORGLYPH**) malloc(sizeof (gcgTENSORGLYPH*) * (totalSize));

    geometry = new gcgOptimizedGridSPH();
    geometry->init((int) fwidth, (int) fheight, (int) fdepth, nParticles, smoothWidth, gridInfluenceRadius, 1);
//    geometry->init((int) fwidth, (int) fheight, (int) fdepth);

    //Creating the tensor grid
    gcgTENSORGRID *grid_entrada = new gcgTENSORGRID();
    grid_entrada->createGrid(fwidth, fheight, fdepth);

    //Opening the input dt-mri file and creating space for the tensor field
    field = (MATRIX3*) malloc(sizeof (MATRIX3) * totalSize);
    confidence = (float*) malloc(sizeof (float) * totalSize);
    norms  = (float*)malloc(sizeof(float)*totalSize);
    normDerivatives  = (VECTOR3*)malloc(sizeof(VECTOR3)*totalSize);

    maxNorm = -9999999999;
    maxDerivN = -9999999999;
    minNorm = 999999999;
    minDerivN = 999999999;
    
    float zVal;
    if(tipoCampo == FIELD2D) zVal = 0;
    else zVal = 999999999999;
    
    gcgSETVECTOR3(maxDerivatives, -9999999999, -9999999999, -1.*zVal);
    gcgSETVECTOR3(minDerivatives, 9999999999, 9999999999, zVal);

    readField();

    //creating glyphs from the tensor field
    bigger = -1;
    createTensorGlyphs(glyphs_original, &noriginals, field, confidence, fwidth, fheight, fdepth, &bigger, true);

    if(tipoCampo == FIELD2D){
        derivateNorms();
    }
    else{
//        derivateNorms3D();
//        derivateNorms3D2();
//        derivateNorms3D3();
        derivateNorms3D4();
    }

    ordSize = 0;

    printf("Max %f | Min %f\n",  maxNorm, minNorm);

    for (unsigned int z = 0; z < fdepth; z++)
        for (unsigned int y = 0; y < fheight; y++)
            for (unsigned int x = 0; x < fwidth; x++) {
                createTensorlineDup(velField, glyphs_original, &noriginals, width, height, depth, x, y, z, bigger, ADDRESSVOLUME(x, y, z, width, height));
            }

    for (int i = 0; i < totalSize; i++) {
       if (glyphs_original[i]->cs < 0.8) {
//           printf("akii");
           tensor_ord[ordSize] = glyphs_original[i];
           ordSize++;
       }
        
    }

    cout << "Total ordenados: " << ordSize - 1 << endl;

    if(ordSize < nParticles) nParticles = ordSize;

//    gcgCOPYVECTOR3(obs, frustum.view);
//    gcgCOPYVECTOR3(fobs, frustum.view);
//    calcDist2Obs(tensor_ord, ordSize);
//    calcConsts(tensor_ord, ordSize);
//    calcK1(tensor_ord, ordSize);
//    calcK2(tensor_ord, ordSize);
//    calcK3(tensor_ord, ordSize);
//    calcKeys(tensor_ord, ordSize, subglyphs);
//    mergesort(tensor_ord, 0, ordSize -1);

//    cout << "-----------------------------------------------------------------------------" << endl;
//    cout << tensor_ord[5]->cl << " " << tensor_ord[6]->cl << " " << tensor_ord[7]->cl << endl;
//    cout << tensor_ord[5]->cp << " " << tensor_ord[6]->cp << " " << tensor_ord[7]->cp << endl;
//    cout << tensor_ord[5]->cs << " " << tensor_ord[6]->cs << " " << tensor_ord[7]->cs << endl;


    float r;
    int p, posX, posY, posZ;
    VECTOR3 pos = {0.0, 0.0, 0.0};
    VECTOR3 mainDirect = {0.0, 0.0, 0.0};
    int n = 0;
    bool restrictedSearch = true;
//    int n;

    /************Criar partículas**************/
//    for (unsigned int i = 0; i < nParticles; i++) {
//
//        criacao
//
//    }

}

void gcgTensor::derivateNorms(){

    VECTOR3 deriv = {0,0,0};
    int v, i;
    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){
                i = ADDRESSGLYPHS(x, y, z, fheight, fheight);
                v = i;
                deriv[0] += -norms[i];
                deriv[1] += norms[i];

                if(((x+1) >= 0 ) || ((x+1) < fwidth)){ ///tirar comparacoes desnecessarias dps
                    i = ADDRESSGLYPHS(x+1, y, z, fheight, fheight);
                    deriv[0] += norms[i];
                    deriv[1] += norms[i];
                    if(((y-1) >= 0 ) || ((y-1) < fheight)){
                        i = ADDRESSGLYPHS(x+1, y-1, z, fheight, fheight);
                        deriv[0] += norms[i];
                        deriv[1] += -norms[i];
                    }
                }
                if(((y-1) >= 0 ) || ((y-1) < fheight)){
                    i = ADDRESSGLYPHS(x, y-1, z, fheight, fheight);
                    deriv[0] += -norms[i];
                    deriv[1] += -norms[i];
                }

                if(deriv[0] > maxDerivatives[0]) maxDerivatives[0] = deriv[0];
                if(deriv[0] < minDerivatives[0]) minDerivatives[0] = deriv[0];
                if(deriv[1] > maxDerivatives[1]) maxDerivatives[1] = deriv[1];
                if(deriv[1] < minDerivatives[1]) minDerivatives[1] = deriv[1];

                
                float tmpV = gcgLENGTHVECTOR3(deriv);
                if(tmpV > maxDerivN) maxDerivN = tmpV;
                if(tmpV < minDerivN) minDerivN = tmpV;
                
                gcgCOPYVECTOR3(normDerivatives[v], deriv);

                gcgSETVECTOR3(deriv, 0., 0., 0.);
            }
        }
    }

}

void gcgTensor::derivateNorms3D(){

    VECTOR3 deriv = {0,0,0};
    int v, i, pos;
    int xi, yj, zk;
    bool xmin = false, xmax = false,ymin = false, ymax = false, zmin = false, zmax = false;
    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){
                v = ADDRESSGLYPHS(x, y, z, fwidth, fheight);

                if(x-1 >= 0) xmin = true;
                if(x+1 < fwidth) xmax = true;
                if(y-1 >= 0) ymin = true;
                if(y+1 < fheight) ymax = true;
                if(z-1 >= 0) zmin = true;
                if(z+1 < fdepth) zmax = true;

                ///Gx e Gy
                if(xmin){
                    pos = ADDRESSGLYPHS(x-1, y, z, fwidth, fheight);
                    deriv[0] -= 2*norms[pos];
                }
                if(xmax){
                    pos = ADDRESSGLYPHS(x+1, y, z, fwidth, fheight);
                    deriv[0] += 2*norms[pos];
                }

                if(ymin){
                    pos = ADDRESSGLYPHS(x, y-1, z, fwidth, fheight);
                    deriv[1] += 2*norms[pos];
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y-1, z, fwidth, fheight);
                        deriv[0] -= norms[pos];
                        deriv[1] += norms[pos];
                    }
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y-1, z, fwidth, fheight);
                        deriv[0] += norms[pos];
                        deriv[1] += norms[pos];
                    }
                }
                if(ymax){
                    pos = ADDRESSGLYPHS(x, y+1, z, fwidth, fheight);
                    deriv[1] -= 2*norms[pos];
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y+1, z, fwidth, fheight);
                        deriv[0] -= norms[pos];
                        deriv[1] -= norms[pos];
                    }
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y+1, z, fwidth, fheight);
                        deriv[0] += norms[pos];
                        deriv[1] -= norms[pos];
                    }
                }

                ///Gz Primeira matriz (z = -1)
                if(zmin){
                    pos = ADDRESSGLYPHS(x, y, z-1, fwidth, fheight);
                    deriv[2] += 4*norms[pos];
                    ///Linha da esquerda
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                    }
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                    }

                    if(ymin){
                        pos = ADDRESSGLYPHS(x, y-1, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                    }
                    if(ymax){
                        pos = ADDRESSGLYPHS(x, y+1, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                    }
                }
                /// Z+1

                if(zmax){
                    pos = ADDRESSGLYPHS(x, y, z+1, fwidth, fheight);
                    deriv[2] -= 4*norms[pos];
                    ///Linha da esquerda
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y, z+1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                    }
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y, z+1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                    }

                    if(ymin){
                        pos = ADDRESSGLYPHS(x, y-1, z+1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z+1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z+1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                    }
                    if(ymax){
                        pos = ADDRESSGLYPHS(x, y+1, z+1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z+1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z+1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                    }
                }

                xmin = false;
                xmax = false;
                ymin = false;
                ymax = false;
                zmin = false;
                zmax = false;


                if(deriv[0] > maxDerivatives[0]) maxDerivatives[0] = deriv[0];
                if(deriv[0] < minDerivatives[0]) minDerivatives[0] = deriv[0];
                if(deriv[1] > maxDerivatives[1]) maxDerivatives[1] = deriv[1];
                if(deriv[1] < minDerivatives[1]) minDerivatives[1] = deriv[1];
                if(deriv[2] > maxDerivatives[2]) maxDerivatives[2] = deriv[2];
                if(deriv[2] < minDerivatives[2]) minDerivatives[2] = deriv[2];

                float tmpV = gcgLENGTHVECTOR3(deriv);
                if(tmpV > maxDerivN) maxDerivN = tmpV;
                if(tmpV < minDerivN) minDerivN = tmpV;

                gcgCOPYVECTOR3(normDerivatives[v], deriv);

                gcgSETVECTOR3(deriv, 0., 0., 0.);
            }
        }
    }

}

void gcgTensor::derivateNorms3D2(){

    VECTOR3 deriv = {0,0,0};
    int v, i, pos;
    int xi, yj, zk;
    bool xmin = false, xmax = false,ymin = false, ymax = false, zmin = false, zmax = false;
    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){
                v = ADDRESSGLYPHS(x, y, z, fwidth, fheight);

                if(x-1 >= 0) xmin = true;
                if(x+1 < fwidth) xmax = true;
                if(y-1 >= 0) ymin = true;
                if(y+1 < fheight) ymax = true;
                if(z-1 >= 0) zmin = true;
                if(z+1 < fdepth) zmax = true;

                ///Primeira matriz
                if(xmin){
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x-1, y-1, z, fwidth, fheight);
                        deriv[1] -= 3*norms[pos];
                        deriv[2] -= 3*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z-1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                            deriv[1] -= norms[pos];
                            deriv[2] -= norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x-1, y-1, z+1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                            deriv[1] += norms[pos];
                            deriv[2] -= norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x-1, y, z, fwidth, fheight);
                    deriv[0] -= 6*norms[pos];
                    if(zmin){
                        pos = ADDRESSGLYPHS(x-1, y, z-1, fwidth, fheight);
                        deriv[0] -= 3*norms[pos];
                        deriv[1] -= 3*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x-1, y, z+1, fwidth, fheight);
                        deriv[0] -= 3*norms[pos];
                        deriv[1] += 3*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x-1, y+1, z, fwidth, fheight);
                        deriv[0] -= 3*norms[pos];
                        deriv[2] += 3*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z-1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                            deriv[1] -= norms[pos];
                            deriv[2] += norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x-1, y+1, z+1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                            deriv[1] += norms[pos];
                            deriv[2] += norms[pos];
                        }
                    }
                }
                ///Segunda Matriz
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x, y-1, z, fwidth, fheight);
                        deriv[2] -= 6*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x, y-1, z-1, fwidth, fheight);
                            deriv[1] -= 3*norms[pos];
                            deriv[2] -= 3*norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x, y-1, z+1, fwidth, fheight);
                            deriv[1] += 3*norms[pos];
                            deriv[2] -= 3*norms[pos];
                        }
                    }
                    ///Linha do meio
                    if(zmin){
                        pos = ADDRESSGLYPHS(x, y, z-1, fwidth, fheight);
                        deriv[1] -= 6*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x, y, z+1, fwidth, fheight);
                        deriv[1] += 6*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x, y+1, z, fwidth, fheight);
                        deriv[2] += 6*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x, y+1, z-1, fwidth, fheight);
                            deriv[1] -= 3*norms[pos];
                            deriv[2] += 3*norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x, y+1, z+1, fwidth, fheight);
                            deriv[1] += 3*norms[pos];
                            deriv[2] += 3*norms[pos];
                        }
                    }
                ///Terceira matriz
                    if(xmax){
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x+1, y-1, z, fwidth, fheight);
                        deriv[0] += 3*norms[pos];
                        deriv[2] -= 3*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x+1, y-1, z-1, fwidth, fheight);
                            deriv[0] += norms[pos];
                            deriv[1] -= norms[pos];
                            deriv[2] -= norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z+1, fwidth, fheight);
                            deriv[0] += norms[pos];
                            deriv[1] += norms[pos];
                            deriv[2] -= norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x+1, y, z, fwidth, fheight);
                    deriv[0] += 6*norms[pos];
                    if(zmin){
                        pos = ADDRESSGLYPHS(x+1, y, z-1, fwidth, fheight);
                        deriv[0] += 3*norms[pos];
                        deriv[1] -= 3*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x+1, y, z+1, fwidth, fheight);
                        deriv[0] += 3*norms[pos];
                        deriv[1] += 3*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x+1, y+1, z, fwidth, fheight);
                        deriv[0] += 3*norms[pos];
                        deriv[2] += 3*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x+1, y+1, z-1, fwidth, fheight);
                            deriv[0] += norms[pos];
                            deriv[1] -= norms[pos];
                            deriv[2] += norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z+1, fwidth, fheight);
                            deriv[0] += norms[pos];
                            deriv[1] += norms[pos];
                            deriv[2] += norms[pos];
                        }
                    }
                }

                xmin = false;
                xmax = false;
                ymin = false;
                ymax = false;
                zmin = false;
                zmax = false;


                if(deriv[0] > maxDerivatives[0]) maxDerivatives[0] = deriv[0];
                if(deriv[0] < minDerivatives[0]) minDerivatives[0] = deriv[0];
                if(deriv[1] > maxDerivatives[1]) maxDerivatives[1] = deriv[1];
                if(deriv[1] < minDerivatives[1]) minDerivatives[1] = deriv[1];
                if(deriv[2] > maxDerivatives[2]) maxDerivatives[2] = deriv[2];
                if(deriv[2] < minDerivatives[2]) minDerivatives[2] = deriv[2];

                float tmpV = gcgLENGTHVECTOR3(deriv);
                if(tmpV > maxDerivN) maxDerivN = tmpV;
                if(tmpV < minDerivN) minDerivN = tmpV;

                gcgCOPYVECTOR3(normDerivatives[v], deriv);

                gcgSETVECTOR3(deriv, 0., 0., 0.);
            }
        }
    }

}

void gcgTensor::derivateNorms3D3(){

    VECTOR3 deriv = {0,0,0};
    int v, i, pos;
    int xi, yj, zk;
    bool xmin = false, xmax = false,ymin = false, ymax = false, zmin = false, zmax = false;
    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){
                v = ADDRESSGLYPHS(x, y, z, fwidth, fheight);

                if(x-1 >= 0) xmin = true;
                if(x+1 < fwidth) xmax = true;
                if(y-1 >= 0) ymin = true;
                if(y+1 < fheight) ymax = true;
                if(z-1 >= 0) zmin = true;
                if(z+1 < fdepth) zmax = true;

                ////////////////////////
                //////////Dx
                ////////////////////////
                if(xmin){
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x-1, y-1, z, fwidth, fheight);
                        deriv[0] += 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z-1, fwidth, fheight);
                            deriv[0] += norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x-1, y-1, z+1, fwidth, fheight);
                            deriv[0] += norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x-1, y, z, fwidth, fheight);
                    deriv[0] += 4*norms[pos];
                    if(zmin){
                        pos = ADDRESSGLYPHS(x-1, y, z-1, fwidth, fheight);
                        deriv[0] += 2*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x-1, y, z+1, fwidth, fheight);
                        deriv[0] += 2*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x-1, y+1, z, fwidth, fheight);
                        deriv[0] += 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z-1, fwidth, fheight);
                            deriv[0] += norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x-1, y+1, z+1, fwidth, fheight);
                            deriv[0] += norms[pos];
                        }
                    }
                }
                if(xmax){
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x+1, y-1, z, fwidth, fheight);
                        deriv[0] -= 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x+1, y-1, z-1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z+1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x+1, y, z, fwidth, fheight);
                    deriv[0] -= 4*norms[pos];
                    if(zmin){
                        pos = ADDRESSGLYPHS(x+1, y, z-1, fwidth, fheight);
                        deriv[0] -= 2*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x+1, y, z+1, fwidth, fheight);
                        deriv[0] -= 2*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x+1, y+1, z, fwidth, fheight);
                        deriv[0] -= 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x+1, y+1, z-1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z+1, fwidth, fheight);
                            deriv[0] -= norms[pos];
                        }
                    }
                }

                /////////////////////
                //Dy
                /////////////////////

                if(ymin){
                    ///Linha da esquerda
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y-1, z, fwidth, fheight);
                        deriv[1] += 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z-1, fwidth, fheight);
                            deriv[1] += norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x-1, y-1, z+1, fwidth, fheight);
                            deriv[1] += norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x, y-1, z, fwidth, fheight);
                    deriv[1] += 4*norms[pos];
                    if(zmin){
                        pos = ADDRESSGLYPHS(x, y-1, z-1, fwidth, fheight);
                        deriv[1] += 2*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x, y-1, z+1, fwidth, fheight);
                        deriv[1] += 2*norms[pos];
                    }
                    ///Linha da direita
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y-1, z, fwidth, fheight);
                        deriv[1] += 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x+1, y-1, z-1, fwidth, fheight);
                            deriv[1] += norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z+1, fwidth, fheight);
                            deriv[1] += norms[pos];
                        }
                    }
                }
                if(ymax){
                    ///Linha da esquerda
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y+1, z, fwidth, fheight);
                        deriv[1] -= 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z-1, fwidth, fheight);
                            deriv[1] -= norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x-1, y+1, z+1, fwidth, fheight);
                            deriv[1] -= norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x, y+1, z, fwidth, fheight);
                    deriv[1] -= 4*norms[pos];
                    if(zmin){
                        pos = ADDRESSGLYPHS(x, y+1, z-1, fwidth, fheight);
                        deriv[1] -= 2*norms[pos];
                    }
                    if(zmax){
                        pos = ADDRESSGLYPHS(x, y+1, z+1, fwidth, fheight);
                        deriv[1] -= 2*norms[pos];
                    }
                    ///Linha da direita
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y+1, z, fwidth, fheight);
                        deriv[1] -= 2*norms[pos];
                        if(zmin){
                            pos = ADDRESSGLYPHS(x+1, y+1, z-1, fwidth, fheight);
                            deriv[1] -= norms[pos];
                        }
                        if(zmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z+1, fwidth, fheight);
                            deriv[1] -= norms[pos];
                        }
                    }
                }

                /////////////////////
                //Dz
                /////////////////////

                if(zmin){
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x, y-1, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x, y, z-1, fwidth, fheight);
                    deriv[2] += 4*norms[pos];
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                    }
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x, y+1, z-1, fwidth, fheight);
                        deriv[2] += 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z-1, fwidth, fheight);
                            deriv[2] += norms[pos];
                        }
                    }
                }
                if(zmax){
                    ///Linha da esquerda
                    if(ymin){
                        pos = ADDRESSGLYPHS(x, y-1, z-1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y-1, z-1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y-1, z-1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                    }
                    ///Linha do meio
                    pos = ADDRESSGLYPHS(x, y, z-1, fwidth, fheight);
                    deriv[2] -= 4*norms[pos];
                    if(xmin){
                        pos = ADDRESSGLYPHS(x-1, y, z-1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                    }
                    if(xmax){
                        pos = ADDRESSGLYPHS(x+1, y, z-1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                    }
                    ///Linha da direita
                    if(ymax){
                        pos = ADDRESSGLYPHS(x, y+1, z-1, fwidth, fheight);
                        deriv[2] -= 2*norms[pos];
                        if(xmin){
                            pos = ADDRESSGLYPHS(x-1, y+1, z-1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                        if(xmax){
                            pos = ADDRESSGLYPHS(x+1, y+1, z-1, fwidth, fheight);
                            deriv[2] -= norms[pos];
                        }
                    }
                }


                //////////////////////
                xmin = false;
                xmax = false;
                ymin = false;
                ymax = false;
                zmin = false;
                zmax = false;


                if(deriv[0] > maxDerivatives[0]) maxDerivatives[0] = deriv[0];
                if(deriv[0] < minDerivatives[0]) minDerivatives[0] = deriv[0];
                if(deriv[1] > maxDerivatives[1]) maxDerivatives[1] = deriv[1];
                if(deriv[1] < minDerivatives[1]) minDerivatives[1] = deriv[1];
                if(deriv[2] > maxDerivatives[2]) maxDerivatives[2] = deriv[2];
                if(deriv[2] < minDerivatives[2]) minDerivatives[2] = deriv[2];

                float tmpV = gcgLENGTHVECTOR3(deriv);
                if(tmpV > maxDerivN) maxDerivN = tmpV;
                if(tmpV < minDerivN) minDerivN = tmpV;

                gcgCOPYVECTOR3(normDerivatives[v], deriv);

                gcgSETVECTOR3(deriv, 0., 0., 0.);
            }
        }
    }

}

void gcgTensor::applyGaussian(){

    printf("---------------------------\nApplying Gaussian!!!\n---------------------\n\n");
    
    MATRIX3 * newField = (MATRIX3*)malloc(sizeof(MATRIX3)*fdepth*fheight*fwidth);
    MATRIX3 tx, ty, tz;
    int v, i, pos;
    int xi, yj, zk;
    int countx;
    int county;
    int countz;
    float maxN = -9999999;
    bool xmin = false, xmax = false,ymin = false, ymax = false, zmin = false, zmax = false;
    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){

                v = ADDRESSGLYPHS(x, y, z, fwidth, fheight);

                gcgSCALEMATRIX3(tx, glyphs_original[v]->tensor, 2.0);
                gcgCOPYMATRIX3(ty, tx);
                gcgCOPYMATRIX3(tz, tx);

                countx = 2;
                county = 2;
                countz = 2;

                if(x-1 >= 0) xmin = true;
                if(x+1 < fwidth) xmax = true;
                if(y-1 >= 0) ymin = true;
                if(y+1 < fheight) ymax = true;
                if(z-1 >= 0) zmin = true;
                if(z+1 < fdepth) zmax = true;

                ////////////////////////
                //////////Dx
                ////////////////////////
                if(xmin){
                    pos = ADDRESSGLYPHS(x-1, y, z, fwidth, fheight);
                    gcgADDMATRIX3(tx, tx, glyphs_original[pos]->tensor);
                    countx++;
                }
                if(xmax){
                    pos = ADDRESSGLYPHS(x+1, y, z, fwidth, fheight);
                    gcgADDMATRIX3(tx, tx, glyphs_original[pos]->tensor);
                    countx++;
                }
                if(ymin){
                    pos = ADDRESSGLYPHS(x, y-1, z, fwidth, fheight);
                    gcgADDMATRIX3(ty, ty, glyphs_original[pos]->tensor);
                    county++;
                }
                if(ymax){
                    pos = ADDRESSGLYPHS(x, y+1, z, fwidth, fheight);
                    gcgADDMATRIX3(ty, ty, glyphs_original[pos]->tensor);
                    county++;
                }
                if(zmin){
                    pos = ADDRESSGLYPHS(x, y, z-1, fwidth, fheight);
                    gcgADDMATRIX3(tz, tz, glyphs_original[pos]->tensor);
                    countz++;
                }
                if(zmax){
                    pos = ADDRESSGLYPHS(x, y, z+1, fwidth, fheight);
                    gcgADDMATRIX3(tz, tz, glyphs_original[pos]->tensor);
                    countz++;
                }

                //////////////////////
                xmin = false;
                xmax = false;
                ymin = false;
                ymax = false;
                zmin = false;
                zmax = false;

                gcgSCALEMATRIX3(tx, tx, 1./countx);
                gcgSCALEMATRIX3(ty, ty, 1./county);
                gcgSCALEMATRIX3(tz, tz, 1./countz);

                gcgADDMATRIX3(tx, tx, ty);
                gcgADDMATRIX3(tx, tx, tz);
                gcgSCALEMATRIX3(tx, tx, 1./3.);
                float n = calculateNorm(tx);
                if(n > maxN) maxN = n;

                gcgCOPYMATRIX3(newField[v], tx);
            }
        }
    }

    ///Normalizing

    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){

                v = ADDRESSGLYPHS(x, y, z, fwidth, fheight);

                MATRIX3 nt;

                gcgSCALEMATRIX3(nt, newField[v], 1./maxN);
                float n = calculateNorm(nt);
                if(n > maxNorm) maxNorm = n;
                if(n < minNorm) minNorm = n;
                glyphs_original[v]->setTensor(nt);
            }
        }
    }

    storeField("smoothhelix", newField);
    exit(0);

    free(newField);

}

void gcgTensor::derivateNorms3D4(){

    VECTOR3 deriv = {0,0,0};
    int v, i, pos;
    int xi, yj, zk;
    bool xmin = false, xmax = false,ymin = false, ymax = false, zmin = false, zmax = false;
    for(int x = 0; x < fwidth; x++){
        for(int y = 0; y < fheight; y++){
            for(int z = 0; z < fdepth; z++){
                v = ADDRESSGLYPHS(x, y, z, fwidth, fheight);

                if(x-1 >= 0) xmin = true;
                if(x+1 < fwidth) xmax = true;
                if(y-1 >= 0) ymin = true;
                if(y+1 < fheight) ymax = true;
                if(z-1 >= 0) zmin = true;
                if(z+1 < fdepth) zmax = true;

                ////////////////////////
                //////////Dx
                ////////////////////////
                if(xmin){
                    pos = ADDRESSGLYPHS(x-1, y, z, fwidth, fheight);
                    deriv[0] -= norms[pos];
                }
                if(xmax){
                    pos = ADDRESSGLYPHS(x+1, y, z, fwidth, fheight);
                    deriv[0] += norms[pos];
                }
                if(ymin){
                    pos = ADDRESSGLYPHS(x, y-1, z, fwidth, fheight);
                    deriv[1] -= norms[pos];
                }
                if(ymax){
                    pos = ADDRESSGLYPHS(x, y+1, z, fwidth, fheight);
                    deriv[1] += norms[pos];
                }
                if(zmin){
                    pos = ADDRESSGLYPHS(x, y, z-1, fwidth, fheight);
                    deriv[2] -= norms[pos];
                }
                if(zmax){
                    pos = ADDRESSGLYPHS(x, y, z+1, fwidth, fheight);
                    deriv[2] += norms[pos];
                }

                //////////////////////
                xmin = false;
                xmax = false;
                ymin = false;
                ymax = false;
                zmin = false;
                zmax = false;


                if(deriv[0] > maxDerivatives[0]) maxDerivatives[0] = deriv[0];
                if(deriv[0] < minDerivatives[0]) minDerivatives[0] = deriv[0];
                if(deriv[1] > maxDerivatives[1]) maxDerivatives[1] = deriv[1];
                if(deriv[1] < minDerivatives[1]) minDerivatives[1] = deriv[1];
                if(deriv[2] > maxDerivatives[2]) maxDerivatives[2] = deriv[2];
                if(deriv[2] < minDerivatives[2]) minDerivatives[2] = deriv[2];


                float tmpV = gcgLENGTHVECTOR3(deriv);
                if(tmpV > maxDerivN) maxDerivN = tmpV;
                if(tmpV < minDerivN) minDerivN = tmpV;

                gcgCOPYVECTOR3(normDerivatives[v], deriv);

                gcgSETVECTOR3(deriv, 0., 0., 0.);
            }
        }
    }

}


gcgTensor::~gcgTensor() {
}

void gcgTensor::readField() {
//    subfields = (MATRIX3 **) malloc (numSubFields*sizeof(MATRIX3*));
//    subFieldDim = (VECTOR3*) malloc (numSubFields*sizeof(VECTOR3));

    if(tipoCampo == CEREBRO){
//        readRawField(totalSize, "smaller.raw", field, confidence);
        readTXTField3(0, "smoothhelix", field, totalSize);
//        readTXTField3(0, "smoothbrainori", field, totalSize);
//        subfields[0] = readTXTField3(0, "smallermenor");
//        subfields[1] = readTXTField3(1, "smallermenor2");

    }
    else if(tipoCampo == HELICE){
        readRawField(totalSize, "dt-helix.raw", field, confidence);
//        subfields[0] = readTXTField3(0, "helicemenor");
//        subfields[1] = readTXTField3(1, "helicemenor2");
    }
    else if(tipoCampo == CEREBRO2){
//        readRawField(totalSize, "dt-helix.raw", field, confidence);
        readTXTField3(0, "smallermenor2", field, totalSize);
//        subfields[1] = readTXTField3(1, "helicemenor2");
    }
    else if(tipoCampo == CEREBRO1){
//        readRawField(totalSize, "dt-helix.raw", field, confidence);
        readTXTField3(0, "smallermenor", field, totalSize);
//        readTXTField3(0, "smoothbrain", field, totalSize);
//        readTXTField3(0, "smoothbrainorinn", field, totalSize);
//        subfields[1] = readTXTField3(1, "helicemenor2");
    }
    else if(tipoCampo == HELICE1){
//        readRawField(totalSize, "dt-helix.raw", field, confidence);
//        readTXTField3(0, "smallermenor", field, totalSize);
        readTXTField3(0, "helicemenor", field, totalSize);
//        subfields[1] = readTXTField3(1, "helicemenor2");
    }
    else if(tipoCampo == PONTAS){
        readTXTField(totalSize, "field1.txt", field);
//        subfields[0] = readTXTField3(0, "field1menor");
//        subfields[1] = readTXTField3(1, "field1menor2");
    }
    else if(tipoCampo == FIELD2D){
        readNewField(totalSize, "field64center.bin", field);
//        readTXTField3(0, "f2dsmooth", field, totalSize);
//        subfields[0] = readTXTField3(0, "field1menor");
//        subfields[1] = readTXTField3(1, "field1menor2");
    }
    else{
        printf("Tipo de campo invalido.");
        exit(0);
    }
}

void gcgTensor::createTensorGlyphs(gcgTENSORGLYPH** glyphs, int *nglyphs, MATRIX3* tensorField, float* confidence, unsigned int width, unsigned int height, unsigned int depth, float *bigger, bool isOriginal) {

    //creating the tensor glyphs and finding the bigger the smaller weights
//    FILE * f = fopen("dbg.txt", "w");
    float smaller = HUGE_VAL;
    for (unsigned int i = 0; i < width; i++)
        for (unsigned int j = 0; j < height; j++)
            for (unsigned int k = 0; k < depth; k++) {
                VECTOR3 v = {(float) i,(float) j,(float) k};
//                VECTOR3 v = {(float) i, (float) j, (float) k};

                //properties of the tensor matrix
                MATRIX3 eigenVectors;
                VECTOR3 eigenValues;
//                if(dbg){
//                    printf("DBG: % d % d %d\n", width, height, depth);
//                    float aaa = (j*width) + i + (k*width*height);
//                    printf("i: %d | j: %d | k: %d | addv2: %d | addv2n %f\n", i, j, k, ADDRESSVOLUME2((i), (j), (k), width, height), aaa);
//                }
                gcgEigenSymmetricMatrix3(tensorField[ADDRESSVOLUME((i), (j), (k), width, height)], eigenVectors, eigenValues);




                int addressG, addressV;
                addressG = ADDRESSGLYPHS(i, j, k, width, height);
                addressV = ADDRESSVOLUME(i, j, k, width, height);

                /******************       Trocado sem wavelet         *********************/
                glyphs[addressG] = new gcgTENSORGLYPH(tensorField[addressV], v, confidence[addressV], isOriginal, 0, NULL, i, j, k, fwidth, fheight, fdepth);

                /******Zerando******/
//                if(glyphs[addressG]->cs > 0.7){
//                    gcgSCALEMATRIX3(glyphs[addressG]->tensor, tensorField[addressV], 0.);
//                    gcgCOPYMATRIX3(glyphs[addressG]->eigenVectors, glyphs[addressG]->tensor);
//                    gcgZEROVECTOR3(glyphs[addressG]->eigenValues);
//                    glyphs[addressG]->cs = 1.;
//                    glyphs[addressG]->cl = 0.;
//                    glyphs[addressG]->cp = 0.;
//                }

                /*******************/

                norms[addressG] = calculateNorm(glyphs[addressG]->tensor);
                if(norms[addressG] < minNorm) minNorm = norms[addressG];
                if(norms[addressG] > maxNorm) maxNorm = norms[addressG];
//                glyphs[addressG] = new gcgTENSORGLYPH(tensorField[addressV], v, confidence[addressV], isOriginal, numSubFields, subglyphs, i, j, k, fwidth, fheight, fdepth);
//                if (glyphs[addressG]->weight > (*bigger)) (*bigger) = glyphs[addressG]->weight;
//                if (glyphs[addressG]->weight < smaller) smaller = glyphs[addressG]->weight;
//
//                (*nglyphs)++;
                
//                if(norms[addressG] >= 1. || norms[addressG] < 0) printf("Norma: %f\n", norms[addressG]);

                
                
            }

    float totalMaxNorm = maxNorm;
    printf("Max: %f | Min: %f\n", maxNorm, minNorm);
//    FILE * f = fopen("dbg.txt", "w");
    if(tipoCampo == CEREBRO || tipoCampo == CEREBRO1){
    for (unsigned int i = 0; i < width; i++)
        for (unsigned int j = 0; j < height; j++)
            for (unsigned int k = 0; k < depth; k++) {
                VECTOR3 v = {(float) i,(float) j,(float) k};
                int addressG = ADDRESSGLYPHS(i, j, k, width, height);
                float bf, af;
                MATRIX3 nt;
                gcgCOPYMATRIX3(nt, glyphs[addressG]->tensor);
//                bf = calculateNorm(nt);
//                fprintf(f, "%f %f %f %f %f %f %f %f %f --- %f\n", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5], nt[6], nt[7], nt[8], bf);
                float fdiv = 1./maxNorm;
                gcgSCALEMATRIX3(nt, nt, fdiv);
//                af = calculateNorm(nt);
//                fprintf(f, "%f %f %f %f %f %f %f %f %f --- %f\n", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5], nt[6], nt[7], nt[8], af);
//                glyphs[addressG]->setTensor(nt);
                
                ///Noa tava normalizando antes
                norms[addressG] = calculateNorm(glyphs[addressG]->tensor);
                if(norms[addressG] < minNorm) minNorm = norms[addressG];
                if(norms[addressG] > maxNorm) maxNorm = norms[addressG];
//                gcgCOPYMATRIX3(nt, glyphs[addressG]->tensor);
//                fprintf(f, "%f %f %f %f %f %f %f %f %f --- %f\n", nt[0], nt[1], nt[2], nt[3], nt[4], nt[5], nt[6], nt[7], nt[8], calculateNorm(nt));
//                if(af < 0 || af > 1) printf("Antes: %f | Dps: %f\n", bf, af);
                
            }
    }
//    fclose(f);
//    exit(0);
    
//    minNorm = 99999999999;
//    maxNorm = -99999999999;
    
//    MATRIX3 aa = {0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.0};
//    int addressG = ADDRESSGLYPHS(30, 30, 0, width, height);
//    glyphs[addressG]->setTensor(aa);
//    addressG = ADDRESSGLYPHS(32, 32, 0, width, height);
//    glyphs[addressG]->setTensor(aa);
//    addressG = ADDRESSGLYPHS(30, 32, 0, width, height);
//    glyphs[addressG]->setTensor(aa);
//    addressG = ADDRESSGLYPHS(32, 47, 0, width, height);
//    printf("%f %f %f %f %f %f %f %f %f\n", glyphs[addressG]->tensor[0], glyphs[addressG]->tensor[1], glyphs[addressG]->tensor[2], glyphs[addressG]->tensor[3], glyphs[addressG]->tensor[4], glyphs[addressG]->tensor[5], glyphs[addressG]->tensor[6], glyphs[addressG]->tensor[7], glyphs[addressG]->tensor[8]);
    
//    applyGaussian();

    printf("weights: raiz(lambda1^2 + lambda2^2 + lambda3^2)\nmaior: %f menor: %f\n", (*bigger), smaller);


//    for (unsigned int i = 0; i < width; i++)
//        for (unsigned int j = 0; j < height; j++)
//            for (unsigned int k = 0; k < depth; k++) {
//                VECTOR3 v = {(float) i,(float) j,(float) k};
////                VECTOR3 v = {(float) i, (float) j, (float) k};
//                int addressG, addressV;
//                addressG = ADDRESSGLYPHS(i, j, k, width, height);
//                addressV = ADDRESSVOLUME(i, j, k, width, height);
//
//                if(FEQUAL(glyphs[addressG]->cl, 1.)){continue;}
//                //properties of the tensor matrix
//                MATRIX3 eigenVectors, newTensor;
//                VECTOR3 eigenValues;
//                gcgSCALEMATRIX3(newTensor, glyphs[addressG]->tensor, 1./totalMaxNorm);
//
//                glyphs[addressG]->setTensor(newTensor);
//
//                /******Zerando******/
////                if(glyphs[addressG]->cs > 0.7){
////                    gcgSCALEMATRIX3(glyphs[addressG]->tensor, tensorField[addressV], 0.);
////                    gcgCOPYMATRIX3(glyphs[addressG]->eigenVectors, glyphs[addressG]->tensor);
////                    gcgZEROVECTOR3(glyphs[addressG]->eigenValues);
////                    glyphs[addressG]->cs = 1.;
////                    glyphs[addressG]->cl = 0.;
////                    glyphs[addressG]->cp = 0.;
////                }
//
//
////                fprintf(f, "%f %f %f %f %f %f %f %f %f\n", glyphs[addressG]->tensor[0], glyphs[addressG]->tensor[1], glyphs[addressG]->tensor[2], glyphs[addressG]->tensor[3], glyphs[addressG]->tensor[4], glyphs[addressG]->tensor[5], glyphs[addressG]->tensor[6], glyphs[addressG]->tensor[7], glyphs[addressG]->tensor[8]);
//                /*******************/
//
//                norms[addressG] = calculateNorm(glyphs[addressG]->tensor);
//                if(norms[addressG] < minNorm) minNorm = norms[addressG];
//                if(norms[addressG] > maxNorm) maxNorm = norms[addressG];
//                if (glyphs[addressG]->weight > (*bigger)) (*bigger) = glyphs[addressG]->weight;
//                if (glyphs[addressG]->weight < smaller) smaller = glyphs[addressG]->weight;
//
//                (*nglyphs)++;
//
//
//            }


    int addressH;

    //adjusting weights
    float diff = (*bigger) - smaller;
    for (unsigned int i = 0; i < width; i++)
        for (unsigned int j = 0; j < height; j++)
            for (unsigned int k = 0; k < depth; k++) {
                addressH = ADDRESSGLYPHS(i, j, k, width, height);
                if (glyphs[addressH] != NULL) {
                    glyphs[addressH]->weight = 1 - (glyphs[addressH]->weight - smaller) / diff;
                }
            }
    printf("bigger: %f - amount of tensor: %d\n", *bigger, (*nglyphs));
//    fclose(f);


}


void gcgTensor::drawField(){

    glPushMatrix();

    glTranslatef(-fwidth / 2, -fheight / 2, -fdepth / 2);
//    printf("fwidth: %d %d %d\n", -fwidth / 2, -fheight / 2, -fdepth / 2);

    VECTOR3 posAVector;
    glDisable(GL_LIGHTING);

    for (int x = 0; x < fwidth; x++)
        for (int y = 0; y < fheight; y++)
            for (int z = 0; z < fdepth; z++) {

                int i = ADDRESSVOLUME(x, y, z, fwidth, fheight);

//                if(glyphs_original[i]->cs < 0.7) {continue;}
                if((tipoCampo != FIELD2D) && (glyphs_original[i]->cs > 0.7)) { continue; }
//                if((tipoCampo != FIELD2D) && ((glyphs_original[i]->cs+glyphs_original[i]->cl+glyphs_original[i]->cp) < 1.)) { continue; }
//                if((tipoCampo != FIELD2D) && (glyphs_original[i]->cs +> 0.8)) { continue; }
//                if(false){}
//                if (glyphs_original[i]->cs > 1.0) {printf("cs errado! %f\n", glyphs_original[i]->cs);}
//                if (glyphs_original[i]->cs <= 1.0 && glyphs_original[i]->cl <= 1.0 && glyphs_original[i]->cp <= 1.0) {}
                else {
//                    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //                        if (original) gcgHeatColor(glyphs_original[i]->cl, color);
    //                        else gcgHeatColor(glyphs_original[i]->equ, color);

                    glColor4f(glyphs_original[i]->FA * fabs(glyphs_original[i]->eigenVectors[0]),glyphs_original[i]->FA * fabs(glyphs_original[i]->eigenVectors[1]),glyphs_original[i]->FA * fabs(glyphs_original[i]->eigenVectors[0]), glyphs_original[i]->FA - 0.2);
                    VECTOR3 color;
                    gcgHeatColor(glyphs_original[i]->cl, color);
//                    glColor4f(glyphs_original[i]->cl * fabs(glyphs_original[i]->eigenVectors[0]),glyphs_original[i]->cl * fabs(glyphs_original[i]->eigenVectors[1]),glyphs_original[i]->cl * fabs(glyphs_original[i]->eigenVectors[0]), 1.);
//                    glColor4f(color[0], color[1], color[2], 1.);

//                    if((x == 30 && y == 30) || (x == 32 && y == 32) || (x == 30 && y == 32)) glColor4f(0.5, 0.5, 0.5, 0.5);
//                    else glColor4f(0., 0.0, 0.0, 1.0);
                    
//                    gcgSETVECTOR3(posAVector, velField[i][0], velField[i][1], velField[i][2]);
                    gcgSETVECTOR3(posAVector, glyphs_original[i]->eigenVectors[0], glyphs_original[i]->eigenVectors[1], glyphs_original[i]->eigenVectors[2]);
                    //                        gcgSETVECTOR3(posAVector, glyphs_original[i]->eigenVectors[0], glyphs_original[i]->eigenVectors[1], glyphs_original[i]->eigenVectors[2]);


                    gcgDrawVectorPyramid(glyphs_original[i]->pos[0], glyphs_original[i]->pos[1], glyphs_original[i]->pos[2], posAVector, 1.0f);
//                    gcgDrawVectorThetrahedron(glyphs_original[i]->pos[0], glyphs_original[i]->pos[1], glyphs_original[i]->pos[2], posAVector, 1.0f);


                    //                        gcgSETVECTOR3(posAVector, velField[i][0], velField[i][1], velField[i][2]);
                    //                        gcgSETVECTOR3(posAVector, glyphs_original[i]->eigenVectors[3], glyphs_original[i]->eigenVectors[4], glyphs_original[i]->eigenVectors[5]);
                    //                        gcgDrawVectorPyramid(glyphs_original[i]->pos[0], glyphs_original[i]->pos[1], glyphs_original[i]->pos[2], posAVector, 1.0f);
                    //                        glColor4f(0.0, 0.0, 1.0, 1.0);
                    //                        gcgSETVECTOR3(posAVector, velField[i][0], velField[i][1], velField[i][2]);
                    //                        gcgSETVECTOR3(posAVector, glyphs_original[i]->eigenVectors[6], glyphs_original[i]->eigenVectors[7], glyphs_original[i]->eigenVectors[8]);
                    //                        gcgDrawVectorPyramid(glyphs_original[i]->pos[0], glyphs_original[i]->pos[1], glyphs_original[i]->pos[2], posAVector, 1.0f);
                }
            }


    glPopMatrix();
}

void gcgTensor::changeDrawing(){
    normShow++;
    normShow = (int) (normShow%6);
}

void interpolateTensorCopy(VECTOR3 position, gcgTensor * teste, VECTOR3 * gl){

    int width = teste->fwidth;
    int height = teste->fheight;
    
    float x1, x2, y1, y2, w, h, dx1, dx2, dy1, dy2;
//    float norm11, norm12, norm21, norm22;
    VECTOR3 norm11, norm12, norm21, norm22;

    x1 = floor(position[0]);

    if(x1 >= width) x1 -= 1;

    x2 = x1+1.;

    if(x2 >= width) x2 -= 1;

    y1 = floor(position[1]);

    if(y1 >= height) y1 -= 1;

    y2 = y1+1.;

    if(y2 >= height) y2 -= 1;

    w = x2 - x1;
    h = y2 - y1;

    dx1 = (float)(x2 - position[0]) / w;
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
    dy1 = (y2 - position[1])/h;
    dy2 = (position[1] - y1)/h;

    gcgSCALEMATRIX3(temp1, R1, dy1);
    gcgSCALEMATRIX3(temp2, R2, dy2);

    gcgADDMATRIX3(P, temp1, temp2);

    VECTOR3 vt;
    vt[2] = 0;
    vt[0] = (dy1*((norm11[0]*dx1)+(norm21[0]*dx2)))+(dy2*((norm12[0]*dx1)+(norm22[0]*dx2)));
    vt[1] = (dy1*((norm11[1]*dx1)+(norm21[1]*dx2)))+(dy2*((norm12[1]*dx1)+(norm22[1]*dx2)));
    
    gcgCOPYVECTOR3(*gl, vt);
    
    
    

}


void gcgTensor::drawFewNorms(){
    
    glPushMatrix();

    glTranslatef(-fwidth/2, -fheight/2, -fdepth/2);
    
    VECTOR3 posAVector;
    MATRIX3 itpTensor;
    VECTOR3 g;
    VECTOR3 ps;
    
    glDisable(GL_LIGHTING);
    glColor4f(0, 0, 0.4, 1.);
    for (int x = 0.5; x < fwidth; x+= 2)
        for (int y = 0.5; y < fheight; y+=2)
            for (int z = 0.5; z < fdepth; z+=2) {
                
                int i = ADDRESSVOLUME(x, y, z, fwidth, fheight);
                gcgSETVECTOR3(ps, x, y, z);
                interpolateTensorCopy(ps, this, &g);
                gcgDrawVectorPyramid(x, y, z, g, 2.);
//                gcgDrawVectorPyramid(glyphs_original[i]->pos[0] - (fwidth/2), glyphs_original[i]->pos[1] - (fheight/2), glyphs_original[i]->pos[2] - (fdepth/2), normDerivatives[i], 2.);
                
            }
    
    glPopMatrix();
}

void gcgTensor::drawNorms(){

    if(normShow == 5){return;}
//    if(normShow == 4){drawFewNorms();return;}

    glPushMatrix();

    glTranslatef(-fwidth/2, -fheight/2, -fdepth/2);
//    printf("fwidth: %d %d %d\n", -fwidth / 2, -fheight / 2, -fdepth / 2);

    VECTOR3 posAVector;
    glDisable(GL_LIGHTING);
    glPointSize(12);
    glBegin(GL_POINTS);
    for (int x = 0; x < fwidth; x++)
        for (int y = 0; y < fheight; y++)
            for (int z = 0; z < fdepth; z++) {

                int i = ADDRESSVOLUME(x, y, z, fwidth, fheight);
                VECTOR3 color;
                if(normShow == 4){
                    float tmpV = gcgLENGTHVECTOR3(normDerivatives[i]);
                    gcgHeatColor((tmpV - minDerivN)/(maxDerivN - minDerivN), color);
//                    glColor3f(0., 0., 1.);
//                    glColor4f(color[0], color[1], color[2], 1.);
                    glColor4f(0, 0, 0.4, 1.);
                    gcgDrawVectorPyramid(glyphs_original[i]->pos[0] - (fwidth/2), glyphs_original[i]->pos[1] - (fheight/2), glyphs_original[i]->pos[2] - (fdepth/2), normDerivatives[i], 1.);
//                    gcgDrawVectorPyramid(glyphs_original[i]->pos[0], glyphs_original[i]->pos[1], glyphs_original[i]->pos[2], posAVector, 1.0f);
                }
                else{

                    switch(normShow){
                        case 0:
                            gcgHeatColor((norms[i] - minNorm)/(maxNorm - minNorm), color);
                            break;
                        case 1:
                            gcgHeatColor((normDerivatives[i][0] - minDerivatives[0])/(maxDerivatives[0] - minDerivatives[0]), color);
                            break;
                        case 2:
                            gcgHeatColor((normDerivatives[i][1] - minDerivatives[1])/(maxDerivatives[1] - minDerivatives[1]), color);
                            break;
                        case 3:
                            gcgHeatColor((normDerivatives[i][2] - minDerivatives[2])/(maxDerivatives[2] - minDerivatives[2]), color);
                            break;
                        default: break;
                    }

                    if(FEQUAL(gcgLENGTHVECTOR3(normDerivatives[i]), 0.0)){continue;}
                    glColor4f(color[0], color[1], color[2], 1.);
                    glVertex3f(glyphs_original[i]->pos[0], glyphs_original[i]->pos[1], glyphs_original[i]->pos[2]);
                }
            }

    glEnd();
    glPointSize(1);
    glPopMatrix();
}

//void gcgTensor::allocateAll(){
//
//
//
//}
