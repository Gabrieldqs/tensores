/*************************************************************************************
    GCG - GROUP FOR COMPUTER GRAPHICS
        Universidade Federal de Juiz de Fora
        Instituto de Ciências Exatas
        Departamento de Ciência da Computação

 **************************************************************************************/

#define NFEQUAL(a, b) (fabs((a)-(b)) < 0.0007)

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "gcg.h"
#include "tensorgenerate.h"
#include "octree.h"
#include "sph.h"
#include "gcgTensor.h"



gcgSph * sph;
extern gcgTensor * teste;
extern int translateW;
extern int translateH;
extern int translateD;
extern int tipoCampo;
extern bool filterParticles;
extern float filterT;
//////////////////////
/// Adicoes Dinamicas
//////////////////////
bool simulationStarted = false;
Octree * _octree;

gcgFRUSTUM frustum;

static gcgRANDOMGAUSSIAN gcgRandom;
void simulacao(void);


int windowW = 640, windowH = 480;

gcgTEXT texto;

//parametros das superquadricas
int slicesphi   = 7;
int slicestheta = 7;


//constantes para a GLUT
int janela;
int largura = 640;
int altura = 480;
int iniX, iniY;
int gira = FALSE;
int translada = FALSE;
float velocidadeFrente = 0;
float velocidadeLado = 0;
float dt = 0;
float bigger = -1;
clock_t tempoant = 0;
bool enableCreation = false;
//int ordSize = 0;

//opcoes da aplicacao
bool drawMesh = false;
bool drawPart = false;
bool drawQuad = false;
bool translationMode = false;
bool original = true;
bool cutVis = false;
bool interpolated_visible = false;
bool drawTensor = false;
bool drawNorms = true;
bool clearGl = true;
bool simula = false;
bool direc = true;
bool oldSimula = false;
VECTOR3 cutBounds1 = {0,0,0}, cutBounds2;
int currentCutPlane = 0;
extern unsigned int simulationStep;


float farDist, farVal = 250;

float alpha = 2.0; //sigma/alpha

//visualization equations paramter
GLfloat gAngle = 0.0;
GLUquadricObj *IDquadric;

/**************************/
/// Coordinate translation
////////////////////////////

void GetOGLPos(int x, int y, VECTOR3 pos)
{
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLfloat winX, winY, winZ;
    GLdouble posX, posY, posZ;

    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );

    winX = (float)x;
    winY = (float)viewport[3] - (float)y;
    glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

//    winZ = 0.0;

    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

    gcgSETVECTOR3(pos, (double)posX, (double)posY, (double)posZ);

    GLdouble pos3D_x, pos3D_y, pos3D_z;

// arrays to hold matrix information

    gluUnProject(largura/2, altura/2, 0.01,
            modelview, projection, viewport,
            &pos3D_x, &pos3D_y, &pos3D_z);

}

/**************************/

/////////////////////////////////////////////////////////////////////////////////////
// Calcula a quantidade de quadros pos segundo dadas as variacoes de tempo. Atualiza
// a cada 500ms.
/////////////////////////////////////////////////////////////////////////////////////

float CalculaQuadrosSegundo(unsigned int dt) {
    static float ultimoqps = 0.0f; // Guarda o ultimo qps calculado
    static unsigned int contaquadro = 0; // Conta o numero de quadros a cada chamada
    static unsigned int tempopassado = 0; // Armazena o tempo do ultimo quadro

    // Incrementa o contador
    contaquadro++;
    tempopassado += dt;

    // Testa se o tempo passado soma 1 segundo
    if (tempopassado > 1000) {
        // Calcula quadros por segundo
        ultimoqps = ((float) contaquadro / (float) tempopassado) * 1000.f;

        // Zera contador de tempo
        tempopassado = 0;

        // Zera contador de quadros
        contaquadro = 0;
    }

    // Retorna quadros por segundo
    return ultimoqps;
}

void ConfiguraOpenGL() {
//    glClearColor(0.8f, 0.8f, 0.8f, 1.0f); //RGBA
    glClearColor(1.f, 1.f, 1.f, 1.0f); //RGBA
    glShadeModel(GL_SMOOTH);

    //teste de buffer
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_CULL_FACE);
    //glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);


}

// Função callback chamada quando o tamanho da janela é alterado
void alteraTamanhoJanela(GLsizei w, GLsizei h) {
    clearGl = true;
    // Para previnir uma divisão por zero
    if (h == 0) h = 1;

    // Especifica o tamanho da viewport
    glViewport(0, 0, w, h);

    largura = w;
    altura = h;

    float totalDist = farVal + sqrt((frustum.x)*(frustum.x))+((frustum.y)*(frustum.y))+((frustum.z)*(frustum.z));

    frustum.setPerspective(40.0, (float) w / (float) h, 0.1, totalDist);

}

void changeCutBounds(int axis, bool add){
//    printf("opaaaa\n"); return;
    switch(axis){
        case 0:{
            if(add){
                cutBounds1[0] += 0.25;
                if(cutBounds1[0] >= cutBounds2[0]) cutBounds1[0] -= 0.25;
            }
            else{
                cutBounds1[0] -= 0.25;
                if(cutBounds1[0] < 0) cutBounds1[0] = 0;
            }
        }
        break;
        case 1:{
            if(add){
                cutBounds2[0] += 0.25;
                if(cutBounds2[0] > sph->width -1) cutBounds2[0] = sph->width - 1 ;
            }
            else{
                cutBounds2[0] -= 0.25;
                if(cutBounds2[0] < cutBounds1[0]) cutBounds2[0] += 0.25;
            }
        }
        break;
        case 2:{
            if(add){
                cutBounds1[1] += 0.25;
                if(cutBounds1[1] >= cutBounds2[1]) cutBounds1[1] -= 0.25;
            }
            else{
                cutBounds1[1] -= 0.25;
                if(cutBounds1[1] < 0) cutBounds1[1] = 0;
            }
        }
        break;
        case 3:{
            if(add){
                cutBounds2[1] += 0.25;
                if(cutBounds2[1] > sph->height -1) cutBounds2[1] = sph->height - 1 ;
            }
            else{
                cutBounds2[1] -= 0.25;
                if(cutBounds2[1] < cutBounds1[1]) cutBounds2[1] += 0.25;
            }
        }
        break;
        case 4:{
            if(add){
                cutBounds1[2] += 0.25;
                if(cutBounds1[2] >= cutBounds2[2]) cutBounds1[2] -= 0.25;
            }
            else{
                cutBounds1[2] -= 0.25;
                if(cutBounds1[2] < 0) cutBounds1[2] = 0;
            }
        }
        break;
        case 5:{
            if(add){
                cutBounds2[2] += 0.25;
                if(cutBounds2[2] > sph->depth -1) cutBounds2[2] = sph->depth - 1 ;
            }
            else{
                cutBounds2[2] -= 0.25;
                if(cutBounds2[2] < cutBounds1[2]) cutBounds2[2] += 0.25;
            }
        }
        case 9:{
            if(add){
                filterT += 0.02;
                if(filterT > 1) filterT = 1;
            }
            else{
                filterT -= 0.02;
                if(filterT < 0) filterT = 0;
            }
        }
        break;
    }
    
//    printf("X: %f-%f | Y: %f-%f | Z: %f-%f\n", cutBounds1[0], cutBounds2[0], cutBounds1[1], cutBounds2[1], cutBounds1[2], cutBounds2[2]);
}

void teclado(unsigned char tecla, int x, int y) {
    clearGl = true;
    switch (tecla) {
//        case 'w':
//        case 'W': velocidadeFrente = sqrt(frustum.x * frustum.x + frustum.y * frustum.y + frustum.z * frustum.z);
//            break;
//        case 's':
//        case 'S': velocidadeFrente = -sqrt(frustum.x * frustum.x + frustum.y * frustum.y + frustum.z * frustum.z);
//            break;
//        case 'a':
//        case 'A': velocidadeLado = -sqrt(frustum.x * frustum.x + frustum.y * frustum.y + frustum.z * frustum.z);
//            break;
//        case 'd':
//        case 'D': velocidadeLado = sqrt(frustum.x * frustum.x + frustum.y * frustum.y + frustum.z * frustum.z);
//        break;
//        case 'b':
//        case 'B': simula = !simula;
//            break;
        case 'h':
        case 'H': sph->particles[1]->position[1] += 0.2;
            break;
        case 'n':
        case 'N': filterParticles = !filterParticles;
            break;
        case 'b':
        case 'B': teste->changeDrawing();
            break;
        case 'm':
        case 'M': sph->particles[1]->position[0] += 0.2;
//
        case 'q':
//        case 'Q': sph->toggleBox();
//        case 'Q': printDebugVals();
        case 'Q': changeTensorApplication();
        break;
//        case 'w':sph->iterate();
        case 'w':
        case 'W':toggleExternals();
            break;
        case 'e':
        case 'E': sph->resetParticles();
            break;
//        case 'a': enableCreation = !enableCreation;
        case 'a': toggleDistortion();
        case 'A':
            break;
        case 's':
        case 'S': { simula = !simula; simulacao(); simula = !simula;}
            break;
        case 'd':
        case 'D': simula = !simula;
            break;
        case 'f':
        case 'F': sph->changeColourType();
            break;
        case 'p':
            glutPostRedisplay(); // Manda redesenhar a janela
            break;
        case ' ': velocidadeFrente = velocidadeLado = 0;
            break;
        case 't':
        case 'T': translationMode = !translationMode;
            break;
        case 'g':
//        case 'G': clearGl = !clearGl;
        case 'G': original = !original;
            break;
//        case 'o':
//        case 'O': original = !original;
//            break;
        case 'c':
        case 'C': cutVis = !cutVis;
            break;
        case 'r':
        case 'R':
            break;
        case 'v':
        case 'V': drawTensor = !drawTensor;
            break;
        case '7': changeCutBounds(currentCutPlane, false);
        break;
        case '8': changeCutBounds(currentCutPlane, true);
        break;
        case '1': currentCutPlane = 0;printf("CurCutPlane: %d\n", currentCutPlane);
        break;
        case '2': currentCutPlane = 1;printf("CurCutPlane: %d\n", currentCutPlane);
        break;
        case '3': currentCutPlane = 2;printf("CurCutPlane: %d\n", currentCutPlane);
        break;
        case '4': currentCutPlane = 3;printf("CurCutPlane: %d\n", currentCutPlane);
        break;
        case '5': currentCutPlane = 4;printf("CurCutPlane: %d\n", currentCutPlane);
        break;
        case '6': currentCutPlane = 5;printf("CurCutPlane: %d\n", currentCutPlane);
        break;
        case '9': currentCutPlane = 9;printf("CurCutPlane: %d\n", currentCutPlane);
        break;

        
        
//        case '4': if (!cutVis) {
//                     if (flag[0] > 0) flag[0]--;
//                  }else{
//                      if ((cutBound[0] > 0) && (cutBound[0] < cutBound[1] )) cutBound[0]++;
//                      else cutBound[0] = 1;
//                  }
//            break;
//        case '6': if (!cutVis){
//                     if (flag[0] < (((int) fheight)-1))  flag[0]++;
//                  }else{
//                    if ((cutBound[1] > cutBound[0]) && (cutBound[1] < (((int) fheight)))) cutBound[1]--;
//                    else cutBound[1] = ((int) fheight)-1;
//                  }
//            break;
//        case '8': if (!cutVis){
//                        if (flag[1] < (((int) fdepth)-1)) flag[1]++;
//                  }else{
//                        if ((cutBound[3] > cutBound[2]) && (cutBound[3] < (((int) fdepth)))) cutBound[3]--;
//                        else cutBound[3] = ((int) fdepth)-1;
//                  }
//            break;
//        case '2': if (!cutVis){
//                        if (flag[1] > 0) flag[1]--;
//                  }else{
//                        if ((cutBound[2] > 0) && (cutBound[2] < cutBound[3] )) cutBound[2]++;
//                        else cutBound[2] =1;
//                  }
//            break;
//        case '9': if (!cutVis){
//                        if (flag[2] > 0) flag[2]--;
//                  }else{
//                        if ((cutBound[4] > 0) && (cutBound[4] < cutBound[5] )) cutBound[4]++;
//                        else cutBound[4] = 1;
//                  }
//            break;
//        case '1': if (!cutVis){
//                        if (flag[2] < (((int) fwidth) -1)) flag[2]++;
//                  }else{
//                        if ((cutBound[5] > cutBound[4]) && (cutBound[5] < (((int) fwidth)))) cutBound[5]--;
//                        else cutBound[5] = ((int) fwidth)-1;
//                  }
//            break;
//       case '5': if (cutVis){
//
//                cutBound[0] = 1;
//                cutBound[1] = (int) (fheight) - 2;
//                cutBound[2] = 1;
//                cutBound[3] = (int) (fdepth) - 2;
//                cutBound[4] = 1;
//                cutBound[5] = (int) (fwidth) - 2;
//
//
//                  }
//            break;
        case 'x':
        case 'X':{
            simulationStarted = !simulationStarted;
//            SIMULA = true;
        }
        break;
        case 'y':
        case 'Y':{
            sph->insertParticles();

        }
        case 'u':
        case 'U':{
            drawPart = !drawPart;
//            changeK(true);
        }
        break;
        case 'i':
        case 'I':{
            sph->saveState();
//            sph->exportState();
//            changeK(false);
        }
        break;
        case 'o':
        case 'O': changeEditMode();
            break;
        case 'j':
        case 'J':{
//             changeMi(true);
        }
        break;
        case 'k':
        case 'K':{
            changeEditValue(true);
//             changeMi(false);
        }
        break;
        case 'l':
        case 'L': changeEditValue(false);
            break;
        default: break;


    }

    glutPostRedisplay();


}


void tecladoEspecial(int tecla, int x, int y) {
//     switch(tecla){
//        case GLUT_KEY_LEFT:
//            changeEditValue(true);
//            break;
//
//        case GLUT_KEY_RIGHT:
//            changeEditValue(false);
//            break;
//
//        case GLUT_KEY_DOWN:
//            break;
//
//        case GLUT_KEY_UP:
//            changeEditMode();
//            break;
//        default:
//            break;
//     }

}



bool isNullVec(VECTOR3 v1){
    if(FEQUAL(v1[0], 0.0) && FEQUAL(v1[1], 0.0) && FEQUAL(v1[2], 0.0)) return true;
    else return false;
}



// Funcao de simulacao de tempo
void simulacao(void) {
   if (simula && simulationStarted){
//       printf("Simula\n");
    glutSetWindow(janela);
    static clock_t tempoini = 0;
    clock_t tempo;


    gcgRANDOMGAUSSIAN gcgRandom;

    tempo = clock();
    dt = ((float) (tempo - tempoini) / (float) CLOCKS_PER_SEC);
    if (dt < EPSILON) return; // Nothing to do
    if (dt > 0.1) dt = 0.1;

    ////////////////////////////////////////////////////
    // Aqui podemos atualizar TODOS os agentes da cena

    // Atualiza a posicao da camera em movimento. Mas só se for necessario.
    //if(fabs(velocidadeFrente) > 0.1 || fabs(velocidadeLado) > 0.1) {
    float x = frustum.x, y = frustum.x, z = frustum.z; // Save

//    frustum.advancePosition((float) velocidadeFrente * dt, 0.0f, (float) velocidadeLado * dt);
//    if ((frustum.x * frustum.x + frustum.y * frustum.y + frustum.z * frustum.z) < (0.6*0.6)) frustum.setPosition(x, y, z); // restore if too close
//    if (!translationMode) frustum.setTarget(0.f, 0.f, 0.f); // Ajusta as direcoes para visualizar o alvo.

    velocidadeFrente -= velocidadeFrente * dt * 2;
    velocidadeLado -= velocidadeLado * dt * 2;
    if(simulationStarted)
        sph->iterate();
    //Calculo dos parametros dependentes da posição


//    for(int i = 0; i < nParticles; i++) {
//        gcgParticleSPH* p = geometry->particles[i];
//        VECTOR3 vt, at;
//        gcgSCALEVECTOR3(vt, p->velocity, TIME_STEP);
//        gcgSCALEVECTOR3(at, p->acceleration, TIME_STEP2*0.5);
//        gcgADDVECTOR3(p->pos, p->pos, vt);
//        gcgADDVECTOR3(p->pos, p->pos, at);
//        gcgCOPYVECTOR3(p->lastAcceleration, p->acceleration);
//        gcgSETVECTOR3(p->acceleration, 0.0, 0.0, 0.0);
//        gcgCOPYVECTOR3(p->lastVelocity, p->velocity);
//    }
//
//    geometry->particles[0]->parentOctree->checkCollisionBruteForce();
//
//    for(int i = 0; i < nParticles; i++) {
//        gcgParticleSPH* p = geometry->particles[i];
//        VECTOR3 at;
//        gcgADDVECTOR3(at, p->lastAcceleration, p->acceleration);
//        gcgSCALEVECTOR3(at, at, 0.5*TIME_STEP);
//        gcgADDVECTOR3(p->velocity, p->lastVelocity, at);
//    }


    glutPostRedisplay(); //manda redesenhar a janela
        tempoini = tempo;
    }
}

// Função callback chamada para gerenciar eventos do mouse

void gerenciaMouse(int button, int state, int x, int y) {
    clearGl = true;
    gira = FALSE;
    translada = FALSE;



    if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON) {
        iniX = x;
        iniY = y;
        gira = TRUE;
//        oldSimula = simula;
//        simula = FALSE;

        if(enableCreation){
            VECTOR3 p;
            GetOGLPos(x, y, p);
            printf("x = %f | y = %f\n", p[0] , p[1]);
        }
    }

    if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON) {
        iniX = x;
        iniY = y;
        translada = TRUE;
        oldSimula = simula;
        simula = FALSE;

    }

    if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) {
//        simula = TRUE;
    }

    if (state == GLUT_UP && button == GLUT_RIGHT_BUTTON) {
//        simula = TRUE;
        simula = oldSimula;

    }


}

// Função callback chamada para gerenciar eventos do mouse

void gerenciaMovimentoMouse(int x, int y) {
    clearGl = true;
    if (gira) {
        frustum.rotateOrbitTrackball(0, 0, 0, (largura - 2.f * iniX) / largura, (altura - 2.f * y) / altura, (largura - 2.f * x) / largura, (altura - 2.f * iniY) / altura);
        if (!translationMode) frustum.setTarget(0.f, 0.f, 0.f); // Ajusta as direcoes para visualizar o alvo.
        glutPostRedisplay();
    }

    if (translada) {
        frustum.advancePosition(.3 * (x - iniX), 0.0, 0.0);
        glutPostRedisplay();
    }
    iniX = x;
    iniY = y;

}


void cleanupQuadric(void) // Properly Kill The Window
{
    gluDeleteQuadric(IDquadric);
    printf("cleanupQuadric completed\n");
}

void desenhaTudo() {

    if (clearGl)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    frustum.setPerspective(40.0, (float) largura / (float) altura, 0.1, farVal + sqrt((frustum.x)*(frustum.x))+((frustum.y)*(frustum.y))+((frustum.z)*(frustum.z)));
    frustum.exportOpenGL();


    // Desenha eixos
//    glDisable(GL_TEXTURE_2D);
//    glBegin(GL_LINES);
//    glColor3f(0.0f, 1.0f, 0.0f);
//    glVertex3f(0.0f, 0.0f, 0.0f);
//    glVertex3f(0.0f, 10.0f, 0.0f);
//    glColor3f(1.0f, 0.0f, 0.0f);
//    glVertex3f(0.0f, 0.0f, 0.0f);
//    glVertex3f(10.0f, 0.0f, 0.0f);
//    glColor3f(0.0f, 0.0f, 1.0f);
//    glVertex3f(0.0f, 0.0f, 0.0f);
//    glVertex3f(0.0f, 0.0f, 10.0f);
//
//
//    glEnd();


    //setting up the light
    VECTOR4 diffuse; //color
    VECTOR4 specular; //bright
    VECTOR4 ambient;
    VECTOR4 emission;
    int materialSpecularity = 40;
    VECTOR4 specularity;
    VECTOR4 pos;

    gcgSETVECTOR4(diffuse, 1.0, 1.0, 1.0, 1.0);
//    gcgSETVECTOR4(specular, 0.24725, 0.1995, 0.0745, 1.0);
    gcgSETVECTOR4(specular, 0.24725, 0.1995, 0.0745, 1.0);
    gcgSETVECTOR4(specularity, 0.628281, 0.655802, 0.366065, 1.0);
//    gcgSETVECTOR4(specularity, 0.628281, 0.655802, 0.366065, 1.0);
    gcgSETVECTOR4(ambient, 0.0, 0.0, 0.0, 1.0);
    gcgSETVECTOR4(pos, frustum.x, frustum.y, frustum.z, 1.0);
    gcgSETVECTOR4(emission, 0.1f, 0.1f, 0.1f, 1.0f);

    if (!drawMesh) {
        glShadeModel(GL_SMOOTH);
        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_LIGHTING);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specularity);
        glMateriali(GL_FRONT, GL_SHININESS, materialSpecularity);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
                glLightfv(GL_LIGHT0, GL_POSITION, pos);
        glEnable(GL_LIGHT0);
    } else glDisable(GL_LIGHTING);



    //setting up alpha
    //glEnable(GL_BLEND); // Turn Blending On
    //glDisable(GL_BLEND); // Turn Blending On
    //glDisable(GL_LIGHTING); // Turn Blending On
    //glDisable(GL_LIGHT0); // Turn Blending On
    //  glDisable(GL_DEPTH_TEST);	// Turn Depth Testing Off

    glDisable(GL_LIGHTING);


//    sph->drawBox();
    if(drawTensor) teste->drawField();
    if(drawNorms) teste->drawNorms();


    glEnable(GL_LIGHTING);
    if(drawPart){
        if(original) {
            sph->drawParticles();
        }
        else{
            sph->drawParticlesSuperQuadricTassio();
        }
    }


    texto.textPosition(0, 0);

    


    glDisable(GL_LIGHTING);
    
    /********Debug de Distorcao*******/
//    glColor3f(1.0,0.,0.);
//    glTranslatef(-32, -32, 0);
////    gcgDrawVectorPyramid(sph->particles[0]->position[0], sph->particles[0]->position[1], sph->particles[0]->position[2], sph->particles[0]->debugDif, 1.0f);
////    gcgDrawVectorPyramid(0, 0,0, v2, 1.0f);
//    glColor3f(0.0,1.,0.);
////    gcgDrawVectorPyramid(sph->particles[0]->position[0], sph->particles[0]->position[1], sph->particles[0]->position[2], sph->particles[0]->debugV, 1.0f);
//
//
//    glPointSize(10.);
//
//    glMatrixMode(GL_MODELVIEW);
//    glPushMatrix();
////    glLoadIdentity();
////    glTranslatef(3, 0., 0.);
////    glTranslatef(-translateW, -translateH, -translateD);
//    gcgParticleSPH * p = sph->particles[0];
//    gcgParticleSPH * p2 = sph->particles[1];
//    glTranslatef(p->position[0], p->position[1], p->position[2]);
//
//    glLineWidth(3.);
//    
//    VECTOR3 dist, dist2;
//    gcgSUBVECTOR3(dist, p2->position, p->position);
//    
//    glBegin(GL_LINES);
//    glColor3f(1.0,0.,0.);
//    glVertex3f(0.,0.,0.);
//    glVertex3f(dist[0],dist[1],dist[2]);
//    
//    glEnd();
//    
//    MATRIX3 mTensor;
//    gcgCOPYMATRIX3(mTensor, p->tensor->tensor);
//    gcgADDMATRIX3(mTensor, mTensor, p2->tensor->tensor);
//    gcgSCALEMATRIX3(mTensor, mTensor, 0.5);
//    gcgAPPLYMATRIX3VECTOR3(dist2, mTensor, dist);
//    
//    glBegin(GL_LINES);
//    glColor3f(0.0,1.,0.);
//    glVertex3f(0.,0.,0.);
//    glVertex3f(dist2[0],dist2[1],dist2[2]);
//    
//    glEnd();
//    
//    MATRIX3 ti1, ti2, mi;
//    gcgIDENTITYMATRIX3(mi);
//    gcgSUBMATRIX3(ti1, mi, p->tensor->tensor);
//    gcgSUBMATRIX3(ti2, mi, p2->tensor->tensor);
//    gcgADDMATRIX3(mTensor, ti1, ti2);
//    gcgSCALEMATRIX3(mTensor, mTensor, 0.5);
//    gcgAPPLYMATRIX3VECTOR3(dist2, mTensor, dist);
//    
//    glBegin(GL_LINES);
//    glColor3f(0.0,0.,1.);
//    glVertex3f(0.,0.,0.);
//    glVertex3f(dist2[0],dist2[1],dist2[2]);
//    
//    glEnd();
//    
//    printf("Dist: %f - %f\n", gcgLENGTHVECTOR3(dist), gcgLENGTHVECTOR3(dist2));
//    glLineWidth(1.);
//    glPopMatrix();

//    gcgDrawVectorPyramid(sph->particles[0]->position[0], sph->particles[0]->position[1], sph->particles[0]->position[2], v2, 1.0f);
/***********************************/



    glColor3f(0.f, 0.f, 0.f);
    texto.gcgprintf("GCG Particle Tracing for Tensor Visualization");
    clock_t tempo = clock();
    texto.gcgprintf("\nQps: %5.2f", CalculaQuadrosSegundo(((tempo - tempoant) * 1000) / CLOCKS_PER_SEC));
    if (translationMode) texto.gcgprintf("\nTranslation Mode(T) ON");
    else texto.gcgprintf("\nTranslation Mode(T) OFF");
    if (cutVis) texto.gcgprintf("\nCute Mode(C) ON");
    else texto.gcgprintf("\nCut Mode(C) OFF");
    texto.gcgprintf("\nTotal time: %f s", sph->simulationTime);
    texto.gcgprintf("\nTotal iterations: %d", simulationStep);
    texto.gcgprintf("\nHeat Colour: %s", printColour());
    texto.gcgprintf("\nActivate Distortion: %s", printDistortion());
    texto.gcgprintf("\nExternal Forces: %s", printUseExternals());
    texto.gcgprintf("\nInverse Tensor: %s\n", printInverseTensor());
    texto.gcgprintf("\n-- Editing: %s --", printEditMode());
    texto.gcgprintf("\nSupport Radius: %f", printRadius());
//    texto.gcgprintf("\nRest Density: %f", printRestDensity());
    texto.gcgprintf("\nStiffness Coef: %f", printStiffness());
    texto.gcgprintf("\nViscosity Coef: %f", printViscosity());
    texto.gcgprintf("\nGradiente Scale: x%.2f", printGradientScale());
    texto.gcgprintf("\nParticles Radius: %.4f", printParticlesRadius());

    tempoant = tempo;

    //Executa os comandos OpenGL e troca os buffers
    glutSwapBuffers();

}

int main(int argv, char** argc) {

    
    ////Se der pa no import eh o inix y z
    switch(tipoCampo){
        case FIELD2D:
            teste = new gcgTensor(64, 64, 1, 20);
//            sph = new gcgSph(900, 30, 30, 1, 64, 64, 1);
//            sph = new gcgSph(2500, 50, 50, 1, 64, 64, 1);
//            sph = new gcgSph(50, 50, 1, 50, 50, 1, 64, 64, 1);
//            sph = new gcgSph(20, 10, 1, 20, 20, 1, 64, 64, 1);
//            sph = new gcgSph(9, 9, 1, 9, 9, 1, 64, 64, 1);
//            sph = new gcgSph(2, 1, 1, 2, 1, 1, 64, 64, 1);
//            sph = new gcgSph(50, 50, 1, 50, 50, 1, 64, 64, 1);
//            sph = new gcgSph(64, 64, 1, 64, 64, 1, 64, 64, 1);
//            sph = new gcgSph(5, 1, 1, 5, 1, 1, 64, 64, 1);
//            sph = new gcgSph(1, 1, 1, 1, 1, 1, 64, 64, 1);
//            printf("\n\n\n%f %f %f\n%f %f %f\n%f %f %f\n\n\n",
//                    sph->particles[0]->tensor->tensor[0], sph->particles[0]->tensor->tensor[1], sph->particles[0]->tensor->tensor[2],
//                    sph->particles[0]->tensor->tensor[3], sph->particles[0]->tensor->tensor[4], sph->particles[0]->tensor->tensor[5],
//                    sph->particles[0]->tensor->tensor[6], sph->particles[0]->tensor->tensor[7], sph->particles[0]->tensor->tensor[8]);
            sph = loadState(737);
            break;
        case HELICE:
            teste = new gcgTensor(38, 39, 40, 2);
//            sph = new gcgSph(4500, 15, 15, 20, 38, 39, 40);
//            sph = new gcgSph(9000, 15, 30, 20, 38, 39, 40);
//            sph = new gcgSph(27000, 15, 15, 15, 38, 39, 40);
//            sph = new gcgSph(15, 15, 20, 15, 15, 20, 38, 39, 40); //9k
//            sph = new gcgSph(15, 20, 30, 10, 15, 20, 38, 39, 40); //9k
//            sph = new gcgSph(10, 25, 20, 8, 20, 20, 38, 39, 40); //5k
//            sph = new gcgSph(30, 30, 30, 20, 20, 20, 38, 39, 40); //27k
//            sph = new gcgSph(30, 50, 40, 25, 30, 28, 38, 39, 40); //60k
//            sph = new gcgSph(3, 3, 3, 3, 3, 3, 38, 39, 40); //60k
//            sph = new gcgSph(1, 1, 1, 1, 1, 1, 38, 39, 40);
//            printf("\n\n\n%f %f %f\n%f %f %f\n%f %f %f\n\n\n",
//                    sph->particles[0]->tensor->tensor[0], sph->particles[0]->tensor->tensor[1], sph->particles[0]->tensor->tensor[2],
//                    sph->particles[0]->tensor->tensor[3], sph->particles[0]->tensor->tensor[4], sph->particles[0]->tensor->tensor[5],
//                    sph->particles[0]->tensor->tensor[6], sph->particles[0]->tensor->tensor[7], sph->particles[0]->tensor->tensor[8]);
//            sph = new gcgSph(27000, 30, 30, 30, 38, 39, 40);
//            sph = loadState(759);
            sph = loadState(7);
//            sph = loadState(1800);
            break;
        case CEREBRO:
            teste = new gcgTensor(74, 95, 80, 20);
//            sph = new gcgSph(1000, 10, 10, 10, 74, 95, 80);
//            sph = new gcgSph(20, 20, 20, 15, 15, 15, 74, 95, 80);
//            sph = new gcgSph(20, 40, 25, 15.25, 25, 20, 74, 95, 80); //20k
//            sph = new gcgSph(25, 50, 28, 15, 35, 25, 74, 95, 80);
//            sph = new gcgSph(95000, 38, 50, 50, 74, 95, 80);
//            sph = new gcgSph(30000, 25, 40,30, 74, 95, 80);
//            sph = loadState(17);
//            sph = loadState(7);
            sph = loadState(2);
//            sph = loadState(0);
            break;
        case CEREBRO2:
            teste = new gcgTensor(19, 24, 20, 20);
            sph = new gcgSph(125, 5, 5, 5, 19, 24, 20);
            break;
        case CEREBRO1:
            teste = new gcgTensor(37, 48, 40, 20);
            sph = new gcgSph(18000, 30, 30, 20, 37, 48, 40);
//            sph = new gcgSph(42000, 30, 40, 35, 37, 48, 40);
//            sph = new gcgSph(1100, 10, 10, 11, 37, 48, 40);
//            sph = loadState(14);
            break;
        case HELICE1:
            teste = new gcgTensor(5, 19, 20, 20);
            sph = new gcgSph(1100, 10, 10, 11, 19, 20, 20);
            break;
        default:
            break;

    }
//    VECTOR3 c1, c2;
//    gcgSETVECTOR3(c1, 0., 0., 0.);
//    gcgSETVECTOR3(c2, fwidth, fheight, fdepth);
//    _octree = new Octree(c1, c2, 0, NULL);

//    sph = new gcgSph(1, 2, 998.29, 10, 10, 10);
//    sph = new gcgSph(1, 10, 998.29, 5, 5, 5, 10, 10, 20);

//        sph = new gcgSph(343, 7, 7, 7, 38, 39, 40);



    teste->bigger = -1.;
    
//    delete sph;
    
    
//    sph = new gcgSph(1, 2, 3, 4, 38, 39, 40);
//    sph = new gcgSph(20, 5, 4, 1, 64, 64, 1);

//    sph = new gcgSph(2, 2, 1, 1, 64, 64, 1);
//    sph = new gcgSph(400, 20, 20, 1, 64, 64, 1);
    gcgSETVECTOR3(cutBounds2, sph->width, sph->height, sph->depth);
    sph->toggleBox();

//    sph = new gcgSph(0.0315, 192, 998.29, 2, 3, 4, 5, 5, 4);
//    teste = new gcgTensor(74, 95, 80, 1920);
//    allocateAll();
//    simulationStarted = true;


    //setting up glut and opengl
    glutInit(&argv, argc);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(largura, altura);

    glutInitWindowPosition(230, 120);
    janela = glutCreateWindow("GCG - DCC/UFJF - Visualizador básico");
    frustum.setPosition(5.0f, 2.f, 5.0f);
    ConfiguraOpenGL(); //inicializa funções do openGL

    glutReshapeFunc(alteraTamanhoJanela); //chamada cada vez que um tamanho de janela � modificado
    glutDisplayFunc(desenhaTudo); //chamada sempre que a janela precisar ser redesenhada
    glutKeyboardFunc(teclado); //chamada cada vez que uma tecla de c�digo ASCII � pressionado
    glutSpecialFunc(tecladoEspecial); //chamada cada vez que uma tecla de c�digo n�o ASCII � pressionado
    glutMotionFunc(gerenciaMovimentoMouse);
    glutMouseFunc(gerenciaMouse);
    glutIdleFunc(simulacao); //chamada a cada ciclo de glutMainLoop()

    sph->changeColourType();
    glutMainLoop(); //mant�m GLUT em loop, considerando eventos acima

    return 0;

}
