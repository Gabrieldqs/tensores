//Molecular Dynamics
//---Codigo Sequencial
//Alessandra Matos Campos
//Joao Paulo Peçanha
//Tassio Knop de Castro


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gl\gl.h>
#include <gl\glu.h>


#include "gcg.h"
#include "volume.h"

#define PARTICLEADDRESS(a, b, c) ((c) * width * height + (b) * width + (a))

#define ISOBJECT(a, b, c) (((a+0.5) - (float) width*0.5) * ((a+0.5) - (float) width*0.5) + ((b+0.5) - (float) height*0.5) * ((b+0.5) - (float) height*0.5) > 10.5  \
                        && ((a+0.5) - (float) width*0.5) * ((a+0.5) - (float) width*0.5) + ((b+0.5) - (float) height*0.5) * ((b+0.5) - (float) height*0.5) < 38.5)


PARTICLE *matrix;                   // volume contendo todas as particulas

unsigned int length, width, height; // Dimensoes do volume de particulas
unsigned int volume;                // Numero inicial de particulas
VECTOR3 flippedSpin;                // Armazena o spin antes de flipar
FLOAT   prevEnergy = INF;           // saves last energy
unsigned int index = 0;             // Last flipped particle
VECTOR3 magneticVector = {0, 0, 0}; // Object's magnetism

FLOAT minenergy = INF, maxenergy = -INF;


unsigned int initMatrix(unsigned int w,unsigned int h, unsigned int l, VECTOR3 *initialSpin) {
    volume = l * w * h;
    matrix = (PARTICLE*) malloc(sizeof(PARTICLE) * volume);
    length = l;
    width  = w;
    height = h;

    if(initialSpin != NULL)
      for(unsigned int i = 0; i < volume; i++) {
        gcgCOPYVECTOR3(matrix[i].spin, *initialSpin);
        matrix[i].energy = 0.0;
      }
    else
      for(unsigned int i = 0; i < volume; i++) {
        gcgRandomVector(matrix[i].spin);
        matrix[i].energy = 0.0;
      }

    unsigned int nparticles = 0;
    for(unsigned int i = 0; i < width; i++)
      for(unsigned int j = 0; j < height; j++)
        for(unsigned int k = 0; k < length; k++)
          if(ISOBJECT(i, j, k)) nparticles++;

    return nparticles;
}

void changeSize(unsigned int newLength,unsigned int newWidth) {
   //arredonda a nova altura para cima, visando uma aproximacao do volume inicial

   unsigned int newHeight = (unsigned int) ceil((float) volume / (newLength * newWidth));
   if(newHeight == 0) newHeight = 1;

   matrix = (PARTICLE*) realloc(matrix, newLength * newWidth * newHeight);
   length = newLength;
   width  = newWidth;
   height = newHeight;
}


// Calcula a energia de todo o sistema, media: media das nIterações sobre o sistema
FLOAT computeEnergy(float A, float J, float D, VECTOR3 magneticField, VECTOR3 magneticVector) {
    unsigned int i, j, k;
    FLOAT totalenergy = 0.0;

    // Comment to get the min/max energy from all iterations
    //minenergy =  INF;
    //maxenergy = -INF;
    gcgZEROVECTOR3(magneticVector);
    for(i = 0; i < width; i++)
      for(j = 0; j < height; j++)
        for(k = 0; k < length; k++)
          if(ISOBJECT(i, j, k)) {
            unsigned int ind = PARTICLEADDRESS(i, j, k);
            matrix[ind].energy = A * computeDipoleDipoleEnergy(i, j, k) - J * computeMagnecticFerFactor(i, j, k) - D * gcgDOTVECTOR3(matrix[ind].spin, magneticField);
            totalenergy += matrix[ind].energy;
            if(matrix[ind].energy < minenergy) minenergy = matrix[ind].energy;
            if(matrix[ind].energy > maxenergy) maxenergy = matrix[ind].energy;
            gcgADDVECTOR3(magneticVector, magneticVector, matrix[ind].spin);
          }

    return totalenergy;
}

FLOAT computeDipoleDipoleEnergy(unsigned int i, unsigned int j, unsigned int k) {
    FLOAT energy = 0.0;
    unsigned int  x = 0, y = 0, z = 0;

    // Current spin
    VECTOR3 *spin = &matrix[PARTICLEADDRESS(i, j, k)].spin;

    for(x = 0; x < width; x++)
      for(y = 0; y < height; y++)
        for(z = 0; z < length; z++)
          if(ISOBJECT(i, j, k) && x != i || y != j || z != k) {
            VECTOR3 *spin2 = &matrix[PARTICLEADDRESS(x, y, z)].spin;
            VECTOR3 dist = {i - x, j - y, k - z};
            FLOAT distance = gcgLENGTHVECTOR3(dist);
            FLOAT cube = (distance * distance * distance);

            FLOAT term1 = gcgDOTVECTOR3(*spin, *spin2) / cube;
            FLOAT term2 = (gcgDOTVECTOR3(*spin, dist) * gcgDOTVECTOR3(*spin2, dist)) / (cube * distance * distance);

            energy += term1 - 3.0 * term2;
          }
    return 0.5 * energy;
}


#define NEIGHBORS 6
FLOAT computeMagnecticFerFactor(unsigned int i, unsigned int j, unsigned int k) {
    static const int neighX[NEIGHBORS] = {-1, 1,   0, 0, 0,  0};
    static const int neighY[NEIGHBORS] = {0,  0,  -1, 1, 0,  0};
    static const int neighZ[NEIGHBORS] = {0,  0,   0, 0, -1, 1};

    VECTOR3  *spinParticle = &matrix[PARTICLEADDRESS(i, j, k)].spin;
    FLOAT     sum = 0.0;

    for(int v = 0; v < NEIGHBORS; v++)
      if(i + neighX[v] >= 0 && i + neighX[v] < width  &&
         j + neighY[v] >= 0 && j + neighY[v] < height &&
         k + neighZ[v] >= 0 && k + neighZ[v] < length)
          if(ISOBJECT(i + neighX[v], j + neighY[v], k + neighZ[v]))
            sum += gcgDOTVECTOR3(*spinParticle, matrix[PARTICLEADDRESS(i + neighX[v], j + neighY[v], k + neighZ[v])].spin);

    return sum;
}


void flipRandomParticle() {
  unsigned int major = MAX(MAX(length, height), width);

  unsigned int x, y, z;

  do
    do {
      x = rand() % major;
      y = rand() % major;
      z = rand() % major;
    } while(x >= width || y >= height || z >= length);
  while(!ISOBJECT(x, y, z));

  // Save spin
  index = PARTICLEADDRESS(x, y, z);
  gcgCOPYVECTOR3(flippedSpin, matrix[index].spin);  // Save old spin
  // Flip spin
  gcgRandomVector(matrix[index].spin);

}

int iterateMetropolis(FLOAT *energy, VECTOR3 resultMag, FLOAT A, FLOAT J, FLOAT D, FLOAT Kt, VECTOR3 magneticField){
  FLOAT curEnergy = computeEnergy(A, J, D, magneticField, resultMag);

  if(curEnergy > prevEnergy)
    if(((float) rand() / (float) RAND_MAX) >= exp(-(curEnergy - prevEnergy) / Kt)) {
      gcgCOPYVECTOR3(matrix[index].spin, flippedSpin); // rolls  back last flip
      *energy = prevEnergy;
      gcgCOPYVECTOR3(resultMag, magneticVector);
      return FALSE;
    }

  // Commit flip
  *energy = prevEnergy = curEnergy;
  gcgCOPYVECTOR3(magneticVector, resultMag);

  return TRUE;
}

void freeMatrix(){
    free(matrix);
}


void gcgHeatColor(float normheat, VECTOR3 color) {
  #define GCG_HEAT_COLORS 6
  const static float colors[GCG_HEAT_COLORS][3] = { {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 0.0},
                                                    {1.0, 1.0, 0.0}, {1.0, 0.5, 0.0}, {1.0, 0.0, 0.0} };
  FLOAT colorpos = normheat * (GCG_HEAT_COLORS - 1);
  int   index = (int) colorpos;
  if(index < 0) {gcgCOPYVECTOR3(color, colors[0]); }
  else if(index >= (GCG_HEAT_COLORS - 1)) { gcgCOPYVECTOR3(color, colors[GCG_HEAT_COLORS - 1]); }
       else {
         FLOAT  s1 = colorpos - (float) index, s0 = 1.0 - s1;
         gcgSETVECTOR3(color, colors[index][0] * s0 + colors[index + 1][0] * s1,
                              colors[index][1] * s0 + colors[index + 1][1] * s1,
                              colors[index][2] * s0 + colors[index + 1][2] * s1 );
       }
}


void gcgComputeAlignMatrix(MATRIX4 matrix, VECTOR3 dir) {
  VECTOR3 v;

  gcgNORMALIZEVECTOR3(v, dir);
  FLOAT nlength = SQR(v[0]) + SQR(v[1]);
  if(nlength + SQR(v[2] - 1.0) > EPSILON * EPSILON) { // diff = V - Z
    VECTOR3 n;
    VECTOR4 quat;

    // Find normal and angle cosine
    // Check normal length to avoid numeric errors
    if(nlength < 10 * EPSILON * EPSILON) {
      gcgSETVECTOR3(n, 0.0, 1.0, 0.0);
    } else gcgSETVECTOR3(n, v[1] * (nlength = 1.0 / sqrt(nlength)), -v[0] * nlength, 0.0); // n = V x Z

    // Check angle: theta = V . Z
    if(v[2] > 1.0) v[2] = 1.0;
    else if(v[2] < -1.0) v[2] = -1.0;

    gcgAXISTOQUATERNION(quat, acos(v[2]), n);
    gcgQUATERNIONTOMATRIX4(matrix, quat);
  } else gcgIDENTITYMATRIX4(matrix);
}


void gcgDrawArrow(FLOAT x, FLOAT y, FLOAT z, VECTOR3 vector, float scale) {
  static unsigned int arrowlist = 0;

  if(arrowlist == 0) {
    arrowlist = glGenLists(1);

    glNewList(arrowlist, GL_COMPILE);
      glBegin(GL_TRIANGLE_FAN);
        glNormal3f(0.0, 0.0, 1.0);
        glVertex3f(0.0,  0.0, 0.5);
        glNormal3f(1.0, 0.0, 0.0);
        glVertex3f(0.11, 0.0, -0.5);
        glNormal3f(0.0, 1.0, 0.0);
        glVertex3f(0.0,  0.11, -0.5);
        glNormal3f(-1.0, 0.0, 0.0);
        glVertex3f(-0.11, 0.0, -0.5);
        glNormal3f(0.0, -1.0, 0.0);
        glVertex3f(0.0,  -0.11, -0.5);
        glNormal3f(1.0, 0.0, 0.0);
        glVertex3f(0.11, 0.0, -0.5);
      glEnd();
    glEndList();
  }

	// Compute rotation matrix
	MATRIX4 matrix;
  gcgComputeAlignMatrix(matrix, vector);

  glListBase(0);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glTranslatef(x, y, z);
  glMultMatrixf(matrix);
  glScalef(scale, scale, scale);
  glCallList(arrowlist);
	glPopMatrix();
}

void drawSpins() {
  // Draw spin vectors
  glDisable(GL_CULL_FACE);
  glPushMatrix();
  glTranslatef((float) width * -0.5f, (float) height * -0.5f, (float) length * -0.5f);

  register unsigned int x, y, z;
  for(x = 0; x < width; x++)
    for(y = 0; y < height; y++)
      for(z = 0; z < length; z++)
        if(ISOBJECT(x, y, z)) {
          VECTOR3 color;
          gcgHeatColor((matrix[PARTICLEADDRESS(x, y, z)].energy - minenergy) / (maxenergy - minenergy), color);
          glColor3f(color[0], color[1], color[2]);
          gcgDrawArrow(x + 0.5, y + 0.5, z + 0.5, matrix[PARTICLEADDRESS(x, y, z)].spin, 0.5);
        }
  glPopMatrix();
}


void drawGrid() {
  // Draw grid
  glDisable(GL_CULL_FACE);
  glPushMatrix();
  glTranslatef((float) width * -0.5f, (float) height * -0.5f, (float) length * -0.5f);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.2, 1.0, 0.2, 0.04);//gcgBlinkingAlpha()*0.08);

  glBegin(GL_LINES);
  for(unsigned int z = 0; z <= length; z++)
    for(unsigned int x = 0; x <= width; x++)
      for(unsigned int y = 0; y <= height; y++) {
          glVertex3f(x, y, z);
          glVertex3f(x, height, z);

          glVertex3f(x, y, z);
          glVertex3f(width, y, z);

          glVertex3f(x, y, z);
          glVertex3f(x, y, length);
      }
  glEnd();

  glPopMatrix();
}



void drawMagneticField(VECTOR3 mag) {
  // Draw magnetic field
  glDisable(GL_CULL_FACE);
  glColor4f(1.0, 1.0, 1.0, gcgBlinkingAlpha() * 0.18);
  gcgDrawArrow(0.0, 0.0, 0.0, mag, 6);

  // Draw magnetic field
  glColor4f(1.0, 0.0, 1.0, gcgBlinkingAlpha() * 0.18);
  gcgDrawArrow(0.0, 0.0, 0.0, magneticVector, 6);

  glDisable(GL_BLEND);

  // Draw text boxes
  static gcgTEXT box;
  MATRIX4 matrix;
  glPushMatrix();
  gcgComputeAlignMatrix(matrix, mag);
  glMultMatrixf(matrix);
  box.enableTextBoxAt3DPos(0, 0, 7, 60, 16);
  box.adjustTextBox(-30, -8, 0, 0);
  box.drawTextBox(1.0, 1.0, 1.0, gcgBlinkingAlpha() * 0.18, 1.0, 1.0, 1.0, gcgBlinkingAlpha() * 0.28, 1.0);
  glColor4f(1.0, 1.0, 1.0, gcgBlinkingAlpha());
  box.setSystemFont(GCG_FONT_SANSSERIF_11_NORMAL);
  box.gcgprintf("External field");
  glColor4f(1.0, 1.0, 1.0, gcgBlinkingAlpha());
  glPopMatrix();

  glPushMatrix();
  gcgComputeAlignMatrix(matrix, magneticVector);
  glMultMatrixf(matrix);
  box.enableTextBoxAt3DPos(0, 0, 7, 68, 16);
  box.adjustTextBox(-34, -8, 0, 0);
  box.drawTextBox(1.0, 0.0, 1.0, gcgBlinkingAlpha() * 0.18, 1.0, 0.0, 1.0, gcgBlinkingAlpha() * 0.28, 1.0);
  glColor4f(1.0, 1.0, 1.0, gcgBlinkingAlpha());
  box.gcgprintf("Magnetization");
  glPopMatrix();
}
