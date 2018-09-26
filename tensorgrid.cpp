//////////////////////////////////////////////////////////////////////////////
// gcgTENSORGRID: class for Discrete Wavelet Transforms in Tensor Fields
//////////////////////////////////////////////////////////////////////////////
// Marcelo Bernardes Vieira
// Virgínia Fernandes Mota
// Tássio Knop de Castro
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tensorgrid.h"
#include <math.h>

#define SOURCE(x,y,z) (z) * this->width * this->height + (y) * this->width + (x)
#define DEST(x,y,z) (z) * dest->width * dest->height + (y) * dest->width + (x)
gcgTENSORGRID::gcgTENSORGRID(){
   data = NULL;
   height = 0;
   width = 0;
   depth = 0;
}

gcgTENSORGRID::~gcgTENSORGRID() {
   destroyGrid();
}

////////////////////////////////////////////////////////////
// Image creation. Alocates all necessary space for the image
// indicated by parameters. Initiates with black color and
// a gray scale palette (if bpp < 16).
////////////////////////////////////////////////////////////
bool gcgTENSORGRID::createGrid(unsigned int iwidth, unsigned int iheight, unsigned int idepth) {

  // Release previous resources
  destroyGrid();

	// Image information
	width   = iwidth;
	height  = iheight;
	depth   = idepth;

 	// Create memory for image data
	this->data = (MATRIX3*) malloc(sizeof(MATRIX3) *depth * height * width);
	if(this->data == NULL) return false;
	memset(this->data, 0, depth * height * width * sizeof(MATRIX3));
  return true;
}

////////////////////////////////////////////////////////////
// Release all allocated resources for this image
////////////////////////////////////////////////////////////
void gcgTENSORGRID::destroyGrid() {
  if(data) delete data;
  width = 0;
  height = 0;
  depth = 0;
  data = NULL;
}

////////////////////////////////////////////////////////////
// Computes the Mean Squared Error. Perfect reconstruction
// gives MSE = 0.
////////////////////////////////////////////////////////////
double gcgTENSORGRID::computeMSEwith(gcgTENSORGRID *tocompare){
  if(tocompare == NULL) return INF;
  if(tocompare->data == NULL || tocompare->width == 0 || tocompare->height == 0 || tocompare->depth == 0) return INF;

  // Check equivalence between both images.
  if(tocompare->width != width || tocompare->height != height || tocompare->depth != depth) return INF;

  // Now we compute the MSE for this case
  double sqsum = 0;

  for (unsigned int z = 0; z < depth; z++)
   for (unsigned int y = 0; y < height; y++)
     for (unsigned int x = 0; x < width; x++){
       MATRIX3 sub;
       gcgSUBMATRIX3(sub, data[SOURCE(x,y,z)], tocompare->data[SOURCE(x,y,z)]);
//       sqsum += SQR(sub);
     }
    return sqsum /= (width * height);
}


// Computes the Peak Signal-to-Reconstructed Ratio. Good results vary between
// 20 e 40 db. Perfect reconstruction gives an infinite PSNR.
double gcgTENSORGRID::computePSNRwith(gcgTENSORGRID *tocompare){
  double mse = computeMSEwith(tocompare);

  // Computes the PSNR in decibels = 20log(base10) (255/sqrt(MSE))
  return (mse < EPSILON)? INF : 20 * log10(255.0 / sqrt(mse));
}


////////////////////////////////////////////////////////////
// Decomposes an tensor field into its wavelet components.
////////////////////////////////////////////////////////////
gcgTENSORGRID *gcgTENSORGRID::DWTDecomposition(int nH, float *H, int nG, float *G) {

    //is the width or the height odd?
    bool addCol = (this->width % 2 == 1);
    bool addRow = (this->height % 2 == 1);
    bool addPlane = (this->depth % 2 == 1);
    //do we have to add one row or collumn?if so, we will create new GRIDs with larger row or collumn to do the convolution
    int w, h, p;
    w = this->width  + ((addCol) ? 1 : 0);
    h = this->height + ((addRow) ? 1 : 0);
    p = this->depth + ((addPlane) ? 1 : 0);

    // Creates space for dwt GRID
    gcgTENSORGRID  *dwt = new gcgTENSORGRID(),
             *dwt2 = new gcgTENSORGRID(),
             *dwt3 = new gcgTENSORGRID(),
             *result = new gcgTENSORGRID();

    dwt->createGrid(w,h,p);
    dwt2->createGrid(w,h,p);
    dwt3->createGrid(w,h,p);
    result->createGrid(this->width, this->height, this->depth);

    this->periodicTransRows(dwt, nH, H, nG, G, addCol, addRow, addPlane);
    dwt->periodicTransCols(dwt2, nH, H, nG, G, addCol, addRow, addPlane);
    dwt2->periodicTransPlanes(dwt3, nH, H, nG, G, addCol, addRow, addPlane);

   // Remove new samples
  // for(unsigned int p =0; p < result->depth; p++)
    //for(unsigned int i =0; i < result->height; i++)
      //memcpy(&result->data[p* result->width * result->height + i * result->width], &dwt3->data[p* dwt3->width * dwt3->height +i * dwt3->width], result->width*sizeof(MATRIX3));


   memcpy(result->data, dwt3->data, dwt3->width*dwt3->height*dwt3->depth * sizeof(MATRIX3));

   // frees first GRID
    delete dwt;
    delete dwt2;
    delete dwt3;

    // Returns dwt coefficients
    return result;
}

////////////////////////////////////////////////////////////
// Do the decomposition in y axis.
// The dimensions of the destiny GRID here are always even. If it was odd,
// we are now creating a new GRID with another row and/or collumn.
////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicTransRows(gcgTENSORGRID *dest, int nH, float *H, int nG, float *G,bool addCol,bool addRow, bool addPlane) {
    const unsigned int centerH = ((nH % 2) ? (nH - 1) : (nH - 2)) / 2;
    const unsigned int centerG = ((nG % 2) ? (nG - 3) : (nG - 2)) / 2;
    const unsigned int windoww = dest->width;
    const unsigned int windowh = (addRow) ? dest->height - 1 : dest->height;//this->height;//(addRow) ? dest->height - 1 : dest->height;
    const unsigned int windowp = (addPlane) ? dest->depth - 1 : dest->depth;
    const bool iseven = windoww % 2 == 0;
    const unsigned int halfw = ((iseven)  ? windoww : windoww + 1) >> 1; // Compute half width
    MATRIX3 s, d, newSample;
    int n = 0;

   //for each depth...
   for (register unsigned int p = 0; p < windowp; p++){
    // For each source row...
    for(register unsigned int y = 0; y < windowh; y++) {
        //s = newSample = 0.0f;
        gcgSETMATRIX3(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        gcgSETMATRIX3(newSample, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        /*if(addCol) {
          //lets discover the new sample of this row: k is the last position on the high-pass convolution
          unsigned int k = windoww - 2;
          for (int m = 0; m < nG; m++) {
            n = (int) k + (int) ((int) m - (int) (centerG));
            n = (n < 0) ? (windoww + n): n;
            if(n == (int) (windoww - 1)) {
              newSample = G[m];
              continue;
            }
            s += G[m] * this->data[SOURCE((n % windoww),y,p)]; // Periodic extension
          }
          newSample = -(s / newSample);
        }*/

        //now we can do the convolution over the rows
        for(register unsigned int k = 0, i = 0; i < halfw; i++, k += 2) { // Makes the diadic convolution
          //s = d = 0.0f;
          gcgSETMATRIX3(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          gcgSETMATRIX3(d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          // Realiza a convolução do filtro G (wavelet)
          for (int m = 0; m < nG; m++) {
            n = (int) k + (int) ((int) m - (int) centerG);
            n = (n < 0) ? (windoww + n) : n;
            if ((n == (int) (windoww - 1) && addCol) || (n > windoww-1)) {
             /* MATRIX3 teste;
              gcgIDENTITYMATRIX3(teste);
              gcgSCALEMATRIX3(dest->data[SOURCE((n % windoww),y,p)], teste, 20);
              gcgADDMATRIX3(d,d,dest->data[SOURCE((n % windoww),y,p)]);*/
              ///d += G[m] * newSample;
              continue;
            }
            //d += G[m] * this->data[SOURCE((n % windoww),y,p)]; // Periodic extension
            //gcgCOPYMATRIX3(newSample, this->data[SOURCE((n % windoww),y,p)]);
            gcgSCALEMATRIX3(dest->data[SOURCE((n % windoww),y,p)], this->data[SOURCE((n % windoww),y,p)], G[m]);
            gcgADDMATRIX3(d,d,dest->data[SOURCE((n % windoww),y,p)]);
          }

          // Store result
          //dest->data[DEST((i+halfw),y,p)] = d;  //high frequency
          gcgCOPYMATRIX3(dest->data[DEST((i+halfw),y,p)], d);
          gcgCOPYMATRIX3(newSample, dest->data[SOURCE((n % windoww),y,p)]);

          // Realiza a convolução do filtro H (escalamento)
          for (int m = 0; m < nH; m++) {
            n = (int) k + (int) ((int) m - (int) centerH);
            n = (n < 0) ? windoww + n : n;
            if ((n == (int) (windoww - 1) && addCol) || (n > windoww-1)) {
              ///s += H[m] * newSample;
              continue;
            }
            //s += H[m] * this->data[SOURCE((n % windoww),y,p)]; // Periodic extension
            gcgCOPYMATRIX3(newSample, dest->data[SOURCE((n % windoww),y,p)]);
            gcgSCALEMATRIX3(dest->data[SOURCE((n % windoww),y,p)], this->data[SOURCE((n % windoww),y,p)], H[m]);
            gcgADDMATRIX3(s,s,dest->data[SOURCE((n % windoww),y,p)]);
          }

          // Store result, Low frequency
          //dest->data[DEST(i,y,p)] = s;
          gcgCOPYMATRIX3(dest->data[DEST(i,y,p)], s);
          gcgCOPYMATRIX3(newSample, dest->data[SOURCE((n % windoww),y,p)]);
        }
    }
   }
}

////////////////////////////////////////////////////////////
// Do the decomposition in x axis.
////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicTransCols(gcgTENSORGRID *dest, int nH, float *H, int nG, float *G, bool addCol, bool addRow, bool addPlane) {
    const unsigned int centerH = ((nH % 2) ? (nH - 1) : (nH - 2)) / 2;
    const unsigned int centerG = ((nG % 2) ? (nG - 3) : (nG - 2)) / 2;
    const unsigned int windoww = (addCol) ? dest->width - 1 : dest->width;;//this->width;//(addCol) ? dest->width - 1 : dest->width;
    const unsigned int windowh = dest->height;
    const unsigned int windowp = (addPlane)? dest->depth - 1: dest->depth;
    const bool iseven = windowh % 2 == 0;
    const unsigned int halfh = ((iseven)  ? windowh : windowh + 1) >> 1; // Compute half height
    MATRIX3 s, d, newSample;
    int n = 0;

   //for each depth...
   for (register unsigned int p = 0; p < windowp; p++){
    // For each source column...
    for(register unsigned int k = 0; k < windoww; k++) {
        //s = newSample = 0.0f;
        gcgSETMATRIX3(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        gcgSETMATRIX3(newSample, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        /*if(addRow) {
          //lets compute the new sample of this collumn
          unsigned int  y = windowh - 2;
          for(int m = 0; m < nG; m++) {
            n = (int) y + (int) ((int) m - (int) centerG);
            n = (n < 0) ? windowh + n : n;
            if (n == (int) windowh - 1) {
                newSample = G[m];
                continue;
            }
            s += G[m] * this->data[SOURCE(k,(n % windowh), p)];//p * dest->width * dest->height + (n % windowh) * this->width + k]; // Periodic extension
          }
          newSample = -(s/newSample);
        }
        */
        //now we continue with the convolution over the collumns
        for(register unsigned int y = 0, i = 0; i < halfh; i++, y += 2) { // Makes the diadic convolution
          //s = d = 0.0f;
          gcgSETMATRIX3(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          gcgSETMATRIX3(d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          // Realiza a convolução do filtro G (wavelet)
          for (int m = 0; m < nG; m++) {
            n = (int) y +(int) ((int) m - (int) centerG);
            n = (n < 0) ? windowh + n : n;
            if ((n == (int) (windoww - 1) && addRow) || (n > windoww-1)) {
              ///d += G[m] * newSample;
              continue;
            }
            //d += G[m] * this->data[SOURCE(k,(n % windowh), p)];
            gcgSCALEMATRIX3(dest->data[SOURCE(k,(n % windowh), p)], this->data[SOURCE(k,(n % windowh), p)], G[m]);
            gcgADDMATRIX3(d,d,dest->data[SOURCE(k,(n % windowh), p)]);
          }

          // Store result
          //dest->data[DEST(k,(i+halfh),p)] = d;
          gcgCOPYMATRIX3(dest->data[DEST(k,(i+halfh),p)], d);

          // Realiza a convolução do filtro H (escalamento)
          for (int m = 0; m < nH; m++) {
            n = (int) y  + (int) ((int) m - (int) centerH);
            n = (n < 0) ? windowh + n : n;
            if ((n == (int) (windoww - 1) && addRow) || (n > windoww-1)) {
             /// s += H[m] * newSample;
              continue;
            }
            //s += H[m] * this->data[SOURCE(k,(n % windowh), p)];//row[p*this->height*this->width + (n % windowh) * this->width + k]; // Periodic extension
            gcgSCALEMATRIX3(dest->data[SOURCE(k,(n % windowh), p)], this->data[SOURCE(k,(n % windowh), p)], H[m]);
            gcgADDMATRIX3(s,s,dest->data[SOURCE(k,(n % windowh), p)]);
          }
          // Store result, Low frequency
          gcgCOPYMATRIX3(dest->data[DEST(k,i,p)], s);
        }
    }
  }
}


////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicTransPlanes(gcgTENSORGRID *dest, int nH, float *H, int nG, float *G, bool addCol, bool addRow, bool addPlane) {
    const unsigned int centerH = ((nH % 2) ? (nH - 1) : (nH - 2)) / 2;
    const unsigned int centerG = ((nG % 2) ? (nG - 3) : (nG - 2)) / 2;
    const unsigned int windoww = (addCol) ? dest->width - 1 : dest->width;
    const unsigned int windowh = (addRow) ? dest->height - 1 : dest->height;
    const unsigned int windowp = dest->depth;//this->depth;
    const bool iseven = windowp % 2 == 0;
    const unsigned int halfp = ((iseven)  ? windowp : windowp + 1) >> 1; // Compute half depth
    MATRIX3 s, d, newSample;
    int n = 0;

   //for each row...
   for (register unsigned int i = 0; i < windowh; i++){
    // For each source column...
    for(register unsigned int k = 0; k < windoww; k++) {
        gcgSETMATRIX3(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        gcgSETMATRIX3(newSample, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        /*s = newSample = 0.0f;
        if (addPlane) {
          //lets compute the new sample of this depth
          unsigned int  lastplane = windowp - 2;
          for(int m = 0; m < nG; m++) {
            n = (int) lastplane + (int) ((int) m - (int) centerG);
            n = (n < 0) ? (windowp + n) : n;
            if (n == (int) (windowp - 1)) {
                newSample = G[m];
                continue;
            }
            s += G[m] * this->data[SOURCE(k,i,(n % windowp))];// * this->width * this->height + i * this->height + k]; // Periodic extension
          }
          newSample = -(s/newSample);
        }*/

        //now we continue with the convolution over the collumns
        for(register unsigned int pd = 0, p = 0; p < halfp; p++, pd += 2) { // Makes the diadic convolution
         // s = d = 0.0f;
            gcgSETMATRIX3(s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            gcgSETMATRIX3(newSample, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          // Realiza a convolução do filtro G (wavelet)
          for (int m = 0; m < nG; m++) {
            n = (int) pd +(int) ((int) m - (int) centerG);
            n = (n < 0) ? (windowp + n) : n;
            if ((n == (int) (windoww - 1) && addPlane) || (n > windoww-1)) {
              ///d += G[m] * newSample;
              continue;
            }
            //d += G[m] * this->data[SOURCE(k,i,(n % windowp))];// * this->height * this->width+ i * this->width + k]; // Periodic extension
            gcgSCALEMATRIX3(dest->data[SOURCE(k,i,(n % windowp))], this->data[SOURCE(k,i,(n % windowp))], G[m]);
            gcgADDMATRIX3(d,d,dest->data[SOURCE(k,i,(n % windowp))]);
          }

          // Store result
          //dest->data[DEST(k,i,(p+halfp))] = d;      // High frequency
          gcgCOPYMATRIX3(dest->data[DEST(k,i,(p+halfp))], d);

          // Realiza a convolução do filtro H (escalamento)
          for (int m = 0; m < nH; m++) {
            n = (int) pd  + (int) ((int) m - (int) centerH);
            n = (n < 0) ? windowp + n : n;
            if ((n == (int) (windoww - 1) && addPlane) || (n > windoww-1)) {
              ///s += H[m] * newSample ;
              continue;
            }
            //s += H[m] * this->data[SOURCE(k,i,(n % windowp))];// * this->height * this->width+ i * this->width + k]; // Periodic extension
            gcgSCALEMATRIX3(dest->data[SOURCE(k,i,(n % windowp))], this->data[SOURCE(k,i,(n % windowp))], H[m]);
            gcgADDMATRIX3(s,s,dest->data[SOURCE(k,i,(n % windowp))]);
          }

          // Store result
          //dest->data[DEST(k,i,p)] = s;  //Low frequency
          gcgCOPYMATRIX3(dest->data[DEST(k,i,p)], s);
        }
    }
  }
}

/*
////////////////////////////////////////////////////////////
// Do the composition in y axis.
////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicSynthRows(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addCol) {
    const int centertH = (int) ((ntH % 2) ? (ntH - 1) : (ntH - 2)) / 2;
    const int centertG = (int) ((ntG % 2) ? (ntG - 3) : (ntG - 2)) / 2;
    const unsigned int windoww = this->width;
    const unsigned int windowh = this->height;
    const unsigned int windowp = this->depth;
    const bool iseven = windoww % 2 == 0;
    const unsigned int halfw = ((iseven)  ? windoww : windoww + 1) >> 1; // Compute half width
    MATRIX3 s2k, s2k1;
    int n;

   //for each depth...
   for (register unsigned int p = 0; p < windowp; p++){
    // For each source row...
    for (register unsigned int y = 0; y < windowh; y++) {
        for (register unsigned int k = 0; k < halfw; k++) {
            s2k = s2k1 = 0.0f;

            // Even positions (low pass)
            for (int m = ntH - 1; m >= 0; m--) {
                register int pos = m - (int) centertH;
                n = (int) k - ((pos % 2) ? pos + 1 : pos) / 2;
                n = (n < 0) ? halfw + n : n;
                if (pos % 2 == 0) s2k += tH[m] * this->data[SOURCE((n%halfw),y,p)]; // Upsample low pass band
            }

            // Even positions (high pass)
            for (int m = ntG - 1; m >= 0; m--) {
                register int pos = m - (int) centertG;
                n = (int) k - ((pos % 2) ? pos + 1 : pos) / 2;
                n = (n < 0) ? halfw + n : n;
                if (pos % 2 == 0) {
                    if (n == (int) (halfw - 1) && addCol) continue;
                    s2k += tG[m] * this->data[SOURCE(((n%halfw) + halfw),y,p)]; // Add detail
                }
            }

            // Store even
            MATRIX3 s2ki  = s2k  + s2k;
            dest->data[DEST((k << 1), y, p)] = s2ki;

            // Odd positions (low pass)
            for (int m = ntH - 1; m >= 0; m--) {
                register int pos = m - (int) centertH;
                n = (int) k + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
                n = (n < 0) ? halfw + n : n;
                if (pos % 2 != 0) s2k1 += tH[m] * this->data[SOURCE((n%halfw),y,p)];  // Upsample low pass band
            }

            // Odd positions (high pass)
            for (int m = ntG - 1; m >= 0; m--) {
                register int pos = m - (int) centertG;
                n = (int) k + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
                n = (n < 0) ? halfw + n: n;
                if (pos % 2 != 0) {
                    if (n == (int) halfw - 1 && addCol) continue;
                    s2k1 += tG[m] * this->data[SOURCE(((n%halfw) + halfw),y,p)]; // Add detail
                }
            }

            // Store odd
            dest->data[DEST(((k << 1) + 1), y, p)] = s2k1 + s2k1;
        }
     }
    }
}


////////////////////////////////////////////////////////////
// Do the composition in x axis.
////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicSynthCols(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addRow) {
  const int centertH = (int) ((ntH % 2) ? (ntH - 1) : (ntH - 2)) / 2;
  const int centertG = (int) ((ntG % 2) ? (ntG - 3) : (ntG - 2)) / 2;
  const unsigned int windoww = this->width;
  const unsigned int windowh = this->height;
  const unsigned int windowp = this->depth;
  const bool iseven = windowh % 2 == 0;
  const unsigned int halfh = ((iseven)  ? windowh : windowh + 1) >> 1; // Compute half height
  MATRIX3 s2y, s2y1;
  int n;

   //for each depth...
  for (register unsigned int p = 0; p < windowp; p++){
  // For each source column...
  for (register unsigned int k = 0; k < windoww; k++){
    for (register unsigned int y = 0; y < halfh; y++) {
      s2y = s2y1 = 0.0f;

      // Even positions(low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) y - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        if (pos % 2 == 0) s2y += tH[m] * this->data[SOURCE(k,(n%halfh), p)]; // Upsample low pass band
      }

      // Even positions(high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) y - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        if (pos % 2 == 0) {
          if (n == (int) (halfh - 1) && addRow) continue;
          s2y += tG[m] * (((iseven || n % halfh < halfh - 1) ? this->data[SOURCE(k,((n%halfh) +halfh), p)] : 0)); // Add detail
        }
      }

      // Store even
      MATRIX3 s2yi  = s2y  + s2y;
      dest->data[DEST(k,(y << 1), p)] = s2yi;

      // Odd positions (low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) y + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        if (pos % 2 != 0) s2y1 += tH[m] * this->data[SOURCE(k,(n%halfh), p)] ;//pixel[0];  // Upsample low pass band
      }

      // Odd positions (high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) y + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        if (pos % 2 != 0) {
          if (n == (int) halfh - 1 && addRow) continue;
          s2y1 +=tG[m] * this->data[SOURCE(k, ((n%halfh) + halfh), p)]; // Add detail
        }
      }

      // Store odd
      dest->data[DEST(k,((y << 1) +1), p)] = s2y1 + s2y1;
    }
  }
  }
}



void gcgTENSORGRID::periodicSynthPlanes(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addPlane) {
  const int centertH = (int) ((ntH % 2) ? (ntH - 1) : (ntH - 2)) / 2;
  const int centertG = (int) ((ntG % 2) ? (ntG - 3) : (ntG - 2)) / 2;
  const unsigned int windoww = this->width;
  const unsigned int windowh = this->height;
  const unsigned int windowp = this->depth;
  const bool iseven = windowp % 2 == 0;
  const unsigned int halfp = ((iseven)  ? windowp : windowp + 1) >> 1; // Compute half depth
  MATRIX3 s2p, s2p1;
  int n;

  //for each source row..
  for (register unsigned int y = 0; y < windowh; y++){
   // For each source column...
   for (register unsigned int k = 0; k < windoww; k++){
    //for each depth...
    for (register unsigned int p = 0; p < halfp; p++) {
      s2p = s2p1 = 0.0f;

      // Even positions(low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) p - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfp + n : n;
        if (pos % 2 == 0) s2p += tH[m] * this->data[SOURCE(k,y,(n % halfp))]; // Upsample low pass band
      }

      // Even positions(high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) p - (int)(((pos % 2) ? pos + 1 : pos) / 2);
        n = (n < 0) ? halfp + n : n;
        if (pos % 2 == 0) {
         if (n == (int) (halfp - 1) && addPlane) continue;
          s2p += tG[m] * ((iseven || (n % halfp) < (halfp - 1)) ? this->data[SOURCE(k,y,(halfp +(n % halfp)))]: 0); // Add detail
        }
      }

      // Store even
      MATRIX3 s2pi  = s2p  + s2p;
      dest->data[DEST(k,y,( p<<1 ))] = s2pi;

      // Odd positions (low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) p + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfp + n : n;
        if (pos % 2 != 0) s2p1 += tH[m] * this->data[SOURCE(k,y,(n % halfp))]; // Upsample low pass band
      }

      // Odd positions (high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) p + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfp + n : n;
        if (pos % 2 != 0) {
          if (n == (int) halfp - 1 && addPlane) continue;
          s2p1 += tG[m] * ((iseven || (n % halfp) < (halfp - 1)) ? this->data[SOURCE(k,y,((n % halfp) + halfp))]: 0);// Add detail
        }
      }
      // Store odd
      dest->data[DEST(k,y,(( p<<1 )+1))] = s2p1 + s2p1;
    }
  }
  }
}


////////////////////////////////////////////////////////////
// Composes an image from its wavelet components.
////////////////////////////////////////////////////////////
gcgTENSORGRID *gcgTENSORGRID::DWTComposition(int nscales, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG) {
  //is the width or the height odd?
  bool addCol = (this->width % 2 == 1);
  bool addRow = (this->height % 2 == 1);
  bool addPlane = (this->depth % 2 == 1);

  //do we have to add one row or collumn?
  int w, h, p;
  w = this->width  + ((addCol) ? 1 : 0);
  h = this->height + ((addRow) ? 1 : 0);
  p = this->depth + ((addPlane) ? 1 : 0);

  // Creates space for dwt GRID
   gcgTENSORGRID  *dwt = new gcgTENSORGRID(),
            *dwt2 = new gcgTENSORGRID(),
            *dwt3 = new gcgTENSORGRID(),
            *result = new gcgTENSORGRID();

   dwt->createGrid(w,h,p);
   dwt2->createGrid(w,h,p);
   dwt3->createGrid(w,h,p);
   result->createGrid(this->width, this->height, this->depth);

   this->periodicSynthPlanes(dwt, ntH, tH, ntG, tG, addPlane);

   gcgIMAGE *codificada = new gcgIMAGE();
   codificada->createImage(dwt->width, dwt->height, 8, false);
   codificada->data = new unsigned char [codificada->rowsize * codificada->height];
   char nome_arquivo[11];
   for (unsigned int z = 0; z < dwt->depth; z++){
     for (unsigned int y = 0; y < dwt->height; y++)
       for (unsigned int x = 0; x <dwt->width; x++){
           codificada->data[y*codificada->rowsize+ x] = (unsigned char) dwt->data[dwt->width*dwt->height*z + y*dwt->width + x];
        }
       sprintf(nome_arquivo, "%s%d%s", "temp/dcoded-p",z, ".bmp");
       codificada->saveBMP(nome_arquivo);
      }

   dwt->periodicSynthCols(dwt2, ntH, tH, ntG, tG, addRow);
   codificada->createImage(dwt2->width, dwt2->height, 8, false);
   codificada->data = new unsigned char [codificada->rowsize * codificada->height];
   for (unsigned int z = 0; z < dwt2->depth; z++){
     for (unsigned int y = 0; y < dwt2->height; y++)
       for (unsigned int x = 0; x <dwt2->width; x++){
           codificada->data[y*codificada->rowsize+ x] = (unsigned char) dwt2->data[dwt2->width*dwt2->height*z + y*dwt2->width + x];
          }
       sprintf(nome_arquivo, "%s%d%s", "temp/dcoded-c",z, ".bmp");
       codificada->saveBMP(nome_arquivo);
       }

   dwt2->periodicSynthRows(dwt3, ntH, tH, ntG, tG, addCol);
   codificada->createImage(dwt3->width, dwt3->height, 8, false);
   codificada->data = new unsigned char [codificada->rowsize * codificada->height];
   for (unsigned int z = 0; z < dwt3->depth; z++){
     for (unsigned int y = 0; y < dwt3->height; y++)
       for (unsigned int x = 0; x <dwt3->width; x++){
           codificada->data[y*codificada->rowsize+ x] = (unsigned char) dwt3->data[dwt3->width*dwt3->height*z + y*dwt3->width + x];
        }
        sprintf(nome_arquivo, "%s%d%s", "temp/dcoded-r",z, ".bmp");
        codificada->saveBMP(nome_arquivo);
   }

  // Remove new samples
  for(unsigned int p = 0; p < result->depth; p++)
    for(unsigned int i =0; i < result->height; i++)
      memcpy(&result->data[p* result->width * result->height + i * result->width], &dwt3->data[p* dwt3->width * dwt3->height + i * dwt3->width], result->width*sizeof(MATRIX3));

  // frees first image
  delete dwt;
  delete dwt2;
  delete dwt3;

  // Returns the image reconstructed from dwt coefficients
  return result;
}
*/

