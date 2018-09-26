//////////////////////////////////////////////////////////////////////////////
// gcgTENSORGRID: class for Discrete Wavelet Transforms in Tensor Fields
//////////////////////////////////////////////////////////////////////////////
// Marcelo Bernardes Vieira
// Virgínia Fernandes Mota
// Tássio Knop de Castro
//////////////////////////////////////////////////////////////////////////////

#ifndef GCGGRID
#define GCGGRID

#include "gcg.h"

class gcgTENSORGRID{
  public:
    // Field info and data
    MATRIX3 *data;
    unsigned int    width;
    unsigned int    height;
    unsigned int    depth;

    //unsigned int    rowsize; // Stride in bytes between one row and the next

    /////////////////////////////////////////////////////////////
    // Field construction, initialization and creation.
    /////////////////////////////////////////////////////////////
     gcgTENSORGRID();
     ~gcgTENSORGRID();

    bool createGrid(unsigned int iwidth, unsigned int iheight, unsigned int idepth);
    void destroyGrid();  // Release all resources used by this grid object.
    /////////////////////////////////////////////////////////////
    // Discrete Wavelet Transform
    /////////////////////////////////////////////////////////////
	gcgTENSORGRID *DWTDecomposition(int nH, float *H, int nG, float *G);
    gcgTENSORGRID *DWTComposition(int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG);

    double computeMSEwith(gcgTENSORGRID *tocompare); // Computes the MSE.
    double computePSNRwith(gcgTENSORGRID *tocompare); // Computes the PSNR in dB.


	private:

    // DWT transform
    void removeNewSamplesf(gcgTENSORGRID* newGrid);

    void periodicTransRows(gcgTENSORGRID *dest, int nH, float *H, int nG, float *G,bool addCol,bool addRow, bool addPlane);
    void periodicTransCols(gcgTENSORGRID *dest, int nH, float *H, int nG, float *G,bool addCol,bool addRow, bool addPlane);
    void periodicTransPlanes(gcgTENSORGRID *dest, int nH, float *H, int nG, float *G,bool addCol,bool addRow, bool addPlane);
    void periodicSynthRows(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addCol);
    void periodicSynthCols(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addRow);
    void periodicSynthPlanes(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addPlane);
};

#endif

//////////////////////////////
/*

void gcgTENSORGRID::periodicTransRows(gcgTENSORGRID *dest, int nH, MATRIX3 *H, int nG, MATRIX3 *G,bool addCol,bool addRow) {
    const unsigned int centerH = ((nH % 2) ? (nH - 1) : (nH - 2)) / 2;
    const unsigned int centerG = ((nG % 2) ? (nG - 3) : (nG - 2)) / 2;
    const unsigned int windoww = dest->width;
    const unsigned int windowh = this->height;//(addRow) ? dest->height - 1 : dest->height;
    const unsigned int windowp = dest->depth;
    const bool iseven = windoww % 2 == 0;
    const unsigned int halfw = ((iseven)  ? windoww : windoww + 1) >> 1; // Compute half width
    MATRIX3 s, d, newSample;
    int n = 0;

   //for each depth...
   for (register unsigned int p = 0; p < windowp; p++){
    // For each source row...
    for(register unsigned int y = 0; y < windowh; y++) {
        s = newSample = 0.0f;
        if(addCol) {
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
        }

        //now we can do the convolution over the rows
        for(register unsigned int k = 0, i = 0; i < halfw; i++, k += 2) { // Makes the diadic convolution
          s = d = 0.0f;

          // Realiza a convolução do filtro G (wavelet)
          for (int m = 0; m < nG; m++) {
            n = (int) k + (int) ((int) m - (int) centerG);
            n = (n < 0) ? (windoww + n) : n;
            if (n == (int) (windoww - 1) && addCol) {
              d += G[m] * newSample;
              continue;
            }
            d += G[m] * this->data[SOURCE((n % windoww),y,p)]; // Periodic extension
          }

          // Store result
          dest->data[DEST((i+halfw),y,p)] = d;  //high frequency

          // Realiza a convolução do filtro H (escalamento)
          for (int m = 0; m < nH; m++) {
            n = (int) k + (int) ((int) m - (int) centerH);
            n = (n < 0) ? windoww + n : n;
            if (n == (int) (windoww - 1) && addCol) {
              s += H[m] * newSample;
              continue;
            }
            s += H[m] * this->data[SOURCE((n % windoww),y,p)]; // Periodic extension
          }

          // Store result, Low frequency
          //coef[0] = s;
          dest->data[DEST(i,y,p)] = s;
        }
    }
   }
}

////////////////////////////////////////////////////////////
// Do the decomposition in x axis.
////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicTransCols(gcgTENSORGRID *dest, int nH, MATRIX3 *H, int nG, MATRIX3 *G, bool addCol, bool addRow) {
    const unsigned int centerH = ((nH % 2) ? (nH - 1) : (nH - 2)) / 2;
    const unsigned int centerG = ((nG % 2) ? (nG - 3) : (nG - 2)) / 2;
    const unsigned int windoww = this->width;//(addCol) ? dest->width - 1 : dest->width;
    const unsigned int windowh = dest->height;
    const bool iseven = windowh % 2 == 0;
    const unsigned int halfh = ((iseven)  ? windowh : windowh + 1) >> 1; // Compute half height
    MATRIX3 s, d, newSample;
    int n = 0;
    const unsigned int windowp = dest->depth;

   //for each depth...
   for (register unsigned int p = 0; p < windowp; p++){
    // For each source column...
    for(register unsigned int k = 0; k < windoww; k++) {
        s = newSample = 0.0f;
        newSample = 0;
        //lets compute the new sample of this collumn
        if (addRow) {
          unsigned int  y = windowh - 2;
          for(int m = 0; m < nG; m++) {
            n = (int) y + (int) ((int) m - (int) centerG);
            n = (n < 0) ? windowh + n : n;
            if (n == (int) windowh - 1) {
                newSample = G[m];
                continue;
            }
            s += G[m] * this->data[DEST(k,(n % windowh), p)];//p * dest->width * dest->height + (n % windowh) * this->width + k]; // Periodic extension
          }
          newSample = -(s/newSample);
        }

        //now we continue with the convolution over the collumns
        for(register unsigned int y = 0, i = 0; i < halfh; i++, y += 2) { // Makes the diadic convolution
          s = d = 0.0f;

          //MATRIX3 *coef = &dest->data[DEST(k,i,p)];//p * dest->width * dest->height + i * dest->width + k];//[i * dest->width + k];
          // Realiza a convolução do filtro G (wavelet)
          for (int m = 0; m < nG; m++) {
            n = (int) y +(int) ((int) m - (int) centerG);
            n = (n < 0) ? windowh + n : n;
            if (n == (int) windowh - 1 && addRow) {
              d += G[m] * newSample;
              continue;
            }
            d += G[m] * this->data[SOURCE(k,(n % windowh), p)];//row[p*dest->height*dest->width + (n % windowh) * dest->width + k]; // Periodic extension
          }

          // Store result
          //coef[halfh * dest->width] = d; // High frequency
          dest->data[DEST(k,(i+halfh),p)] = d;
          // Realiza a convolução do filtro H (escalamento)
          for (int m = 0; m < nH; m++) {
            n = (int) y  + (int) ((int) m - (int) centerH);
            n = (n < 0) ? windowh + n : n;
            if (n == (int) (windowh - 1) && addRow) {
              s += H[m] * newSample ;
              continue;
            }
            s += H[m] * this->data[SOURCE(k,(n % windowh), p)];//row[p*this->height*this->width + (n % windowh) * this->width + k]; // Periodic extension
          }

          // Store result
          //coef[0] = s; // Low frequency
          dest->data[DEST(k,i,p)] = s;
        }
    }
  }
}

////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicTransPlanes(gcgTENSORGRID *dest, int nH, MATRIX3 *H, int nG, MATRIX3 *G, bool addCol, bool addRow, bool addPlane) {
    const unsigned int centerH = ((nH % 2) ? (nH - 1) : (nH - 2)) / 2;
    const unsigned int centerG = ((nG % 2) ? (nG - 3) : (nG - 2)) / 2;
    const unsigned int windoww = this->width;//(addCol) ? dest->width - 1 : dest->width;
    const unsigned int windowh = dest->height;//(addRow) ? dest->height - 1 : dest->height;
    const unsigned int windowp = dest->depth;
    const bool iseven = windowp % 2 == 0;
    MATRIX3 s, d, newSample;
    int n = 0;
    const unsigned int halfp = ((iseven)  ? windowp : windowp + 1) >> 1; // Compute half depth

   //for each depth...
   for (register unsigned int i = 0; i < windowh; i++){
    // For each source column...
    for(register unsigned int k = 0; k < windoww; k++) {
        s = newSample = 0.0f;
        newSample = 0;
        //lets compute the new sample of this depth
        if (addPlane) {
          unsigned int  lastplane = windowp - 2;
          for(int m = 0; m < nG; m++) {
            n = (int) lastplane + (int) ((int) m - (int) centerG);
            n = (n < 0) ? windowp + n : n;
            if (n == (int) windowp - 1) {
                newSample = G[m];
                continue;
            }
            s += G[m] * this->data[SOURCE(k,i,(n % windowp))];// * this->width * this->height + i * this->height + k]; // Periodic extension
          }
          newSample = -(s/newSample);
        }

        //now we continue with the convolution over the collumns
        for(register unsigned int pd = 0, p = 0; p < halfp; p++, pd += 2) { // Makes the diadic convolution
          s = d = 0.0f;

         // MATRIX3 *coef = &dest->data[p * dest->width * dest->height + i * dest->width + k];
          // Realiza a convolução do filtro G (wavelet)
          for (int m = 0; m < nG; m++) {
            n = (int) pd +(int) ((int) m - (int) centerG);
            n = (n < 0) ? windowp + n : n;
            if (n == (int) (windowp - 1) && addPlane) {
              d += G[m] * newSample;
              continue;
            }
            d += G[m] * this->data[SOURCE(k,i,(n % windowp))];// * this->height * this->width+ i * this->width + k]; // Periodic extension
          }

          // Store result
          dest->data[DEST(k,i,(p+halfp))] = d;// * dest->width * dest->height + i * dest->width + k] = d;//dest
          //coef[halfp * dest->height * dest->width] = d; // High frequency

          // Realiza a convolução do filtro H (escalamento)
          for (int m = 0; m < nH; m++) {
            n = (int) pd  + (int) ((int) m - (int) centerH);
            n = (n < 0) ? windowp + n : n;
            if (n == (int) (windowp - 1) && addPlane) {
              s += H[m] * newSample ;
              continue;
            }
            s += H[m] * this->data[SOURCE(k,i,(n % windowp))];// * this->height * this->width+ i * this->width + k]; // Periodic extension
          }

          // Store result
          dest->data[DEST(k,i,p)] = s;// * dest->width * dest->height + i * dest->width + k] = s;// Low frequency
        }
    }
  }
}

////////////////////////////////////////////////////////////
// Do the composition in y axis.
////////////////////////////////////////////////////////////
void gcgTENSORGRID::periodicSynthRows(gcgTENSORGRID *dest, int ntH, MATRIX3 *tH, int ntG, MATRIX3 *tG, bool addCol) {
    const int centertH = (int) ((ntH % 2) ? (ntH - 1) : (ntH - 2)) / 2;
    const int centertG = (int) ((ntG % 2) ? (ntG - 3) : (ntG - 2)) / 2;
    const unsigned int windoww = this->width;
    const unsigned int windowh = dest->height;
    const unsigned int windowp = dest->depth;
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
  const unsigned int windoww = dest->width;
  const unsigned int windowh = this->height;
  const bool iseven = windowh % 2 == 0;
  const unsigned int halfh = ((iseven)  ? windowh : windowh + 1) >> 1; // Compute half height
  MATRIX3 s2y, s2y1;
  int n;

  const unsigned int windowp = dest->depth;

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
        //MATRIX3 *pixel = &this->data[p * this->width * this->height + (n % halfh) * width + k]; // Periodic extension
        if (pos % 2 == 0) s2y += tH[m] * this->data[SOURCE(k,(n%halfh), p)];//pixel[0]; // Upsample low pass band
      }

      // Even positions(high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) y - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        //MATRIX3 *pixel = &this->data[p * this->width * this->height + (n % halfh) * width + k]; // Periodic extension
        if (pos % 2 == 0) {
          if (n == (int) halfh - 1 && addRow) continue;
          s2y += tG[m] * (((iseven || n % halfh < halfh - 1) ? this->data[SOURCE(k,((n%halfh) +halfh), p)] : 0));//this->data[p * dest->width * dest->height + ((n % (halfh - 1)) + halfh) * width + k])); // Add detail
        }
      }

      // Store even
      MATRIX3 s2yi  = s2y  + s2y;
      //MATRIX3 *coef = &dest->data[p * dest->width * dest->height + (y << 1) * dest->width + k];
      dest->data[DEST(k,(y << 1), p)] = s2yi;

      // Odd positions (low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) y + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        //MATRIX3 *pixel = &this->data[p * this->width * this->height + (n % halfh) * width + k]; // Periodic extension
        if (pos % 2 != 0) s2y1 += tH[m] * this->data[SOURCE(k,(n%halfh), p)] ;//pixel[0];  // Upsample low pass band
      }

      // Odd positions (high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) y + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfh + n : n;
        //MATRIX3 *pixel = &this->data[p * this->width * this->height + (n % halfh) * width + k]; // Periodic extension
        if (pos % 2 != 0) {
          if (n == (int) halfh - 1 && addRow) continue;
          s2y1 +=tG[m] * this->data[SOURCE(k, ((n%halfh) + halfh), p)];//((int) pixel[halfh * width]); // Add detail
        }
      }

      // Store odd
      //coef[width] = s2y1 + s2y1;
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
  const unsigned int windowp = dest->depth;
  const bool iseven = windowp % 2 == 0;
  const unsigned int halfp = ((iseven)  ? windowp : windowp + 1) >> 1; // Compute half depth
  MATRIX3 s2y, s2y1;
  int n;

for (register unsigned int p = 0; p < halfp; p++) {
   //for each depth...
  for (register unsigned int y = 0; y < windowh; y++){
  // For each source column...
  for (register unsigned int k = 0; k < windoww; k++){

      s2y = s2y1 = 0.0f;

      // Even positions(low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) p - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfp + n : n;
        //MATRIX3 *pixel = &this->data[(n % halfp) * dest->width * dest->height + y* width + k]; // Periodic extension
        if (pos % 2 == 0) s2y += tH[m] * this->data[SOURCE(k,y,(n % halfp))];// * dest->width * dest->height + y* width + k];//pixel[0]; // Upsample low pass band
      }

      // Even positions(high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) p - (int)(((pos % 2) ? pos + 1 : pos) / 2);
        n = (n < 0) ? halfp + n : n;
        //MATRIX3 *pixel = &this->data[SOURCE(k,y,(n % halfp))];// * dest->width * dest->height + y * width + k]; // Periodic extension
        if (pos % 2 == 0) {
          //if (n == (int) halfp - 1 && addPlane) continue;
         int aux = SOURCE(k,y,((n%halfp) +halfp));
         int aux2 = (n%halfp) +halfp;
         //printf("%d ", aux2);

         // s2y += tG[m] * ((iseven || (n % halfp) < (halfp - 1)) ? this->data[SOURCE(k,y,(halfp +(n % halfp)))]: 0);//this->data[SOURCE(k,y,(n % (halfp - 1) + halfp))];//]: this->data[SOURCE(k,y,((n % (halfp - 1)) + halfp))];//  * dest->width * dest->height + y * width + k])); // Add detail
         // s2y += tG[m] * this->data[SOURCE(k,y,(halfp + (n % halfp)))];
         ///O PROBLEMA EH AQUI!!! VAGABUNDO, CANALHA, SEM VERGONHA!
         if (iseven && !(y%2)) s2y += tG[m] * this->data[SOURCE(k,y,((n%halfp) +halfp))];//this->data[p * dest->width * dest->height + ((n % (halfh - 1)) + halfh) * width + k])); // Add detail
         //else if (n > (halfp-2)) s2y += 0;//tG[m] * (((iseven || n % halfp < halfp - 1) ? this->data[SOURCE(k,y,((n%halfp) +halfp))] : 0));//this->data[p * dest->width * dest->height + ((n % (halfh - 1)) + halfh) * width + k])); // Add detail
        }// pixel[halfp * dest->width * dest->height ]
        //((iseven || n % halfh < halfh - 1) ? pixel[halfh * width] : this->data[p * dest->width * dest->height + ((n % (halfh - 1)) + halfh) * width + k])); // Add detail
      }

      // Store even
      MATRIX3 s2yi  = s2y  + s2y;
      //MATRIX3 *coef = &dest->data[( p<<1 ) * dest->width * dest->height + y * dest->width + k];
      dest->data[DEST(k,y,( p<<1 ))] = s2yi;// * dest->width * dest->height + y * dest->width + k][0] = s2yi;

      // Odd positions (low pass)
      for (int m = ntH - 1; m >= 0; m--) {
        register int pos = m - (int) centertH;
        n = (int) p + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfp + n : n;
        //MATRIX3 *pixel = &this->data[(n % halfp) * dest->width * dest->height + y * width + k]; // Periodic extension
        if (pos % 2 != 0) s2y1 += tH[m] * this->data[SOURCE(k,y,(n % halfp))];// * dest->width * dest->height + y * width + k];//pixel[0];  // Upsample low pass band
      }

      // Odd positions (high pass)
      for (int m = ntG - 1; m >= 0; m--) {
        register int pos = m - (int) centertG;
        n = (int) p + 1 - ((pos % 2) ? pos + 1 : pos) / 2;
        n = (n < 0) ? halfp + n : n;
        //MATRIX3 *pixel = &this->data[SOURCE(k,y,(n % halfp))];// * dest->width * dest->height + y * width + k]; // Periodic extension
        if (pos % 2 != 0) {
          if (n == (int) halfp - 1 && addPlane) continue;
          s2y1 += tG[m] * this->data[SOURCE(k,y,((n % halfp) + halfp))];//pixel[halfp * dest->width * dest->height]); // Add detail
        }

      }

      // Store odd
      //coef[dest->height * dest->width] = s2y1 + s2y1;
      dest->data[DEST(k,y,(( p<<1 )+1))] = s2y1 + s2y1;
    }
  }
  }
}
*/
