#ifndef GCGTENSORGENERATE
#define GCGTENSORGENERATE

#include <stdio.h>
#include <stdlib.h>
#include "gcg.h"
#include "tensorgrid.h"


class gcgGCGTENSORGENERATE{

  public:

      gcgGCGTENSORGENERATE(int width, int height, int depth );
     ~gcgGCGTENSORGENERATE();

     gcgTENSORGRID *grid;



  private:

};

#endif

