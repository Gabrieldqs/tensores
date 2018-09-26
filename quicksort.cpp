#include "quicksort.h"

/////////

void mergesort(gcgTENSORGLYPH** a, int low, int high) {
 int i = 0;
 int length = high - low + 1;
 int pivot  = 0;
 int merge1 = 0;
 int merge2 = 0;
 gcgTENSORGLYPH * working[length];

 if(low == high)
  return;

 pivot  = (low + high) / 2;

 mergesort(a, low, pivot);
 mergesort(a, pivot + 1, high);
 
 for(i = 0; i < length; i++)
  working[i] = a[low + i];

 merge1 = 0;
 merge2 = pivot - low + 1;

 for(i = 0; i < length; i++) {
  if(merge2 <= high - low)
   if(merge1 <= pivot - low)
    if(key(working[merge1]) < key(working[merge2]))
     a[i + low] = working[merge2++];
    else
     a[i + low] = working[merge1++];
   else
    a[i + low] = working[merge2++];
  else
   a[i + low] = working[merge1++];
 }
}
//////////////////////////////////


/* sort everything inbetween `low' <-> `high' */
void nqs(gcgTENSORGLYPH** arr, int low, int high, bool verbose) {
 int i = low;
 int j = high;
 gcgTENSORGLYPH* y;
 /* compare value */
 float z = key(arr[(low + high) / 2]);

 /* partition */
 do {
  /* find member above ... */
  while(key(arr[i]) < z) i++;

  /* find element below ... */
  while(key(arr[j]) > z) j--;

  if(i <= j) {
      if(verbose)printf("Trocou nqs!\n");
   /* swap two elements */
   y = arr[i];
   arr[i] = arr[j]; 
   arr[j] = y;
   i++; 
   j--;
  }
 } while(i <= j);

 /* recurse */
 if(low < j) 
  quicksort(arr, low, j);

 if(i < high) 
  quicksort(arr, i, high); 
}

int partition32(Item a[], int l, int r, float* O /* = 0 */)
  { int i = l-1, j = r; Item v = a[r];
    for (;;)
      {
        ++i;
        while (less(a[i], v)) ++i;
        --j;
        while (less(v, a[j])){ --j; if (j == l) break;}
        if (i >= j) break;
        exch(a[i], a[j]);
      }
    exch(a[i], a[r]);
    return i;
  }



void bubbleSort(gcgTENSORGLYPH** arr, int n, bool verbose) {

      bool swapped = true;

      int j = 0;

      gcgTENSORGLYPH* tmp;

      while (swapped) {

            swapped = false;

            j++;

            for (int i = 0; i < n - j; i++) {

                  if (key(arr[i]) > key(arr[i + 1])) {
                      if(verbose)printf("Trocou!\n");
                        tmp = arr[i];

                        arr[i] = arr[i + 1];

                        arr[i + 1] = tmp;

                        swapped = true;

                  }

            }

      }

}

void NquickSort(gcgTENSORGLYPH** arr, int left, int right, float * O) {

      int i = left, j = right;
      gcgTENSORGLYPH* tmp;
      float pivot = key(arr[(left + right) / 2]);
      /* partition */
      while (i <= j) {
            while (key(arr[i]) < pivot)
                  i++;
            while (key(arr[j]) > pivot)
                  j--;
            if (i <= j) {
                printf("Trocou quick!\n");
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      };

      /* recursion */
      if (left < j)
            NquickSort(arr, left, j, O);
      if (i < right)
            NquickSort(arr, i, right, O);
}


void quicksortM32(Item a[], int l, int r, float* O /* = 0 */)
  { int i;
    if (r-l <= M) return;
    exch(a[(l+r)/2], a[r-1]);
    compexch(a[l], a[r-1]);
      compexch(a[l], a[r]);
        compexch(a[r-1], a[r]);
    i = partition32(a, l+1, r-1, O);
    quicksortM32(a, l, i-1, O);
    quicksortM32(a, i+1, r, O);
  }
void sortM32(Item a[], int l, int r)
  {
    quicksortM32(a, l, r, 0);
    insertionSort(a, r);
  }


void partitionQK(gcgTENSORGLYPH** A, int left, int right, int* i, int* j) {
    gcgTENSORGLYPH* x;
    gcgTENSORGLYPH* w;
    *i = left;
    *j = right;
        int divi = (int) ((*i + *j) / 2);
    x = A[divi];
    

    do {
        while (x->dist2Obs > A[*i]->dist2Obs) *i = *i + 1;
        while (x->dist2Obs > A[*j]->dist2Obs) *j = *j - 1;
        if (*i <= *j) {
            w = A[*i];
            A[*i] = A[*j];
            A[*j] = w;
            *i = *i + 1;
            *j = *j - 1;
        }
    } while (*i <= *j);
}

void sortQK(gcgTENSORGLYPH** A, int left, int right) {
    int i, j;
    partitionQK(A, left, right, &i, &j);
    if (left < j) sortQK(A, left, j);
    if (i < right) sortQK(A, i, right);
    return;
}
int partition(gcgTENSORGLYPH **v, int inicio, int fim) {
    int i, j;
    gcgTENSORGLYPH *x, *t = NULL;

    x = v[inicio];
    i = inicio - 1;
    j = fim + 1;


    

//    for (;;) {
//
//        do {
//            j--;
//        } while (v[j]->cp > x->cp);
//        do {
//            i++;
//        } while (v[i]->cp < x->cp);
//        if (i < j) {
//            t = v[i];
//            v[i] = v[j];
//            v[j] = t;
//        } else
//            return j;
//
//    }

        for (;;) {
//menor
        do {
            j--;            
        } while (v[j]->k3 > x->k3);
//        } while (v[j]->dist2Obs > x->dist2Obs);
        do {
            i++;
            
        } while (v[i]->k3 > x->k3);
//} while (v[i]->dist2Obs < x->dist2Obs);
        if (i < j) {
            t = v[i];
            v[i] = v[j];
            v[j] = t;
        } else
            return j;

    }

}

void quicksort(gcgTENSORGLYPH **v, int inicio, int fim) {
    int q;

    if (inicio < fim) {
        q = partition(v, inicio, fim);
        quicksort(v, inicio, q);
        quicksort(v, q + 1, fim);
    }



}

void quickSort(gcgTENSORGLYPH** A, int n) {
    sortQK(A, 0, n - 1);
}


void swap(gcgTENSORGLYPH **a, gcgTENSORGLYPH **b){
     gcgTENSORGLYPH *t = (*a);
     a = b;
     b = &t;
}

gcgTENSORGLYPH* median3( gcgTENSORGLYPH** a, int left, int right )
{
  int center;
  center = (left + right) / 2;
  if (a[left]->k3   > a[center]->k3) swap( &a[left],   &a[center] );
  if (a[left]->k3   > a[right]->k3)  swap( &a[left],   &a[right] );
  if (a[center]->k3 > a[right]->k3)  swap( &a[center], &a[right] );

  /* invariant: a[left] <= a[center] <= a[right] */
     swap( &a[center], &a[right-1] );  /* why? */
     return( a[right-1] );

}

void quicksortM3(gcgTENSORGLYPH** a, int left, int right ) {
  int i, j;
  gcgTENSORGLYPH* pivot;

    pivot = median3( a, left, right );
    i = left; j = right-1;              /* why? */

    for (;;) {
      while (a[++i]->k3 < pivot->k3);
      while (a[--j]->k3 > pivot->k3);
      if (i < j)
        swap( &a[i], &a[j] );
      else
        break;
    }
    swap( &a[i], &a[right-1] );        /* why? */
     quicksortM3( a, left, i-1 );
     quicksortM3( a, i+1, right);

}





void remakeHeap(int left, int right, gcgTENSORGLYPH** A, char type) {
    int i = left;
    int j;
    gcgTENSORGLYPH* x;
    j = i * 2;
    x = A[i];
    while (j <= right) {
        switch (type) {
            case 's':
            case 'S':
                if (j < right) {
                    if (A[j]->cs < A[j + 1]->cs)
                        j++;
                }
                if (x->cs >= A[j]->cs) goto L999;

                break;
            case 'l':
            case 'L':
                if (j < right) {
                    if (A[j]->cl < A[j + 1]->cl)
                        j++;
                }
                if (x->cl >= A[j]->cl) goto L999;
                break;
            case 'p':
            case 'P':
                if (j < right) {
                    if (A[j]->cp < A[j + 1]->cp)
                        j++;
                }
                if (x->cp >= A[j]->cp) goto L999;
                break;
            default:
                if (j < right) {
                    if (A[j]->cp < A[j + 1]->cp)
                        j++;
                }
                if (x->cp >= A[j]->cp) goto L999;
                break;
        }
        A[i] = A[j];
        i = j;
        j = i * 2;
    }
L999:
    A[i] = x;
}

void createHeap(gcgTENSORGLYPH** A, int n, char type) {
    int left;
    left = (n / 2) + 1;
    while (left > 0) {
        left--;
        remakeHeap(left, n, A, type);
    }
}


void swap3(gcgTENSORGLYPH **a, gcgTENSORGLYPH **b) {
  gcgTENSORGLYPH *tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

int partition3(gcgTENSORGLYPH **vec, int left, int right, float obX, float obY, float obZ) {
  int i, j;
  float distJ, distL = 0.0;
//  distJ = sqrt(((obX- vec[j]->pos[0])*(obX - vec[j]->pos[0]))+((obY - vec[j]->pos[1])*(obY - vec[j]->pos[1]))+((obZ -  vec[j]->pos[2])*(obZ -  vec[j]->pos[2])));
//  distL = sqrt(((obX- vec[left]->pos[0])*(obX - vec[left]->pos[0]))+((obY - vec[left]->pos[1])*(obY - vec[left]->pos[1]))+((obZ -  vec[left]->pos[2])*(obZ -  vec[left]->pos[2])));
  i = left;
  for (j = left + 1; j <= right; ++j) {
     distJ = sqrt(((obX- vec[j]->pos[0])*(obX - vec[j]->pos[0]))+((obY - vec[j]->pos[1])*(obY - vec[j]->pos[1]))+((obZ -  vec[j]->pos[2])*(obZ -  vec[j]->pos[2])));
     distL = sqrt(((obX- vec[left]->pos[0])*(obX - vec[left]->pos[0]))+((obY - vec[left]->pos[1])*(obY - vec[left]->pos[1]))+((obZ -  vec[left]->pos[2])*(obZ -  vec[left]->pos[2])));

    if (distJ < distL) {
      ++i;
      swap3(&vec[i], &vec[j]);
    }
  }
  swap3(&vec[left], &vec[i]);

  return i;
}

void quickSort3(gcgTENSORGLYPH **vec, int left, int right, float obX, float obY, float obZ) {
  int r;

  if (right > left) {
    r = partition3(vec, left, right,  obX,  obY,  obZ);
    quickSort3(vec, left, r - 1,  obX,  obY,  obZ);
    quickSort3(vec, r + 1, right,  obX,  obY,  obZ);
  }
}


void swapbubble( gcgTENSORGLYPH **v, int i)
{

    gcgTENSORGLYPH * aux;

   aux=v[i];
   v[i] = v[i+1];
   v[i+1] = aux;

}

  void bubble(gcgTENSORGLYPH **v, int qtd)
{
	int i;
	int trocou;

	do
	{
		qtd--;
		trocou = 0;

		for(i = 0; i < qtd; i++){
                  	if(v[i]->dist2Obs >  v[i+1]->dist2Obs)
			{

				swapbubble(v, i);
				trocou = 1;

			}
                }
	}while(trocou);
}
  void insertionSort(gcgTENSORGLYPH **v, int n)
{
	int i, j;
        gcgTENSORGLYPH *chave;
        int upL = (int) (n * 0.22);
	for(j=1; j<n; j++)
//        for(j=(n-upL); j<n; j++)
	{
		chave = v[j];
		i = j-1;
//		while((i >= 0) && (v[i]->dist2Obs > chave->dist2Obs))
                while((i >= 0) && (v[i]->k3 < chave->k3))
		{
			v[i+1] = v[i];
			i--;
		}
		v[i+1] = chave;
	}
}