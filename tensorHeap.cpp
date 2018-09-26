#include "tensorHeap.h"


void Refaz(int Esq, int Dir, gcgTENSORGLYPH** A)
{ 
  char i = Esq;
  int j;
  gcgTENSORGLYPH* x;
  j = i * 2;
  x = A[i];
  while (j <= Dir)
    { if (j < Dir)
       
      { if (A[j]->cs < A[j+1]->cs)
        j++;
      }
      if (x->cs >= A[j]->cs) goto L999;
      A[i] = A[j];
      i = j; j = i * 2;
    }
  L999: A[i] = x;
}
void constroiHeap(gcgTENSORGLYPH** A, int n)
{ int Esq;
  Esq = n / 2 + 1;
  while (Esq > 1)
    { Esq--;
      Refaz(Esq, n, A);
    }
}

//
//Item Max(Item *A)
//{ return (A[1]);
//}
//
//Item RetiraMax(Item *A, Indice *n)
//{ Item Maximo;
//  if (*n < 1)
//  printf(" Erro: heap vazio \n");
//  else
//    { Maximo = A[1];
//      A[1]   = A[*n];
//      (*n)--;
//      Refaz(1, *n, A);
//    }
//  return Maximo;
//}
//void AumentaChave(Indice i, ChaveTipo ChaveNova, Item *A)
//{ Item x;
//  if (ChaveNova < A[i].Chave)
//  { printf("Erro: ChaveNova menor que a chave atual\n");
//    return;
//  }
//  A[i].Chave = ChaveNova;
//  while (i > 1 && A[i / 2].Chave < A[i].Chave)
//    { x = A[i / 2];
//      A[i / 2] = A[i];
//      A[i] = x;
//      i /= 2;
//    }
//}
//
//void Insere(Item *x, Item *A, Indice *n)
//{(*n)++;
//  A[*n] = *x;
//  A[*n].Chave = INT_MIN;
//  AumentaChave(*n, x->Chave, A);
//}
//
//void Imprime(Item *V, Indice *n)
//{ for (i = 1; i <= *n; i++)
//    printf("%d ", V[i].Chave);
//  putchar('\n');
//}
//
//
//int main(int argc, char *argv[])
//{ Item TEMP;
//  n = 7;
//  for (i = 1; i <= n; i++)
//    scanf("%d", &A[i].Chave);
//  /* Teste: 20 15 18 10 12 9 13 */
//  printf("Desordenado: ");
//  Imprime(A, &n);
//
//  printf("Constroi   : ");
//  Constroi(A, &n);
//  Imprime(A, &n);
//
//  printf("Aumenta chave posicao 6 para 25: ");
//  AumentaChave(6, 25, A);
//  Imprime(A, &n);
//
//  x.Chave = 13;
//  printf("Insere%3d: ", x.Chave);
//  Insere(&x, A, &n);
//  Imprime(A, &n);
//
//  TEMP = Max(A);
//  //printf("Max:%3ld\n", TEMP.Chave);
//
//  x = RetiraMax(A, &n);
//  printf("Retira%3d: ", x.Chave);
//  Imprime(A, &n);
//  return(0);
//}
