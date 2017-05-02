/*
 *
 * Originally Tony's code.
 *
 * S.H.Lee modified it.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#undef  uchar
#define uchar   unsigned char
#undef  ushort
#define ushort  unsigned short
#undef  uint
#define uint    unsigned int
#undef  ulong
#define ulong   unsigned long


/* vector gives a one dimensional array, matrix a two dimensional array
   and tensorx a x dimensional array. */

void *Svector(size_t nitem, size_t size)
{
  int	nbyte = nitem * size;
  char *address;

  if(nbyte == 0) nbyte = 1;

  address = malloc(nbyte);

  if(address == NULL) {
    printf("Svector(%d, %d)", (int)nitem, (int)size);
    exit(1) ;
  }

  return(address);
}

/* this routine relies on sizeof(char) == 1 */
void **Smatrix(size_t nitem0, size_t nitem1, size_t size)
{
  size_t tabsize, nbyte;
  void **array0;
  char  *array1;
  int    i;

  /* raw table size for pointers */
  tabsize = nitem0 * sizeof(void*);

  /* align to nearest item boundary */
  if(tabsize % size != 0)
    tabsize += size - tabsize % size;
  
  /* compute the total space needed */
  nbyte = tabsize + nitem0 * nitem1 * size;
  if(nbyte == 0) nbyte = 1;

  /* grab all the space */
  array0 = malloc(nbyte);

  /* exit if problems */
  if(array0 == NULL) {
    printf("Smatrix(%d, %d, %d)", (int)nitem0, (int)nitem1, (int)size);
    exit(1) ;
  }

  /* compute data area */
  array1 = (char*) array0 + tabsize;

  /* make table entries point into data */
  for(i = 0; i < nitem0; i++)
    array0[i] = array1 + i * nitem1 * size;

  return(array0);
}

/* cubic */

void ***Scubic(size_t nitem0, size_t nitem1, size_t nitem2, size_t size)
{
  size_t  nbyte1   ;
  size_t  nbyte    ;
  size_t  tabsize  ;
  size_t  tabsize1 ;
  void ***array0   ;
  char   *array1   ;
  char   *array2   ;
  int     i ,j     ;

  tabsize1 = nitem1 * sizeof(void*);

  if(tabsize1 % size != 0)
    tabsize1 += size - tabsize1 % size ;
  
  nbyte1 = tabsize1 + nitem1 * nitem2 * size;
  if(nbyte1 == 0) nbyte1 = 1;

  tabsize  = nitem0 * sizeof(void **) ; 

  if(tabsize % size != 0)
    tabsize += size - tabsize % size ;

  nbyte = tabsize + nbyte1 * nitem0 ; 
  if(nbyte == 0) nbyte = 1;

  array0 = (void ***) malloc(nbyte);

  /* exit if problems */
  if(array0 == NULL) {
    printf("Scubic(%d, %d, %d, %d)", (int)nitem0, (int)nitem1, (int)nitem2, (int)size);
    exit(1) ;
  }

  /* compute data area */
  array1 = (char*) array0 + tabsize;

  for(i = 0; i < nitem0; i++)
    array0[i] = (void **) (array1 + i * nbyte1) ;

  for (i = 0 ; i < nitem0 ; i++) {
    array2 = array1 + i * nbyte1 + tabsize1 ;

    for (j = 0 ; j < nitem1 ; j++) {
     *((void **)(array1 + i * nbyte1 + j * sizeof(void *))) = 
       array2 + j * nitem2 * size ;
    }
  }

  return(array0);
}


/* Here are some useful derived functions */
char* SvectorChar(size_t nitem)
{
  return( (char *) Svector(nitem, sizeof(char)));
}

uchar* SvectorUchar(size_t nitem)
{
  return( (uchar *) Svector(nitem, sizeof(unsigned char)));
}

short* SvectorShort(size_t nitem)
{
  return( (short *) Svector(nitem, sizeof(short)));
}

ushort* SvectorUshort(size_t nitem)
{
  return( (ushort *) Svector(nitem, sizeof(ushort)));
}

int* SvectorInt(size_t nitem)
{
  return( (int *) Svector(nitem, sizeof(int)));
}

uint* SvectorUint(size_t nitem)
{
  return( (uint *) Svector(nitem, sizeof(uint)));
}

float* SvectorFloat(size_t nitem)
{
  return( (float *) Svector(nitem, sizeof(float)));
}

double* SvectorDouble(size_t nitem)
{
  return( (double *) Svector(nitem, sizeof(double)));
}

char** SmatrixChar(size_t nitem0, size_t nitem1)
{
  return((char**) Smatrix(nitem0, nitem1, sizeof(char)));
}

uchar** SmatrixUchar(size_t nitem0, size_t nitem1)
{
  return((uchar**) Smatrix(nitem0, nitem1, sizeof(uchar)));
}

short** SmatrixShort(size_t nitem0, size_t nitem1)
{
  return((short**) Smatrix(nitem0, nitem1, sizeof(short)));
}

ushort** SmatrixUshort(size_t nitem0, size_t nitem1)
{
  return((ushort**) Smatrix(nitem0, nitem1, sizeof(ushort)));
}

int** SmatrixInt(size_t nitem0, size_t nitem1)
{
  return((int**) Smatrix(nitem0, nitem1, sizeof(int)));
}

uint** SmatrixUint(size_t nitem0, size_t nitem1)
{
  return((uint**) Smatrix(nitem0, nitem1, sizeof(uint)));
}

float** SmatrixFloat(size_t nitem0, size_t nitem1)
{
  return((float**) Smatrix(nitem0, nitem1, sizeof(float)));
}

double** SmatrixDouble(size_t nitem0, size_t nitem1)
{
  return((double**) Smatrix(nitem0, nitem1, sizeof(double)));
}

char*** ScubicChar(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((char***) Scubic(nitem0, nitem1, nitem2, sizeof(char)));
}

uchar*** ScubicUchar(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((uchar***) Scubic(nitem0, nitem1, nitem2, sizeof(uchar)));
}

short*** ScubicShort(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((short***) Scubic(nitem0, nitem1, nitem2, sizeof(short)));
}

ushort*** ScubicUshort(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((ushort***) Scubic(nitem0, nitem1, nitem2, sizeof(ushort)));
}

int*** ScubicInt(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((int***) Scubic(nitem0, nitem1, nitem2, sizeof(int)));
}

uint*** ScubicUint(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((uint***) Scubic(nitem0, nitem1, nitem2, sizeof(uint)));
}

float*** ScubicFloat(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((float***) Scubic(nitem0, nitem1, nitem2, sizeof(float)));
}

double*** ScubicDouble(size_t nitem0, size_t nitem1, size_t nitem2)
{
  return((double***) Scubic(nitem0, nitem1, nitem2, sizeof(double)));
}


void P2free (char **dptr)
{
	if(dptr==NULL) return;
	free((char *)dptr);
}


void	**alloc2d(size_t dim1, size_t dim2, size_t size)
{
	int		i;		
	size_t	nelem;	

	char	*p;	
	void	**pp;

	nelem = dim1 * dim2;

	p = (char *)calloc(nelem, size);

	pp = (void **) calloc(dim1, sizeof(char *));

	for (i = 0; i < dim1; i++)
		pp[i] = p + (i * dim2 * size);

	return (pp);
}


int 	free2d(void **p)
{

	if ( p != NULL && *p != NULL) 
		free( (void *) *p);
	if ( p != NULL)
		free( (void *) p);
	return(1);
}

