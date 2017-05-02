/******************************************************************************
*                                                                             *
*       Copyright (C) 1993,1994 Tony Robinson				      *
*                                                                             *
*       See the file SLICENSE for conditions on distribution and usage	      *
*                                                                             *
******************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>


#undef  uchar
#define uchar   unsigned char
#undef  ushort
#define ushort  unsigned short
#undef  uint
#define uint    unsigned int
#undef  ulong
#define ulong   unsigned long

#ifdef __STDC__
extern void *Svector(size_t nitem, size_t size)  ; 

/* this routine relies on sizeof(char) == 1 */

extern void **Smatrix(size_t nitem0, size_t nitem1, size_t size)  ; 
extern void **alloc2d(size_t dim1, size_t dim2, size_t size);
extern int 	free2d(void **p);
extern void P2free(void **ptr);

/* Here are some useful derived functions */

extern char* SvectorChar(size_t nitem)  ;

extern uchar* SvectorUchar(size_t nitem) ;

extern short* SvectorShort(size_t nitem) ;

extern ushort* SvectorUshort(size_t nitem)  ;

extern int* SvectorInt(size_t nitem) ;

extern uint* SvectorUint(size_t nitem) ;

extern float* SvectorFloat(size_t nitem) ;

extern double* SvectorDouble(size_t nitem) ;

extern char** SmatrixChar(size_t nitem0, size_t nitem1) ;

extern uchar** SmatrixUchar(size_t nitem0, size_t nitem1) ;

extern short** SmatrixShort(size_t nitem0, size_t nitem1) ;

extern ushort** SmatrixUshort(size_t nitem0, size_t nitem1) ;

extern int** SmatrixInt(size_t nitem0, size_t nitem1) ;

extern uint** SmatrixUint(size_t nitem0, size_t nitem1) ;

extern float** SmatrixFloat(size_t nitem0, size_t nitem1) ;

extern double** SmatrixDouble(size_t nitem0, size_t nitem1) ;

extern char*** ScubicChar(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern uchar*** ScubicUchar(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern short*** ScubicShort(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern ushort*** ScubicUshort(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern int*** ScubicInt(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern uint*** ScubicUint(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern float*** ScubicFloat(size_t nitem0, size_t nitem1, size_t nitem2) ;

extern double*** ScubicDouble(size_t nitem0, size_t nitem1, size_t nitem2) ;

#else

extern void *Svector()  ; 

/* this routine relies on sizeof(char) == 1 */

extern void **Smatrix()  ; 

/* Here are some useful derived functions */

extern char* SvectorChar()  ;

extern uchar* SvectorUchar() ;

extern short* SvectorShort() ;

extern ushort* SvectorUshort()  ;

extern int* SvectorInt() ;

extern uint* SvectorUint() ;

extern float* SvectorFloat() ;

extern double* SvectorDouble() ;

extern char** SmatrixChar() ;

extern uchar** SmatrixUchar() ;

extern short** SmatrixShort() ;

extern ushort** SmatrixUshort() ;

extern int** SmatrixInt() ;

extern uint** SmatrixUint() ;

extern float** SmatrixFloat() ;

extern double** SmatrixDouble() ;

extern char*** ScubicChar() ;

extern uchar*** ScubicUchar() ;

extern short*** ScubicShort() ;

extern ushort*** ScubicUshort() ;

extern int*** ScubicInt() ;

extern uint*** ScubicUint() ;

extern float*** ScubicFloat() ;

extern double*** ScubicDouble() ;


extern void **alloc2d();
extern int 	free2d();
extern void P2free();


#endif

