#ifndef __BTUTIL_H__
#define __BTUIIL_H__

extern void sort_matrix(double **in_matrix, int height, int width) ;
extern int GetPairMatrixSize(double **in_matrix, int height, int width) ;
extern int MakePairMatrix(pair_matrix_t **pairMatrix, double **matrix, int nRecord, int nCol);
extern int *GetIdxToCases(double *target, int nRecord1, int nUniq);
extern double *GetUniqTargets(double **matrix, double *target, int nRecord1, int nFeat1, int *nUniq);
extern int *GetIdxToCases(double *target, int nRecord1, int nUniq);


#endif
