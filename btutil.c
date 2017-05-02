/******************************************************
 *                                                    *
 * btutil.c : Bradley-Terry Model utility.c           *
 *                                                    *
 * Sangho Lee                                         *
 *                                                    *
 ******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include "Smemory.h"

#include "pair_matrix.h"

void sort_matrix(double **in_matrix, int height, int width) ;
int GetPairMatrixSize(double **in_matrix, int height, int width) ;
int MakePairMatrix(pair_matrix_t **pairMatrix, double **matrix, int nRecord, int nCol);
int *GetIdxToCases(double *target, int nRecord1, int nUniq);
double *GetUniqTargets(double **matrix, double *target, int nRecord1, int nFeat1, int *nUniq);
int *GetIdxToCases(double *target, int nRecord1, int nUniq);
//////////////////////////////////////////////////////////////

// sort matrix by the first column

void sort_matrix(double **in_matrix, int height, int width) 
{
	int      i , j , idx ; 
	double  **workmatrix ;
    double  *workvector  ; // for sorting target vector
	int     *workvectIdx ;
	extern void merge_sortidx() ;

	workmatrix  = (double **) SmatrixDouble(height, width)    ;
	workvectIdx = (int *)     malloc(sizeof(int) * height)    ;
	workvector  = (double *)  malloc(sizeof(double) * height) ;

	// copy the matrix to workmatrix
	for (i = 0 ; i < height ; ++i) {
	    for (j = 0 ; j < width ; ++j)
		workmatrix[i][j] = in_matrix[i][j]  ;
	    workvector[i]    = workmatrix[i][0] ; // List Set Index
	    workvectIdx[i]   = i ;
	}

	// quick sorting workvectIdx
    merge_sortidx(workvector, workvectIdx, 0, height - 1) ;

	// copy workmatrix to in_matrix
	for (i = 0 ; i < height ; ++i) {
	    idx = workvectIdx[i] ;
	    for (j = 0 ; j < width ; ++j)
            in_matrix[i][j] = workmatrix[idx][j] ; 
	}

	// free memory
	free(workmatrix)  ;
	free(workvectIdx) ;
	free(workvector)  ;

#ifdef DEBUG
    for (i = 0 ; i < height ; ++i)  {
        for (j = 0 ; j < width ; ++j)
            printf("%f ", in_matrix[i][j]) ;
        printf("\n") ;
    }
#endif

}

//////////////////////////////////////////////////////////////

int GetPairMatrixSize(double **in_matrix, int height, int width) 
{
    int matrixSize = 0, listSize ; 
    int idx ;
    int prevIdx ;

    prevIdx = (int) in_matrix[0][0] ;

    listSize = 1 ; 
    for (idx = 1 ; idx < height ; ++idx) {
        if ((int) in_matrix[idx][0] == prevIdx) ++listSize ;
        else {
            matrixSize += (listSize-1) * listSize ;
            listSize = 1 ; 
            prevIdx = (int) in_matrix[idx][0] ;
        }
    }
    matrixSize += (listSize-1) * listSize ;

    return matrixSize ;
}

//////////////////////////////////////////////////////////////


int MakePairMatrix(pair_matrix_t **pairMatrix, double **matrix, int nRecord, int nCol)
{
    double *workVector ;
    int    *workVectorIdx, *workIdx ;
	int     numunique, nPairs ;
	int     i , idx, tmpIdx, uniqIdx ; 
    int     beg, end, leftidx, rightidx ;
    int     begleft, endleft;
    int     realIdx1, realIdx2 ;
    double  prevValue, prevValue2 ;
	extern void merge_sortidx() ;
#ifdef DEBUG
	int		k;
#endif


    nPairs = 0 ; 
    beg = 0 ; prevValue = matrix[beg][0] ;
    while (beg < nRecord) {

        // find the list of the same Index
        for (end = beg ; end < nRecord && matrix[end][0] == prevValue ; ++end) ;
        
        // beg .. end - 1
        workVectorIdx = SvectorInt   (end - beg)     ;
        workVector    = SvectorDouble(end - beg)     ;
        workIdx       = SvectorInt   (end - beg + 1) ;
	    for (i = 0 ; i < end - beg ; ++i) {
	        workVector[i]    = matrix[beg+i][nCol-1] ; // target value
	        workVectorIdx[i] = i                     ;
        }

        // Sort the same list
        merge_sortidx(workVector, workVectorIdx, 0, end - beg - 1) ;

        // Adjust the Index
        for (i = 0 ; i < end - beg ; ++i)
            workVectorIdx[i] += beg ;

        // Compute the number of unique targets
        prevValue2  = matrix[workVectorIdx[0]][nCol-1] ;
        numunique  = 0 ;
        workIdx[0] = 0 ;
	    for (idx = 1 ; idx < end - beg ; ++idx) {
            tmpIdx = workVectorIdx[idx] ;
            if (matrix[tmpIdx][nCol-1] != prevValue2) {
               prevValue2 = matrix[tmpIdx][nCol-1] ;
               workIdx[++numunique] = idx ;
            } 
        }
        workIdx[++numunique] = idx ;

        // Write the result
        for (uniqIdx = 0 ; uniqIdx < numunique - 1 ; ++uniqIdx) {
            begleft = workIdx[uniqIdx]     ; 
            endleft = workIdx[uniqIdx + 1] ;

            for (leftidx = begleft ; leftidx < endleft ; ++leftidx) {
                realIdx1 = workVectorIdx[leftidx] ;

                for (rightidx = endleft ; rightidx < end - beg ; ++rightidx) {
                    realIdx2 = workVectorIdx[rightidx] ; 
					
					pairMatrix[ nPairs ]->x = matrix[realIdx2];
					pairMatrix[ nPairs ]->y = matrix[realIdx1];

					//pairMatrix[ nPairs + 1 ]->y = matrix[realIdx1];
					//pairMatrix[ nPairs + 1 ]->x = matrix[realIdx2];


                    //pairTarget[nPairs]   = 0.0 ; // lose
                    //pairTarget[nPairs+1] = 1.0 ; // win
                    //nPairs += 2 ;
					nPairs++;
                } // rightidx
            } // leftidx
        } // uniqIdx

        beg = end ; 
        if (beg < nRecord) prevValue = matrix[beg][0] ;

        free(workVectorIdx) ;
        free(workVector)    ;
        free(workIdx)       ;
    }

#ifdef DEBUG
    for (i = 0 ; i < nPairs ; ++i) {
        for (k = 1 ; k < nCol -1  ; ++k){
			//printf("%lf\t", pairMatrix[ i ]->x[ k ] - pairMatrix[ i ]->y[k]);
			printf("%lf\t", pairMatrix[ i ]->x[ k ] );
		}
		printf("%d\n", 1);
        for (k = 1 ; k < nCol -1  ; ++k){
			//printf("%lf\t", pairMatrix[ i ]->y[ k ] - pairMatrix[ i ]->x[k]);
			printf("%lf\t", pairMatrix[ i ]->y[ k ]);
		}
		printf("%d\n", -1);
    }
#endif
   	fprintf(stderr,"Num of compare pair %d\n", nPairs) ;

    return nPairs ;

}

//////////////////////////////////////////////////////////////

double *GetUniqTargets(double **matrix, double *target, int nRecord1, int nFeat1, int *nUniq)
{
    double **workMatrix, *workVector ;
    double *tmpVector, *tmpVector2 ;
    int    *workVectorIdx ;
	int     numunique ;
	int     i , j , idx ; 
	extern void merge_sortidx() ;

    workMatrix    = SmatrixDouble(nRecord1, nFeat1) ;
    workVectorIdx = SvectorInt   (nRecord1) ;
    workVector    = SvectorDouble(nRecord1) ;
    tmpVector     = SvectorDouble(nRecord1) ; // 0 .. nRecord + 1 

	// copy the matrix to workmatrix 
	for (i = 0 ; i < nRecord1 ; ++i) {
	    for (j = 0 ; j < nFeat1 ; ++j)
		    workMatrix[i][j] = matrix[i][j] ;
	    workVector[i]    = target[i]    ;
	    workVectorIdx[i] = i            ;
	}

	// quick sorting workvectIdx
    // workVectorIdx holds the sorted Index to workVector.
    merge_sortidx(workVector, workVectorIdx, 1, nRecord1 - 1) ;

	// copy workmatrix to in_matrix
	for (i = 1 ; i < nRecord1 ; ++i) {
	    idx = workVectorIdx[i] ;
	    for (j = 0 ; j < nFeat1 ; ++j) 
            matrix[i][j] = workMatrix[idx][j] ; 
        tmpVector[i-1] = target[i] = workVector[idx] ;
	}

    // tmpVector : 0 .. nRecord1-2 
	numunique = 0 ; 
	for (idx = 1 ; idx < nRecord1 - 1 ; ++idx) 
	    if (tmpVector[numunique] != tmpVector[idx]) 
		    tmpVector[++numunique] = tmpVector[idx] ;

    tmpVector2 = SvectorDouble(numunique+1) ; 
    for (idx = 0 ; idx < numunique+1 ; ++idx)
        tmpVector2[idx] = tmpVector[idx] ;

    free(tmpVector)     ;
    free(workMatrix)    ;
    free(workVectorIdx) ;
    free(workVector)    ;

    printf("# of unique target levels : %d\n", numunique+1) ;
    for (idx = 0 ; idx < numunique+1 ; ++idx)
        printf("%d : %lf\n",idx, tmpVector2[idx]) ; 

    *nUniq = numunique+1 ;
    return tmpVector2 ;
}

//////////////////////////////////////////////////////////////

int *GetIdxToCases(double *target, int nRecord1, int nUniq)
{
    int idx, prevIdx, prevValue ;
    int *tmpVector ;

    tmpVector = SvectorInt(nUniq+1) ;
    tmpVector[0] = 1 ; // tUniq[0] --> target[1]
    prevIdx   = 0 ;
    prevValue = target[ tmpVector[prevIdx] ] ; 

    for (idx = 2 ; idx < nRecord1 ; ++idx) {
        if (target[idx] != prevValue) {
           prevIdx++ ; 
           tmpVector[prevIdx] = idx ;
           prevValue = target[idx] ;
        }
    }
    tmpVector[++prevIdx] = idx ; // Idx to final data 

    if (prevIdx != nUniq+1) {
        fprintf(stderr,"Error in GetIdxToCases : %d, %d\n", prevIdx, nUniq+1) ;
        exit(1) ;
    }

    return tmpVector ;
}

//////////////////////////////////////////////////////////////

