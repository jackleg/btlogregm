/*****************************************************
 *                                                    *
 * linmdrvr.c : linear machine driver.c               *
 *                                                    *
 * Sangho Lee                                         *
 *                                                    *
 * Korea Advanced Institute of Science and Technology *
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
#include "logregm.h"
#include "btutil.h"

#ifndef __STDC__
#define __DATE__ "1996. 12. 25"
#define __TIME__ "00:00:00"
#endif

//extern int Sgetargs  () ;  
extern int Sgetargs(int argc, char **argv,...);
extern void sort_matrix(double **in_matrix, int height, int width) ;

//////////////////////////////////////////////////////////////////////

void BT_ReadMatrix(double **matrix, double *target, 
                char inputfile[], int nRecord, int nFeats)
{
    FILE *fptr ; 
    int idx, idx2 ;

    if ((fptr = fopen(inputfile,"rt")) == NULL) {
        fprintf(stderr,"Error in opening %s\n", inputfile) ;
        exit(1) ;
    }

    for (idx = 0 ; idx < nRecord ; ++idx) {
        for (idx2 = 0 ; idx2 < nFeats ; ++idx2)
            fscanf(fptr,"%lf", &matrix[idx][idx2]) ;
        fscanf(fptr,"%lf", &target[idx]) ;
        matrix[idx][idx2] = target[idx] ;
    }

    fclose(fptr) ;
}

////////////////////////////////////////////////////////////////////////

void WriteBeta(double beta[], int nFeats, char parafile[])
{
    FILE *fptr ;
    int idx ;

    if ((fptr = fopen(parafile,"wt")) == NULL) {
        fprintf(stderr,"Can't write %s\n", parafile) ;
        exit(1) ;
    }

    for (idx = 0 ; idx < nFeats+1 ; ++idx) 
        fprintf(fptr,"%f\n", beta[idx]) ;

    fclose(fptr) ;
}
    
////////////////////////////////////////////////////////////////////////

void WriteBTBeta(double beta[], int nFeats, char parafile[])
{
    FILE *fptr ;
    int idx ;

    if ((fptr = fopen(parafile,"wt")) == NULL) {
        fprintf(stderr,"Can't write %s\n", parafile) ;
        exit(1) ;
    }

    for (idx = 0 ; idx < nFeats ; ++idx) 
        fprintf(fptr,"%.10lf\n", beta[idx]) ;

    fclose(fptr) ;
}
    
////////////////////////////////////////////////////////////////////////

void ReadBeta(double beta[], int nFeats, char parafile[])
{
    FILE *fptr ;
    int idx ;

    if ((fptr = fopen(parafile,"rt")) == NULL) {
        fprintf(stderr,"Error in opening %s\n", parafile) ;
        exit(1) ;
    }

    for (idx = 0 ; idx < nFeats ; ++idx) 
        fscanf(fptr,"%lf", &beta[idx]) ;

    fclose(fptr) ;
}

////////////////////////////////////////////////////////////////////////

void DumpBetaScreen(double *beta, double *se_beta, int nFeats)
{
    int idx ;

    printf("----------------------------------------------\n") ;
    printf("                Coeff.     Std.Err     Z score\n") ;
    printf("----------------------------------------------\n") ;
    for (idx = 1 ; idx < nFeats ; ++idx) {
		if( se_beta[idx] >= (1.0e-10) ){
		   	printf("   %-2d var.  %10f  %10f  %10f\n", idx, beta[idx], se_beta[idx], beta[idx]/se_beta[idx]) ;
		}
		else {
		   	printf("   %-2d var.  %10f  %10f    Inf\n", idx, beta[idx], se_beta[idx]) ;
		}
               //beta[idx], se_beta[idx], beta[idx]) ;
    }
    printf("----------------------------------------------\n") ;

}

////////////////////////////////////////////////////////////////////////

void WriteResult(char outputfile[], double *target, double *esti, int nRecord)
{
    FILE *fptr ;
    int idx ;

    if ((fptr = fopen(outputfile,"wt")) == NULL) {
        fprintf(stderr,"Error in writing %s\n", outputfile) ;
        exit(1) ;
    }

    for (idx = 1 ; idx < nRecord+1 ; ++idx) 
        fprintf(fptr,"%f %f\n", target[idx], esti[idx]) ;

    fclose(fptr) ;
}

////////////////////////////////////////////////////////////////////////

void WriteBTResult(char outputfile[], double **matrix, int nCol, 
            double *esti, int nRecord)
{
    FILE *fptr ;
    int idx ;

    if ((fptr = fopen(outputfile,"wt")) == NULL) {
        fprintf(stderr,"Error in writing %s\n", outputfile) ;
        exit(1) ;
    }

    for (idx = 0 ; idx < nRecord ; ++idx) 
        fprintf(fptr,"%d %f %f\n", (int) matrix[idx][0], matrix[idx][nCol-1], esti[idx]) ;

    fclose(fptr) ;
}

////////////////////////////////////////////////////////////////////////

void cut_tail(char *str)
{
	int idx ; 
	for (idx = strlen(str) - 1 ; idx >= 0 ; --idx) {
		if (str[idx] == '\n' || str[idx] == '\r' || 
                str[idx] == ' ') str[idx] = '\0' ;
		else break ;
	}
}

////////////////////////////////////////////////////////////////////////

//
// Given the name of a file of numeric data, this routine counts
// the numbers of rows and columns. It's assumed that the number
// of entries is the same in each row, and an error is flagged if this
// is not the case.
//

void file_size (char *file_name, int *width, int *height)
{
    FILE *f;
    int buf_size = 600, _width = 0;
    char *buffer, *ptr;

    *width = *height = 0;

    buffer = (char*) malloc (buf_size * sizeof (char));
    if (buffer == NULL) {
        fprintf(stderr,"Malloc Error %s %d\n",__FILE__, __LINE__) ;
        exit(1) ;
    }

    /* Open price file - abort if filename invalid */
    f = fopen(file_name, "r");
    if (f == NULL) {
        fprintf(stderr,"\n File not found : %s\n", file_name) ;
        fprintf(stderr,"Error : %s %d\n", __FILE__, __LINE__) ;
        exit(1) ;
    }

    /* Get number of entries in first row */
    if (fgets(buffer, buf_size, f) != NULL) {
        cut_tail(buffer) ;

        ++*height;
        ptr = strtok (buffer, "\t ,");
        while (ptr != NULL) {
            ++*width;
            ptr = strtok (NULL, "\t ,");
        }
    }

    /* Count numbers of subsequent rows */
    while (!feof(f)) {
        if (fgets(buffer, buf_size, f) != NULL) {
            cut_tail(buffer) ;
            if (strlen(buffer) > strlen("\n")) { /* if line is more than a NL char */
                if (buffer[0] == '#')
                    continue;
            
                ++*height;
                _width = 0;
                ptr = strtok (buffer, "\t ,");
                while (ptr != NULL) {
                    ++_width;
                    ptr = strtok (NULL, "\t ,");
                }

                if (*width != _width) {
                    printf("\n Number of entries in file %s did not agree", file_name);
                    fprintf(stderr,"\n Error : %s %d\n", __FILE__, __LINE__) ;
                    exit(1) ;
                }
            }
        }
    }
    fclose(f) ; 
    free (buffer);
}

////////////////////////////////////////////////////////////////////////////

void HowTo()
{
    printf("HowTo Page.\n") ;
    printf("the first column should have an index to list.\n") ;
    printf("the first feature begins at the second column.\n") ;
}

////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[])
{
    double **matrix, *target ;
    double *beta, *se_beta ;
    //double *tUniq, **pairMatrix, *pairTarget ;
    //double *pairTarget ;
    int     approxSize, howto = 0 ;
    int     pairSize, noexp=1 ;
    char   *parafile   = "beta.btm"  ;
    //char   *inputfile  = "datafname" ;
	char	*inputfile="list_work/tmp.dat";
    char   *outputfile = "outfile"   ;
    int     nRecord = 1 ;
    int     nFeats  = 1 ;
    int     n_vars, n_samples, do_sort = 0 ;

    int     trainflag  = 1            ; 


	pair_matrix_t	**pairMatrix;



    Sgetargs(argc, argv, "-h", "-help", NULL,
		"-train",   NULL, &trainflag,   "if set, train a bradley-terry logistic reg, if not test it",
	    "-mdl",     "%s", &parafile ,   "output bradley-terry LogReg file name (human-readable)",
	    "-input",   "%s", &inputfile,   "input data file name",
        "-output",  "%s", &outputfile,  "output file name, when testing on input file",
        "-how",     NULL, &howto,       "How-To page",
        "-noexp",   NULL, &noexp,       "output weighted sum only",
        "-sort",    NULL, &do_sort,     "do sort the data by the first column",
	         NULL,
	         NULL);

    if (howto) { HowTo() ; exit(1) ; } ;


#if 1
    // Read File Size
    file_size (inputfile, &n_vars, &n_samples) ;


    nRecord = n_samples  ;
    nFeats  = n_vars - 1 ;

    fprintf(stderr,"# of records  : %d\n", nRecord)  ;
    fprintf(stderr,"# of features : %d\n", nFeats+1) ;
    fprintf(stderr,"Program runs on %s mode.\n", (trainflag) ? "TRAIN" : "TEST") ;

    matrix  = SmatrixDouble(nRecord, nFeats+1) ; // 0 .. nFeats
    target  = SvectorDouble(nRecord)           ;
    //esti    = SvectorDouble(nRecord)           ; // 1 .. pairSize
    beta    = SvectorDouble(nFeats+1)          ;
    se_beta = SvectorDouble(nFeats+1)          ;

    BT_ReadMatrix(matrix, target, inputfile, nRecord, nFeats) ;
    if (do_sort) sort_matrix(matrix, nRecord, nFeats+1) ;

    if (trainflag==1) {

        approxSize = GetPairMatrixSize(matrix, nRecord, nFeats+1) ;

        // now make a pair matrix
        pairMatrix = (pair_matrix_t **)Smatrix( approxSize + 1, 1 , sizeof( pair_matrix_t));
        pairSize   = MakePairMatrix(pairMatrix, matrix, nRecord, nFeats+1) ;

        // 1 .. pairSize-1 , 1 .. nFeats - 1
        LearnLogisticRegression(pairSize, pairMatrix, nFeats - 1, beta, se_beta, TRAIN_CG) ;
        WriteBTBeta(beta, nFeats, parafile) ; // 1 .. nFeats - 1
        EstiBTStrength(nRecord, matrix, nFeats, beta, noexp, outputfile) ; // 1 .. nFeats -1
        //EstiBTStrength(nRecord, matrix, nFeats, esti, beta, noexp) ; // 1 .. nFeats -1
        //WriteBTResult(outputfile, matrix, nFeats+1, esti, nRecord) ;
        DumpBetaScreen(beta, se_beta, nFeats) ;

		free(pairMatrix);
        //free(pairTarget) ;

    } else {
        ReadBeta(beta, nFeats, parafile) ;
        //EstiBTStrength(nRecord, matrix, nFeats, esti, beta, noexp) ;
        EstiBTStrength(nRecord, matrix, nFeats, beta, noexp, outputfile) ; // 1 .. nFeats -1
        //WriteBTResult(outputfile, matrix, nFeats+1, esti, nRecord) ;
    }
#endif

    free(matrix)  ;
    free(target)  ;
    //free(esti)    ;
    free(beta)    ;
    free(se_beta) ;

    return 0 ;
}

/*-------------------------------------------------------------------*/



