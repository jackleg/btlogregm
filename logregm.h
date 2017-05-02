#ifndef __LOGREGM_H__
#define __LOGREGM_H__

extern int LearnLogisticRegression(int ngroups, pair_matrix_t **matrix, int nfeats, double *beta, double *se_beta, int nFlag);
extern void EstiLogisticRegression(int ngroups, double **matrix, int nfeats, double *estimated, double *beta);
//extern void EstiBTStrength(int ngroups, double **matrix, int nfeats, double *estimated, double *beta, int noexp);
extern void EstiBTStrength(int ngroups, double **matrix, int nfeats, double *beta, int noexp, char *outputfile);


enum {
	TRAIN_CG = 0,
	TRAIN_BFGS
};


#endif

