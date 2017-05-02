#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Smemory.h"
#include "pair_matrix.h"
#include "logregm.h"


#define NUM_ITER	100
#define	EPSILON	1.0e-6


static int train_cg(int nSample, pair_matrix_t **pair_data, int nDim, double *wgt, double *se_beta);
static int train_bfgs(int nSample, pair_matrix_t **pair_data, int nDim, double *wgt, double *se_beta);

int LearnLogisticRegression(int nSample, pair_matrix_t **pair_data, int nDim, double *wgt, double *se_beta, int nFlag)
{

	if( nFlag == TRAIN_CG ){
	   	train_cg(nSample, pair_data, nDim, wgt, se_beta);
	}
	else if ( nFlag == TRAIN_BFGS ){
	   	train_bfgs(nSample, pair_data, nDim, wgt, se_beta);
	}
	else{
		fprintf(stderr,"Wrong Training method\n");
		exit(1);
	}

	return 0;
}



static int train_bfgs(int nSample, pair_matrix_t **pair_data, int nDim, double *wgt, double *se_beta)
{
	int	iter;
	int	num_iter;
	int	idx_n, idx_d;

	double	*grad;
	double	*old_grad;
	double	*old_wgt;
	double	*prev_wgt;
	double	*dir_u;
	double	*dw;
	double	*dg;
	double	dwdg;

	double	sum;
	double	sigmoid;
	double	lambda;

	double	ag;
	double	aw;
	double	b;

	double	ug;
	double	ux;
	double	uhu;
	double	dev;
	double	old_dev;
	double	tmp_f;

	double	*delta;
	double	*data_x;
	double	*data_y;
	double	converge;

	double	*wgt_sum;
	double	*wgt_sq_sum;



	delta = SvectorDouble( nDim );
	grad = SvectorDouble( nDim );
	old_grad = SvectorDouble( nDim );
	old_wgt = SvectorDouble( nDim + 1);
	prev_wgt = SvectorDouble( nDim + 1);
	dw = SvectorDouble( nDim + 1);
	dg = SvectorDouble( nDim );
	dir_u = SvectorDouble( nDim );

	wgt_sum = SvectorDouble( nDim + 1);
	wgt_sq_sum = SvectorDouble( nDim + 1);
	for( idx_d = 0; idx_d <= nDim; idx_d++){
		wgt_sum[ idx_d ] = 0.0;
		wgt_sq_sum[ idx_d ] = 0.0;
	}

	num_iter = NUM_ITER;
	old_dev = -100000.0;
	lambda = 0.2;




	for( iter = 1; iter < num_iter ; iter++){

		
		for( idx_d = 0; idx_d < nDim; idx_d++){
			grad[ idx_d ] = 0.0;
			old_wgt[ idx_d + 1] = wgt[ idx_d + 1];
		}


		for( idx_n = 0; idx_n < nSample; idx_n++){

			data_x = pair_data[ idx_n ]->x + 1;
			data_y = pair_data[ idx_n ]->y + 1;

		   	sum = 0.0;
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
			   	sum += (wgt[ idx_d + 1] * ( data_x[ idx_d ] - data_y[idx_d]));
			}
			sigmoid = 1.0 / ( 1 + exp( sum ));
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
				grad[ idx_d ] += (sigmoid * (data_x[ idx_d ] - data_y[ idx_d  ]));
			}

		}

		for( idx_d = 0; idx_d < nDim; idx_d++){
			grad[ idx_d ] = 2.0*grad[ idx_d ] - lambda * wgt[ idx_d + 1];
		}
		



		if( iter == 1 ){
			for( idx_d = 0; idx_d < nDim; idx_d++){
				dir_u[ idx_d ] = -grad[ idx_d ] ;
			}
		}
		else{


			for( idx_d = 0; idx_d < nDim; idx_d++){
				dw[ idx_d + 1] = wgt[ idx_d + 1] - prev_wgt[ idx_d + 1 ];
				dg[ idx_d ] = grad[ idx_d ] - old_grad[ idx_d ];
			}

			dwdg = 0.0;
			
			b = 0.0;
			ag = 0.0;
			aw = 0.0;
			for( idx_d = 0; idx_d < nDim; idx_d++){
				dwdg += (dw[ idx_d + 1] * dg[ idx_d ]);
				b += ( dg[ idx_d ] * dg[ idx_d ] );

				ag += dw[ idx_d + 1] * grad[ idx_d ];
				aw += dg[ idx_d ] * grad[ idx_d];
			}

			b = 1.0 + b / dwdg ;

			ag = ag / dwdg;
			aw = aw / dwdg - b * ag;



			for( idx_d = 0; idx_d < nDim; idx_d++){
				dir_u [ idx_d ] = -grad [ idx_d ]+ aw * dw[ idx_d + 1] + ag * dg[ idx_d ];
			}

		}

		for( idx_d = 0; idx_d < nDim; idx_d++){
			prev_wgt [ idx_d + 1] = wgt[ idx_d + 1];
		}


		// line search along dir_u

		ug = 0.0;
		uhu = 0.0;
		for( idx_d = 0; idx_d < nDim; idx_d++){
			ug += ( dir_u[ idx_d ] * grad[ idx_d] );
			uhu += ( dir_u[ idx_d ] * dir_u[ idx_d ] );
		}
		uhu *= lambda;

		for( idx_n = 0; idx_n < nSample; idx_n++){

			data_x = pair_data[ idx_n ]->x + 1;
			data_y = pair_data[ idx_n ]->y + 1;

		   	ux = 0.0;
			sum = 0.0;
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
				ux += (dir_u[ idx_d ] * (data_x[ idx_d ] - data_y[ idx_d]));
			   	sum += (wgt[ idx_d + 1] * (data_x[ idx_d ] - data_y[idx_d]));
			}
			tmp_f = exp( sum );
			uhu += (( ux * ux ) / ( 2 + tmp_f + 1.0/tmp_f )); 

		}


		tmp_f = ug / (2.0 * uhu) ;
		for( idx_d = 0; idx_d < nDim; idx_d++){
			wgt[ idx_d + 1] = wgt[ idx_d + 1] + tmp_f * dir_u[ idx_d ];
			old_grad[ idx_d ] = grad[ idx_d ];

			//printf("wgt[%d] = %lf %lf\n", idx_d, wgt[ idx_d ], tmp_f);
		}

		for( idx_d = 0; idx_d <= nDim; idx_d++){
			tmp_f = wgt[ idx_d ]  ;
			wgt_sum[ idx_d ] += tmp_f;
			wgt_sq_sum[ idx_d ] += (tmp_f * tmp_f );
		}

#if 0
		// converge check
		dev = 0.0;
		for( idx_d = 0; idx_d <= nDim; idx_d++){
			tmp_f = ( old_wgt[ idx_d ] - wgt[ idx_d ] ) ;
			dev += ( tmp_f * tmp_f);
			se_beta[ idx_d ] = fabs(tmp_f);
		}
		converge = sqrt(dev / nDim);
#endif


#if 1
		dev = 0.0;
		for( idx_n = 0; idx_n < nSample; idx_n++){

			data_x = pair_data[ idx_n ]->x + 1;
			data_y = pair_data[ idx_n ]->y + 1;

			sum = 0.0;
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
			   	sum += (wgt[ idx_d + 1] * (data_x[ idx_d ] - data_y[ idx_d ]));
			}

			dev += log( 1.0 + exp(-sum));

		}

		sum = 0.0;
		for( idx_d = 0; idx_d < nDim; idx_d++){
			sum += (wgt[ idx_d + 1] * wgt[ idx_d + 1]);
		}
		dev = 2.0*dev + 0.5 * lambda * sum ;


		converge = fabs( 1.0 - old_dev / dev );
		//printf("iter[%d]:  %lf %lf %lf\n", iter, dev, old_dev, fabs( dev - old_dev));
#endif
		if( converge < EPSILON ) break;

		if( (iter + 1 ) % 5 == 0 ){
		   	fprintf(stderr, " iter[%d]:\t%lf\t%lf\n", iter, dev, converge);
		}
		/////////////////////


		old_dev = dev;



	}
   
	fprintf(stderr, " iter[%d]:\t%lf\t%lf\n", iter, dev, converge);
#if 0

	for( idx_n = 0; idx_n < nSample; idx_n++){

		data_x = pair_data[ idx_n ]->x + 1;
		data_y = pair_data[ idx_n ]->y + 1;

		for( idx_d = 0; idx_d < nDim; idx_d++){
			tmp_f = (data_x[ idx_d ] - data_y[ idx_d]);
			tmp_f = tmp_f * wgt[idx_d + 1];
			wgt_sum[idx_d + 1] += (tmp_f);
			wgt_sq_sum[idx_d + 1] += (tmp_f * tmp_f);
		}
	}
#endif

	for( idx_d = 0; idx_d <= nDim; idx_d++){
		wgt_sum[ idx_d ] /= (double)iter;
		wgt_sq_sum[ idx_d ] /= (double)iter;

		tmp_f = wgt_sq_sum[idx_d] - (wgt_sum[idx_d] * wgt_sum[idx_d]);

		if( tmp_f > 0.0 ){
		   	se_beta[ idx_d ] = sqrt(tmp_f / iter);
		}
		else 
		   	se_beta[ idx_d ] = 0.0;
	}

	free(delta );
	free(grad );
	free(old_grad );
	free(old_wgt );
	free(dw );
	free(dg );
	free(prev_wgt );
	free(dir_u );

	free( wgt_sum );
	free( wgt_sq_sum );



	return 0;
}



static int train_cg(int nSample, pair_matrix_t **pair_data, int nDim, double *wgt, double *se_beta)
{
	int	iter;
	int	num_iter;
	int	idx_n, idx_d;

	double	*grad;
	double	*old_grad;
	double	*old_wgt;
	double	*dir_u;

	double	sum;
	double	sigmoid;
	double	lambda;

	double	ug;
	double	ux;
	double	uhu;
	double	dev;
	double	old_dev;
	double	tmp_f;

	double	beta;
	double	*delta;
	double	gd;
	double	ud;
	double	*data_x;
	double	*data_y;


	double	*wgt_sum;
	double	*wgt_sq_sum;

	double	converge;



	delta = SvectorDouble( nDim );
	grad = SvectorDouble( nDim );
	old_grad = SvectorDouble( nDim );
	old_wgt = SvectorDouble( nDim + 1);
	dir_u = SvectorDouble( nDim );

	wgt_sum = SvectorDouble( nDim + 1);
	wgt_sq_sum = SvectorDouble( nDim + 1);
	for( idx_d = 0; idx_d <= nDim; idx_d++){
		wgt_sum[ idx_d ] = 0.0;
		wgt_sq_sum[ idx_d ] = 0.0;
	}

	num_iter = NUM_ITER;
	old_dev = -100000.0;
	lambda = 0.2; 
	


	for( iter = 1; iter < num_iter ; iter++){

		
		for( idx_d = 0; idx_d < nDim; idx_d++){
			grad[ idx_d ] = 0.0;
			old_wgt[ idx_d + 1] = wgt[ idx_d + 1];
		}


		for( idx_n = 0; idx_n < nSample; idx_n++){

			data_x = pair_data[ idx_n ]->x + 1;
			data_y = pair_data[ idx_n ]->y + 1;

		   	sum = 0.0;
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
			   	sum += (wgt[ idx_d + 1] * ( data_x[ idx_d ] - data_y[idx_d]));
			}
			sigmoid = 1.0 / ( 1 + exp( sum ));
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
				grad[ idx_d ] += (sigmoid * (data_x[ idx_d ] - data_y[ idx_d  ]));
			}

		}

		for( idx_d = 0; idx_d < nDim; idx_d++){
			grad[ idx_d ] = 2.0*grad[ idx_d ] - lambda * wgt[ idx_d + 1];
		}
		



		if( iter == 1 ){
			memcpy( dir_u, grad, nDim * sizeof( double ));
		}
		else{

			//dir_u = cg_dir ( dir_u, grad, old_grad, nDim);

			for( idx_d = 0; idx_d < nDim; idx_d++){
				delta[ idx_d ] = grad[ idx_d ] - old_grad[ idx_d ];
			}

			gd = 0.0;
			ud = 0.0;
			for( idx_d = 0; idx_d < nDim; idx_d++){
				gd += ( grad[ idx_d ] * delta[ idx_d ]);
				ud += ( dir_u[ idx_d ] * delta[ idx_d ]);
			}

			beta = gd / ud;

			for( idx_d = 0; idx_d < nDim; idx_d++){
				dir_u [ idx_d ] = grad [ idx_d ] - beta * dir_u[ idx_d ];
			}

		}


		// line search along dir_u

		ug = 0.0;
		uhu = 0.0;
		for( idx_d = 0; idx_d < nDim; idx_d++){
			ug += ( dir_u[ idx_d ] * grad[ idx_d] );
			uhu += ( dir_u[ idx_d ] * dir_u[ idx_d ] );
		}
		uhu *= lambda;

		for( idx_n = 0; idx_n < nSample; idx_n++){

			data_x = pair_data[ idx_n ]->x + 1;
			data_y = pair_data[ idx_n ]->y + 1;

		   	ux = 0.0;
			sum = 0.0;
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
				ux += (dir_u[ idx_d ] * (data_x[ idx_d ] - data_y[ idx_d]));
			   	sum += (wgt[ idx_d + 1] * (data_x[ idx_d ] - data_y[idx_d]));
			}
			tmp_f = exp( sum );
			uhu += (( ux * ux ) / ( 2 + tmp_f + 1.0/tmp_f )); 

		}


		tmp_f = ug / (2.0 * uhu) ;
		for( idx_d = 0; idx_d < nDim; idx_d++){
			wgt[ idx_d + 1] = wgt[ idx_d + 1] + tmp_f * dir_u[ idx_d ];
			old_grad[ idx_d ] = grad[ idx_d ];
		}

		for( idx_d = 0; idx_d <= nDim; idx_d++){
			tmp_f = wgt[ idx_d ]  ;
			wgt_sum[ idx_d ] += tmp_f;
			wgt_sq_sum[ idx_d ] += (tmp_f * tmp_f );
		}

#if 0
		dev = 0.0;
		for( idx_d = 0; idx_d <= nDim; idx_d++){
			tmp_f = ( old_wgt[ idx_d ] - wgt[ idx_d ] ) ;
			dev += ( tmp_f * tmp_f);
			se_beta[ idx_d ] = fabs(tmp_f);

			fprintf(stderr,"%d %lf %lf %lf\n", idx_d, old_wgt[ idx_d ], wgt[idx_d ], tmp_f);
		}
		converge = sqrt(dev / nDim);
#endif

#if 1
		// converge check
		dev = 0.0;
		for( idx_n = 0; idx_n < nSample; idx_n++){

			data_x = pair_data[ idx_n ]->x + 1;
			data_y = pair_data[ idx_n ]->y + 1;

			// class +1
			sum = 0.0;
		   	for( idx_d = 0; idx_d < nDim; idx_d++){
			   	sum += (wgt[ idx_d + 1] * (data_x[ idx_d ] - data_y[ idx_d ]));
			}
			dev += log( 1.0 + exp(-sum));

		}

		sum = 0.0;
		for( idx_d = 0; idx_d < nDim; idx_d++){
			sum += (wgt[ idx_d + 1] * wgt[ idx_d + 1]);
		}
		dev = 2.0*dev + 0.5 * lambda * sum ;


		converge = fabs( 1.0 - old_dev / dev );
#endif
		if( converge < EPSILON) break;
		/////////////////////
		if( (iter + 1 ) % 5  == 0){
		   	fprintf(stderr, " iter[%d]:\t%lf\t%lf\n", iter, dev, converge);
		}


		old_dev = dev;



	}
   
	fprintf(stderr, " iter[%d]:\t%lf\t%lf\n", iter, dev, converge);

#if 0
	for( idx_n = 0; idx_n < nSample; idx_n++){

		data_x = pair_data[ idx_n ]->x + 1;
		data_y = pair_data[ idx_n ]->y + 1;

		for( idx_d = 0; idx_d < nDim; idx_d++){
			tmp_f = (data_x[ idx_d ] - data_y[ idx_d]);
			tmp_f = tmp_f * wgt[idx_d + 1];
			wgt_sum[idx_d + 1] += (tmp_f);
			wgt_sq_sum[idx_d + 1] += (tmp_f * tmp_f);
		}
	}

#endif
	for( idx_d = 0; idx_d <= nDim; idx_d++){
		wgt_sum[ idx_d ] /= (double)iter;
		wgt_sq_sum[ idx_d ] /= (double)iter;

		tmp_f = wgt_sq_sum[idx_d] - (wgt_sum[idx_d] * wgt_sum[idx_d]);

		if( tmp_f > 0.0 ){
		   	se_beta[ idx_d ] = sqrt(tmp_f / iter);
		}
		else 
		   	se_beta[ idx_d ] = 0.0;
	}

	free(delta );
	free(grad );
	free(old_grad );
	free(old_wgt );
	free(dir_u );

	free( wgt_sum );
	free( wgt_sq_sum );



	return 0;
}


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))


//////////////////////////////////////////////////////////////////////

void EstiLogisticRegression(int ngroups, double **matrix, int nfeats, double *estimated, double *beta)
{
    int    idx, idx2    ;
    double esti, expsum ;

    for (idx = 1 ; idx < ngroups+1 ; ++idx) {
        expsum = -beta[0] ;
        for (idx2 = 1 ; idx2 < nfeats+1 ; ++idx2)
            expsum += (-beta[idx2]*matrix[idx][idx2]) ;
        esti = 1.0/(1.0+exp(expsum)) ;
        estimated[idx] = esti ;
    }
}

//////////////////////////////////////////////////////////////////////

//void EstiBTStrength(int ngroups, double **matrix, int nfeats, double *estimated, double *beta, int noexp)
void EstiBTStrength(int ngroups, double **matrix, int nfeats, double *beta, int noexp, char *outputfile)
{
    int    idx, idx2    ;
    double expsum ;

	FILE *fptr ;

	if ((fptr = fopen(outputfile,"wt")) == NULL) {
		fprintf(stderr,"Error in writing %s\n", outputfile) ;
		exit(1) ;
	}

    if (noexp) {

        for (idx = 0 ; idx < ngroups ; ++idx) {
            expsum = 0 ; 
            for (idx2 = 1 ; idx2 < nfeats ; ++idx2)
                expsum += beta[idx2]*matrix[idx][idx2] ;
            //estimated[idx] = expsum ;
		   	fprintf(fptr,"%d %f %f\n", (int) matrix[idx][0], matrix[idx][nfeats], expsum) ;
        }

    } else {

        for (idx = 0 ; idx < ngroups ; ++idx) {
            expsum = 0 ; 
            for (idx2 = 1 ; idx2 < nfeats ; ++idx2)
                expsum += beta[idx2]*matrix[idx][idx2] ;
		   	fprintf(fptr,"%d %f %f\n", (int) matrix[idx][0], matrix[idx][nfeats], expsum) ;
        }
    }

	fclose(fptr);
}

