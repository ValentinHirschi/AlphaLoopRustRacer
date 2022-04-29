#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <mp++/complex128.hpp>

using namespace mppp;
using namespace mppp::literals;

void load_inputs(double lm[%(n_lm)d][4][%(n_inputs)d], double sglm[%(n_sglm)d][4][%(n_inputs)d], double extm[%(n_extm)d][4][%(n_inputs)d], double sgextm[%(n_sgextm)d][4][%(n_inputs)d], double target_res[%(n_inputs)d], int *n_lm, int *n_ext, int *n_sg_ext, int *n_inputs ) {

    FILE *myfile;
    double buff[(%(n_lm)d+%(n_extm)d+%(n_sgextm)d)*4+1];
    int i,j,k;

    myfile=fopen("./%(input_file_name)s", "r");

    fscanf(myfile,"%%d",n_lm);
    fscanf(myfile,"%%d",n_ext);
    fscanf(myfile,"%%d",n_sg_ext);
    fscanf(myfile,"%%d",n_inputs);
    
    for(i = 0; i < (*n_inputs); i++)
    {
        for (j = 0 ; j <= (%(n_lm)d+%(n_extm)d+%(n_sgextm)d)*4; j++)
        {
            fscanf(myfile,"%%lf",&buff[j]);
        }
        for (j = 0; j <= (*n_lm)-1; j++) {
            for (k = 0; k < 4; k++) {
                lm[j][k][i] = buff[j*4+k];
            }
        }
        for (j = 0; j <= (*n_lm)+(*n_ext)-1; j++) {
            for (k = 0; k < 4; k++) {
                sglm[j][k][i] = buff[j*4+k];
            }
        }
        for (j = (*n_lm); j <= (*n_lm)+(*n_ext)-1; j++) {
            for (k = 0; k < 4; k++) {
                extm[j-(*n_lm)][k][i] = buff[j*4+k];
            }
        }
        for (j = (*n_lm)+(*n_ext); j <= (*n_lm)+(*n_ext)+(*n_sg_ext)-1; j++) {
            for (k = 0; k < 4; k++) {
                sgextm[j-(*n_lm)-(*n_ext)][k][i] = buff[j*4+k];
            }
        }
        target_res[i] = buff[((*n_lm)+(*n_ext)+(*n_sg_ext))*4];
    }

    fclose(myfile);
}

void race(real128 lm[%(n_lm)d][4], real128 sglm[%(n_sglm)d][4], real128 extm[%(n_extm)d][4], real128 sgextm[%(n_sgextm)d][4], double *out_res) {

    real128 em[%(n_edges)d][4];
    real128 sd[%(n_edges)d][%(n_edges)d];
    real128 em_osE[%(n_edges)d];
    real128 num, denom;
    const double pi=acos(-1.);
    const bool DBG=false;
    
    real128 res;
    int i,j;

    if (DBG) {
        for (i = 0; i < %(n_lm)d; i++) {
            printf("lm[%%d]=(%%.16e,%%.16e,%%.16e,%%.16e)\n",i, ((double)lm[i][0]), ((double)lm[i][1]), ((double)lm[i][2]), ((double)lm[i][3]));
        }
        for (i = 0; i < %(n_sglm)d; i++) {
            printf("sglm[%%d]=(%%.16e,%%.16e,%%.16e,%%.16e)\n",i, ((double)sglm[i][0]), ((double)sglm[i][1]), ((double)sglm[i][2]), ((double)sglm[i][3]));
        }
        for (i = 0; i < %(n_extm)d; i++) {
            printf("extm[%%d]=(%%.16e,%%.16e,%%.16e,%%.16e)\n",i, ((double)extm[i][0]), ((double)extm[i][1]), ((double)extm[i][2]), ((double)extm[i][3]));
        }
        for (i = 0; i < %(n_sgextm)d; i++) {
            printf("sgextm[%%d]=(%%.16e,%%.16e,%%.16e,%%.16e)\n",i, ((double)sgextm[i][0]), ((double)sgextm[i][1]), ((double)sgextm[i][2]), ((double)sgextm[i][3]));
        }
    }
    
    res = real128(0.q);

//           Evaluae all common quantities
%(warmup_code)s
    
//           Now evaluate all cuts
%(evaluate_ltd_cut)s

    *out_res = ((double) res);

//           Finalise computation
%(wrapup_code)s

    if (DBG) {
        printf("Final result: %%.16e\n",((double)*out_res));
    }
}

int main() {
    double lm[%(n_lm)d][4][%(n_inputs)d], sglm[%(n_sglm)d][4][%(n_inputs)d], extm[%(n_extm)d][4][%(n_inputs)d], sgextm[%(n_sgextm)d][4][%(n_inputs)d], target_res[%(n_inputs)d];
    double res[%(n_inputs)d];
    int n_lm, n_ext, n_sg_ext, n_inputs;
    load_inputs(lm, sglm, extm, sgextm, target_res, &n_lm, &n_ext, &n_sg_ext, &n_inputs);

    real128 a_lm[3][4], a_sglm[4][4], a_extm[1][4], a_sgextm[1][4];
    int i,j,k,l;
    int n_runs=1;

    struct timespec start, end;

    clock_gettime(CLOCK_REALTIME, &start);
    for (i = 0; i < n_runs; i++) {
        for (j = 0; j < n_inputs; j++) {
            for (k = 0; k < n_lm; k++) {
                for (l = 0; l < 4; l++) {
                    a_lm[k][l]=real128(lm[k][l][j]);    
                }
            }
            for (k = 0; k < n_lm+n_ext; k++) {
                for (l = 0; l < 4; l++) {
                    a_sglm[k][l]=real128(sglm[k][l][j]);
                }
            }
            for (k = 0; k < n_ext; k++) {
                for (l = 0; l < 4; l++) {
                    a_extm[k][l]=(extm[k][l][j]); 
                }
            }
            for (k = 0; k < n_sg_ext; k++) {
                for (l = 0; l < 4; l++) {
                    a_sgextm[k][l]=real128(sgextm[k][l][j]);
                }
            }
            race(a_lm, a_sglm, a_extm, a_sgextm, &res[j]);
        }
    }
    clock_gettime(CLOCK_REALTIME, &end);

    double t_eval = ((end.tv_sec - start.tv_sec)*1.0e6 + (end.tv_nsec - start.tv_nsec) / 1.0e3)/((double) n_runs*n_inputs);
    printf("%%.10e\n", t_eval);

    double max_diff = -1.0;
    int max_diff_index = 0;
    for (j = 0; j < n_inputs; j++) {
        if (fabs((res[j]-target_res[j])/target_res[j])>max_diff) {
            max_diff = fabs((res[j]-target_res[j])/target_res[j]);
            max_diff_index = j;
        }
    }
    printf("%%.3e\n", max_diff*100.);
 }