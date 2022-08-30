#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

typedef struct DEF_tables {
    char *eosName;
    int N_mA;
    int N_alpha0;
    int N_beta0;
    double *alpha0_table;
    double *beta0_table;
    double ***mA_table;
    double ***alphaA_table;
    double ***betaA_table;
    double ***kA_table;
} DEF_tables;

typedef struct DEF_factors {
    double alpha0;
    double beta0;
    double mA;
    double alphaA;
    double betaA;
    double kA;
    char *companion_type;
    double mB;
    double alphaB;
    double betaB;
    double kB;
} DEF_factors;

typedef struct DEF_PPK_derivatives {
    double dk_dm;
    double dk_dm2;
    double dgamma_dm;
    double dgamma_dm2;
    double dsi_dm;
    double dsi_dm2;
    double ddr_dm;
    double ddr_dm2;
    double ddth_dm;
    double ddth_dm2;
    double dpbdot_dm;
    double dpbdot_dm2;
} DEF_PPK_derivatives;

void DEF_tables_load(DEF_tables *tables);
void DEF_tables_free(DEF_tables *tables);

void allocate_DEF_arrays(pulsar *psr, DEF_tables *tables);
void interpolate_DEF_arrays(pulsar *psr, DEF_tables *tables);
void interpolate_DEF_params(pulsar *psr, double mA, double *alphaA, double *betaA, double *kA);

void interpolate_DEF_factors(pulsar *psr, DEF_factors *factors);

double*** allocate_3Darray(int dim1, int dim2, int dim3);
void read_3Darray(char *fileName, double ***array, int dim1, int dim2, int dim3);
void free_3Darray(double ***array, int dim1, int dim2, int dim3);
int locate_in_array(double value, double *array, int N);
double bilinear_interpolation(double f00, double f10, double f01, double f11, double x, double y);
double interpolate_linear(double x_value, int ind, double *x_array, double *y_array);