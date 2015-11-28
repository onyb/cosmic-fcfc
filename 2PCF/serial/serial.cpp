/* Header files */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "utils.hpp"

using namespace std;


/* Binning Functions */

#define BIN_SIZE 256

/* This function counts the number of pairs in the data file */
/* We will use this kernel to calculate real-real pairs and random-random pairs */

void binning(float *xd,float *yd,float *zd,float *ZZ,int number_lines,int points_per_degree, int number_of_degrees)
{
    float angle;

    float x,y,z; //MCM
    float xx,yy,zz; //MCM

    /* We start the counting */

    for (int i=0;i<number_lines;i++)
    {
        x = xd[i];//MCM
        y = yd[i];//MCM
        z = zd[i];//MCM

        /* Beginning the sequential calculation case (CPU). We use "while" rather than "if" as recommended in the book "Cuda by Example" */

        for(int j = 0; j < number_lines; ++j)
        {
            xx = xd[j];//MCM
            yy = yd[j];//MCM
            zz = zd[j];//MCM

            /* We make the dot product */
            angle = x * xx + y * yy + z * zz;//MCM

            /* Sometimes "angle" is higher than one, due to numnerical precision, to solve it we use the next sentence */

            angle=fminf(angle,1.0);
            angle=acosf(angle)*180.0/M_PI;

            /* We finally count the number of pairs separated at an angular distance "angle", always in shared memory */

            if(angle < number_of_degrees)
            {
                ZZ[int(angle*points_per_degree)] += 1.0;
            }
        }
    }
}

/* This function counts the number of pairs that there are between two data groups */
/* We will use this kernel to calculate real-random pairs and real_1-real_2 pairs (cross-correlation) */
/* NOTE that this kernel has NOT been merged with 'binning' above: this is for speed optimization, we avoid passing extra variables to the GPU */


void binning_mix(float *xd_real, float *yd_real, float *zd_real, float *xd_sim, float *yd_sim, float *zd_sim, float *ZY, int lines_number_1, int lines_number_2, int points_per_degree, int number_of_degrees)
{

    /* We define variables (arrays) in shared memory */

    float angle;
    float x,y,z; //MCM
    float xx,yy,zz; //MCM

    /* We start the counting */

    for (int i=0;i<lines_number_1;i++)
    {
        x = xd_real[i];//MCM
        y = yd_real[i];//MCM
        z = zd_real[i];//MCM

        /* Beginning the sequential calculation case (CPU). */

        for(int j = 0; j< lines_number_2; ++j)
        {
            xx = xd_sim[j];//MCM
            yy = yd_sim[j];//MCM
            zz = zd_sim[j];//MCM

            /* We make the dot product */
            angle = x * xx + y * yy + z * zz;//MCM

            /* Sometimes "angle" is higher than one, due to numnerical precision, to solve it we use the next sentence */

            angle=fminf(angle,1.0);
            angle=acosf(angle)*180.0/M_PI;

            /* We finally count the number of pairs separated an angular distance "angle", always in shared memory */

            if(angle < number_of_degrees)
            {
                ZY[int(angle*points_per_degree)] += 1.0;
            }
        }
    }
}


/* Main Function*/

int main(int argc, char *argv[])
{
    /* Checking if the input files and call to script meet the requirements */

    const int mode = verification(argc, argv);
    if (mode == 0)
    {
        exit(1);
    }

    /* Definition of variables, these variables are the same for the auto-correlation and the cross-correlation */

    int points_per_degree;
    int number_of_degrees;
    char *input_real_file_1, *input_real_file_2, *input_random_file, *output_file;
    double W, angle_theta, poissonian_error, norm_cost_1, norm_cost_2;

    int real_lines_number_1, real_lines_number_2, random_lines_number, max_lines;

    float *xd_real_1,*yd_real_1,*zd_real_1;
    float *xd_real_2,*yd_real_2,*zd_real_2;
    float *xd_rand,*yd_rand,*zd_rand;

    /* Assignment of Variables with inputs */

    input_real_file_1 = argv[1];
    input_real_file_2 = argv[2];
    input_random_file = argv[3];
    points_per_degree = atoi(argv[4]);
    number_of_degrees = int(256/float(points_per_degree));
    output_file=argv[5];

    /* Counting lines in every input file */

    real_lines_number_1 = count_lines(input_real_file_1);
    if(real_lines_number_1 == -1)
        std::cerr << "Incorrectly formatted file: " << input_real_file_1 << std::endl;

    /* const int mode = (std::string(argv[1]) == std::string(argv[2])) ? AUTO : CROSS; */

    if(mode == CROSS){
        real_lines_number_2 = count_lines(input_real_file_2);
        if(real_lines_number_2 == -1)
            std::cerr << "Incorrectly formatted file: " << input_real_file_2 << std::endl;
    }

    random_lines_number = count_lines(input_random_file);
    if(random_lines_number == -1)
        std::cerr << "Incorrectly formatted file: " << input_random_file << std::endl;

    /* We define variables to store the real,random data */

    xd_real_1 = (float *)malloc(real_lines_number_1 * sizeof (float));
    yd_real_1 = (float *)malloc(real_lines_number_1 * sizeof (float));
    zd_real_1 = (float *)malloc(real_lines_number_1 * sizeof (float));
    xd_rand = (float *)malloc(random_lines_number * sizeof (float));
    yd_rand = (float *)malloc(random_lines_number * sizeof (float));
    zd_rand = (float *)malloc(random_lines_number * sizeof (float));

    if(mode == CROSS)
    {
        xd_real_2 = (float *)malloc(real_lines_number_2 * sizeof (float));
        yd_real_2 = (float *)malloc(real_lines_number_2 * sizeof (float));
        zd_real_2 = (float *)malloc(real_lines_number_2 * sizeof (float));
    }

    /* Opening the first input file */

    eq2cart(input_real_file_1,real_lines_number_1,xd_real_1,yd_real_1,zd_real_1);

    /* Opening the second input file */

    if(mode == CROSS)
    {
        eq2cart(input_real_file_2,real_lines_number_2,xd_real_2,yd_real_2,zd_real_2);
    }

    /* Opening the third input file */

    eq2cart(input_random_file,random_lines_number,xd_rand,yd_rand,zd_rand);


    /* We define variables to store the pairs between real data (DD), random data (RR) and both together (DR) */

    /* on CPU */
    float *DD;
    float *DR;
    float *RR;
    float *D1D2;
    float *D1R;
    float *D2R;
    RR = (float *)malloc(BIN_SIZE*sizeof(float));
    if(mode == CROSS)
    {
        D1D2 = (float *)malloc(BIN_SIZE*sizeof(float));
        D1R = (float *)malloc(BIN_SIZE*sizeof(float));
        D2R = (float *)malloc(BIN_SIZE*sizeof(float));
        for (int i = 0; i < BIN_SIZE; i++)
        {
            D1D2[i] = 0.0;
            RR[i] = 0.0;
            D1R[i] = 0.0;
            D2R[i] = 0.0;
        }
    }
    else
    {
        DD = (float *)malloc(BIN_SIZE*sizeof(float));
        DR = (float *)malloc(BIN_SIZE*sizeof(float));
        for (int i=0; i< real_lines_number_1; i++)
        {
           DD[i] = 0.0;
           RR[i] = 0.0;
           DR[i] = 0.0;
        }
    }

    /* We determine which is the maximum number of lines */
    max_lines = max(real_lines_number_1,random_lines_number);
    if(mode == CROSS)
    {
        max_lines = max(max_lines,real_lines_number_2);
    }

    if(mode == CROSS)
    {
        /* Loop Part */
        binning(xd_rand, yd_rand, zd_rand, RR, random_lines_number, points_per_degree, number_of_degrees);
        binning_mix(xd_real_1, yd_real_1, zd_real_1, xd_real_2, yd_real_2, zd_real_2, D1D2, real_lines_number_1, real_lines_number_2, points_per_degree, number_of_degrees);
        binning_mix(xd_real_1, yd_real_1, zd_real_1, xd_rand, yd_rand, zd_rand, D1R, real_lines_number_1, random_lines_number, points_per_degree, number_of_degrees);
        binning_mix(xd_real_2, yd_real_2, zd_real_2, xd_rand, yd_rand, zd_rand, D2R, real_lines_number_2, random_lines_number, points_per_degree, number_of_degrees);
    }
    else
    {
        /* Loop Part */
        binning(xd_real_1, yd_real_1, zd_real_1, DD, real_lines_number_1, points_per_degree, number_of_degrees);
        binning(xd_rand, yd_rand, zd_rand, RR, random_lines_number, points_per_degree, number_of_degrees);
        binning_mix(xd_real_1, yd_real_1, zd_real_1, xd_rand, yd_rand, zd_rand, DR, real_lines_number_1, random_lines_number, points_per_degree, number_of_degrees);
    }

    /* Opening the output file */

    std::ofstream f_out(output_file);

   /* We calculate the normalization factor */

   norm_cost_1=float(random_lines_number)/float(real_lines_number_1);
   if(mode == CROSS)
   {
       norm_cost_2=float(random_lines_number)/float(real_lines_number_2);
   }

   if(mode == CROSS)
   {

        for (int i=1;i<real_lines_number_1;i++)
        {
            /* The angle corresponding to the W value */

            angle_theta=(1.0/points_per_degree)/2.0+(i*(1.0/points_per_degree));

            /* We are calculating the Landy & Szalay estimator */

            W=(norm_cost_1)*(norm_cost_2)*(D1D2[i]/RR[i])-(norm_cost_1)*(D1R[i]/RR[i])-(norm_cost_2)*(D2R[i]/RR[i])+1.0;
            poissonian_error=(1.0+W)/sqrt(D1D2[i]);
            f_out<<angle_theta<<"\t"<<W<<"\t"<<poissonian_error<<"\t"<<D1D2[i]<<"\t"<<D1R[i]<<"\t"<<D2R[i]<<"\t"<<RR[i]<<endl;
        }
    }
    else
    {
        for (int i=0;i<real_lines_number_1;i++)
        {
            /* The angle corresponding to the W value */

            angle_theta=(1.0/points_per_degree)/2.0+(i*(1.0/points_per_degree));
            W=(((pow(norm_cost_1,2)*DD[i])-(2*norm_cost_1*DR[i]))/RR[i])+1.0;
            poissonian_error=(1.0+W)/sqrt(DD[i]);
            f_out<<angle_theta<<"\t"<<W<<"\t"<<poissonian_error<<"\t"<<DD[i]<<"\t"<<DR[i]<<"\t"<<RR[i]<<endl;
        }
    }

    /* Closing output files */

    f_out.close();

    /* Freeing memory on the CPU */

    free(xd_real_1);
    free(yd_real_1);
    free(zd_real_1);
    free(xd_rand);
    free(yd_rand);
    free(zd_rand);

    if(mode == CROSS)
    {
        free(xd_real_2);
        free(yd_real_2);
        free(zd_real_2);
        free(D1D2);
        free(D1R);
        free(D2R);
    }
    else
    {
        free(DD);
        free(DR);
        free(RR);
    }

    return(0);
}

