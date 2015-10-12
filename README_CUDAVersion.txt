
/* Authors: Miguel Cardenas-Montes (1)
            Rafael Ponce (2)
            Ignacio Sevilla (3)

(1): miguel.cardenas@ciemat.es
(2): rafael.ponce@ciemat.es
(3): ignacio.sevilla@ciemat.es

NOTE: This code has been successfully tested on NVIDIA GTX295, C1060, C2050, C2070 and GT555M. The latter is a gaming GPU, for this reason does not accept too many points */

               /* TO COMPILE: IT'S CLEAR THAT COMPILATION DEPENDS ON YOUR PLATFORM */

/* For Ubuntu with Cuda Libraries you have to do once:

    user@user:yourpath$ nvcc -c -arch=sm_20 GP2PCF_GPU_v0.2.cu
    user@user:yourpath$ g++ -lcudart GP2PCF_GPU_v0.2.o -o GP2PCF_GPU

    NOTE: Sometimes could be necesary to do: (m value depend on your architecture, 32 bits or 64 bits)
        nvcc -c -m=32 -arch=sm_20 GP2PCF_GPU_v0.2.cu   or
        nvcc -c -m=64 -arch=sm_20 GP2PCF_GPU_v0.2.cu


        If you have to indicate where the libraries are you will have to do something like that:

        nvcc -arch=sm_20 GP2PCF_GPU_v0.2.cu
        g++ -lcudart -lpthread -L/usr/local/cuda/lib64 GP2PCF_GPU_v0.2.o -o GP2PCF_GPU.cu


Now, if everything went well, you just created a binary file called "GP2PCF_GPU" */


                                            /* TO RUN */

/* This program calculates the 2pt-angular auto-correlation and cross-correlation functions by brute force using a GPU
The program receives as input parameters the following:

    Catalog_1: The path of your galaxy catalog, usually is a galaxy catalog for a specific redshift bin
    Catalog_2: The path of the same galaxy catalog you put as the first input (if you want to do the auto-correlation function) or another galaxy catalog (if you want to do the cross-correlation function between two redshift bins) 
    Random_catalog: The path corresponding to the random catalog
    Number_of_points: It's the number of bins per degree that you want
    Output_correlation_file: The path of your output file 

And the program will return inside the output file six or seven columns, depend on whether you did the auto-correlation or the cross-correlation function.

The input files must have only two columns:

    RA: Right Ascension in degrees 
    DEC: Declination also in degrees

The output will be:

    For the auto-correlation:

    angle(degrees) 
    W(corr-func.) 
    Poissonian Error 
    Real-Real Pairs 
    Real-Rand Pairs 
    Rand-Rand Pairs

    For the cross-correlation:

    angle(degrees) 
    W(corr-func.) 
    Poissonian Error 
    Real_1-Real_2 Pairs 
    Real_1-Rand Pairs 
    Real_2-Rand Pairs
    Rand-Rand Pairs

To run the program:

If you want to do calculate the auto-correlation you must do:

        ./GP2PCF_GPU Catalog_1.txt Catalog_1.txt Random_cat.txt points output.txt

If you want to do an cross-correlation you must do:

        ./GP2PCF_GPU Catalog_1.txt Catalog_2.txt Random_cat.txt points output.txt

EXAMPLES TO CALL THE FUNCTION:

   To do the correlation with 16 points per degree

       ./GP2PCF_GPU test.cat test.cat test_rand.cat 16 W_auto_test.txt

   To do cross-correlation with 32 points per degree

       ./GP2PCF_GPU test.cat test2.cat test_rand.cat 32 W_cross_test.txt    


IMPORTANT: 

1.- Note that if the first and second file have the same name, the program will calculate the auto-correlation function. If not, the program will calculate the cross-correlation function

2.- For memory reasons the points per degree should be between 4 points and 64 points and it must be multiple of 4.

3.- For hardware limitations, there is a relationship between the points per degree and the degree you can reach. You can see this relationship in the next table

Points	Angle
4	64
8	32
16	16
32	8
64	4 */
