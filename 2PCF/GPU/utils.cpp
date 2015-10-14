/*
 * This project is dual licensed. You may license this software under one of the
   following licences:

   + Creative Commons Attribution-Share Alike 3.0 Unported License
     http://creativecommons.org/licenses/by-nc-sa/3.0/

   + GNU GENERAL PUBLIC LICENSE v3, designated as a "BY-SA Compatible License"
     as defined in BY-SA 4.0 on 8 October 2015

 * See the LICENSE file in the root directory of this source tree for full
   copyright disclosure, and other details.
 */

 /* This function counts the number of columns in a file */
/* NOTE that this is done only checking the first row */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>

#include "utils.hpp"

int cols_number(char *input_file)
{

    /* Definition of Variables */

    char row[650];
    char * pch;
    int columns=0;

    /* Opening the input file */

    std::ifstream f_in(input_file, std::ios::in);

    /* Reading the first line */

    f_in.getline(row,650,'\n');

    /* Closing file */

    f_in.close();

    /* Counting columns */

    pch = strtok (row,"'\t'' '");
    while (pch != NULL)
    {
        columns++;
        pch = strtok (NULL, "'\t'' '");
    }
    return (columns);
}

/* This function counts the number of rows in a file */

int counting_lines(char *input_file)
{

    /* Definition of variables */

    int lines=0;
    char line[650];

    /* Opening the input file */

    std::ifstream f_in;
    f_in.open(input_file, std::ios::in);

    /* Counting lines */

    while(!f_in.eof())
    {
        if(f_in.getline(line, 650, '\n' ) != NULL)
        {
            lines=lines+1;
        }
    }
    lines=lines-1;

    /* Closing file */

    f_in.close();
    
    return(lines);
}

/* This function checks the input data */

int verification(int argc, char *argv[])
{

    /* Definition of variables */
   
    int points_per_degree;
    int columns;
    int f;

    /* Opening, checking and closing the first input file */

    std::ifstream fileinputdata;
    for(f=1;f<4;f++){
        fileinputdata.open(argv[f], std::ios::in);
        if (fileinputdata == NULL)
        {
            printf("Error opening data file number %d \n",f);
            return(0);
        }
        fileinputdata.close();
    }

    /* Checking other input parameters */

    if (argv[4]==NULL)
    {
        printf("You must introduce a number of points that you want per degree");
        return(0);
    }

    points_per_degree=atoi(argv[4]);
 

    if (points_per_degree<4 or points_per_degree%4!=0 or points_per_degree>64)
    {
        printf("The points per degree must be [4,8,16,32,64]");
        return(0);
    }

    if (argv[5]==NULL)
    {
        printf("You must introduce a name for the output file");
        return(0);
    }

    /* Checking cols number in every input file */

    for(f=1;f<4;f++){
        columns=cols_number(argv[f]);
        if (columns != 2 )
        {
            printf("Number of columns in file number %d must be exactly 2 and the first one has to be the right ascension and the second one the declination, both in degrees",f);
            return(0);
        }
    }

    return(1);
}

/* Equatorial to cartesian conversion */ 

int eq2cart(char *filename, int nlines, float *xd, float *yd, float *zd){
        std::ifstream infile(filename);
        int n;
        float ra, dec, phi, cth, th;

        /* We will store the data in cartesian coordinates, that's why we pass the ecuatorial coordinates to spheric coordinates and then to cartesian coordinates */

        for (n=0;n<nlines;n++)
        {
            infile>>ra>>dec; /* reading ecuatorial coordinates */
            phi=ra*M_PI/180.0;
            cth=cos((90.0-dec)*M_PI/180.0);
            th=acos(cth);

            /* to cartesian coordinates */
            
            xd[n]=cos(phi)*sin(th);
            yd[n]=sin(phi)*sin(th);
            zd[n]=cth;    
        }

        infile.close();
        
        return(0);
}
