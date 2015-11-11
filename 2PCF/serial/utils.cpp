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
#include <sstream>

#include "utils.hpp"

int count_lines(const char *file)
{
    std::ifstream infile(file);
    std::string line;
    int c = 0;
    while (std::getline(infile, line))
    {
        std::istringstream buffer(line);
        float a, b;
        buffer >> a >> b;
        if (!buffer || !buffer.eof())
            return -1;
        c++;
    }
    return c;
}

/* This function checks the input data */

int verification(int argc, char *argv[])
{

    /* Definition of variables */

    int points_per_degree;
    int columns;
    int f;
    int mode;

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

    /* If all input files are okay, determine and set the mode */

    mode = (std::string(argv[1]) == std::string(argv[2])) ? AUTO : CROSS;

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

    return(mode);
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
