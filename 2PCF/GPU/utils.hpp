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

#define AUTO 1
#define CROSS 2

int count_lines(const char *file);

int cols_number(char *input_file);

int counting_lines(char *input_file);

int verification(int argc, char *argv[]);

int eq2cart(char *filename, int nlines, float *xd, float *yd, float *zd);
