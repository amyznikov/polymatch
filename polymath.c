/*
 * polymath.c
 *
 *  Created on: Oct 18, 2016
 *      Authors:
 *        Andrey Myznikov <andrey.myznikov@gmail.com>
 *        Igor Uspeniev <igor.uspeniev@gmail.com>
 *
 *  http://www.studfiles.ru/preview/4466389/
 */

#define __GNU_SOURCE

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <malloc.h>
#include "polymath.h"

#define MAX_POINTS    100
#define MAX_POLYGONS ((MAX_POINTS)*(MAX_POINTS)*(MAX_POINTS)*(MAX_POINTS))

#define PDBG(...)\
    fprintf(stderr,"%s(): %d ", __func__, __LINE__), \
    fprintf(stderr,__VA_ARGS__), \
    fputs( "\n",stderr)


////////////////////////////////////////////////////////////
struct poly {
  uint16_t idx[4];
  float d1, d2;
  //double cx, cy;
};





////////////////////////////////////////////////////////////



size_t load_points_from_file(FILE * fp, float x[], float y[], size_t maxsize)
{
  size_t size = 0;
  while ( size < maxsize && fscanf(fp, "%f %f", &x[size], &y[size]) == 2 ) {
    ++size;
  }

  if ( !feof(fp) ) {
    PDBG("WARNING: EOF NOT REACHED: There are still some unparesed data");
  }

  return size;
}


static void dump_points_to_file(FILE * fp, const float x[], const float y[], size_t n)
{
  for ( size_t i = 0; i < n; ++i ) {
    fprintf(fp, "%3zu\t%8g\t%8g\n", i, x[i], y[i]);
  }
}


static size_t load_points(const char * filename, float x[], float y[], size_t maxsize)
{
  size_t size = 0;
  FILE * fp = NULL;

  if ( !(fp = fopen(filename, "r")) ) {
    PDBG("fopen(%s) fails: %s", filename, strerror(errno));
  }
  else if ( !(size = load_points_from_file(fp, x, y, maxsize)) ) {
    PDBG("load_points_from_file(%s) fails: %s", filename, strerror(errno));
  }

  if ( fp ) {
    fclose(fp);
  }

  return size;
}




static void dump_polygons_to_file(FILE * fp, const struct poly poly[], size_t n)
{
  size_t i0, i1, i2, i3;

  for ( size_t i = 0; i < n; ++i ) {
    const struct poly * p = poly + i;

    i0 = p->idx[0];
    i1 = p->idx[1];
    i2 = p->idx[2];
    i3 = p->idx[3];

    fprintf(fp, "%4zu\t"
        "[%4zu\t%4zu\t%4zu\t%4zu]\t"


        "[%9.3e\t%9.3e]\n",
        i,
        i0, i1, i2, i3,


        p->d1, p->d2

        );
  }

  //"[%9.3e\t%9.3e\t%9.3e\t%9.3e]\t"
  // p->e1, p->e2, 1.0/p->e1, 1.0/p->e2

  // "[%9.3e\t%9.3e]\t"
  // p->cx, p->cy,
}




static void dump_polygons(const char * filename, const struct poly poly[], size_t n)
{
  FILE * fp = NULL;

  if ( strcmp(filename, "stdout") == 0 ) {
    fp = stdout;
  }
  else if ( strcmp(filename, "stderr") == 0 ) {
    fp = stderr;
  }
  else if ( !(fp = fopen(filename, "w") ) ) {
    PDBG("fopen(%s) fails: %s", filename, strerror(errno));
  }


  if ( fp ) {

    dump_polygons_to_file(fp, poly, n);

    if ( fp != stdout && fp != stderr ) {
      fclose(fp);
    }
  }

}

static void calc_center(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double * restrict cx, double * restrict cy)
{
  double c1 = (x2 - x0) / (y2 - y0);
  double c2 = (x3 - x1) / (y3 - y1);
  *cx = (c1 * c2 * y1 - c1 * c2 * y0 - c1 * x1 + c2 * x0) / (c2 - c1);
  *cy = (c2 * y1 - c1 * y0 - x1 + x0) / (c2 - c1);
}


static size_t create_polygons(const float x[], const float y[], size_t N, struct poly poly[], size_t nmax)
{
  size_t n = 0;

  double a00 = 0, a01 = 0, a10 = 0, a11 = 0, b0 = 0, b1 = 0, det = 0, detx = 0, dety = 0;

  // now we are here
  for ( size_t i0 = 0; i0 < N-3 && n < nmax; ++i0 ) {

    const double x0 = x[i0];
    const double y0 = y[i0];

    for ( size_t i1 = i0 + 1; i1 < N-2 && n < nmax; ++i1 ) {

      const double x1 = x[i1];
      const double y1 = y[i1];

      for ( size_t i2 = i1 + 1; i2 < N-1 && n < nmax; ++i2 ) {

        for ( size_t i3 = i2 + 1; i3 < N && n < nmax; ++i3 ) {

          double x2, y2, x3, y3;
          double d1 = 0, d2 = 0;
          double cx, cy;

          bool swap = false;


          x2 = x[i2];
          y2 = y[i2];
          x3 = x[i3];
          y3 = y[i3];
          calc_center(x0, y0, x1, y1, x2, y2, x3, y3, &cx, &cy);

          if ( isnan(cx) || !isnan(cy) ) {
            swap = true;
          }
          else {
            d1 = hypot(cx - x0, cy - y0) / hypot(cx - x2, cy - y2);
            d2 = hypot(cx - x1, cy - y1) / hypot(cx - x3, cy - y3);
            if ( isnan(d1) || isnan(d2) || isinf(d2) || isinf(d2) || d1 < FLT_EPSILON || d2 < FLT_EPSILON ) {
              swap = true;
            }
          }


          if ( swap ) {

            x2 = x[i3];
            y2 = y[i3];
            x3 = x[i2];
            y3 = y[i2];
            calc_center(x0, y0, x1, y1, x2, y2, x3, y3, &cx, &cy);

            if ( isnan(cx) || isnan(cy) ) {
              // bad poly
              continue;
            }
            else {
              d1 = hypot(cx - x0, cy - y0) / hypot(cx - x2, cy - y2);
              d2 = hypot(cx - x1, cy - y1) / hypot(cx - x3, cy - y3);
              if ( isnan(d1) || isnan(d2) || isinf(d2) || isinf(d2) || d1 < FLT_EPSILON || d2 < FLT_EPSILON ) {
                // bad poly
                continue;
              }
            }
          }

          // save polygon

//          poly [n].cx = cx;
//          poly [n].cy = cy;
          poly [n].idx[0] = i0;
          poly [n].idx[1] = i1;
          poly[n].d1 = d1 < 1 ? 1 / d1 : d1;
          poly[n].d2 = d2 < 1 ? 1 / d2 : d2;

          if ( swap ) {
            poly [n].idx[2] = i3;
            poly [n].idx[3] = i2;
          }
          else {
            poly [n].idx[2] = i2;
            poly [n].idx[3] = i3;
          }

          if ( poly[n].d1 < poly[n].d2 ) {
            double tmp = poly[n].d1;
            poly[n].d1 = poly[n].d2;
            poly[n].d2 = tmp;
          }

          ++n;
        }

      }

    }

  }

  return n;
}


static inline bool eq(double d1, double d2, double tol)
{
  return fabs(d1 - d2) / fabs(d1 + d2) <= tol / 2;
}

/**
 * fill the match[] array with indexes of matching p2[] items
 */
static size_t match_polygons(double tol, const struct poly poly1[], size_t n1, const struct poly poly2[], size_t n2, uint16_t m1[], uint16_t m2[])
{
  const struct poly * p1, *p2;
  size_t n = 0;

  for ( size_t i = 0; i < n1; ++i ) {

    p1 = poly1 + i;

    for ( size_t j = 0; j < n2; ++j ) {

      p2 = poly2 + j;

      if ( (eq(p1->d1, p2->d1, tol) && eq(p1->d2, p2->d2, tol))
          || (eq(p1->d1, p2->d2, tol) && eq(p1->d2, p2->d1, tol)) ) {

        m1[n] = i;
        m2[n] = j;
        ++n;
      }

    }
  }

  return n;
}


static void dump_matches_to_file(FILE * fp,
    const float x1[], const float y1[], size_t N1,
    const float x2[], const float y2[], size_t N2,
    const struct poly poly1[], size_t np1,
    const struct poly poly2[], size_t np2,
    const uint16_t m1[], const uint16_t m2[], size_t nm)
{

  for ( size_t im = 0; im < nm; ++im ) {

    // i-th element of poly1[] matches to j-th element of poly2[]
    uint16_t i = m1[im];
    uint16_t j = m2[im];

    const struct poly * p1 = poly1 + i;
    const struct poly * p2 = poly2 + j;


    fprintf(fp,
        "%zu\t"
        "[%3u\t%3u\t%3u\t%3u]\t"
        "[%3u\t%3u\t%3u\t%3u]\t"
        "[%9.6e\t%9.6e]\t"
        "[%9.6e\t%9.6e]\t"
        "{%g %g\t%g %g\t%g %g\t%g %g}\t"
        "{%g %g\t%g %g\t%g %g\t%g %g}\n"
        ,

        im,
        p1->idx[0], p1->idx[1], p1->idx[2], p1->idx[3],
        p2->idx[0], p2->idx[1], p2->idx[2], p2->idx[3],
        p1->d1, p1->d2,
        p2->d1, p2->d2,
        x1[p1->idx[0]], y1[p1->idx[0]], x1[p1->idx[1]], y1[p1->idx[1]], x1[p1->idx[2]], y1[p1->idx[2]], x1[p1->idx[3]], y1[p1->idx[3]],
        x2[p2->idx[0]], y2[p2->idx[0]], x2[p2->idx[1]], y2[p2->idx[1]], x2[p2->idx[2]], y2[p2->idx[2]], x2[p2->idx[3]], y2[p2->idx[3]]
        );

//
//    "{%7g\t%7g\t%7g\t%7g}\n"
//
//
  }

}


static void affine_transform(float x[], float y[], size_t n, double x0, double y0, double s, double alpha)
{
  for ( size_t i = 0; i < n; ++i ) {
    double xx = x0 + s * (x[i] * cos(alpha) + y[i] * sin(alpha));
    double yy = y0 + s * (-x[i] * sin(alpha) + y[i] * cos(alpha));
    x[i] = xx, y[i] = yy;
  }

}




/*
 * Usage:
 *  polymath <first-file> <second-file>
 */

int main(int argc, char *argv[])
{
  const char * first_input_file_name = NULL;
  const char * second_input_file_name = NULL;


  static float x1[MAX_POINTS], x2[MAX_POINTS];
  static float y1[MAX_POINTS], y2[MAX_POINTS];
  size_t N1 = 0, N2 = 0;


  struct poly * poly1 = NULL, * poly2 = NULL;
  size_t n1 = 0, n2 = 0;

  size_t nmatch = 0, maxmatch = 0;
  uint16_t * m1 = NULL;
  uint16_t * m2 = NULL;


  // affine transform for second data
  double X0 = 0, Y0 = 0, S = 1, Alpha = 0;

  double tol = FLT_EPSILON;



  // Parse command line args

  for ( int i = 1; i < argc; ++i ) {

    if ( strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0  ){\
      printf("Usage:\n");
      printf(" polymath <first-file> <second-file> [x0=<x0>] [y0=<y0>] [s=<s>] [a=<a>] []tol=<match-tolerance-to-d>\n");
      return 0;
    }

    if ( strncmp(argv[i], "x0=", 3) == 0 ) {
      if ( sscanf(argv[i] + 3, "%lf", &X0) != 1 ) {
        fprintf(stderr, "Syntax Error in argument: %s\n", argv[i]);
        return 1;
      }
    }
    else  if ( strncmp(argv[i], "y0=", 3) == 0 ) {
      if ( sscanf(argv[i] + 3, "%lf", &Y0) != 1 ) {
        fprintf(stderr, "Syntax Error in argument: %s\n", argv[i]);
        return 1;
      }
    }
    else if ( strncmp(argv[i], "s=", 2) == 0 ) {
      if ( sscanf(argv[i] + 2, "%lf", &S) != 1 ) {
        fprintf(stderr, "Syntax Error in argument: %s\n", argv[i]);
        return 1;
      }
    }
    else if ( strncmp(argv[i], "a=", 2) == 0 ) {
      if ( sscanf(argv[i] + 2, "%lf", &Alpha) != 1 ) {
        fprintf(stderr, "Syntax Error in argument: %s\n", argv[i]);
        return 1;
      }
    }
    else if ( strncmp(argv[i], "tol=", 4) == 0 ) {
      if ( sscanf(argv[i] + 4, "%lf", &tol) != 1 ) {
        fprintf(stderr, "Syntax Error in argument: %s\n", argv[i]);
        return 1;
      }
    }
    else if ( !first_input_file_name ) {
      first_input_file_name = argv[i];
    }
    else if ( !second_input_file_name ) {
      second_input_file_name = argv[i];
    }
    else {
      fprintf(stderr, "Invalid argument %s\n", argv[i]);
      return 1;
    }

  }


  if ( !first_input_file_name || !second_input_file_name ) {
    fprintf(stderr, "Input file name not provided, try --help\n");
    return 1;
  }


  //
  // allocate work memory
  //

  if ( !(poly1 = malloc(sizeof(*poly1) * MAX_POLYGONS)) ) {
    fprintf(stderr, "malloc(poly1) fails: %s\n", strerror(errno));
    return 1;
  }

  if ( !(poly2 = malloc(sizeof(*poly2) * MAX_POLYGONS)) ) {
    fprintf(stderr, "malloc(poly2) fails: %s\n", strerror(errno));
    return 1;
  }

  //
  // Load data from files
  //

  if ( (N1 = load_points(first_input_file_name, x1, y1, MAX_POINTS)) < 4 ) {
    fprintf(stderr, "load_points(%s) fails: %s\n", first_input_file_name, strerror(errno));
    return 1;
  }

  if ( (N2 = load_points(second_input_file_name, x2, y2, MAX_POINTS)) < 4 ) {
    fprintf(stderr, "load_points(%s) fails: %s\n", second_input_file_name, strerror(errno));
    return 1;
  }


  if ( X0 != 0 || Y0 != 0 || S != 1 || Alpha != 0 ) {
    affine_transform(x2, y2, N2, X0, Y0, S, Alpha);
  }





  //
  // Dump points
  //

  printf("Initial %s:\n", first_input_file_name);
  dump_points_to_file(stdout, x1, y1, N1);

  printf("Initial %s:\n", first_input_file_name);
  dump_points_to_file(stdout, x2, y2, N2);

  //
  // create two sets of polygons
  //

  if ( !(n1 = create_polygons(x1, y1, N1, poly1, MAX_POLYGONS)) ) {
    fprintf(stderr, "create_polygons(%s) fails: %s\n", first_input_file_name, strerror(errno));
    return 1;
  }

  if ( !(n2 = create_polygons(x2, y2, N2, poly2, MAX_POLYGONS)) ) {
    fprintf(stderr, "create_polygons(%s) fails: %s\n", first_input_file_name, strerror(errno));
    return 1;
  }


  //
  // Dump polygons
  //
  printf("%s polygons:\n", first_input_file_name );
  dump_polygons_to_file(stdout, poly1, n1);

  printf("%s polygons:\n", second_input_file_name);
  dump_polygons_to_file(stdout, poly2, n2);


  //
  // Match polygons
  //
  maxmatch = n1 < n2 ? n1 : n2;

  if ( !(m1 = malloc(maxmatch * sizeof(*m1))) ) {
    fprintf(stderr, "malloc(m1) fails: %s\n", strerror(errno));
    return 1;
  }

  if ( !(m2 = malloc(maxmatch * sizeof(*m2))) ) {
    fprintf(stderr, "malloc(m2) fails: %s\n", strerror(errno));
    return 1;
  }


  nmatch = match_polygons(tol, poly1, n1, poly2, n2, m1, m2);
  if ( nmatch < 1 ) {
    fprintf(stderr, "NO MATCHES FOUND\n");
    return 0;
  }



  //
  // Dump matches
  //
  printf("MATCHES:\n");
  dump_matches_to_file(stderr, x1, y1, N1, x2, y2, N2,
      poly1, n1, poly2, n2,
      m1, m2, nmatch);

  return 0;
}

