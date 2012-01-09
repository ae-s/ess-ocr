/* gcc deskew.c -o deskew -llept -lm -ljpeg -lgif -ltiff -lpng */

/*====================================================================*
 -  Copyright (C) 2005 La Monte H.P. Yarroll
 -  Copyright (C) 2001 Leptonica.  All rights reserved.
 -  This software is distributed in the hope that it will be
 -  useful, but with NO WARRANTY OF ANY KIND.
 -  No author or distributor accepts responsibility to anyone for the
 -  consequences of using this software, or for whether it serves any
 -  particular purpose or works at all, unless he or she says so in
 -  writing.  Everyone is granted permission to copy, modify and
 -  redistribute this source code, for commercial or non-commercial
 -  purposes, with the following restrictions: (1) the origin of this
 -  source code must not be misrepresented; (2) modified versions must
 -  be plainly marked as such; and (3) this notice may not be removed
 -  or altered from any source or modified source distribution.
 *====================================================================*/


/*
 * deskew.c
 *
 *     Deskew the input (binary) image using the best method in Leptonica.
 *     The best version does a linear sweep followed by a binary
 *     (angle-splitting) search.  The basic method is to find the vertical
 *     shear angle such that the differential variance of ON pixels between
 *     each line and it's neighbor, when summed over all lines, is maximized.
 *     Use the best area mapping rotation to produce a new image.
 */

#include <stdio.h>
#include <stdlib.h>
#include <leptonica/allheaders.h>

    /* Default sweep angle parameters for pixFindSkew() */
static const l_float32  DEFAULT_SWEEP_RANGE = 5.;    /* degrees */
static const l_float32  DEFAULT_SWEEP_DELTA = 1.;    /* degrees */

    /* Default final angle difference parameter for binary
     * search in pixFindSkew().  The expected accuracy is
     * not better than the inverse image width in pixels,
     * say, 1/2000 radians, or about 0.03 degrees. */
static const l_float32  DEFAULT_MINBS_DELTA = 0.01;  /* degrees */

    /* Default scale factors for pixFindSkew() */
static const l_int32  DEFAULT_SWEEP_REDUCTION = 4;  /* sweep part; 4 is good */
static const l_int32  DEFAULT_BS_REDUCTION = 2;  /* binary search part */

    /* Minimum angle for deskewing in pixDeskew() */
static const l_float32  MIN_DESKEW_ANGLE = 0.1;  /* degree */

    /* Minimum allowed confidence (ratio) for deskewing in pixDeskew() */
static const l_float32  MIN_ALLOWED_CONFIDENCE = 3.0;

    /* Minimum allowed maxscore to give nonzero confidence */
static const l_int32  MIN_VALID_MAXSCORE = 10000;

    /* Constant setting threshold for minimum allowed minscore
     * to give nonzero confidence; multiply this constant by
     *  (height * width^2) */
static const l_float32  MINSCORE_THRESHOLD_CONSTANT = 0.000002;


#ifndef  NO_CONSOLE_IO
#define  DEBUG_PRINT_SCORES     0
#define  DEBUG_PRINT_SWEEP      0
#define  DEBUG_PRINT_BINARY     0
#define  DEBUG_THRESHOLD        0
#define  DEBUG_PLOT_SCORES      0
#endif  /* ~NO_CONSOLE_IO */


    /* deskew */
#define   DESKEW_REDUCTION      2      /* 1, 2 or 4 */

    /* sweep only */
#define   SWEEP_RANGE           5.     /* degrees */
#define   SWEEP_DELTA           0.2    /* degrees */
#define   SWEEP_REDUCTION       2      /* 1, 2, 4 or 8 */

    /* sweep and search */
#define   SWEEP_RANGE2          5.     /* degrees */
#define   SWEEP_DELTA2          1.     /* degrees */
#define   SWEEP_REDUCTION2      2      /* 1, 2, 4 or 8 */
#define   SEARCH_REDUCTION      2      /* 1, 2, 4 or 8 */
#define   SEARCH_MIN_DELTA      0.01   /* degrees */


/*!
 *  deskew()
 *
 *      Input:  pixs
 *              redsearch  (for binary search: reduction factor = 1, 2 or 4)
 *      Return: deskewed pix, or NULL on error
 */
PIX *
deskew(PIX     *pixs,
       l_int32  redsearch)
{
l_float32  angle, conf, deg2rad;
PIX       *pixg;  /* gray version */
PIX       *pixb; /* binary version */
PIX       *pixd;  /* destination image */

    PROCNAME("deskew");

    if (!pixs)
	return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);

    /* Calculate a skew angle.  We may need to make a binary version of the
     * image for this calculation.
     */
    if (pixGetDepth(pixs) != 1) {
	/* FIX ME:  We should probably pick a threshold value with more care.  */
	/* Create a grayscale image if we need one.  */
	if (pixGetDepth(pixs) >= 24) {
	    pixg = pixConvertRGBToGray(pixs, 0.0, 0.0, 0.0);
	} else {
	    pixg = pixs;
	}
	    
	pixb = pixThresholdToBinary(pixg, 127);
	if (pixg != pixs) {
	    pixDestroy(&pixg);
	}
	/* Assert:  We are done with any gray image.  */
    } else {
	pixb = pixs;
    }
    /* Assert: We have a valid binary image.  */
    if (redsearch != 1 && redsearch != 2 && redsearch != 4)
	return (PIX *)ERROR_PTR("redsearch not in {1,2,4}", procName, NULL);

    deg2rad = 3.1415926535 / 180.;
    if (pixFindSkewSweepAndSearch(pixb, &angle, &conf,
				  DEFAULT_SWEEP_REDUCTION, redsearch,
				  DEFAULT_SWEEP_RANGE, DEFAULT_SWEEP_DELTA,
				  DEFAULT_MINBS_DELTA)) {
	pixd = pixClone(pixs);
	goto finish;
    }
	
    if (L_ABS(angle) < MIN_DESKEW_ANGLE || conf < MIN_ALLOWED_CONFIDENCE) {
	pixd = pixClone(pixs);
	goto finish;
    }

    /* If the pixel depth of pixs is 1, we need to use a bit-depth
     * independent rotate instead of the more accurate area mapping rotate.
     */
    if (pixGetDepth(pixs) == 1) {
	if ((pixd = pixRotateShear(pixs, 0, 0, deg2rad * angle, 0xffffff00)) == NULL) {
	    pixd = pixClone(pixs);
	}
    } else {
#if defined(COLOR_ROTATE)
	if ((pixd = pixRotateAMColorFast(pixs, deg2rad * angle)) == NULL) {
	    pixd = pixClone(pixs);
	}
#else
	if ((pixd = pixRotateAM(pixs, deg2rad * angle, 0xffffff00)) == NULL) {
	    pixd = pixClone(pixs);
	}
#endif
    }

   finish:
    if (pixb != pixs) {
	pixDestroy(&pixb);
    }
    return pixd;
}

main(int    argc,
     char **argv)
{
char        *filein, *fileout;
l_float32    deg2rad;
l_float32    angle, conf;
PIX         *pixs, *pixd;
static char  mainName[] = "deskew";

    if (argc != 3)
	exit(ERROR_INT(" Syntax:  deskew filein fileout", mainName, 1));

    filein = argv[1];
    fileout = argv[2];

    deg2rad = 3.1415926535 / 180.;

    if ((pixs = pixRead(filein)) == NULL)
	exit(ERROR_INT("pixs not made", mainName, 1));

    pixd = deskew(pixs, DESKEW_REDUCTION);
    pixWrite(fileout, pixd, IFF_PNG);
    pixDestroy(&pixd);
    pixDestroy(&pixs);
    exit(0);
}
