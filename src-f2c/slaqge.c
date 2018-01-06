#line 1 "slaqge.f"
/* slaqge.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "slaqge.f"
/* > \brief \b SLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sg
eequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqge.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqge.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqge.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED */
/*       INTEGER            LDA, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQGE equilibrates a general M by N matrix A using the row and */
/* > column scaling factors in the vectors R and C. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M by N matrix A. */
/* >          On exit, the equilibrated matrix.  See EQUED for the form of */
/* >          the equilibrated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[in] R */
/* > \verbatim */
/* >          R is REAL array, dimension (M) */
/* >          The row scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >          The column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] ROWCND */
/* > \verbatim */
/* >          ROWCND is REAL */
/* >          Ratio of the smallest R(i) to the largest R(i). */
/* > \endverbatim */
/* > */
/* > \param[in] COLCND */
/* > \verbatim */
/* >          COLCND is REAL */
/* >          Ratio of the smallest C(i) to the largest C(i). */
/* > \endverbatim */
/* > */
/* > \param[in] AMAX */
/* > \verbatim */
/* >          AMAX is REAL */
/* >          Absolute value of largest matrix entry. */
/* > \endverbatim */
/* > */
/* > \param[out] EQUED */
/* > \verbatim */
/* >          EQUED is CHARACTER*1 */
/* >          Specifies the form of equilibration that was done. */
/* >          = 'N':  No equilibration */
/* >          = 'R':  Row equilibration, i.e., A has been premultiplied by */
/* >                  diag(R). */
/* >          = 'C':  Column equilibration, i.e., A has been postmultiplied */
/* >                  by diag(C). */
/* >          = 'B':  Both row and column equilibration, i.e., A has been */
/* >                  replaced by diag(R) * A * diag(C). */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  THRESH is a threshold value used to decide if row or column scaling */
/* >  should be done based on the ratio of the row or column scaling */
/* >  factors.  If ROWCND < THRESH, row scaling is done, and if */
/* >  COLCND < THRESH, column scaling is done. */
/* > */
/* >  LARGE and SMALL are threshold values used to decide if row scaling */
/* >  should be done based on the absolute size of the largest matrix */
/* >  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGEauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaqge_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, char *equed, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal cj, large, small;
    extern doublereal slamch_(char *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 177 "slaqge.f"
    /* Parameter adjustments */
#line 177 "slaqge.f"
    a_dim1 = *lda;
#line 177 "slaqge.f"
    a_offset = 1 + a_dim1;
#line 177 "slaqge.f"
    a -= a_offset;
#line 177 "slaqge.f"
    --r__;
#line 177 "slaqge.f"
    --c__;
#line 177 "slaqge.f"

#line 177 "slaqge.f"
    /* Function Body */
#line 177 "slaqge.f"
    if (*m <= 0 || *n <= 0) {
#line 178 "slaqge.f"
	*(unsigned char *)equed = 'N';
#line 179 "slaqge.f"
	return 0;
#line 180 "slaqge.f"
    }

/*     Initialize LARGE and SMALL. */

#line 184 "slaqge.f"
    small = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 185 "slaqge.f"
    large = 1. / small;

#line 187 "slaqge.f"
    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

#line 192 "slaqge.f"
	if (*colcnd >= .1) {

/*           No column scaling */

#line 196 "slaqge.f"
	    *(unsigned char *)equed = 'N';
#line 197 "slaqge.f"
	} else {

/*           Column scaling */

#line 201 "slaqge.f"
	    i__1 = *n;
#line 201 "slaqge.f"
	    for (j = 1; j <= i__1; ++j) {
#line 202 "slaqge.f"
		cj = c__[j];
#line 203 "slaqge.f"
		i__2 = *m;
#line 203 "slaqge.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 204 "slaqge.f"
		    a[i__ + j * a_dim1] = cj * a[i__ + j * a_dim1];
#line 205 "slaqge.f"
/* L10: */
#line 205 "slaqge.f"
		}
#line 206 "slaqge.f"
/* L20: */
#line 206 "slaqge.f"
	    }
#line 207 "slaqge.f"
	    *(unsigned char *)equed = 'C';
#line 208 "slaqge.f"
	}
#line 209 "slaqge.f"
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

#line 213 "slaqge.f"
	i__1 = *n;
#line 213 "slaqge.f"
	for (j = 1; j <= i__1; ++j) {
#line 214 "slaqge.f"
	    i__2 = *m;
#line 214 "slaqge.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "slaqge.f"
		a[i__ + j * a_dim1] = r__[i__] * a[i__ + j * a_dim1];
#line 216 "slaqge.f"
/* L30: */
#line 216 "slaqge.f"
	    }
#line 217 "slaqge.f"
/* L40: */
#line 217 "slaqge.f"
	}
#line 218 "slaqge.f"
	*(unsigned char *)equed = 'R';
#line 219 "slaqge.f"
    } else {

/*        Row and column scaling */

#line 223 "slaqge.f"
	i__1 = *n;
#line 223 "slaqge.f"
	for (j = 1; j <= i__1; ++j) {
#line 224 "slaqge.f"
	    cj = c__[j];
#line 225 "slaqge.f"
	    i__2 = *m;
#line 225 "slaqge.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 226 "slaqge.f"
		a[i__ + j * a_dim1] = cj * r__[i__] * a[i__ + j * a_dim1];
#line 227 "slaqge.f"
/* L50: */
#line 227 "slaqge.f"
	    }
#line 228 "slaqge.f"
/* L60: */
#line 228 "slaqge.f"
	}
#line 229 "slaqge.f"
	*(unsigned char *)equed = 'B';
#line 230 "slaqge.f"
    }

#line 232 "slaqge.f"
    return 0;

/*     End of SLAQGE */

} /* slaqge_ */

