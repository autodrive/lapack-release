#line 1 "dlaqge.f"
/* dlaqge.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqge.f"
/* > \brief \b DLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sg
eequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqge.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqge.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqge.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED */
/*       INTEGER            LDA, M, N */
/*       DOUBLE PRECISION   AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQGE equilibrates a general M by N matrix A using the row and */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          R is DOUBLE PRECISION array, dimension (M) */
/* >          The row scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N) */
/* >          The column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] ROWCND */
/* > \verbatim */
/* >          ROWCND is DOUBLE PRECISION */
/* >          Ratio of the smallest R(i) to the largest R(i). */
/* > \endverbatim */
/* > */
/* > \param[in] COLCND */
/* > \verbatim */
/* >          COLCND is DOUBLE PRECISION */
/* >          Ratio of the smallest C(i) to the largest C(i). */
/* > \endverbatim */
/* > */
/* > \param[in] AMAX */
/* > \verbatim */
/* >          AMAX is DOUBLE PRECISION */
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

/* > \ingroup doubleGEauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaqge_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, char *equed, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal cj, large, small;
    extern doublereal dlamch_(char *, ftnlen);


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

#line 177 "dlaqge.f"
    /* Parameter adjustments */
#line 177 "dlaqge.f"
    a_dim1 = *lda;
#line 177 "dlaqge.f"
    a_offset = 1 + a_dim1;
#line 177 "dlaqge.f"
    a -= a_offset;
#line 177 "dlaqge.f"
    --r__;
#line 177 "dlaqge.f"
    --c__;
#line 177 "dlaqge.f"

#line 177 "dlaqge.f"
    /* Function Body */
#line 177 "dlaqge.f"
    if (*m <= 0 || *n <= 0) {
#line 178 "dlaqge.f"
	*(unsigned char *)equed = 'N';
#line 179 "dlaqge.f"
	return 0;
#line 180 "dlaqge.f"
    }

/*     Initialize LARGE and SMALL. */

#line 184 "dlaqge.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 185 "dlaqge.f"
    large = 1. / small;

#line 187 "dlaqge.f"
    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

#line 192 "dlaqge.f"
	if (*colcnd >= .1) {

/*           No column scaling */

#line 196 "dlaqge.f"
	    *(unsigned char *)equed = 'N';
#line 197 "dlaqge.f"
	} else {

/*           Column scaling */

#line 201 "dlaqge.f"
	    i__1 = *n;
#line 201 "dlaqge.f"
	    for (j = 1; j <= i__1; ++j) {
#line 202 "dlaqge.f"
		cj = c__[j];
#line 203 "dlaqge.f"
		i__2 = *m;
#line 203 "dlaqge.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 204 "dlaqge.f"
		    a[i__ + j * a_dim1] = cj * a[i__ + j * a_dim1];
#line 205 "dlaqge.f"
/* L10: */
#line 205 "dlaqge.f"
		}
#line 206 "dlaqge.f"
/* L20: */
#line 206 "dlaqge.f"
	    }
#line 207 "dlaqge.f"
	    *(unsigned char *)equed = 'C';
#line 208 "dlaqge.f"
	}
#line 209 "dlaqge.f"
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

#line 213 "dlaqge.f"
	i__1 = *n;
#line 213 "dlaqge.f"
	for (j = 1; j <= i__1; ++j) {
#line 214 "dlaqge.f"
	    i__2 = *m;
#line 214 "dlaqge.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "dlaqge.f"
		a[i__ + j * a_dim1] = r__[i__] * a[i__ + j * a_dim1];
#line 216 "dlaqge.f"
/* L30: */
#line 216 "dlaqge.f"
	    }
#line 217 "dlaqge.f"
/* L40: */
#line 217 "dlaqge.f"
	}
#line 218 "dlaqge.f"
	*(unsigned char *)equed = 'R';
#line 219 "dlaqge.f"
    } else {

/*        Row and column scaling */

#line 223 "dlaqge.f"
	i__1 = *n;
#line 223 "dlaqge.f"
	for (j = 1; j <= i__1; ++j) {
#line 224 "dlaqge.f"
	    cj = c__[j];
#line 225 "dlaqge.f"
	    i__2 = *m;
#line 225 "dlaqge.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 226 "dlaqge.f"
		a[i__ + j * a_dim1] = cj * r__[i__] * a[i__ + j * a_dim1];
#line 227 "dlaqge.f"
/* L50: */
#line 227 "dlaqge.f"
	    }
#line 228 "dlaqge.f"
/* L60: */
#line 228 "dlaqge.f"
	}
#line 229 "dlaqge.f"
	*(unsigned char *)equed = 'B';
#line 230 "dlaqge.f"
    }

#line 232 "dlaqge.f"
    return 0;

/*     End of DLAQGE */

} /* dlaqge_ */

