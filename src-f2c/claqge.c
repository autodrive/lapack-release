#line 1 "claqge.f"
/* claqge.f -- translated by f2c (version 20100827).
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

#line 1 "claqge.f"
/* > \brief \b CLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sg
eequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAQGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqge.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqge.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqge.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED */
/*       INTEGER            LDA, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ), R( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQGE equilibrates a general M by N matrix A using the row and */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \ingroup complexGEauxiliary */

/*  ===================================================================== */
/* Subroutine */ int claqge_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, char *equed, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1;

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

#line 179 "claqge.f"
    /* Parameter adjustments */
#line 179 "claqge.f"
    a_dim1 = *lda;
#line 179 "claqge.f"
    a_offset = 1 + a_dim1;
#line 179 "claqge.f"
    a -= a_offset;
#line 179 "claqge.f"
    --r__;
#line 179 "claqge.f"
    --c__;
#line 179 "claqge.f"

#line 179 "claqge.f"
    /* Function Body */
#line 179 "claqge.f"
    if (*m <= 0 || *n <= 0) {
#line 180 "claqge.f"
	*(unsigned char *)equed = 'N';
#line 181 "claqge.f"
	return 0;
#line 182 "claqge.f"
    }

/*     Initialize LARGE and SMALL. */

#line 186 "claqge.f"
    small = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 187 "claqge.f"
    large = 1. / small;

#line 189 "claqge.f"
    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

#line 194 "claqge.f"
	if (*colcnd >= .1) {

/*           No column scaling */

#line 198 "claqge.f"
	    *(unsigned char *)equed = 'N';
#line 199 "claqge.f"
	} else {

/*           Column scaling */

#line 203 "claqge.f"
	    i__1 = *n;
#line 203 "claqge.f"
	    for (j = 1; j <= i__1; ++j) {
#line 204 "claqge.f"
		cj = c__[j];
#line 205 "claqge.f"
		i__2 = *m;
#line 205 "claqge.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 206 "claqge.f"
		    i__3 = i__ + j * a_dim1;
#line 206 "claqge.f"
		    i__4 = i__ + j * a_dim1;
#line 206 "claqge.f"
		    z__1.r = cj * a[i__4].r, z__1.i = cj * a[i__4].i;
#line 206 "claqge.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 207 "claqge.f"
/* L10: */
#line 207 "claqge.f"
		}
#line 208 "claqge.f"
/* L20: */
#line 208 "claqge.f"
	    }
#line 209 "claqge.f"
	    *(unsigned char *)equed = 'C';
#line 210 "claqge.f"
	}
#line 211 "claqge.f"
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

#line 215 "claqge.f"
	i__1 = *n;
#line 215 "claqge.f"
	for (j = 1; j <= i__1; ++j) {
#line 216 "claqge.f"
	    i__2 = *m;
#line 216 "claqge.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 217 "claqge.f"
		i__3 = i__ + j * a_dim1;
#line 217 "claqge.f"
		i__4 = i__;
#line 217 "claqge.f"
		i__5 = i__ + j * a_dim1;
#line 217 "claqge.f"
		z__1.r = r__[i__4] * a[i__5].r, z__1.i = r__[i__4] * a[i__5]
			.i;
#line 217 "claqge.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 218 "claqge.f"
/* L30: */
#line 218 "claqge.f"
	    }
#line 219 "claqge.f"
/* L40: */
#line 219 "claqge.f"
	}
#line 220 "claqge.f"
	*(unsigned char *)equed = 'R';
#line 221 "claqge.f"
    } else {

/*        Row and column scaling */

#line 225 "claqge.f"
	i__1 = *n;
#line 225 "claqge.f"
	for (j = 1; j <= i__1; ++j) {
#line 226 "claqge.f"
	    cj = c__[j];
#line 227 "claqge.f"
	    i__2 = *m;
#line 227 "claqge.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 228 "claqge.f"
		i__3 = i__ + j * a_dim1;
#line 228 "claqge.f"
		d__1 = cj * r__[i__];
#line 228 "claqge.f"
		i__4 = i__ + j * a_dim1;
#line 228 "claqge.f"
		z__1.r = d__1 * a[i__4].r, z__1.i = d__1 * a[i__4].i;
#line 228 "claqge.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 229 "claqge.f"
/* L50: */
#line 229 "claqge.f"
	    }
#line 230 "claqge.f"
/* L60: */
#line 230 "claqge.f"
	}
#line 231 "claqge.f"
	*(unsigned char *)equed = 'B';
#line 232 "claqge.f"
    }

#line 234 "claqge.f"
    return 0;

/*     End of CLAQGE */

} /* claqge_ */

