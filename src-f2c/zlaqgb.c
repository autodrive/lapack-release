#line 1 "zlaqgb.f"
/* zlaqgb.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqgb.f"
/* > \brief \b ZLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqgb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqgb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqgb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                          AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED */
/*       INTEGER            KL, KU, LDAB, M, N */
/*       DOUBLE PRECISION   AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), R( * ) */
/*       COMPLEX*16         AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQGB equilibrates a general M by N band matrix A with KL */
/* > subdiagonals and KU superdiagonals using the row and scaling factors */
/* > in the vectors R and C. */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl) */
/* > */
/* >          On exit, the equilibrated matrix, in the same storage format */
/* >          as A.  See EQUED for the form of the equilibrated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDA >= KL+KU+1. */
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

/* > \date September 2012 */

/* > \ingroup complex16GBauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed,
	 ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublereal cj, large, small;
    extern doublereal dlamch_(char *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 199 "zlaqgb.f"
    /* Parameter adjustments */
#line 199 "zlaqgb.f"
    ab_dim1 = *ldab;
#line 199 "zlaqgb.f"
    ab_offset = 1 + ab_dim1;
#line 199 "zlaqgb.f"
    ab -= ab_offset;
#line 199 "zlaqgb.f"
    --r__;
#line 199 "zlaqgb.f"
    --c__;
#line 199 "zlaqgb.f"

#line 199 "zlaqgb.f"
    /* Function Body */
#line 199 "zlaqgb.f"
    if (*m <= 0 || *n <= 0) {
#line 200 "zlaqgb.f"
	*(unsigned char *)equed = 'N';
#line 201 "zlaqgb.f"
	return 0;
#line 202 "zlaqgb.f"
    }

/*     Initialize LARGE and SMALL. */

#line 206 "zlaqgb.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 207 "zlaqgb.f"
    large = 1. / small;

#line 209 "zlaqgb.f"
    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

#line 214 "zlaqgb.f"
	if (*colcnd >= .1) {

/*           No column scaling */

#line 218 "zlaqgb.f"
	    *(unsigned char *)equed = 'N';
#line 219 "zlaqgb.f"
	} else {

/*           Column scaling */

#line 223 "zlaqgb.f"
	    i__1 = *n;
#line 223 "zlaqgb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 224 "zlaqgb.f"
		cj = c__[j];
/* Computing MAX */
#line 225 "zlaqgb.f"
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 225 "zlaqgb.f"
		i__5 = *m, i__6 = j + *kl;
#line 225 "zlaqgb.f"
		i__4 = min(i__5,i__6);
#line 225 "zlaqgb.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 226 "zlaqgb.f"
		    i__2 = *ku + 1 + i__ - j + j * ab_dim1;
#line 226 "zlaqgb.f"
		    i__3 = *ku + 1 + i__ - j + j * ab_dim1;
#line 226 "zlaqgb.f"
		    z__1.r = cj * ab[i__3].r, z__1.i = cj * ab[i__3].i;
#line 226 "zlaqgb.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 227 "zlaqgb.f"
/* L10: */
#line 227 "zlaqgb.f"
		}
#line 228 "zlaqgb.f"
/* L20: */
#line 228 "zlaqgb.f"
	    }
#line 229 "zlaqgb.f"
	    *(unsigned char *)equed = 'C';
#line 230 "zlaqgb.f"
	}
#line 231 "zlaqgb.f"
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

#line 235 "zlaqgb.f"
	i__1 = *n;
#line 235 "zlaqgb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 236 "zlaqgb.f"
	    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 236 "zlaqgb.f"
	    i__5 = *m, i__6 = j + *kl;
#line 236 "zlaqgb.f"
	    i__3 = min(i__5,i__6);
#line 236 "zlaqgb.f"
	    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 237 "zlaqgb.f"
		i__4 = *ku + 1 + i__ - j + j * ab_dim1;
#line 237 "zlaqgb.f"
		i__2 = i__;
#line 237 "zlaqgb.f"
		i__5 = *ku + 1 + i__ - j + j * ab_dim1;
#line 237 "zlaqgb.f"
		z__1.r = r__[i__2] * ab[i__5].r, z__1.i = r__[i__2] * ab[i__5]
			.i;
#line 237 "zlaqgb.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 238 "zlaqgb.f"
/* L30: */
#line 238 "zlaqgb.f"
	    }
#line 239 "zlaqgb.f"
/* L40: */
#line 239 "zlaqgb.f"
	}
#line 240 "zlaqgb.f"
	*(unsigned char *)equed = 'R';
#line 241 "zlaqgb.f"
    } else {

/*        Row and column scaling */

#line 245 "zlaqgb.f"
	i__1 = *n;
#line 245 "zlaqgb.f"
	for (j = 1; j <= i__1; ++j) {
#line 246 "zlaqgb.f"
	    cj = c__[j];
/* Computing MAX */
#line 247 "zlaqgb.f"
	    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 247 "zlaqgb.f"
	    i__5 = *m, i__6 = j + *kl;
#line 247 "zlaqgb.f"
	    i__2 = min(i__5,i__6);
#line 247 "zlaqgb.f"
	    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 248 "zlaqgb.f"
		i__3 = *ku + 1 + i__ - j + j * ab_dim1;
#line 248 "zlaqgb.f"
		d__1 = cj * r__[i__];
#line 248 "zlaqgb.f"
		i__4 = *ku + 1 + i__ - j + j * ab_dim1;
#line 248 "zlaqgb.f"
		z__1.r = d__1 * ab[i__4].r, z__1.i = d__1 * ab[i__4].i;
#line 248 "zlaqgb.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 249 "zlaqgb.f"
/* L50: */
#line 249 "zlaqgb.f"
	    }
#line 250 "zlaqgb.f"
/* L60: */
#line 250 "zlaqgb.f"
	}
#line 251 "zlaqgb.f"
	*(unsigned char *)equed = 'B';
#line 252 "zlaqgb.f"
    }

#line 254 "zlaqgb.f"
    return 0;

/*     End of ZLAQGB */

} /* zlaqgb_ */

