#line 1 "slaqgb.f"
/* slaqgb.f -- translated by f2c (version 20100827).
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

#line 1 "slaqgb.f"
/* > \brief \b SLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqgb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqgb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqgb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                          AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED */
/*       INTEGER            KL, KU, LDAB, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQGB equilibrates a general M by N band matrix A with KL */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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

/* > \date September 2012 */

/* > \ingroup realGBauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed,
	 ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j;
    static doublereal cj, large, small;
    extern doublereal slamch_(char *, ftnlen);


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

#line 197 "slaqgb.f"
    /* Parameter adjustments */
#line 197 "slaqgb.f"
    ab_dim1 = *ldab;
#line 197 "slaqgb.f"
    ab_offset = 1 + ab_dim1;
#line 197 "slaqgb.f"
    ab -= ab_offset;
#line 197 "slaqgb.f"
    --r__;
#line 197 "slaqgb.f"
    --c__;
#line 197 "slaqgb.f"

#line 197 "slaqgb.f"
    /* Function Body */
#line 197 "slaqgb.f"
    if (*m <= 0 || *n <= 0) {
#line 198 "slaqgb.f"
	*(unsigned char *)equed = 'N';
#line 199 "slaqgb.f"
	return 0;
#line 200 "slaqgb.f"
    }

/*     Initialize LARGE and SMALL. */

#line 204 "slaqgb.f"
    small = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 205 "slaqgb.f"
    large = 1. / small;

#line 207 "slaqgb.f"
    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

#line 212 "slaqgb.f"
	if (*colcnd >= .1) {

/*           No column scaling */

#line 216 "slaqgb.f"
	    *(unsigned char *)equed = 'N';
#line 217 "slaqgb.f"
	} else {

/*           Column scaling */

#line 221 "slaqgb.f"
	    i__1 = *n;
#line 221 "slaqgb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 222 "slaqgb.f"
		cj = c__[j];
/* Computing MAX */
#line 223 "slaqgb.f"
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 223 "slaqgb.f"
		i__5 = *m, i__6 = j + *kl;
#line 223 "slaqgb.f"
		i__4 = min(i__5,i__6);
#line 223 "slaqgb.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 224 "slaqgb.f"
		    ab[*ku + 1 + i__ - j + j * ab_dim1] = cj * ab[*ku + 1 + 
			    i__ - j + j * ab_dim1];
#line 225 "slaqgb.f"
/* L10: */
#line 225 "slaqgb.f"
		}
#line 226 "slaqgb.f"
/* L20: */
#line 226 "slaqgb.f"
	    }
#line 227 "slaqgb.f"
	    *(unsigned char *)equed = 'C';
#line 228 "slaqgb.f"
	}
#line 229 "slaqgb.f"
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

#line 233 "slaqgb.f"
	i__1 = *n;
#line 233 "slaqgb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 234 "slaqgb.f"
	    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 234 "slaqgb.f"
	    i__5 = *m, i__6 = j + *kl;
#line 234 "slaqgb.f"
	    i__3 = min(i__5,i__6);
#line 234 "slaqgb.f"
	    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 235 "slaqgb.f"
		ab[*ku + 1 + i__ - j + j * ab_dim1] = r__[i__] * ab[*ku + 1 + 
			i__ - j + j * ab_dim1];
#line 236 "slaqgb.f"
/* L30: */
#line 236 "slaqgb.f"
	    }
#line 237 "slaqgb.f"
/* L40: */
#line 237 "slaqgb.f"
	}
#line 238 "slaqgb.f"
	*(unsigned char *)equed = 'R';
#line 239 "slaqgb.f"
    } else {

/*        Row and column scaling */

#line 243 "slaqgb.f"
	i__1 = *n;
#line 243 "slaqgb.f"
	for (j = 1; j <= i__1; ++j) {
#line 244 "slaqgb.f"
	    cj = c__[j];
/* Computing MAX */
#line 245 "slaqgb.f"
	    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 245 "slaqgb.f"
	    i__5 = *m, i__6 = j + *kl;
#line 245 "slaqgb.f"
	    i__2 = min(i__5,i__6);
#line 245 "slaqgb.f"
	    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 246 "slaqgb.f"
		ab[*ku + 1 + i__ - j + j * ab_dim1] = cj * r__[i__] * ab[*ku 
			+ 1 + i__ - j + j * ab_dim1];
#line 247 "slaqgb.f"
/* L50: */
#line 247 "slaqgb.f"
	    }
#line 248 "slaqgb.f"
/* L60: */
#line 248 "slaqgb.f"
	}
#line 249 "slaqgb.f"
	*(unsigned char *)equed = 'B';
#line 250 "slaqgb.f"
    }

#line 252 "slaqgb.f"
    return 0;

/*     End of SLAQGB */

} /* slaqgb_ */

