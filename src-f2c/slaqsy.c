#line 1 "slaqsy.f"
/* slaqsy.f -- translated by f2c (version 20100827).
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

#line 1 "slaqsy.f"
/* > \brief \b SLAQSY scales a symmetric/Hermitian matrix, using scaling factors computed by spoequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqsy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqsy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqsy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            LDA, N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQSY equilibrates a symmetric matrix A using the scaling factors */
/* > in the vector S. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          n by n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if EQUED = 'Y', the equilibrated matrix: */
/* >          diag(S) * A * diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL array, dimension (N) */
/* >          The scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] SCOND */
/* > \verbatim */
/* >          SCOND is REAL */
/* >          Ratio of the smallest S(i) to the largest S(i). */
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
/* >          Specifies whether or not equilibration was done. */
/* >          = 'N':  No equilibration. */
/* >          = 'Y':  Equilibration was done, i.e., A has been replaced by */
/* >                  diag(S) * A * diag(S). */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  THRESH is a threshold value used to decide if scaling should be done */
/* >  based on the ratio of the scaling factors.  If SCOND < THRESH, */
/* >  scaling is done. */
/* > */
/* >  LARGE and SMALL are threshold values used to decide if scaling should */
/* >  be done based on the absolute size of the largest matrix element. */
/* >  If AMAX > LARGE or AMAX < SMALL, scaling is done. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realSYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaqsy_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed, 
	ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal cj, large;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small;
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
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 169 "slaqsy.f"
    /* Parameter adjustments */
#line 169 "slaqsy.f"
    a_dim1 = *lda;
#line 169 "slaqsy.f"
    a_offset = 1 + a_dim1;
#line 169 "slaqsy.f"
    a -= a_offset;
#line 169 "slaqsy.f"
    --s;
#line 169 "slaqsy.f"

#line 169 "slaqsy.f"
    /* Function Body */
#line 169 "slaqsy.f"
    if (*n <= 0) {
#line 170 "slaqsy.f"
	*(unsigned char *)equed = 'N';
#line 171 "slaqsy.f"
	return 0;
#line 172 "slaqsy.f"
    }

/*     Initialize LARGE and SMALL. */

#line 176 "slaqsy.f"
    small = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 177 "slaqsy.f"
    large = 1. / small;

#line 179 "slaqsy.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 183 "slaqsy.f"
	*(unsigned char *)equed = 'N';
#line 184 "slaqsy.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 188 "slaqsy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored. */

#line 192 "slaqsy.f"
	    i__1 = *n;
#line 192 "slaqsy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 193 "slaqsy.f"
		cj = s[j];
#line 194 "slaqsy.f"
		i__2 = j;
#line 194 "slaqsy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 195 "slaqsy.f"
		    a[i__ + j * a_dim1] = cj * s[i__] * a[i__ + j * a_dim1];
#line 196 "slaqsy.f"
/* L10: */
#line 196 "slaqsy.f"
		}
#line 197 "slaqsy.f"
/* L20: */
#line 197 "slaqsy.f"
	    }
#line 198 "slaqsy.f"
	} else {

/*           Lower triangle of A is stored. */

#line 202 "slaqsy.f"
	    i__1 = *n;
#line 202 "slaqsy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 203 "slaqsy.f"
		cj = s[j];
#line 204 "slaqsy.f"
		i__2 = *n;
#line 204 "slaqsy.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 205 "slaqsy.f"
		    a[i__ + j * a_dim1] = cj * s[i__] * a[i__ + j * a_dim1];
#line 206 "slaqsy.f"
/* L30: */
#line 206 "slaqsy.f"
		}
#line 207 "slaqsy.f"
/* L40: */
#line 207 "slaqsy.f"
	    }
#line 208 "slaqsy.f"
	}
#line 209 "slaqsy.f"
	*(unsigned char *)equed = 'Y';
#line 210 "slaqsy.f"
    }

#line 212 "slaqsy.f"
    return 0;

/*     End of SLAQSY */

} /* slaqsy_ */

