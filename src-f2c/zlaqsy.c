#line 1 "zlaqsy.f"
/* zlaqsy.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqsy.f"
/* > \brief \b ZLAQSY scales a symmetric/Hermitian matrix, using scaling factors computed by spoequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqsy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqsy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqsy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            LDA, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   S( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQSY equilibrates a symmetric matrix A using the scaling factors */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          S is DOUBLE PRECISION array, dimension (N) */
/* >          The scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] SCOND */
/* > \verbatim */
/* >          SCOND is DOUBLE PRECISION */
/* >          Ratio of the smallest S(i) to the largest S(i). */
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

/* > \ingroup complex16SYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaqsy_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublereal cj, large;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small;
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
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 171 "zlaqsy.f"
    /* Parameter adjustments */
#line 171 "zlaqsy.f"
    a_dim1 = *lda;
#line 171 "zlaqsy.f"
    a_offset = 1 + a_dim1;
#line 171 "zlaqsy.f"
    a -= a_offset;
#line 171 "zlaqsy.f"
    --s;
#line 171 "zlaqsy.f"

#line 171 "zlaqsy.f"
    /* Function Body */
#line 171 "zlaqsy.f"
    if (*n <= 0) {
#line 172 "zlaqsy.f"
	*(unsigned char *)equed = 'N';
#line 173 "zlaqsy.f"
	return 0;
#line 174 "zlaqsy.f"
    }

/*     Initialize LARGE and SMALL. */

#line 178 "zlaqsy.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 179 "zlaqsy.f"
    large = 1. / small;

#line 181 "zlaqsy.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 185 "zlaqsy.f"
	*(unsigned char *)equed = 'N';
#line 186 "zlaqsy.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 190 "zlaqsy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored. */

#line 194 "zlaqsy.f"
	    i__1 = *n;
#line 194 "zlaqsy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 195 "zlaqsy.f"
		cj = s[j];
#line 196 "zlaqsy.f"
		i__2 = j;
#line 196 "zlaqsy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 197 "zlaqsy.f"
		    i__3 = i__ + j * a_dim1;
#line 197 "zlaqsy.f"
		    d__1 = cj * s[i__];
#line 197 "zlaqsy.f"
		    i__4 = i__ + j * a_dim1;
#line 197 "zlaqsy.f"
		    z__1.r = d__1 * a[i__4].r, z__1.i = d__1 * a[i__4].i;
#line 197 "zlaqsy.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 198 "zlaqsy.f"
/* L10: */
#line 198 "zlaqsy.f"
		}
#line 199 "zlaqsy.f"
/* L20: */
#line 199 "zlaqsy.f"
	    }
#line 200 "zlaqsy.f"
	} else {

/*           Lower triangle of A is stored. */

#line 204 "zlaqsy.f"
	    i__1 = *n;
#line 204 "zlaqsy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "zlaqsy.f"
		cj = s[j];
#line 206 "zlaqsy.f"
		i__2 = *n;
#line 206 "zlaqsy.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 207 "zlaqsy.f"
		    i__3 = i__ + j * a_dim1;
#line 207 "zlaqsy.f"
		    d__1 = cj * s[i__];
#line 207 "zlaqsy.f"
		    i__4 = i__ + j * a_dim1;
#line 207 "zlaqsy.f"
		    z__1.r = d__1 * a[i__4].r, z__1.i = d__1 * a[i__4].i;
#line 207 "zlaqsy.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 208 "zlaqsy.f"
/* L30: */
#line 208 "zlaqsy.f"
		}
#line 209 "zlaqsy.f"
/* L40: */
#line 209 "zlaqsy.f"
	    }
#line 210 "zlaqsy.f"
	}
#line 211 "zlaqsy.f"
	*(unsigned char *)equed = 'Y';
#line 212 "zlaqsy.f"
    }

#line 214 "zlaqsy.f"
    return 0;

/*     End of ZLAQSY */

} /* zlaqsy_ */

