#line 1 "zlaqhe.f"
/* zlaqhe.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqhe.f"
/* > \brief \b ZLAQHE scales a Hermitian matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQHE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqhe.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqhe.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqhe.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) */

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
/* > ZLAQHE equilibrates a Hermitian matrix A using the scaling factors */
/* > in the vector S. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          Hermitian matrix A is stored. */
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
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
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

/* > \ingroup complex16HEauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaqhe_(char *uplo, integer *n, doublecomplex *a, 
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 174 "zlaqhe.f"
    /* Parameter adjustments */
#line 174 "zlaqhe.f"
    a_dim1 = *lda;
#line 174 "zlaqhe.f"
    a_offset = 1 + a_dim1;
#line 174 "zlaqhe.f"
    a -= a_offset;
#line 174 "zlaqhe.f"
    --s;
#line 174 "zlaqhe.f"

#line 174 "zlaqhe.f"
    /* Function Body */
#line 174 "zlaqhe.f"
    if (*n <= 0) {
#line 175 "zlaqhe.f"
	*(unsigned char *)equed = 'N';
#line 176 "zlaqhe.f"
	return 0;
#line 177 "zlaqhe.f"
    }

/*     Initialize LARGE and SMALL. */

#line 181 "zlaqhe.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 182 "zlaqhe.f"
    large = 1. / small;

#line 184 "zlaqhe.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 188 "zlaqhe.f"
	*(unsigned char *)equed = 'N';
#line 189 "zlaqhe.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 193 "zlaqhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored. */

#line 197 "zlaqhe.f"
	    i__1 = *n;
#line 197 "zlaqhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 198 "zlaqhe.f"
		cj = s[j];
#line 199 "zlaqhe.f"
		i__2 = j - 1;
#line 199 "zlaqhe.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 200 "zlaqhe.f"
		    i__3 = i__ + j * a_dim1;
#line 200 "zlaqhe.f"
		    d__1 = cj * s[i__];
#line 200 "zlaqhe.f"
		    i__4 = i__ + j * a_dim1;
#line 200 "zlaqhe.f"
		    z__1.r = d__1 * a[i__4].r, z__1.i = d__1 * a[i__4].i;
#line 200 "zlaqhe.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 201 "zlaqhe.f"
/* L10: */
#line 201 "zlaqhe.f"
		}
#line 202 "zlaqhe.f"
		i__2 = j + j * a_dim1;
#line 202 "zlaqhe.f"
		i__3 = j + j * a_dim1;
#line 202 "zlaqhe.f"
		d__1 = cj * cj * a[i__3].r;
#line 202 "zlaqhe.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 203 "zlaqhe.f"
/* L20: */
#line 203 "zlaqhe.f"
	    }
#line 204 "zlaqhe.f"
	} else {

/*           Lower triangle of A is stored. */

#line 208 "zlaqhe.f"
	    i__1 = *n;
#line 208 "zlaqhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 209 "zlaqhe.f"
		cj = s[j];
#line 210 "zlaqhe.f"
		i__2 = j + j * a_dim1;
#line 210 "zlaqhe.f"
		i__3 = j + j * a_dim1;
#line 210 "zlaqhe.f"
		d__1 = cj * cj * a[i__3].r;
#line 210 "zlaqhe.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 211 "zlaqhe.f"
		i__2 = *n;
#line 211 "zlaqhe.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 212 "zlaqhe.f"
		    i__3 = i__ + j * a_dim1;
#line 212 "zlaqhe.f"
		    d__1 = cj * s[i__];
#line 212 "zlaqhe.f"
		    i__4 = i__ + j * a_dim1;
#line 212 "zlaqhe.f"
		    z__1.r = d__1 * a[i__4].r, z__1.i = d__1 * a[i__4].i;
#line 212 "zlaqhe.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 213 "zlaqhe.f"
/* L30: */
#line 213 "zlaqhe.f"
		}
#line 214 "zlaqhe.f"
/* L40: */
#line 214 "zlaqhe.f"
	    }
#line 215 "zlaqhe.f"
	}
#line 216 "zlaqhe.f"
	*(unsigned char *)equed = 'Y';
#line 217 "zlaqhe.f"
    }

#line 219 "zlaqhe.f"
    return 0;

/*     End of ZLAQHE */

} /* zlaqhe_ */

