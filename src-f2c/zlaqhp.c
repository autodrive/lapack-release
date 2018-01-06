#line 1 "zlaqhp.f"
/* zlaqhp.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqhp.f"
/* > \brief \b ZLAQHP scales a Hermitian matrix stored in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQHP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqhp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqhp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqhp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   S( * ) */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQHP equilibrates a Hermitian matrix A using the scaling factors */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in */
/* >          the same storage format as A. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaqhp_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed, 
	ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, jc;
    static doublereal cj, large;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small;
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 166 "zlaqhp.f"
    /* Parameter adjustments */
#line 166 "zlaqhp.f"
    --s;
#line 166 "zlaqhp.f"
    --ap;
#line 166 "zlaqhp.f"

#line 166 "zlaqhp.f"
    /* Function Body */
#line 166 "zlaqhp.f"
    if (*n <= 0) {
#line 167 "zlaqhp.f"
	*(unsigned char *)equed = 'N';
#line 168 "zlaqhp.f"
	return 0;
#line 169 "zlaqhp.f"
    }

/*     Initialize LARGE and SMALL. */

#line 173 "zlaqhp.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 174 "zlaqhp.f"
    large = 1. / small;

#line 176 "zlaqhp.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 180 "zlaqhp.f"
	*(unsigned char *)equed = 'N';
#line 181 "zlaqhp.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 185 "zlaqhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored. */

#line 189 "zlaqhp.f"
	    jc = 1;
#line 190 "zlaqhp.f"
	    i__1 = *n;
#line 190 "zlaqhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 191 "zlaqhp.f"
		cj = s[j];
#line 192 "zlaqhp.f"
		i__2 = j - 1;
#line 192 "zlaqhp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 193 "zlaqhp.f"
		    i__3 = jc + i__ - 1;
#line 193 "zlaqhp.f"
		    d__1 = cj * s[i__];
#line 193 "zlaqhp.f"
		    i__4 = jc + i__ - 1;
#line 193 "zlaqhp.f"
		    z__1.r = d__1 * ap[i__4].r, z__1.i = d__1 * ap[i__4].i;
#line 193 "zlaqhp.f"
		    ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 194 "zlaqhp.f"
/* L10: */
#line 194 "zlaqhp.f"
		}
#line 195 "zlaqhp.f"
		i__2 = jc + j - 1;
#line 195 "zlaqhp.f"
		i__3 = jc + j - 1;
#line 195 "zlaqhp.f"
		d__1 = cj * cj * ap[i__3].r;
#line 195 "zlaqhp.f"
		ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 196 "zlaqhp.f"
		jc += j;
#line 197 "zlaqhp.f"
/* L20: */
#line 197 "zlaqhp.f"
	    }
#line 198 "zlaqhp.f"
	} else {

/*           Lower triangle of A is stored. */

#line 202 "zlaqhp.f"
	    jc = 1;
#line 203 "zlaqhp.f"
	    i__1 = *n;
#line 203 "zlaqhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 204 "zlaqhp.f"
		cj = s[j];
#line 205 "zlaqhp.f"
		i__2 = jc;
#line 205 "zlaqhp.f"
		i__3 = jc;
#line 205 "zlaqhp.f"
		d__1 = cj * cj * ap[i__3].r;
#line 205 "zlaqhp.f"
		ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 206 "zlaqhp.f"
		i__2 = *n;
#line 206 "zlaqhp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 207 "zlaqhp.f"
		    i__3 = jc + i__ - j;
#line 207 "zlaqhp.f"
		    d__1 = cj * s[i__];
#line 207 "zlaqhp.f"
		    i__4 = jc + i__ - j;
#line 207 "zlaqhp.f"
		    z__1.r = d__1 * ap[i__4].r, z__1.i = d__1 * ap[i__4].i;
#line 207 "zlaqhp.f"
		    ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 208 "zlaqhp.f"
/* L30: */
#line 208 "zlaqhp.f"
		}
#line 209 "zlaqhp.f"
		jc = jc + *n - j + 1;
#line 210 "zlaqhp.f"
/* L40: */
#line 210 "zlaqhp.f"
	    }
#line 211 "zlaqhp.f"
	}
#line 212 "zlaqhp.f"
	*(unsigned char *)equed = 'Y';
#line 213 "zlaqhp.f"
    }

#line 215 "zlaqhp.f"
    return 0;

/*     End of ZLAQHP */

} /* zlaqhp_ */

