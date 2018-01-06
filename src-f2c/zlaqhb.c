#line 1 "zlaqhb.f"
/* zlaqhb.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqhb.f"
/* > \brief \b ZLAQHB scales a Hermitian band matrix, using scaling factors computed by cpbequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQHB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqhb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqhb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqhb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            KD, LDAB, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   S( * ) */
/*       COMPLEX*16         AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAQHB equilibrates a Hermitian band matrix A */
/* > using the scaling factors in the vector S. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of super-diagonals of the matrix A if UPLO = 'U', */
/* >          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**H *U or A = L*L**H of the band */
/* >          matrix A, in the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
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

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaqhb_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
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

#line 181 "zlaqhb.f"
    /* Parameter adjustments */
#line 181 "zlaqhb.f"
    ab_dim1 = *ldab;
#line 181 "zlaqhb.f"
    ab_offset = 1 + ab_dim1;
#line 181 "zlaqhb.f"
    ab -= ab_offset;
#line 181 "zlaqhb.f"
    --s;
#line 181 "zlaqhb.f"

#line 181 "zlaqhb.f"
    /* Function Body */
#line 181 "zlaqhb.f"
    if (*n <= 0) {
#line 182 "zlaqhb.f"
	*(unsigned char *)equed = 'N';
#line 183 "zlaqhb.f"
	return 0;
#line 184 "zlaqhb.f"
    }

/*     Initialize LARGE and SMALL. */

#line 188 "zlaqhb.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 189 "zlaqhb.f"
    large = 1. / small;

#line 191 "zlaqhb.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 195 "zlaqhb.f"
	*(unsigned char *)equed = 'N';
#line 196 "zlaqhb.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 200 "zlaqhb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored in band format. */

#line 204 "zlaqhb.f"
	    i__1 = *n;
#line 204 "zlaqhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "zlaqhb.f"
		cj = s[j];
/* Computing MAX */
#line 206 "zlaqhb.f"
		i__2 = 1, i__3 = j - *kd;
#line 206 "zlaqhb.f"
		i__4 = j - 1;
#line 206 "zlaqhb.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 207 "zlaqhb.f"
		    i__2 = *kd + 1 + i__ - j + j * ab_dim1;
#line 207 "zlaqhb.f"
		    d__1 = cj * s[i__];
#line 207 "zlaqhb.f"
		    i__3 = *kd + 1 + i__ - j + j * ab_dim1;
#line 207 "zlaqhb.f"
		    z__1.r = d__1 * ab[i__3].r, z__1.i = d__1 * ab[i__3].i;
#line 207 "zlaqhb.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 208 "zlaqhb.f"
/* L10: */
#line 208 "zlaqhb.f"
		}
#line 209 "zlaqhb.f"
		i__4 = *kd + 1 + j * ab_dim1;
#line 209 "zlaqhb.f"
		i__2 = *kd + 1 + j * ab_dim1;
#line 209 "zlaqhb.f"
		d__1 = cj * cj * ab[i__2].r;
#line 209 "zlaqhb.f"
		ab[i__4].r = d__1, ab[i__4].i = 0.;
#line 210 "zlaqhb.f"
/* L20: */
#line 210 "zlaqhb.f"
	    }
#line 211 "zlaqhb.f"
	} else {

/*           Lower triangle of A is stored. */

#line 215 "zlaqhb.f"
	    i__1 = *n;
#line 215 "zlaqhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 216 "zlaqhb.f"
		cj = s[j];
#line 217 "zlaqhb.f"
		i__4 = j * ab_dim1 + 1;
#line 217 "zlaqhb.f"
		i__2 = j * ab_dim1 + 1;
#line 217 "zlaqhb.f"
		d__1 = cj * cj * ab[i__2].r;
#line 217 "zlaqhb.f"
		ab[i__4].r = d__1, ab[i__4].i = 0.;
/* Computing MIN */
#line 218 "zlaqhb.f"
		i__2 = *n, i__3 = j + *kd;
#line 218 "zlaqhb.f"
		i__4 = min(i__2,i__3);
#line 218 "zlaqhb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 219 "zlaqhb.f"
		    i__2 = i__ + 1 - j + j * ab_dim1;
#line 219 "zlaqhb.f"
		    d__1 = cj * s[i__];
#line 219 "zlaqhb.f"
		    i__3 = i__ + 1 - j + j * ab_dim1;
#line 219 "zlaqhb.f"
		    z__1.r = d__1 * ab[i__3].r, z__1.i = d__1 * ab[i__3].i;
#line 219 "zlaqhb.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 220 "zlaqhb.f"
/* L30: */
#line 220 "zlaqhb.f"
		}
#line 221 "zlaqhb.f"
/* L40: */
#line 221 "zlaqhb.f"
	    }
#line 222 "zlaqhb.f"
	}
#line 223 "zlaqhb.f"
	*(unsigned char *)equed = 'Y';
#line 224 "zlaqhb.f"
    }

#line 226 "zlaqhb.f"
    return 0;

/*     End of ZLAQHB */

} /* zlaqhb_ */

