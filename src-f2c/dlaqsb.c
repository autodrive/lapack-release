#line 1 "dlaqsb.f"
/* dlaqsb.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqsb.f"
/* > \brief \b DLAQSB scales a symmetric/Hermitian band matrix, using scaling factors computed by spbequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQSB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqsb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqsb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqsb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQSB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            KD, LDAB, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQSB equilibrates a symmetric band matrix A using the scaling */
/* > factors in the vector S. */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**T*U or A = L*L**T of the band */
/* >          matrix A, in the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaqsb_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	 char *equed, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j;
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

#line 179 "dlaqsb.f"
    /* Parameter adjustments */
#line 179 "dlaqsb.f"
    ab_dim1 = *ldab;
#line 179 "dlaqsb.f"
    ab_offset = 1 + ab_dim1;
#line 179 "dlaqsb.f"
    ab -= ab_offset;
#line 179 "dlaqsb.f"
    --s;
#line 179 "dlaqsb.f"

#line 179 "dlaqsb.f"
    /* Function Body */
#line 179 "dlaqsb.f"
    if (*n <= 0) {
#line 180 "dlaqsb.f"
	*(unsigned char *)equed = 'N';
#line 181 "dlaqsb.f"
	return 0;
#line 182 "dlaqsb.f"
    }

/*     Initialize LARGE and SMALL. */

#line 186 "dlaqsb.f"
    small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 187 "dlaqsb.f"
    large = 1. / small;

#line 189 "dlaqsb.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 193 "dlaqsb.f"
	*(unsigned char *)equed = 'N';
#line 194 "dlaqsb.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 198 "dlaqsb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored in band format. */

#line 202 "dlaqsb.f"
	    i__1 = *n;
#line 202 "dlaqsb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 203 "dlaqsb.f"
		cj = s[j];
/* Computing MAX */
#line 204 "dlaqsb.f"
		i__2 = 1, i__3 = j - *kd;
#line 204 "dlaqsb.f"
		i__4 = j;
#line 204 "dlaqsb.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 205 "dlaqsb.f"
		    ab[*kd + 1 + i__ - j + j * ab_dim1] = cj * s[i__] * ab[*
			    kd + 1 + i__ - j + j * ab_dim1];
#line 206 "dlaqsb.f"
/* L10: */
#line 206 "dlaqsb.f"
		}
#line 207 "dlaqsb.f"
/* L20: */
#line 207 "dlaqsb.f"
	    }
#line 208 "dlaqsb.f"
	} else {

/*           Lower triangle of A is stored. */

#line 212 "dlaqsb.f"
	    i__1 = *n;
#line 212 "dlaqsb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 213 "dlaqsb.f"
		cj = s[j];
/* Computing MIN */
#line 214 "dlaqsb.f"
		i__2 = *n, i__3 = j + *kd;
#line 214 "dlaqsb.f"
		i__4 = min(i__2,i__3);
#line 214 "dlaqsb.f"
		for (i__ = j; i__ <= i__4; ++i__) {
#line 215 "dlaqsb.f"
		    ab[i__ + 1 - j + j * ab_dim1] = cj * s[i__] * ab[i__ + 1 
			    - j + j * ab_dim1];
#line 216 "dlaqsb.f"
/* L30: */
#line 216 "dlaqsb.f"
		}
#line 217 "dlaqsb.f"
/* L40: */
#line 217 "dlaqsb.f"
	    }
#line 218 "dlaqsb.f"
	}
#line 219 "dlaqsb.f"
	*(unsigned char *)equed = 'Y';
#line 220 "dlaqsb.f"
    }

#line 222 "dlaqsb.f"
    return 0;

/*     End of DLAQSB */

} /* dlaqsb_ */

