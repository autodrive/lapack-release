#line 1 "slaqsp.f"
/* slaqsp.f -- translated by f2c (version 20100827).
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

#line 1 "slaqsp.f"
/* > \brief \b SLAQSP scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by
 sppequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqsp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqsp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqsp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQSP equilibrates a symmetric matrix A using the scaling factors */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaqsp_(char *uplo, integer *n, doublereal *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed, 
	ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jc;
    static doublereal cj, large;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small;
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

#line 161 "slaqsp.f"
    /* Parameter adjustments */
#line 161 "slaqsp.f"
    --s;
#line 161 "slaqsp.f"
    --ap;
#line 161 "slaqsp.f"

#line 161 "slaqsp.f"
    /* Function Body */
#line 161 "slaqsp.f"
    if (*n <= 0) {
#line 162 "slaqsp.f"
	*(unsigned char *)equed = 'N';
#line 163 "slaqsp.f"
	return 0;
#line 164 "slaqsp.f"
    }

/*     Initialize LARGE and SMALL. */

#line 168 "slaqsp.f"
    small = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 169 "slaqsp.f"
    large = 1. / small;

#line 171 "slaqsp.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 175 "slaqsp.f"
	*(unsigned char *)equed = 'N';
#line 176 "slaqsp.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 180 "slaqsp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored. */

#line 184 "slaqsp.f"
	    jc = 1;
#line 185 "slaqsp.f"
	    i__1 = *n;
#line 185 "slaqsp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 186 "slaqsp.f"
		cj = s[j];
#line 187 "slaqsp.f"
		i__2 = j;
#line 187 "slaqsp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 188 "slaqsp.f"
		    ap[jc + i__ - 1] = cj * s[i__] * ap[jc + i__ - 1];
#line 189 "slaqsp.f"
/* L10: */
#line 189 "slaqsp.f"
		}
#line 190 "slaqsp.f"
		jc += j;
#line 191 "slaqsp.f"
/* L20: */
#line 191 "slaqsp.f"
	    }
#line 192 "slaqsp.f"
	} else {

/*           Lower triangle of A is stored. */

#line 196 "slaqsp.f"
	    jc = 1;
#line 197 "slaqsp.f"
	    i__1 = *n;
#line 197 "slaqsp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 198 "slaqsp.f"
		cj = s[j];
#line 199 "slaqsp.f"
		i__2 = *n;
#line 199 "slaqsp.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 200 "slaqsp.f"
		    ap[jc + i__ - j] = cj * s[i__] * ap[jc + i__ - j];
#line 201 "slaqsp.f"
/* L30: */
#line 201 "slaqsp.f"
		}
#line 202 "slaqsp.f"
		jc = jc + *n - j + 1;
#line 203 "slaqsp.f"
/* L40: */
#line 203 "slaqsp.f"
	    }
#line 204 "slaqsp.f"
	}
#line 205 "slaqsp.f"
	*(unsigned char *)equed = 'Y';
#line 206 "slaqsp.f"
    }

#line 208 "slaqsp.f"
    return 0;

/*     End of SLAQSP */

} /* slaqsp_ */

