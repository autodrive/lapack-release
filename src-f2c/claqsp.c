#line 1 "claqsp.f"
/* claqsp.f -- translated by f2c (version 20100827).
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

#line 1 "claqsp.f"
/* > \brief \b CLAQSP scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by
 sppequ. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAQSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqsp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqsp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqsp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, UPLO */
/*       INTEGER            N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               S( * ) */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQSP equilibrates a symmetric matrix A using the scaling factors */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int claqsp_(char *uplo, integer *n, doublecomplex *ap, 
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

#line 163 "claqsp.f"
    /* Parameter adjustments */
#line 163 "claqsp.f"
    --s;
#line 163 "claqsp.f"
    --ap;
#line 163 "claqsp.f"

#line 163 "claqsp.f"
    /* Function Body */
#line 163 "claqsp.f"
    if (*n <= 0) {
#line 164 "claqsp.f"
	*(unsigned char *)equed = 'N';
#line 165 "claqsp.f"
	return 0;
#line 166 "claqsp.f"
    }

/*     Initialize LARGE and SMALL. */

#line 170 "claqsp.f"
    small = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 171 "claqsp.f"
    large = 1. / small;

#line 173 "claqsp.f"
    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

#line 177 "claqsp.f"
	*(unsigned char *)equed = 'N';
#line 178 "claqsp.f"
    } else {

/*        Replace A by diag(S) * A * diag(S). */

#line 182 "claqsp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Upper triangle of A is stored. */

#line 186 "claqsp.f"
	    jc = 1;
#line 187 "claqsp.f"
	    i__1 = *n;
#line 187 "claqsp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 188 "claqsp.f"
		cj = s[j];
#line 189 "claqsp.f"
		i__2 = j;
#line 189 "claqsp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 190 "claqsp.f"
		    i__3 = jc + i__ - 1;
#line 190 "claqsp.f"
		    d__1 = cj * s[i__];
#line 190 "claqsp.f"
		    i__4 = jc + i__ - 1;
#line 190 "claqsp.f"
		    z__1.r = d__1 * ap[i__4].r, z__1.i = d__1 * ap[i__4].i;
#line 190 "claqsp.f"
		    ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 191 "claqsp.f"
/* L10: */
#line 191 "claqsp.f"
		}
#line 192 "claqsp.f"
		jc += j;
#line 193 "claqsp.f"
/* L20: */
#line 193 "claqsp.f"
	    }
#line 194 "claqsp.f"
	} else {

/*           Lower triangle of A is stored. */

#line 198 "claqsp.f"
	    jc = 1;
#line 199 "claqsp.f"
	    i__1 = *n;
#line 199 "claqsp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 200 "claqsp.f"
		cj = s[j];
#line 201 "claqsp.f"
		i__2 = *n;
#line 201 "claqsp.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 202 "claqsp.f"
		    i__3 = jc + i__ - j;
#line 202 "claqsp.f"
		    d__1 = cj * s[i__];
#line 202 "claqsp.f"
		    i__4 = jc + i__ - j;
#line 202 "claqsp.f"
		    z__1.r = d__1 * ap[i__4].r, z__1.i = d__1 * ap[i__4].i;
#line 202 "claqsp.f"
		    ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 203 "claqsp.f"
/* L30: */
#line 203 "claqsp.f"
		}
#line 204 "claqsp.f"
		jc = jc + *n - j + 1;
#line 205 "claqsp.f"
/* L40: */
#line 205 "claqsp.f"
	    }
#line 206 "claqsp.f"
	}
#line 207 "claqsp.f"
	*(unsigned char *)equed = 'Y';
#line 208 "claqsp.f"
    }

#line 210 "claqsp.f"
    return 0;

/*     End of CLAQSP */

} /* claqsp_ */

