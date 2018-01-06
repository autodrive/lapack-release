#line 1 "spptrf.f"
/* spptrf.f -- translated by f2c (version 20100827).
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

#line 1 "spptrf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b16 = -1.;

/* > \brief \b SPPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPPTRF( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPPTRF computes the Cholesky factorization of a real symmetric */
/* > positive definite matrix A stored in packed format. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
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
/* >          See below for further details. */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**T*U or A = L*L**T, in the same */
/* >          storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i is not */
/* >                positive definite, and the factorization could not be */
/* >                completed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The packed storage scheme is illustrated by the following example */
/* >  when N = 4, UPLO = 'U': */
/* > */
/* >  Two-dimensional storage of the symmetric matrix A: */
/* > */
/* >     a11 a12 a13 a14 */
/* >         a22 a23 a24 */
/* >             a33 a34     (aij = aji) */
/* >                 a44 */
/* > */
/* >  Packed storage of the upper triangle of A: */
/* > */
/* >  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ] */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int spptrf_(char *uplo, integer *n, doublereal *ap, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int sspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int stpsv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 161 "spptrf.f"
    /* Parameter adjustments */
#line 161 "spptrf.f"
    --ap;
#line 161 "spptrf.f"

#line 161 "spptrf.f"
    /* Function Body */
#line 161 "spptrf.f"
    *info = 0;
#line 162 "spptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "spptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "spptrf.f"
	*info = -1;
#line 165 "spptrf.f"
    } else if (*n < 0) {
#line 166 "spptrf.f"
	*info = -2;
#line 167 "spptrf.f"
    }
#line 168 "spptrf.f"
    if (*info != 0) {
#line 169 "spptrf.f"
	i__1 = -(*info);
#line 169 "spptrf.f"
	xerbla_("SPPTRF", &i__1, (ftnlen)6);
#line 170 "spptrf.f"
	return 0;
#line 171 "spptrf.f"
    }

/*     Quick return if possible */

#line 175 "spptrf.f"
    if (*n == 0) {
#line 175 "spptrf.f"
	return 0;
#line 175 "spptrf.f"
    }

#line 178 "spptrf.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**T*U. */

#line 182 "spptrf.f"
	jj = 0;
#line 183 "spptrf.f"
	i__1 = *n;
#line 183 "spptrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "spptrf.f"
	    jc = jj + 1;
#line 185 "spptrf.f"
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

#line 189 "spptrf.f"
	    if (j > 1) {
#line 189 "spptrf.f"
		i__2 = j - 1;
#line 189 "spptrf.f"
		stpsv_("Upper", "Transpose", "Non-unit", &i__2, &ap[1], &ap[
			jc], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 189 "spptrf.f"
	    }

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 195 "spptrf.f"
	    i__2 = j - 1;
#line 195 "spptrf.f"
	    ajj = ap[jj] - sdot_(&i__2, &ap[jc], &c__1, &ap[jc], &c__1);
#line 196 "spptrf.f"
	    if (ajj <= 0.) {
#line 197 "spptrf.f"
		ap[jj] = ajj;
#line 198 "spptrf.f"
		goto L30;
#line 199 "spptrf.f"
	    }
#line 200 "spptrf.f"
	    ap[jj] = sqrt(ajj);
#line 201 "spptrf.f"
/* L10: */
#line 201 "spptrf.f"
	}
#line 202 "spptrf.f"
    } else {

/*        Compute the Cholesky factorization A = L*L**T. */

#line 206 "spptrf.f"
	jj = 1;
#line 207 "spptrf.f"
	i__1 = *n;
#line 207 "spptrf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 211 "spptrf.f"
	    ajj = ap[jj];
#line 212 "spptrf.f"
	    if (ajj <= 0.) {
#line 213 "spptrf.f"
		ap[jj] = ajj;
#line 214 "spptrf.f"
		goto L30;
#line 215 "spptrf.f"
	    }
#line 216 "spptrf.f"
	    ajj = sqrt(ajj);
#line 217 "spptrf.f"
	    ap[jj] = ajj;

/*           Compute elements J+1:N of column J and update the trailing */
/*           submatrix. */

#line 222 "spptrf.f"
	    if (j < *n) {
#line 223 "spptrf.f"
		i__2 = *n - j;
#line 223 "spptrf.f"
		d__1 = 1. / ajj;
#line 223 "spptrf.f"
		sscal_(&i__2, &d__1, &ap[jj + 1], &c__1);
#line 224 "spptrf.f"
		i__2 = *n - j;
#line 224 "spptrf.f"
		sspr_("Lower", &i__2, &c_b16, &ap[jj + 1], &c__1, &ap[jj + *n 
			- j + 1], (ftnlen)5);
#line 226 "spptrf.f"
		jj = jj + *n - j + 1;
#line 227 "spptrf.f"
	    }
#line 228 "spptrf.f"
/* L20: */
#line 228 "spptrf.f"
	}
#line 229 "spptrf.f"
    }
#line 230 "spptrf.f"
    goto L40;

#line 232 "spptrf.f"
L30:
#line 233 "spptrf.f"
    *info = j;

#line 235 "spptrf.f"
L40:
#line 236 "spptrf.f"
    return 0;

/*     End of SPPTRF */

} /* spptrf_ */

