#line 1 "dpptrf.f"
/* dpptrf.f -- translated by f2c (version 20100827).
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

#line 1 "dpptrf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b16 = -1.;

/* > \brief \b DPPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPPTRF( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPPTRF computes the Cholesky factorization of a real symmetric */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dpptrf_(char *uplo, integer *n, doublereal *ap, integer *
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
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dtpsv_(char *, char *, char *, integer *, 
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

#line 161 "dpptrf.f"
    /* Parameter adjustments */
#line 161 "dpptrf.f"
    --ap;
#line 161 "dpptrf.f"

#line 161 "dpptrf.f"
    /* Function Body */
#line 161 "dpptrf.f"
    *info = 0;
#line 162 "dpptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "dpptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "dpptrf.f"
	*info = -1;
#line 165 "dpptrf.f"
    } else if (*n < 0) {
#line 166 "dpptrf.f"
	*info = -2;
#line 167 "dpptrf.f"
    }
#line 168 "dpptrf.f"
    if (*info != 0) {
#line 169 "dpptrf.f"
	i__1 = -(*info);
#line 169 "dpptrf.f"
	xerbla_("DPPTRF", &i__1, (ftnlen)6);
#line 170 "dpptrf.f"
	return 0;
#line 171 "dpptrf.f"
    }

/*     Quick return if possible */

#line 175 "dpptrf.f"
    if (*n == 0) {
#line 175 "dpptrf.f"
	return 0;
#line 175 "dpptrf.f"
    }

#line 178 "dpptrf.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**T*U. */

#line 182 "dpptrf.f"
	jj = 0;
#line 183 "dpptrf.f"
	i__1 = *n;
#line 183 "dpptrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "dpptrf.f"
	    jc = jj + 1;
#line 185 "dpptrf.f"
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

#line 189 "dpptrf.f"
	    if (j > 1) {
#line 189 "dpptrf.f"
		i__2 = j - 1;
#line 189 "dpptrf.f"
		dtpsv_("Upper", "Transpose", "Non-unit", &i__2, &ap[1], &ap[
			jc], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 189 "dpptrf.f"
	    }

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 195 "dpptrf.f"
	    i__2 = j - 1;
#line 195 "dpptrf.f"
	    ajj = ap[jj] - ddot_(&i__2, &ap[jc], &c__1, &ap[jc], &c__1);
#line 196 "dpptrf.f"
	    if (ajj <= 0.) {
#line 197 "dpptrf.f"
		ap[jj] = ajj;
#line 198 "dpptrf.f"
		goto L30;
#line 199 "dpptrf.f"
	    }
#line 200 "dpptrf.f"
	    ap[jj] = sqrt(ajj);
#line 201 "dpptrf.f"
/* L10: */
#line 201 "dpptrf.f"
	}
#line 202 "dpptrf.f"
    } else {

/*        Compute the Cholesky factorization A = L*L**T. */

#line 206 "dpptrf.f"
	jj = 1;
#line 207 "dpptrf.f"
	i__1 = *n;
#line 207 "dpptrf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 211 "dpptrf.f"
	    ajj = ap[jj];
#line 212 "dpptrf.f"
	    if (ajj <= 0.) {
#line 213 "dpptrf.f"
		ap[jj] = ajj;
#line 214 "dpptrf.f"
		goto L30;
#line 215 "dpptrf.f"
	    }
#line 216 "dpptrf.f"
	    ajj = sqrt(ajj);
#line 217 "dpptrf.f"
	    ap[jj] = ajj;

/*           Compute elements J+1:N of column J and update the trailing */
/*           submatrix. */

#line 222 "dpptrf.f"
	    if (j < *n) {
#line 223 "dpptrf.f"
		i__2 = *n - j;
#line 223 "dpptrf.f"
		d__1 = 1. / ajj;
#line 223 "dpptrf.f"
		dscal_(&i__2, &d__1, &ap[jj + 1], &c__1);
#line 224 "dpptrf.f"
		i__2 = *n - j;
#line 224 "dpptrf.f"
		dspr_("Lower", &i__2, &c_b16, &ap[jj + 1], &c__1, &ap[jj + *n 
			- j + 1], (ftnlen)5);
#line 226 "dpptrf.f"
		jj = jj + *n - j + 1;
#line 227 "dpptrf.f"
	    }
#line 228 "dpptrf.f"
/* L20: */
#line 228 "dpptrf.f"
	}
#line 229 "dpptrf.f"
    }
#line 230 "dpptrf.f"
    goto L40;

#line 232 "dpptrf.f"
L30:
#line 233 "dpptrf.f"
    *info = j;

#line 235 "dpptrf.f"
L40:
#line 236 "dpptrf.f"
    return 0;

/*     End of DPPTRF */

} /* dpptrf_ */

