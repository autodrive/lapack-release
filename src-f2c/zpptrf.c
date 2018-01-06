#line 1 "zpptrf.f"
/* zpptrf.f -- translated by f2c (version 20100827).
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

#line 1 "zpptrf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b16 = -1.;

/* > \brief \b ZPPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPPTRF( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPTRF computes the Cholesky factorization of a complex Hermitian */
/* > positive definite matrix A stored in packed format. */
/* > */
/* > The factorization has the form */
/* >    A = U**H * U,  if UPLO = 'U', or */
/* >    A = L  * L**H,  if UPLO = 'L', */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* >          See below for further details. */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**H*U or A = L*L**H, in the same */
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

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The packed storage scheme is illustrated by the following example */
/* >  when N = 4, UPLO = 'U': */
/* > */
/* >  Two-dimensional storage of the Hermitian matrix A: */
/* > */
/* >     a11 a12 a13 a14 */
/* >         a22 a23 a24 */
/* >             a33 a34     (aij = conjg(aji)) */
/* >                 a44 */
/* > */
/* >  Packed storage of the upper triangle of A: */
/* > */
/* >  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ] */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zpptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    extern /* Subroutine */ int zhpr_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ztpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), zdscal_(integer *, 
	    doublereal *, doublecomplex *, integer *);


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

#line 161 "zpptrf.f"
    /* Parameter adjustments */
#line 161 "zpptrf.f"
    --ap;
#line 161 "zpptrf.f"

#line 161 "zpptrf.f"
    /* Function Body */
#line 161 "zpptrf.f"
    *info = 0;
#line 162 "zpptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "zpptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "zpptrf.f"
	*info = -1;
#line 165 "zpptrf.f"
    } else if (*n < 0) {
#line 166 "zpptrf.f"
	*info = -2;
#line 167 "zpptrf.f"
    }
#line 168 "zpptrf.f"
    if (*info != 0) {
#line 169 "zpptrf.f"
	i__1 = -(*info);
#line 169 "zpptrf.f"
	xerbla_("ZPPTRF", &i__1, (ftnlen)6);
#line 170 "zpptrf.f"
	return 0;
#line 171 "zpptrf.f"
    }

/*     Quick return if possible */

#line 175 "zpptrf.f"
    if (*n == 0) {
#line 175 "zpptrf.f"
	return 0;
#line 175 "zpptrf.f"
    }

#line 178 "zpptrf.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**H * U. */

#line 182 "zpptrf.f"
	jj = 0;
#line 183 "zpptrf.f"
	i__1 = *n;
#line 183 "zpptrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "zpptrf.f"
	    jc = jj + 1;
#line 185 "zpptrf.f"
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

#line 189 "zpptrf.f"
	    if (j > 1) {
#line 189 "zpptrf.f"
		i__2 = j - 1;
#line 189 "zpptrf.f"
		ztpsv_("Upper", "Conjugate transpose", "Non-unit", &i__2, &ap[
			1], &ap[jc], &c__1, (ftnlen)5, (ftnlen)19, (ftnlen)8);
#line 189 "zpptrf.f"
	    }

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 195 "zpptrf.f"
	    i__2 = jj;
#line 195 "zpptrf.f"
	    d__1 = ap[i__2].r;
#line 195 "zpptrf.f"
	    i__3 = j - 1;
#line 195 "zpptrf.f"
	    zdotc_(&z__2, &i__3, &ap[jc], &c__1, &ap[jc], &c__1);
#line 195 "zpptrf.f"
	    z__1.r = d__1 - z__2.r, z__1.i = -z__2.i;
#line 195 "zpptrf.f"
	    ajj = z__1.r;
#line 197 "zpptrf.f"
	    if (ajj <= 0.) {
#line 198 "zpptrf.f"
		i__2 = jj;
#line 198 "zpptrf.f"
		ap[i__2].r = ajj, ap[i__2].i = 0.;
#line 199 "zpptrf.f"
		goto L30;
#line 200 "zpptrf.f"
	    }
#line 201 "zpptrf.f"
	    i__2 = jj;
#line 201 "zpptrf.f"
	    d__1 = sqrt(ajj);
#line 201 "zpptrf.f"
	    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 202 "zpptrf.f"
/* L10: */
#line 202 "zpptrf.f"
	}
#line 203 "zpptrf.f"
    } else {

/*        Compute the Cholesky factorization A = L * L**H. */

#line 207 "zpptrf.f"
	jj = 1;
#line 208 "zpptrf.f"
	i__1 = *n;
#line 208 "zpptrf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 212 "zpptrf.f"
	    i__2 = jj;
#line 212 "zpptrf.f"
	    ajj = ap[i__2].r;
#line 213 "zpptrf.f"
	    if (ajj <= 0.) {
#line 214 "zpptrf.f"
		i__2 = jj;
#line 214 "zpptrf.f"
		ap[i__2].r = ajj, ap[i__2].i = 0.;
#line 215 "zpptrf.f"
		goto L30;
#line 216 "zpptrf.f"
	    }
#line 217 "zpptrf.f"
	    ajj = sqrt(ajj);
#line 218 "zpptrf.f"
	    i__2 = jj;
#line 218 "zpptrf.f"
	    ap[i__2].r = ajj, ap[i__2].i = 0.;

/*           Compute elements J+1:N of column J and update the trailing */
/*           submatrix. */

#line 223 "zpptrf.f"
	    if (j < *n) {
#line 224 "zpptrf.f"
		i__2 = *n - j;
#line 224 "zpptrf.f"
		d__1 = 1. / ajj;
#line 224 "zpptrf.f"
		zdscal_(&i__2, &d__1, &ap[jj + 1], &c__1);
#line 225 "zpptrf.f"
		i__2 = *n - j;
#line 225 "zpptrf.f"
		zhpr_("Lower", &i__2, &c_b16, &ap[jj + 1], &c__1, &ap[jj + *n 
			- j + 1], (ftnlen)5);
#line 227 "zpptrf.f"
		jj = jj + *n - j + 1;
#line 228 "zpptrf.f"
	    }
#line 229 "zpptrf.f"
/* L20: */
#line 229 "zpptrf.f"
	}
#line 230 "zpptrf.f"
    }
#line 231 "zpptrf.f"
    goto L40;

#line 233 "zpptrf.f"
L30:
#line 234 "zpptrf.f"
    *info = j;

#line 236 "zpptrf.f"
L40:
#line 237 "zpptrf.f"
    return 0;

/*     End of ZPPTRF */

} /* zpptrf_ */

