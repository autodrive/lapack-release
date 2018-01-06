#line 1 "cpptrf.f"
/* cpptrf.f -- translated by f2c (version 20100827).
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

#line 1 "cpptrf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b16 = -1.;

/* > \brief \b CPPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPPTRF( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPTRF computes the Cholesky factorization of a complex Hermitian */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int cpptrf_(char *uplo, integer *n, doublecomplex *ap, 
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
    extern /* Subroutine */ int chpr_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), csscal_(integer *, doublereal *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 161 "cpptrf.f"
    /* Parameter adjustments */
#line 161 "cpptrf.f"
    --ap;
#line 161 "cpptrf.f"

#line 161 "cpptrf.f"
    /* Function Body */
#line 161 "cpptrf.f"
    *info = 0;
#line 162 "cpptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "cpptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "cpptrf.f"
	*info = -1;
#line 165 "cpptrf.f"
    } else if (*n < 0) {
#line 166 "cpptrf.f"
	*info = -2;
#line 167 "cpptrf.f"
    }
#line 168 "cpptrf.f"
    if (*info != 0) {
#line 169 "cpptrf.f"
	i__1 = -(*info);
#line 169 "cpptrf.f"
	xerbla_("CPPTRF", &i__1, (ftnlen)6);
#line 170 "cpptrf.f"
	return 0;
#line 171 "cpptrf.f"
    }

/*     Quick return if possible */

#line 175 "cpptrf.f"
    if (*n == 0) {
#line 175 "cpptrf.f"
	return 0;
#line 175 "cpptrf.f"
    }

#line 178 "cpptrf.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**H * U. */

#line 182 "cpptrf.f"
	jj = 0;
#line 183 "cpptrf.f"
	i__1 = *n;
#line 183 "cpptrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "cpptrf.f"
	    jc = jj + 1;
#line 185 "cpptrf.f"
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

#line 189 "cpptrf.f"
	    if (j > 1) {
#line 189 "cpptrf.f"
		i__2 = j - 1;
#line 189 "cpptrf.f"
		ctpsv_("Upper", "Conjugate transpose", "Non-unit", &i__2, &ap[
			1], &ap[jc], &c__1, (ftnlen)5, (ftnlen)19, (ftnlen)8);
#line 189 "cpptrf.f"
	    }

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 195 "cpptrf.f"
	    i__2 = jj;
#line 195 "cpptrf.f"
	    d__1 = ap[i__2].r;
#line 195 "cpptrf.f"
	    i__3 = j - 1;
#line 195 "cpptrf.f"
	    cdotc_(&z__2, &i__3, &ap[jc], &c__1, &ap[jc], &c__1);
#line 195 "cpptrf.f"
	    z__1.r = d__1 - z__2.r, z__1.i = -z__2.i;
#line 195 "cpptrf.f"
	    ajj = z__1.r;
#line 197 "cpptrf.f"
	    if (ajj <= 0.) {
#line 198 "cpptrf.f"
		i__2 = jj;
#line 198 "cpptrf.f"
		ap[i__2].r = ajj, ap[i__2].i = 0.;
#line 199 "cpptrf.f"
		goto L30;
#line 200 "cpptrf.f"
	    }
#line 201 "cpptrf.f"
	    i__2 = jj;
#line 201 "cpptrf.f"
	    d__1 = sqrt(ajj);
#line 201 "cpptrf.f"
	    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 202 "cpptrf.f"
/* L10: */
#line 202 "cpptrf.f"
	}
#line 203 "cpptrf.f"
    } else {

/*        Compute the Cholesky factorization A = L * L**H. */

#line 207 "cpptrf.f"
	jj = 1;
#line 208 "cpptrf.f"
	i__1 = *n;
#line 208 "cpptrf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 212 "cpptrf.f"
	    i__2 = jj;
#line 212 "cpptrf.f"
	    ajj = ap[i__2].r;
#line 213 "cpptrf.f"
	    if (ajj <= 0.) {
#line 214 "cpptrf.f"
		i__2 = jj;
#line 214 "cpptrf.f"
		ap[i__2].r = ajj, ap[i__2].i = 0.;
#line 215 "cpptrf.f"
		goto L30;
#line 216 "cpptrf.f"
	    }
#line 217 "cpptrf.f"
	    ajj = sqrt(ajj);
#line 218 "cpptrf.f"
	    i__2 = jj;
#line 218 "cpptrf.f"
	    ap[i__2].r = ajj, ap[i__2].i = 0.;

/*           Compute elements J+1:N of column J and update the trailing */
/*           submatrix. */

#line 223 "cpptrf.f"
	    if (j < *n) {
#line 224 "cpptrf.f"
		i__2 = *n - j;
#line 224 "cpptrf.f"
		d__1 = 1. / ajj;
#line 224 "cpptrf.f"
		csscal_(&i__2, &d__1, &ap[jj + 1], &c__1);
#line 225 "cpptrf.f"
		i__2 = *n - j;
#line 225 "cpptrf.f"
		chpr_("Lower", &i__2, &c_b16, &ap[jj + 1], &c__1, &ap[jj + *n 
			- j + 1], (ftnlen)5);
#line 227 "cpptrf.f"
		jj = jj + *n - j + 1;
#line 228 "cpptrf.f"
	    }
#line 229 "cpptrf.f"
/* L20: */
#line 229 "cpptrf.f"
	}
#line 230 "cpptrf.f"
    }
#line 231 "cpptrf.f"
    goto L40;

#line 233 "cpptrf.f"
L30:
#line 234 "cpptrf.f"
    *info = j;

#line 236 "cpptrf.f"
L40:
#line 237 "cpptrf.f"
    return 0;

/*     End of CPPTRF */

} /* cpptrf_ */

