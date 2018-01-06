#line 1 "dspgv.f"
/* dspgv.f -- translated by f2c (version 20100827).
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

#line 1 "dspgv.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSPGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspgv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspgv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspgv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPGV computes all the eigenvalues and, optionally, the eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be symmetric, stored in packed format, */
/* > and B is also positive definite. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          Specifies the problem type to be solved: */
/* >          = 1:  A*x = (lambda)*B*x */
/* >          = 2:  A*B*x = (lambda)*x */
/* >          = 3:  B*A*x = (lambda)*x */
/* > \endverbatim */
/* > */
/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only; */
/* >          = 'V':  Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangles of A and B are stored; */
/* >          = 'L':  Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension */
/* >                            (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the contents of AP are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BP */
/* > \verbatim */
/* >          BP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          B, packed columnwise in a linear array.  The j-th column of B */
/* >          is stored in the array BP as follows: */
/* >          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the triangular factor U or L from the Cholesky */
/* >          factorization B = U**T*U or B = L*L**T, in the same storage */
/* >          format as B. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors.  The eigenvectors are normalized as follows: */
/* >          if ITYPE = 1 or 2, Z**T*B*Z = I; */
/* >          if ITYPE = 3, Z**T*inv(B)*Z = I. */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  DPPTRF or DSPEV returned an error code: */
/* >             <= N:  if INFO = i, DSPEV failed to converge; */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero. */
/* >             > N:   if INFO = n + i, for 1 <= i <= n, then the leading */
/* >                    minor of order i of B is not positive definite. */
/* >                    The factorization of B could not be completed and */
/* >                    no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dspgv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *info, ftnlen jobz_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;

    /* Local variables */
    static integer j, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dspev_(char *, char *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static char trans[1];
    static logical upper;
    extern /* Subroutine */ int dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dtpsv_(char *, char *, char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical wantz;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpptrf_(
	    char *, integer *, doublereal *, integer *, ftnlen), dspgst_(
	    integer *, char *, integer *, doublereal *, doublereal *, integer 
	    *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 196 "dspgv.f"
    /* Parameter adjustments */
#line 196 "dspgv.f"
    --ap;
#line 196 "dspgv.f"
    --bp;
#line 196 "dspgv.f"
    --w;
#line 196 "dspgv.f"
    z_dim1 = *ldz;
#line 196 "dspgv.f"
    z_offset = 1 + z_dim1;
#line 196 "dspgv.f"
    z__ -= z_offset;
#line 196 "dspgv.f"
    --work;
#line 196 "dspgv.f"

#line 196 "dspgv.f"
    /* Function Body */
#line 196 "dspgv.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 197 "dspgv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 199 "dspgv.f"
    *info = 0;
#line 200 "dspgv.f"
    if (*itype < 1 || *itype > 3) {
#line 201 "dspgv.f"
	*info = -1;
#line 202 "dspgv.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 203 "dspgv.f"
	*info = -2;
#line 204 "dspgv.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 205 "dspgv.f"
	*info = -3;
#line 206 "dspgv.f"
    } else if (*n < 0) {
#line 207 "dspgv.f"
	*info = -4;
#line 208 "dspgv.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 209 "dspgv.f"
	*info = -9;
#line 210 "dspgv.f"
    }
#line 211 "dspgv.f"
    if (*info != 0) {
#line 212 "dspgv.f"
	i__1 = -(*info);
#line 212 "dspgv.f"
	xerbla_("DSPGV ", &i__1, (ftnlen)6);
#line 213 "dspgv.f"
	return 0;
#line 214 "dspgv.f"
    }

/*     Quick return if possible */

#line 218 "dspgv.f"
    if (*n == 0) {
#line 218 "dspgv.f"
	return 0;
#line 218 "dspgv.f"
    }

/*     Form a Cholesky factorization of B. */

#line 223 "dspgv.f"
    dpptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 224 "dspgv.f"
    if (*info != 0) {
#line 225 "dspgv.f"
	*info = *n + *info;
#line 226 "dspgv.f"
	return 0;
#line 227 "dspgv.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 231 "dspgv.f"
    dspgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 232 "dspgv.f"
    dspev_(jobz, uplo, n, &ap[1], &w[1], &z__[z_offset], ldz, &work[1], info, 
	    (ftnlen)1, (ftnlen)1);

#line 234 "dspgv.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 238 "dspgv.f"
	neig = *n;
#line 239 "dspgv.f"
	if (*info > 0) {
#line 239 "dspgv.f"
	    neig = *info - 1;
#line 239 "dspgv.f"
	}
#line 241 "dspgv.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 246 "dspgv.f"
	    if (upper) {
#line 247 "dspgv.f"
		*(unsigned char *)trans = 'N';
#line 248 "dspgv.f"
	    } else {
#line 249 "dspgv.f"
		*(unsigned char *)trans = 'T';
#line 250 "dspgv.f"
	    }

#line 252 "dspgv.f"
	    i__1 = neig;
#line 252 "dspgv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 253 "dspgv.f"
		dtpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 255 "dspgv.f"
/* L10: */
#line 255 "dspgv.f"
	    }

#line 257 "dspgv.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 262 "dspgv.f"
	    if (upper) {
#line 263 "dspgv.f"
		*(unsigned char *)trans = 'T';
#line 264 "dspgv.f"
	    } else {
#line 265 "dspgv.f"
		*(unsigned char *)trans = 'N';
#line 266 "dspgv.f"
	    }

#line 268 "dspgv.f"
	    i__1 = neig;
#line 268 "dspgv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 269 "dspgv.f"
		dtpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 271 "dspgv.f"
/* L20: */
#line 271 "dspgv.f"
	    }
#line 272 "dspgv.f"
	}
#line 273 "dspgv.f"
    }
#line 274 "dspgv.f"
    return 0;

/*     End of DSPGV */

} /* dspgv_ */

