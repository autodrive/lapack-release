#line 1 "chpgv.f"
/* chpgv.f -- translated by f2c (version 20100827).
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

#line 1 "chpgv.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHPGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpgv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpgv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpgv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, */
/*                         RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AP( * ), BP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPGV computes all the eigenvalues and, optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be Hermitian, stored in packed format, */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
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
/* >          BP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          B, packed columnwise in a linear array.  The j-th column of B */
/* >          is stored in the array BP as follows: */
/* >          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the triangular factor U or L from the Cholesky */
/* >          factorization B = U**H*U or B = L*L**H, in the same storage */
/* >          format as B. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors.  The eigenvectors are normalized as follows: */
/* >          if ITYPE = 1 or 2, Z**H*B*Z = I; */
/* >          if ITYPE = 3, Z**H*inv(B)*Z = I. */
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
/* >          WORK is COMPLEX array, dimension (max(1, 2*N-1)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (max(1, 3*N-2)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  CPPTRF or CHPEV returned an error code: */
/* >             <= N:  if INFO = i, CHPEV failed to converge; */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not convergeto zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= n, then the leading */
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

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int chpgv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;

    /* Local variables */
    static integer j, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int chpev_(char *, char *, integer *, 
	    doublecomplex *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, doublereal *, integer *, ftnlen, ftnlen);
    static char trans[1];
    extern /* Subroutine */ int ctpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical wantz;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), chpgst_(
	    integer *, char *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), cpptrf_(char *, integer *, doublecomplex *, 
	    integer *, ftnlen);


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

#line 200 "chpgv.f"
    /* Parameter adjustments */
#line 200 "chpgv.f"
    --ap;
#line 200 "chpgv.f"
    --bp;
#line 200 "chpgv.f"
    --w;
#line 200 "chpgv.f"
    z_dim1 = *ldz;
#line 200 "chpgv.f"
    z_offset = 1 + z_dim1;
#line 200 "chpgv.f"
    z__ -= z_offset;
#line 200 "chpgv.f"
    --work;
#line 200 "chpgv.f"
    --rwork;
#line 200 "chpgv.f"

#line 200 "chpgv.f"
    /* Function Body */
#line 200 "chpgv.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 201 "chpgv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 203 "chpgv.f"
    *info = 0;
#line 204 "chpgv.f"
    if (*itype < 1 || *itype > 3) {
#line 205 "chpgv.f"
	*info = -1;
#line 206 "chpgv.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 207 "chpgv.f"
	*info = -2;
#line 208 "chpgv.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 209 "chpgv.f"
	*info = -3;
#line 210 "chpgv.f"
    } else if (*n < 0) {
#line 211 "chpgv.f"
	*info = -4;
#line 212 "chpgv.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 213 "chpgv.f"
	*info = -9;
#line 214 "chpgv.f"
    }
#line 215 "chpgv.f"
    if (*info != 0) {
#line 216 "chpgv.f"
	i__1 = -(*info);
#line 216 "chpgv.f"
	xerbla_("CHPGV ", &i__1, (ftnlen)6);
#line 217 "chpgv.f"
	return 0;
#line 218 "chpgv.f"
    }

/*     Quick return if possible */

#line 222 "chpgv.f"
    if (*n == 0) {
#line 222 "chpgv.f"
	return 0;
#line 222 "chpgv.f"
    }

/*     Form a Cholesky factorization of B. */

#line 227 "chpgv.f"
    cpptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 228 "chpgv.f"
    if (*info != 0) {
#line 229 "chpgv.f"
	*info = *n + *info;
#line 230 "chpgv.f"
	return 0;
#line 231 "chpgv.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 235 "chpgv.f"
    chpgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 236 "chpgv.f"
    chpev_(jobz, uplo, n, &ap[1], &w[1], &z__[z_offset], ldz, &work[1], &
	    rwork[1], info, (ftnlen)1, (ftnlen)1);

#line 238 "chpgv.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 242 "chpgv.f"
	neig = *n;
#line 243 "chpgv.f"
	if (*info > 0) {
#line 243 "chpgv.f"
	    neig = *info - 1;
#line 243 "chpgv.f"
	}
#line 245 "chpgv.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y */

#line 250 "chpgv.f"
	    if (upper) {
#line 251 "chpgv.f"
		*(unsigned char *)trans = 'N';
#line 252 "chpgv.f"
	    } else {
#line 253 "chpgv.f"
		*(unsigned char *)trans = 'C';
#line 254 "chpgv.f"
	    }

#line 256 "chpgv.f"
	    i__1 = neig;
#line 256 "chpgv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 257 "chpgv.f"
		ctpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 259 "chpgv.f"
/* L10: */
#line 259 "chpgv.f"
	    }

#line 261 "chpgv.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H*y */

#line 266 "chpgv.f"
	    if (upper) {
#line 267 "chpgv.f"
		*(unsigned char *)trans = 'C';
#line 268 "chpgv.f"
	    } else {
#line 269 "chpgv.f"
		*(unsigned char *)trans = 'N';
#line 270 "chpgv.f"
	    }

#line 272 "chpgv.f"
	    i__1 = neig;
#line 272 "chpgv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 273 "chpgv.f"
		ctpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 275 "chpgv.f"
/* L20: */
#line 275 "chpgv.f"
	    }
#line 276 "chpgv.f"
	}
#line 277 "chpgv.f"
    }
#line 278 "chpgv.f"
    return 0;

/*     End of CHPGV */

} /* chpgv_ */

