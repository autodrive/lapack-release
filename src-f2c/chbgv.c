#line 1 "chbgv.f"
/* chbgv.f -- translated by f2c (version 20100827).
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

#line 1 "chbgv.f"
/* > \brief \b CHBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbgv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbgv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbgv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, */
/*                         LDZ, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBGV computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian */
/* > and banded, and B is also positive definite. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* > \param[in] KA */
/* > \verbatim */
/* >          KA is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'. KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of superdiagonals of the matrix B if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'. KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first ka+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). */
/* > */
/* >          On exit, the contents of AB are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BB */
/* > \verbatim */
/* >          BB is COMPLEX array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**H*S, as returned by CPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* >          LDBB is INTEGER */
/* >          The leading dimension of the array BB.  LDBB >= KB+1. */
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
/* >          eigenvectors, with the i-th column of Z holding the */
/* >          eigenvector associated with W(i). The eigenvectors are */
/* >          normalized so that Z**H*B*Z = I. */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, and i is: */
/* >             <= N:  the algorithm failed to converge: */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then CPBSTF */
/* >                    returned INFO = i: B is not positive definite. */
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
/* Subroutine */ int chbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;

    /* Local variables */
    static integer inde;
    static char vect[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper, wantz;
    extern /* Subroutine */ int chbtrd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), chbgst_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    cpbstf_(char *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), ssterf_(integer *, doublereal *, doublereal *, integer *
	    );


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

#line 219 "chbgv.f"
    /* Parameter adjustments */
#line 219 "chbgv.f"
    ab_dim1 = *ldab;
#line 219 "chbgv.f"
    ab_offset = 1 + ab_dim1;
#line 219 "chbgv.f"
    ab -= ab_offset;
#line 219 "chbgv.f"
    bb_dim1 = *ldbb;
#line 219 "chbgv.f"
    bb_offset = 1 + bb_dim1;
#line 219 "chbgv.f"
    bb -= bb_offset;
#line 219 "chbgv.f"
    --w;
#line 219 "chbgv.f"
    z_dim1 = *ldz;
#line 219 "chbgv.f"
    z_offset = 1 + z_dim1;
#line 219 "chbgv.f"
    z__ -= z_offset;
#line 219 "chbgv.f"
    --work;
#line 219 "chbgv.f"
    --rwork;
#line 219 "chbgv.f"

#line 219 "chbgv.f"
    /* Function Body */
#line 219 "chbgv.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 220 "chbgv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 222 "chbgv.f"
    *info = 0;
#line 223 "chbgv.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 224 "chbgv.f"
	*info = -1;
#line 225 "chbgv.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 226 "chbgv.f"
	*info = -2;
#line 227 "chbgv.f"
    } else if (*n < 0) {
#line 228 "chbgv.f"
	*info = -3;
#line 229 "chbgv.f"
    } else if (*ka < 0) {
#line 230 "chbgv.f"
	*info = -4;
#line 231 "chbgv.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 232 "chbgv.f"
	*info = -5;
#line 233 "chbgv.f"
    } else if (*ldab < *ka + 1) {
#line 234 "chbgv.f"
	*info = -7;
#line 235 "chbgv.f"
    } else if (*ldbb < *kb + 1) {
#line 236 "chbgv.f"
	*info = -9;
#line 237 "chbgv.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 238 "chbgv.f"
	*info = -12;
#line 239 "chbgv.f"
    }
#line 240 "chbgv.f"
    if (*info != 0) {
#line 241 "chbgv.f"
	i__1 = -(*info);
#line 241 "chbgv.f"
	xerbla_("CHBGV ", &i__1, (ftnlen)6);
#line 242 "chbgv.f"
	return 0;
#line 243 "chbgv.f"
    }

/*     Quick return if possible */

#line 247 "chbgv.f"
    if (*n == 0) {
#line 247 "chbgv.f"
	return 0;
#line 247 "chbgv.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 252 "chbgv.f"
    cpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 253 "chbgv.f"
    if (*info != 0) {
#line 254 "chbgv.f"
	*info = *n + *info;
#line 255 "chbgv.f"
	return 0;
#line 256 "chbgv.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 260 "chbgv.f"
    inde = 1;
#line 261 "chbgv.f"
    indwrk = inde + *n;
#line 262 "chbgv.f"
    chbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &z__[z_offset], ldz, &work[1], &rwork[indwrk], &iinfo, (ftnlen)1,
	     (ftnlen)1);

/*     Reduce to tridiagonal form. */

#line 267 "chbgv.f"
    if (wantz) {
#line 268 "chbgv.f"
	*(unsigned char *)vect = 'U';
#line 269 "chbgv.f"
    } else {
#line 270 "chbgv.f"
	*(unsigned char *)vect = 'N';
#line 271 "chbgv.f"
    }
#line 272 "chbgv.f"
    chbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &w[1], &rwork[inde], &
	    z__[z_offset], ldz, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR. */

#line 277 "chbgv.f"
    if (! wantz) {
#line 278 "chbgv.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 279 "chbgv.f"
    } else {
#line 280 "chbgv.f"
	csteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[
		indwrk], info, (ftnlen)1);
#line 282 "chbgv.f"
    }
#line 283 "chbgv.f"
    return 0;

/*     End of CHBGV */

} /* chbgv_ */

