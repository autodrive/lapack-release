#line 1 "dsbgv.f"
/* dsbgv.f -- translated by f2c (version 20100827).
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

#line 1 "dsbgv.f"
/* > \brief \b DSBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, */
/*                         LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), W( * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBGV computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
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
/* >          BB is DOUBLE PRECISION array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**T*S, as returned by DPBSTF. */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors, with the i-th column of Z holding the */
/* >          eigenvector associated with W(i). The eigenvectors are */
/* >          normalized so that Z**T*B*Z = I. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
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
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then DPBSTF */
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

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dsbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;

    /* Local variables */
    static integer inde;
    static char vect[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper, wantz;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpbstf_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dsbtrd_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen), dsbgst_(char *, char *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen), dsterf_(integer *, doublereal *, 
	    doublereal *, integer *);
    static integer indwrk;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);


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

#line 212 "dsbgv.f"
    /* Parameter adjustments */
#line 212 "dsbgv.f"
    ab_dim1 = *ldab;
#line 212 "dsbgv.f"
    ab_offset = 1 + ab_dim1;
#line 212 "dsbgv.f"
    ab -= ab_offset;
#line 212 "dsbgv.f"
    bb_dim1 = *ldbb;
#line 212 "dsbgv.f"
    bb_offset = 1 + bb_dim1;
#line 212 "dsbgv.f"
    bb -= bb_offset;
#line 212 "dsbgv.f"
    --w;
#line 212 "dsbgv.f"
    z_dim1 = *ldz;
#line 212 "dsbgv.f"
    z_offset = 1 + z_dim1;
#line 212 "dsbgv.f"
    z__ -= z_offset;
#line 212 "dsbgv.f"
    --work;
#line 212 "dsbgv.f"

#line 212 "dsbgv.f"
    /* Function Body */
#line 212 "dsbgv.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 213 "dsbgv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 215 "dsbgv.f"
    *info = 0;
#line 216 "dsbgv.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 217 "dsbgv.f"
	*info = -1;
#line 218 "dsbgv.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 219 "dsbgv.f"
	*info = -2;
#line 220 "dsbgv.f"
    } else if (*n < 0) {
#line 221 "dsbgv.f"
	*info = -3;
#line 222 "dsbgv.f"
    } else if (*ka < 0) {
#line 223 "dsbgv.f"
	*info = -4;
#line 224 "dsbgv.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 225 "dsbgv.f"
	*info = -5;
#line 226 "dsbgv.f"
    } else if (*ldab < *ka + 1) {
#line 227 "dsbgv.f"
	*info = -7;
#line 228 "dsbgv.f"
    } else if (*ldbb < *kb + 1) {
#line 229 "dsbgv.f"
	*info = -9;
#line 230 "dsbgv.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 231 "dsbgv.f"
	*info = -12;
#line 232 "dsbgv.f"
    }
#line 233 "dsbgv.f"
    if (*info != 0) {
#line 234 "dsbgv.f"
	i__1 = -(*info);
#line 234 "dsbgv.f"
	xerbla_("DSBGV ", &i__1, (ftnlen)6);
#line 235 "dsbgv.f"
	return 0;
#line 236 "dsbgv.f"
    }

/*     Quick return if possible */

#line 240 "dsbgv.f"
    if (*n == 0) {
#line 240 "dsbgv.f"
	return 0;
#line 240 "dsbgv.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 245 "dsbgv.f"
    dpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 246 "dsbgv.f"
    if (*info != 0) {
#line 247 "dsbgv.f"
	*info = *n + *info;
#line 248 "dsbgv.f"
	return 0;
#line 249 "dsbgv.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 253 "dsbgv.f"
    inde = 1;
#line 254 "dsbgv.f"
    indwrk = inde + *n;
#line 255 "dsbgv.f"
    dsbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &z__[z_offset], ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1)
	    ;

/*     Reduce to tridiagonal form. */

#line 260 "dsbgv.f"
    if (wantz) {
#line 261 "dsbgv.f"
	*(unsigned char *)vect = 'U';
#line 262 "dsbgv.f"
    } else {
#line 263 "dsbgv.f"
	*(unsigned char *)vect = 'N';
#line 264 "dsbgv.f"
    }
#line 265 "dsbgv.f"
    dsbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[
	    z_offset], ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR. */

#line 270 "dsbgv.f"
    if (! wantz) {
#line 271 "dsbgv.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 272 "dsbgv.f"
    } else {
#line 273 "dsbgv.f"
	dsteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[
		indwrk], info, (ftnlen)1);
#line 275 "dsbgv.f"
    }
#line 276 "dsbgv.f"
    return 0;

/*     End of DSBGV */

} /* dsbgv_ */

