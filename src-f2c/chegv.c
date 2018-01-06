#line 1 "chegv.f"
/* chegv.f -- translated by f2c (version 20100827).
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

#line 1 "chegv.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CHEGV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, */
/*                         LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEGV computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be Hermitian and B is also */
/* > positive definite. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* > */
/* >          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* >          matrix Z of eigenvectors.  The eigenvectors are normalized */
/* >          as follows: */
/* >          if ITYPE = 1 or 2, Z**H*B*Z = I; */
/* >          if ITYPE = 3, Z**H*inv(B)*Z = I. */
/* >          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U') */
/* >          or the lower triangle (if UPLO='L') of A, including the */
/* >          diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB, N) */
/* >          On entry, the Hermitian positive definite matrix B. */
/* >          If UPLO = 'U', the leading N-by-N upper triangular part of B */
/* >          contains the upper triangular part of the matrix B. */
/* >          If UPLO = 'L', the leading N-by-N lower triangular part of B */
/* >          contains the lower triangular part of the matrix B. */
/* > */
/* >          On exit, if INFO <= N, the part of B containing the matrix is */
/* >          overwritten by the triangular factor U or L from the Cholesky */
/* >          factorization B = U**H*U or B = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,2*N-1). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the blocksize for CHETRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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
/* >          > 0:  CPOTRF or CHEEV returned an error code: */
/* >             <= N:  if INFO = i, CHEEV failed to converge; */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then the leading */
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

/* > \date December 2016 */

/* > \ingroup complexHEeigen */

/*  ===================================================================== */
/* Subroutine */ int chegv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
	 integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer nb, neig;
    extern /* Subroutine */ int cheev_(char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static char trans[1];
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper, wantz;
    extern /* Subroutine */ int chegst_(integer *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cpotrf_(
	    char *, integer *, doublecomplex *, integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.7.0) -- */
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

#line 224 "chegv.f"
    /* Parameter adjustments */
#line 224 "chegv.f"
    a_dim1 = *lda;
#line 224 "chegv.f"
    a_offset = 1 + a_dim1;
#line 224 "chegv.f"
    a -= a_offset;
#line 224 "chegv.f"
    b_dim1 = *ldb;
#line 224 "chegv.f"
    b_offset = 1 + b_dim1;
#line 224 "chegv.f"
    b -= b_offset;
#line 224 "chegv.f"
    --w;
#line 224 "chegv.f"
    --work;
#line 224 "chegv.f"
    --rwork;
#line 224 "chegv.f"

#line 224 "chegv.f"
    /* Function Body */
#line 224 "chegv.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 225 "chegv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 226 "chegv.f"
    lquery = *lwork == -1;

#line 228 "chegv.f"
    *info = 0;
#line 229 "chegv.f"
    if (*itype < 1 || *itype > 3) {
#line 230 "chegv.f"
	*info = -1;
#line 231 "chegv.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 232 "chegv.f"
	*info = -2;
#line 233 "chegv.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 234 "chegv.f"
	*info = -3;
#line 235 "chegv.f"
    } else if (*n < 0) {
#line 236 "chegv.f"
	*info = -4;
#line 237 "chegv.f"
    } else if (*lda < max(1,*n)) {
#line 238 "chegv.f"
	*info = -6;
#line 239 "chegv.f"
    } else if (*ldb < max(1,*n)) {
#line 240 "chegv.f"
	*info = -8;
#line 241 "chegv.f"
    }

#line 243 "chegv.f"
    if (*info == 0) {
#line 244 "chegv.f"
	nb = ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 245 "chegv.f"
	i__1 = 1, i__2 = (nb + 1) * *n;
#line 245 "chegv.f"
	lwkopt = max(i__1,i__2);
#line 246 "chegv.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/* Computing MAX */
#line 248 "chegv.f"
	i__1 = 1, i__2 = (*n << 1) - 1;
#line 248 "chegv.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 249 "chegv.f"
	    *info = -11;
#line 250 "chegv.f"
	}
#line 251 "chegv.f"
    }

#line 253 "chegv.f"
    if (*info != 0) {
#line 254 "chegv.f"
	i__1 = -(*info);
#line 254 "chegv.f"
	xerbla_("CHEGV ", &i__1, (ftnlen)6);
#line 255 "chegv.f"
	return 0;
#line 256 "chegv.f"
    } else if (lquery) {
#line 257 "chegv.f"
	return 0;
#line 258 "chegv.f"
    }

/*     Quick return if possible */

#line 262 "chegv.f"
    if (*n == 0) {
#line 262 "chegv.f"
	return 0;
#line 262 "chegv.f"
    }

/*     Form a Cholesky factorization of B. */

#line 267 "chegv.f"
    cpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 268 "chegv.f"
    if (*info != 0) {
#line 269 "chegv.f"
	*info = *n + *info;
#line 270 "chegv.f"
	return 0;
#line 271 "chegv.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 275 "chegv.f"
    chegst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 276 "chegv.f"
    cheev_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, &rwork[1]
	    , info, (ftnlen)1, (ftnlen)1);

#line 278 "chegv.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 282 "chegv.f"
	neig = *n;
#line 283 "chegv.f"
	if (*info > 0) {
#line 283 "chegv.f"
	    neig = *info - 1;
#line 283 "chegv.f"
	}
#line 285 "chegv.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y */

#line 290 "chegv.f"
	    if (upper) {
#line 291 "chegv.f"
		*(unsigned char *)trans = 'N';
#line 292 "chegv.f"
	    } else {
#line 293 "chegv.f"
		*(unsigned char *)trans = 'C';
#line 294 "chegv.f"
	    }

#line 296 "chegv.f"
	    ctrsm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b1, &b[
		    b_offset], ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);

#line 299 "chegv.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H*y */

#line 304 "chegv.f"
	    if (upper) {
#line 305 "chegv.f"
		*(unsigned char *)trans = 'C';
#line 306 "chegv.f"
	    } else {
#line 307 "chegv.f"
		*(unsigned char *)trans = 'N';
#line 308 "chegv.f"
	    }

#line 310 "chegv.f"
	    ctrmm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b1, &b[
		    b_offset], ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);
#line 312 "chegv.f"
	}
#line 313 "chegv.f"
    }

#line 315 "chegv.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 317 "chegv.f"
    return 0;

/*     End of CHEGV */

} /* chegv_ */

