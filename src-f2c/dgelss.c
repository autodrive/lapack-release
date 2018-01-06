#line 1 "dgelss.f"
/* dgelss.f -- translated by f2c (version 20100827).
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

#line 1 "dgelss.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b46 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;

/* > \brief <b> DGELSS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGELSS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelss.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelss.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelss.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGELSS computes the minimum norm solution to a real linear least */
/* > squares problem: */
/* > */
/* > Minimize 2-norm(| b - A*x |). */
/* > */
/* > using the singular value decomposition (SVD) of A. A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix */
/* > X. */
/* > */
/* > The effective rank of A is determined by treating as zero those */
/* > singular values which are less than RCOND times the largest singular */
/* > value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the first min(m,n) rows of A are overwritten with */
/* >          its right singular vectors, stored rowwise. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, B is overwritten by the N-by-NRHS solution */
/* >          matrix X.  If m >= n and RANK = n, the residual */
/* >          sum-of-squares for the solution in the i-th column is given */
/* >          by the sum of squares of elements n+1:m in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,max(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A in decreasing order. */
/* >          The condition number of A in the 2-norm = S(1)/S(min(m,n)). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          RCOND is used to determine the effective rank of A. */
/* >          Singular values S(i) <= RCOND*S(1) are treated as zero. */
/* >          If RCOND < 0, machine precision is used instead. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* >          RANK is INTEGER */
/* >          The effective rank of A, i.e., the number of singular values */
/* >          which are greater than RCOND*S(1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= 1, and also: */
/* >          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS ) */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  the algorithm for computing the SVD failed to converge; */
/* >                if INFO = i, i off-diagonal elements of an intermediate */
/* >                bidiagonal form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleGEsolve */

/*  ===================================================================== */
/* Subroutine */ int dgelss_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, bl, ie, il, mm;
    static doublereal dum[1], eps, thr, anrm, bnrm;
    static integer itau, lwork_dgebrd__, lwork_dgelqf__, lwork_dgeqrf__, 
	    lwork_dorgbr__, lwork_dormbr__, lwork_dormlq__, lwork_dormqr__;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer iascl, ibscl;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), drscl_(integer *, 
	    doublereal *, doublereal *, integer *);
    static integer chunk;
    static doublereal sfmin;
    static integer minmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxmn, itaup, itauq, mnthr, iwork;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebrd_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static integer bdspac;
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dlacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dbdsqr_(char *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dorgbr_(char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);
    static doublereal bignum;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), dormlq_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 224 "dgelss.f"
    /* Parameter adjustments */
#line 224 "dgelss.f"
    a_dim1 = *lda;
#line 224 "dgelss.f"
    a_offset = 1 + a_dim1;
#line 224 "dgelss.f"
    a -= a_offset;
#line 224 "dgelss.f"
    b_dim1 = *ldb;
#line 224 "dgelss.f"
    b_offset = 1 + b_dim1;
#line 224 "dgelss.f"
    b -= b_offset;
#line 224 "dgelss.f"
    --s;
#line 224 "dgelss.f"
    --work;
#line 224 "dgelss.f"

#line 224 "dgelss.f"
    /* Function Body */
#line 224 "dgelss.f"
    *info = 0;
#line 225 "dgelss.f"
    minmn = min(*m,*n);
#line 226 "dgelss.f"
    maxmn = max(*m,*n);
#line 227 "dgelss.f"
    lquery = *lwork == -1;
#line 228 "dgelss.f"
    if (*m < 0) {
#line 229 "dgelss.f"
	*info = -1;
#line 230 "dgelss.f"
    } else if (*n < 0) {
#line 231 "dgelss.f"
	*info = -2;
#line 232 "dgelss.f"
    } else if (*nrhs < 0) {
#line 233 "dgelss.f"
	*info = -3;
#line 234 "dgelss.f"
    } else if (*lda < max(1,*m)) {
#line 235 "dgelss.f"
	*info = -5;
#line 236 "dgelss.f"
    } else if (*ldb < max(1,maxmn)) {
#line 237 "dgelss.f"
	*info = -7;
#line 238 "dgelss.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 247 "dgelss.f"
    if (*info == 0) {
#line 248 "dgelss.f"
	minwrk = 1;
#line 249 "dgelss.f"
	maxwrk = 1;
#line 250 "dgelss.f"
	if (minmn > 0) {
#line 251 "dgelss.f"
	    mm = *m;
#line 252 "dgelss.f"
	    mnthr = ilaenv_(&c__6, "DGELSS", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
#line 253 "dgelss.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns */

/*              Compute space needed for DGEQRF */
#line 259 "dgelss.f"
		dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 260 "dgelss.f"
		lwork_dgeqrf__ = (integer) dum[0];
/*              Compute space needed for DORMQR */
#line 262 "dgelss.f"
		dormqr_("L", "T", m, nrhs, n, &a[a_offset], lda, dum, &b[
			b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (ftnlen)
			1);
#line 264 "dgelss.f"
		lwork_dormqr__ = (integer) dum[0];
#line 265 "dgelss.f"
		mm = *n;
/* Computing MAX */
#line 266 "dgelss.f"
		i__1 = maxwrk, i__2 = *n + lwork_dgeqrf__;
#line 266 "dgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 267 "dgelss.f"
		i__1 = maxwrk, i__2 = *n + lwork_dormqr__;
#line 267 "dgelss.f"
		maxwrk = max(i__1,i__2);
#line 268 "dgelss.f"
	    }
#line 269 "dgelss.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined */

/*              Compute workspace needed for DBDSQR */

/* Computing MAX */
#line 275 "dgelss.f"
		i__1 = 1, i__2 = *n * 5;
#line 275 "dgelss.f"
		bdspac = max(i__1,i__2);
/*              Compute space needed for DGEBRD */
#line 277 "dgelss.f"
		dgebrd_(&mm, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, 
			&c_n1, info);
#line 279 "dgelss.f"
		lwork_dgebrd__ = (integer) dum[0];
/*              Compute space needed for DORMBR */
#line 281 "dgelss.f"
		dormbr_("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, dum, &
			b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 283 "dgelss.f"
		lwork_dormbr__ = (integer) dum[0];
/*              Compute space needed for DORGBR */
#line 285 "dgelss.f"
		dorgbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			info, (ftnlen)1);
#line 287 "dgelss.f"
		lwork_dorgbr__ = (integer) dum[0];
/*              Compute total workspace needed */
/* Computing MAX */
#line 289 "dgelss.f"
		i__1 = maxwrk, i__2 = *n * 3 + lwork_dgebrd__;
#line 289 "dgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 290 "dgelss.f"
		i__1 = maxwrk, i__2 = *n * 3 + lwork_dormbr__;
#line 290 "dgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 291 "dgelss.f"
		i__1 = maxwrk, i__2 = *n * 3 + lwork_dorgbr__;
#line 291 "dgelss.f"
		maxwrk = max(i__1,i__2);
#line 292 "dgelss.f"
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 293 "dgelss.f"
		i__1 = maxwrk, i__2 = *n * *nrhs;
#line 293 "dgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 294 "dgelss.f"
		i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1,
			i__2);
#line 294 "dgelss.f"
		minwrk = max(i__1,bdspac);
#line 295 "dgelss.f"
		maxwrk = max(minwrk,maxwrk);
#line 296 "dgelss.f"
	    }
#line 297 "dgelss.f"
	    if (*n > *m) {

/*              Compute workspace needed for DBDSQR */

/* Computing MAX */
#line 301 "dgelss.f"
		i__1 = 1, i__2 = *m * 5;
#line 301 "dgelss.f"
		bdspac = max(i__1,i__2);
/* Computing MAX */
#line 302 "dgelss.f"
		i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *n, i__1 = max(i__1,
			i__2);
#line 302 "dgelss.f"
		minwrk = max(i__1,bdspac);
#line 303 "dgelss.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                 than rows */

/*                 Compute space needed for DGELQF */
#line 309 "dgelss.f"
		    dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 311 "dgelss.f"
		    lwork_dgelqf__ = (integer) dum[0];
/*                 Compute space needed for DGEBRD */
#line 313 "dgelss.f"
		    dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, 
			    dum, &c_n1, info);
#line 315 "dgelss.f"
		    lwork_dgebrd__ = (integer) dum[0];
/*                 Compute space needed for DORMBR */
#line 317 "dgelss.f"
		    dormbr_("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 319 "dgelss.f"
		    lwork_dormbr__ = (integer) dum[0];
/*                 Compute space needed for DORGBR */
#line 321 "dgelss.f"
		    dorgbr_("P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 323 "dgelss.f"
		    lwork_dorgbr__ = (integer) dum[0];
/*                 Compute space needed for DORMLQ */
#line 325 "dgelss.f"
		    dormlq_("L", "T", n, nrhs, m, &a[a_offset], lda, dum, &b[
			    b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1);
#line 327 "dgelss.f"
		    lwork_dormlq__ = (integer) dum[0];
/*                 Compute total workspace needed */
#line 329 "dgelss.f"
		    maxwrk = *m + lwork_dgelqf__;
/* Computing MAX */
#line 330 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + 
			    lwork_dgebrd__;
#line 330 "dgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 331 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + 
			    lwork_dormbr__;
#line 331 "dgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 332 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + 
			    lwork_dorgbr__;
#line 332 "dgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 333 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + *m + bdspac;
#line 333 "dgelss.f"
		    maxwrk = max(i__1,i__2);
#line 334 "dgelss.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 335 "dgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 335 "dgelss.f"
			maxwrk = max(i__1,i__2);
#line 336 "dgelss.f"
		    } else {
/* Computing MAX */
#line 337 "dgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 337 "dgelss.f"
			maxwrk = max(i__1,i__2);
#line 338 "dgelss.f"
		    }
/* Computing MAX */
#line 339 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m + lwork_dormlq__;
#line 339 "dgelss.f"
		    maxwrk = max(i__1,i__2);
#line 340 "dgelss.f"
		} else {

/*                 Path 2 - underdetermined */

/*                 Compute space needed for DGEBRD */
#line 345 "dgelss.f"
		    dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, 
			    dum, &c_n1, info);
#line 347 "dgelss.f"
		    lwork_dgebrd__ = (integer) dum[0];
/*                 Compute space needed for DORMBR */
#line 349 "dgelss.f"
		    dormbr_("Q", "L", "T", m, nrhs, m, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 351 "dgelss.f"
		    lwork_dormbr__ = (integer) dum[0];
/*                 Compute space needed for DORGBR */
#line 353 "dgelss.f"
		    dorgbr_("P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 355 "dgelss.f"
		    lwork_dorgbr__ = (integer) dum[0];
#line 356 "dgelss.f"
		    maxwrk = *m * 3 + lwork_dgebrd__;
/* Computing MAX */
#line 357 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + lwork_dormbr__;
#line 357 "dgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 358 "dgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + lwork_dorgbr__;
#line 358 "dgelss.f"
		    maxwrk = max(i__1,i__2);
#line 359 "dgelss.f"
		    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 360 "dgelss.f"
		    i__1 = maxwrk, i__2 = *n * *nrhs;
#line 360 "dgelss.f"
		    maxwrk = max(i__1,i__2);
#line 361 "dgelss.f"
		}
#line 362 "dgelss.f"
	    }
#line 363 "dgelss.f"
	    maxwrk = max(minwrk,maxwrk);
#line 364 "dgelss.f"
	}
#line 365 "dgelss.f"
	work[1] = (doublereal) maxwrk;

#line 367 "dgelss.f"
	if (*lwork < minwrk && ! lquery) {
#line 367 "dgelss.f"
	    *info = -12;
#line 367 "dgelss.f"
	}
#line 369 "dgelss.f"
    }

#line 371 "dgelss.f"
    if (*info != 0) {
#line 372 "dgelss.f"
	i__1 = -(*info);
#line 372 "dgelss.f"
	xerbla_("DGELSS", &i__1, (ftnlen)6);
#line 373 "dgelss.f"
	return 0;
#line 374 "dgelss.f"
    } else if (lquery) {
#line 375 "dgelss.f"
	return 0;
#line 376 "dgelss.f"
    }

/*     Quick return if possible */

#line 380 "dgelss.f"
    if (*m == 0 || *n == 0) {
#line 381 "dgelss.f"
	*rank = 0;
#line 382 "dgelss.f"
	return 0;
#line 383 "dgelss.f"
    }

/*     Get machine parameters */

#line 387 "dgelss.f"
    eps = dlamch_("P", (ftnlen)1);
#line 388 "dgelss.f"
    sfmin = dlamch_("S", (ftnlen)1);
#line 389 "dgelss.f"
    smlnum = sfmin / eps;
#line 390 "dgelss.f"
    bignum = 1. / smlnum;
#line 391 "dgelss.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 395 "dgelss.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 396 "dgelss.f"
    iascl = 0;
#line 397 "dgelss.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 401 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 402 "dgelss.f"
	iascl = 1;
#line 403 "dgelss.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 407 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 408 "dgelss.f"
	iascl = 2;
#line 409 "dgelss.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 413 "dgelss.f"
	i__1 = max(*m,*n);
#line 413 "dgelss.f"
	dlaset_("F", &i__1, nrhs, &c_b46, &c_b46, &b[b_offset], ldb, (ftnlen)
		1);
#line 414 "dgelss.f"
	dlaset_("F", &minmn, &c__1, &c_b46, &c_b46, &s[1], &minmn, (ftnlen)1);
#line 415 "dgelss.f"
	*rank = 0;
#line 416 "dgelss.f"
	goto L70;
#line 417 "dgelss.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 421 "dgelss.f"
    bnrm = dlange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 422 "dgelss.f"
    ibscl = 0;
#line 423 "dgelss.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 427 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 428 "dgelss.f"
	ibscl = 1;
#line 429 "dgelss.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 433 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 434 "dgelss.f"
	ibscl = 2;
#line 435 "dgelss.f"
    }

/*     Overdetermined case */

#line 439 "dgelss.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

#line 443 "dgelss.f"
	mm = *m;
#line 444 "dgelss.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

#line 448 "dgelss.f"
	    mm = *n;
#line 449 "dgelss.f"
	    itau = 1;
#line 450 "dgelss.f"
	    iwork = itau + *n;

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 455 "dgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 455 "dgelss.f"
	    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1,
		     info);

/*           Multiply B by transpose(Q) */
/*           (Workspace: need N+NRHS, prefer N+NRHS*NB) */

#line 461 "dgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 461 "dgelss.f"
	    dormqr_("L", "T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R */

#line 466 "dgelss.f"
	    if (*n > 1) {
#line 466 "dgelss.f"
		i__1 = *n - 1;
#line 466 "dgelss.f"
		i__2 = *n - 1;
#line 466 "dgelss.f"
		dlaset_("L", &i__1, &i__2, &c_b46, &c_b46, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 466 "dgelss.f"
	    }
#line 468 "dgelss.f"
	}

#line 470 "dgelss.f"
	ie = 1;
#line 471 "dgelss.f"
	itauq = ie + *n;
#line 472 "dgelss.f"
	itaup = itauq + *n;
#line 473 "dgelss.f"
	iwork = itaup + *n;

/*        Bidiagonalize R in A */
/*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

#line 478 "dgelss.f"
	i__1 = *lwork - iwork + 1;
#line 478 "dgelss.f"
	dgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R */
/*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

#line 485 "dgelss.f"
	i__1 = *lwork - iwork + 1;
#line 485 "dgelss.f"
	dormbr_("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in A */
/*        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 491 "dgelss.f"
	i__1 = *lwork - iwork + 1;
#line 491 "dgelss.f"
	dorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &
		i__1, info, (ftnlen)1);
#line 493 "dgelss.f"
	iwork = ie + *n;

/*        Perform bidiagonal QR iteration */
/*          multiply B by transpose of left singular vectors */
/*          compute right singular vectors in A */
/*        (Workspace: need BDSPAC) */

#line 500 "dgelss.f"
	dbdsqr_("U", n, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda, 
		dum, &c__1, &b[b_offset], ldb, &work[iwork], info, (ftnlen)1);
#line 502 "dgelss.f"
	if (*info != 0) {
#line 502 "dgelss.f"
	    goto L70;
#line 502 "dgelss.f"
	}

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 507 "dgelss.f"
	d__1 = *rcond * s[1];
#line 507 "dgelss.f"
	thr = max(d__1,sfmin);
#line 508 "dgelss.f"
	if (*rcond < 0.) {
/* Computing MAX */
#line 508 "dgelss.f"
	    d__1 = eps * s[1];
#line 508 "dgelss.f"
	    thr = max(d__1,sfmin);
#line 508 "dgelss.f"
	}
#line 510 "dgelss.f"
	*rank = 0;
#line 511 "dgelss.f"
	i__1 = *n;
#line 511 "dgelss.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 512 "dgelss.f"
	    if (s[i__] > thr) {
#line 513 "dgelss.f"
		drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 514 "dgelss.f"
		++(*rank);
#line 515 "dgelss.f"
	    } else {
#line 516 "dgelss.f"
		dlaset_("F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1], 
			ldb, (ftnlen)1);
#line 517 "dgelss.f"
	    }
#line 518 "dgelss.f"
/* L10: */
#line 518 "dgelss.f"
	}

/*        Multiply B by right singular vectors */
/*        (Workspace: need N, prefer N*NRHS) */

#line 523 "dgelss.f"
	if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 524 "dgelss.f"
	    dgemm_("T", "N", n, nrhs, n, &c_b79, &a[a_offset], lda, &b[
		    b_offset], ldb, &c_b46, &work[1], ldb, (ftnlen)1, (ftnlen)
		    1);
#line 526 "dgelss.f"
	    dlacpy_("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1)
		    ;
#line 527 "dgelss.f"
	} else if (*nrhs > 1) {
#line 528 "dgelss.f"
	    chunk = *lwork / *n;
#line 529 "dgelss.f"
	    i__1 = *nrhs;
#line 529 "dgelss.f"
	    i__2 = chunk;
#line 529 "dgelss.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 530 "dgelss.f"
		i__3 = *nrhs - i__ + 1;
#line 530 "dgelss.f"
		bl = min(i__3,chunk);
#line 531 "dgelss.f"
		dgemm_("T", "N", n, &bl, n, &c_b79, &a[a_offset], lda, &b[i__ 
			* b_dim1 + 1], ldb, &c_b46, &work[1], n, (ftnlen)1, (
			ftnlen)1);
#line 533 "dgelss.f"
		dlacpy_("G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (
			ftnlen)1);
#line 534 "dgelss.f"
/* L20: */
#line 534 "dgelss.f"
	    }
#line 535 "dgelss.f"
	} else {
#line 536 "dgelss.f"
	    dgemv_("T", n, n, &c_b79, &a[a_offset], lda, &b[b_offset], &c__1, 
		    &c_b46, &work[1], &c__1, (ftnlen)1);
#line 537 "dgelss.f"
	    dcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 538 "dgelss.f"
	}

#line 540 "dgelss.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 540 "dgelss.f"
	i__2 = *m, i__1 = (*m << 1) - 4, i__2 = max(i__2,i__1), i__2 = max(
		i__2,*nrhs), i__1 = *n - *m * 3;
#line 540 "dgelss.f"
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__2,i__1)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm */

#line 546 "dgelss.f"
	    ldwork = *m;
/* Computing MAX */
/* Computing MAX */
#line 547 "dgelss.f"
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 547 "dgelss.f"
	    i__2 = (*m << 2) + *m * *lda + max(i__3,i__4), i__1 = *m * *lda + 
		    *m + *m * *nrhs;
#line 547 "dgelss.f"
	    if (*lwork >= max(i__2,i__1)) {
#line 547 "dgelss.f"
		ldwork = *lda;
#line 547 "dgelss.f"
	    }
#line 549 "dgelss.f"
	    itau = 1;
#line 550 "dgelss.f"
	    iwork = *m + 1;

/*        Compute A=L*Q */
/*        (Workspace: need 2*M, prefer M+M*NB) */

#line 555 "dgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 555 "dgelss.f"
	    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
		     info);
#line 557 "dgelss.f"
	    il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */

#line 561 "dgelss.f"
	    dlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 562 "dgelss.f"
	    i__2 = *m - 1;
#line 562 "dgelss.f"
	    i__1 = *m - 1;
#line 562 "dgelss.f"
	    dlaset_("U", &i__2, &i__1, &c_b46, &c_b46, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 564 "dgelss.f"
	    ie = il + ldwork * *m;
#line 565 "dgelss.f"
	    itauq = ie + *m;
#line 566 "dgelss.f"
	    itaup = itauq + *m;
#line 567 "dgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL) */
/*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

#line 572 "dgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 572 "dgelss.f"
	    dgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L */
/*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

#line 579 "dgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 579 "dgelss.f"
	    dormbr_("Q", "L", "T", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[iwork], &i__2, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) */
/*        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB) */

#line 586 "dgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 586 "dgelss.f"
	    dorgbr_("P", m, m, m, &work[il], &ldwork, &work[itaup], &work[
		    iwork], &i__2, info, (ftnlen)1);
#line 588 "dgelss.f"
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of L in WORK(IL) and */
/*           multiplying B by transpose of left singular vectors */
/*        (Workspace: need M*M+M+BDSPAC) */

#line 595 "dgelss.f"
	    dbdsqr_("U", m, m, &c__0, nrhs, &s[1], &work[ie], &work[il], &
		    ldwork, &a[a_offset], lda, &b[b_offset], ldb, &work[iwork]
		    , info, (ftnlen)1);
#line 597 "dgelss.f"
	    if (*info != 0) {
#line 597 "dgelss.f"
		goto L70;
#line 597 "dgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 602 "dgelss.f"
	    d__1 = *rcond * s[1];
#line 602 "dgelss.f"
	    thr = max(d__1,sfmin);
#line 603 "dgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 603 "dgelss.f"
		d__1 = eps * s[1];
#line 603 "dgelss.f"
		thr = max(d__1,sfmin);
#line 603 "dgelss.f"
	    }
#line 605 "dgelss.f"
	    *rank = 0;
#line 606 "dgelss.f"
	    i__2 = *m;
#line 606 "dgelss.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 607 "dgelss.f"
		if (s[i__] > thr) {
#line 608 "dgelss.f"
		    drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 609 "dgelss.f"
		    ++(*rank);
#line 610 "dgelss.f"
		} else {
#line 611 "dgelss.f"
		    dlaset_("F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1]
			    , ldb, (ftnlen)1);
#line 612 "dgelss.f"
		}
#line 613 "dgelss.f"
/* L30: */
#line 613 "dgelss.f"
	    }
#line 614 "dgelss.f"
	    iwork = ie;

/*        Multiply B by right singular vectors of L in WORK(IL) */
/*        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS) */

#line 619 "dgelss.f"
	    if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
#line 620 "dgelss.f"
		dgemm_("T", "N", m, nrhs, m, &c_b79, &work[il], &ldwork, &b[
			b_offset], ldb, &c_b46, &work[iwork], ldb, (ftnlen)1, 
			(ftnlen)1);
#line 622 "dgelss.f"
		dlacpy_("G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 623 "dgelss.f"
	    } else if (*nrhs > 1) {
#line 624 "dgelss.f"
		chunk = (*lwork - iwork + 1) / *m;
#line 625 "dgelss.f"
		i__2 = *nrhs;
#line 625 "dgelss.f"
		i__1 = chunk;
#line 625 "dgelss.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 626 "dgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 626 "dgelss.f"
		    bl = min(i__3,chunk);
#line 627 "dgelss.f"
		    dgemm_("T", "N", m, &bl, m, &c_b79, &work[il], &ldwork, &
			    b[i__ * b_dim1 + 1], ldb, &c_b46, &work[iwork], m,
			     (ftnlen)1, (ftnlen)1);
#line 629 "dgelss.f"
		    dlacpy_("G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1]
			    , ldb, (ftnlen)1);
#line 631 "dgelss.f"
/* L40: */
#line 631 "dgelss.f"
		}
#line 632 "dgelss.f"
	    } else {
#line 633 "dgelss.f"
		dgemv_("T", m, m, &c_b79, &work[il], &ldwork, &b[b_dim1 + 1], 
			&c__1, &c_b46, &work[iwork], &c__1, (ftnlen)1);
#line 635 "dgelss.f"
		dcopy_(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
#line 636 "dgelss.f"
	    }

/*        Zero out below first M rows of B */

#line 640 "dgelss.f"
	    i__1 = *n - *m;
#line 640 "dgelss.f"
	    dlaset_("F", &i__1, nrhs, &c_b46, &c_b46, &b[*m + 1 + b_dim1], 
		    ldb, (ftnlen)1);
#line 641 "dgelss.f"
	    iwork = itau + *m;

/*        Multiply transpose(Q) by B */
/*        (Workspace: need M+NRHS, prefer M+NRHS*NB) */

#line 646 "dgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 646 "dgelss.f"
	    dormlq_("L", "T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 649 "dgelss.f"
	} else {

/*        Path 2 - remaining underdetermined cases */

#line 653 "dgelss.f"
	    ie = 1;
#line 654 "dgelss.f"
	    itauq = ie + *m;
#line 655 "dgelss.f"
	    itaup = itauq + *m;
#line 656 "dgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize A */
/*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 661 "dgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 661 "dgelss.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors */
/*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

#line 668 "dgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 668 "dgelss.f"
	    dormbr_("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors in A */
/*        (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 674 "dgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 674 "dgelss.f"
	    dorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
		    iwork], &i__1, info, (ftnlen)1);
#line 676 "dgelss.f"
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of A in A and */
/*           multiplying B by transpose of left singular vectors */
/*        (Workspace: need BDSPAC) */

#line 683 "dgelss.f"
	    dbdsqr_("L", m, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], 
		    lda, dum, &c__1, &b[b_offset], ldb, &work[iwork], info, (
		    ftnlen)1);
#line 685 "dgelss.f"
	    if (*info != 0) {
#line 685 "dgelss.f"
		goto L70;
#line 685 "dgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 690 "dgelss.f"
	    d__1 = *rcond * s[1];
#line 690 "dgelss.f"
	    thr = max(d__1,sfmin);
#line 691 "dgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 691 "dgelss.f"
		d__1 = eps * s[1];
#line 691 "dgelss.f"
		thr = max(d__1,sfmin);
#line 691 "dgelss.f"
	    }
#line 693 "dgelss.f"
	    *rank = 0;
#line 694 "dgelss.f"
	    i__1 = *m;
#line 694 "dgelss.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 695 "dgelss.f"
		if (s[i__] > thr) {
#line 696 "dgelss.f"
		    drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 697 "dgelss.f"
		    ++(*rank);
#line 698 "dgelss.f"
		} else {
#line 699 "dgelss.f"
		    dlaset_("F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1]
			    , ldb, (ftnlen)1);
#line 700 "dgelss.f"
		}
#line 701 "dgelss.f"
/* L50: */
#line 701 "dgelss.f"
	    }

/*        Multiply B by right singular vectors of A */
/*        (Workspace: need N, prefer N*NRHS) */

#line 706 "dgelss.f"
	    if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 707 "dgelss.f"
		dgemm_("T", "N", n, nrhs, m, &c_b79, &a[a_offset], lda, &b[
			b_offset], ldb, &c_b46, &work[1], ldb, (ftnlen)1, (
			ftnlen)1);
#line 709 "dgelss.f"
		dlacpy_("F", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 710 "dgelss.f"
	    } else if (*nrhs > 1) {
#line 711 "dgelss.f"
		chunk = *lwork / *n;
#line 712 "dgelss.f"
		i__1 = *nrhs;
#line 712 "dgelss.f"
		i__2 = chunk;
#line 712 "dgelss.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 713 "dgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 713 "dgelss.f"
		    bl = min(i__3,chunk);
#line 714 "dgelss.f"
		    dgemm_("T", "N", n, &bl, m, &c_b79, &a[a_offset], lda, &b[
			    i__ * b_dim1 + 1], ldb, &c_b46, &work[1], n, (
			    ftnlen)1, (ftnlen)1);
#line 716 "dgelss.f"
		    dlacpy_("F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], 
			    ldb, (ftnlen)1);
#line 717 "dgelss.f"
/* L60: */
#line 717 "dgelss.f"
		}
#line 718 "dgelss.f"
	    } else {
#line 719 "dgelss.f"
		dgemv_("T", m, n, &c_b79, &a[a_offset], lda, &b[b_offset], &
			c__1, &c_b46, &work[1], &c__1, (ftnlen)1);
#line 720 "dgelss.f"
		dcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 721 "dgelss.f"
	    }
#line 722 "dgelss.f"
	}
#line 722 "dgelss.f"
    }

/*     Undo scaling */

#line 726 "dgelss.f"
    if (iascl == 1) {
#line 727 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 728 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 730 "dgelss.f"
    } else if (iascl == 2) {
#line 731 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 732 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 734 "dgelss.f"
    }
#line 735 "dgelss.f"
    if (ibscl == 1) {
#line 736 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 737 "dgelss.f"
    } else if (ibscl == 2) {
#line 738 "dgelss.f"
	dlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 739 "dgelss.f"
    }

#line 741 "dgelss.f"
L70:
#line 742 "dgelss.f"
    work[1] = (doublereal) maxwrk;
#line 743 "dgelss.f"
    return 0;

/*     End of DGELSS */

} /* dgelss_ */

