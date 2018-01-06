#line 1 "sgelss.f"
/* sgelss.f -- translated by f2c (version 20100827).
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

#line 1 "sgelss.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b50 = 0.;
static doublereal c_b84 = 1.;

/* > \brief <b> SGELSS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGELSS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelss.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelss.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelss.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGELSS computes the minimum norm solution to a real linear least */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          S is REAL array, dimension (min(M,N)) */
/* >          The singular values of A in decreasing order. */
/* >          The condition number of A in the 2-norm = S(1)/S(min(m,n)). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup realGEsolve */

/*  ===================================================================== */
/* Subroutine */ int sgelss_(integer *m, integer *n, integer *nrhs, 
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
    static integer itau, lwork_sgebrd__, lwork_sgeqrf__, lwork_sorgbr__, 
	    lwork_sormbr__, lwork_sormlq__, iascl, ibscl, lwork_sormqr__, 
	    chunk;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal sfmin;
    static integer minmn, maxmn;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer itaup, itauq;
    extern /* Subroutine */ int srscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer mnthr, iwork;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slabad_(doublereal *, doublereal *);
    static integer bdspac;
    extern /* Subroutine */ int sgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     sgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), slacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    sbdsqr_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), sorgbr_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int sormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int sormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 223 "sgelss.f"
    /* Parameter adjustments */
#line 223 "sgelss.f"
    a_dim1 = *lda;
#line 223 "sgelss.f"
    a_offset = 1 + a_dim1;
#line 223 "sgelss.f"
    a -= a_offset;
#line 223 "sgelss.f"
    b_dim1 = *ldb;
#line 223 "sgelss.f"
    b_offset = 1 + b_dim1;
#line 223 "sgelss.f"
    b -= b_offset;
#line 223 "sgelss.f"
    --s;
#line 223 "sgelss.f"
    --work;
#line 223 "sgelss.f"

#line 223 "sgelss.f"
    /* Function Body */
#line 223 "sgelss.f"
    *info = 0;
#line 224 "sgelss.f"
    minmn = min(*m,*n);
#line 225 "sgelss.f"
    maxmn = max(*m,*n);
#line 226 "sgelss.f"
    lquery = *lwork == -1;
#line 227 "sgelss.f"
    if (*m < 0) {
#line 228 "sgelss.f"
	*info = -1;
#line 229 "sgelss.f"
    } else if (*n < 0) {
#line 230 "sgelss.f"
	*info = -2;
#line 231 "sgelss.f"
    } else if (*nrhs < 0) {
#line 232 "sgelss.f"
	*info = -3;
#line 233 "sgelss.f"
    } else if (*lda < max(1,*m)) {
#line 234 "sgelss.f"
	*info = -5;
#line 235 "sgelss.f"
    } else if (*ldb < max(1,maxmn)) {
#line 236 "sgelss.f"
	*info = -7;
#line 237 "sgelss.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 246 "sgelss.f"
    if (*info == 0) {
#line 247 "sgelss.f"
	minwrk = 1;
#line 248 "sgelss.f"
	maxwrk = 1;
#line 249 "sgelss.f"
	if (minmn > 0) {
#line 250 "sgelss.f"
	    mm = *m;
#line 251 "sgelss.f"
	    mnthr = ilaenv_(&c__6, "SGELSS", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
#line 252 "sgelss.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns */

/*              Compute space needed for SGEQRF */
#line 258 "sgelss.f"
		sgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 259 "sgelss.f"
		lwork_sgeqrf__ = (integer) dum[0];
/*              Compute space needed for SORMQR */
#line 261 "sgelss.f"
		sormqr_("L", "T", m, nrhs, n, &a[a_offset], lda, dum, &b[
			b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (ftnlen)
			1);
#line 263 "sgelss.f"
		lwork_sormqr__ = (integer) dum[0];
#line 264 "sgelss.f"
		mm = *n;
/* Computing MAX */
#line 265 "sgelss.f"
		i__1 = maxwrk, i__2 = *n + lwork_sgeqrf__;
#line 265 "sgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 266 "sgelss.f"
		i__1 = maxwrk, i__2 = *n + lwork_sormqr__;
#line 266 "sgelss.f"
		maxwrk = max(i__1,i__2);
#line 267 "sgelss.f"
	    }
#line 268 "sgelss.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined */

/*              Compute workspace needed for SBDSQR */

/* Computing MAX */
#line 274 "sgelss.f"
		i__1 = 1, i__2 = *n * 5;
#line 274 "sgelss.f"
		bdspac = max(i__1,i__2);
/*              Compute space needed for SGEBRD */
#line 276 "sgelss.f"
		sgebrd_(&mm, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, 
			&c_n1, info);
#line 278 "sgelss.f"
		lwork_sgebrd__ = (integer) dum[0];
/*              Compute space needed for SORMBR */
#line 280 "sgelss.f"
		sormbr_("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, dum, &
			b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 282 "sgelss.f"
		lwork_sormbr__ = (integer) dum[0];
/*              Compute space needed for SORGBR */
#line 284 "sgelss.f"
		sorgbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			info, (ftnlen)1);
#line 286 "sgelss.f"
		lwork_sorgbr__ = (integer) dum[0];
/*              Compute total workspace needed */
/* Computing MAX */
#line 288 "sgelss.f"
		i__1 = maxwrk, i__2 = *n * 3 + lwork_sgebrd__;
#line 288 "sgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 289 "sgelss.f"
		i__1 = maxwrk, i__2 = *n * 3 + lwork_sormbr__;
#line 289 "sgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 290 "sgelss.f"
		i__1 = maxwrk, i__2 = *n * 3 + lwork_sorgbr__;
#line 290 "sgelss.f"
		maxwrk = max(i__1,i__2);
#line 291 "sgelss.f"
		maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 292 "sgelss.f"
		i__1 = maxwrk, i__2 = *n * *nrhs;
#line 292 "sgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 293 "sgelss.f"
		i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1,
			i__2);
#line 293 "sgelss.f"
		minwrk = max(i__1,bdspac);
#line 294 "sgelss.f"
		maxwrk = max(minwrk,maxwrk);
#line 295 "sgelss.f"
	    }
#line 296 "sgelss.f"
	    if (*n > *m) {

/*              Compute workspace needed for SBDSQR */

/* Computing MAX */
#line 300 "sgelss.f"
		i__1 = 1, i__2 = *m * 5;
#line 300 "sgelss.f"
		bdspac = max(i__1,i__2);
/* Computing MAX */
#line 301 "sgelss.f"
		i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *n, i__1 = max(i__1,
			i__2);
#line 301 "sgelss.f"
		minwrk = max(i__1,bdspac);
#line 302 "sgelss.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                 than rows */

/*                 Compute space needed for SGEBRD */
#line 308 "sgelss.f"
		    sgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, 
			    dum, &c_n1, info);
#line 310 "sgelss.f"
		    lwork_sgebrd__ = (integer) dum[0];
/*                 Compute space needed for SORMBR */
#line 312 "sgelss.f"
		    sormbr_("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 314 "sgelss.f"
		    lwork_sormbr__ = (integer) dum[0];
/*                 Compute space needed for SORGBR */
#line 316 "sgelss.f"
		    sorgbr_("P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 318 "sgelss.f"
		    lwork_sorgbr__ = (integer) dum[0];
/*                 Compute space needed for SORMLQ */
#line 320 "sgelss.f"
		    sormlq_("L", "T", n, nrhs, m, &a[a_offset], lda, dum, &b[
			    b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1);
#line 322 "sgelss.f"
		    lwork_sormlq__ = (integer) dum[0];
/*                 Compute total workspace needed */
#line 324 "sgelss.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 326 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + 
			    lwork_sgebrd__;
#line 326 "sgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 327 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + 
			    lwork_sormbr__;
#line 327 "sgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 328 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + 
			    lwork_sorgbr__;
#line 328 "sgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 329 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m * *m + *m + bdspac;
#line 329 "sgelss.f"
		    maxwrk = max(i__1,i__2);
#line 330 "sgelss.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 331 "sgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 331 "sgelss.f"
			maxwrk = max(i__1,i__2);
#line 332 "sgelss.f"
		    } else {
/* Computing MAX */
#line 333 "sgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 333 "sgelss.f"
			maxwrk = max(i__1,i__2);
#line 334 "sgelss.f"
		    }
/* Computing MAX */
#line 335 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m + lwork_sormlq__;
#line 335 "sgelss.f"
		    maxwrk = max(i__1,i__2);
#line 336 "sgelss.f"
		} else {

/*                 Path 2 - underdetermined */

/*                 Compute space needed for SGEBRD */
#line 341 "sgelss.f"
		    sgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, 
			    dum, &c_n1, info);
#line 343 "sgelss.f"
		    lwork_sgebrd__ = (integer) dum[0];
/*                 Compute space needed for SORMBR */
#line 345 "sgelss.f"
		    sormbr_("Q", "L", "T", m, nrhs, m, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 347 "sgelss.f"
		    lwork_sormbr__ = (integer) dum[0];
/*                 Compute space needed for SORGBR */
#line 349 "sgelss.f"
		    sorgbr_("P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 351 "sgelss.f"
		    lwork_sorgbr__ = (integer) dum[0];
#line 352 "sgelss.f"
		    maxwrk = *m * 3 + lwork_sgebrd__;
/* Computing MAX */
#line 353 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + lwork_sormbr__;
#line 353 "sgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 354 "sgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + lwork_sorgbr__;
#line 354 "sgelss.f"
		    maxwrk = max(i__1,i__2);
#line 355 "sgelss.f"
		    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
#line 356 "sgelss.f"
		    i__1 = maxwrk, i__2 = *n * *nrhs;
#line 356 "sgelss.f"
		    maxwrk = max(i__1,i__2);
#line 357 "sgelss.f"
		}
#line 358 "sgelss.f"
	    }
#line 359 "sgelss.f"
	    maxwrk = max(minwrk,maxwrk);
#line 360 "sgelss.f"
	}
#line 361 "sgelss.f"
	work[1] = (doublereal) maxwrk;

#line 363 "sgelss.f"
	if (*lwork < minwrk && ! lquery) {
#line 363 "sgelss.f"
	    *info = -12;
#line 363 "sgelss.f"
	}
#line 365 "sgelss.f"
    }

#line 367 "sgelss.f"
    if (*info != 0) {
#line 368 "sgelss.f"
	i__1 = -(*info);
#line 368 "sgelss.f"
	xerbla_("SGELSS", &i__1, (ftnlen)6);
#line 369 "sgelss.f"
	return 0;
#line 370 "sgelss.f"
    } else if (lquery) {
#line 371 "sgelss.f"
	return 0;
#line 372 "sgelss.f"
    }

/*     Quick return if possible */

#line 376 "sgelss.f"
    if (*m == 0 || *n == 0) {
#line 377 "sgelss.f"
	*rank = 0;
#line 378 "sgelss.f"
	return 0;
#line 379 "sgelss.f"
    }

/*     Get machine parameters */

#line 383 "sgelss.f"
    eps = slamch_("P", (ftnlen)1);
#line 384 "sgelss.f"
    sfmin = slamch_("S", (ftnlen)1);
#line 385 "sgelss.f"
    smlnum = sfmin / eps;
#line 386 "sgelss.f"
    bignum = 1. / smlnum;
#line 387 "sgelss.f"
    slabad_(&smlnum, &bignum);

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 391 "sgelss.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 392 "sgelss.f"
    iascl = 0;
#line 393 "sgelss.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 397 "sgelss.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 398 "sgelss.f"
	iascl = 1;
#line 399 "sgelss.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 403 "sgelss.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 404 "sgelss.f"
	iascl = 2;
#line 405 "sgelss.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 409 "sgelss.f"
	i__1 = max(*m,*n);
#line 409 "sgelss.f"
	slaset_("F", &i__1, nrhs, &c_b50, &c_b50, &b[b_offset], ldb, (ftnlen)
		1);
#line 410 "sgelss.f"
	slaset_("F", &minmn, &c__1, &c_b50, &c_b50, &s[1], &c__1, (ftnlen)1);
#line 411 "sgelss.f"
	*rank = 0;
#line 412 "sgelss.f"
	goto L70;
#line 413 "sgelss.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 417 "sgelss.f"
    bnrm = slange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 418 "sgelss.f"
    ibscl = 0;
#line 419 "sgelss.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 423 "sgelss.f"
	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 424 "sgelss.f"
	ibscl = 1;
#line 425 "sgelss.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 429 "sgelss.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 430 "sgelss.f"
	ibscl = 2;
#line 431 "sgelss.f"
    }

/*     Overdetermined case */

#line 435 "sgelss.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

#line 439 "sgelss.f"
	mm = *m;
#line 440 "sgelss.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

#line 444 "sgelss.f"
	    mm = *n;
#line 445 "sgelss.f"
	    itau = 1;
#line 446 "sgelss.f"
	    iwork = itau + *n;

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 451 "sgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 451 "sgelss.f"
	    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1,
		     info);

/*           Multiply B by transpose(Q) */
/*           (Workspace: need N+NRHS, prefer N+NRHS*NB) */

#line 457 "sgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 457 "sgelss.f"
	    sormqr_("L", "T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R */

#line 462 "sgelss.f"
	    if (*n > 1) {
#line 462 "sgelss.f"
		i__1 = *n - 1;
#line 462 "sgelss.f"
		i__2 = *n - 1;
#line 462 "sgelss.f"
		slaset_("L", &i__1, &i__2, &c_b50, &c_b50, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 462 "sgelss.f"
	    }
#line 464 "sgelss.f"
	}

#line 466 "sgelss.f"
	ie = 1;
#line 467 "sgelss.f"
	itauq = ie + *n;
#line 468 "sgelss.f"
	itaup = itauq + *n;
#line 469 "sgelss.f"
	iwork = itaup + *n;

/*        Bidiagonalize R in A */
/*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

#line 474 "sgelss.f"
	i__1 = *lwork - iwork + 1;
#line 474 "sgelss.f"
	sgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R */
/*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

#line 481 "sgelss.f"
	i__1 = *lwork - iwork + 1;
#line 481 "sgelss.f"
	sormbr_("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in A */
/*        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

#line 487 "sgelss.f"
	i__1 = *lwork - iwork + 1;
#line 487 "sgelss.f"
	sorgbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &
		i__1, info, (ftnlen)1);
#line 489 "sgelss.f"
	iwork = ie + *n;

/*        Perform bidiagonal QR iteration */
/*          multiply B by transpose of left singular vectors */
/*          compute right singular vectors in A */
/*        (Workspace: need BDSPAC) */

#line 496 "sgelss.f"
	sbdsqr_("U", n, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda, 
		dum, &c__1, &b[b_offset], ldb, &work[iwork], info, (ftnlen)1);
#line 498 "sgelss.f"
	if (*info != 0) {
#line 498 "sgelss.f"
	    goto L70;
#line 498 "sgelss.f"
	}

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 503 "sgelss.f"
	d__1 = *rcond * s[1];
#line 503 "sgelss.f"
	thr = max(d__1,sfmin);
#line 504 "sgelss.f"
	if (*rcond < 0.) {
/* Computing MAX */
#line 504 "sgelss.f"
	    d__1 = eps * s[1];
#line 504 "sgelss.f"
	    thr = max(d__1,sfmin);
#line 504 "sgelss.f"
	}
#line 506 "sgelss.f"
	*rank = 0;
#line 507 "sgelss.f"
	i__1 = *n;
#line 507 "sgelss.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 508 "sgelss.f"
	    if (s[i__] > thr) {
#line 509 "sgelss.f"
		srscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 510 "sgelss.f"
		++(*rank);
#line 511 "sgelss.f"
	    } else {
#line 512 "sgelss.f"
		slaset_("F", &c__1, nrhs, &c_b50, &c_b50, &b[i__ + b_dim1], 
			ldb, (ftnlen)1);
#line 513 "sgelss.f"
	    }
#line 514 "sgelss.f"
/* L10: */
#line 514 "sgelss.f"
	}

/*        Multiply B by right singular vectors */
/*        (Workspace: need N, prefer N*NRHS) */

#line 519 "sgelss.f"
	if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 520 "sgelss.f"
	    sgemm_("T", "N", n, nrhs, n, &c_b84, &a[a_offset], lda, &b[
		    b_offset], ldb, &c_b50, &work[1], ldb, (ftnlen)1, (ftnlen)
		    1);
#line 522 "sgelss.f"
	    slacpy_("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1)
		    ;
#line 523 "sgelss.f"
	} else if (*nrhs > 1) {
#line 524 "sgelss.f"
	    chunk = *lwork / *n;
#line 525 "sgelss.f"
	    i__1 = *nrhs;
#line 525 "sgelss.f"
	    i__2 = chunk;
#line 525 "sgelss.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 526 "sgelss.f"
		i__3 = *nrhs - i__ + 1;
#line 526 "sgelss.f"
		bl = min(i__3,chunk);
#line 527 "sgelss.f"
		sgemm_("T", "N", n, &bl, n, &c_b84, &a[a_offset], lda, &b[i__ 
			* b_dim1 + 1], ldb, &c_b50, &work[1], n, (ftnlen)1, (
			ftnlen)1);
#line 529 "sgelss.f"
		slacpy_("G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (
			ftnlen)1);
#line 530 "sgelss.f"
/* L20: */
#line 530 "sgelss.f"
	    }
#line 531 "sgelss.f"
	} else {
#line 532 "sgelss.f"
	    sgemv_("T", n, n, &c_b84, &a[a_offset], lda, &b[b_offset], &c__1, 
		    &c_b50, &work[1], &c__1, (ftnlen)1);
#line 533 "sgelss.f"
	    scopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 534 "sgelss.f"
	}

#line 536 "sgelss.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 536 "sgelss.f"
	i__2 = *m, i__1 = (*m << 1) - 4, i__2 = max(i__2,i__1), i__2 = max(
		i__2,*nrhs), i__1 = *n - *m * 3;
#line 536 "sgelss.f"
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__2,i__1)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm */

#line 542 "sgelss.f"
	    ldwork = *m;
/* Computing MAX */
/* Computing MAX */
#line 543 "sgelss.f"
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 543 "sgelss.f"
	    i__2 = (*m << 2) + *m * *lda + max(i__3,i__4), i__1 = *m * *lda + 
		    *m + *m * *nrhs;
#line 543 "sgelss.f"
	    if (*lwork >= max(i__2,i__1)) {
#line 543 "sgelss.f"
		ldwork = *lda;
#line 543 "sgelss.f"
	    }
#line 545 "sgelss.f"
	    itau = 1;
#line 546 "sgelss.f"
	    iwork = *m + 1;

/*        Compute A=L*Q */
/*        (Workspace: need 2*M, prefer M+M*NB) */

#line 551 "sgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 551 "sgelss.f"
	    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
		     info);
#line 553 "sgelss.f"
	    il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */

#line 557 "sgelss.f"
	    slacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 558 "sgelss.f"
	    i__2 = *m - 1;
#line 558 "sgelss.f"
	    i__1 = *m - 1;
#line 558 "sgelss.f"
	    slaset_("U", &i__2, &i__1, &c_b50, &c_b50, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 560 "sgelss.f"
	    ie = il + ldwork * *m;
#line 561 "sgelss.f"
	    itauq = ie + *m;
#line 562 "sgelss.f"
	    itaup = itauq + *m;
#line 563 "sgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL) */
/*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

#line 568 "sgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 568 "sgelss.f"
	    sgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L */
/*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

#line 575 "sgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 575 "sgelss.f"
	    sormbr_("Q", "L", "T", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[iwork], &i__2, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) */
/*        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB) */

#line 582 "sgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 582 "sgelss.f"
	    sorgbr_("P", m, m, m, &work[il], &ldwork, &work[itaup], &work[
		    iwork], &i__2, info, (ftnlen)1);
#line 584 "sgelss.f"
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of L in WORK(IL) and */
/*           multiplying B by transpose of left singular vectors */
/*        (Workspace: need M*M+M+BDSPAC) */

#line 591 "sgelss.f"
	    sbdsqr_("U", m, m, &c__0, nrhs, &s[1], &work[ie], &work[il], &
		    ldwork, &a[a_offset], lda, &b[b_offset], ldb, &work[iwork]
		    , info, (ftnlen)1);
#line 593 "sgelss.f"
	    if (*info != 0) {
#line 593 "sgelss.f"
		goto L70;
#line 593 "sgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 598 "sgelss.f"
	    d__1 = *rcond * s[1];
#line 598 "sgelss.f"
	    thr = max(d__1,sfmin);
#line 599 "sgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 599 "sgelss.f"
		d__1 = eps * s[1];
#line 599 "sgelss.f"
		thr = max(d__1,sfmin);
#line 599 "sgelss.f"
	    }
#line 601 "sgelss.f"
	    *rank = 0;
#line 602 "sgelss.f"
	    i__2 = *m;
#line 602 "sgelss.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 603 "sgelss.f"
		if (s[i__] > thr) {
#line 604 "sgelss.f"
		    srscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 605 "sgelss.f"
		    ++(*rank);
#line 606 "sgelss.f"
		} else {
#line 607 "sgelss.f"
		    slaset_("F", &c__1, nrhs, &c_b50, &c_b50, &b[i__ + b_dim1]
			    , ldb, (ftnlen)1);
#line 608 "sgelss.f"
		}
#line 609 "sgelss.f"
/* L30: */
#line 609 "sgelss.f"
	    }
#line 610 "sgelss.f"
	    iwork = ie;

/*        Multiply B by right singular vectors of L in WORK(IL) */
/*        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS) */

#line 615 "sgelss.f"
	    if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
#line 616 "sgelss.f"
		sgemm_("T", "N", m, nrhs, m, &c_b84, &work[il], &ldwork, &b[
			b_offset], ldb, &c_b50, &work[iwork], ldb, (ftnlen)1, 
			(ftnlen)1);
#line 618 "sgelss.f"
		slacpy_("G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 619 "sgelss.f"
	    } else if (*nrhs > 1) {
#line 620 "sgelss.f"
		chunk = (*lwork - iwork + 1) / *m;
#line 621 "sgelss.f"
		i__2 = *nrhs;
#line 621 "sgelss.f"
		i__1 = chunk;
#line 621 "sgelss.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 622 "sgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 622 "sgelss.f"
		    bl = min(i__3,chunk);
#line 623 "sgelss.f"
		    sgemm_("T", "N", m, &bl, m, &c_b84, &work[il], &ldwork, &
			    b[i__ * b_dim1 + 1], ldb, &c_b50, &work[iwork], m,
			     (ftnlen)1, (ftnlen)1);
#line 625 "sgelss.f"
		    slacpy_("G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1]
			    , ldb, (ftnlen)1);
#line 627 "sgelss.f"
/* L40: */
#line 627 "sgelss.f"
		}
#line 628 "sgelss.f"
	    } else {
#line 629 "sgelss.f"
		sgemv_("T", m, m, &c_b84, &work[il], &ldwork, &b[b_dim1 + 1], 
			&c__1, &c_b50, &work[iwork], &c__1, (ftnlen)1);
#line 631 "sgelss.f"
		scopy_(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
#line 632 "sgelss.f"
	    }

/*        Zero out below first M rows of B */

#line 636 "sgelss.f"
	    i__1 = *n - *m;
#line 636 "sgelss.f"
	    slaset_("F", &i__1, nrhs, &c_b50, &c_b50, &b[*m + 1 + b_dim1], 
		    ldb, (ftnlen)1);
#line 637 "sgelss.f"
	    iwork = itau + *m;

/*        Multiply transpose(Q) by B */
/*        (Workspace: need M+NRHS, prefer M+NRHS*NB) */

#line 642 "sgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 642 "sgelss.f"
	    sormlq_("L", "T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 645 "sgelss.f"
	} else {

/*        Path 2 - remaining underdetermined cases */

#line 649 "sgelss.f"
	    ie = 1;
#line 650 "sgelss.f"
	    itauq = ie + *m;
#line 651 "sgelss.f"
	    itaup = itauq + *m;
#line 652 "sgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize A */
/*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 657 "sgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 657 "sgelss.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors */
/*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

#line 664 "sgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 664 "sgelss.f"
	    sormbr_("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors in A */
/*        (Workspace: need 4*M, prefer 3*M+M*NB) */

#line 670 "sgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 670 "sgelss.f"
	    sorgbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
		    iwork], &i__1, info, (ftnlen)1);
#line 672 "sgelss.f"
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of A in A and */
/*           multiplying B by transpose of left singular vectors */
/*        (Workspace: need BDSPAC) */

#line 679 "sgelss.f"
	    sbdsqr_("L", m, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], 
		    lda, dum, &c__1, &b[b_offset], ldb, &work[iwork], info, (
		    ftnlen)1);
#line 681 "sgelss.f"
	    if (*info != 0) {
#line 681 "sgelss.f"
		goto L70;
#line 681 "sgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 686 "sgelss.f"
	    d__1 = *rcond * s[1];
#line 686 "sgelss.f"
	    thr = max(d__1,sfmin);
#line 687 "sgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 687 "sgelss.f"
		d__1 = eps * s[1];
#line 687 "sgelss.f"
		thr = max(d__1,sfmin);
#line 687 "sgelss.f"
	    }
#line 689 "sgelss.f"
	    *rank = 0;
#line 690 "sgelss.f"
	    i__1 = *m;
#line 690 "sgelss.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 691 "sgelss.f"
		if (s[i__] > thr) {
#line 692 "sgelss.f"
		    srscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 693 "sgelss.f"
		    ++(*rank);
#line 694 "sgelss.f"
		} else {
#line 695 "sgelss.f"
		    slaset_("F", &c__1, nrhs, &c_b50, &c_b50, &b[i__ + b_dim1]
			    , ldb, (ftnlen)1);
#line 696 "sgelss.f"
		}
#line 697 "sgelss.f"
/* L50: */
#line 697 "sgelss.f"
	    }

/*        Multiply B by right singular vectors of A */
/*        (Workspace: need N, prefer N*NRHS) */

#line 702 "sgelss.f"
	    if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 703 "sgelss.f"
		sgemm_("T", "N", n, nrhs, m, &c_b84, &a[a_offset], lda, &b[
			b_offset], ldb, &c_b50, &work[1], ldb, (ftnlen)1, (
			ftnlen)1);
#line 705 "sgelss.f"
		slacpy_("F", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 706 "sgelss.f"
	    } else if (*nrhs > 1) {
#line 707 "sgelss.f"
		chunk = *lwork / *n;
#line 708 "sgelss.f"
		i__1 = *nrhs;
#line 708 "sgelss.f"
		i__2 = chunk;
#line 708 "sgelss.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 709 "sgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 709 "sgelss.f"
		    bl = min(i__3,chunk);
#line 710 "sgelss.f"
		    sgemm_("T", "N", n, &bl, m, &c_b84, &a[a_offset], lda, &b[
			    i__ * b_dim1 + 1], ldb, &c_b50, &work[1], n, (
			    ftnlen)1, (ftnlen)1);
#line 712 "sgelss.f"
		    slacpy_("F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], 
			    ldb, (ftnlen)1);
#line 713 "sgelss.f"
/* L60: */
#line 713 "sgelss.f"
		}
#line 714 "sgelss.f"
	    } else {
#line 715 "sgelss.f"
		sgemv_("T", m, n, &c_b84, &a[a_offset], lda, &b[b_offset], &
			c__1, &c_b50, &work[1], &c__1, (ftnlen)1);
#line 716 "sgelss.f"
		scopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 717 "sgelss.f"
	    }
#line 718 "sgelss.f"
	}
#line 718 "sgelss.f"
    }

/*     Undo scaling */

#line 722 "sgelss.f"
    if (iascl == 1) {
#line 723 "sgelss.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 724 "sgelss.f"
	slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 726 "sgelss.f"
    } else if (iascl == 2) {
#line 727 "sgelss.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 728 "sgelss.f"
	slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 730 "sgelss.f"
    }
#line 731 "sgelss.f"
    if (ibscl == 1) {
#line 732 "sgelss.f"
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 733 "sgelss.f"
    } else if (ibscl == 2) {
#line 734 "sgelss.f"
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 735 "sgelss.f"
    }

#line 737 "sgelss.f"
L70:
#line 738 "sgelss.f"
    work[1] = (doublereal) maxwrk;
#line 739 "sgelss.f"
    return 0;

/*     End of SGELSS */

} /* sgelss_ */

