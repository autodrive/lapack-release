#line 1 "cgelss.f"
/* cgelss.f -- translated by f2c (version 20100827).
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

#line 1 "cgelss.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b59 = 0.;

/* > \brief <b> CGELSS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGELSS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelss.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelss.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelss.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), S( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGELSS computes the minimum norm solution to a complex linear */
/* > least squares problem: */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the first min(m,n) rows of A are overwritten with */
/* >          its right singular vectors, stored rowwise. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, B is overwritten by the N-by-NRHS solution matrix X. */
/* >          If m >= n and RANK = n, the residual sum-of-squares for */
/* >          the solution in the i-th column is given by the sum of */
/* >          squares of the modulus of elements n+1:m in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M,N). */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= 1, and also: */
/* >          LWORK >=  2*min(M,N) + max(M,N,NRHS) */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (5*min(M,N)) */
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

/* > \date June 2016 */

/* > \ingroup complexGEsolve */

/*  ===================================================================== */
/* Subroutine */ int cgelss_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, bl, ie, il, mm;
    static doublecomplex dum[1];
    static doublereal eps, thr, anrm, bnrm;
    static integer itau, lwork_cgebrd__, lwork_cgelqf__, lwork_cgeqrf__, 
	    lwork_cungbr__, lwork_cunmbr__, lwork_cunmlq__, lwork_cunmqr__;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer iascl, ibscl;
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer chunk;
    static doublereal sfmin;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer minmn, maxmn, itaup, itauq, mnthr, iwork;
    extern /* Subroutine */ int cgebrd_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), slabad_(
	    doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), cgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen), cbdsqr_(char *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int cungbr_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen), slascl_(char *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     integer *, integer *, ftnlen), cunmbr_(char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), csrscl_(integer *, 
	    doublereal *, doublecomplex *, integer *), slaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), cunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static integer irwork;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 235 "cgelss.f"
    /* Parameter adjustments */
#line 235 "cgelss.f"
    a_dim1 = *lda;
#line 235 "cgelss.f"
    a_offset = 1 + a_dim1;
#line 235 "cgelss.f"
    a -= a_offset;
#line 235 "cgelss.f"
    b_dim1 = *ldb;
#line 235 "cgelss.f"
    b_offset = 1 + b_dim1;
#line 235 "cgelss.f"
    b -= b_offset;
#line 235 "cgelss.f"
    --s;
#line 235 "cgelss.f"
    --work;
#line 235 "cgelss.f"
    --rwork;
#line 235 "cgelss.f"

#line 235 "cgelss.f"
    /* Function Body */
#line 235 "cgelss.f"
    *info = 0;
#line 236 "cgelss.f"
    minmn = min(*m,*n);
#line 237 "cgelss.f"
    maxmn = max(*m,*n);
#line 238 "cgelss.f"
    lquery = *lwork == -1;
#line 239 "cgelss.f"
    if (*m < 0) {
#line 240 "cgelss.f"
	*info = -1;
#line 241 "cgelss.f"
    } else if (*n < 0) {
#line 242 "cgelss.f"
	*info = -2;
#line 243 "cgelss.f"
    } else if (*nrhs < 0) {
#line 244 "cgelss.f"
	*info = -3;
#line 245 "cgelss.f"
    } else if (*lda < max(1,*m)) {
#line 246 "cgelss.f"
	*info = -5;
#line 247 "cgelss.f"
    } else if (*ldb < max(1,maxmn)) {
#line 248 "cgelss.f"
	*info = -7;
#line 249 "cgelss.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace refers */
/*       to real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 259 "cgelss.f"
    if (*info == 0) {
#line 260 "cgelss.f"
	minwrk = 1;
#line 261 "cgelss.f"
	maxwrk = 1;
#line 262 "cgelss.f"
	if (minmn > 0) {
#line 263 "cgelss.f"
	    mm = *m;
#line 264 "cgelss.f"
	    mnthr = ilaenv_(&c__6, "CGELSS", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
#line 265 "cgelss.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns */

/*              Compute space needed for CGEQRF */
#line 271 "cgelss.f"
		cgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 272 "cgelss.f"
		lwork_cgeqrf__ = (integer) dum[0].r;
/*              Compute space needed for CUNMQR */
#line 274 "cgelss.f"
		cunmqr_("L", "C", m, nrhs, n, &a[a_offset], lda, dum, &b[
			b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (ftnlen)
			1);
#line 276 "cgelss.f"
		lwork_cunmqr__ = (integer) dum[0].r;
#line 277 "cgelss.f"
		mm = *n;
/* Computing MAX */
#line 278 "cgelss.f"
		i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "CGEQRF", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 278 "cgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 280 "cgelss.f"
		i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, "CUNMQR", 
			"LC", m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 280 "cgelss.f"
		maxwrk = max(i__1,i__2);
#line 282 "cgelss.f"
	    }
#line 283 "cgelss.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined */

/*              Compute space needed for CGEBRD */
#line 288 "cgelss.f"
		cgebrd_(&mm, n, &a[a_offset], lda, &s[1], &s[1], dum, dum, 
			dum, &c_n1, info);
#line 290 "cgelss.f"
		lwork_cgebrd__ = (integer) dum[0].r;
/*              Compute space needed for CUNMBR */
#line 292 "cgelss.f"
		cunmbr_("Q", "L", "C", &mm, nrhs, n, &a[a_offset], lda, dum, &
			b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 294 "cgelss.f"
		lwork_cunmbr__ = (integer) dum[0].r;
/*              Compute space needed for CUNGBR */
#line 296 "cgelss.f"
		cungbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			info, (ftnlen)1);
#line 298 "cgelss.f"
		lwork_cungbr__ = (integer) dum[0].r;
/*              Compute total workspace needed */
/* Computing MAX */
#line 300 "cgelss.f"
		i__1 = maxwrk, i__2 = (*n << 1) + lwork_cgebrd__;
#line 300 "cgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 301 "cgelss.f"
		i__1 = maxwrk, i__2 = (*n << 1) + lwork_cunmbr__;
#line 301 "cgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 302 "cgelss.f"
		i__1 = maxwrk, i__2 = (*n << 1) + lwork_cungbr__;
#line 302 "cgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 303 "cgelss.f"
		i__1 = maxwrk, i__2 = *n * *nrhs;
#line 303 "cgelss.f"
		maxwrk = max(i__1,i__2);
#line 304 "cgelss.f"
		minwrk = (*n << 1) + max(*nrhs,*m);
#line 305 "cgelss.f"
	    }
#line 306 "cgelss.f"
	    if (*n > *m) {
#line 307 "cgelss.f"
		minwrk = (*m << 1) + max(*nrhs,*n);
#line 308 "cgelss.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                 than rows */

/*                 Compute space needed for CGELQF */
#line 314 "cgelss.f"
		    cgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 316 "cgelss.f"
		    lwork_cgelqf__ = (integer) dum[0].r;
/*                 Compute space needed for CGEBRD */
#line 318 "cgelss.f"
		    cgebrd_(m, m, &a[a_offset], lda, &s[1], &s[1], dum, dum, 
			    dum, &c_n1, info);
#line 320 "cgelss.f"
		    lwork_cgebrd__ = (integer) dum[0].r;
/*                 Compute space needed for CUNMBR */
#line 322 "cgelss.f"
		    cunmbr_("Q", "L", "C", m, nrhs, n, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 324 "cgelss.f"
		    lwork_cunmbr__ = (integer) dum[0].r;
/*                 Compute space needed for CUNGBR */
#line 326 "cgelss.f"
		    cungbr_("P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 328 "cgelss.f"
		    lwork_cungbr__ = (integer) dum[0].r;
/*                 Compute space needed for CUNMLQ */
#line 330 "cgelss.f"
		    cunmlq_("L", "C", n, nrhs, m, &a[a_offset], lda, dum, &b[
			    b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1);
#line 332 "cgelss.f"
		    lwork_cunmlq__ = (integer) dum[0].r;
/*                 Compute total workspace needed */
#line 334 "cgelss.f"
		    maxwrk = *m + lwork_cgelqf__;
/* Computing MAX */
#line 335 "cgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * *m + lwork_cgebrd__;
#line 335 "cgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 336 "cgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * *m + lwork_cunmbr__;
#line 336 "cgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 337 "cgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * *m + lwork_cungbr__;
#line 337 "cgelss.f"
		    maxwrk = max(i__1,i__2);
#line 338 "cgelss.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 339 "cgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 339 "cgelss.f"
			maxwrk = max(i__1,i__2);
#line 340 "cgelss.f"
		    } else {
/* Computing MAX */
#line 341 "cgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 341 "cgelss.f"
			maxwrk = max(i__1,i__2);
#line 342 "cgelss.f"
		    }
/* Computing MAX */
#line 343 "cgelss.f"
		    i__1 = maxwrk, i__2 = *m + lwork_cunmlq__;
#line 343 "cgelss.f"
		    maxwrk = max(i__1,i__2);
#line 344 "cgelss.f"
		} else {

/*                 Path 2 - underdetermined */

/*                 Compute space needed for CGEBRD */
#line 349 "cgelss.f"
		    cgebrd_(m, n, &a[a_offset], lda, &s[1], &s[1], dum, dum, 
			    dum, &c_n1, info);
#line 351 "cgelss.f"
		    lwork_cgebrd__ = (integer) dum[0].r;
/*                 Compute space needed for CUNMBR */
#line 353 "cgelss.f"
		    cunmbr_("Q", "L", "C", m, nrhs, m, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 355 "cgelss.f"
		    lwork_cunmbr__ = (integer) dum[0].r;
/*                 Compute space needed for CUNGBR */
#line 357 "cgelss.f"
		    cungbr_("P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 359 "cgelss.f"
		    lwork_cungbr__ = (integer) dum[0].r;
#line 360 "cgelss.f"
		    maxwrk = (*m << 1) + lwork_cgebrd__;
/* Computing MAX */
#line 361 "cgelss.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cunmbr__;
#line 361 "cgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 362 "cgelss.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_cungbr__;
#line 362 "cgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 363 "cgelss.f"
		    i__1 = maxwrk, i__2 = *n * *nrhs;
#line 363 "cgelss.f"
		    maxwrk = max(i__1,i__2);
#line 364 "cgelss.f"
		}
#line 365 "cgelss.f"
	    }
#line 366 "cgelss.f"
	    maxwrk = max(minwrk,maxwrk);
#line 367 "cgelss.f"
	}
#line 368 "cgelss.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 370 "cgelss.f"
	if (*lwork < minwrk && ! lquery) {
#line 370 "cgelss.f"
	    *info = -12;
#line 370 "cgelss.f"
	}
#line 372 "cgelss.f"
    }

#line 374 "cgelss.f"
    if (*info != 0) {
#line 375 "cgelss.f"
	i__1 = -(*info);
#line 375 "cgelss.f"
	xerbla_("CGELSS", &i__1, (ftnlen)6);
#line 376 "cgelss.f"
	return 0;
#line 377 "cgelss.f"
    } else if (lquery) {
#line 378 "cgelss.f"
	return 0;
#line 379 "cgelss.f"
    }

/*     Quick return if possible */

#line 383 "cgelss.f"
    if (*m == 0 || *n == 0) {
#line 384 "cgelss.f"
	*rank = 0;
#line 385 "cgelss.f"
	return 0;
#line 386 "cgelss.f"
    }

/*     Get machine parameters */

#line 390 "cgelss.f"
    eps = slamch_("P", (ftnlen)1);
#line 391 "cgelss.f"
    sfmin = slamch_("S", (ftnlen)1);
#line 392 "cgelss.f"
    smlnum = sfmin / eps;
#line 393 "cgelss.f"
    bignum = 1. / smlnum;
#line 394 "cgelss.f"
    slabad_(&smlnum, &bignum);

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 398 "cgelss.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 399 "cgelss.f"
    iascl = 0;
#line 400 "cgelss.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 404 "cgelss.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 405 "cgelss.f"
	iascl = 1;
#line 406 "cgelss.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 410 "cgelss.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 411 "cgelss.f"
	iascl = 2;
#line 412 "cgelss.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 416 "cgelss.f"
	i__1 = max(*m,*n);
#line 416 "cgelss.f"
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 417 "cgelss.f"
	slaset_("F", &minmn, &c__1, &c_b59, &c_b59, &s[1], &minmn, (ftnlen)1);
#line 418 "cgelss.f"
	*rank = 0;
#line 419 "cgelss.f"
	goto L70;
#line 420 "cgelss.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 424 "cgelss.f"
    bnrm = clange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 425 "cgelss.f"
    ibscl = 0;
#line 426 "cgelss.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 430 "cgelss.f"
	clascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 431 "cgelss.f"
	ibscl = 1;
#line 432 "cgelss.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 436 "cgelss.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 437 "cgelss.f"
	ibscl = 2;
#line 438 "cgelss.f"
    }

/*     Overdetermined case */

#line 442 "cgelss.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

#line 446 "cgelss.f"
	mm = *m;
#line 447 "cgelss.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

#line 451 "cgelss.f"
	    mm = *n;
#line 452 "cgelss.f"
	    itau = 1;
#line 453 "cgelss.f"
	    iwork = itau + *n;

/*           Compute A=Q*R */
/*           (CWorkspace: need 2*N, prefer N+N*NB) */
/*           (RWorkspace: none) */

#line 459 "cgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 459 "cgelss.f"
	    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1,
		     info);

/*           Multiply B by transpose(Q) */
/*           (CWorkspace: need N+NRHS, prefer N+NRHS*NB) */
/*           (RWorkspace: none) */

#line 466 "cgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 466 "cgelss.f"
	    cunmqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R */

#line 471 "cgelss.f"
	    if (*n > 1) {
#line 471 "cgelss.f"
		i__1 = *n - 1;
#line 471 "cgelss.f"
		i__2 = *n - 1;
#line 471 "cgelss.f"
		claset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 471 "cgelss.f"
	    }
#line 474 "cgelss.f"
	}

#line 476 "cgelss.f"
	ie = 1;
#line 477 "cgelss.f"
	itauq = 1;
#line 478 "cgelss.f"
	itaup = itauq + *n;
#line 479 "cgelss.f"
	iwork = itaup + *n;

/*        Bidiagonalize R in A */
/*        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB) */
/*        (RWorkspace: need N) */

#line 485 "cgelss.f"
	i__1 = *lwork - iwork + 1;
#line 485 "cgelss.f"
	cgebrd_(&mm, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], &
		work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R */
/*        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB) */
/*        (RWorkspace: none) */

#line 493 "cgelss.f"
	i__1 = *lwork - iwork + 1;
#line 493 "cgelss.f"
	cunmbr_("Q", "L", "C", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in A */
/*        (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 500 "cgelss.f"
	i__1 = *lwork - iwork + 1;
#line 500 "cgelss.f"
	cungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &
		i__1, info, (ftnlen)1);
#line 502 "cgelss.f"
	irwork = ie + *n;

/*        Perform bidiagonal QR iteration */
/*          multiply B by transpose of left singular vectors */
/*          compute right singular vectors in A */
/*        (CWorkspace: none) */
/*        (RWorkspace: need BDSPAC) */

#line 510 "cgelss.f"
	cbdsqr_("U", n, n, &c__0, nrhs, &s[1], &rwork[ie], &a[a_offset], lda, 
		dum, &c__1, &b[b_offset], ldb, &rwork[irwork], info, (ftnlen)
		1);
#line 512 "cgelss.f"
	if (*info != 0) {
#line 512 "cgelss.f"
	    goto L70;
#line 512 "cgelss.f"
	}

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 517 "cgelss.f"
	d__1 = *rcond * s[1];
#line 517 "cgelss.f"
	thr = max(d__1,sfmin);
#line 518 "cgelss.f"
	if (*rcond < 0.) {
/* Computing MAX */
#line 518 "cgelss.f"
	    d__1 = eps * s[1];
#line 518 "cgelss.f"
	    thr = max(d__1,sfmin);
#line 518 "cgelss.f"
	}
#line 520 "cgelss.f"
	*rank = 0;
#line 521 "cgelss.f"
	i__1 = *n;
#line 521 "cgelss.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 522 "cgelss.f"
	    if (s[i__] > thr) {
#line 523 "cgelss.f"
		csrscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 524 "cgelss.f"
		++(*rank);
#line 525 "cgelss.f"
	    } else {
#line 526 "cgelss.f"
		claset_("F", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 527 "cgelss.f"
	    }
#line 528 "cgelss.f"
/* L10: */
#line 528 "cgelss.f"
	}

/*        Multiply B by right singular vectors */
/*        (CWorkspace: need N, prefer N*NRHS) */
/*        (RWorkspace: none) */

#line 534 "cgelss.f"
	if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 535 "cgelss.f"
	    cgemm_("C", "N", n, nrhs, n, &c_b2, &a[a_offset], lda, &b[
		    b_offset], ldb, &c_b1, &work[1], ldb, (ftnlen)1, (ftnlen)
		    1);
#line 537 "cgelss.f"
	    clacpy_("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1)
		    ;
#line 538 "cgelss.f"
	} else if (*nrhs > 1) {
#line 539 "cgelss.f"
	    chunk = *lwork / *n;
#line 540 "cgelss.f"
	    i__1 = *nrhs;
#line 540 "cgelss.f"
	    i__2 = chunk;
#line 540 "cgelss.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 541 "cgelss.f"
		i__3 = *nrhs - i__ + 1;
#line 541 "cgelss.f"
		bl = min(i__3,chunk);
#line 542 "cgelss.f"
		cgemm_("C", "N", n, &bl, n, &c_b2, &a[a_offset], lda, &b[i__ *
			 b_dim1 + 1], ldb, &c_b1, &work[1], n, (ftnlen)1, (
			ftnlen)1);
#line 544 "cgelss.f"
		clacpy_("G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (
			ftnlen)1);
#line 545 "cgelss.f"
/* L20: */
#line 545 "cgelss.f"
	    }
#line 546 "cgelss.f"
	} else {
#line 547 "cgelss.f"
	    cgemv_("C", n, n, &c_b2, &a[a_offset], lda, &b[b_offset], &c__1, &
		    c_b1, &work[1], &c__1, (ftnlen)1);
#line 548 "cgelss.f"
	    ccopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 549 "cgelss.f"
	}

#line 551 "cgelss.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 551 "cgelss.f"
	i__2 = max(*m,*nrhs), i__1 = *n - (*m << 1);
#line 551 "cgelss.f"
	if (*n >= mnthr && *lwork >= *m * 3 + *m * *m + max(i__2,i__1)) {

/*        Underdetermined case, M much less than N */

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm */

#line 559 "cgelss.f"
	    ldwork = *m;
/* Computing MAX */
#line 560 "cgelss.f"
	    i__2 = max(*m,*nrhs), i__1 = *n - (*m << 1);
#line 560 "cgelss.f"
	    if (*lwork >= *m * 3 + *m * *lda + max(i__2,i__1)) {
#line 560 "cgelss.f"
		ldwork = *lda;
#line 560 "cgelss.f"
	    }
#line 562 "cgelss.f"
	    itau = 1;
#line 563 "cgelss.f"
	    iwork = *m + 1;

/*        Compute A=L*Q */
/*        (CWorkspace: need 2*M, prefer M+M*NB) */
/*        (RWorkspace: none) */

#line 569 "cgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 569 "cgelss.f"
	    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
		     info);
#line 571 "cgelss.f"
	    il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */

#line 575 "cgelss.f"
	    clacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 576 "cgelss.f"
	    i__2 = *m - 1;
#line 576 "cgelss.f"
	    i__1 = *m - 1;
#line 576 "cgelss.f"
	    claset_("U", &i__2, &i__1, &c_b1, &c_b1, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 578 "cgelss.f"
	    ie = 1;
#line 579 "cgelss.f"
	    itauq = il + ldwork * *m;
#line 580 "cgelss.f"
	    itaup = itauq + *m;
#line 581 "cgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL) */
/*        (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */
/*        (RWorkspace: need M) */

#line 587 "cgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 587 "cgelss.f"
	    cgebrd_(m, m, &work[il], &ldwork, &s[1], &rwork[ie], &work[itauq],
		     &work[itaup], &work[iwork], &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L */
/*        (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB) */
/*        (RWorkspace: none) */

#line 595 "cgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 595 "cgelss.f"
	    cunmbr_("Q", "L", "C", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[iwork], &i__2, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) */
/*        (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */
/*        (RWorkspace: none) */

#line 603 "cgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 603 "cgelss.f"
	    cungbr_("P", m, m, m, &work[il], &ldwork, &work[itaup], &work[
		    iwork], &i__2, info, (ftnlen)1);
#line 605 "cgelss.f"
	    irwork = ie + *m;

/*        Perform bidiagonal QR iteration, computing right singular */
/*        vectors of L in WORK(IL) and multiplying B by transpose of */
/*        left singular vectors */
/*        (CWorkspace: need M*M) */
/*        (RWorkspace: need BDSPAC) */

#line 613 "cgelss.f"
	    cbdsqr_("U", m, m, &c__0, nrhs, &s[1], &rwork[ie], &work[il], &
		    ldwork, &a[a_offset], lda, &b[b_offset], ldb, &rwork[
		    irwork], info, (ftnlen)1);
#line 615 "cgelss.f"
	    if (*info != 0) {
#line 615 "cgelss.f"
		goto L70;
#line 615 "cgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 620 "cgelss.f"
	    d__1 = *rcond * s[1];
#line 620 "cgelss.f"
	    thr = max(d__1,sfmin);
#line 621 "cgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 621 "cgelss.f"
		d__1 = eps * s[1];
#line 621 "cgelss.f"
		thr = max(d__1,sfmin);
#line 621 "cgelss.f"
	    }
#line 623 "cgelss.f"
	    *rank = 0;
#line 624 "cgelss.f"
	    i__2 = *m;
#line 624 "cgelss.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 625 "cgelss.f"
		if (s[i__] > thr) {
#line 626 "cgelss.f"
		    csrscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 627 "cgelss.f"
		    ++(*rank);
#line 628 "cgelss.f"
		} else {
#line 629 "cgelss.f"
		    claset_("F", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], 
			    ldb, (ftnlen)1);
#line 630 "cgelss.f"
		}
#line 631 "cgelss.f"
/* L30: */
#line 631 "cgelss.f"
	    }
#line 632 "cgelss.f"
	    iwork = il + *m * ldwork;

/*        Multiply B by right singular vectors of L in WORK(IL) */
/*        (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS) */
/*        (RWorkspace: none) */

#line 638 "cgelss.f"
	    if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
#line 639 "cgelss.f"
		cgemm_("C", "N", m, nrhs, m, &c_b2, &work[il], &ldwork, &b[
			b_offset], ldb, &c_b1, &work[iwork], ldb, (ftnlen)1, (
			ftnlen)1);
#line 641 "cgelss.f"
		clacpy_("G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 642 "cgelss.f"
	    } else if (*nrhs > 1) {
#line 643 "cgelss.f"
		chunk = (*lwork - iwork + 1) / *m;
#line 644 "cgelss.f"
		i__2 = *nrhs;
#line 644 "cgelss.f"
		i__1 = chunk;
#line 644 "cgelss.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 645 "cgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 645 "cgelss.f"
		    bl = min(i__3,chunk);
#line 646 "cgelss.f"
		    cgemm_("C", "N", m, &bl, m, &c_b2, &work[il], &ldwork, &b[
			    i__ * b_dim1 + 1], ldb, &c_b1, &work[iwork], m, (
			    ftnlen)1, (ftnlen)1);
#line 648 "cgelss.f"
		    clacpy_("G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1]
			    , ldb, (ftnlen)1);
#line 650 "cgelss.f"
/* L40: */
#line 650 "cgelss.f"
		}
#line 651 "cgelss.f"
	    } else {
#line 652 "cgelss.f"
		cgemv_("C", m, m, &c_b2, &work[il], &ldwork, &b[b_dim1 + 1], &
			c__1, &c_b1, &work[iwork], &c__1, (ftnlen)1);
#line 654 "cgelss.f"
		ccopy_(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
#line 655 "cgelss.f"
	    }

/*        Zero out below first M rows of B */

#line 659 "cgelss.f"
	    i__1 = *n - *m;
#line 659 "cgelss.f"
	    claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[*m + 1 + b_dim1], ldb, 
		    (ftnlen)1);
#line 660 "cgelss.f"
	    iwork = itau + *m;

/*        Multiply transpose(Q) by B */
/*        (CWorkspace: need M+NRHS, prefer M+NHRS*NB) */
/*        (RWorkspace: none) */

#line 666 "cgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 666 "cgelss.f"
	    cunmlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 669 "cgelss.f"
	} else {

/*        Path 2 - remaining underdetermined cases */

#line 673 "cgelss.f"
	    ie = 1;
#line 674 "cgelss.f"
	    itauq = 1;
#line 675 "cgelss.f"
	    itaup = itauq + *m;
#line 676 "cgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize A */
/*        (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB) */
/*        (RWorkspace: need N) */

#line 682 "cgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 682 "cgelss.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors */
/*        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB) */
/*        (RWorkspace: none) */

#line 690 "cgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 690 "cgelss.f"
	    cunmbr_("Q", "L", "C", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors in A */
/*        (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*        (RWorkspace: none) */

#line 697 "cgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 697 "cgelss.f"
	    cungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
		    iwork], &i__1, info, (ftnlen)1);
#line 699 "cgelss.f"
	    irwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of A in A and */
/*           multiplying B by transpose of left singular vectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need BDSPAC) */

#line 707 "cgelss.f"
	    cbdsqr_("L", m, n, &c__0, nrhs, &s[1], &rwork[ie], &a[a_offset], 
		    lda, dum, &c__1, &b[b_offset], ldb, &rwork[irwork], info, 
		    (ftnlen)1);
#line 709 "cgelss.f"
	    if (*info != 0) {
#line 709 "cgelss.f"
		goto L70;
#line 709 "cgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 714 "cgelss.f"
	    d__1 = *rcond * s[1];
#line 714 "cgelss.f"
	    thr = max(d__1,sfmin);
#line 715 "cgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 715 "cgelss.f"
		d__1 = eps * s[1];
#line 715 "cgelss.f"
		thr = max(d__1,sfmin);
#line 715 "cgelss.f"
	    }
#line 717 "cgelss.f"
	    *rank = 0;
#line 718 "cgelss.f"
	    i__1 = *m;
#line 718 "cgelss.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 719 "cgelss.f"
		if (s[i__] > thr) {
#line 720 "cgelss.f"
		    csrscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 721 "cgelss.f"
		    ++(*rank);
#line 722 "cgelss.f"
		} else {
#line 723 "cgelss.f"
		    claset_("F", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], 
			    ldb, (ftnlen)1);
#line 724 "cgelss.f"
		}
#line 725 "cgelss.f"
/* L50: */
#line 725 "cgelss.f"
	    }

/*        Multiply B by right singular vectors of A */
/*        (CWorkspace: need N, prefer N*NRHS) */
/*        (RWorkspace: none) */

#line 731 "cgelss.f"
	    if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 732 "cgelss.f"
		cgemm_("C", "N", n, nrhs, m, &c_b2, &a[a_offset], lda, &b[
			b_offset], ldb, &c_b1, &work[1], ldb, (ftnlen)1, (
			ftnlen)1);
#line 734 "cgelss.f"
		clacpy_("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 735 "cgelss.f"
	    } else if (*nrhs > 1) {
#line 736 "cgelss.f"
		chunk = *lwork / *n;
#line 737 "cgelss.f"
		i__1 = *nrhs;
#line 737 "cgelss.f"
		i__2 = chunk;
#line 737 "cgelss.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 738 "cgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 738 "cgelss.f"
		    bl = min(i__3,chunk);
#line 739 "cgelss.f"
		    cgemm_("C", "N", n, &bl, m, &c_b2, &a[a_offset], lda, &b[
			    i__ * b_dim1 + 1], ldb, &c_b1, &work[1], n, (
			    ftnlen)1, (ftnlen)1);
#line 741 "cgelss.f"
		    clacpy_("F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], 
			    ldb, (ftnlen)1);
#line 742 "cgelss.f"
/* L60: */
#line 742 "cgelss.f"
		}
#line 743 "cgelss.f"
	    } else {
#line 744 "cgelss.f"
		cgemv_("C", m, n, &c_b2, &a[a_offset], lda, &b[b_offset], &
			c__1, &c_b1, &work[1], &c__1, (ftnlen)1);
#line 745 "cgelss.f"
		ccopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 746 "cgelss.f"
	    }
#line 747 "cgelss.f"
	}
#line 747 "cgelss.f"
    }

/*     Undo scaling */

#line 751 "cgelss.f"
    if (iascl == 1) {
#line 752 "cgelss.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 753 "cgelss.f"
	slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 755 "cgelss.f"
    } else if (iascl == 2) {
#line 756 "cgelss.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 757 "cgelss.f"
	slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 759 "cgelss.f"
    }
#line 760 "cgelss.f"
    if (ibscl == 1) {
#line 761 "cgelss.f"
	clascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 762 "cgelss.f"
    } else if (ibscl == 2) {
#line 763 "cgelss.f"
	clascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 764 "cgelss.f"
    }
#line 765 "cgelss.f"
L70:
#line 766 "cgelss.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 767 "cgelss.f"
    return 0;

/*     End of CGELSS */

} /* cgelss_ */

