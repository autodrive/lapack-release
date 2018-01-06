#line 1 "zgelsy.f"
/* zgelsy.f -- translated by f2c (version 20100827).
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

#line 1 "zgelsy.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static integer c__2 = 2;

/* > \brief <b> ZGELSY solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGELSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelsy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelsy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelsy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/*                          WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGELSY computes the minimum-norm solution to a complex linear least */
/* > squares problem: */
/* >     minimize || A * X - B || */
/* > using a complete orthogonal factorization of A.  A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > */
/* > The routine first computes a QR factorization with column pivoting: */
/* >     A * P = Q * [ R11 R12 ] */
/* >                 [  0  R22 ] */
/* > with R11 defined as the largest leading submatrix whose estimated */
/* > condition number is less than 1/RCOND.  The order of R11, RANK, */
/* > is the effective rank of A. */
/* > */
/* > Then, R22 is considered to be negligible, and R12 is annihilated */
/* > by unitary transformations from the right, arriving at the */
/* > complete orthogonal factorization: */
/* >    A * P = Q * [ T11 0 ] * Z */
/* >                [  0  0 ] */
/* > The minimum-norm solution is then */
/* >    X = P * Z**H [ inv(T11)*Q1**H*B ] */
/* >                 [        0         ] */
/* > where Q1 consists of the first RANK columns of Q. */
/* > */
/* > This routine is basically identical to the original xGELSX except */
/* > three differences: */
/* >   o The permutation of matrix B (the right hand side) is faster and */
/* >     more simple. */
/* >   o The call to the subroutine xGEQPF has been substituted by the */
/* >     the call to the subroutine xGEQP3. This subroutine is a Blas-3 */
/* >     version of the QR factorization with column pivoting. */
/* >   o Matrix B (the right hand side) is updated with Blas-3. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of */
/* >          columns of matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A has been overwritten by details of its */
/* >          complete orthogonal factorization. */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,M,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* >          JPVT is INTEGER array, dimension (N) */
/* >          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/* >          to the front of AP, otherwise column i is a free column. */
/* >          On exit, if JPVT(i) = k, then the i-th column of A*P */
/* >          was the k-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          RCOND is used to determine the effective rank of A, which */
/* >          is defined as the order of the largest leading triangular */
/* >          submatrix R11 in the QR factorization with pivoting of A, */
/* >          whose estimated condition number < 1/RCOND. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* >          RANK is INTEGER */
/* >          The effective rank of A, i.e., the order of the submatrix */
/* >          R11.  This is the same as the order of the submatrix T11 */
/* >          in the complete orthogonal factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          The unblocked strategy requires that: */
/* >            LWORK >= MN + MAX( 2*MN, N+1, MN+NRHS ) */
/* >          where MN = min(M,N). */
/* >          The block algorithm requires that: */
/* >            LWORK >= MN + MAX( 2*MN, NB*(N+1), MN+MN*NB, MN+NB*NRHS ) */
/* >          where NB is an upper bound on the blocksize returned */
/* >          by ILAENV for the routines ZGEQP3, ZTZRZF, CTZRQF, ZUNMQR, */
/* >          and ZUNMRZ. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GEsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n */
/* >    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgelsy_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex c1, c2, s1, s2;
    static integer nb, mn, nb1, nb2, nb3, nb4;
    static doublereal anrm, bnrm, smin, smax;
    static integer iascl, ibscl, ismin, ismax;
    static doublereal wsize;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    zlaic1_(integer *, integer *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *), dlabad_(doublereal *, doublereal *), zgeqp3_(
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), zlaset_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen);
    static doublereal sminpr, smaxpr, smlnum;
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), ztzrzf_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 261 "zgelsy.f"
    /* Parameter adjustments */
#line 261 "zgelsy.f"
    a_dim1 = *lda;
#line 261 "zgelsy.f"
    a_offset = 1 + a_dim1;
#line 261 "zgelsy.f"
    a -= a_offset;
#line 261 "zgelsy.f"
    b_dim1 = *ldb;
#line 261 "zgelsy.f"
    b_offset = 1 + b_dim1;
#line 261 "zgelsy.f"
    b -= b_offset;
#line 261 "zgelsy.f"
    --jpvt;
#line 261 "zgelsy.f"
    --work;
#line 261 "zgelsy.f"
    --rwork;
#line 261 "zgelsy.f"

#line 261 "zgelsy.f"
    /* Function Body */
#line 261 "zgelsy.f"
    mn = min(*m,*n);
#line 262 "zgelsy.f"
    ismin = mn + 1;
#line 263 "zgelsy.f"
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

#line 267 "zgelsy.f"
    *info = 0;
#line 268 "zgelsy.f"
    nb1 = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 269 "zgelsy.f"
    nb2 = ilaenv_(&c__1, "ZGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 270 "zgelsy.f"
    nb3 = ilaenv_(&c__1, "ZUNMQR", " ", m, n, nrhs, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 271 "zgelsy.f"
    nb4 = ilaenv_(&c__1, "ZUNMRQ", " ", m, n, nrhs, &c_n1, (ftnlen)6, (ftnlen)
	    1);
/* Computing MAX */
#line 272 "zgelsy.f"
    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 272 "zgelsy.f"
    nb = max(i__1,nb4);
/* Computing MAX */
#line 273 "zgelsy.f"
    i__1 = 1, i__2 = mn + (*n << 1) + nb * (*n + 1), i__1 = max(i__1,i__2), 
	    i__2 = (mn << 1) + nb * *nrhs;
#line 273 "zgelsy.f"
    lwkopt = max(i__1,i__2);
#line 274 "zgelsy.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 274 "zgelsy.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 275 "zgelsy.f"
    lquery = *lwork == -1;
#line 276 "zgelsy.f"
    if (*m < 0) {
#line 277 "zgelsy.f"
	*info = -1;
#line 278 "zgelsy.f"
    } else if (*n < 0) {
#line 279 "zgelsy.f"
	*info = -2;
#line 280 "zgelsy.f"
    } else if (*nrhs < 0) {
#line 281 "zgelsy.f"
	*info = -3;
#line 282 "zgelsy.f"
    } else if (*lda < max(1,*m)) {
#line 283 "zgelsy.f"
	*info = -5;
#line 284 "zgelsy.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 284 "zgelsy.f"
	i__1 = max(1,*m);
#line 284 "zgelsy.f"
	if (*ldb < max(i__1,*n)) {
#line 285 "zgelsy.f"
	    *info = -7;
#line 286 "zgelsy.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 286 "zgelsy.f"
	    i__1 = mn << 1, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = mn + 
		    *nrhs;
#line 286 "zgelsy.f"
	    if (*lwork < mn + max(i__1,i__2) && ! lquery) {
#line 288 "zgelsy.f"
		*info = -12;
#line 289 "zgelsy.f"
	    }
#line 289 "zgelsy.f"
	}
#line 289 "zgelsy.f"
    }

#line 291 "zgelsy.f"
    if (*info != 0) {
#line 292 "zgelsy.f"
	i__1 = -(*info);
#line 292 "zgelsy.f"
	xerbla_("ZGELSY", &i__1, (ftnlen)6);
#line 293 "zgelsy.f"
	return 0;
#line 294 "zgelsy.f"
    } else if (lquery) {
#line 295 "zgelsy.f"
	return 0;
#line 296 "zgelsy.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 300 "zgelsy.f"
    i__1 = min(*m,*n);
#line 300 "zgelsy.f"
    if (min(i__1,*nrhs) == 0) {
#line 301 "zgelsy.f"
	*rank = 0;
#line 302 "zgelsy.f"
	return 0;
#line 303 "zgelsy.f"
    }

/*     Get machine parameters */

#line 307 "zgelsy.f"
    smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 308 "zgelsy.f"
    bignum = 1. / smlnum;
#line 309 "zgelsy.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max entries outside range [SMLNUM,BIGNUM] */

#line 313 "zgelsy.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 314 "zgelsy.f"
    iascl = 0;
#line 315 "zgelsy.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 319 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 320 "zgelsy.f"
	iascl = 1;
#line 321 "zgelsy.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 325 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 326 "zgelsy.f"
	iascl = 2;
#line 327 "zgelsy.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 331 "zgelsy.f"
	i__1 = max(*m,*n);
#line 331 "zgelsy.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 332 "zgelsy.f"
	*rank = 0;
#line 333 "zgelsy.f"
	goto L70;
#line 334 "zgelsy.f"
    }

#line 336 "zgelsy.f"
    bnrm = zlange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 337 "zgelsy.f"
    ibscl = 0;
#line 338 "zgelsy.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 342 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 343 "zgelsy.f"
	ibscl = 1;
#line 344 "zgelsy.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 348 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 349 "zgelsy.f"
	ibscl = 2;
#line 350 "zgelsy.f"
    }

/*     Compute QR factorization with column pivoting of A: */
/*        A * P = Q * R */

#line 355 "zgelsy.f"
    i__1 = *lwork - mn;
#line 355 "zgelsy.f"
    zgeqp3_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], &i__1,
	     &rwork[1], info);
#line 357 "zgelsy.f"
    i__1 = mn + 1;
#line 357 "zgelsy.f"
    wsize = mn + work[i__1].r;

/*     complex workspace: MN+NB*(N+1). real workspace 2*N. */
/*     Details of Householder rotations stored in WORK(1:MN). */

/*     Determine RANK using incremental condition estimation */

#line 364 "zgelsy.f"
    i__1 = ismin;
#line 364 "zgelsy.f"
    work[i__1].r = 1., work[i__1].i = 0.;
#line 365 "zgelsy.f"
    i__1 = ismax;
#line 365 "zgelsy.f"
    work[i__1].r = 1., work[i__1].i = 0.;
#line 366 "zgelsy.f"
    smax = z_abs(&a[a_dim1 + 1]);
#line 367 "zgelsy.f"
    smin = smax;
#line 368 "zgelsy.f"
    if (z_abs(&a[a_dim1 + 1]) == 0.) {
#line 369 "zgelsy.f"
	*rank = 0;
#line 370 "zgelsy.f"
	i__1 = max(*m,*n);
#line 370 "zgelsy.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 371 "zgelsy.f"
	goto L70;
#line 372 "zgelsy.f"
    } else {
#line 373 "zgelsy.f"
	*rank = 1;
#line 374 "zgelsy.f"
    }

#line 376 "zgelsy.f"
L10:
#line 377 "zgelsy.f"
    if (*rank < mn) {
#line 378 "zgelsy.f"
	i__ = *rank + 1;
#line 379 "zgelsy.f"
	zlaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 381 "zgelsy.f"
	zlaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 384 "zgelsy.f"
	if (smaxpr * *rcond <= sminpr) {
#line 385 "zgelsy.f"
	    i__1 = *rank;
#line 385 "zgelsy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "zgelsy.f"
		i__2 = ismin + i__ - 1;
#line 386 "zgelsy.f"
		i__3 = ismin + i__ - 1;
#line 386 "zgelsy.f"
		z__1.r = s1.r * work[i__3].r - s1.i * work[i__3].i, z__1.i = 
			s1.r * work[i__3].i + s1.i * work[i__3].r;
#line 386 "zgelsy.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 387 "zgelsy.f"
		i__2 = ismax + i__ - 1;
#line 387 "zgelsy.f"
		i__3 = ismax + i__ - 1;
#line 387 "zgelsy.f"
		z__1.r = s2.r * work[i__3].r - s2.i * work[i__3].i, z__1.i = 
			s2.r * work[i__3].i + s2.i * work[i__3].r;
#line 387 "zgelsy.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 388 "zgelsy.f"
/* L20: */
#line 388 "zgelsy.f"
	    }
#line 389 "zgelsy.f"
	    i__1 = ismin + *rank;
#line 389 "zgelsy.f"
	    work[i__1].r = c1.r, work[i__1].i = c1.i;
#line 390 "zgelsy.f"
	    i__1 = ismax + *rank;
#line 390 "zgelsy.f"
	    work[i__1].r = c2.r, work[i__1].i = c2.i;
#line 391 "zgelsy.f"
	    smin = sminpr;
#line 392 "zgelsy.f"
	    smax = smaxpr;
#line 393 "zgelsy.f"
	    ++(*rank);
#line 394 "zgelsy.f"
	    goto L10;
#line 395 "zgelsy.f"
	}
#line 396 "zgelsy.f"
    }

/*     complex workspace: 3*MN. */

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ] */
/*     where R11 = R(1:RANK,1:RANK) */

/*     [R11,R12] = [ T11, 0 ] * Y */

#line 406 "zgelsy.f"
    if (*rank < *n) {
#line 406 "zgelsy.f"
	i__1 = *lwork - (mn << 1);
#line 406 "zgelsy.f"
	ztzrzf_(rank, n, &a[a_offset], lda, &work[mn + 1], &work[(mn << 1) + 
		1], &i__1, info);
#line 406 "zgelsy.f"
    }

/*     complex workspace: 2*MN. */
/*     Details of Householder rotations stored in WORK(MN+1:2*MN) */

/*     B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */

#line 415 "zgelsy.f"
    i__1 = *lwork - (mn << 1);
#line 415 "zgelsy.f"
    zunmqr_("Left", "Conjugate transpose", m, nrhs, &mn, &a[a_offset], lda, &
	    work[1], &b[b_offset], ldb, &work[(mn << 1) + 1], &i__1, info, (
	    ftnlen)4, (ftnlen)19);
/* Computing MAX */
#line 417 "zgelsy.f"
    i__1 = (mn << 1) + 1;
#line 417 "zgelsy.f"
    d__1 = wsize, d__2 = (mn << 1) + work[i__1].r;
#line 417 "zgelsy.f"
    wsize = max(d__1,d__2);

/*     complex workspace: 2*MN+NB*NRHS. */

/*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

#line 423 "zgelsy.f"
    ztrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b2, &a[
	    a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);

#line 426 "zgelsy.f"
    i__1 = *nrhs;
#line 426 "zgelsy.f"
    for (j = 1; j <= i__1; ++j) {
#line 427 "zgelsy.f"
	i__2 = *n;
#line 427 "zgelsy.f"
	for (i__ = *rank + 1; i__ <= i__2; ++i__) {
#line 428 "zgelsy.f"
	    i__3 = i__ + j * b_dim1;
#line 428 "zgelsy.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 429 "zgelsy.f"
/* L30: */
#line 429 "zgelsy.f"
	}
#line 430 "zgelsy.f"
/* L40: */
#line 430 "zgelsy.f"
    }

/*     B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS) */

#line 434 "zgelsy.f"
    if (*rank < *n) {
#line 435 "zgelsy.f"
	i__1 = *n - *rank;
#line 435 "zgelsy.f"
	i__2 = *lwork - (mn << 1);
#line 435 "zgelsy.f"
	zunmrz_("Left", "Conjugate transpose", n, nrhs, rank, &i__1, &a[
		a_offset], lda, &work[mn + 1], &b[b_offset], ldb, &work[(mn <<
		 1) + 1], &i__2, info, (ftnlen)4, (ftnlen)19);
#line 438 "zgelsy.f"
    }

/*     complex workspace: 2*MN+NRHS. */

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

#line 444 "zgelsy.f"
    i__1 = *nrhs;
#line 444 "zgelsy.f"
    for (j = 1; j <= i__1; ++j) {
#line 445 "zgelsy.f"
	i__2 = *n;
#line 445 "zgelsy.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 446 "zgelsy.f"
	    i__3 = jpvt[i__];
#line 446 "zgelsy.f"
	    i__4 = i__ + j * b_dim1;
#line 446 "zgelsy.f"
	    work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 447 "zgelsy.f"
/* L50: */
#line 447 "zgelsy.f"
	}
#line 448 "zgelsy.f"
	zcopy_(n, &work[1], &c__1, &b[j * b_dim1 + 1], &c__1);
#line 449 "zgelsy.f"
/* L60: */
#line 449 "zgelsy.f"
    }

/*     complex workspace: N. */

/*     Undo scaling */

#line 455 "zgelsy.f"
    if (iascl == 1) {
#line 456 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 457 "zgelsy.f"
	zlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 459 "zgelsy.f"
    } else if (iascl == 2) {
#line 460 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 461 "zgelsy.f"
	zlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 463 "zgelsy.f"
    }
#line 464 "zgelsy.f"
    if (ibscl == 1) {
#line 465 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 466 "zgelsy.f"
    } else if (ibscl == 2) {
#line 467 "zgelsy.f"
	zlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 468 "zgelsy.f"
    }

#line 470 "zgelsy.f"
L70:
#line 471 "zgelsy.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 471 "zgelsy.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 473 "zgelsy.f"
    return 0;

/*     End of ZGELSY */

} /* zgelsy_ */

