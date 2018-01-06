#line 1 "zgels.f"
/* zgels.f -- translated by f2c (version 20100827).
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

#line 1 "zgels.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;

/* > \brief <b> ZGELS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGELS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgels.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgels.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgels.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGELS solves overdetermined or underdetermined complex linear systems */
/* > involving an M-by-N matrix A, or its conjugate-transpose, using a QR */
/* > or LQ factorization of A.  It is assumed that A has full rank. */
/* > */
/* > The following options are provided: */
/* > */
/* > 1. If TRANS = 'N' and m >= n:  find the least squares solution of */
/* >    an overdetermined system, i.e., solve the least squares problem */
/* >                 minimize || B - A*X ||. */
/* > */
/* > 2. If TRANS = 'N' and m < n:  find the minimum norm solution of */
/* >    an underdetermined system A * X = B. */
/* > */
/* > 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of */
/* >    an underdetermined system A**H * X = B. */
/* > */
/* > 4. If TRANS = 'C' and m < n:  find the least squares solution of */
/* >    an overdetermined system, i.e., solve the least squares problem */
/* >                 minimize || B - A**H * X ||. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': the linear system involves A; */
/* >          = 'C': the linear system involves A**H. */
/* > \endverbatim */
/* > */
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
/* >          columns of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >            if M >= N, A is overwritten by details of its QR */
/* >                       factorization as returned by ZGEQRF; */
/* >            if M <  N, A is overwritten by details of its LQ */
/* >                       factorization as returned by ZGELQF. */
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
/* >          On entry, the matrix B of right hand side vectors, stored */
/* >          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* >          if TRANS = 'C'. */
/* >          On exit, if INFO = 0, B is overwritten by the solution */
/* >          vectors, stored columnwise: */
/* >          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* >          squares solution vectors; the residual sum of squares for the */
/* >          solution in each column is given by the sum of squares of the */
/* >          modulus of elements N+1 to M in that column; */
/* >          if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'C' and m >= n, rows 1 to M of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'C' and m < n, rows 1 to M of B contain the */
/* >          least squares solution vectors; the residual sum of squares */
/* >          for the solution in each column is given by the sum of */
/* >          squares of the modulus of elements M+1 to N in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= MAX(1,M,N). */
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
/* >          LWORK >= max( 1, MN + max( MN, NRHS ) ). */
/* >          For optimal performance, */
/* >          LWORK >= max( 1, MN + max( MN, NRHS )*NB ). */
/* >          where MN = min(M,N) and NB is the optimum block size. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO =  i, the i-th diagonal element of the */
/* >                triangular factor of A is zero, so that A does not have */
/* >                full rank; the least squares solution could not be */
/* >                computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GEsolve */

/*  ===================================================================== */
/* Subroutine */ int zgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, nb, mn;
    static doublereal anrm, bnrm;
    static integer brow;
    static logical tpsd;
    static integer iascl, ibscl;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer wsize;
    static doublereal rwork[1];
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer scllen;
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), zgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *), zlaset_(
	    char *, integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen);
    static doublereal smlnum;
    static logical lquery;
    extern /* Subroutine */ int zunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), ztrtrs_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 231 "zgels.f"
    /* Parameter adjustments */
#line 231 "zgels.f"
    a_dim1 = *lda;
#line 231 "zgels.f"
    a_offset = 1 + a_dim1;
#line 231 "zgels.f"
    a -= a_offset;
#line 231 "zgels.f"
    b_dim1 = *ldb;
#line 231 "zgels.f"
    b_offset = 1 + b_dim1;
#line 231 "zgels.f"
    b -= b_offset;
#line 231 "zgels.f"
    --work;
#line 231 "zgels.f"

#line 231 "zgels.f"
    /* Function Body */
#line 231 "zgels.f"
    *info = 0;
#line 232 "zgels.f"
    mn = min(*m,*n);
#line 233 "zgels.f"
    lquery = *lwork == -1;
#line 234 "zgels.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1))) {
#line 235 "zgels.f"
	*info = -1;
#line 236 "zgels.f"
    } else if (*m < 0) {
#line 237 "zgels.f"
	*info = -2;
#line 238 "zgels.f"
    } else if (*n < 0) {
#line 239 "zgels.f"
	*info = -3;
#line 240 "zgels.f"
    } else if (*nrhs < 0) {
#line 241 "zgels.f"
	*info = -4;
#line 242 "zgels.f"
    } else if (*lda < max(1,*m)) {
#line 243 "zgels.f"
	*info = -6;
#line 244 "zgels.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 244 "zgels.f"
	i__1 = max(1,*m);
#line 244 "zgels.f"
	if (*ldb < max(i__1,*n)) {
#line 245 "zgels.f"
	    *info = -8;
#line 246 "zgels.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 246 "zgels.f"
	    i__1 = 1, i__2 = mn + max(mn,*nrhs);
#line 246 "zgels.f"
	    if (*lwork < max(i__1,i__2) && ! lquery) {
#line 248 "zgels.f"
		*info = -10;
#line 249 "zgels.f"
	    }
#line 249 "zgels.f"
	}
#line 249 "zgels.f"
    }

/*     Figure out optimal block size */

#line 253 "zgels.f"
    if (*info == 0 || *info == -10) {

#line 255 "zgels.f"
	tpsd = TRUE_;
#line 256 "zgels.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 256 "zgels.f"
	    tpsd = FALSE_;
#line 256 "zgels.f"
	}

#line 259 "zgels.f"
	if (*m >= *n) {
#line 260 "zgels.f"
	    nb = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 261 "zgels.f"
	    if (tpsd) {
/* Computing MAX */
#line 262 "zgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "ZUNMQR", "LN", m, nrhs, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 262 "zgels.f"
		nb = max(i__1,i__2);
#line 264 "zgels.f"
	    } else {
/* Computing MAX */
#line 265 "zgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", m, nrhs, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 265 "zgels.f"
		nb = max(i__1,i__2);
#line 267 "zgels.f"
	    }
#line 268 "zgels.f"
	} else {
#line 269 "zgels.f"
	    nb = ilaenv_(&c__1, "ZGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 270 "zgels.f"
	    if (tpsd) {
/* Computing MAX */
#line 271 "zgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "ZUNMLQ", "LC", n, nrhs, m, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 271 "zgels.f"
		nb = max(i__1,i__2);
#line 273 "zgels.f"
	    } else {
/* Computing MAX */
#line 274 "zgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "ZUNMLQ", "LN", n, nrhs, m, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 274 "zgels.f"
		nb = max(i__1,i__2);
#line 276 "zgels.f"
	    }
#line 277 "zgels.f"
	}

/* Computing MAX */
#line 279 "zgels.f"
	i__1 = 1, i__2 = mn + max(mn,*nrhs) * nb;
#line 279 "zgels.f"
	wsize = max(i__1,i__2);
#line 280 "zgels.f"
	d__1 = (doublereal) wsize;
#line 280 "zgels.f"
	work[1].r = d__1, work[1].i = 0.;

#line 282 "zgels.f"
    }

#line 284 "zgels.f"
    if (*info != 0) {
#line 285 "zgels.f"
	i__1 = -(*info);
#line 285 "zgels.f"
	xerbla_("ZGELS ", &i__1, (ftnlen)6);
#line 286 "zgels.f"
	return 0;
#line 287 "zgels.f"
    } else if (lquery) {
#line 288 "zgels.f"
	return 0;
#line 289 "zgels.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 293 "zgels.f"
    i__1 = min(*m,*n);
#line 293 "zgels.f"
    if (min(i__1,*nrhs) == 0) {
#line 294 "zgels.f"
	i__1 = max(*m,*n);
#line 294 "zgels.f"
	zlaset_("Full", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)
		4);
#line 295 "zgels.f"
	return 0;
#line 296 "zgels.f"
    }

/*     Get machine parameters */

#line 300 "zgels.f"
    smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 301 "zgels.f"
    bignum = 1. / smlnum;
#line 302 "zgels.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 306 "zgels.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, rwork, (ftnlen)1);
#line 307 "zgels.f"
    iascl = 0;
#line 308 "zgels.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 312 "zgels.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 313 "zgels.f"
	iascl = 1;
#line 314 "zgels.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 318 "zgels.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 319 "zgels.f"
	iascl = 2;
#line 320 "zgels.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 324 "zgels.f"
	i__1 = max(*m,*n);
#line 324 "zgels.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 325 "zgels.f"
	goto L50;
#line 326 "zgels.f"
    }

#line 328 "zgels.f"
    brow = *m;
#line 329 "zgels.f"
    if (tpsd) {
#line 329 "zgels.f"
	brow = *n;
#line 329 "zgels.f"
    }
#line 331 "zgels.f"
    bnrm = zlange_("M", &brow, nrhs, &b[b_offset], ldb, rwork, (ftnlen)1);
#line 332 "zgels.f"
    ibscl = 0;
#line 333 "zgels.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 337 "zgels.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 339 "zgels.f"
	ibscl = 1;
#line 340 "zgels.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 344 "zgels.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 346 "zgels.f"
	ibscl = 2;
#line 347 "zgels.f"
    }

#line 349 "zgels.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 353 "zgels.f"
	i__1 = *lwork - mn;
#line 353 "zgels.f"
	zgeqrf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least N, optimally N*NB */

#line 358 "zgels.f"
	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */

#line 364 "zgels.f"
	    i__1 = *lwork - mn;
#line 364 "zgels.f"
	    zunmqr_("Left", "Conjugate transpose", m, nrhs, n, &a[a_offset], 
		    lda, &work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, 
		    info, (ftnlen)4, (ftnlen)19);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 372 "zgels.f"
	    ztrtrs_("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset]
		    , lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

#line 375 "zgels.f"
	    if (*info > 0) {
#line 376 "zgels.f"
		return 0;
#line 377 "zgels.f"
	    }

#line 379 "zgels.f"
	    scllen = *n;

#line 381 "zgels.f"
	} else {

/*           Underdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS) */

#line 387 "zgels.f"
	    ztrtrs_("Upper", "Conjugate transpose", "Non-unit", n, nrhs, &a[
		    a_offset], lda, &b[b_offset], ldb, info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

#line 390 "zgels.f"
	    if (*info > 0) {
#line 391 "zgels.f"
		return 0;
#line 392 "zgels.f"
	    }

/*           B(N+1:M,1:NRHS) = ZERO */

#line 396 "zgels.f"
	    i__1 = *nrhs;
#line 396 "zgels.f"
	    for (j = 1; j <= i__1; ++j) {
#line 397 "zgels.f"
		i__2 = *m;
#line 397 "zgels.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 398 "zgels.f"
		    i__3 = i__ + j * b_dim1;
#line 398 "zgels.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 399 "zgels.f"
/* L10: */
#line 399 "zgels.f"
		}
#line 400 "zgels.f"
/* L20: */
#line 400 "zgels.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 404 "zgels.f"
	    i__1 = *lwork - mn;
#line 404 "zgels.f"
	    zunmqr_("Left", "No transpose", m, nrhs, n, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)12);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 410 "zgels.f"
	    scllen = *m;

#line 412 "zgels.f"
	}

#line 414 "zgels.f"
    } else {

/*        Compute LQ factorization of A */

#line 418 "zgels.f"
	i__1 = *lwork - mn;
#line 418 "zgels.f"
	zgelqf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least M, optimally M*NB. */

#line 423 "zgels.f"
	if (! tpsd) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 429 "zgels.f"
	    ztrtrs_("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset]
		    , lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

#line 432 "zgels.f"
	    if (*info > 0) {
#line 433 "zgels.f"
		return 0;
#line 434 "zgels.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 438 "zgels.f"
	    i__1 = *nrhs;
#line 438 "zgels.f"
	    for (j = 1; j <= i__1; ++j) {
#line 439 "zgels.f"
		i__2 = *n;
#line 439 "zgels.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 440 "zgels.f"
		    i__3 = i__ + j * b_dim1;
#line 440 "zgels.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 441 "zgels.f"
/* L30: */
#line 441 "zgels.f"
		}
#line 442 "zgels.f"
/* L40: */
#line 442 "zgels.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS) */

#line 446 "zgels.f"
	    i__1 = *lwork - mn;
#line 446 "zgels.f"
	    zunmlq_("Left", "Conjugate transpose", n, nrhs, m, &a[a_offset], 
		    lda, &work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, 
		    info, (ftnlen)4, (ftnlen)19);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 452 "zgels.f"
	    scllen = *n;

#line 454 "zgels.f"
	} else {

/*           overdetermined system min || A**H * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 460 "zgels.f"
	    i__1 = *lwork - mn;
#line 460 "zgels.f"
	    zunmlq_("Left", "No transpose", n, nrhs, m, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)12);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS) */

#line 468 "zgels.f"
	    ztrtrs_("Lower", "Conjugate transpose", "Non-unit", m, nrhs, &a[
		    a_offset], lda, &b[b_offset], ldb, info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

#line 471 "zgels.f"
	    if (*info > 0) {
#line 472 "zgels.f"
		return 0;
#line 473 "zgels.f"
	    }

#line 475 "zgels.f"
	    scllen = *m;

#line 477 "zgels.f"
	}

#line 479 "zgels.f"
    }

/*     Undo scaling */

#line 483 "zgels.f"
    if (iascl == 1) {
#line 484 "zgels.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 486 "zgels.f"
    } else if (iascl == 2) {
#line 487 "zgels.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 489 "zgels.f"
    }
#line 490 "zgels.f"
    if (ibscl == 1) {
#line 491 "zgels.f"
	zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 493 "zgels.f"
    } else if (ibscl == 2) {
#line 494 "zgels.f"
	zlascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 496 "zgels.f"
    }

#line 498 "zgels.f"
L50:
#line 499 "zgels.f"
    d__1 = (doublereal) wsize;
#line 499 "zgels.f"
    work[1].r = d__1, work[1].i = 0.;

#line 501 "zgels.f"
    return 0;

/*     End of ZGELS */

} /* zgels_ */

