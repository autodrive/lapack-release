#line 1 "cgels.f"
/* cgels.f -- translated by f2c (version 20100827).
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

#line 1 "cgels.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;

/* > \brief <b> CGELS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGELS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgels.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgels.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgels.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGELS solves overdetermined or underdetermined complex linear systems */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >            if M >= N, A is overwritten by details of its QR */
/* >                       factorization as returned by CGEQRF; */
/* >            if M <  N, A is overwritten by details of its LQ */
/* >                       factorization as returned by CGELQF. */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complexGEsolve */

/*  ===================================================================== */
/* Subroutine */ int cgels_(char *trans, integer *m, integer *n, integer *
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
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int cgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer scllen;
    static doublereal bignum;
    extern /* Subroutine */ int cunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal smlnum;
    static logical lquery;
    extern /* Subroutine */ int ctrtrs_(char *, char *, char *, integer *, 
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

#line 231 "cgels.f"
    /* Parameter adjustments */
#line 231 "cgels.f"
    a_dim1 = *lda;
#line 231 "cgels.f"
    a_offset = 1 + a_dim1;
#line 231 "cgels.f"
    a -= a_offset;
#line 231 "cgels.f"
    b_dim1 = *ldb;
#line 231 "cgels.f"
    b_offset = 1 + b_dim1;
#line 231 "cgels.f"
    b -= b_offset;
#line 231 "cgels.f"
    --work;
#line 231 "cgels.f"

#line 231 "cgels.f"
    /* Function Body */
#line 231 "cgels.f"
    *info = 0;
#line 232 "cgels.f"
    mn = min(*m,*n);
#line 233 "cgels.f"
    lquery = *lwork == -1;
#line 234 "cgels.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1))) {
#line 235 "cgels.f"
	*info = -1;
#line 236 "cgels.f"
    } else if (*m < 0) {
#line 237 "cgels.f"
	*info = -2;
#line 238 "cgels.f"
    } else if (*n < 0) {
#line 239 "cgels.f"
	*info = -3;
#line 240 "cgels.f"
    } else if (*nrhs < 0) {
#line 241 "cgels.f"
	*info = -4;
#line 242 "cgels.f"
    } else if (*lda < max(1,*m)) {
#line 243 "cgels.f"
	*info = -6;
#line 244 "cgels.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 244 "cgels.f"
	i__1 = max(1,*m);
#line 244 "cgels.f"
	if (*ldb < max(i__1,*n)) {
#line 245 "cgels.f"
	    *info = -8;
#line 246 "cgels.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 246 "cgels.f"
	    i__1 = 1, i__2 = mn + max(mn,*nrhs);
#line 246 "cgels.f"
	    if (*lwork < max(i__1,i__2) && ! lquery) {
#line 248 "cgels.f"
		*info = -10;
#line 249 "cgels.f"
	    }
#line 249 "cgels.f"
	}
#line 249 "cgels.f"
    }

/*     Figure out optimal block size */

#line 253 "cgels.f"
    if (*info == 0 || *info == -10) {

#line 255 "cgels.f"
	tpsd = TRUE_;
#line 256 "cgels.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 256 "cgels.f"
	    tpsd = FALSE_;
#line 256 "cgels.f"
	}

#line 259 "cgels.f"
	if (*m >= *n) {
#line 260 "cgels.f"
	    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 261 "cgels.f"
	    if (tpsd) {
/* Computing MAX */
#line 262 "cgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMQR", "LN", m, nrhs, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 262 "cgels.f"
		nb = max(i__1,i__2);
#line 264 "cgels.f"
	    } else {
/* Computing MAX */
#line 265 "cgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMQR", "LC", m, nrhs, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 265 "cgels.f"
		nb = max(i__1,i__2);
#line 267 "cgels.f"
	    }
#line 268 "cgels.f"
	} else {
#line 269 "cgels.f"
	    nb = ilaenv_(&c__1, "CGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 270 "cgels.f"
	    if (tpsd) {
/* Computing MAX */
#line 271 "cgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMLQ", "LC", n, nrhs, m, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 271 "cgels.f"
		nb = max(i__1,i__2);
#line 273 "cgels.f"
	    } else {
/* Computing MAX */
#line 274 "cgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "CUNMLQ", "LN", n, nrhs, m, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 274 "cgels.f"
		nb = max(i__1,i__2);
#line 276 "cgels.f"
	    }
#line 277 "cgels.f"
	}

/* Computing MAX */
#line 279 "cgels.f"
	i__1 = 1, i__2 = mn + max(mn,*nrhs) * nb;
#line 279 "cgels.f"
	wsize = max(i__1,i__2);
#line 280 "cgels.f"
	d__1 = (doublereal) wsize;
#line 280 "cgels.f"
	work[1].r = d__1, work[1].i = 0.;

#line 282 "cgels.f"
    }

#line 284 "cgels.f"
    if (*info != 0) {
#line 285 "cgels.f"
	i__1 = -(*info);
#line 285 "cgels.f"
	xerbla_("CGELS ", &i__1, (ftnlen)6);
#line 286 "cgels.f"
	return 0;
#line 287 "cgels.f"
    } else if (lquery) {
#line 288 "cgels.f"
	return 0;
#line 289 "cgels.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 293 "cgels.f"
    i__1 = min(*m,*n);
#line 293 "cgels.f"
    if (min(i__1,*nrhs) == 0) {
#line 294 "cgels.f"
	i__1 = max(*m,*n);
#line 294 "cgels.f"
	claset_("Full", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)
		4);
#line 295 "cgels.f"
	return 0;
#line 296 "cgels.f"
    }

/*     Get machine parameters */

#line 300 "cgels.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 301 "cgels.f"
    bignum = 1. / smlnum;
#line 302 "cgels.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 306 "cgels.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, rwork, (ftnlen)1);
#line 307 "cgels.f"
    iascl = 0;
#line 308 "cgels.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 312 "cgels.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 313 "cgels.f"
	iascl = 1;
#line 314 "cgels.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 318 "cgels.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 319 "cgels.f"
	iascl = 2;
#line 320 "cgels.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 324 "cgels.f"
	i__1 = max(*m,*n);
#line 324 "cgels.f"
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 325 "cgels.f"
	goto L50;
#line 326 "cgels.f"
    }

#line 328 "cgels.f"
    brow = *m;
#line 329 "cgels.f"
    if (tpsd) {
#line 329 "cgels.f"
	brow = *n;
#line 329 "cgels.f"
    }
#line 331 "cgels.f"
    bnrm = clange_("M", &brow, nrhs, &b[b_offset], ldb, rwork, (ftnlen)1);
#line 332 "cgels.f"
    ibscl = 0;
#line 333 "cgels.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 337 "cgels.f"
	clascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 339 "cgels.f"
	ibscl = 1;
#line 340 "cgels.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 344 "cgels.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 346 "cgels.f"
	ibscl = 2;
#line 347 "cgels.f"
    }

#line 349 "cgels.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 353 "cgels.f"
	i__1 = *lwork - mn;
#line 353 "cgels.f"
	cgeqrf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least N, optimally N*NB */

#line 358 "cgels.f"
	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */

#line 364 "cgels.f"
	    i__1 = *lwork - mn;
#line 364 "cgels.f"
	    cunmqr_("Left", "Conjugate transpose", m, nrhs, n, &a[a_offset], 
		    lda, &work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, 
		    info, (ftnlen)4, (ftnlen)19);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 372 "cgels.f"
	    ctrtrs_("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset]
		    , lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

#line 375 "cgels.f"
	    if (*info > 0) {
#line 376 "cgels.f"
		return 0;
#line 377 "cgels.f"
	    }

#line 379 "cgels.f"
	    scllen = *n;

#line 381 "cgels.f"
	} else {

/*           Underdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS) */

#line 387 "cgels.f"
	    ctrtrs_("Upper", "Conjugate transpose", "Non-unit", n, nrhs, &a[
		    a_offset], lda, &b[b_offset], ldb, info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

#line 390 "cgels.f"
	    if (*info > 0) {
#line 391 "cgels.f"
		return 0;
#line 392 "cgels.f"
	    }

/*           B(N+1:M,1:NRHS) = ZERO */

#line 396 "cgels.f"
	    i__1 = *nrhs;
#line 396 "cgels.f"
	    for (j = 1; j <= i__1; ++j) {
#line 397 "cgels.f"
		i__2 = *m;
#line 397 "cgels.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 398 "cgels.f"
		    i__3 = i__ + j * b_dim1;
#line 398 "cgels.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 399 "cgels.f"
/* L10: */
#line 399 "cgels.f"
		}
#line 400 "cgels.f"
/* L20: */
#line 400 "cgels.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 404 "cgels.f"
	    i__1 = *lwork - mn;
#line 404 "cgels.f"
	    cunmqr_("Left", "No transpose", m, nrhs, n, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)12);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 410 "cgels.f"
	    scllen = *m;

#line 412 "cgels.f"
	}

#line 414 "cgels.f"
    } else {

/*        Compute LQ factorization of A */

#line 418 "cgels.f"
	i__1 = *lwork - mn;
#line 418 "cgels.f"
	cgelqf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least M, optimally M*NB. */

#line 423 "cgels.f"
	if (! tpsd) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 429 "cgels.f"
	    ctrtrs_("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset]
		    , lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

#line 432 "cgels.f"
	    if (*info > 0) {
#line 433 "cgels.f"
		return 0;
#line 434 "cgels.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 438 "cgels.f"
	    i__1 = *nrhs;
#line 438 "cgels.f"
	    for (j = 1; j <= i__1; ++j) {
#line 439 "cgels.f"
		i__2 = *n;
#line 439 "cgels.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 440 "cgels.f"
		    i__3 = i__ + j * b_dim1;
#line 440 "cgels.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 441 "cgels.f"
/* L30: */
#line 441 "cgels.f"
		}
#line 442 "cgels.f"
/* L40: */
#line 442 "cgels.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS) */

#line 446 "cgels.f"
	    i__1 = *lwork - mn;
#line 446 "cgels.f"
	    cunmlq_("Left", "Conjugate transpose", n, nrhs, m, &a[a_offset], 
		    lda, &work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, 
		    info, (ftnlen)4, (ftnlen)19);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 452 "cgels.f"
	    scllen = *n;

#line 454 "cgels.f"
	} else {

/*           overdetermined system min || A**H * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 460 "cgels.f"
	    i__1 = *lwork - mn;
#line 460 "cgels.f"
	    cunmlq_("Left", "No transpose", n, nrhs, m, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)12);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS) */

#line 468 "cgels.f"
	    ctrtrs_("Lower", "Conjugate transpose", "Non-unit", m, nrhs, &a[
		    a_offset], lda, &b[b_offset], ldb, info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

#line 471 "cgels.f"
	    if (*info > 0) {
#line 472 "cgels.f"
		return 0;
#line 473 "cgels.f"
	    }

#line 475 "cgels.f"
	    scllen = *m;

#line 477 "cgels.f"
	}

#line 479 "cgels.f"
    }

/*     Undo scaling */

#line 483 "cgels.f"
    if (iascl == 1) {
#line 484 "cgels.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 486 "cgels.f"
    } else if (iascl == 2) {
#line 487 "cgels.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 489 "cgels.f"
    }
#line 490 "cgels.f"
    if (ibscl == 1) {
#line 491 "cgels.f"
	clascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 493 "cgels.f"
    } else if (ibscl == 2) {
#line 494 "cgels.f"
	clascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 496 "cgels.f"
    }

#line 498 "cgels.f"
L50:
#line 499 "cgels.f"
    d__1 = (doublereal) wsize;
#line 499 "cgels.f"
    work[1].r = d__1, work[1].i = 0.;

#line 501 "cgels.f"
    return 0;

/*     End of CGELS */

} /* cgels_ */

