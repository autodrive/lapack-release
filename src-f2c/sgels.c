#line 1 "sgels.f"
/* sgels.f -- translated by f2c (version 20100827).
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

#line 1 "sgels.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b33 = 0.;
static integer c__0 = 0;

/* > \brief <b> SGELS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGELS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgels.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgels.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgels.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGELS solves overdetermined or underdetermined real linear systems */
/* > involving an M-by-N matrix A, or its transpose, using a QR or LQ */
/* > factorization of A.  It is assumed that A has full rank. */
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
/* > 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of */
/* >    an undetermined system A**T * X = B. */
/* > */
/* > 4. If TRANS = 'T' and m < n:  find the least squares solution of */
/* >    an overdetermined system, i.e., solve the least squares problem */
/* >                 minimize || B - A**T * X ||. */
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
/* >          = 'T': the linear system involves A**T. */
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
/* >          columns of the matrices B and X. NRHS >=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >            if M >= N, A is overwritten by details of its QR */
/* >                       factorization as returned by SGEQRF; */
/* >            if M <  N, A is overwritten by details of its LQ */
/* >                       factorization as returned by SGELQF. */
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
/* >          On entry, the matrix B of right hand side vectors, stored */
/* >          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* >          if TRANS = 'T'. */
/* >          On exit, if INFO = 0, B is overwritten by the solution */
/* >          vectors, stored columnwise: */
/* >          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* >          squares solution vectors; the residual sum of squares for the */
/* >          solution in each column is given by the sum of squares of */
/* >          elements N+1 to M in that column; */
/* >          if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'T' and m >= n, rows 1 to M of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'T' and m < n, rows 1 to M of B contain the */
/* >          least squares solution vectors; the residual sum of squares */
/* >          for the solution in each column is given by the sum of */
/* >          squares of elements M+1 to N in that column. */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \date November 2011 */

/* > \ingroup realGEsolve */

/*  ===================================================================== */
/* Subroutine */ int sgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

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
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer scllen;
    static doublereal bignum;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     sgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), slaset_(char *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int sormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    strtrs_(char *, char *, char *, integer *, integer *, doublereal *
	    , integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
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

#line 230 "sgels.f"
    /* Parameter adjustments */
#line 230 "sgels.f"
    a_dim1 = *lda;
#line 230 "sgels.f"
    a_offset = 1 + a_dim1;
#line 230 "sgels.f"
    a -= a_offset;
#line 230 "sgels.f"
    b_dim1 = *ldb;
#line 230 "sgels.f"
    b_offset = 1 + b_dim1;
#line 230 "sgels.f"
    b -= b_offset;
#line 230 "sgels.f"
    --work;
#line 230 "sgels.f"

#line 230 "sgels.f"
    /* Function Body */
#line 230 "sgels.f"
    *info = 0;
#line 231 "sgels.f"
    mn = min(*m,*n);
#line 232 "sgels.f"
    lquery = *lwork == -1;
#line 233 "sgels.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1))) {
#line 234 "sgels.f"
	*info = -1;
#line 235 "sgels.f"
    } else if (*m < 0) {
#line 236 "sgels.f"
	*info = -2;
#line 237 "sgels.f"
    } else if (*n < 0) {
#line 238 "sgels.f"
	*info = -3;
#line 239 "sgels.f"
    } else if (*nrhs < 0) {
#line 240 "sgels.f"
	*info = -4;
#line 241 "sgels.f"
    } else if (*lda < max(1,*m)) {
#line 242 "sgels.f"
	*info = -6;
#line 243 "sgels.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 243 "sgels.f"
	i__1 = max(1,*m);
#line 243 "sgels.f"
	if (*ldb < max(i__1,*n)) {
#line 244 "sgels.f"
	    *info = -8;
#line 245 "sgels.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 245 "sgels.f"
	    i__1 = 1, i__2 = mn + max(mn,*nrhs);
#line 245 "sgels.f"
	    if (*lwork < max(i__1,i__2) && ! lquery) {
#line 247 "sgels.f"
		*info = -10;
#line 248 "sgels.f"
	    }
#line 248 "sgels.f"
	}
#line 248 "sgels.f"
    }

/*     Figure out optimal block size */

#line 252 "sgels.f"
    if (*info == 0 || *info == -10) {

#line 254 "sgels.f"
	tpsd = TRUE_;
#line 255 "sgels.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 255 "sgels.f"
	    tpsd = FALSE_;
#line 255 "sgels.f"
	}

#line 258 "sgels.f"
	if (*m >= *n) {
#line 259 "sgels.f"
	    nb = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 260 "sgels.f"
	    if (tpsd) {
/* Computing MAX */
#line 261 "sgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "SORMQR", "LN", m, nrhs, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 261 "sgels.f"
		nb = max(i__1,i__2);
#line 263 "sgels.f"
	    } else {
/* Computing MAX */
#line 264 "sgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "SORMQR", "LT", m, nrhs, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 264 "sgels.f"
		nb = max(i__1,i__2);
#line 266 "sgels.f"
	    }
#line 267 "sgels.f"
	} else {
#line 268 "sgels.f"
	    nb = ilaenv_(&c__1, "SGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 269 "sgels.f"
	    if (tpsd) {
/* Computing MAX */
#line 270 "sgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "SORMLQ", "LT", n, nrhs, m, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 270 "sgels.f"
		nb = max(i__1,i__2);
#line 272 "sgels.f"
	    } else {
/* Computing MAX */
#line 273 "sgels.f"
		i__1 = nb, i__2 = ilaenv_(&c__1, "SORMLQ", "LN", n, nrhs, m, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 273 "sgels.f"
		nb = max(i__1,i__2);
#line 275 "sgels.f"
	    }
#line 276 "sgels.f"
	}

/* Computing MAX */
#line 278 "sgels.f"
	i__1 = 1, i__2 = mn + max(mn,*nrhs) * nb;
#line 278 "sgels.f"
	wsize = max(i__1,i__2);
#line 279 "sgels.f"
	work[1] = (doublereal) wsize;

#line 281 "sgels.f"
    }

#line 283 "sgels.f"
    if (*info != 0) {
#line 284 "sgels.f"
	i__1 = -(*info);
#line 284 "sgels.f"
	xerbla_("SGELS ", &i__1, (ftnlen)6);
#line 285 "sgels.f"
	return 0;
#line 286 "sgels.f"
    } else if (lquery) {
#line 287 "sgels.f"
	return 0;
#line 288 "sgels.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 292 "sgels.f"
    i__1 = min(*m,*n);
#line 292 "sgels.f"
    if (min(i__1,*nrhs) == 0) {
#line 293 "sgels.f"
	i__1 = max(*m,*n);
#line 293 "sgels.f"
	slaset_("Full", &i__1, nrhs, &c_b33, &c_b33, &b[b_offset], ldb, (
		ftnlen)4);
#line 294 "sgels.f"
	return 0;
#line 295 "sgels.f"
    }

/*     Get machine parameters */

#line 299 "sgels.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 300 "sgels.f"
    bignum = 1. / smlnum;
#line 301 "sgels.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 305 "sgels.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, rwork, (ftnlen)1);
#line 306 "sgels.f"
    iascl = 0;
#line 307 "sgels.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 311 "sgels.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 312 "sgels.f"
	iascl = 1;
#line 313 "sgels.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 317 "sgels.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 318 "sgels.f"
	iascl = 2;
#line 319 "sgels.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 323 "sgels.f"
	i__1 = max(*m,*n);
#line 323 "sgels.f"
	slaset_("F", &i__1, nrhs, &c_b33, &c_b33, &b[b_offset], ldb, (ftnlen)
		1);
#line 324 "sgels.f"
	goto L50;
#line 325 "sgels.f"
    }

#line 327 "sgels.f"
    brow = *m;
#line 328 "sgels.f"
    if (tpsd) {
#line 328 "sgels.f"
	brow = *n;
#line 328 "sgels.f"
    }
#line 330 "sgels.f"
    bnrm = slange_("M", &brow, nrhs, &b[b_offset], ldb, rwork, (ftnlen)1);
#line 331 "sgels.f"
    ibscl = 0;
#line 332 "sgels.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 336 "sgels.f"
	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 338 "sgels.f"
	ibscl = 1;
#line 339 "sgels.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 343 "sgels.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 345 "sgels.f"
	ibscl = 2;
#line 346 "sgels.f"
    }

#line 348 "sgels.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 352 "sgels.f"
	i__1 = *lwork - mn;
#line 352 "sgels.f"
	sgeqrf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least N, optimally N*NB */

#line 357 "sgels.f"
	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 363 "sgels.f"
	    i__1 = *lwork - mn;
#line 363 "sgels.f"
	    sormqr_("Left", "Transpose", m, nrhs, n, &a[a_offset], lda, &work[
		    1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)9);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 371 "sgels.f"
	    strtrs_("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset]
		    , lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

#line 374 "sgels.f"
	    if (*info > 0) {
#line 375 "sgels.f"
		return 0;
#line 376 "sgels.f"
	    }

#line 378 "sgels.f"
	    scllen = *n;

#line 380 "sgels.f"
	} else {

/*           Overdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */

#line 386 "sgels.f"
	    strtrs_("Upper", "Transpose", "Non-unit", n, nrhs, &a[a_offset], 
		    lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);

#line 389 "sgels.f"
	    if (*info > 0) {
#line 390 "sgels.f"
		return 0;
#line 391 "sgels.f"
	    }

/*           B(N+1:M,1:NRHS) = ZERO */

#line 395 "sgels.f"
	    i__1 = *nrhs;
#line 395 "sgels.f"
	    for (j = 1; j <= i__1; ++j) {
#line 396 "sgels.f"
		i__2 = *m;
#line 396 "sgels.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 397 "sgels.f"
		    b[i__ + j * b_dim1] = 0.;
#line 398 "sgels.f"
/* L10: */
#line 398 "sgels.f"
		}
#line 399 "sgels.f"
/* L20: */
#line 399 "sgels.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 403 "sgels.f"
	    i__1 = *lwork - mn;
#line 403 "sgels.f"
	    sormqr_("Left", "No transpose", m, nrhs, n, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)12);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 409 "sgels.f"
	    scllen = *m;

#line 411 "sgels.f"
	}

#line 413 "sgels.f"
    } else {

/*        Compute LQ factorization of A */

#line 417 "sgels.f"
	i__1 = *lwork - mn;
#line 417 "sgels.f"
	sgelqf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least M, optimally M*NB. */

#line 422 "sgels.f"
	if (! tpsd) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 428 "sgels.f"
	    strtrs_("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset]
		    , lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

#line 431 "sgels.f"
	    if (*info > 0) {
#line 432 "sgels.f"
		return 0;
#line 433 "sgels.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 437 "sgels.f"
	    i__1 = *nrhs;
#line 437 "sgels.f"
	    for (j = 1; j <= i__1; ++j) {
#line 438 "sgels.f"
		i__2 = *n;
#line 438 "sgels.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 439 "sgels.f"
		    b[i__ + j * b_dim1] = 0.;
#line 440 "sgels.f"
/* L30: */
#line 440 "sgels.f"
		}
#line 441 "sgels.f"
/* L40: */
#line 441 "sgels.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */

#line 445 "sgels.f"
	    i__1 = *lwork - mn;
#line 445 "sgels.f"
	    sormlq_("Left", "Transpose", n, nrhs, m, &a[a_offset], lda, &work[
		    1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)9);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 451 "sgels.f"
	    scllen = *n;

#line 453 "sgels.f"
	} else {

/*           overdetermined system min || A**T * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 459 "sgels.f"
	    i__1 = *lwork - mn;
#line 459 "sgels.f"
	    sormlq_("Left", "No transpose", n, nrhs, m, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info, (
		    ftnlen)4, (ftnlen)12);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */

#line 467 "sgels.f"
	    strtrs_("Lower", "Transpose", "Non-unit", m, nrhs, &a[a_offset], 
		    lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);

#line 470 "sgels.f"
	    if (*info > 0) {
#line 471 "sgels.f"
		return 0;
#line 472 "sgels.f"
	    }

#line 474 "sgels.f"
	    scllen = *m;

#line 476 "sgels.f"
	}

#line 478 "sgels.f"
    }

/*     Undo scaling */

#line 482 "sgels.f"
    if (iascl == 1) {
#line 483 "sgels.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 485 "sgels.f"
    } else if (iascl == 2) {
#line 486 "sgels.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 488 "sgels.f"
    }
#line 489 "sgels.f"
    if (ibscl == 1) {
#line 490 "sgels.f"
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 492 "sgels.f"
    } else if (ibscl == 2) {
#line 493 "sgels.f"
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 495 "sgels.f"
    }

#line 497 "sgels.f"
L50:
#line 498 "sgels.f"
    work[1] = (doublereal) wsize;

#line 500 "sgels.f"
    return 0;

/*     End of SGELS */

} /* sgels_ */

