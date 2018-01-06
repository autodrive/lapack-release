#line 1 "dgetsls.f"
/* dgetsls.f -- translated by f2c (version 20100827).
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

#line 1 "dgetsls.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c_n2 = -2;
static doublereal c_b23 = 0.;
static integer c__0 = 0;

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, */
/*     $                     WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGETSLS solves overdetermined or underdetermined real linear systems */
/* > involving an M-by-N matrix A, using a tall skinny QR or short wide LQ */
/* > factorization of A.  It is assumed that A has full rank. */
/* > */
/* > */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          A is overwritten by details of its QR or LQ */
/* >          factorization as returned by DGEQR or DGELQ. */
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
/* >          On entry, the matrix B of right hand side vectors, stored */
/* >          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* >          if TRANS = 'T'. */
/* >          On exit, if INFO = 0, B is overwritten by the solution */
/* >          vectors, stored columnwise: */
/* >          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* >          squares solution vectors. */
/* >          if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'T' and m >= n, rows 1 to M of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'T' and m < n, rows 1 to M of B contain the */
/* >          least squares solution vectors. */
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
/* >          (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal */
/* >          or optimal, if query was assumed) LWORK. */
/* >          See LWORK for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If LWORK = -1 or -2, then a workspace query is assumed. */
/* >          If LWORK = -1, the routine calculates optimal size of WORK for the */
/* >          optimal performance and returns this value in WORK(1). */
/* >          If LWORK = -2, the routine calculates minimal size of WORK and */
/* >          returns this value in WORK(1). */
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

/* > \ingroup doubleGEsolve */

/*  ===================================================================== */
/* Subroutine */ int dgetsls_(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal tq[5];
    static integer lw1, lw2, mnk, lwm, lwo;
    static doublereal anrm, bnrm;
    static logical tran;
    static integer brow, tszm, tszo, info2, iascl, ibscl;
    extern /* Subroutine */ int dgelq_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgeqr_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static integer minmn, maxmn;
    static doublereal workq;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgemlq_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dgemqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer scllen;
    static doublereal bignum, smlnum;
    static integer wsizem, wsizeo;
    static logical lquery;
    extern /* Subroutine */ int dtrtrs_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 207 "dgetsls.f"
    /* Parameter adjustments */
#line 207 "dgetsls.f"
    a_dim1 = *lda;
#line 207 "dgetsls.f"
    a_offset = 1 + a_dim1;
#line 207 "dgetsls.f"
    a -= a_offset;
#line 207 "dgetsls.f"
    b_dim1 = *ldb;
#line 207 "dgetsls.f"
    b_offset = 1 + b_dim1;
#line 207 "dgetsls.f"
    b -= b_offset;
#line 207 "dgetsls.f"
    --work;
#line 207 "dgetsls.f"

#line 207 "dgetsls.f"
    /* Function Body */
#line 207 "dgetsls.f"
    *info = 0;
#line 208 "dgetsls.f"
    minmn = min(*m,*n);
#line 209 "dgetsls.f"
    maxmn = max(*m,*n);
#line 210 "dgetsls.f"
    mnk = max(minmn,*nrhs);
#line 211 "dgetsls.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

#line 213 "dgetsls.f"
    lquery = *lwork == -1 || *lwork == -2;
#line 214 "dgetsls.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1))) {
#line 216 "dgetsls.f"
	*info = -1;
#line 217 "dgetsls.f"
    } else if (*m < 0) {
#line 218 "dgetsls.f"
	*info = -2;
#line 219 "dgetsls.f"
    } else if (*n < 0) {
#line 220 "dgetsls.f"
	*info = -3;
#line 221 "dgetsls.f"
    } else if (*nrhs < 0) {
#line 222 "dgetsls.f"
	*info = -4;
#line 223 "dgetsls.f"
    } else if (*lda < max(1,*m)) {
#line 224 "dgetsls.f"
	*info = -6;
#line 225 "dgetsls.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 225 "dgetsls.f"
	i__1 = max(1,*m);
#line 225 "dgetsls.f"
	if (*ldb < max(i__1,*n)) {
#line 226 "dgetsls.f"
	    *info = -8;
#line 227 "dgetsls.f"
	}
#line 227 "dgetsls.f"
    }

#line 229 "dgetsls.f"
    if (*info == 0) {

/*     Determine the block size and minimum LWORK */

#line 233 "dgetsls.f"
	if (*m >= *n) {
#line 234 "dgetsls.f"
	    dgeqr_(m, n, &a[a_offset], lda, tq, &c_n1, &workq, &c_n1, &info2);
#line 235 "dgetsls.f"
	    tszo = (integer) tq[0];
#line 236 "dgetsls.f"
	    lwo = (integer) workq;
#line 237 "dgetsls.f"
	    dgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 239 "dgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq;
#line 239 "dgetsls.f"
	    lwo = max(i__1,i__2);
#line 240 "dgetsls.f"
	    dgeqr_(m, n, &a[a_offset], lda, tq, &c_n2, &workq, &c_n2, &info2);
#line 241 "dgetsls.f"
	    tszm = (integer) tq[0];
#line 242 "dgetsls.f"
	    lwm = (integer) workq;
#line 243 "dgetsls.f"
	    dgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszm, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 245 "dgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq;
#line 245 "dgetsls.f"
	    lwm = max(i__1,i__2);
#line 246 "dgetsls.f"
	    wsizeo = tszo + lwo;
#line 247 "dgetsls.f"
	    wsizem = tszm + lwm;
#line 248 "dgetsls.f"
	} else {
#line 249 "dgetsls.f"
	    dgelq_(m, n, &a[a_offset], lda, tq, &c_n1, &workq, &c_n1, &info2);
#line 250 "dgetsls.f"
	    tszo = (integer) tq[0];
#line 251 "dgetsls.f"
	    lwo = (integer) workq;
#line 252 "dgetsls.f"
	    dgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 254 "dgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq;
#line 254 "dgetsls.f"
	    lwo = max(i__1,i__2);
#line 255 "dgetsls.f"
	    dgelq_(m, n, &a[a_offset], lda, tq, &c_n2, &workq, &c_n2, &info2);
#line 256 "dgetsls.f"
	    tszm = (integer) tq[0];
#line 257 "dgetsls.f"
	    lwm = (integer) workq;
#line 258 "dgetsls.f"
	    dgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 260 "dgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq;
#line 260 "dgetsls.f"
	    lwm = max(i__1,i__2);
#line 261 "dgetsls.f"
	    wsizeo = tszo + lwo;
#line 262 "dgetsls.f"
	    wsizem = tszm + lwm;
#line 263 "dgetsls.f"
	}

#line 265 "dgetsls.f"
	if (*lwork < wsizem && ! lquery) {
#line 266 "dgetsls.f"
	    *info = -10;
#line 267 "dgetsls.f"
	}

#line 269 "dgetsls.f"
    }

#line 271 "dgetsls.f"
    if (*info != 0) {
#line 272 "dgetsls.f"
	i__1 = -(*info);
#line 272 "dgetsls.f"
	xerbla_("DGETSLS", &i__1, (ftnlen)7);
#line 273 "dgetsls.f"
	work[1] = (doublereal) wsizeo;
#line 274 "dgetsls.f"
	return 0;
#line 275 "dgetsls.f"
    }
#line 276 "dgetsls.f"
    if (lquery) {
#line 277 "dgetsls.f"
	if (*lwork == -1) {
#line 277 "dgetsls.f"
	    work[1] = (doublereal) wsizeo;
#line 277 "dgetsls.f"
	}
#line 278 "dgetsls.f"
	if (*lwork == -2) {
#line 278 "dgetsls.f"
	    work[1] = (doublereal) wsizem;
#line 278 "dgetsls.f"
	}
#line 279 "dgetsls.f"
	return 0;
#line 280 "dgetsls.f"
    }
#line 281 "dgetsls.f"
    if (*lwork < wsizeo) {
#line 282 "dgetsls.f"
	lw1 = tszm;
#line 283 "dgetsls.f"
	lw2 = lwm;
#line 284 "dgetsls.f"
    } else {
#line 285 "dgetsls.f"
	lw1 = tszo;
#line 286 "dgetsls.f"
	lw2 = lwo;
#line 287 "dgetsls.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 291 "dgetsls.f"
    i__1 = min(*m,*n);
#line 291 "dgetsls.f"
    if (min(i__1,*nrhs) == 0) {
#line 292 "dgetsls.f"
	i__1 = max(*m,*n);
#line 292 "dgetsls.f"
	dlaset_("FULL", &i__1, nrhs, &c_b23, &c_b23, &b[b_offset], ldb, (
		ftnlen)4);
#line 294 "dgetsls.f"
	return 0;
#line 295 "dgetsls.f"
    }

/*     Get machine parameters */

#line 299 "dgetsls.f"
    smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 300 "dgetsls.f"
    bignum = 1. / smlnum;
#line 301 "dgetsls.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 305 "dgetsls.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 306 "dgetsls.f"
    iascl = 0;
#line 307 "dgetsls.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 311 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 312 "dgetsls.f"
	iascl = 1;
#line 313 "dgetsls.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 317 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 318 "dgetsls.f"
	iascl = 2;
#line 319 "dgetsls.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 323 "dgetsls.f"
	dlaset_("F", &maxmn, nrhs, &c_b23, &c_b23, &b[b_offset], ldb, (ftnlen)
		1);
#line 324 "dgetsls.f"
	goto L50;
#line 325 "dgetsls.f"
    }

#line 327 "dgetsls.f"
    brow = *m;
#line 328 "dgetsls.f"
    if (tran) {
#line 329 "dgetsls.f"
	brow = *n;
#line 330 "dgetsls.f"
    }
#line 331 "dgetsls.f"
    bnrm = dlange_("M", &brow, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 332 "dgetsls.f"
    ibscl = 0;
#line 333 "dgetsls.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 337 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 339 "dgetsls.f"
	ibscl = 1;
#line 340 "dgetsls.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 344 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 346 "dgetsls.f"
	ibscl = 2;
#line 347 "dgetsls.f"
    }

#line 349 "dgetsls.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 353 "dgetsls.f"
	dgeqr_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);
#line 355 "dgetsls.f"
	if (! tran) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 361 "dgetsls.f"
	    dgemqr_("L", "T", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 367 "dgetsls.f"
	    dtrtrs_("U", "N", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 369 "dgetsls.f"
	    if (*info > 0) {
#line 370 "dgetsls.f"
		return 0;
#line 371 "dgetsls.f"
	    }
#line 372 "dgetsls.f"
	    scllen = *n;
#line 373 "dgetsls.f"
	} else {

/*           Overdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */

#line 379 "dgetsls.f"
	    dtrtrs_("U", "T", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 382 "dgetsls.f"
	    if (*info > 0) {
#line 383 "dgetsls.f"
		return 0;
#line 384 "dgetsls.f"
	    }

/*           B(N+1:M,1:NRHS) = ZERO */

#line 388 "dgetsls.f"
	    i__1 = *nrhs;
#line 388 "dgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 389 "dgetsls.f"
		i__2 = *m;
#line 389 "dgetsls.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 390 "dgetsls.f"
		    b[i__ + j * b_dim1] = 0.;
#line 391 "dgetsls.f"
/* L10: */
#line 391 "dgetsls.f"
		}
#line 392 "dgetsls.f"
/* L20: */
#line 392 "dgetsls.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 396 "dgetsls.f"
	    dgemqr_("L", "N", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

#line 400 "dgetsls.f"
	    scllen = *m;

#line 402 "dgetsls.f"
	}

#line 404 "dgetsls.f"
    } else {

/*        Compute LQ factorization of A */

#line 408 "dgetsls.f"
	dgelq_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);

/*        workspace at least M, optimally M*NB. */

#line 413 "dgetsls.f"
	if (! tran) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 419 "dgetsls.f"
	    dtrtrs_("L", "N", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 422 "dgetsls.f"
	    if (*info > 0) {
#line 423 "dgetsls.f"
		return 0;
#line 424 "dgetsls.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 428 "dgetsls.f"
	    i__1 = *nrhs;
#line 428 "dgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 429 "dgetsls.f"
		i__2 = *n;
#line 429 "dgetsls.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 430 "dgetsls.f"
		    b[i__ + j * b_dim1] = 0.;
#line 431 "dgetsls.f"
/* L30: */
#line 431 "dgetsls.f"
		}
#line 432 "dgetsls.f"
/* L40: */
#line 432 "dgetsls.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */

#line 436 "dgetsls.f"
	    dgemlq_("L", "T", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 442 "dgetsls.f"
	    scllen = *n;

#line 444 "dgetsls.f"
	} else {

/*           overdetermined system min || A**T * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 450 "dgetsls.f"
	    dgemlq_("L", "N", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */

#line 458 "dgetsls.f"
	    dtrtrs_("Lower", "Transpose", "Non-unit", m, nrhs, &a[a_offset], 
		    lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);

#line 461 "dgetsls.f"
	    if (*info > 0) {
#line 462 "dgetsls.f"
		return 0;
#line 463 "dgetsls.f"
	    }

#line 465 "dgetsls.f"
	    scllen = *m;

#line 467 "dgetsls.f"
	}

#line 469 "dgetsls.f"
    }

/*     Undo scaling */

#line 473 "dgetsls.f"
    if (iascl == 1) {
#line 474 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 476 "dgetsls.f"
    } else if (iascl == 2) {
#line 477 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 479 "dgetsls.f"
    }
#line 480 "dgetsls.f"
    if (ibscl == 1) {
#line 481 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 483 "dgetsls.f"
    } else if (ibscl == 2) {
#line 484 "dgetsls.f"
	dlascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 486 "dgetsls.f"
    }

#line 488 "dgetsls.f"
L50:
#line 489 "dgetsls.f"
    work[1] = (doublereal) (tszo + lwo);
#line 490 "dgetsls.f"
    return 0;

/*     End of DGETSLS */

} /* dgetsls_ */

