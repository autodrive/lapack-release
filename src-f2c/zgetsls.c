#line 1 "zgetsls.f"
/* zgetsls.f -- translated by f2c (version 20100827).
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

#line 1 "zgetsls.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c_n1 = -1;
static integer c_n2 = -2;
static integer c__0 = 0;

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, */
/*     $                     WORK, LWORK, INFO ) */

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
/* > ZGETSLS solves overdetermined or underdetermined complex linear systems */
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
/* > 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of */
/* >    an undetermined system A**T * X = B. */
/* > */
/* > 4. If TRANS = 'C' and m < n:  find the least squares solution of */
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
/* >          = 'C': the linear system involves A**C. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          A is overwritten by details of its QR or LQ */
/* >          factorization as returned by ZGEQR or ZGELQ. */
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
/* >          squares solution vectors. */
/* >          if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'C' and m >= n, rows 1 to M of B contain the */
/* >          minimum norm solution vectors; */
/* >          if TRANS = 'C' and m < n, rows 1 to M of B contain the */
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
/* >          (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complex16GEsolve */

/*  ===================================================================== */
/* Subroutine */ int zgetsls_(char *trans, integer *m, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublecomplex tq[5];
    static integer lw1, lw2, mnk, lwm, lwo;
    static doublereal anrm, bnrm;
    static logical tran;
    static integer brow, tszm, tszo, info2, iascl, ibscl;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer minmn, maxmn;
    extern /* Subroutine */ int zgelq_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *), zgeqr_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *);
    static doublecomplex workq;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer scllen;
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), zgemlq_(char *, char *, integer *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *, integer *, ftnlen, ftnlen), zlaset_(char *, integer *, integer 
	    *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), zgemqr_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal smlnum;
    static integer wsizem, wsizeo;
    static logical lquery;
    extern /* Subroutine */ int ztrtrs_(char *, char *, char *, integer *, 
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 210 "zgetsls.f"
    /* Parameter adjustments */
#line 210 "zgetsls.f"
    a_dim1 = *lda;
#line 210 "zgetsls.f"
    a_offset = 1 + a_dim1;
#line 210 "zgetsls.f"
    a -= a_offset;
#line 210 "zgetsls.f"
    b_dim1 = *ldb;
#line 210 "zgetsls.f"
    b_offset = 1 + b_dim1;
#line 210 "zgetsls.f"
    b -= b_offset;
#line 210 "zgetsls.f"
    --work;
#line 210 "zgetsls.f"

#line 210 "zgetsls.f"
    /* Function Body */
#line 210 "zgetsls.f"
    *info = 0;
#line 211 "zgetsls.f"
    minmn = min(*m,*n);
#line 212 "zgetsls.f"
    maxmn = max(*m,*n);
#line 213 "zgetsls.f"
    mnk = max(minmn,*nrhs);
#line 214 "zgetsls.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);

#line 216 "zgetsls.f"
    lquery = *lwork == -1 || *lwork == -2;
#line 217 "zgetsls.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1))) {
#line 219 "zgetsls.f"
	*info = -1;
#line 220 "zgetsls.f"
    } else if (*m < 0) {
#line 221 "zgetsls.f"
	*info = -2;
#line 222 "zgetsls.f"
    } else if (*n < 0) {
#line 223 "zgetsls.f"
	*info = -3;
#line 224 "zgetsls.f"
    } else if (*nrhs < 0) {
#line 225 "zgetsls.f"
	*info = -4;
#line 226 "zgetsls.f"
    } else if (*lda < max(1,*m)) {
#line 227 "zgetsls.f"
	*info = -6;
#line 228 "zgetsls.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 228 "zgetsls.f"
	i__1 = max(1,*m);
#line 228 "zgetsls.f"
	if (*ldb < max(i__1,*n)) {
#line 229 "zgetsls.f"
	    *info = -8;
#line 230 "zgetsls.f"
	}
#line 230 "zgetsls.f"
    }

#line 232 "zgetsls.f"
    if (*info == 0) {

/*     Determine the block size and minimum LWORK */

#line 236 "zgetsls.f"
	if (*m >= *n) {
#line 237 "zgetsls.f"
	    zgeqr_(m, n, &a[a_offset], lda, tq, &c_n1, &workq, &c_n1, &info2);
#line 238 "zgetsls.f"
	    tszo = (integer) tq[0].r;
#line 239 "zgetsls.f"
	    lwo = (integer) workq.r;
#line 240 "zgetsls.f"
	    zgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 242 "zgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq.r;
#line 242 "zgetsls.f"
	    lwo = max(i__1,i__2);
#line 243 "zgetsls.f"
	    zgeqr_(m, n, &a[a_offset], lda, tq, &c_n2, &workq, &c_n2, &info2);
#line 244 "zgetsls.f"
	    tszm = (integer) tq[0].r;
#line 245 "zgetsls.f"
	    lwm = (integer) workq.r;
#line 246 "zgetsls.f"
	    zgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszm, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 248 "zgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq.r;
#line 248 "zgetsls.f"
	    lwm = max(i__1,i__2);
#line 249 "zgetsls.f"
	    wsizeo = tszo + lwo;
#line 250 "zgetsls.f"
	    wsizem = tszm + lwm;
#line 251 "zgetsls.f"
	} else {
#line 252 "zgetsls.f"
	    zgelq_(m, n, &a[a_offset], lda, tq, &c_n1, &workq, &c_n1, &info2);
#line 253 "zgetsls.f"
	    tszo = (integer) tq[0].r;
#line 254 "zgetsls.f"
	    lwo = (integer) workq.r;
#line 255 "zgetsls.f"
	    zgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 257 "zgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq.r;
#line 257 "zgetsls.f"
	    lwo = max(i__1,i__2);
#line 258 "zgetsls.f"
	    zgelq_(m, n, &a[a_offset], lda, tq, &c_n2, &workq, &c_n2, &info2);
#line 259 "zgetsls.f"
	    tszm = (integer) tq[0].r;
#line 260 "zgetsls.f"
	    lwm = (integer) workq.r;
#line 261 "zgetsls.f"
	    zgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 263 "zgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq.r;
#line 263 "zgetsls.f"
	    lwm = max(i__1,i__2);
#line 264 "zgetsls.f"
	    wsizeo = tszo + lwo;
#line 265 "zgetsls.f"
	    wsizem = tszm + lwm;
#line 266 "zgetsls.f"
	}

#line 268 "zgetsls.f"
	if (*lwork < wsizem && ! lquery) {
#line 269 "zgetsls.f"
	    *info = -10;
#line 270 "zgetsls.f"
	}

#line 272 "zgetsls.f"
    }

#line 274 "zgetsls.f"
    if (*info != 0) {
#line 275 "zgetsls.f"
	i__1 = -(*info);
#line 275 "zgetsls.f"
	xerbla_("ZGETSLS", &i__1, (ftnlen)7);
#line 276 "zgetsls.f"
	d__1 = (doublereal) wsizeo;
#line 276 "zgetsls.f"
	work[1].r = d__1, work[1].i = 0.;
#line 277 "zgetsls.f"
	return 0;
#line 278 "zgetsls.f"
    }
#line 279 "zgetsls.f"
    if (lquery) {
#line 280 "zgetsls.f"
	if (*lwork == -1) {
#line 280 "zgetsls.f"
	    d__1 = (doublereal) wsizeo;
#line 280 "zgetsls.f"
	    work[1].r = d__1, work[1].i = 0.;
#line 280 "zgetsls.f"
	}
#line 281 "zgetsls.f"
	if (*lwork == -2) {
#line 281 "zgetsls.f"
	    d__1 = (doublereal) wsizem;
#line 281 "zgetsls.f"
	    work[1].r = d__1, work[1].i = 0.;
#line 281 "zgetsls.f"
	}
#line 282 "zgetsls.f"
	return 0;
#line 283 "zgetsls.f"
    }
#line 284 "zgetsls.f"
    if (*lwork < wsizeo) {
#line 285 "zgetsls.f"
	lw1 = tszm;
#line 286 "zgetsls.f"
	lw2 = lwm;
#line 287 "zgetsls.f"
    } else {
#line 288 "zgetsls.f"
	lw1 = tszo;
#line 289 "zgetsls.f"
	lw2 = lwo;
#line 290 "zgetsls.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 294 "zgetsls.f"
    i__1 = min(*m,*n);
#line 294 "zgetsls.f"
    if (min(i__1,*nrhs) == 0) {
#line 295 "zgetsls.f"
	i__1 = max(*m,*n);
#line 295 "zgetsls.f"
	zlaset_("FULL", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)
		4);
#line 297 "zgetsls.f"
	return 0;
#line 298 "zgetsls.f"
    }

/*     Get machine parameters */

#line 302 "zgetsls.f"
    smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 303 "zgetsls.f"
    bignum = 1. / smlnum;
#line 304 "zgetsls.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 308 "zgetsls.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 309 "zgetsls.f"
    iascl = 0;
#line 310 "zgetsls.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 314 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 315 "zgetsls.f"
	iascl = 1;
#line 316 "zgetsls.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 320 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 321 "zgetsls.f"
	iascl = 2;
#line 322 "zgetsls.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 326 "zgetsls.f"
	zlaset_("F", &maxmn, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1)
		;
#line 327 "zgetsls.f"
	goto L50;
#line 328 "zgetsls.f"
    }

#line 330 "zgetsls.f"
    brow = *m;
#line 331 "zgetsls.f"
    if (tran) {
#line 332 "zgetsls.f"
	brow = *n;
#line 333 "zgetsls.f"
    }
#line 334 "zgetsls.f"
    bnrm = zlange_("M", &brow, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 335 "zgetsls.f"
    ibscl = 0;
#line 336 "zgetsls.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 340 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 342 "zgetsls.f"
	ibscl = 1;
#line 343 "zgetsls.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 347 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 349 "zgetsls.f"
	ibscl = 2;
#line 350 "zgetsls.f"
    }

#line 352 "zgetsls.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 356 "zgetsls.f"
	zgeqr_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);
#line 358 "zgetsls.f"
	if (! tran) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 364 "zgetsls.f"
	    zgemqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 370 "zgetsls.f"
	    ztrtrs_("U", "N", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 372 "zgetsls.f"
	    if (*info > 0) {
#line 373 "zgetsls.f"
		return 0;
#line 374 "zgetsls.f"
	    }
#line 375 "zgetsls.f"
	    scllen = *n;
#line 376 "zgetsls.f"
	} else {

/*           Overdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */

#line 382 "zgetsls.f"
	    ztrtrs_("U", "C", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 385 "zgetsls.f"
	    if (*info > 0) {
#line 386 "zgetsls.f"
		return 0;
#line 387 "zgetsls.f"
	    }

/*           B(N+1:M,1:NRHS) = CZERO */

#line 391 "zgetsls.f"
	    i__1 = *nrhs;
#line 391 "zgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 392 "zgetsls.f"
		i__2 = *m;
#line 392 "zgetsls.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 393 "zgetsls.f"
		    i__3 = i__ + j * b_dim1;
#line 393 "zgetsls.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 394 "zgetsls.f"
/* L10: */
#line 394 "zgetsls.f"
		}
#line 395 "zgetsls.f"
/* L20: */
#line 395 "zgetsls.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 399 "zgetsls.f"
	    zgemqr_("L", "N", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

#line 403 "zgetsls.f"
	    scllen = *m;

#line 405 "zgetsls.f"
	}

#line 407 "zgetsls.f"
    } else {

/*        Compute LQ factorization of A */

#line 411 "zgetsls.f"
	zgelq_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);

/*        workspace at least M, optimally M*NB. */

#line 416 "zgetsls.f"
	if (! tran) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 422 "zgetsls.f"
	    ztrtrs_("L", "N", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 425 "zgetsls.f"
	    if (*info > 0) {
#line 426 "zgetsls.f"
		return 0;
#line 427 "zgetsls.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 431 "zgetsls.f"
	    i__1 = *nrhs;
#line 431 "zgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 432 "zgetsls.f"
		i__2 = *n;
#line 432 "zgetsls.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 433 "zgetsls.f"
		    i__3 = i__ + j * b_dim1;
#line 433 "zgetsls.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 434 "zgetsls.f"
/* L30: */
#line 434 "zgetsls.f"
		}
#line 435 "zgetsls.f"
/* L40: */
#line 435 "zgetsls.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */

#line 439 "zgetsls.f"
	    zgemlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 445 "zgetsls.f"
	    scllen = *n;

#line 447 "zgetsls.f"
	} else {

/*           overdetermined system min || A**T * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 453 "zgetsls.f"
	    zgemlq_("L", "N", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */

#line 461 "zgetsls.f"
	    ztrtrs_("L", "C", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 464 "zgetsls.f"
	    if (*info > 0) {
#line 465 "zgetsls.f"
		return 0;
#line 466 "zgetsls.f"
	    }

#line 468 "zgetsls.f"
	    scllen = *m;

#line 470 "zgetsls.f"
	}

#line 472 "zgetsls.f"
    }

/*     Undo scaling */

#line 476 "zgetsls.f"
    if (iascl == 1) {
#line 477 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 479 "zgetsls.f"
    } else if (iascl == 2) {
#line 480 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 482 "zgetsls.f"
    }
#line 483 "zgetsls.f"
    if (ibscl == 1) {
#line 484 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 486 "zgetsls.f"
    } else if (ibscl == 2) {
#line 487 "zgetsls.f"
	zlascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 489 "zgetsls.f"
    }

#line 491 "zgetsls.f"
L50:
#line 492 "zgetsls.f"
    d__1 = (doublereal) (tszo + lwo);
#line 492 "zgetsls.f"
    work[1].r = d__1, work[1].i = 0.;
#line 493 "zgetsls.f"
    return 0;

/*     End of ZGETSLS */

} /* zgetsls_ */

