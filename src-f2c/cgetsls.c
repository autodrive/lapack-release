#line 1 "cgetsls.f"
/* cgetsls.f -- translated by f2c (version 20100827).
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

#line 1 "cgetsls.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c_n1 = -1;
static integer c_n2 = -2;
static integer c__0 = 0;

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, */
/*     $                     WORK, LWORK, INFO ) */

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
/* > CGETSLS solves overdetermined or underdetermined complex linear systems */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          A is overwritten by details of its QR or LQ */
/* >          factorization as returned by CGEQR or CGELQ. */
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
/* >          (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complexGEsolve */

/*  ===================================================================== */
/* Subroutine */ int cgetsls_(char *trans, integer *m, integer *n, integer *
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
    extern /* Subroutine */ int cgelq_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgeqr_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *);
    static integer minmn, maxmn;
    static doublecomplex workq;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int cgemlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), cgemqr_(char *, char 
	    *, integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    static integer scllen;
    static doublereal bignum, smlnum;
    static integer wsizem, wsizeo;
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 210 "cgetsls.f"
    /* Parameter adjustments */
#line 210 "cgetsls.f"
    a_dim1 = *lda;
#line 210 "cgetsls.f"
    a_offset = 1 + a_dim1;
#line 210 "cgetsls.f"
    a -= a_offset;
#line 210 "cgetsls.f"
    b_dim1 = *ldb;
#line 210 "cgetsls.f"
    b_offset = 1 + b_dim1;
#line 210 "cgetsls.f"
    b -= b_offset;
#line 210 "cgetsls.f"
    --work;
#line 210 "cgetsls.f"

#line 210 "cgetsls.f"
    /* Function Body */
#line 210 "cgetsls.f"
    *info = 0;
#line 211 "cgetsls.f"
    minmn = min(*m,*n);
#line 212 "cgetsls.f"
    maxmn = max(*m,*n);
#line 213 "cgetsls.f"
    mnk = max(minmn,*nrhs);
#line 214 "cgetsls.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);

#line 216 "cgetsls.f"
    lquery = *lwork == -1 || *lwork == -2;
#line 217 "cgetsls.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1))) {
#line 219 "cgetsls.f"
	*info = -1;
#line 220 "cgetsls.f"
    } else if (*m < 0) {
#line 221 "cgetsls.f"
	*info = -2;
#line 222 "cgetsls.f"
    } else if (*n < 0) {
#line 223 "cgetsls.f"
	*info = -3;
#line 224 "cgetsls.f"
    } else if (*nrhs < 0) {
#line 225 "cgetsls.f"
	*info = -4;
#line 226 "cgetsls.f"
    } else if (*lda < max(1,*m)) {
#line 227 "cgetsls.f"
	*info = -6;
#line 228 "cgetsls.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 228 "cgetsls.f"
	i__1 = max(1,*m);
#line 228 "cgetsls.f"
	if (*ldb < max(i__1,*n)) {
#line 229 "cgetsls.f"
	    *info = -8;
#line 230 "cgetsls.f"
	}
#line 230 "cgetsls.f"
    }

#line 232 "cgetsls.f"
    if (*info == 0) {

/*     Determine the block size and minimum LWORK */

#line 236 "cgetsls.f"
	if (*m >= *n) {
#line 237 "cgetsls.f"
	    cgeqr_(m, n, &a[a_offset], lda, tq, &c_n1, &workq, &c_n1, &info2);
#line 238 "cgetsls.f"
	    tszo = (integer) tq[0].r;
#line 239 "cgetsls.f"
	    lwo = (integer) workq.r;
#line 240 "cgetsls.f"
	    cgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 242 "cgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq.r;
#line 242 "cgetsls.f"
	    lwo = max(i__1,i__2);
#line 243 "cgetsls.f"
	    cgeqr_(m, n, &a[a_offset], lda, tq, &c_n2, &workq, &c_n2, &info2);
#line 244 "cgetsls.f"
	    tszm = (integer) tq[0].r;
#line 245 "cgetsls.f"
	    lwm = (integer) workq.r;
#line 246 "cgetsls.f"
	    cgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszm, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 248 "cgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq.r;
#line 248 "cgetsls.f"
	    lwm = max(i__1,i__2);
#line 249 "cgetsls.f"
	    wsizeo = tszo + lwo;
#line 250 "cgetsls.f"
	    wsizem = tszm + lwm;
#line 251 "cgetsls.f"
	} else {
#line 252 "cgetsls.f"
	    cgelq_(m, n, &a[a_offset], lda, tq, &c_n1, &workq, &c_n1, &info2);
#line 253 "cgetsls.f"
	    tszo = (integer) tq[0].r;
#line 254 "cgetsls.f"
	    lwo = (integer) workq.r;
#line 255 "cgetsls.f"
	    cgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 257 "cgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq.r;
#line 257 "cgetsls.f"
	    lwo = max(i__1,i__2);
#line 258 "cgetsls.f"
	    cgelq_(m, n, &a[a_offset], lda, tq, &c_n2, &workq, &c_n2, &info2);
#line 259 "cgetsls.f"
	    tszm = (integer) tq[0].r;
#line 260 "cgetsls.f"
	    lwm = (integer) workq.r;
#line 261 "cgetsls.f"
	    cgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, &workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 263 "cgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq.r;
#line 263 "cgetsls.f"
	    lwm = max(i__1,i__2);
#line 264 "cgetsls.f"
	    wsizeo = tszo + lwo;
#line 265 "cgetsls.f"
	    wsizem = tszm + lwm;
#line 266 "cgetsls.f"
	}

#line 268 "cgetsls.f"
	if (*lwork < wsizem && ! lquery) {
#line 269 "cgetsls.f"
	    *info = -10;
#line 270 "cgetsls.f"
	}

#line 272 "cgetsls.f"
    }

#line 274 "cgetsls.f"
    if (*info != 0) {
#line 275 "cgetsls.f"
	i__1 = -(*info);
#line 275 "cgetsls.f"
	xerbla_("CGETSLS", &i__1, (ftnlen)7);
#line 276 "cgetsls.f"
	d__1 = (doublereal) wsizeo;
#line 276 "cgetsls.f"
	work[1].r = d__1, work[1].i = 0.;
#line 277 "cgetsls.f"
	return 0;
#line 278 "cgetsls.f"
    }
#line 279 "cgetsls.f"
    if (lquery) {
#line 280 "cgetsls.f"
	if (*lwork == -1) {
#line 280 "cgetsls.f"
	    d__1 = (doublereal) wsizeo;
#line 280 "cgetsls.f"
	    work[1].r = d__1, work[1].i = 0.;
#line 280 "cgetsls.f"
	}
#line 281 "cgetsls.f"
	if (*lwork == -2) {
#line 281 "cgetsls.f"
	    d__1 = (doublereal) wsizem;
#line 281 "cgetsls.f"
	    work[1].r = d__1, work[1].i = 0.;
#line 281 "cgetsls.f"
	}
#line 282 "cgetsls.f"
	return 0;
#line 283 "cgetsls.f"
    }
#line 284 "cgetsls.f"
    if (*lwork < wsizeo) {
#line 285 "cgetsls.f"
	lw1 = tszm;
#line 286 "cgetsls.f"
	lw2 = lwm;
#line 287 "cgetsls.f"
    } else {
#line 288 "cgetsls.f"
	lw1 = tszo;
#line 289 "cgetsls.f"
	lw2 = lwo;
#line 290 "cgetsls.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 294 "cgetsls.f"
    i__1 = min(*m,*n);
#line 294 "cgetsls.f"
    if (min(i__1,*nrhs) == 0) {
#line 295 "cgetsls.f"
	i__1 = max(*m,*n);
#line 295 "cgetsls.f"
	claset_("FULL", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)
		4);
#line 297 "cgetsls.f"
	return 0;
#line 298 "cgetsls.f"
    }

/*     Get machine parameters */

#line 302 "cgetsls.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 303 "cgetsls.f"
    bignum = 1. / smlnum;
#line 304 "cgetsls.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 308 "cgetsls.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 309 "cgetsls.f"
    iascl = 0;
#line 310 "cgetsls.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 314 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 315 "cgetsls.f"
	iascl = 1;
#line 316 "cgetsls.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 320 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 321 "cgetsls.f"
	iascl = 2;
#line 322 "cgetsls.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 326 "cgetsls.f"
	claset_("F", &maxmn, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1)
		;
#line 327 "cgetsls.f"
	goto L50;
#line 328 "cgetsls.f"
    }

#line 330 "cgetsls.f"
    brow = *m;
#line 331 "cgetsls.f"
    if (tran) {
#line 332 "cgetsls.f"
	brow = *n;
#line 333 "cgetsls.f"
    }
#line 334 "cgetsls.f"
    bnrm = clange_("M", &brow, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 335 "cgetsls.f"
    ibscl = 0;
#line 336 "cgetsls.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 340 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 342 "cgetsls.f"
	ibscl = 1;
#line 343 "cgetsls.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 347 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 349 "cgetsls.f"
	ibscl = 2;
#line 350 "cgetsls.f"
    }

#line 352 "cgetsls.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 356 "cgetsls.f"
	cgeqr_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);
#line 358 "cgetsls.f"
	if (! tran) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 364 "cgetsls.f"
	    cgemqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 370 "cgetsls.f"
	    ctrtrs_("U", "N", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 372 "cgetsls.f"
	    if (*info > 0) {
#line 373 "cgetsls.f"
		return 0;
#line 374 "cgetsls.f"
	    }
#line 375 "cgetsls.f"
	    scllen = *n;
#line 376 "cgetsls.f"
	} else {

/*           Overdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */

#line 382 "cgetsls.f"
	    ctrtrs_("U", "C", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 385 "cgetsls.f"
	    if (*info > 0) {
#line 386 "cgetsls.f"
		return 0;
#line 387 "cgetsls.f"
	    }

/*           B(N+1:M,1:NRHS) = CZERO */

#line 391 "cgetsls.f"
	    i__1 = *nrhs;
#line 391 "cgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 392 "cgetsls.f"
		i__2 = *m;
#line 392 "cgetsls.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 393 "cgetsls.f"
		    i__3 = i__ + j * b_dim1;
#line 393 "cgetsls.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 394 "cgetsls.f"
/* L10: */
#line 394 "cgetsls.f"
		}
#line 395 "cgetsls.f"
/* L20: */
#line 395 "cgetsls.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 399 "cgetsls.f"
	    cgemqr_("L", "N", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

#line 403 "cgetsls.f"
	    scllen = *m;

#line 405 "cgetsls.f"
	}

#line 407 "cgetsls.f"
    } else {

/*        Compute LQ factorization of A */

#line 411 "cgetsls.f"
	cgelq_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);

/*        workspace at least M, optimally M*NB. */

#line 416 "cgetsls.f"
	if (! tran) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 422 "cgetsls.f"
	    ctrtrs_("L", "N", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 425 "cgetsls.f"
	    if (*info > 0) {
#line 426 "cgetsls.f"
		return 0;
#line 427 "cgetsls.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 431 "cgetsls.f"
	    i__1 = *nrhs;
#line 431 "cgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 432 "cgetsls.f"
		i__2 = *n;
#line 432 "cgetsls.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 433 "cgetsls.f"
		    i__3 = i__ + j * b_dim1;
#line 433 "cgetsls.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 434 "cgetsls.f"
/* L30: */
#line 434 "cgetsls.f"
		}
#line 435 "cgetsls.f"
/* L40: */
#line 435 "cgetsls.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */

#line 439 "cgetsls.f"
	    cgemlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 445 "cgetsls.f"
	    scllen = *n;

#line 447 "cgetsls.f"
	} else {

/*           overdetermined system min || A**T * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 453 "cgetsls.f"
	    cgemlq_("L", "N", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */

#line 461 "cgetsls.f"
	    ctrtrs_("L", "C", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 464 "cgetsls.f"
	    if (*info > 0) {
#line 465 "cgetsls.f"
		return 0;
#line 466 "cgetsls.f"
	    }

#line 468 "cgetsls.f"
	    scllen = *m;

#line 470 "cgetsls.f"
	}

#line 472 "cgetsls.f"
    }

/*     Undo scaling */

#line 476 "cgetsls.f"
    if (iascl == 1) {
#line 477 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 479 "cgetsls.f"
    } else if (iascl == 2) {
#line 480 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 482 "cgetsls.f"
    }
#line 483 "cgetsls.f"
    if (ibscl == 1) {
#line 484 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 486 "cgetsls.f"
    } else if (ibscl == 2) {
#line 487 "cgetsls.f"
	clascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 489 "cgetsls.f"
    }

#line 491 "cgetsls.f"
L50:
#line 492 "cgetsls.f"
    d__1 = (doublereal) (tszo + lwo);
#line 492 "cgetsls.f"
    work[1].r = d__1, work[1].i = 0.;
#line 493 "cgetsls.f"
    return 0;

/*     End of ZGETSLS */

} /* cgetsls_ */

