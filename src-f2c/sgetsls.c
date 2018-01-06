#line 1 "sgetsls.f"
/* sgetsls.f -- translated by f2c (version 20100827).
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

#line 1 "sgetsls.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c_n2 = -2;
static doublereal c_b23 = 0.;
static integer c__0 = 0;

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, */
/*     $                     WORK, LWORK, INFO ) */

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
/* > SGETSLS solves overdetermined or underdetermined real linear systems */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          A is overwritten by details of its QR or LQ */
/* >          factorization as returned by SGEQR or SGELQ. */
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
/* >          (workspace) REAL array, dimension (MAX(1,LWORK)) */
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

/* > \date June 2017 */

/* > \ingroup doubleGEsolve */

/*  ===================================================================== */
/* Subroutine */ int sgetsls_(char *trans, integer *m, integer *n, integer *
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgelq_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static integer minmn, maxmn;
    extern /* Subroutine */ int sgeqr_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static doublereal workq[1];
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer scllen;
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), sgemlq_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    sgemqr_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal smlnum;
    static integer wsizem, wsizeo;
    static logical lquery;
    extern /* Subroutine */ int strtrs_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 207 "sgetsls.f"
    /* Parameter adjustments */
#line 207 "sgetsls.f"
    a_dim1 = *lda;
#line 207 "sgetsls.f"
    a_offset = 1 + a_dim1;
#line 207 "sgetsls.f"
    a -= a_offset;
#line 207 "sgetsls.f"
    b_dim1 = *ldb;
#line 207 "sgetsls.f"
    b_offset = 1 + b_dim1;
#line 207 "sgetsls.f"
    b -= b_offset;
#line 207 "sgetsls.f"
    --work;
#line 207 "sgetsls.f"

#line 207 "sgetsls.f"
    /* Function Body */
#line 207 "sgetsls.f"
    *info = 0;
#line 208 "sgetsls.f"
    minmn = min(*m,*n);
#line 209 "sgetsls.f"
    maxmn = max(*m,*n);
#line 210 "sgetsls.f"
    mnk = max(minmn,*nrhs);
#line 211 "sgetsls.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

#line 213 "sgetsls.f"
    lquery = *lwork == -1 || *lwork == -2;
#line 214 "sgetsls.f"
    if (! (lsame_(trans, "N", (ftnlen)1, (ftnlen)1) || lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1))) {
#line 216 "sgetsls.f"
	*info = -1;
#line 217 "sgetsls.f"
    } else if (*m < 0) {
#line 218 "sgetsls.f"
	*info = -2;
#line 219 "sgetsls.f"
    } else if (*n < 0) {
#line 220 "sgetsls.f"
	*info = -3;
#line 221 "sgetsls.f"
    } else if (*nrhs < 0) {
#line 222 "sgetsls.f"
	*info = -4;
#line 223 "sgetsls.f"
    } else if (*lda < max(1,*m)) {
#line 224 "sgetsls.f"
	*info = -6;
#line 225 "sgetsls.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 225 "sgetsls.f"
	i__1 = max(1,*m);
#line 225 "sgetsls.f"
	if (*ldb < max(i__1,*n)) {
#line 226 "sgetsls.f"
	    *info = -8;
#line 227 "sgetsls.f"
	}
#line 227 "sgetsls.f"
    }

#line 229 "sgetsls.f"
    if (*info == 0) {

/*     Determine the block size and minimum LWORK */

#line 233 "sgetsls.f"
	if (*m >= *n) {
#line 234 "sgetsls.f"
	    sgeqr_(m, n, &a[a_offset], lda, tq, &c_n1, workq, &c_n1, &info2);
#line 235 "sgetsls.f"
	    tszo = (integer) tq[0];
#line 236 "sgetsls.f"
	    lwo = (integer) workq[0];
#line 237 "sgetsls.f"
	    sgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 239 "sgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq[0];
#line 239 "sgetsls.f"
	    lwo = max(i__1,i__2);
#line 240 "sgetsls.f"
	    sgeqr_(m, n, &a[a_offset], lda, tq, &c_n2, workq, &c_n2, &info2);
#line 241 "sgetsls.f"
	    tszm = (integer) tq[0];
#line 242 "sgetsls.f"
	    lwm = (integer) workq[0];
#line 243 "sgetsls.f"
	    sgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszm, &b[
		    b_offset], ldb, workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 245 "sgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq[0];
#line 245 "sgetsls.f"
	    lwm = max(i__1,i__2);
#line 246 "sgetsls.f"
	    wsizeo = tszo + lwo;
#line 247 "sgetsls.f"
	    wsizem = tszm + lwm;
#line 248 "sgetsls.f"
	} else {
#line 249 "sgetsls.f"
	    sgelq_(m, n, &a[a_offset], lda, tq, &c_n1, workq, &c_n1, &info2);
#line 250 "sgetsls.f"
	    tszo = (integer) tq[0];
#line 251 "sgetsls.f"
	    lwo = (integer) workq[0];
#line 252 "sgetsls.f"
	    sgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 254 "sgetsls.f"
	    i__1 = lwo, i__2 = (integer) workq[0];
#line 254 "sgetsls.f"
	    lwo = max(i__1,i__2);
#line 255 "sgetsls.f"
	    sgelq_(m, n, &a[a_offset], lda, tq, &c_n2, workq, &c_n2, &info2);
#line 256 "sgetsls.f"
	    tszm = (integer) tq[0];
#line 257 "sgetsls.f"
	    lwm = (integer) workq[0];
#line 258 "sgetsls.f"
	    sgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[
		    b_offset], ldb, workq, &c_n1, &info2, (ftnlen)1, (ftnlen)
		    1);
/* Computing MAX */
#line 260 "sgetsls.f"
	    i__1 = lwm, i__2 = (integer) workq[0];
#line 260 "sgetsls.f"
	    lwm = max(i__1,i__2);
#line 261 "sgetsls.f"
	    wsizeo = tszo + lwo;
#line 262 "sgetsls.f"
	    wsizem = tszm + lwm;
#line 263 "sgetsls.f"
	}

#line 265 "sgetsls.f"
	if (*lwork < wsizem && ! lquery) {
#line 266 "sgetsls.f"
	    *info = -10;
#line 267 "sgetsls.f"
	}

#line 269 "sgetsls.f"
    }

#line 271 "sgetsls.f"
    if (*info != 0) {
#line 272 "sgetsls.f"
	i__1 = -(*info);
#line 272 "sgetsls.f"
	xerbla_("SGETSLS", &i__1, (ftnlen)7);
#line 273 "sgetsls.f"
	work[1] = (doublereal) wsizeo;
#line 274 "sgetsls.f"
	return 0;
#line 275 "sgetsls.f"
    }
#line 276 "sgetsls.f"
    if (lquery) {
#line 277 "sgetsls.f"
	if (*lwork == -1) {
#line 277 "sgetsls.f"
	    work[1] = (doublereal) wsizeo;
#line 277 "sgetsls.f"
	}
#line 278 "sgetsls.f"
	if (*lwork == -2) {
#line 278 "sgetsls.f"
	    work[1] = (doublereal) wsizem;
#line 278 "sgetsls.f"
	}
#line 279 "sgetsls.f"
	return 0;
#line 280 "sgetsls.f"
    }
#line 281 "sgetsls.f"
    if (*lwork < wsizeo) {
#line 282 "sgetsls.f"
	lw1 = tszm;
#line 283 "sgetsls.f"
	lw2 = lwm;
#line 284 "sgetsls.f"
    } else {
#line 285 "sgetsls.f"
	lw1 = tszo;
#line 286 "sgetsls.f"
	lw2 = lwo;
#line 287 "sgetsls.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 291 "sgetsls.f"
    i__1 = min(*m,*n);
#line 291 "sgetsls.f"
    if (min(i__1,*nrhs) == 0) {
#line 292 "sgetsls.f"
	i__1 = max(*m,*n);
#line 292 "sgetsls.f"
	slaset_("FULL", &i__1, nrhs, &c_b23, &c_b23, &b[b_offset], ldb, (
		ftnlen)4);
#line 294 "sgetsls.f"
	return 0;
#line 295 "sgetsls.f"
    }

/*     Get machine parameters */

#line 299 "sgetsls.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 300 "sgetsls.f"
    bignum = 1. / smlnum;
#line 301 "sgetsls.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

#line 305 "sgetsls.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 306 "sgetsls.f"
    iascl = 0;
#line 307 "sgetsls.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 311 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 312 "sgetsls.f"
	iascl = 1;
#line 313 "sgetsls.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 317 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 318 "sgetsls.f"
	iascl = 2;
#line 319 "sgetsls.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 323 "sgetsls.f"
	slaset_("F", &maxmn, nrhs, &c_b23, &c_b23, &b[b_offset], ldb, (ftnlen)
		1);
#line 324 "sgetsls.f"
	goto L50;
#line 325 "sgetsls.f"
    }

#line 327 "sgetsls.f"
    brow = *m;
#line 328 "sgetsls.f"
    if (tran) {
#line 329 "sgetsls.f"
	brow = *n;
#line 330 "sgetsls.f"
    }
#line 331 "sgetsls.f"
    bnrm = slange_("M", &brow, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 332 "sgetsls.f"
    ibscl = 0;
#line 333 "sgetsls.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 337 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 339 "sgetsls.f"
	ibscl = 1;
#line 340 "sgetsls.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 344 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);
#line 346 "sgetsls.f"
	ibscl = 2;
#line 347 "sgetsls.f"
    }

#line 349 "sgetsls.f"
    if (*m >= *n) {

/*        compute QR factorization of A */

#line 353 "sgetsls.f"
	sgeqr_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);
#line 355 "sgetsls.f"
	if (! tran) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 361 "sgetsls.f"
	    sgemqr_("L", "T", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

#line 367 "sgetsls.f"
	    strtrs_("U", "N", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 369 "sgetsls.f"
	    if (*info > 0) {
#line 370 "sgetsls.f"
		return 0;
#line 371 "sgetsls.f"
	    }
#line 372 "sgetsls.f"
	    scllen = *n;
#line 373 "sgetsls.f"
	} else {

/*           Overdetermined system of equations A**T * X = B */

/*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */

#line 379 "sgetsls.f"
	    strtrs_("U", "T", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 382 "sgetsls.f"
	    if (*info > 0) {
#line 383 "sgetsls.f"
		return 0;
#line 384 "sgetsls.f"
	    }

/*           B(N+1:M,1:NRHS) = ZERO */

#line 388 "sgetsls.f"
	    i__1 = *nrhs;
#line 388 "sgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 389 "sgetsls.f"
		i__2 = *m;
#line 389 "sgetsls.f"
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
#line 390 "sgetsls.f"
		    b[i__ + j * b_dim1] = 0.;
#line 391 "sgetsls.f"
/* L10: */
#line 391 "sgetsls.f"
		}
#line 392 "sgetsls.f"
/* L20: */
#line 392 "sgetsls.f"
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

#line 396 "sgetsls.f"
	    sgemqr_("L", "N", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

#line 400 "sgetsls.f"
	    scllen = *m;

#line 402 "sgetsls.f"
	}

#line 404 "sgetsls.f"
    } else {

/*        Compute LQ factorization of A */

#line 408 "sgetsls.f"
	sgelq_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, 
		info);

/*        workspace at least M, optimally M*NB. */

#line 413 "sgetsls.f"
	if (! tran) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

#line 419 "sgetsls.f"
	    strtrs_("L", "N", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], 
		    ldb, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 422 "sgetsls.f"
	    if (*info > 0) {
#line 423 "sgetsls.f"
		return 0;
#line 424 "sgetsls.f"
	    }

/*           B(M+1:N,1:NRHS) = 0 */

#line 428 "sgetsls.f"
	    i__1 = *nrhs;
#line 428 "sgetsls.f"
	    for (j = 1; j <= i__1; ++j) {
#line 429 "sgetsls.f"
		i__2 = *n;
#line 429 "sgetsls.f"
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 430 "sgetsls.f"
		    b[i__ + j * b_dim1] = 0.;
#line 431 "sgetsls.f"
/* L30: */
#line 431 "sgetsls.f"
		}
#line 432 "sgetsls.f"
/* L40: */
#line 432 "sgetsls.f"
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */

#line 436 "sgetsls.f"
	    sgemlq_("L", "T", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

#line 442 "sgetsls.f"
	    scllen = *n;

#line 444 "sgetsls.f"
	} else {

/*           overdetermined system min || A**T * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

#line 450 "sgetsls.f"
	    sgemlq_("L", "N", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &
		    lw1, &b[b_offset], ldb, &work[1], &lw2, info, (ftnlen)1, (
		    ftnlen)1);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */

#line 458 "sgetsls.f"
	    strtrs_("Lower", "Transpose", "Non-unit", m, nrhs, &a[a_offset], 
		    lda, &b[b_offset], ldb, info, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);

#line 461 "sgetsls.f"
	    if (*info > 0) {
#line 462 "sgetsls.f"
		return 0;
#line 463 "sgetsls.f"
	    }

#line 465 "sgetsls.f"
	    scllen = *m;

#line 467 "sgetsls.f"
	}

#line 469 "sgetsls.f"
    }

/*     Undo scaling */

#line 473 "sgetsls.f"
    if (iascl == 1) {
#line 474 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 476 "sgetsls.f"
    } else if (iascl == 2) {
#line 477 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 479 "sgetsls.f"
    }
#line 480 "sgetsls.f"
    if (ibscl == 1) {
#line 481 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 483 "sgetsls.f"
    } else if (ibscl == 2) {
#line 484 "sgetsls.f"
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
		, ldb, info, (ftnlen)1);
#line 486 "sgetsls.f"
    }

#line 488 "sgetsls.f"
L50:
#line 489 "sgetsls.f"
    work[1] = (doublereal) (tszo + lwo);
#line 490 "sgetsls.f"
    return 0;

/*     End of SGETSLS */

} /* sgetsls_ */

