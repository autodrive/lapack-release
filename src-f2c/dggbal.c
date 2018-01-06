#line 1 "dggbal.f"
/* dggbal.f -- translated by f2c (version 20100827).
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

#line 1 "dggbal.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b35 = 10.;
static doublereal c_b71 = .5;

/* > \brief \b DGGBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, */
/*                          RSCALE, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ), */
/*      $                   RSCALE( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGBAL balances a pair of general real matrices (A,B).  This */
/* > involves, first, permuting A and B by similarity transformations to */
/* > isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N */
/* > elements on the diagonal; and second, applying a diagonal similarity */
/* > transformation to rows and columns ILO to IHI to make the rows */
/* > and columns as close in norm as possible. Both steps are optional. */
/* > */
/* > Balancing may reduce the 1-norm of the matrices, and improve the */
/* > accuracy of the computed eigenvalues and/or eigenvectors in the */
/* > generalized eigenvalue problem A*x = lambda*B*x. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies the operations to be performed on A and B: */
/* >          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0 */
/* >                  and RSCALE(I) = 1.0 for i = 1,...,N. */
/* >          = 'P':  permute only; */
/* >          = 'S':  scale only; */
/* >          = 'B':  both permute and scale. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the input matrix A. */
/* >          On exit,  A is overwritten by the balanced matrix. */
/* >          If JOB = 'N', A is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the input matrix B. */
/* >          On exit,  B is overwritten by the balanced matrix. */
/* >          If JOB = 'N', B is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          ILO and IHI are set to integers such that on exit */
/* >          A(i,j) = 0 and B(i,j) = 0 if i > j and */
/* >          j = 1,...,ILO-1 or i = IHI+1,...,N. */
/* >          If JOB = 'N' or 'S', ILO = 1 and IHI = N. */
/* > \endverbatim */
/* > */
/* > \param[out] LSCALE */
/* > \verbatim */
/* >          LSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and scaling factors applied */
/* >          to the left side of A and B.  If P(j) is the index of the */
/* >          row interchanged with row j, and D(j) */
/* >          is the scaling factor applied to row j, then */
/* >            LSCALE(j) = P(j)    for J = 1,...,ILO-1 */
/* >                      = D(j)    for J = ILO,...,IHI */
/* >                      = P(j)    for J = IHI+1,...,N. */
/* >          The order in which the interchanges are made is N to IHI+1, */
/* >          then 1 to ILO-1. */
/* > \endverbatim */
/* > */
/* > \param[out] RSCALE */
/* > \verbatim */
/* >          RSCALE is DOUBLE PRECISION array, dimension (N) */
/* >          Details of the permutations and scaling factors applied */
/* >          to the right side of A and B.  If P(j) is the index of the */
/* >          column interchanged with column j, and D(j) */
/* >          is the scaling factor applied to column j, then */
/* >            LSCALE(j) = P(j)    for J = 1,...,ILO-1 */
/* >                      = D(j)    for J = ILO,...,IHI */
/* >                      = P(j)    for J = IHI+1,...,N. */
/* >          The order in which the interchanges are made is N to IHI+1, */
/* >          then 1 to ILO-1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (lwork) */
/* >          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and */
/* >          at least 1 when JOB = 'N' or 'P'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleGBcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  See R.C. WARD, Balancing the generalized eigenvalue problem, */
/* >                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dggbal_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_lg10(doublereal *), d_sign(doublereal *, doublereal *), pow_di(
	    doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal t;
    static integer jc;
    static doublereal ta, tb, tc;
    static integer ir;
    static doublereal ew;
    static integer it, nr, ip1, jp1, lm1;
    static doublereal cab, rab, ewc, cor, sum;
    static integer nrp2, icab, lcab;
    static doublereal beta, coef;
    static integer irab, lrab;
    static doublereal basl, cmax;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal coef2, coef5, gamma, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal sfmin, sfmax;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer iflow;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer kount;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal pgamma;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lsfmin, lsfmax;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

/*     Test the input parameters */

#line 226 "dggbal.f"
    /* Parameter adjustments */
#line 226 "dggbal.f"
    a_dim1 = *lda;
#line 226 "dggbal.f"
    a_offset = 1 + a_dim1;
#line 226 "dggbal.f"
    a -= a_offset;
#line 226 "dggbal.f"
    b_dim1 = *ldb;
#line 226 "dggbal.f"
    b_offset = 1 + b_dim1;
#line 226 "dggbal.f"
    b -= b_offset;
#line 226 "dggbal.f"
    --lscale;
#line 226 "dggbal.f"
    --rscale;
#line 226 "dggbal.f"
    --work;
#line 226 "dggbal.f"

#line 226 "dggbal.f"
    /* Function Body */
#line 226 "dggbal.f"
    *info = 0;
#line 227 "dggbal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 229 "dggbal.f"
	*info = -1;
#line 230 "dggbal.f"
    } else if (*n < 0) {
#line 231 "dggbal.f"
	*info = -2;
#line 232 "dggbal.f"
    } else if (*lda < max(1,*n)) {
#line 233 "dggbal.f"
	*info = -4;
#line 234 "dggbal.f"
    } else if (*ldb < max(1,*n)) {
#line 235 "dggbal.f"
	*info = -6;
#line 236 "dggbal.f"
    }
#line 237 "dggbal.f"
    if (*info != 0) {
#line 238 "dggbal.f"
	i__1 = -(*info);
#line 238 "dggbal.f"
	xerbla_("DGGBAL", &i__1, (ftnlen)6);
#line 239 "dggbal.f"
	return 0;
#line 240 "dggbal.f"
    }

/*     Quick return if possible */

#line 244 "dggbal.f"
    if (*n == 0) {
#line 245 "dggbal.f"
	*ilo = 1;
#line 246 "dggbal.f"
	*ihi = *n;
#line 247 "dggbal.f"
	return 0;
#line 248 "dggbal.f"
    }

#line 250 "dggbal.f"
    if (*n == 1) {
#line 251 "dggbal.f"
	*ilo = 1;
#line 252 "dggbal.f"
	*ihi = *n;
#line 253 "dggbal.f"
	lscale[1] = 1.;
#line 254 "dggbal.f"
	rscale[1] = 1.;
#line 255 "dggbal.f"
	return 0;
#line 256 "dggbal.f"
    }

#line 258 "dggbal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 259 "dggbal.f"
	*ilo = 1;
#line 260 "dggbal.f"
	*ihi = *n;
#line 261 "dggbal.f"
	i__1 = *n;
#line 261 "dggbal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "dggbal.f"
	    lscale[i__] = 1.;
#line 263 "dggbal.f"
	    rscale[i__] = 1.;
#line 264 "dggbal.f"
/* L10: */
#line 264 "dggbal.f"
	}
#line 265 "dggbal.f"
	return 0;
#line 266 "dggbal.f"
    }

#line 268 "dggbal.f"
    k = 1;
#line 269 "dggbal.f"
    l = *n;
#line 270 "dggbal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 270 "dggbal.f"
	goto L190;
#line 270 "dggbal.f"
    }

#line 273 "dggbal.f"
    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues. */

/*     Find row with one nonzero in columns 1 through L */

#line 279 "dggbal.f"
L20:
#line 280 "dggbal.f"
    l = lm1;
#line 281 "dggbal.f"
    if (l != 1) {
#line 281 "dggbal.f"
	goto L30;
#line 281 "dggbal.f"
    }

#line 284 "dggbal.f"
    rscale[1] = 1.;
#line 285 "dggbal.f"
    lscale[1] = 1.;
#line 286 "dggbal.f"
    goto L190;

#line 288 "dggbal.f"
L30:
#line 289 "dggbal.f"
    lm1 = l - 1;
#line 290 "dggbal.f"
    for (i__ = l; i__ >= 1; --i__) {
#line 291 "dggbal.f"
	i__1 = lm1;
#line 291 "dggbal.f"
	for (j = 1; j <= i__1; ++j) {
#line 292 "dggbal.f"
	    jp1 = j + 1;
#line 293 "dggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 293 "dggbal.f"
		goto L50;
#line 293 "dggbal.f"
	    }
#line 295 "dggbal.f"
/* L40: */
#line 295 "dggbal.f"
	}
#line 296 "dggbal.f"
	j = l;
#line 297 "dggbal.f"
	goto L70;

#line 299 "dggbal.f"
L50:
#line 300 "dggbal.f"
	i__1 = l;
#line 300 "dggbal.f"
	for (j = jp1; j <= i__1; ++j) {
#line 301 "dggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 301 "dggbal.f"
		goto L80;
#line 301 "dggbal.f"
	    }
#line 303 "dggbal.f"
/* L60: */
#line 303 "dggbal.f"
	}
#line 304 "dggbal.f"
	j = jp1 - 1;

#line 306 "dggbal.f"
L70:
#line 307 "dggbal.f"
	m = l;
#line 308 "dggbal.f"
	iflow = 1;
#line 309 "dggbal.f"
	goto L160;
#line 310 "dggbal.f"
L80:
#line 310 "dggbal.f"
	;
#line 310 "dggbal.f"
    }
#line 311 "dggbal.f"
    goto L100;

/*     Find column with one nonzero in rows K through N */

#line 315 "dggbal.f"
L90:
#line 316 "dggbal.f"
    ++k;

#line 318 "dggbal.f"
L100:
#line 319 "dggbal.f"
    i__1 = l;
#line 319 "dggbal.f"
    for (j = k; j <= i__1; ++j) {
#line 320 "dggbal.f"
	i__2 = lm1;
#line 320 "dggbal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 321 "dggbal.f"
	    ip1 = i__ + 1;
#line 322 "dggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 322 "dggbal.f"
		goto L120;
#line 322 "dggbal.f"
	    }
#line 324 "dggbal.f"
/* L110: */
#line 324 "dggbal.f"
	}
#line 325 "dggbal.f"
	i__ = l;
#line 326 "dggbal.f"
	goto L140;
#line 327 "dggbal.f"
L120:
#line 328 "dggbal.f"
	i__2 = l;
#line 328 "dggbal.f"
	for (i__ = ip1; i__ <= i__2; ++i__) {
#line 329 "dggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 329 "dggbal.f"
		goto L150;
#line 329 "dggbal.f"
	    }
#line 331 "dggbal.f"
/* L130: */
#line 331 "dggbal.f"
	}
#line 332 "dggbal.f"
	i__ = ip1 - 1;
#line 333 "dggbal.f"
L140:
#line 334 "dggbal.f"
	m = k;
#line 335 "dggbal.f"
	iflow = 2;
#line 336 "dggbal.f"
	goto L160;
#line 337 "dggbal.f"
L150:
#line 337 "dggbal.f"
	;
#line 337 "dggbal.f"
    }
#line 338 "dggbal.f"
    goto L190;

/*     Permute rows M and I */

#line 342 "dggbal.f"
L160:
#line 343 "dggbal.f"
    lscale[m] = (doublereal) i__;
#line 344 "dggbal.f"
    if (i__ == m) {
#line 344 "dggbal.f"
	goto L170;
#line 344 "dggbal.f"
    }
#line 346 "dggbal.f"
    i__1 = *n - k + 1;
#line 346 "dggbal.f"
    dswap_(&i__1, &a[i__ + k * a_dim1], lda, &a[m + k * a_dim1], lda);
#line 347 "dggbal.f"
    i__1 = *n - k + 1;
#line 347 "dggbal.f"
    dswap_(&i__1, &b[i__ + k * b_dim1], ldb, &b[m + k * b_dim1], ldb);

/*     Permute columns M and J */

#line 351 "dggbal.f"
L170:
#line 352 "dggbal.f"
    rscale[m] = (doublereal) j;
#line 353 "dggbal.f"
    if (j == m) {
#line 353 "dggbal.f"
	goto L180;
#line 353 "dggbal.f"
    }
#line 355 "dggbal.f"
    dswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 356 "dggbal.f"
    dswap_(&l, &b[j * b_dim1 + 1], &c__1, &b[m * b_dim1 + 1], &c__1);

#line 358 "dggbal.f"
L180:
#line 359 "dggbal.f"
    switch (iflow) {
#line 359 "dggbal.f"
	case 1:  goto L20;
#line 359 "dggbal.f"
	case 2:  goto L90;
#line 359 "dggbal.f"
    }

#line 361 "dggbal.f"
L190:
#line 362 "dggbal.f"
    *ilo = k;
#line 363 "dggbal.f"
    *ihi = l;

#line 365 "dggbal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 366 "dggbal.f"
	i__1 = *ihi;
#line 366 "dggbal.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 367 "dggbal.f"
	    lscale[i__] = 1.;
#line 368 "dggbal.f"
	    rscale[i__] = 1.;
#line 369 "dggbal.f"
/* L195: */
#line 369 "dggbal.f"
	}
#line 370 "dggbal.f"
	return 0;
#line 371 "dggbal.f"
    }

#line 373 "dggbal.f"
    if (*ilo == *ihi) {
#line 373 "dggbal.f"
	return 0;
#line 373 "dggbal.f"
    }

/*     Balance the submatrix in rows ILO to IHI. */

#line 378 "dggbal.f"
    nr = *ihi - *ilo + 1;
#line 379 "dggbal.f"
    i__1 = *ihi;
#line 379 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 380 "dggbal.f"
	rscale[i__] = 0.;
#line 381 "dggbal.f"
	lscale[i__] = 0.;

#line 383 "dggbal.f"
	work[i__] = 0.;
#line 384 "dggbal.f"
	work[i__ + *n] = 0.;
#line 385 "dggbal.f"
	work[i__ + (*n << 1)] = 0.;
#line 386 "dggbal.f"
	work[i__ + *n * 3] = 0.;
#line 387 "dggbal.f"
	work[i__ + (*n << 2)] = 0.;
#line 388 "dggbal.f"
	work[i__ + *n * 5] = 0.;
#line 389 "dggbal.f"
/* L200: */
#line 389 "dggbal.f"
    }

/*     Compute right side vector in resulting linear equations */

#line 393 "dggbal.f"
    basl = d_lg10(&c_b35);
#line 394 "dggbal.f"
    i__1 = *ihi;
#line 394 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 395 "dggbal.f"
	i__2 = *ihi;
#line 395 "dggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 396 "dggbal.f"
	    tb = b[i__ + j * b_dim1];
#line 397 "dggbal.f"
	    ta = a[i__ + j * a_dim1];
#line 398 "dggbal.f"
	    if (ta == 0.) {
#line 398 "dggbal.f"
		goto L210;
#line 398 "dggbal.f"
	    }
#line 400 "dggbal.f"
	    d__1 = abs(ta);
#line 400 "dggbal.f"
	    ta = d_lg10(&d__1) / basl;
#line 401 "dggbal.f"
L210:
#line 402 "dggbal.f"
	    if (tb == 0.) {
#line 402 "dggbal.f"
		goto L220;
#line 402 "dggbal.f"
	    }
#line 404 "dggbal.f"
	    d__1 = abs(tb);
#line 404 "dggbal.f"
	    tb = d_lg10(&d__1) / basl;
#line 405 "dggbal.f"
L220:
#line 406 "dggbal.f"
	    work[i__ + (*n << 2)] = work[i__ + (*n << 2)] - ta - tb;
#line 407 "dggbal.f"
	    work[j + *n * 5] = work[j + *n * 5] - ta - tb;
#line 408 "dggbal.f"
/* L230: */
#line 408 "dggbal.f"
	}
#line 409 "dggbal.f"
/* L240: */
#line 409 "dggbal.f"
    }

#line 411 "dggbal.f"
    coef = 1. / (doublereal) (nr << 1);
#line 412 "dggbal.f"
    coef2 = coef * coef;
#line 413 "dggbal.f"
    coef5 = coef2 * .5;
#line 414 "dggbal.f"
    nrp2 = nr + 2;
#line 415 "dggbal.f"
    beta = 0.;
#line 416 "dggbal.f"
    it = 1;

/*     Start generalized conjugate gradient iteration */

#line 420 "dggbal.f"
L250:

#line 422 "dggbal.f"
    gamma = ddot_(&nr, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1) + ddot_(&nr, &work[*ilo + *n * 5], &c__1, &work[*ilo + *
	    n * 5], &c__1);

#line 425 "dggbal.f"
    ew = 0.;
#line 426 "dggbal.f"
    ewc = 0.;
#line 427 "dggbal.f"
    i__1 = *ihi;
#line 427 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 428 "dggbal.f"
	ew += work[i__ + (*n << 2)];
#line 429 "dggbal.f"
	ewc += work[i__ + *n * 5];
#line 430 "dggbal.f"
/* L260: */
#line 430 "dggbal.f"
    }

/* Computing 2nd power */
#line 432 "dggbal.f"
    d__1 = ew;
/* Computing 2nd power */
#line 432 "dggbal.f"
    d__2 = ewc;
/* Computing 2nd power */
#line 432 "dggbal.f"
    d__3 = ew - ewc;
#line 432 "dggbal.f"
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
#line 433 "dggbal.f"
    if (gamma == 0.) {
#line 433 "dggbal.f"
	goto L350;
#line 433 "dggbal.f"
    }
#line 435 "dggbal.f"
    if (it != 1) {
#line 435 "dggbal.f"
	beta = gamma / pgamma;
#line 435 "dggbal.f"
    }
#line 437 "dggbal.f"
    t = coef5 * (ewc - ew * 3.);
#line 438 "dggbal.f"
    tc = coef5 * (ew - ewc * 3.);

#line 440 "dggbal.f"
    dscal_(&nr, &beta, &work[*ilo], &c__1);
#line 441 "dggbal.f"
    dscal_(&nr, &beta, &work[*ilo + *n], &c__1);

#line 443 "dggbal.f"
    daxpy_(&nr, &coef, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + *n], &
	    c__1);
#line 444 "dggbal.f"
    daxpy_(&nr, &coef, &work[*ilo + *n * 5], &c__1, &work[*ilo], &c__1);

#line 446 "dggbal.f"
    i__1 = *ihi;
#line 446 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 447 "dggbal.f"
	work[i__] += tc;
#line 448 "dggbal.f"
	work[i__ + *n] += t;
#line 449 "dggbal.f"
/* L270: */
#line 449 "dggbal.f"
    }

/*     Apply matrix to vector */

#line 453 "dggbal.f"
    i__1 = *ihi;
#line 453 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 454 "dggbal.f"
	kount = 0;
#line 455 "dggbal.f"
	sum = 0.;
#line 456 "dggbal.f"
	i__2 = *ihi;
#line 456 "dggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 457 "dggbal.f"
	    if (a[i__ + j * a_dim1] == 0.) {
#line 457 "dggbal.f"
		goto L280;
#line 457 "dggbal.f"
	    }
#line 459 "dggbal.f"
	    ++kount;
#line 460 "dggbal.f"
	    sum += work[j];
#line 461 "dggbal.f"
L280:
#line 462 "dggbal.f"
	    if (b[i__ + j * b_dim1] == 0.) {
#line 462 "dggbal.f"
		goto L290;
#line 462 "dggbal.f"
	    }
#line 464 "dggbal.f"
	    ++kount;
#line 465 "dggbal.f"
	    sum += work[j];
#line 466 "dggbal.f"
L290:
#line 466 "dggbal.f"
	    ;
#line 466 "dggbal.f"
	}
#line 467 "dggbal.f"
	work[i__ + (*n << 1)] = (doublereal) kount * work[i__ + *n] + sum;
#line 468 "dggbal.f"
/* L300: */
#line 468 "dggbal.f"
    }

#line 470 "dggbal.f"
    i__1 = *ihi;
#line 470 "dggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 471 "dggbal.f"
	kount = 0;
#line 472 "dggbal.f"
	sum = 0.;
#line 473 "dggbal.f"
	i__2 = *ihi;
#line 473 "dggbal.f"
	for (i__ = *ilo; i__ <= i__2; ++i__) {
#line 474 "dggbal.f"
	    if (a[i__ + j * a_dim1] == 0.) {
#line 474 "dggbal.f"
		goto L310;
#line 474 "dggbal.f"
	    }
#line 476 "dggbal.f"
	    ++kount;
#line 477 "dggbal.f"
	    sum += work[i__ + *n];
#line 478 "dggbal.f"
L310:
#line 479 "dggbal.f"
	    if (b[i__ + j * b_dim1] == 0.) {
#line 479 "dggbal.f"
		goto L320;
#line 479 "dggbal.f"
	    }
#line 481 "dggbal.f"
	    ++kount;
#line 482 "dggbal.f"
	    sum += work[i__ + *n];
#line 483 "dggbal.f"
L320:
#line 483 "dggbal.f"
	    ;
#line 483 "dggbal.f"
	}
#line 484 "dggbal.f"
	work[j + *n * 3] = (doublereal) kount * work[j] + sum;
#line 485 "dggbal.f"
/* L330: */
#line 485 "dggbal.f"
    }

#line 487 "dggbal.f"
    sum = ddot_(&nr, &work[*ilo + *n], &c__1, &work[*ilo + (*n << 1)], &c__1) 
	    + ddot_(&nr, &work[*ilo], &c__1, &work[*ilo + *n * 3], &c__1);
#line 489 "dggbal.f"
    alpha = gamma / sum;

/*     Determine correction to current iteration */

#line 493 "dggbal.f"
    cmax = 0.;
#line 494 "dggbal.f"
    i__1 = *ihi;
#line 494 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 495 "dggbal.f"
	cor = alpha * work[i__ + *n];
#line 496 "dggbal.f"
	if (abs(cor) > cmax) {
#line 496 "dggbal.f"
	    cmax = abs(cor);
#line 496 "dggbal.f"
	}
#line 498 "dggbal.f"
	lscale[i__] += cor;
#line 499 "dggbal.f"
	cor = alpha * work[i__];
#line 500 "dggbal.f"
	if (abs(cor) > cmax) {
#line 500 "dggbal.f"
	    cmax = abs(cor);
#line 500 "dggbal.f"
	}
#line 502 "dggbal.f"
	rscale[i__] += cor;
#line 503 "dggbal.f"
/* L340: */
#line 503 "dggbal.f"
    }
#line 504 "dggbal.f"
    if (cmax < .5) {
#line 504 "dggbal.f"
	goto L350;
#line 504 "dggbal.f"
    }

#line 507 "dggbal.f"
    d__1 = -alpha;
#line 507 "dggbal.f"
    daxpy_(&nr, &d__1, &work[*ilo + (*n << 1)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1);
#line 508 "dggbal.f"
    d__1 = -alpha;
#line 508 "dggbal.f"
    daxpy_(&nr, &d__1, &work[*ilo + *n * 3], &c__1, &work[*ilo + *n * 5], &
	    c__1);

#line 510 "dggbal.f"
    pgamma = gamma;
#line 511 "dggbal.f"
    ++it;
#line 512 "dggbal.f"
    if (it <= nrp2) {
#line 512 "dggbal.f"
	goto L250;
#line 512 "dggbal.f"
    }

/*     End generalized conjugate gradient iteration */

#line 517 "dggbal.f"
L350:
#line 518 "dggbal.f"
    sfmin = dlamch_("S", (ftnlen)1);
#line 519 "dggbal.f"
    sfmax = 1. / sfmin;
#line 520 "dggbal.f"
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
#line 521 "dggbal.f"
    lsfmax = (integer) (d_lg10(&sfmax) / basl);
#line 522 "dggbal.f"
    i__1 = *ihi;
#line 522 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 523 "dggbal.f"
	i__2 = *n - *ilo + 1;
#line 523 "dggbal.f"
	irab = idamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
#line 524 "dggbal.f"
	rab = (d__1 = a[i__ + (irab + *ilo - 1) * a_dim1], abs(d__1));
#line 525 "dggbal.f"
	i__2 = *n - *ilo + 1;
#line 525 "dggbal.f"
	irab = idamax_(&i__2, &b[i__ + *ilo * b_dim1], ldb);
/* Computing MAX */
#line 526 "dggbal.f"
	d__2 = rab, d__3 = (d__1 = b[i__ + (irab + *ilo - 1) * b_dim1], abs(
		d__1));
#line 526 "dggbal.f"
	rab = max(d__2,d__3);
#line 527 "dggbal.f"
	d__1 = rab + sfmin;
#line 527 "dggbal.f"
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 528 "dggbal.f"
	ir = (integer) (lscale[i__] + d_sign(&c_b71, &lscale[i__]));
/* Computing MIN */
#line 529 "dggbal.f"
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
#line 529 "dggbal.f"
	ir = min(i__2,i__3);
#line 530 "dggbal.f"
	lscale[i__] = pow_di(&c_b35, &ir);
#line 531 "dggbal.f"
	icab = idamax_(ihi, &a[i__ * a_dim1 + 1], &c__1);
#line 532 "dggbal.f"
	cab = (d__1 = a[icab + i__ * a_dim1], abs(d__1));
#line 533 "dggbal.f"
	icab = idamax_(ihi, &b[i__ * b_dim1 + 1], &c__1);
/* Computing MAX */
#line 534 "dggbal.f"
	d__2 = cab, d__3 = (d__1 = b[icab + i__ * b_dim1], abs(d__1));
#line 534 "dggbal.f"
	cab = max(d__2,d__3);
#line 535 "dggbal.f"
	d__1 = cab + sfmin;
#line 535 "dggbal.f"
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 536 "dggbal.f"
	jc = (integer) (rscale[i__] + d_sign(&c_b71, &rscale[i__]));
/* Computing MIN */
#line 537 "dggbal.f"
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
#line 537 "dggbal.f"
	jc = min(i__2,i__3);
#line 538 "dggbal.f"
	rscale[i__] = pow_di(&c_b35, &jc);
#line 539 "dggbal.f"
/* L360: */
#line 539 "dggbal.f"
    }

/*     Row scaling of matrices A and B */

#line 543 "dggbal.f"
    i__1 = *ihi;
#line 543 "dggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 544 "dggbal.f"
	i__2 = *n - *ilo + 1;
#line 544 "dggbal.f"
	dscal_(&i__2, &lscale[i__], &a[i__ + *ilo * a_dim1], lda);
#line 545 "dggbal.f"
	i__2 = *n - *ilo + 1;
#line 545 "dggbal.f"
	dscal_(&i__2, &lscale[i__], &b[i__ + *ilo * b_dim1], ldb);
#line 546 "dggbal.f"
/* L370: */
#line 546 "dggbal.f"
    }

/*     Column scaling of matrices A and B */

#line 550 "dggbal.f"
    i__1 = *ihi;
#line 550 "dggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 551 "dggbal.f"
	dscal_(ihi, &rscale[j], &a[j * a_dim1 + 1], &c__1);
#line 552 "dggbal.f"
	dscal_(ihi, &rscale[j], &b[j * b_dim1 + 1], &c__1);
#line 553 "dggbal.f"
/* L380: */
#line 553 "dggbal.f"
    }

#line 555 "dggbal.f"
    return 0;

/*     End of DGGBAL */

} /* dggbal_ */

