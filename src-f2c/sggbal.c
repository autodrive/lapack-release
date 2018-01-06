#line 1 "sggbal.f"
/* sggbal.f -- translated by f2c (version 20100827).
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

#line 1 "sggbal.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b35 = 10.;
static doublereal c_b71 = .5;

/* > \brief \b SGGBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggbal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggbal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggbal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, */
/*                          RSCALE, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), LSCALE( * ), */
/*      $                   RSCALE( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGBAL balances a pair of general real matrices (A,B).  This */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          B is REAL array, dimension (LDB,N) */
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
/* >          LSCALE is REAL array, dimension (N) */
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
/* >          RSCALE is REAL array, dimension (N) */
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
/* >          WORK is REAL array, dimension (lwork) */
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

/* > \date November 2011 */

/* > \ingroup realGBcomputational */

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
/* Subroutine */ int sggbal_(char *job, integer *n, doublereal *a, integer *
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
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal coef2, coef5, gamma, alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sfmin, sfmax;
    static integer iflow;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer kount;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal pgamma;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static integer lsfmin, lsfmax;


/*  -- LAPACK computational routine (version 3.4.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 226 "sggbal.f"
    /* Parameter adjustments */
#line 226 "sggbal.f"
    a_dim1 = *lda;
#line 226 "sggbal.f"
    a_offset = 1 + a_dim1;
#line 226 "sggbal.f"
    a -= a_offset;
#line 226 "sggbal.f"
    b_dim1 = *ldb;
#line 226 "sggbal.f"
    b_offset = 1 + b_dim1;
#line 226 "sggbal.f"
    b -= b_offset;
#line 226 "sggbal.f"
    --lscale;
#line 226 "sggbal.f"
    --rscale;
#line 226 "sggbal.f"
    --work;
#line 226 "sggbal.f"

#line 226 "sggbal.f"
    /* Function Body */
#line 226 "sggbal.f"
    *info = 0;
#line 227 "sggbal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 229 "sggbal.f"
	*info = -1;
#line 230 "sggbal.f"
    } else if (*n < 0) {
#line 231 "sggbal.f"
	*info = -2;
#line 232 "sggbal.f"
    } else if (*lda < max(1,*n)) {
#line 233 "sggbal.f"
	*info = -4;
#line 234 "sggbal.f"
    } else if (*ldb < max(1,*n)) {
#line 235 "sggbal.f"
	*info = -6;
#line 236 "sggbal.f"
    }
#line 237 "sggbal.f"
    if (*info != 0) {
#line 238 "sggbal.f"
	i__1 = -(*info);
#line 238 "sggbal.f"
	xerbla_("SGGBAL", &i__1, (ftnlen)6);
#line 239 "sggbal.f"
	return 0;
#line 240 "sggbal.f"
    }

/*     Quick return if possible */

#line 244 "sggbal.f"
    if (*n == 0) {
#line 245 "sggbal.f"
	*ilo = 1;
#line 246 "sggbal.f"
	*ihi = *n;
#line 247 "sggbal.f"
	return 0;
#line 248 "sggbal.f"
    }

#line 250 "sggbal.f"
    if (*n == 1) {
#line 251 "sggbal.f"
	*ilo = 1;
#line 252 "sggbal.f"
	*ihi = *n;
#line 253 "sggbal.f"
	lscale[1] = 1.;
#line 254 "sggbal.f"
	rscale[1] = 1.;
#line 255 "sggbal.f"
	return 0;
#line 256 "sggbal.f"
    }

#line 258 "sggbal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 259 "sggbal.f"
	*ilo = 1;
#line 260 "sggbal.f"
	*ihi = *n;
#line 261 "sggbal.f"
	i__1 = *n;
#line 261 "sggbal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "sggbal.f"
	    lscale[i__] = 1.;
#line 263 "sggbal.f"
	    rscale[i__] = 1.;
#line 264 "sggbal.f"
/* L10: */
#line 264 "sggbal.f"
	}
#line 265 "sggbal.f"
	return 0;
#line 266 "sggbal.f"
    }

#line 268 "sggbal.f"
    k = 1;
#line 269 "sggbal.f"
    l = *n;
#line 270 "sggbal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 270 "sggbal.f"
	goto L190;
#line 270 "sggbal.f"
    }

#line 273 "sggbal.f"
    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues. */

/*     Find row with one nonzero in columns 1 through L */

#line 279 "sggbal.f"
L20:
#line 280 "sggbal.f"
    l = lm1;
#line 281 "sggbal.f"
    if (l != 1) {
#line 281 "sggbal.f"
	goto L30;
#line 281 "sggbal.f"
    }

#line 284 "sggbal.f"
    rscale[1] = 1.;
#line 285 "sggbal.f"
    lscale[1] = 1.;
#line 286 "sggbal.f"
    goto L190;

#line 288 "sggbal.f"
L30:
#line 289 "sggbal.f"
    lm1 = l - 1;
#line 290 "sggbal.f"
    for (i__ = l; i__ >= 1; --i__) {
#line 291 "sggbal.f"
	i__1 = lm1;
#line 291 "sggbal.f"
	for (j = 1; j <= i__1; ++j) {
#line 292 "sggbal.f"
	    jp1 = j + 1;
#line 293 "sggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 293 "sggbal.f"
		goto L50;
#line 293 "sggbal.f"
	    }
#line 295 "sggbal.f"
/* L40: */
#line 295 "sggbal.f"
	}
#line 296 "sggbal.f"
	j = l;
#line 297 "sggbal.f"
	goto L70;

#line 299 "sggbal.f"
L50:
#line 300 "sggbal.f"
	i__1 = l;
#line 300 "sggbal.f"
	for (j = jp1; j <= i__1; ++j) {
#line 301 "sggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 301 "sggbal.f"
		goto L80;
#line 301 "sggbal.f"
	    }
#line 303 "sggbal.f"
/* L60: */
#line 303 "sggbal.f"
	}
#line 304 "sggbal.f"
	j = jp1 - 1;

#line 306 "sggbal.f"
L70:
#line 307 "sggbal.f"
	m = l;
#line 308 "sggbal.f"
	iflow = 1;
#line 309 "sggbal.f"
	goto L160;
#line 310 "sggbal.f"
L80:
#line 310 "sggbal.f"
	;
#line 310 "sggbal.f"
    }
#line 311 "sggbal.f"
    goto L100;

/*     Find column with one nonzero in rows K through N */

#line 315 "sggbal.f"
L90:
#line 316 "sggbal.f"
    ++k;

#line 318 "sggbal.f"
L100:
#line 319 "sggbal.f"
    i__1 = l;
#line 319 "sggbal.f"
    for (j = k; j <= i__1; ++j) {
#line 320 "sggbal.f"
	i__2 = lm1;
#line 320 "sggbal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 321 "sggbal.f"
	    ip1 = i__ + 1;
#line 322 "sggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 322 "sggbal.f"
		goto L120;
#line 322 "sggbal.f"
	    }
#line 324 "sggbal.f"
/* L110: */
#line 324 "sggbal.f"
	}
#line 325 "sggbal.f"
	i__ = l;
#line 326 "sggbal.f"
	goto L140;
#line 327 "sggbal.f"
L120:
#line 328 "sggbal.f"
	i__2 = l;
#line 328 "sggbal.f"
	for (i__ = ip1; i__ <= i__2; ++i__) {
#line 329 "sggbal.f"
	    if (a[i__ + j * a_dim1] != 0. || b[i__ + j * b_dim1] != 0.) {
#line 329 "sggbal.f"
		goto L150;
#line 329 "sggbal.f"
	    }
#line 331 "sggbal.f"
/* L130: */
#line 331 "sggbal.f"
	}
#line 332 "sggbal.f"
	i__ = ip1 - 1;
#line 333 "sggbal.f"
L140:
#line 334 "sggbal.f"
	m = k;
#line 335 "sggbal.f"
	iflow = 2;
#line 336 "sggbal.f"
	goto L160;
#line 337 "sggbal.f"
L150:
#line 337 "sggbal.f"
	;
#line 337 "sggbal.f"
    }
#line 338 "sggbal.f"
    goto L190;

/*     Permute rows M and I */

#line 342 "sggbal.f"
L160:
#line 343 "sggbal.f"
    lscale[m] = (doublereal) i__;
#line 344 "sggbal.f"
    if (i__ == m) {
#line 344 "sggbal.f"
	goto L170;
#line 344 "sggbal.f"
    }
#line 346 "sggbal.f"
    i__1 = *n - k + 1;
#line 346 "sggbal.f"
    sswap_(&i__1, &a[i__ + k * a_dim1], lda, &a[m + k * a_dim1], lda);
#line 347 "sggbal.f"
    i__1 = *n - k + 1;
#line 347 "sggbal.f"
    sswap_(&i__1, &b[i__ + k * b_dim1], ldb, &b[m + k * b_dim1], ldb);

/*     Permute columns M and J */

#line 351 "sggbal.f"
L170:
#line 352 "sggbal.f"
    rscale[m] = (doublereal) j;
#line 353 "sggbal.f"
    if (j == m) {
#line 353 "sggbal.f"
	goto L180;
#line 353 "sggbal.f"
    }
#line 355 "sggbal.f"
    sswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 356 "sggbal.f"
    sswap_(&l, &b[j * b_dim1 + 1], &c__1, &b[m * b_dim1 + 1], &c__1);

#line 358 "sggbal.f"
L180:
#line 359 "sggbal.f"
    switch (iflow) {
#line 359 "sggbal.f"
	case 1:  goto L20;
#line 359 "sggbal.f"
	case 2:  goto L90;
#line 359 "sggbal.f"
    }

#line 361 "sggbal.f"
L190:
#line 362 "sggbal.f"
    *ilo = k;
#line 363 "sggbal.f"
    *ihi = l;

#line 365 "sggbal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 366 "sggbal.f"
	i__1 = *ihi;
#line 366 "sggbal.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 367 "sggbal.f"
	    lscale[i__] = 1.;
#line 368 "sggbal.f"
	    rscale[i__] = 1.;
#line 369 "sggbal.f"
/* L195: */
#line 369 "sggbal.f"
	}
#line 370 "sggbal.f"
	return 0;
#line 371 "sggbal.f"
    }

#line 373 "sggbal.f"
    if (*ilo == *ihi) {
#line 373 "sggbal.f"
	return 0;
#line 373 "sggbal.f"
    }

/*     Balance the submatrix in rows ILO to IHI. */

#line 378 "sggbal.f"
    nr = *ihi - *ilo + 1;
#line 379 "sggbal.f"
    i__1 = *ihi;
#line 379 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 380 "sggbal.f"
	rscale[i__] = 0.;
#line 381 "sggbal.f"
	lscale[i__] = 0.;

#line 383 "sggbal.f"
	work[i__] = 0.;
#line 384 "sggbal.f"
	work[i__ + *n] = 0.;
#line 385 "sggbal.f"
	work[i__ + (*n << 1)] = 0.;
#line 386 "sggbal.f"
	work[i__ + *n * 3] = 0.;
#line 387 "sggbal.f"
	work[i__ + (*n << 2)] = 0.;
#line 388 "sggbal.f"
	work[i__ + *n * 5] = 0.;
#line 389 "sggbal.f"
/* L200: */
#line 389 "sggbal.f"
    }

/*     Compute right side vector in resulting linear equations */

#line 393 "sggbal.f"
    basl = d_lg10(&c_b35);
#line 394 "sggbal.f"
    i__1 = *ihi;
#line 394 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 395 "sggbal.f"
	i__2 = *ihi;
#line 395 "sggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 396 "sggbal.f"
	    tb = b[i__ + j * b_dim1];
#line 397 "sggbal.f"
	    ta = a[i__ + j * a_dim1];
#line 398 "sggbal.f"
	    if (ta == 0.) {
#line 398 "sggbal.f"
		goto L210;
#line 398 "sggbal.f"
	    }
#line 400 "sggbal.f"
	    d__1 = abs(ta);
#line 400 "sggbal.f"
	    ta = d_lg10(&d__1) / basl;
#line 401 "sggbal.f"
L210:
#line 402 "sggbal.f"
	    if (tb == 0.) {
#line 402 "sggbal.f"
		goto L220;
#line 402 "sggbal.f"
	    }
#line 404 "sggbal.f"
	    d__1 = abs(tb);
#line 404 "sggbal.f"
	    tb = d_lg10(&d__1) / basl;
#line 405 "sggbal.f"
L220:
#line 406 "sggbal.f"
	    work[i__ + (*n << 2)] = work[i__ + (*n << 2)] - ta - tb;
#line 407 "sggbal.f"
	    work[j + *n * 5] = work[j + *n * 5] - ta - tb;
#line 408 "sggbal.f"
/* L230: */
#line 408 "sggbal.f"
	}
#line 409 "sggbal.f"
/* L240: */
#line 409 "sggbal.f"
    }

#line 411 "sggbal.f"
    coef = 1. / (doublereal) (nr << 1);
#line 412 "sggbal.f"
    coef2 = coef * coef;
#line 413 "sggbal.f"
    coef5 = coef2 * .5;
#line 414 "sggbal.f"
    nrp2 = nr + 2;
#line 415 "sggbal.f"
    beta = 0.;
#line 416 "sggbal.f"
    it = 1;

/*     Start generalized conjugate gradient iteration */

#line 420 "sggbal.f"
L250:

#line 422 "sggbal.f"
    gamma = sdot_(&nr, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1) + sdot_(&nr, &work[*ilo + *n * 5], &c__1, &work[*ilo + *
	    n * 5], &c__1);

#line 425 "sggbal.f"
    ew = 0.;
#line 426 "sggbal.f"
    ewc = 0.;
#line 427 "sggbal.f"
    i__1 = *ihi;
#line 427 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 428 "sggbal.f"
	ew += work[i__ + (*n << 2)];
#line 429 "sggbal.f"
	ewc += work[i__ + *n * 5];
#line 430 "sggbal.f"
/* L260: */
#line 430 "sggbal.f"
    }

/* Computing 2nd power */
#line 432 "sggbal.f"
    d__1 = ew;
/* Computing 2nd power */
#line 432 "sggbal.f"
    d__2 = ewc;
/* Computing 2nd power */
#line 432 "sggbal.f"
    d__3 = ew - ewc;
#line 432 "sggbal.f"
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
#line 433 "sggbal.f"
    if (gamma == 0.) {
#line 433 "sggbal.f"
	goto L350;
#line 433 "sggbal.f"
    }
#line 435 "sggbal.f"
    if (it != 1) {
#line 435 "sggbal.f"
	beta = gamma / pgamma;
#line 435 "sggbal.f"
    }
#line 437 "sggbal.f"
    t = coef5 * (ewc - ew * 3.);
#line 438 "sggbal.f"
    tc = coef5 * (ew - ewc * 3.);

#line 440 "sggbal.f"
    sscal_(&nr, &beta, &work[*ilo], &c__1);
#line 441 "sggbal.f"
    sscal_(&nr, &beta, &work[*ilo + *n], &c__1);

#line 443 "sggbal.f"
    saxpy_(&nr, &coef, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + *n], &
	    c__1);
#line 444 "sggbal.f"
    saxpy_(&nr, &coef, &work[*ilo + *n * 5], &c__1, &work[*ilo], &c__1);

#line 446 "sggbal.f"
    i__1 = *ihi;
#line 446 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 447 "sggbal.f"
	work[i__] += tc;
#line 448 "sggbal.f"
	work[i__ + *n] += t;
#line 449 "sggbal.f"
/* L270: */
#line 449 "sggbal.f"
    }

/*     Apply matrix to vector */

#line 453 "sggbal.f"
    i__1 = *ihi;
#line 453 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 454 "sggbal.f"
	kount = 0;
#line 455 "sggbal.f"
	sum = 0.;
#line 456 "sggbal.f"
	i__2 = *ihi;
#line 456 "sggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 457 "sggbal.f"
	    if (a[i__ + j * a_dim1] == 0.) {
#line 457 "sggbal.f"
		goto L280;
#line 457 "sggbal.f"
	    }
#line 459 "sggbal.f"
	    ++kount;
#line 460 "sggbal.f"
	    sum += work[j];
#line 461 "sggbal.f"
L280:
#line 462 "sggbal.f"
	    if (b[i__ + j * b_dim1] == 0.) {
#line 462 "sggbal.f"
		goto L290;
#line 462 "sggbal.f"
	    }
#line 464 "sggbal.f"
	    ++kount;
#line 465 "sggbal.f"
	    sum += work[j];
#line 466 "sggbal.f"
L290:
#line 466 "sggbal.f"
	    ;
#line 466 "sggbal.f"
	}
#line 467 "sggbal.f"
	work[i__ + (*n << 1)] = (doublereal) kount * work[i__ + *n] + sum;
#line 468 "sggbal.f"
/* L300: */
#line 468 "sggbal.f"
    }

#line 470 "sggbal.f"
    i__1 = *ihi;
#line 470 "sggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 471 "sggbal.f"
	kount = 0;
#line 472 "sggbal.f"
	sum = 0.;
#line 473 "sggbal.f"
	i__2 = *ihi;
#line 473 "sggbal.f"
	for (i__ = *ilo; i__ <= i__2; ++i__) {
#line 474 "sggbal.f"
	    if (a[i__ + j * a_dim1] == 0.) {
#line 474 "sggbal.f"
		goto L310;
#line 474 "sggbal.f"
	    }
#line 476 "sggbal.f"
	    ++kount;
#line 477 "sggbal.f"
	    sum += work[i__ + *n];
#line 478 "sggbal.f"
L310:
#line 479 "sggbal.f"
	    if (b[i__ + j * b_dim1] == 0.) {
#line 479 "sggbal.f"
		goto L320;
#line 479 "sggbal.f"
	    }
#line 481 "sggbal.f"
	    ++kount;
#line 482 "sggbal.f"
	    sum += work[i__ + *n];
#line 483 "sggbal.f"
L320:
#line 483 "sggbal.f"
	    ;
#line 483 "sggbal.f"
	}
#line 484 "sggbal.f"
	work[j + *n * 3] = (doublereal) kount * work[j] + sum;
#line 485 "sggbal.f"
/* L330: */
#line 485 "sggbal.f"
    }

#line 487 "sggbal.f"
    sum = sdot_(&nr, &work[*ilo + *n], &c__1, &work[*ilo + (*n << 1)], &c__1) 
	    + sdot_(&nr, &work[*ilo], &c__1, &work[*ilo + *n * 3], &c__1);
#line 489 "sggbal.f"
    alpha = gamma / sum;

/*     Determine correction to current iteration */

#line 493 "sggbal.f"
    cmax = 0.;
#line 494 "sggbal.f"
    i__1 = *ihi;
#line 494 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 495 "sggbal.f"
	cor = alpha * work[i__ + *n];
#line 496 "sggbal.f"
	if (abs(cor) > cmax) {
#line 496 "sggbal.f"
	    cmax = abs(cor);
#line 496 "sggbal.f"
	}
#line 498 "sggbal.f"
	lscale[i__] += cor;
#line 499 "sggbal.f"
	cor = alpha * work[i__];
#line 500 "sggbal.f"
	if (abs(cor) > cmax) {
#line 500 "sggbal.f"
	    cmax = abs(cor);
#line 500 "sggbal.f"
	}
#line 502 "sggbal.f"
	rscale[i__] += cor;
#line 503 "sggbal.f"
/* L340: */
#line 503 "sggbal.f"
    }
#line 504 "sggbal.f"
    if (cmax < .5) {
#line 504 "sggbal.f"
	goto L350;
#line 504 "sggbal.f"
    }

#line 507 "sggbal.f"
    d__1 = -alpha;
#line 507 "sggbal.f"
    saxpy_(&nr, &d__1, &work[*ilo + (*n << 1)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1);
#line 508 "sggbal.f"
    d__1 = -alpha;
#line 508 "sggbal.f"
    saxpy_(&nr, &d__1, &work[*ilo + *n * 3], &c__1, &work[*ilo + *n * 5], &
	    c__1);

#line 510 "sggbal.f"
    pgamma = gamma;
#line 511 "sggbal.f"
    ++it;
#line 512 "sggbal.f"
    if (it <= nrp2) {
#line 512 "sggbal.f"
	goto L250;
#line 512 "sggbal.f"
    }

/*     End generalized conjugate gradient iteration */

#line 517 "sggbal.f"
L350:
#line 518 "sggbal.f"
    sfmin = slamch_("S", (ftnlen)1);
#line 519 "sggbal.f"
    sfmax = 1. / sfmin;
#line 520 "sggbal.f"
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
#line 521 "sggbal.f"
    lsfmax = (integer) (d_lg10(&sfmax) / basl);
#line 522 "sggbal.f"
    i__1 = *ihi;
#line 522 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 523 "sggbal.f"
	i__2 = *n - *ilo + 1;
#line 523 "sggbal.f"
	irab = isamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
#line 524 "sggbal.f"
	rab = (d__1 = a[i__ + (irab + *ilo - 1) * a_dim1], abs(d__1));
#line 525 "sggbal.f"
	i__2 = *n - *ilo + 1;
#line 525 "sggbal.f"
	irab = isamax_(&i__2, &b[i__ + *ilo * b_dim1], ldb);
/* Computing MAX */
#line 526 "sggbal.f"
	d__2 = rab, d__3 = (d__1 = b[i__ + (irab + *ilo - 1) * b_dim1], abs(
		d__1));
#line 526 "sggbal.f"
	rab = max(d__2,d__3);
#line 527 "sggbal.f"
	d__1 = rab + sfmin;
#line 527 "sggbal.f"
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 528 "sggbal.f"
	ir = (integer) (lscale[i__] + d_sign(&c_b71, &lscale[i__]));
/* Computing MIN */
#line 529 "sggbal.f"
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
#line 529 "sggbal.f"
	ir = min(i__2,i__3);
#line 530 "sggbal.f"
	lscale[i__] = pow_di(&c_b35, &ir);
#line 531 "sggbal.f"
	icab = isamax_(ihi, &a[i__ * a_dim1 + 1], &c__1);
#line 532 "sggbal.f"
	cab = (d__1 = a[icab + i__ * a_dim1], abs(d__1));
#line 533 "sggbal.f"
	icab = isamax_(ihi, &b[i__ * b_dim1 + 1], &c__1);
/* Computing MAX */
#line 534 "sggbal.f"
	d__2 = cab, d__3 = (d__1 = b[icab + i__ * b_dim1], abs(d__1));
#line 534 "sggbal.f"
	cab = max(d__2,d__3);
#line 535 "sggbal.f"
	d__1 = cab + sfmin;
#line 535 "sggbal.f"
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 536 "sggbal.f"
	jc = (integer) (rscale[i__] + d_sign(&c_b71, &rscale[i__]));
/* Computing MIN */
#line 537 "sggbal.f"
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
#line 537 "sggbal.f"
	jc = min(i__2,i__3);
#line 538 "sggbal.f"
	rscale[i__] = pow_di(&c_b35, &jc);
#line 539 "sggbal.f"
/* L360: */
#line 539 "sggbal.f"
    }

/*     Row scaling of matrices A and B */

#line 543 "sggbal.f"
    i__1 = *ihi;
#line 543 "sggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 544 "sggbal.f"
	i__2 = *n - *ilo + 1;
#line 544 "sggbal.f"
	sscal_(&i__2, &lscale[i__], &a[i__ + *ilo * a_dim1], lda);
#line 545 "sggbal.f"
	i__2 = *n - *ilo + 1;
#line 545 "sggbal.f"
	sscal_(&i__2, &lscale[i__], &b[i__ + *ilo * b_dim1], ldb);
#line 546 "sggbal.f"
/* L370: */
#line 546 "sggbal.f"
    }

/*     Column scaling of matrices A and B */

#line 550 "sggbal.f"
    i__1 = *ihi;
#line 550 "sggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 551 "sggbal.f"
	sscal_(ihi, &rscale[j], &a[j * a_dim1 + 1], &c__1);
#line 552 "sggbal.f"
	sscal_(ihi, &rscale[j], &b[j * b_dim1 + 1], &c__1);
#line 553 "sggbal.f"
/* L380: */
#line 553 "sggbal.f"
    }

#line 555 "sggbal.f"
    return 0;

/*     End of SGGBAL */

} /* sggbal_ */

