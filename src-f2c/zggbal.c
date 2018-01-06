#line 1 "zggbal.f"
/* zggbal.f -- translated by f2c (version 20100827).
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

#line 1 "zggbal.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = 10.;
static doublereal c_b72 = .5;

/* > \brief \b ZGGBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggbal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggbal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggbal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, */
/*                          RSCALE, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), WORK( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGBAL balances a pair of general complex matrices (A,B).  This */
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
/* >                  and RSCALE(I) = 1.0 for i=1,...,N; */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the input matrix A. */
/* >          On exit, A is overwritten by the balanced matrix. */
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
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On entry, the input matrix B. */
/* >          On exit, B is overwritten by the balanced matrix. */
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
/* >          row interchanged with row j, and D(j) is the scaling factor */
/* >          applied to row j, then */
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
/* >          column interchanged with column j, and D(j) is the scaling */
/* >          factor applied to column j, then */
/* >            RSCALE(j) = P(j)    for J = 1,...,ILO-1 */
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

/* > \date June 2016 */

/* > \ingroup complex16GBcomputational */

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
/* Subroutine */ int zggbal_(char *job, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_lg10(doublereal *), d_imag(doublecomplex *), z_abs(doublecomplex 
	    *), d_sign(doublereal *, doublereal *), pow_di(doublereal *, 
	    integer *);

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
    static integer iflow;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer kount;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal pgamma;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static integer lsfmin;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static integer lsfmax;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 235 "zggbal.f"
    /* Parameter adjustments */
#line 235 "zggbal.f"
    a_dim1 = *lda;
#line 235 "zggbal.f"
    a_offset = 1 + a_dim1;
#line 235 "zggbal.f"
    a -= a_offset;
#line 235 "zggbal.f"
    b_dim1 = *ldb;
#line 235 "zggbal.f"
    b_offset = 1 + b_dim1;
#line 235 "zggbal.f"
    b -= b_offset;
#line 235 "zggbal.f"
    --lscale;
#line 235 "zggbal.f"
    --rscale;
#line 235 "zggbal.f"
    --work;
#line 235 "zggbal.f"

#line 235 "zggbal.f"
    /* Function Body */
#line 235 "zggbal.f"
    *info = 0;
#line 236 "zggbal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 238 "zggbal.f"
	*info = -1;
#line 239 "zggbal.f"
    } else if (*n < 0) {
#line 240 "zggbal.f"
	*info = -2;
#line 241 "zggbal.f"
    } else if (*lda < max(1,*n)) {
#line 242 "zggbal.f"
	*info = -4;
#line 243 "zggbal.f"
    } else if (*ldb < max(1,*n)) {
#line 244 "zggbal.f"
	*info = -6;
#line 245 "zggbal.f"
    }
#line 246 "zggbal.f"
    if (*info != 0) {
#line 247 "zggbal.f"
	i__1 = -(*info);
#line 247 "zggbal.f"
	xerbla_("ZGGBAL", &i__1, (ftnlen)6);
#line 248 "zggbal.f"
	return 0;
#line 249 "zggbal.f"
    }

/*     Quick return if possible */

#line 253 "zggbal.f"
    if (*n == 0) {
#line 254 "zggbal.f"
	*ilo = 1;
#line 255 "zggbal.f"
	*ihi = *n;
#line 256 "zggbal.f"
	return 0;
#line 257 "zggbal.f"
    }

#line 259 "zggbal.f"
    if (*n == 1) {
#line 260 "zggbal.f"
	*ilo = 1;
#line 261 "zggbal.f"
	*ihi = *n;
#line 262 "zggbal.f"
	lscale[1] = 1.;
#line 263 "zggbal.f"
	rscale[1] = 1.;
#line 264 "zggbal.f"
	return 0;
#line 265 "zggbal.f"
    }

#line 267 "zggbal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 268 "zggbal.f"
	*ilo = 1;
#line 269 "zggbal.f"
	*ihi = *n;
#line 270 "zggbal.f"
	i__1 = *n;
#line 270 "zggbal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "zggbal.f"
	    lscale[i__] = 1.;
#line 272 "zggbal.f"
	    rscale[i__] = 1.;
#line 273 "zggbal.f"
/* L10: */
#line 273 "zggbal.f"
	}
#line 274 "zggbal.f"
	return 0;
#line 275 "zggbal.f"
    }

#line 277 "zggbal.f"
    k = 1;
#line 278 "zggbal.f"
    l = *n;
#line 279 "zggbal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 279 "zggbal.f"
	goto L190;
#line 279 "zggbal.f"
    }

#line 282 "zggbal.f"
    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues. */

/*     Find row with one nonzero in columns 1 through L */

#line 288 "zggbal.f"
L20:
#line 289 "zggbal.f"
    l = lm1;
#line 290 "zggbal.f"
    if (l != 1) {
#line 290 "zggbal.f"
	goto L30;
#line 290 "zggbal.f"
    }

#line 293 "zggbal.f"
    rscale[1] = 1.;
#line 294 "zggbal.f"
    lscale[1] = 1.;
#line 295 "zggbal.f"
    goto L190;

#line 297 "zggbal.f"
L30:
#line 298 "zggbal.f"
    lm1 = l - 1;
#line 299 "zggbal.f"
    for (i__ = l; i__ >= 1; --i__) {
#line 300 "zggbal.f"
	i__1 = lm1;
#line 300 "zggbal.f"
	for (j = 1; j <= i__1; ++j) {
#line 301 "zggbal.f"
	    jp1 = j + 1;
#line 302 "zggbal.f"
	    i__2 = i__ + j * a_dim1;
#line 302 "zggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 302 "zggbal.f"
	    if (a[i__2].r != 0. || a[i__2].i != 0. || (b[i__3].r != 0. || b[
		    i__3].i != 0.)) {
#line 302 "zggbal.f"
		goto L50;
#line 302 "zggbal.f"
	    }
#line 304 "zggbal.f"
/* L40: */
#line 304 "zggbal.f"
	}
#line 305 "zggbal.f"
	j = l;
#line 306 "zggbal.f"
	goto L70;

#line 308 "zggbal.f"
L50:
#line 309 "zggbal.f"
	i__1 = l;
#line 309 "zggbal.f"
	for (j = jp1; j <= i__1; ++j) {
#line 310 "zggbal.f"
	    i__2 = i__ + j * a_dim1;
#line 310 "zggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 310 "zggbal.f"
	    if (a[i__2].r != 0. || a[i__2].i != 0. || (b[i__3].r != 0. || b[
		    i__3].i != 0.)) {
#line 310 "zggbal.f"
		goto L80;
#line 310 "zggbal.f"
	    }
#line 312 "zggbal.f"
/* L60: */
#line 312 "zggbal.f"
	}
#line 313 "zggbal.f"
	j = jp1 - 1;

#line 315 "zggbal.f"
L70:
#line 316 "zggbal.f"
	m = l;
#line 317 "zggbal.f"
	iflow = 1;
#line 318 "zggbal.f"
	goto L160;
#line 319 "zggbal.f"
L80:
#line 319 "zggbal.f"
	;
#line 319 "zggbal.f"
    }
#line 320 "zggbal.f"
    goto L100;

/*     Find column with one nonzero in rows K through N */

#line 324 "zggbal.f"
L90:
#line 325 "zggbal.f"
    ++k;

#line 327 "zggbal.f"
L100:
#line 328 "zggbal.f"
    i__1 = l;
#line 328 "zggbal.f"
    for (j = k; j <= i__1; ++j) {
#line 329 "zggbal.f"
	i__2 = lm1;
#line 329 "zggbal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 330 "zggbal.f"
	    ip1 = i__ + 1;
#line 331 "zggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 331 "zggbal.f"
	    i__4 = i__ + j * b_dim1;
#line 331 "zggbal.f"
	    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 0. || b[
		    i__4].i != 0.)) {
#line 331 "zggbal.f"
		goto L120;
#line 331 "zggbal.f"
	    }
#line 333 "zggbal.f"
/* L110: */
#line 333 "zggbal.f"
	}
#line 334 "zggbal.f"
	i__ = l;
#line 335 "zggbal.f"
	goto L140;
#line 336 "zggbal.f"
L120:
#line 337 "zggbal.f"
	i__2 = l;
#line 337 "zggbal.f"
	for (i__ = ip1; i__ <= i__2; ++i__) {
#line 338 "zggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 338 "zggbal.f"
	    i__4 = i__ + j * b_dim1;
#line 338 "zggbal.f"
	    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 0. || b[
		    i__4].i != 0.)) {
#line 338 "zggbal.f"
		goto L150;
#line 338 "zggbal.f"
	    }
#line 340 "zggbal.f"
/* L130: */
#line 340 "zggbal.f"
	}
#line 341 "zggbal.f"
	i__ = ip1 - 1;
#line 342 "zggbal.f"
L140:
#line 343 "zggbal.f"
	m = k;
#line 344 "zggbal.f"
	iflow = 2;
#line 345 "zggbal.f"
	goto L160;
#line 346 "zggbal.f"
L150:
#line 346 "zggbal.f"
	;
#line 346 "zggbal.f"
    }
#line 347 "zggbal.f"
    goto L190;

/*     Permute rows M and I */

#line 351 "zggbal.f"
L160:
#line 352 "zggbal.f"
    lscale[m] = (doublereal) i__;
#line 353 "zggbal.f"
    if (i__ == m) {
#line 353 "zggbal.f"
	goto L170;
#line 353 "zggbal.f"
    }
#line 355 "zggbal.f"
    i__1 = *n - k + 1;
#line 355 "zggbal.f"
    zswap_(&i__1, &a[i__ + k * a_dim1], lda, &a[m + k * a_dim1], lda);
#line 356 "zggbal.f"
    i__1 = *n - k + 1;
#line 356 "zggbal.f"
    zswap_(&i__1, &b[i__ + k * b_dim1], ldb, &b[m + k * b_dim1], ldb);

/*     Permute columns M and J */

#line 360 "zggbal.f"
L170:
#line 361 "zggbal.f"
    rscale[m] = (doublereal) j;
#line 362 "zggbal.f"
    if (j == m) {
#line 362 "zggbal.f"
	goto L180;
#line 362 "zggbal.f"
    }
#line 364 "zggbal.f"
    zswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 365 "zggbal.f"
    zswap_(&l, &b[j * b_dim1 + 1], &c__1, &b[m * b_dim1 + 1], &c__1);

#line 367 "zggbal.f"
L180:
#line 368 "zggbal.f"
    switch (iflow) {
#line 368 "zggbal.f"
	case 1:  goto L20;
#line 368 "zggbal.f"
	case 2:  goto L90;
#line 368 "zggbal.f"
    }

#line 370 "zggbal.f"
L190:
#line 371 "zggbal.f"
    *ilo = k;
#line 372 "zggbal.f"
    *ihi = l;

#line 374 "zggbal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 375 "zggbal.f"
	i__1 = *ihi;
#line 375 "zggbal.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 376 "zggbal.f"
	    lscale[i__] = 1.;
#line 377 "zggbal.f"
	    rscale[i__] = 1.;
#line 378 "zggbal.f"
/* L195: */
#line 378 "zggbal.f"
	}
#line 379 "zggbal.f"
	return 0;
#line 380 "zggbal.f"
    }

#line 382 "zggbal.f"
    if (*ilo == *ihi) {
#line 382 "zggbal.f"
	return 0;
#line 382 "zggbal.f"
    }

/*     Balance the submatrix in rows ILO to IHI. */

#line 387 "zggbal.f"
    nr = *ihi - *ilo + 1;
#line 388 "zggbal.f"
    i__1 = *ihi;
#line 388 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 389 "zggbal.f"
	rscale[i__] = 0.;
#line 390 "zggbal.f"
	lscale[i__] = 0.;

#line 392 "zggbal.f"
	work[i__] = 0.;
#line 393 "zggbal.f"
	work[i__ + *n] = 0.;
#line 394 "zggbal.f"
	work[i__ + (*n << 1)] = 0.;
#line 395 "zggbal.f"
	work[i__ + *n * 3] = 0.;
#line 396 "zggbal.f"
	work[i__ + (*n << 2)] = 0.;
#line 397 "zggbal.f"
	work[i__ + *n * 5] = 0.;
#line 398 "zggbal.f"
/* L200: */
#line 398 "zggbal.f"
    }

/*     Compute right side vector in resulting linear equations */

#line 402 "zggbal.f"
    basl = d_lg10(&c_b36);
#line 403 "zggbal.f"
    i__1 = *ihi;
#line 403 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 404 "zggbal.f"
	i__2 = *ihi;
#line 404 "zggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 405 "zggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 405 "zggbal.f"
	    if (a[i__3].r == 0. && a[i__3].i == 0.) {
#line 406 "zggbal.f"
		ta = 0.;
#line 407 "zggbal.f"
		goto L210;
#line 408 "zggbal.f"
	    }
#line 409 "zggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 409 "zggbal.f"
	    d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j *
		     a_dim1]), abs(d__2));
#line 409 "zggbal.f"
	    ta = d_lg10(&d__3) / basl;

#line 411 "zggbal.f"
L210:
#line 412 "zggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 412 "zggbal.f"
	    if (b[i__3].r == 0. && b[i__3].i == 0.) {
#line 413 "zggbal.f"
		tb = 0.;
#line 414 "zggbal.f"
		goto L220;
#line 415 "zggbal.f"
	    }
#line 416 "zggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 416 "zggbal.f"
	    d__3 = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + j *
		     b_dim1]), abs(d__2));
#line 416 "zggbal.f"
	    tb = d_lg10(&d__3) / basl;

#line 418 "zggbal.f"
L220:
#line 419 "zggbal.f"
	    work[i__ + (*n << 2)] = work[i__ + (*n << 2)] - ta - tb;
#line 420 "zggbal.f"
	    work[j + *n * 5] = work[j + *n * 5] - ta - tb;
#line 421 "zggbal.f"
/* L230: */
#line 421 "zggbal.f"
	}
#line 422 "zggbal.f"
/* L240: */
#line 422 "zggbal.f"
    }

#line 424 "zggbal.f"
    coef = 1. / (doublereal) (nr << 1);
#line 425 "zggbal.f"
    coef2 = coef * coef;
#line 426 "zggbal.f"
    coef5 = coef2 * .5;
#line 427 "zggbal.f"
    nrp2 = nr + 2;
#line 428 "zggbal.f"
    beta = 0.;
#line 429 "zggbal.f"
    it = 1;

/*     Start generalized conjugate gradient iteration */

#line 433 "zggbal.f"
L250:

#line 435 "zggbal.f"
    gamma = ddot_(&nr, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1) + ddot_(&nr, &work[*ilo + *n * 5], &c__1, &work[*ilo + *
	    n * 5], &c__1);

#line 438 "zggbal.f"
    ew = 0.;
#line 439 "zggbal.f"
    ewc = 0.;
#line 440 "zggbal.f"
    i__1 = *ihi;
#line 440 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 441 "zggbal.f"
	ew += work[i__ + (*n << 2)];
#line 442 "zggbal.f"
	ewc += work[i__ + *n * 5];
#line 443 "zggbal.f"
/* L260: */
#line 443 "zggbal.f"
    }

/* Computing 2nd power */
#line 445 "zggbal.f"
    d__1 = ew;
/* Computing 2nd power */
#line 445 "zggbal.f"
    d__2 = ewc;
/* Computing 2nd power */
#line 445 "zggbal.f"
    d__3 = ew - ewc;
#line 445 "zggbal.f"
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
#line 446 "zggbal.f"
    if (gamma == 0.) {
#line 446 "zggbal.f"
	goto L350;
#line 446 "zggbal.f"
    }
#line 448 "zggbal.f"
    if (it != 1) {
#line 448 "zggbal.f"
	beta = gamma / pgamma;
#line 448 "zggbal.f"
    }
#line 450 "zggbal.f"
    t = coef5 * (ewc - ew * 3.);
#line 451 "zggbal.f"
    tc = coef5 * (ew - ewc * 3.);

#line 453 "zggbal.f"
    dscal_(&nr, &beta, &work[*ilo], &c__1);
#line 454 "zggbal.f"
    dscal_(&nr, &beta, &work[*ilo + *n], &c__1);

#line 456 "zggbal.f"
    daxpy_(&nr, &coef, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + *n], &
	    c__1);
#line 457 "zggbal.f"
    daxpy_(&nr, &coef, &work[*ilo + *n * 5], &c__1, &work[*ilo], &c__1);

#line 459 "zggbal.f"
    i__1 = *ihi;
#line 459 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 460 "zggbal.f"
	work[i__] += tc;
#line 461 "zggbal.f"
	work[i__ + *n] += t;
#line 462 "zggbal.f"
/* L270: */
#line 462 "zggbal.f"
    }

/*     Apply matrix to vector */

#line 466 "zggbal.f"
    i__1 = *ihi;
#line 466 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 467 "zggbal.f"
	kount = 0;
#line 468 "zggbal.f"
	sum = 0.;
#line 469 "zggbal.f"
	i__2 = *ihi;
#line 469 "zggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 470 "zggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 470 "zggbal.f"
	    if (a[i__3].r == 0. && a[i__3].i == 0.) {
#line 470 "zggbal.f"
		goto L280;
#line 470 "zggbal.f"
	    }
#line 472 "zggbal.f"
	    ++kount;
#line 473 "zggbal.f"
	    sum += work[j];
#line 474 "zggbal.f"
L280:
#line 475 "zggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 475 "zggbal.f"
	    if (b[i__3].r == 0. && b[i__3].i == 0.) {
#line 475 "zggbal.f"
		goto L290;
#line 475 "zggbal.f"
	    }
#line 477 "zggbal.f"
	    ++kount;
#line 478 "zggbal.f"
	    sum += work[j];
#line 479 "zggbal.f"
L290:
#line 479 "zggbal.f"
	    ;
#line 479 "zggbal.f"
	}
#line 480 "zggbal.f"
	work[i__ + (*n << 1)] = (doublereal) kount * work[i__ + *n] + sum;
#line 481 "zggbal.f"
/* L300: */
#line 481 "zggbal.f"
    }

#line 483 "zggbal.f"
    i__1 = *ihi;
#line 483 "zggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 484 "zggbal.f"
	kount = 0;
#line 485 "zggbal.f"
	sum = 0.;
#line 486 "zggbal.f"
	i__2 = *ihi;
#line 486 "zggbal.f"
	for (i__ = *ilo; i__ <= i__2; ++i__) {
#line 487 "zggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 487 "zggbal.f"
	    if (a[i__3].r == 0. && a[i__3].i == 0.) {
#line 487 "zggbal.f"
		goto L310;
#line 487 "zggbal.f"
	    }
#line 489 "zggbal.f"
	    ++kount;
#line 490 "zggbal.f"
	    sum += work[i__ + *n];
#line 491 "zggbal.f"
L310:
#line 492 "zggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 492 "zggbal.f"
	    if (b[i__3].r == 0. && b[i__3].i == 0.) {
#line 492 "zggbal.f"
		goto L320;
#line 492 "zggbal.f"
	    }
#line 494 "zggbal.f"
	    ++kount;
#line 495 "zggbal.f"
	    sum += work[i__ + *n];
#line 496 "zggbal.f"
L320:
#line 496 "zggbal.f"
	    ;
#line 496 "zggbal.f"
	}
#line 497 "zggbal.f"
	work[j + *n * 3] = (doublereal) kount * work[j] + sum;
#line 498 "zggbal.f"
/* L330: */
#line 498 "zggbal.f"
    }

#line 500 "zggbal.f"
    sum = ddot_(&nr, &work[*ilo + *n], &c__1, &work[*ilo + (*n << 1)], &c__1) 
	    + ddot_(&nr, &work[*ilo], &c__1, &work[*ilo + *n * 3], &c__1);
#line 502 "zggbal.f"
    alpha = gamma / sum;

/*     Determine correction to current iteration */

#line 506 "zggbal.f"
    cmax = 0.;
#line 507 "zggbal.f"
    i__1 = *ihi;
#line 507 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 508 "zggbal.f"
	cor = alpha * work[i__ + *n];
#line 509 "zggbal.f"
	if (abs(cor) > cmax) {
#line 509 "zggbal.f"
	    cmax = abs(cor);
#line 509 "zggbal.f"
	}
#line 511 "zggbal.f"
	lscale[i__] += cor;
#line 512 "zggbal.f"
	cor = alpha * work[i__];
#line 513 "zggbal.f"
	if (abs(cor) > cmax) {
#line 513 "zggbal.f"
	    cmax = abs(cor);
#line 513 "zggbal.f"
	}
#line 515 "zggbal.f"
	rscale[i__] += cor;
#line 516 "zggbal.f"
/* L340: */
#line 516 "zggbal.f"
    }
#line 517 "zggbal.f"
    if (cmax < .5) {
#line 517 "zggbal.f"
	goto L350;
#line 517 "zggbal.f"
    }

#line 520 "zggbal.f"
    d__1 = -alpha;
#line 520 "zggbal.f"
    daxpy_(&nr, &d__1, &work[*ilo + (*n << 1)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1);
#line 521 "zggbal.f"
    d__1 = -alpha;
#line 521 "zggbal.f"
    daxpy_(&nr, &d__1, &work[*ilo + *n * 3], &c__1, &work[*ilo + *n * 5], &
	    c__1);

#line 523 "zggbal.f"
    pgamma = gamma;
#line 524 "zggbal.f"
    ++it;
#line 525 "zggbal.f"
    if (it <= nrp2) {
#line 525 "zggbal.f"
	goto L250;
#line 525 "zggbal.f"
    }

/*     End generalized conjugate gradient iteration */

#line 530 "zggbal.f"
L350:
#line 531 "zggbal.f"
    sfmin = dlamch_("S", (ftnlen)1);
#line 532 "zggbal.f"
    sfmax = 1. / sfmin;
#line 533 "zggbal.f"
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
#line 534 "zggbal.f"
    lsfmax = (integer) (d_lg10(&sfmax) / basl);
#line 535 "zggbal.f"
    i__1 = *ihi;
#line 535 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 536 "zggbal.f"
	i__2 = *n - *ilo + 1;
#line 536 "zggbal.f"
	irab = izamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
#line 537 "zggbal.f"
	rab = z_abs(&a[i__ + (irab + *ilo - 1) * a_dim1]);
#line 538 "zggbal.f"
	i__2 = *n - *ilo + 1;
#line 538 "zggbal.f"
	irab = izamax_(&i__2, &b[i__ + *ilo * b_dim1], ldb);
/* Computing MAX */
#line 539 "zggbal.f"
	d__1 = rab, d__2 = z_abs(&b[i__ + (irab + *ilo - 1) * b_dim1]);
#line 539 "zggbal.f"
	rab = max(d__1,d__2);
#line 540 "zggbal.f"
	d__1 = rab + sfmin;
#line 540 "zggbal.f"
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 541 "zggbal.f"
	ir = (integer) (lscale[i__] + d_sign(&c_b72, &lscale[i__]));
/* Computing MIN */
#line 542 "zggbal.f"
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
#line 542 "zggbal.f"
	ir = min(i__2,i__3);
#line 543 "zggbal.f"
	lscale[i__] = pow_di(&c_b36, &ir);
#line 544 "zggbal.f"
	icab = izamax_(ihi, &a[i__ * a_dim1 + 1], &c__1);
#line 545 "zggbal.f"
	cab = z_abs(&a[icab + i__ * a_dim1]);
#line 546 "zggbal.f"
	icab = izamax_(ihi, &b[i__ * b_dim1 + 1], &c__1);
/* Computing MAX */
#line 547 "zggbal.f"
	d__1 = cab, d__2 = z_abs(&b[icab + i__ * b_dim1]);
#line 547 "zggbal.f"
	cab = max(d__1,d__2);
#line 548 "zggbal.f"
	d__1 = cab + sfmin;
#line 548 "zggbal.f"
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 549 "zggbal.f"
	jc = (integer) (rscale[i__] + d_sign(&c_b72, &rscale[i__]));
/* Computing MIN */
#line 550 "zggbal.f"
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
#line 550 "zggbal.f"
	jc = min(i__2,i__3);
#line 551 "zggbal.f"
	rscale[i__] = pow_di(&c_b36, &jc);
#line 552 "zggbal.f"
/* L360: */
#line 552 "zggbal.f"
    }

/*     Row scaling of matrices A and B */

#line 556 "zggbal.f"
    i__1 = *ihi;
#line 556 "zggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 557 "zggbal.f"
	i__2 = *n - *ilo + 1;
#line 557 "zggbal.f"
	zdscal_(&i__2, &lscale[i__], &a[i__ + *ilo * a_dim1], lda);
#line 558 "zggbal.f"
	i__2 = *n - *ilo + 1;
#line 558 "zggbal.f"
	zdscal_(&i__2, &lscale[i__], &b[i__ + *ilo * b_dim1], ldb);
#line 559 "zggbal.f"
/* L370: */
#line 559 "zggbal.f"
    }

/*     Column scaling of matrices A and B */

#line 563 "zggbal.f"
    i__1 = *ihi;
#line 563 "zggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 564 "zggbal.f"
	zdscal_(ihi, &rscale[j], &a[j * a_dim1 + 1], &c__1);
#line 565 "zggbal.f"
	zdscal_(ihi, &rscale[j], &b[j * b_dim1 + 1], &c__1);
#line 566 "zggbal.f"
/* L380: */
#line 566 "zggbal.f"
    }

#line 568 "zggbal.f"
    return 0;

/*     End of ZGGBAL */

} /* zggbal_ */

