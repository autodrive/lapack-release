#line 1 "cggbal.f"
/* cggbal.f -- translated by f2c (version 20100827).
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

#line 1 "cggbal.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = 10.;
static doublereal c_b72 = .5;

/* > \brief \b CGGBAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggbal.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggbal.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggbal.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, */
/*                          RSCALE, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               LSCALE( * ), RSCALE( * ), WORK( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGBAL balances a pair of general complex matrices (A,B).  This */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          LSCALE is REAL array, dimension (N) */
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
/* >          RSCALE is REAL array, dimension (N) */
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

/* > \date December 2016 */

/* > \ingroup complexGBcomputational */

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
/* Subroutine */ int cggbal_(char *job, integer *n, doublecomplex *a, integer 
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
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal coef2, coef5, gamma, alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sfmin;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal sfmax;
    static integer iflow, kount;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal pgamma;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 235 "cggbal.f"
    /* Parameter adjustments */
#line 235 "cggbal.f"
    a_dim1 = *lda;
#line 235 "cggbal.f"
    a_offset = 1 + a_dim1;
#line 235 "cggbal.f"
    a -= a_offset;
#line 235 "cggbal.f"
    b_dim1 = *ldb;
#line 235 "cggbal.f"
    b_offset = 1 + b_dim1;
#line 235 "cggbal.f"
    b -= b_offset;
#line 235 "cggbal.f"
    --lscale;
#line 235 "cggbal.f"
    --rscale;
#line 235 "cggbal.f"
    --work;
#line 235 "cggbal.f"

#line 235 "cggbal.f"
    /* Function Body */
#line 235 "cggbal.f"
    *info = 0;
#line 236 "cggbal.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 238 "cggbal.f"
	*info = -1;
#line 239 "cggbal.f"
    } else if (*n < 0) {
#line 240 "cggbal.f"
	*info = -2;
#line 241 "cggbal.f"
    } else if (*lda < max(1,*n)) {
#line 242 "cggbal.f"
	*info = -4;
#line 243 "cggbal.f"
    } else if (*ldb < max(1,*n)) {
#line 244 "cggbal.f"
	*info = -6;
#line 245 "cggbal.f"
    }
#line 246 "cggbal.f"
    if (*info != 0) {
#line 247 "cggbal.f"
	i__1 = -(*info);
#line 247 "cggbal.f"
	xerbla_("CGGBAL", &i__1, (ftnlen)6);
#line 248 "cggbal.f"
	return 0;
#line 249 "cggbal.f"
    }

/*     Quick return if possible */

#line 253 "cggbal.f"
    if (*n == 0) {
#line 254 "cggbal.f"
	*ilo = 1;
#line 255 "cggbal.f"
	*ihi = *n;
#line 256 "cggbal.f"
	return 0;
#line 257 "cggbal.f"
    }

#line 259 "cggbal.f"
    if (*n == 1) {
#line 260 "cggbal.f"
	*ilo = 1;
#line 261 "cggbal.f"
	*ihi = *n;
#line 262 "cggbal.f"
	lscale[1] = 1.;
#line 263 "cggbal.f"
	rscale[1] = 1.;
#line 264 "cggbal.f"
	return 0;
#line 265 "cggbal.f"
    }

#line 267 "cggbal.f"
    if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 268 "cggbal.f"
	*ilo = 1;
#line 269 "cggbal.f"
	*ihi = *n;
#line 270 "cggbal.f"
	i__1 = *n;
#line 270 "cggbal.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "cggbal.f"
	    lscale[i__] = 1.;
#line 272 "cggbal.f"
	    rscale[i__] = 1.;
#line 273 "cggbal.f"
/* L10: */
#line 273 "cggbal.f"
	}
#line 274 "cggbal.f"
	return 0;
#line 275 "cggbal.f"
    }

#line 277 "cggbal.f"
    k = 1;
#line 278 "cggbal.f"
    l = *n;
#line 279 "cggbal.f"
    if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 279 "cggbal.f"
	goto L190;
#line 279 "cggbal.f"
    }

#line 282 "cggbal.f"
    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues. */

/*     Find row with one nonzero in columns 1 through L */

#line 288 "cggbal.f"
L20:
#line 289 "cggbal.f"
    l = lm1;
#line 290 "cggbal.f"
    if (l != 1) {
#line 290 "cggbal.f"
	goto L30;
#line 290 "cggbal.f"
    }

#line 293 "cggbal.f"
    rscale[1] = 1.;
#line 294 "cggbal.f"
    lscale[1] = 1.;
#line 295 "cggbal.f"
    goto L190;

#line 297 "cggbal.f"
L30:
#line 298 "cggbal.f"
    lm1 = l - 1;
#line 299 "cggbal.f"
    for (i__ = l; i__ >= 1; --i__) {
#line 300 "cggbal.f"
	i__1 = lm1;
#line 300 "cggbal.f"
	for (j = 1; j <= i__1; ++j) {
#line 301 "cggbal.f"
	    jp1 = j + 1;
#line 302 "cggbal.f"
	    i__2 = i__ + j * a_dim1;
#line 302 "cggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 302 "cggbal.f"
	    if (a[i__2].r != 0. || a[i__2].i != 0. || (b[i__3].r != 0. || b[
		    i__3].i != 0.)) {
#line 302 "cggbal.f"
		goto L50;
#line 302 "cggbal.f"
	    }
#line 304 "cggbal.f"
/* L40: */
#line 304 "cggbal.f"
	}
#line 305 "cggbal.f"
	j = l;
#line 306 "cggbal.f"
	goto L70;

#line 308 "cggbal.f"
L50:
#line 309 "cggbal.f"
	i__1 = l;
#line 309 "cggbal.f"
	for (j = jp1; j <= i__1; ++j) {
#line 310 "cggbal.f"
	    i__2 = i__ + j * a_dim1;
#line 310 "cggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 310 "cggbal.f"
	    if (a[i__2].r != 0. || a[i__2].i != 0. || (b[i__3].r != 0. || b[
		    i__3].i != 0.)) {
#line 310 "cggbal.f"
		goto L80;
#line 310 "cggbal.f"
	    }
#line 312 "cggbal.f"
/* L60: */
#line 312 "cggbal.f"
	}
#line 313 "cggbal.f"
	j = jp1 - 1;

#line 315 "cggbal.f"
L70:
#line 316 "cggbal.f"
	m = l;
#line 317 "cggbal.f"
	iflow = 1;
#line 318 "cggbal.f"
	goto L160;
#line 319 "cggbal.f"
L80:
#line 319 "cggbal.f"
	;
#line 319 "cggbal.f"
    }
#line 320 "cggbal.f"
    goto L100;

/*     Find column with one nonzero in rows K through N */

#line 324 "cggbal.f"
L90:
#line 325 "cggbal.f"
    ++k;

#line 327 "cggbal.f"
L100:
#line 328 "cggbal.f"
    i__1 = l;
#line 328 "cggbal.f"
    for (j = k; j <= i__1; ++j) {
#line 329 "cggbal.f"
	i__2 = lm1;
#line 329 "cggbal.f"
	for (i__ = k; i__ <= i__2; ++i__) {
#line 330 "cggbal.f"
	    ip1 = i__ + 1;
#line 331 "cggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 331 "cggbal.f"
	    i__4 = i__ + j * b_dim1;
#line 331 "cggbal.f"
	    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 0. || b[
		    i__4].i != 0.)) {
#line 331 "cggbal.f"
		goto L120;
#line 331 "cggbal.f"
	    }
#line 333 "cggbal.f"
/* L110: */
#line 333 "cggbal.f"
	}
#line 334 "cggbal.f"
	i__ = l;
#line 335 "cggbal.f"
	goto L140;
#line 336 "cggbal.f"
L120:
#line 337 "cggbal.f"
	i__2 = l;
#line 337 "cggbal.f"
	for (i__ = ip1; i__ <= i__2; ++i__) {
#line 338 "cggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 338 "cggbal.f"
	    i__4 = i__ + j * b_dim1;
#line 338 "cggbal.f"
	    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 0. || b[
		    i__4].i != 0.)) {
#line 338 "cggbal.f"
		goto L150;
#line 338 "cggbal.f"
	    }
#line 340 "cggbal.f"
/* L130: */
#line 340 "cggbal.f"
	}
#line 341 "cggbal.f"
	i__ = ip1 - 1;
#line 342 "cggbal.f"
L140:
#line 343 "cggbal.f"
	m = k;
#line 344 "cggbal.f"
	iflow = 2;
#line 345 "cggbal.f"
	goto L160;
#line 346 "cggbal.f"
L150:
#line 346 "cggbal.f"
	;
#line 346 "cggbal.f"
    }
#line 347 "cggbal.f"
    goto L190;

/*     Permute rows M and I */

#line 351 "cggbal.f"
L160:
#line 352 "cggbal.f"
    lscale[m] = (doublereal) i__;
#line 353 "cggbal.f"
    if (i__ == m) {
#line 353 "cggbal.f"
	goto L170;
#line 353 "cggbal.f"
    }
#line 355 "cggbal.f"
    i__1 = *n - k + 1;
#line 355 "cggbal.f"
    cswap_(&i__1, &a[i__ + k * a_dim1], lda, &a[m + k * a_dim1], lda);
#line 356 "cggbal.f"
    i__1 = *n - k + 1;
#line 356 "cggbal.f"
    cswap_(&i__1, &b[i__ + k * b_dim1], ldb, &b[m + k * b_dim1], ldb);

/*     Permute columns M and J */

#line 360 "cggbal.f"
L170:
#line 361 "cggbal.f"
    rscale[m] = (doublereal) j;
#line 362 "cggbal.f"
    if (j == m) {
#line 362 "cggbal.f"
	goto L180;
#line 362 "cggbal.f"
    }
#line 364 "cggbal.f"
    cswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
#line 365 "cggbal.f"
    cswap_(&l, &b[j * b_dim1 + 1], &c__1, &b[m * b_dim1 + 1], &c__1);

#line 367 "cggbal.f"
L180:
#line 368 "cggbal.f"
    switch (iflow) {
#line 368 "cggbal.f"
	case 1:  goto L20;
#line 368 "cggbal.f"
	case 2:  goto L90;
#line 368 "cggbal.f"
    }

#line 370 "cggbal.f"
L190:
#line 371 "cggbal.f"
    *ilo = k;
#line 372 "cggbal.f"
    *ihi = l;

#line 374 "cggbal.f"
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
#line 375 "cggbal.f"
	i__1 = *ihi;
#line 375 "cggbal.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 376 "cggbal.f"
	    lscale[i__] = 1.;
#line 377 "cggbal.f"
	    rscale[i__] = 1.;
#line 378 "cggbal.f"
/* L195: */
#line 378 "cggbal.f"
	}
#line 379 "cggbal.f"
	return 0;
#line 380 "cggbal.f"
    }

#line 382 "cggbal.f"
    if (*ilo == *ihi) {
#line 382 "cggbal.f"
	return 0;
#line 382 "cggbal.f"
    }

/*     Balance the submatrix in rows ILO to IHI. */

#line 387 "cggbal.f"
    nr = *ihi - *ilo + 1;
#line 388 "cggbal.f"
    i__1 = *ihi;
#line 388 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 389 "cggbal.f"
	rscale[i__] = 0.;
#line 390 "cggbal.f"
	lscale[i__] = 0.;

#line 392 "cggbal.f"
	work[i__] = 0.;
#line 393 "cggbal.f"
	work[i__ + *n] = 0.;
#line 394 "cggbal.f"
	work[i__ + (*n << 1)] = 0.;
#line 395 "cggbal.f"
	work[i__ + *n * 3] = 0.;
#line 396 "cggbal.f"
	work[i__ + (*n << 2)] = 0.;
#line 397 "cggbal.f"
	work[i__ + *n * 5] = 0.;
#line 398 "cggbal.f"
/* L200: */
#line 398 "cggbal.f"
    }

/*     Compute right side vector in resulting linear equations */

#line 402 "cggbal.f"
    basl = d_lg10(&c_b36);
#line 403 "cggbal.f"
    i__1 = *ihi;
#line 403 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 404 "cggbal.f"
	i__2 = *ihi;
#line 404 "cggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 405 "cggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 405 "cggbal.f"
	    if (a[i__3].r == 0. && a[i__3].i == 0.) {
#line 406 "cggbal.f"
		ta = 0.;
#line 407 "cggbal.f"
		goto L210;
#line 408 "cggbal.f"
	    }
#line 409 "cggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 409 "cggbal.f"
	    d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j *
		     a_dim1]), abs(d__2));
#line 409 "cggbal.f"
	    ta = d_lg10(&d__3) / basl;

#line 411 "cggbal.f"
L210:
#line 412 "cggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 412 "cggbal.f"
	    if (b[i__3].r == 0. && b[i__3].i == 0.) {
#line 413 "cggbal.f"
		tb = 0.;
#line 414 "cggbal.f"
		goto L220;
#line 415 "cggbal.f"
	    }
#line 416 "cggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 416 "cggbal.f"
	    d__3 = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + j *
		     b_dim1]), abs(d__2));
#line 416 "cggbal.f"
	    tb = d_lg10(&d__3) / basl;

#line 418 "cggbal.f"
L220:
#line 419 "cggbal.f"
	    work[i__ + (*n << 2)] = work[i__ + (*n << 2)] - ta - tb;
#line 420 "cggbal.f"
	    work[j + *n * 5] = work[j + *n * 5] - ta - tb;
#line 421 "cggbal.f"
/* L230: */
#line 421 "cggbal.f"
	}
#line 422 "cggbal.f"
/* L240: */
#line 422 "cggbal.f"
    }

#line 424 "cggbal.f"
    coef = 1. / (doublereal) (nr << 1);
#line 425 "cggbal.f"
    coef2 = coef * coef;
#line 426 "cggbal.f"
    coef5 = coef2 * .5;
#line 427 "cggbal.f"
    nrp2 = nr + 2;
#line 428 "cggbal.f"
    beta = 0.;
#line 429 "cggbal.f"
    it = 1;

/*     Start generalized conjugate gradient iteration */

#line 433 "cggbal.f"
L250:

#line 435 "cggbal.f"
    gamma = sdot_(&nr, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1) + sdot_(&nr, &work[*ilo + *n * 5], &c__1, &work[*ilo + *
	    n * 5], &c__1);

#line 438 "cggbal.f"
    ew = 0.;
#line 439 "cggbal.f"
    ewc = 0.;
#line 440 "cggbal.f"
    i__1 = *ihi;
#line 440 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 441 "cggbal.f"
	ew += work[i__ + (*n << 2)];
#line 442 "cggbal.f"
	ewc += work[i__ + *n * 5];
#line 443 "cggbal.f"
/* L260: */
#line 443 "cggbal.f"
    }

/* Computing 2nd power */
#line 445 "cggbal.f"
    d__1 = ew;
/* Computing 2nd power */
#line 445 "cggbal.f"
    d__2 = ewc;
/* Computing 2nd power */
#line 445 "cggbal.f"
    d__3 = ew - ewc;
#line 445 "cggbal.f"
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
#line 446 "cggbal.f"
    if (gamma == 0.) {
#line 446 "cggbal.f"
	goto L350;
#line 446 "cggbal.f"
    }
#line 448 "cggbal.f"
    if (it != 1) {
#line 448 "cggbal.f"
	beta = gamma / pgamma;
#line 448 "cggbal.f"
    }
#line 450 "cggbal.f"
    t = coef5 * (ewc - ew * 3.);
#line 451 "cggbal.f"
    tc = coef5 * (ew - ewc * 3.);

#line 453 "cggbal.f"
    sscal_(&nr, &beta, &work[*ilo], &c__1);
#line 454 "cggbal.f"
    sscal_(&nr, &beta, &work[*ilo + *n], &c__1);

#line 456 "cggbal.f"
    saxpy_(&nr, &coef, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + *n], &
	    c__1);
#line 457 "cggbal.f"
    saxpy_(&nr, &coef, &work[*ilo + *n * 5], &c__1, &work[*ilo], &c__1);

#line 459 "cggbal.f"
    i__1 = *ihi;
#line 459 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 460 "cggbal.f"
	work[i__] += tc;
#line 461 "cggbal.f"
	work[i__ + *n] += t;
#line 462 "cggbal.f"
/* L270: */
#line 462 "cggbal.f"
    }

/*     Apply matrix to vector */

#line 466 "cggbal.f"
    i__1 = *ihi;
#line 466 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 467 "cggbal.f"
	kount = 0;
#line 468 "cggbal.f"
	sum = 0.;
#line 469 "cggbal.f"
	i__2 = *ihi;
#line 469 "cggbal.f"
	for (j = *ilo; j <= i__2; ++j) {
#line 470 "cggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 470 "cggbal.f"
	    if (a[i__3].r == 0. && a[i__3].i == 0.) {
#line 470 "cggbal.f"
		goto L280;
#line 470 "cggbal.f"
	    }
#line 472 "cggbal.f"
	    ++kount;
#line 473 "cggbal.f"
	    sum += work[j];
#line 474 "cggbal.f"
L280:
#line 475 "cggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 475 "cggbal.f"
	    if (b[i__3].r == 0. && b[i__3].i == 0.) {
#line 475 "cggbal.f"
		goto L290;
#line 475 "cggbal.f"
	    }
#line 477 "cggbal.f"
	    ++kount;
#line 478 "cggbal.f"
	    sum += work[j];
#line 479 "cggbal.f"
L290:
#line 479 "cggbal.f"
	    ;
#line 479 "cggbal.f"
	}
#line 480 "cggbal.f"
	work[i__ + (*n << 1)] = (doublereal) kount * work[i__ + *n] + sum;
#line 481 "cggbal.f"
/* L300: */
#line 481 "cggbal.f"
    }

#line 483 "cggbal.f"
    i__1 = *ihi;
#line 483 "cggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 484 "cggbal.f"
	kount = 0;
#line 485 "cggbal.f"
	sum = 0.;
#line 486 "cggbal.f"
	i__2 = *ihi;
#line 486 "cggbal.f"
	for (i__ = *ilo; i__ <= i__2; ++i__) {
#line 487 "cggbal.f"
	    i__3 = i__ + j * a_dim1;
#line 487 "cggbal.f"
	    if (a[i__3].r == 0. && a[i__3].i == 0.) {
#line 487 "cggbal.f"
		goto L310;
#line 487 "cggbal.f"
	    }
#line 489 "cggbal.f"
	    ++kount;
#line 490 "cggbal.f"
	    sum += work[i__ + *n];
#line 491 "cggbal.f"
L310:
#line 492 "cggbal.f"
	    i__3 = i__ + j * b_dim1;
#line 492 "cggbal.f"
	    if (b[i__3].r == 0. && b[i__3].i == 0.) {
#line 492 "cggbal.f"
		goto L320;
#line 492 "cggbal.f"
	    }
#line 494 "cggbal.f"
	    ++kount;
#line 495 "cggbal.f"
	    sum += work[i__ + *n];
#line 496 "cggbal.f"
L320:
#line 496 "cggbal.f"
	    ;
#line 496 "cggbal.f"
	}
#line 497 "cggbal.f"
	work[j + *n * 3] = (doublereal) kount * work[j] + sum;
#line 498 "cggbal.f"
/* L330: */
#line 498 "cggbal.f"
    }

#line 500 "cggbal.f"
    sum = sdot_(&nr, &work[*ilo + *n], &c__1, &work[*ilo + (*n << 1)], &c__1) 
	    + sdot_(&nr, &work[*ilo], &c__1, &work[*ilo + *n * 3], &c__1);
#line 502 "cggbal.f"
    alpha = gamma / sum;

/*     Determine correction to current iteration */

#line 506 "cggbal.f"
    cmax = 0.;
#line 507 "cggbal.f"
    i__1 = *ihi;
#line 507 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 508 "cggbal.f"
	cor = alpha * work[i__ + *n];
#line 509 "cggbal.f"
	if (abs(cor) > cmax) {
#line 509 "cggbal.f"
	    cmax = abs(cor);
#line 509 "cggbal.f"
	}
#line 511 "cggbal.f"
	lscale[i__] += cor;
#line 512 "cggbal.f"
	cor = alpha * work[i__];
#line 513 "cggbal.f"
	if (abs(cor) > cmax) {
#line 513 "cggbal.f"
	    cmax = abs(cor);
#line 513 "cggbal.f"
	}
#line 515 "cggbal.f"
	rscale[i__] += cor;
#line 516 "cggbal.f"
/* L340: */
#line 516 "cggbal.f"
    }
#line 517 "cggbal.f"
    if (cmax < .5) {
#line 517 "cggbal.f"
	goto L350;
#line 517 "cggbal.f"
    }

#line 520 "cggbal.f"
    d__1 = -alpha;
#line 520 "cggbal.f"
    saxpy_(&nr, &d__1, &work[*ilo + (*n << 1)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1);
#line 521 "cggbal.f"
    d__1 = -alpha;
#line 521 "cggbal.f"
    saxpy_(&nr, &d__1, &work[*ilo + *n * 3], &c__1, &work[*ilo + *n * 5], &
	    c__1);

#line 523 "cggbal.f"
    pgamma = gamma;
#line 524 "cggbal.f"
    ++it;
#line 525 "cggbal.f"
    if (it <= nrp2) {
#line 525 "cggbal.f"
	goto L250;
#line 525 "cggbal.f"
    }

/*     End generalized conjugate gradient iteration */

#line 530 "cggbal.f"
L350:
#line 531 "cggbal.f"
    sfmin = slamch_("S", (ftnlen)1);
#line 532 "cggbal.f"
    sfmax = 1. / sfmin;
#line 533 "cggbal.f"
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
#line 534 "cggbal.f"
    lsfmax = (integer) (d_lg10(&sfmax) / basl);
#line 535 "cggbal.f"
    i__1 = *ihi;
#line 535 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 536 "cggbal.f"
	i__2 = *n - *ilo + 1;
#line 536 "cggbal.f"
	irab = icamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
#line 537 "cggbal.f"
	rab = z_abs(&a[i__ + (irab + *ilo - 1) * a_dim1]);
#line 538 "cggbal.f"
	i__2 = *n - *ilo + 1;
#line 538 "cggbal.f"
	irab = icamax_(&i__2, &b[i__ + *ilo * b_dim1], ldb);
/* Computing MAX */
#line 539 "cggbal.f"
	d__1 = rab, d__2 = z_abs(&b[i__ + (irab + *ilo - 1) * b_dim1]);
#line 539 "cggbal.f"
	rab = max(d__1,d__2);
#line 540 "cggbal.f"
	d__1 = rab + sfmin;
#line 540 "cggbal.f"
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 541 "cggbal.f"
	ir = (integer) (lscale[i__] + d_sign(&c_b72, &lscale[i__]));
/* Computing MIN */
#line 542 "cggbal.f"
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
#line 542 "cggbal.f"
	ir = min(i__2,i__3);
#line 543 "cggbal.f"
	lscale[i__] = pow_di(&c_b36, &ir);
#line 544 "cggbal.f"
	icab = icamax_(ihi, &a[i__ * a_dim1 + 1], &c__1);
#line 545 "cggbal.f"
	cab = z_abs(&a[icab + i__ * a_dim1]);
#line 546 "cggbal.f"
	icab = icamax_(ihi, &b[i__ * b_dim1 + 1], &c__1);
/* Computing MAX */
#line 547 "cggbal.f"
	d__1 = cab, d__2 = z_abs(&b[icab + i__ * b_dim1]);
#line 547 "cggbal.f"
	cab = max(d__1,d__2);
#line 548 "cggbal.f"
	d__1 = cab + sfmin;
#line 548 "cggbal.f"
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 549 "cggbal.f"
	jc = (integer) (rscale[i__] + d_sign(&c_b72, &rscale[i__]));
/* Computing MIN */
#line 550 "cggbal.f"
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
#line 550 "cggbal.f"
	jc = min(i__2,i__3);
#line 551 "cggbal.f"
	rscale[i__] = pow_di(&c_b36, &jc);
#line 552 "cggbal.f"
/* L360: */
#line 552 "cggbal.f"
    }

/*     Row scaling of matrices A and B */

#line 556 "cggbal.f"
    i__1 = *ihi;
#line 556 "cggbal.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 557 "cggbal.f"
	i__2 = *n - *ilo + 1;
#line 557 "cggbal.f"
	csscal_(&i__2, &lscale[i__], &a[i__ + *ilo * a_dim1], lda);
#line 558 "cggbal.f"
	i__2 = *n - *ilo + 1;
#line 558 "cggbal.f"
	csscal_(&i__2, &lscale[i__], &b[i__ + *ilo * b_dim1], ldb);
#line 559 "cggbal.f"
/* L370: */
#line 559 "cggbal.f"
    }

/*     Column scaling of matrices A and B */

#line 563 "cggbal.f"
    i__1 = *ihi;
#line 563 "cggbal.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 564 "cggbal.f"
	csscal_(ihi, &rscale[j], &a[j * a_dim1 + 1], &c__1);
#line 565 "cggbal.f"
	csscal_(ihi, &rscale[j], &b[j * b_dim1 + 1], &c__1);
#line 566 "cggbal.f"
/* L380: */
#line 566 "cggbal.f"
    }

#line 568 "cggbal.f"
    return 0;

/*     End of CGGBAL */

} /* cggbal_ */

