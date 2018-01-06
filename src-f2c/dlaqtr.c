#line 1 "dlaqtr.f"
/* dlaqtr.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqtr.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b21 = 1.;
static doublereal c_b25 = 0.;
static logical c_true = TRUE_;

/* > \brief \b DLAQTR solves a real quasi-triangular system of equations, or a complex quasi-triangular system
 of special form, in real arithmetic. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            LREAL, LTRAN */
/*       INTEGER            INFO, LDT, N */
/*       DOUBLE PRECISION   SCALE, W */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( * ), T( LDT, * ), WORK( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQTR solves the real quasi-triangular system */
/* > */
/* >              op(T)*p = scale*c,               if LREAL = .TRUE. */
/* > */
/* > or the complex quasi-triangular systems */
/* > */
/* >            op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE. */
/* > */
/* > in real arithmetic, where T is upper quasi-triangular. */
/* > If LREAL = .FALSE., then the first diagonal block of T must be */
/* > 1 by 1, B is the specially structured matrix */
/* > */
/* >                B = [ b(1) b(2) ... b(n) ] */
/* >                    [       w            ] */
/* >                    [           w        ] */
/* >                    [              .     ] */
/* >                    [                 w  ] */
/* > */
/* > op(A) = A or A**T, A**T denotes the transpose of */
/* > matrix A. */
/* > */
/* > On input, X = [ c ].  On output, X = [ p ]. */
/* >               [ d ]                  [ q ] */
/* > */
/* > This subroutine is designed for the condition number estimation */
/* > in routine DTRSNA. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] LTRAN */
/* > \verbatim */
/* >          LTRAN is LOGICAL */
/* >          On entry, LTRAN specifies the option of conjugate transpose: */
/* >             = .FALSE.,    op(T+i*B) = T+i*B, */
/* >             = .TRUE.,     op(T+i*B) = (T+i*B)**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LREAL */
/* > \verbatim */
/* >          LREAL is LOGICAL */
/* >          On entry, LREAL specifies the input matrix structure: */
/* >             = .FALSE.,    the input is complex */
/* >             = .TRUE.,     the input is real */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          On entry, N specifies the order of T+i*B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
/* >          On entry, T contains a matrix in Schur canonical form. */
/* >          If LREAL = .FALSE., then the first diagonal block of T mu */
/* >          be 1 by 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the matrix T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, B contains the elements to form the matrix */
/* >          B as described above. */
/* >          If LREAL = .TRUE., B is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION */
/* >          On entry, W is the diagonal element of the matrix B. */
/* >          If LREAL = .TRUE., W is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          On exit, SCALE is the scale factor. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (2*N) */
/* >          On entry, X contains the right hand side of the system. */
/* >          On exit, X is overwritten by the solution. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          On exit, INFO is set to */
/* >             0: successful exit. */
/* >               1: the some diagonal 1 by 1 block has been perturbed by */
/* >                  a small number SMIN to keep nonsingularity. */
/* >               2: the some diagonal 2 by 2 block has been perturbed by */
/* >                  a small number in DLALN2 to keep nonsingularity. */
/* >          NOTE: In the interests of speed, this routine does not */
/* >                check the inputs for errors. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaqtr_(logical *ltran, logical *lreal, integer *n, 
	doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal 
	*scale, doublereal *x, doublereal *work, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static doublereal d__[4]	/* was [2][2] */;
    static integer i__, j, k;
    static doublereal v[4]	/* was [2][2] */, z__;
    static integer j1, j2, n1, n2;
    static doublereal si, xj, sr, rec, eps, tjj, tmp;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin, xmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer jnext;
    static doublereal sminw, xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bignum;
    static logical notran;
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

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

/*     Do not test the input parameters for errors */

#line 212 "dlaqtr.f"
    /* Parameter adjustments */
#line 212 "dlaqtr.f"
    t_dim1 = *ldt;
#line 212 "dlaqtr.f"
    t_offset = 1 + t_dim1;
#line 212 "dlaqtr.f"
    t -= t_offset;
#line 212 "dlaqtr.f"
    --b;
#line 212 "dlaqtr.f"
    --x;
#line 212 "dlaqtr.f"
    --work;
#line 212 "dlaqtr.f"

#line 212 "dlaqtr.f"
    /* Function Body */
#line 212 "dlaqtr.f"
    notran = ! (*ltran);
#line 213 "dlaqtr.f"
    *info = 0;

/*     Quick return if possible */

#line 217 "dlaqtr.f"
    if (*n == 0) {
#line 217 "dlaqtr.f"
	return 0;
#line 217 "dlaqtr.f"
    }

/*     Set constants to control overflow */

#line 222 "dlaqtr.f"
    eps = dlamch_("P", (ftnlen)1);
#line 223 "dlaqtr.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 224 "dlaqtr.f"
    bignum = 1. / smlnum;

#line 226 "dlaqtr.f"
    xnorm = dlange_("M", n, n, &t[t_offset], ldt, d__, (ftnlen)1);
#line 227 "dlaqtr.f"
    if (! (*lreal)) {
/* Computing MAX */
#line 227 "dlaqtr.f"
	d__1 = xnorm, d__2 = abs(*w), d__1 = max(d__1,d__2), d__2 = dlange_(
		"M", n, &c__1, &b[1], n, d__, (ftnlen)1);
#line 227 "dlaqtr.f"
	xnorm = max(d__1,d__2);
#line 227 "dlaqtr.f"
    }
/* Computing MAX */
#line 229 "dlaqtr.f"
    d__1 = smlnum, d__2 = eps * xnorm;
#line 229 "dlaqtr.f"
    smin = max(d__1,d__2);

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 234 "dlaqtr.f"
    work[1] = 0.;
#line 235 "dlaqtr.f"
    i__1 = *n;
#line 235 "dlaqtr.f"
    for (j = 2; j <= i__1; ++j) {
#line 236 "dlaqtr.f"
	i__2 = j - 1;
#line 236 "dlaqtr.f"
	work[j] = dasum_(&i__2, &t[j * t_dim1 + 1], &c__1);
#line 237 "dlaqtr.f"
/* L10: */
#line 237 "dlaqtr.f"
    }

#line 239 "dlaqtr.f"
    if (! (*lreal)) {
#line 240 "dlaqtr.f"
	i__1 = *n;
#line 240 "dlaqtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 241 "dlaqtr.f"
	    work[i__] += (d__1 = b[i__], abs(d__1));
#line 242 "dlaqtr.f"
/* L20: */
#line 242 "dlaqtr.f"
	}
#line 243 "dlaqtr.f"
    }

#line 245 "dlaqtr.f"
    n2 = *n << 1;
#line 246 "dlaqtr.f"
    n1 = *n;
#line 247 "dlaqtr.f"
    if (! (*lreal)) {
#line 247 "dlaqtr.f"
	n1 = n2;
#line 247 "dlaqtr.f"
    }
#line 249 "dlaqtr.f"
    k = idamax_(&n1, &x[1], &c__1);
#line 250 "dlaqtr.f"
    xmax = (d__1 = x[k], abs(d__1));
#line 251 "dlaqtr.f"
    *scale = 1.;

#line 253 "dlaqtr.f"
    if (xmax > bignum) {
#line 254 "dlaqtr.f"
	*scale = bignum / xmax;
#line 255 "dlaqtr.f"
	dscal_(&n1, scale, &x[1], &c__1);
#line 256 "dlaqtr.f"
	xmax = bignum;
#line 257 "dlaqtr.f"
    }

#line 259 "dlaqtr.f"
    if (*lreal) {

#line 261 "dlaqtr.f"
	if (notran) {

/*           Solve T*p = scale*c */

#line 265 "dlaqtr.f"
	    jnext = *n;
#line 266 "dlaqtr.f"
	    for (j = *n; j >= 1; --j) {
#line 267 "dlaqtr.f"
		if (j > jnext) {
#line 267 "dlaqtr.f"
		    goto L30;
#line 267 "dlaqtr.f"
		}
#line 269 "dlaqtr.f"
		j1 = j;
#line 270 "dlaqtr.f"
		j2 = j;
#line 271 "dlaqtr.f"
		jnext = j - 1;
#line 272 "dlaqtr.f"
		if (j > 1) {
#line 273 "dlaqtr.f"
		    if (t[j + (j - 1) * t_dim1] != 0.) {
#line 274 "dlaqtr.f"
			j1 = j - 1;
#line 275 "dlaqtr.f"
			jnext = j - 2;
#line 276 "dlaqtr.f"
		    }
#line 277 "dlaqtr.f"
		}

#line 279 "dlaqtr.f"
		if (j1 == j2) {

/*                 Meet 1 by 1 diagonal block */

/*                 Scale to avoid overflow when computing */
/*                     x(j) = b(j)/T(j,j) */

#line 286 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));
#line 287 "dlaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1));
#line 288 "dlaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 289 "dlaqtr.f"
		    if (tjj < smin) {
#line 290 "dlaqtr.f"
			tmp = smin;
#line 291 "dlaqtr.f"
			tjj = smin;
#line 292 "dlaqtr.f"
			*info = 1;
#line 293 "dlaqtr.f"
		    }

#line 295 "dlaqtr.f"
		    if (xj == 0.) {
#line 295 "dlaqtr.f"
			goto L30;
#line 295 "dlaqtr.f"
		    }

#line 298 "dlaqtr.f"
		    if (tjj < 1.) {
#line 299 "dlaqtr.f"
			if (xj > bignum * tjj) {
#line 300 "dlaqtr.f"
			    rec = 1. / xj;
#line 301 "dlaqtr.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 302 "dlaqtr.f"
			    *scale *= rec;
#line 303 "dlaqtr.f"
			    xmax *= rec;
#line 304 "dlaqtr.f"
			}
#line 305 "dlaqtr.f"
		    }
#line 306 "dlaqtr.f"
		    x[j1] /= tmp;
#line 307 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));

/*                 Scale x if necessary to avoid overflow when adding a */
/*                 multiple of column j1 of T. */

#line 312 "dlaqtr.f"
		    if (xj > 1.) {
#line 313 "dlaqtr.f"
			rec = 1. / xj;
#line 314 "dlaqtr.f"
			if (work[j1] > (bignum - xmax) * rec) {
#line 315 "dlaqtr.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 316 "dlaqtr.f"
			    *scale *= rec;
#line 317 "dlaqtr.f"
			}
#line 318 "dlaqtr.f"
		    }
#line 319 "dlaqtr.f"
		    if (j1 > 1) {
#line 320 "dlaqtr.f"
			i__1 = j1 - 1;
#line 320 "dlaqtr.f"
			d__1 = -x[j1];
#line 320 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 321 "dlaqtr.f"
			i__1 = j1 - 1;
#line 321 "dlaqtr.f"
			k = idamax_(&i__1, &x[1], &c__1);
#line 322 "dlaqtr.f"
			xmax = (d__1 = x[k], abs(d__1));
#line 323 "dlaqtr.f"
		    }

#line 325 "dlaqtr.f"
		} else {

/*                 Meet 2 by 2 diagonal block */

/*                 Call 2 by 2 linear system solve, to take */
/*                 care of possible overflow by scaling factor. */

#line 332 "dlaqtr.f"
		    d__[0] = x[j1];
#line 333 "dlaqtr.f"
		    d__[1] = x[j2];
#line 334 "dlaqtr.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b21, &t[j1 + j1 
			    * t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &
			    c_b25, &c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
#line 337 "dlaqtr.f"
		    if (ierr != 0) {
#line 337 "dlaqtr.f"
			*info = 2;
#line 337 "dlaqtr.f"
		    }

#line 340 "dlaqtr.f"
		    if (scaloc != 1.) {
#line 341 "dlaqtr.f"
			dscal_(n, &scaloc, &x[1], &c__1);
#line 342 "dlaqtr.f"
			*scale *= scaloc;
#line 343 "dlaqtr.f"
		    }
#line 344 "dlaqtr.f"
		    x[j1] = v[0];
#line 345 "dlaqtr.f"
		    x[j2] = v[1];

/*                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2)) */
/*                 to avoid overflow in updating right-hand side. */

/* Computing MAX */
#line 350 "dlaqtr.f"
		    d__1 = abs(v[0]), d__2 = abs(v[1]);
#line 350 "dlaqtr.f"
		    xj = max(d__1,d__2);
#line 351 "dlaqtr.f"
		    if (xj > 1.) {
#line 352 "dlaqtr.f"
			rec = 1. / xj;
/* Computing MAX */
#line 353 "dlaqtr.f"
			d__1 = work[j1], d__2 = work[j2];
#line 353 "dlaqtr.f"
			if (max(d__1,d__2) > (bignum - xmax) * rec) {
#line 355 "dlaqtr.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 356 "dlaqtr.f"
			    *scale *= rec;
#line 357 "dlaqtr.f"
			}
#line 358 "dlaqtr.f"
		    }

/*                 Update right-hand side */

#line 362 "dlaqtr.f"
		    if (j1 > 1) {
#line 363 "dlaqtr.f"
			i__1 = j1 - 1;
#line 363 "dlaqtr.f"
			d__1 = -x[j1];
#line 363 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 364 "dlaqtr.f"
			i__1 = j1 - 1;
#line 364 "dlaqtr.f"
			d__1 = -x[j2];
#line 364 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 365 "dlaqtr.f"
			i__1 = j1 - 1;
#line 365 "dlaqtr.f"
			k = idamax_(&i__1, &x[1], &c__1);
#line 366 "dlaqtr.f"
			xmax = (d__1 = x[k], abs(d__1));
#line 367 "dlaqtr.f"
		    }

#line 369 "dlaqtr.f"
		}

#line 371 "dlaqtr.f"
L30:
#line 371 "dlaqtr.f"
		;
#line 371 "dlaqtr.f"
	    }

#line 373 "dlaqtr.f"
	} else {

/*           Solve T**T*p = scale*c */

#line 377 "dlaqtr.f"
	    jnext = 1;
#line 378 "dlaqtr.f"
	    i__1 = *n;
#line 378 "dlaqtr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 379 "dlaqtr.f"
		if (j < jnext) {
#line 379 "dlaqtr.f"
		    goto L40;
#line 379 "dlaqtr.f"
		}
#line 381 "dlaqtr.f"
		j1 = j;
#line 382 "dlaqtr.f"
		j2 = j;
#line 383 "dlaqtr.f"
		jnext = j + 1;
#line 384 "dlaqtr.f"
		if (j < *n) {
#line 385 "dlaqtr.f"
		    if (t[j + 1 + j * t_dim1] != 0.) {
#line 386 "dlaqtr.f"
			j2 = j + 1;
#line 387 "dlaqtr.f"
			jnext = j + 2;
#line 388 "dlaqtr.f"
		    }
#line 389 "dlaqtr.f"
		}

#line 391 "dlaqtr.f"
		if (j1 == j2) {

/*                 1 by 1 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side element by inner product. */

#line 398 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));
#line 399 "dlaqtr.f"
		    if (xmax > 1.) {
#line 400 "dlaqtr.f"
			rec = 1. / xmax;
#line 401 "dlaqtr.f"
			if (work[j1] > (bignum - xj) * rec) {
#line 402 "dlaqtr.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 403 "dlaqtr.f"
			    *scale *= rec;
#line 404 "dlaqtr.f"
			    xmax *= rec;
#line 405 "dlaqtr.f"
			}
#line 406 "dlaqtr.f"
		    }

#line 408 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 408 "dlaqtr.f"
		    x[j1] -= ddot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &
			    c__1);

#line 410 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));
#line 411 "dlaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1));
#line 412 "dlaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 413 "dlaqtr.f"
		    if (tjj < smin) {
#line 414 "dlaqtr.f"
			tmp = smin;
#line 415 "dlaqtr.f"
			tjj = smin;
#line 416 "dlaqtr.f"
			*info = 1;
#line 417 "dlaqtr.f"
		    }

#line 419 "dlaqtr.f"
		    if (tjj < 1.) {
#line 420 "dlaqtr.f"
			if (xj > bignum * tjj) {
#line 421 "dlaqtr.f"
			    rec = 1. / xj;
#line 422 "dlaqtr.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 423 "dlaqtr.f"
			    *scale *= rec;
#line 424 "dlaqtr.f"
			    xmax *= rec;
#line 425 "dlaqtr.f"
			}
#line 426 "dlaqtr.f"
		    }
#line 427 "dlaqtr.f"
		    x[j1] /= tmp;
/* Computing MAX */
#line 428 "dlaqtr.f"
		    d__2 = xmax, d__3 = (d__1 = x[j1], abs(d__1));
#line 428 "dlaqtr.f"
		    xmax = max(d__2,d__3);

#line 430 "dlaqtr.f"
		} else {

/*                 2 by 2 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side elements by inner product. */

/* Computing MAX */
#line 437 "dlaqtr.f"
		    d__3 = (d__1 = x[j1], abs(d__1)), d__4 = (d__2 = x[j2], 
			    abs(d__2));
#line 437 "dlaqtr.f"
		    xj = max(d__3,d__4);
#line 438 "dlaqtr.f"
		    if (xmax > 1.) {
#line 439 "dlaqtr.f"
			rec = 1. / xmax;
/* Computing MAX */
#line 440 "dlaqtr.f"
			d__1 = work[j2], d__2 = work[j1];
#line 440 "dlaqtr.f"
			if (max(d__1,d__2) > (bignum - xj) * rec) {
#line 442 "dlaqtr.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 443 "dlaqtr.f"
			    *scale *= rec;
#line 444 "dlaqtr.f"
			    xmax *= rec;
#line 445 "dlaqtr.f"
			}
#line 446 "dlaqtr.f"
		    }

#line 448 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 448 "dlaqtr.f"
		    d__[0] = x[j1] - ddot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);
#line 450 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 450 "dlaqtr.f"
		    d__[1] = x[j2] - ddot_(&i__2, &t[j2 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);

#line 453 "dlaqtr.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b21, &t[j1 + j1 *
			     t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &c_b25,
			     &c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
#line 456 "dlaqtr.f"
		    if (ierr != 0) {
#line 456 "dlaqtr.f"
			*info = 2;
#line 456 "dlaqtr.f"
		    }

#line 459 "dlaqtr.f"
		    if (scaloc != 1.) {
#line 460 "dlaqtr.f"
			dscal_(n, &scaloc, &x[1], &c__1);
#line 461 "dlaqtr.f"
			*scale *= scaloc;
#line 462 "dlaqtr.f"
		    }
#line 463 "dlaqtr.f"
		    x[j1] = v[0];
#line 464 "dlaqtr.f"
		    x[j2] = v[1];
/* Computing MAX */
#line 465 "dlaqtr.f"
		    d__3 = (d__1 = x[j1], abs(d__1)), d__4 = (d__2 = x[j2], 
			    abs(d__2)), d__3 = max(d__3,d__4);
#line 465 "dlaqtr.f"
		    xmax = max(d__3,xmax);

#line 467 "dlaqtr.f"
		}
#line 468 "dlaqtr.f"
L40:
#line 468 "dlaqtr.f"
		;
#line 468 "dlaqtr.f"
	    }
#line 469 "dlaqtr.f"
	}

#line 471 "dlaqtr.f"
    } else {

/* Computing MAX */
#line 473 "dlaqtr.f"
	d__1 = eps * abs(*w);
#line 473 "dlaqtr.f"
	sminw = max(d__1,smin);
#line 474 "dlaqtr.f"
	if (notran) {

/*           Solve (T + iB)*(p+iq) = c+id */

#line 478 "dlaqtr.f"
	    jnext = *n;
#line 479 "dlaqtr.f"
	    for (j = *n; j >= 1; --j) {
#line 480 "dlaqtr.f"
		if (j > jnext) {
#line 480 "dlaqtr.f"
		    goto L70;
#line 480 "dlaqtr.f"
		}
#line 482 "dlaqtr.f"
		j1 = j;
#line 483 "dlaqtr.f"
		j2 = j;
#line 484 "dlaqtr.f"
		jnext = j - 1;
#line 485 "dlaqtr.f"
		if (j > 1) {
#line 486 "dlaqtr.f"
		    if (t[j + (j - 1) * t_dim1] != 0.) {
#line 487 "dlaqtr.f"
			j1 = j - 1;
#line 488 "dlaqtr.f"
			jnext = j - 2;
#line 489 "dlaqtr.f"
		    }
#line 490 "dlaqtr.f"
		}

#line 492 "dlaqtr.f"
		if (j1 == j2) {

/*                 1 by 1 diagonal block */

/*                 Scale if necessary to avoid overflow in division */

#line 498 "dlaqtr.f"
		    z__ = *w;
#line 499 "dlaqtr.f"
		    if (j1 == 1) {
#line 499 "dlaqtr.f"
			z__ = b[1];
#line 499 "dlaqtr.f"
		    }
#line 501 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], abs(
			    d__2));
#line 502 "dlaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1)) + abs(z__);
#line 503 "dlaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 504 "dlaqtr.f"
		    if (tjj < sminw) {
#line 505 "dlaqtr.f"
			tmp = sminw;
#line 506 "dlaqtr.f"
			tjj = sminw;
#line 507 "dlaqtr.f"
			*info = 1;
#line 508 "dlaqtr.f"
		    }

#line 510 "dlaqtr.f"
		    if (xj == 0.) {
#line 510 "dlaqtr.f"
			goto L70;
#line 510 "dlaqtr.f"
		    }

#line 513 "dlaqtr.f"
		    if (tjj < 1.) {
#line 514 "dlaqtr.f"
			if (xj > bignum * tjj) {
#line 515 "dlaqtr.f"
			    rec = 1. / xj;
#line 516 "dlaqtr.f"
			    dscal_(&n2, &rec, &x[1], &c__1);
#line 517 "dlaqtr.f"
			    *scale *= rec;
#line 518 "dlaqtr.f"
			    xmax *= rec;
#line 519 "dlaqtr.f"
			}
#line 520 "dlaqtr.f"
		    }
#line 521 "dlaqtr.f"
		    dladiv_(&x[j1], &x[*n + j1], &tmp, &z__, &sr, &si);
#line 522 "dlaqtr.f"
		    x[j1] = sr;
#line 523 "dlaqtr.f"
		    x[*n + j1] = si;
#line 524 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], abs(
			    d__2));

/*                 Scale x if necessary to avoid overflow when adding a */
/*                 multiple of column j1 of T. */

#line 529 "dlaqtr.f"
		    if (xj > 1.) {
#line 530 "dlaqtr.f"
			rec = 1. / xj;
#line 531 "dlaqtr.f"
			if (work[j1] > (bignum - xmax) * rec) {
#line 532 "dlaqtr.f"
			    dscal_(&n2, &rec, &x[1], &c__1);
#line 533 "dlaqtr.f"
			    *scale *= rec;
#line 534 "dlaqtr.f"
			}
#line 535 "dlaqtr.f"
		    }

#line 537 "dlaqtr.f"
		    if (j1 > 1) {
#line 538 "dlaqtr.f"
			i__1 = j1 - 1;
#line 538 "dlaqtr.f"
			d__1 = -x[j1];
#line 538 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 539 "dlaqtr.f"
			i__1 = j1 - 1;
#line 539 "dlaqtr.f"
			d__1 = -x[*n + j1];
#line 539 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[*
				n + 1], &c__1);

#line 542 "dlaqtr.f"
			x[1] += b[j1] * x[*n + j1];
#line 543 "dlaqtr.f"
			x[*n + 1] -= b[j1] * x[j1];

#line 545 "dlaqtr.f"
			xmax = 0.;
#line 546 "dlaqtr.f"
			i__1 = j1 - 1;
#line 546 "dlaqtr.f"
			for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 547 "dlaqtr.f"
			    d__3 = xmax, d__4 = (d__1 = x[k], abs(d__1)) + (
				    d__2 = x[k + *n], abs(d__2));
#line 547 "dlaqtr.f"
			    xmax = max(d__3,d__4);
#line 549 "dlaqtr.f"
/* L50: */
#line 549 "dlaqtr.f"
			}
#line 550 "dlaqtr.f"
		    }

#line 552 "dlaqtr.f"
		} else {

/*                 Meet 2 by 2 diagonal block */

#line 556 "dlaqtr.f"
		    d__[0] = x[j1];
#line 557 "dlaqtr.f"
		    d__[1] = x[j2];
#line 558 "dlaqtr.f"
		    d__[2] = x[*n + j1];
#line 559 "dlaqtr.f"
		    d__[3] = x[*n + j2];
#line 560 "dlaqtr.f"
		    d__1 = -(*w);
#line 560 "dlaqtr.f"
		    dlaln2_(&c_false, &c__2, &c__2, &sminw, &c_b21, &t[j1 + 
			    j1 * t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &
			    c_b25, &d__1, v, &c__2, &scaloc, &xnorm, &ierr);
#line 563 "dlaqtr.f"
		    if (ierr != 0) {
#line 563 "dlaqtr.f"
			*info = 2;
#line 563 "dlaqtr.f"
		    }

#line 566 "dlaqtr.f"
		    if (scaloc != 1.) {
#line 567 "dlaqtr.f"
			i__1 = *n << 1;
#line 567 "dlaqtr.f"
			dscal_(&i__1, &scaloc, &x[1], &c__1);
#line 568 "dlaqtr.f"
			*scale = scaloc * *scale;
#line 569 "dlaqtr.f"
		    }
#line 570 "dlaqtr.f"
		    x[j1] = v[0];
#line 571 "dlaqtr.f"
		    x[j2] = v[1];
#line 572 "dlaqtr.f"
		    x[*n + j1] = v[2];
#line 573 "dlaqtr.f"
		    x[*n + j2] = v[3];

/*                 Scale X(J1), .... to avoid overflow in */
/*                 updating right hand side. */

/* Computing MAX */
#line 578 "dlaqtr.f"
		    d__1 = abs(v[0]) + abs(v[2]), d__2 = abs(v[1]) + abs(v[3])
			    ;
#line 578 "dlaqtr.f"
		    xj = max(d__1,d__2);
#line 580 "dlaqtr.f"
		    if (xj > 1.) {
#line 581 "dlaqtr.f"
			rec = 1. / xj;
/* Computing MAX */
#line 582 "dlaqtr.f"
			d__1 = work[j1], d__2 = work[j2];
#line 582 "dlaqtr.f"
			if (max(d__1,d__2) > (bignum - xmax) * rec) {
#line 584 "dlaqtr.f"
			    dscal_(&n2, &rec, &x[1], &c__1);
#line 585 "dlaqtr.f"
			    *scale *= rec;
#line 586 "dlaqtr.f"
			}
#line 587 "dlaqtr.f"
		    }

/*                 Update the right-hand side. */

#line 591 "dlaqtr.f"
		    if (j1 > 1) {
#line 592 "dlaqtr.f"
			i__1 = j1 - 1;
#line 592 "dlaqtr.f"
			d__1 = -x[j1];
#line 592 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 593 "dlaqtr.f"
			i__1 = j1 - 1;
#line 593 "dlaqtr.f"
			d__1 = -x[j2];
#line 593 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);

#line 595 "dlaqtr.f"
			i__1 = j1 - 1;
#line 595 "dlaqtr.f"
			d__1 = -x[*n + j1];
#line 595 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[*
				n + 1], &c__1);
#line 597 "dlaqtr.f"
			i__1 = j1 - 1;
#line 597 "dlaqtr.f"
			d__1 = -x[*n + j2];
#line 597 "dlaqtr.f"
			daxpy_(&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[*
				n + 1], &c__1);

#line 600 "dlaqtr.f"
			x[1] = x[1] + b[j1] * x[*n + j1] + b[j2] * x[*n + j2];
#line 602 "dlaqtr.f"
			x[*n + 1] = x[*n + 1] - b[j1] * x[j1] - b[j2] * x[j2];

#line 605 "dlaqtr.f"
			xmax = 0.;
#line 606 "dlaqtr.f"
			i__1 = j1 - 1;
#line 606 "dlaqtr.f"
			for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 607 "dlaqtr.f"
			    d__3 = (d__1 = x[k], abs(d__1)) + (d__2 = x[k + *
				    n], abs(d__2));
#line 607 "dlaqtr.f"
			    xmax = max(d__3,xmax);
#line 609 "dlaqtr.f"
/* L60: */
#line 609 "dlaqtr.f"
			}
#line 610 "dlaqtr.f"
		    }

#line 612 "dlaqtr.f"
		}
#line 613 "dlaqtr.f"
L70:
#line 613 "dlaqtr.f"
		;
#line 613 "dlaqtr.f"
	    }

#line 615 "dlaqtr.f"
	} else {

/*           Solve (T + iB)**T*(p+iq) = c+id */

#line 619 "dlaqtr.f"
	    jnext = 1;
#line 620 "dlaqtr.f"
	    i__1 = *n;
#line 620 "dlaqtr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 621 "dlaqtr.f"
		if (j < jnext) {
#line 621 "dlaqtr.f"
		    goto L80;
#line 621 "dlaqtr.f"
		}
#line 623 "dlaqtr.f"
		j1 = j;
#line 624 "dlaqtr.f"
		j2 = j;
#line 625 "dlaqtr.f"
		jnext = j + 1;
#line 626 "dlaqtr.f"
		if (j < *n) {
#line 627 "dlaqtr.f"
		    if (t[j + 1 + j * t_dim1] != 0.) {
#line 628 "dlaqtr.f"
			j2 = j + 1;
#line 629 "dlaqtr.f"
			jnext = j + 2;
#line 630 "dlaqtr.f"
		    }
#line 631 "dlaqtr.f"
		}

#line 633 "dlaqtr.f"
		if (j1 == j2) {

/*                 1 by 1 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side element by inner product. */

#line 640 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[j1 + *n], abs(
			    d__2));
#line 641 "dlaqtr.f"
		    if (xmax > 1.) {
#line 642 "dlaqtr.f"
			rec = 1. / xmax;
#line 643 "dlaqtr.f"
			if (work[j1] > (bignum - xj) * rec) {
#line 644 "dlaqtr.f"
			    dscal_(&n2, &rec, &x[1], &c__1);
#line 645 "dlaqtr.f"
			    *scale *= rec;
#line 646 "dlaqtr.f"
			    xmax *= rec;
#line 647 "dlaqtr.f"
			}
#line 648 "dlaqtr.f"
		    }

#line 650 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 650 "dlaqtr.f"
		    x[j1] -= ddot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &
			    c__1);
#line 651 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 651 "dlaqtr.f"
		    x[*n + j1] -= ddot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[
			    *n + 1], &c__1);
#line 653 "dlaqtr.f"
		    if (j1 > 1) {
#line 654 "dlaqtr.f"
			x[j1] -= b[j1] * x[*n + 1];
#line 655 "dlaqtr.f"
			x[*n + j1] += b[j1] * x[1];
#line 656 "dlaqtr.f"
		    }
#line 657 "dlaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[j1 + *n], abs(
			    d__2));

#line 659 "dlaqtr.f"
		    z__ = *w;
#line 660 "dlaqtr.f"
		    if (j1 == 1) {
#line 660 "dlaqtr.f"
			z__ = b[1];
#line 660 "dlaqtr.f"
		    }

/*                 Scale if necessary to avoid overflow in */
/*                 complex division */

#line 666 "dlaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1)) + abs(z__);
#line 667 "dlaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 668 "dlaqtr.f"
		    if (tjj < sminw) {
#line 669 "dlaqtr.f"
			tmp = sminw;
#line 670 "dlaqtr.f"
			tjj = sminw;
#line 671 "dlaqtr.f"
			*info = 1;
#line 672 "dlaqtr.f"
		    }

#line 674 "dlaqtr.f"
		    if (tjj < 1.) {
#line 675 "dlaqtr.f"
			if (xj > bignum * tjj) {
#line 676 "dlaqtr.f"
			    rec = 1. / xj;
#line 677 "dlaqtr.f"
			    dscal_(&n2, &rec, &x[1], &c__1);
#line 678 "dlaqtr.f"
			    *scale *= rec;
#line 679 "dlaqtr.f"
			    xmax *= rec;
#line 680 "dlaqtr.f"
			}
#line 681 "dlaqtr.f"
		    }
#line 682 "dlaqtr.f"
		    d__1 = -z__;
#line 682 "dlaqtr.f"
		    dladiv_(&x[j1], &x[*n + j1], &tmp, &d__1, &sr, &si);
#line 683 "dlaqtr.f"
		    x[j1] = sr;
#line 684 "dlaqtr.f"
		    x[j1 + *n] = si;
/* Computing MAX */
#line 685 "dlaqtr.f"
		    d__3 = (d__1 = x[j1], abs(d__1)) + (d__2 = x[j1 + *n], 
			    abs(d__2));
#line 685 "dlaqtr.f"
		    xmax = max(d__3,xmax);

#line 687 "dlaqtr.f"
		} else {

/*                 2 by 2 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side element by inner product. */

/* Computing MAX */
#line 694 "dlaqtr.f"
		    d__5 = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], 
			    abs(d__2)), d__6 = (d__3 = x[j2], abs(d__3)) + (
			    d__4 = x[*n + j2], abs(d__4));
#line 694 "dlaqtr.f"
		    xj = max(d__5,d__6);
#line 696 "dlaqtr.f"
		    if (xmax > 1.) {
#line 697 "dlaqtr.f"
			rec = 1. / xmax;
/* Computing MAX */
#line 698 "dlaqtr.f"
			d__1 = work[j1], d__2 = work[j2];
#line 698 "dlaqtr.f"
			if (max(d__1,d__2) > (bignum - xj) / xmax) {
#line 700 "dlaqtr.f"
			    dscal_(&n2, &rec, &x[1], &c__1);
#line 701 "dlaqtr.f"
			    *scale *= rec;
#line 702 "dlaqtr.f"
			    xmax *= rec;
#line 703 "dlaqtr.f"
			}
#line 704 "dlaqtr.f"
		    }

#line 706 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 706 "dlaqtr.f"
		    d__[0] = x[j1] - ddot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);
#line 708 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 708 "dlaqtr.f"
		    d__[1] = x[j2] - ddot_(&i__2, &t[j2 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);
#line 710 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 710 "dlaqtr.f"
		    d__[2] = x[*n + j1] - ddot_(&i__2, &t[j1 * t_dim1 + 1], &
			    c__1, &x[*n + 1], &c__1);
#line 712 "dlaqtr.f"
		    i__2 = j1 - 1;
#line 712 "dlaqtr.f"
		    d__[3] = x[*n + j2] - ddot_(&i__2, &t[j2 * t_dim1 + 1], &
			    c__1, &x[*n + 1], &c__1);
#line 714 "dlaqtr.f"
		    d__[0] -= b[j1] * x[*n + 1];
#line 715 "dlaqtr.f"
		    d__[1] -= b[j2] * x[*n + 1];
#line 716 "dlaqtr.f"
		    d__[2] += b[j1] * x[1];
#line 717 "dlaqtr.f"
		    d__[3] += b[j2] * x[1];

#line 719 "dlaqtr.f"
		    dlaln2_(&c_true, &c__2, &c__2, &sminw, &c_b21, &t[j1 + j1 
			    * t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &
			    c_b25, w, v, &c__2, &scaloc, &xnorm, &ierr);
#line 722 "dlaqtr.f"
		    if (ierr != 0) {
#line 722 "dlaqtr.f"
			*info = 2;
#line 722 "dlaqtr.f"
		    }

#line 725 "dlaqtr.f"
		    if (scaloc != 1.) {
#line 726 "dlaqtr.f"
			dscal_(&n2, &scaloc, &x[1], &c__1);
#line 727 "dlaqtr.f"
			*scale = scaloc * *scale;
#line 728 "dlaqtr.f"
		    }
#line 729 "dlaqtr.f"
		    x[j1] = v[0];
#line 730 "dlaqtr.f"
		    x[j2] = v[1];
#line 731 "dlaqtr.f"
		    x[*n + j1] = v[2];
#line 732 "dlaqtr.f"
		    x[*n + j2] = v[3];
/* Computing MAX */
#line 733 "dlaqtr.f"
		    d__5 = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], 
			    abs(d__2)), d__6 = (d__3 = x[j2], abs(d__3)) + (
			    d__4 = x[*n + j2], abs(d__4)), d__5 = max(d__5,
			    d__6);
#line 733 "dlaqtr.f"
		    xmax = max(d__5,xmax);

#line 736 "dlaqtr.f"
		}

#line 738 "dlaqtr.f"
L80:
#line 738 "dlaqtr.f"
		;
#line 738 "dlaqtr.f"
	    }

#line 740 "dlaqtr.f"
	}

#line 742 "dlaqtr.f"
    }

#line 744 "dlaqtr.f"
    return 0;

/*     End of DLAQTR */

} /* dlaqtr_ */

