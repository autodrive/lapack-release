#line 1 "slaqtr.f"
/* slaqtr.f -- translated by f2c (version 20100827).
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

#line 1 "slaqtr.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b21 = 1.;
static doublereal c_b25 = 0.;
static logical c_true = TRUE_;

/* > \brief \b SLAQTR solves a real quasi-triangular system of equations, or a complex quasi-triangular system
 of special form, in real arithmetic. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            LREAL, LTRAN */
/*       INTEGER            INFO, LDT, N */
/*       REAL               SCALE, W */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               B( * ), T( LDT, * ), WORK( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQTR solves the real quasi-triangular system */
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
/* > in routine STRSNA. */
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
/* >          T is REAL array, dimension (LDT,N) */
/* >          On entry, T contains a matrix in Schur canonical form. */
/* >          If LREAL = .FALSE., then the first diagonal block of T must */
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
/* >          B is REAL array, dimension (N) */
/* >          On entry, B contains the elements to form the matrix */
/* >          B as described above. */
/* >          If LREAL = .TRUE., B is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is REAL */
/* >          On entry, W is the diagonal element of the matrix B. */
/* >          If LREAL = .TRUE., W is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          On exit, SCALE is the scale factor. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (2*N) */
/* >          On entry, X contains the right hand side of the system. */
/* >          On exit, X is overwritten by the solution. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
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
/* >                  a small number in SLALN2 to keep nonsingularity. */
/* >          NOTE: In the interests of speed, this routine does not */
/* >                check the inputs for errors. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaqtr_(logical *ltran, logical *lreal, integer *n, 
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
    static integer ierr;
    static doublereal smin;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal xmax;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer jnext;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal sminw, xnorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), slaln2_(logical *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , integer *);
    static doublereal scaloc;
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int sladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static logical notran;
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 212 "slaqtr.f"
    /* Parameter adjustments */
#line 212 "slaqtr.f"
    t_dim1 = *ldt;
#line 212 "slaqtr.f"
    t_offset = 1 + t_dim1;
#line 212 "slaqtr.f"
    t -= t_offset;
#line 212 "slaqtr.f"
    --b;
#line 212 "slaqtr.f"
    --x;
#line 212 "slaqtr.f"
    --work;
#line 212 "slaqtr.f"

#line 212 "slaqtr.f"
    /* Function Body */
#line 212 "slaqtr.f"
    notran = ! (*ltran);
#line 213 "slaqtr.f"
    *info = 0;

/*     Quick return if possible */

#line 217 "slaqtr.f"
    if (*n == 0) {
#line 217 "slaqtr.f"
	return 0;
#line 217 "slaqtr.f"
    }

/*     Set constants to control overflow */

#line 222 "slaqtr.f"
    eps = slamch_("P", (ftnlen)1);
#line 223 "slaqtr.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 224 "slaqtr.f"
    bignum = 1. / smlnum;

#line 226 "slaqtr.f"
    xnorm = slange_("M", n, n, &t[t_offset], ldt, d__, (ftnlen)1);
#line 227 "slaqtr.f"
    if (! (*lreal)) {
/* Computing MAX */
#line 227 "slaqtr.f"
	d__1 = xnorm, d__2 = abs(*w), d__1 = max(d__1,d__2), d__2 = slange_(
		"M", n, &c__1, &b[1], n, d__, (ftnlen)1);
#line 227 "slaqtr.f"
	xnorm = max(d__1,d__2);
#line 227 "slaqtr.f"
    }
/* Computing MAX */
#line 229 "slaqtr.f"
    d__1 = smlnum, d__2 = eps * xnorm;
#line 229 "slaqtr.f"
    smin = max(d__1,d__2);

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 234 "slaqtr.f"
    work[1] = 0.;
#line 235 "slaqtr.f"
    i__1 = *n;
#line 235 "slaqtr.f"
    for (j = 2; j <= i__1; ++j) {
#line 236 "slaqtr.f"
	i__2 = j - 1;
#line 236 "slaqtr.f"
	work[j] = sasum_(&i__2, &t[j * t_dim1 + 1], &c__1);
#line 237 "slaqtr.f"
/* L10: */
#line 237 "slaqtr.f"
    }

#line 239 "slaqtr.f"
    if (! (*lreal)) {
#line 240 "slaqtr.f"
	i__1 = *n;
#line 240 "slaqtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 241 "slaqtr.f"
	    work[i__] += (d__1 = b[i__], abs(d__1));
#line 242 "slaqtr.f"
/* L20: */
#line 242 "slaqtr.f"
	}
#line 243 "slaqtr.f"
    }

#line 245 "slaqtr.f"
    n2 = *n << 1;
#line 246 "slaqtr.f"
    n1 = *n;
#line 247 "slaqtr.f"
    if (! (*lreal)) {
#line 247 "slaqtr.f"
	n1 = n2;
#line 247 "slaqtr.f"
    }
#line 249 "slaqtr.f"
    k = isamax_(&n1, &x[1], &c__1);
#line 250 "slaqtr.f"
    xmax = (d__1 = x[k], abs(d__1));
#line 251 "slaqtr.f"
    *scale = 1.;

#line 253 "slaqtr.f"
    if (xmax > bignum) {
#line 254 "slaqtr.f"
	*scale = bignum / xmax;
#line 255 "slaqtr.f"
	sscal_(&n1, scale, &x[1], &c__1);
#line 256 "slaqtr.f"
	xmax = bignum;
#line 257 "slaqtr.f"
    }

#line 259 "slaqtr.f"
    if (*lreal) {

#line 261 "slaqtr.f"
	if (notran) {

/*           Solve T*p = scale*c */

#line 265 "slaqtr.f"
	    jnext = *n;
#line 266 "slaqtr.f"
	    for (j = *n; j >= 1; --j) {
#line 267 "slaqtr.f"
		if (j > jnext) {
#line 267 "slaqtr.f"
		    goto L30;
#line 267 "slaqtr.f"
		}
#line 269 "slaqtr.f"
		j1 = j;
#line 270 "slaqtr.f"
		j2 = j;
#line 271 "slaqtr.f"
		jnext = j - 1;
#line 272 "slaqtr.f"
		if (j > 1) {
#line 273 "slaqtr.f"
		    if (t[j + (j - 1) * t_dim1] != 0.) {
#line 274 "slaqtr.f"
			j1 = j - 1;
#line 275 "slaqtr.f"
			jnext = j - 2;
#line 276 "slaqtr.f"
		    }
#line 277 "slaqtr.f"
		}

#line 279 "slaqtr.f"
		if (j1 == j2) {

/*                 Meet 1 by 1 diagonal block */

/*                 Scale to avoid overflow when computing */
/*                     x(j) = b(j)/T(j,j) */

#line 286 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));
#line 287 "slaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1));
#line 288 "slaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 289 "slaqtr.f"
		    if (tjj < smin) {
#line 290 "slaqtr.f"
			tmp = smin;
#line 291 "slaqtr.f"
			tjj = smin;
#line 292 "slaqtr.f"
			*info = 1;
#line 293 "slaqtr.f"
		    }

#line 295 "slaqtr.f"
		    if (xj == 0.) {
#line 295 "slaqtr.f"
			goto L30;
#line 295 "slaqtr.f"
		    }

#line 298 "slaqtr.f"
		    if (tjj < 1.) {
#line 299 "slaqtr.f"
			if (xj > bignum * tjj) {
#line 300 "slaqtr.f"
			    rec = 1. / xj;
#line 301 "slaqtr.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 302 "slaqtr.f"
			    *scale *= rec;
#line 303 "slaqtr.f"
			    xmax *= rec;
#line 304 "slaqtr.f"
			}
#line 305 "slaqtr.f"
		    }
#line 306 "slaqtr.f"
		    x[j1] /= tmp;
#line 307 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));

/*                 Scale x if necessary to avoid overflow when adding a */
/*                 multiple of column j1 of T. */

#line 312 "slaqtr.f"
		    if (xj > 1.) {
#line 313 "slaqtr.f"
			rec = 1. / xj;
#line 314 "slaqtr.f"
			if (work[j1] > (bignum - xmax) * rec) {
#line 315 "slaqtr.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 316 "slaqtr.f"
			    *scale *= rec;
#line 317 "slaqtr.f"
			}
#line 318 "slaqtr.f"
		    }
#line 319 "slaqtr.f"
		    if (j1 > 1) {
#line 320 "slaqtr.f"
			i__1 = j1 - 1;
#line 320 "slaqtr.f"
			d__1 = -x[j1];
#line 320 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 321 "slaqtr.f"
			i__1 = j1 - 1;
#line 321 "slaqtr.f"
			k = isamax_(&i__1, &x[1], &c__1);
#line 322 "slaqtr.f"
			xmax = (d__1 = x[k], abs(d__1));
#line 323 "slaqtr.f"
		    }

#line 325 "slaqtr.f"
		} else {

/*                 Meet 2 by 2 diagonal block */

/*                 Call 2 by 2 linear system solve, to take */
/*                 care of possible overflow by scaling factor. */

#line 332 "slaqtr.f"
		    d__[0] = x[j1];
#line 333 "slaqtr.f"
		    d__[1] = x[j2];
#line 334 "slaqtr.f"
		    slaln2_(&c_false, &c__2, &c__1, &smin, &c_b21, &t[j1 + j1 
			    * t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &
			    c_b25, &c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
#line 337 "slaqtr.f"
		    if (ierr != 0) {
#line 337 "slaqtr.f"
			*info = 2;
#line 337 "slaqtr.f"
		    }

#line 340 "slaqtr.f"
		    if (scaloc != 1.) {
#line 341 "slaqtr.f"
			sscal_(n, &scaloc, &x[1], &c__1);
#line 342 "slaqtr.f"
			*scale *= scaloc;
#line 343 "slaqtr.f"
		    }
#line 344 "slaqtr.f"
		    x[j1] = v[0];
#line 345 "slaqtr.f"
		    x[j2] = v[1];

/*                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2)) */
/*                 to avoid overflow in updating right-hand side. */

/* Computing MAX */
#line 350 "slaqtr.f"
		    d__1 = abs(v[0]), d__2 = abs(v[1]);
#line 350 "slaqtr.f"
		    xj = max(d__1,d__2);
#line 351 "slaqtr.f"
		    if (xj > 1.) {
#line 352 "slaqtr.f"
			rec = 1. / xj;
/* Computing MAX */
#line 353 "slaqtr.f"
			d__1 = work[j1], d__2 = work[j2];
#line 353 "slaqtr.f"
			if (max(d__1,d__2) > (bignum - xmax) * rec) {
#line 355 "slaqtr.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 356 "slaqtr.f"
			    *scale *= rec;
#line 357 "slaqtr.f"
			}
#line 358 "slaqtr.f"
		    }

/*                 Update right-hand side */

#line 362 "slaqtr.f"
		    if (j1 > 1) {
#line 363 "slaqtr.f"
			i__1 = j1 - 1;
#line 363 "slaqtr.f"
			d__1 = -x[j1];
#line 363 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 364 "slaqtr.f"
			i__1 = j1 - 1;
#line 364 "slaqtr.f"
			d__1 = -x[j2];
#line 364 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 365 "slaqtr.f"
			i__1 = j1 - 1;
#line 365 "slaqtr.f"
			k = isamax_(&i__1, &x[1], &c__1);
#line 366 "slaqtr.f"
			xmax = (d__1 = x[k], abs(d__1));
#line 367 "slaqtr.f"
		    }

#line 369 "slaqtr.f"
		}

#line 371 "slaqtr.f"
L30:
#line 371 "slaqtr.f"
		;
#line 371 "slaqtr.f"
	    }

#line 373 "slaqtr.f"
	} else {

/*           Solve T**T*p = scale*c */

#line 377 "slaqtr.f"
	    jnext = 1;
#line 378 "slaqtr.f"
	    i__1 = *n;
#line 378 "slaqtr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 379 "slaqtr.f"
		if (j < jnext) {
#line 379 "slaqtr.f"
		    goto L40;
#line 379 "slaqtr.f"
		}
#line 381 "slaqtr.f"
		j1 = j;
#line 382 "slaqtr.f"
		j2 = j;
#line 383 "slaqtr.f"
		jnext = j + 1;
#line 384 "slaqtr.f"
		if (j < *n) {
#line 385 "slaqtr.f"
		    if (t[j + 1 + j * t_dim1] != 0.) {
#line 386 "slaqtr.f"
			j2 = j + 1;
#line 387 "slaqtr.f"
			jnext = j + 2;
#line 388 "slaqtr.f"
		    }
#line 389 "slaqtr.f"
		}

#line 391 "slaqtr.f"
		if (j1 == j2) {

/*                 1 by 1 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side element by inner product. */

#line 398 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));
#line 399 "slaqtr.f"
		    if (xmax > 1.) {
#line 400 "slaqtr.f"
			rec = 1. / xmax;
#line 401 "slaqtr.f"
			if (work[j1] > (bignum - xj) * rec) {
#line 402 "slaqtr.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 403 "slaqtr.f"
			    *scale *= rec;
#line 404 "slaqtr.f"
			    xmax *= rec;
#line 405 "slaqtr.f"
			}
#line 406 "slaqtr.f"
		    }

#line 408 "slaqtr.f"
		    i__2 = j1 - 1;
#line 408 "slaqtr.f"
		    x[j1] -= sdot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &
			    c__1);

#line 410 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1));
#line 411 "slaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1));
#line 412 "slaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 413 "slaqtr.f"
		    if (tjj < smin) {
#line 414 "slaqtr.f"
			tmp = smin;
#line 415 "slaqtr.f"
			tjj = smin;
#line 416 "slaqtr.f"
			*info = 1;
#line 417 "slaqtr.f"
		    }

#line 419 "slaqtr.f"
		    if (tjj < 1.) {
#line 420 "slaqtr.f"
			if (xj > bignum * tjj) {
#line 421 "slaqtr.f"
			    rec = 1. / xj;
#line 422 "slaqtr.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 423 "slaqtr.f"
			    *scale *= rec;
#line 424 "slaqtr.f"
			    xmax *= rec;
#line 425 "slaqtr.f"
			}
#line 426 "slaqtr.f"
		    }
#line 427 "slaqtr.f"
		    x[j1] /= tmp;
/* Computing MAX */
#line 428 "slaqtr.f"
		    d__2 = xmax, d__3 = (d__1 = x[j1], abs(d__1));
#line 428 "slaqtr.f"
		    xmax = max(d__2,d__3);

#line 430 "slaqtr.f"
		} else {

/*                 2 by 2 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side elements by inner product. */

/* Computing MAX */
#line 437 "slaqtr.f"
		    d__3 = (d__1 = x[j1], abs(d__1)), d__4 = (d__2 = x[j2], 
			    abs(d__2));
#line 437 "slaqtr.f"
		    xj = max(d__3,d__4);
#line 438 "slaqtr.f"
		    if (xmax > 1.) {
#line 439 "slaqtr.f"
			rec = 1. / xmax;
/* Computing MAX */
#line 440 "slaqtr.f"
			d__1 = work[j2], d__2 = work[j1];
#line 440 "slaqtr.f"
			if (max(d__1,d__2) > (bignum - xj) * rec) {
#line 442 "slaqtr.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 443 "slaqtr.f"
			    *scale *= rec;
#line 444 "slaqtr.f"
			    xmax *= rec;
#line 445 "slaqtr.f"
			}
#line 446 "slaqtr.f"
		    }

#line 448 "slaqtr.f"
		    i__2 = j1 - 1;
#line 448 "slaqtr.f"
		    d__[0] = x[j1] - sdot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);
#line 450 "slaqtr.f"
		    i__2 = j1 - 1;
#line 450 "slaqtr.f"
		    d__[1] = x[j2] - sdot_(&i__2, &t[j2 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);

#line 453 "slaqtr.f"
		    slaln2_(&c_true, &c__2, &c__1, &smin, &c_b21, &t[j1 + j1 *
			     t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &c_b25,
			     &c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
#line 456 "slaqtr.f"
		    if (ierr != 0) {
#line 456 "slaqtr.f"
			*info = 2;
#line 456 "slaqtr.f"
		    }

#line 459 "slaqtr.f"
		    if (scaloc != 1.) {
#line 460 "slaqtr.f"
			sscal_(n, &scaloc, &x[1], &c__1);
#line 461 "slaqtr.f"
			*scale *= scaloc;
#line 462 "slaqtr.f"
		    }
#line 463 "slaqtr.f"
		    x[j1] = v[0];
#line 464 "slaqtr.f"
		    x[j2] = v[1];
/* Computing MAX */
#line 465 "slaqtr.f"
		    d__3 = (d__1 = x[j1], abs(d__1)), d__4 = (d__2 = x[j2], 
			    abs(d__2)), d__3 = max(d__3,d__4);
#line 465 "slaqtr.f"
		    xmax = max(d__3,xmax);

#line 467 "slaqtr.f"
		}
#line 468 "slaqtr.f"
L40:
#line 468 "slaqtr.f"
		;
#line 468 "slaqtr.f"
	    }
#line 469 "slaqtr.f"
	}

#line 471 "slaqtr.f"
    } else {

/* Computing MAX */
#line 473 "slaqtr.f"
	d__1 = eps * abs(*w);
#line 473 "slaqtr.f"
	sminw = max(d__1,smin);
#line 474 "slaqtr.f"
	if (notran) {

/*           Solve (T + iB)*(p+iq) = c+id */

#line 478 "slaqtr.f"
	    jnext = *n;
#line 479 "slaqtr.f"
	    for (j = *n; j >= 1; --j) {
#line 480 "slaqtr.f"
		if (j > jnext) {
#line 480 "slaqtr.f"
		    goto L70;
#line 480 "slaqtr.f"
		}
#line 482 "slaqtr.f"
		j1 = j;
#line 483 "slaqtr.f"
		j2 = j;
#line 484 "slaqtr.f"
		jnext = j - 1;
#line 485 "slaqtr.f"
		if (j > 1) {
#line 486 "slaqtr.f"
		    if (t[j + (j - 1) * t_dim1] != 0.) {
#line 487 "slaqtr.f"
			j1 = j - 1;
#line 488 "slaqtr.f"
			jnext = j - 2;
#line 489 "slaqtr.f"
		    }
#line 490 "slaqtr.f"
		}

#line 492 "slaqtr.f"
		if (j1 == j2) {

/*                 1 by 1 diagonal block */

/*                 Scale if necessary to avoid overflow in division */

#line 498 "slaqtr.f"
		    z__ = *w;
#line 499 "slaqtr.f"
		    if (j1 == 1) {
#line 499 "slaqtr.f"
			z__ = b[1];
#line 499 "slaqtr.f"
		    }
#line 501 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], abs(
			    d__2));
#line 502 "slaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1)) + abs(z__);
#line 503 "slaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 504 "slaqtr.f"
		    if (tjj < sminw) {
#line 505 "slaqtr.f"
			tmp = sminw;
#line 506 "slaqtr.f"
			tjj = sminw;
#line 507 "slaqtr.f"
			*info = 1;
#line 508 "slaqtr.f"
		    }

#line 510 "slaqtr.f"
		    if (xj == 0.) {
#line 510 "slaqtr.f"
			goto L70;
#line 510 "slaqtr.f"
		    }

#line 513 "slaqtr.f"
		    if (tjj < 1.) {
#line 514 "slaqtr.f"
			if (xj > bignum * tjj) {
#line 515 "slaqtr.f"
			    rec = 1. / xj;
#line 516 "slaqtr.f"
			    sscal_(&n2, &rec, &x[1], &c__1);
#line 517 "slaqtr.f"
			    *scale *= rec;
#line 518 "slaqtr.f"
			    xmax *= rec;
#line 519 "slaqtr.f"
			}
#line 520 "slaqtr.f"
		    }
#line 521 "slaqtr.f"
		    sladiv_(&x[j1], &x[*n + j1], &tmp, &z__, &sr, &si);
#line 522 "slaqtr.f"
		    x[j1] = sr;
#line 523 "slaqtr.f"
		    x[*n + j1] = si;
#line 524 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], abs(
			    d__2));

/*                 Scale x if necessary to avoid overflow when adding a */
/*                 multiple of column j1 of T. */

#line 529 "slaqtr.f"
		    if (xj > 1.) {
#line 530 "slaqtr.f"
			rec = 1. / xj;
#line 531 "slaqtr.f"
			if (work[j1] > (bignum - xmax) * rec) {
#line 532 "slaqtr.f"
			    sscal_(&n2, &rec, &x[1], &c__1);
#line 533 "slaqtr.f"
			    *scale *= rec;
#line 534 "slaqtr.f"
			}
#line 535 "slaqtr.f"
		    }

#line 537 "slaqtr.f"
		    if (j1 > 1) {
#line 538 "slaqtr.f"
			i__1 = j1 - 1;
#line 538 "slaqtr.f"
			d__1 = -x[j1];
#line 538 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 539 "slaqtr.f"
			i__1 = j1 - 1;
#line 539 "slaqtr.f"
			d__1 = -x[*n + j1];
#line 539 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[*
				n + 1], &c__1);

#line 542 "slaqtr.f"
			x[1] += b[j1] * x[*n + j1];
#line 543 "slaqtr.f"
			x[*n + 1] -= b[j1] * x[j1];

#line 545 "slaqtr.f"
			xmax = 0.;
#line 546 "slaqtr.f"
			i__1 = j1 - 1;
#line 546 "slaqtr.f"
			for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 547 "slaqtr.f"
			    d__3 = xmax, d__4 = (d__1 = x[k], abs(d__1)) + (
				    d__2 = x[k + *n], abs(d__2));
#line 547 "slaqtr.f"
			    xmax = max(d__3,d__4);
#line 549 "slaqtr.f"
/* L50: */
#line 549 "slaqtr.f"
			}
#line 550 "slaqtr.f"
		    }

#line 552 "slaqtr.f"
		} else {

/*                 Meet 2 by 2 diagonal block */

#line 556 "slaqtr.f"
		    d__[0] = x[j1];
#line 557 "slaqtr.f"
		    d__[1] = x[j2];
#line 558 "slaqtr.f"
		    d__[2] = x[*n + j1];
#line 559 "slaqtr.f"
		    d__[3] = x[*n + j2];
#line 560 "slaqtr.f"
		    d__1 = -(*w);
#line 560 "slaqtr.f"
		    slaln2_(&c_false, &c__2, &c__2, &sminw, &c_b21, &t[j1 + 
			    j1 * t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &
			    c_b25, &d__1, v, &c__2, &scaloc, &xnorm, &ierr);
#line 563 "slaqtr.f"
		    if (ierr != 0) {
#line 563 "slaqtr.f"
			*info = 2;
#line 563 "slaqtr.f"
		    }

#line 566 "slaqtr.f"
		    if (scaloc != 1.) {
#line 567 "slaqtr.f"
			i__1 = *n << 1;
#line 567 "slaqtr.f"
			sscal_(&i__1, &scaloc, &x[1], &c__1);
#line 568 "slaqtr.f"
			*scale = scaloc * *scale;
#line 569 "slaqtr.f"
		    }
#line 570 "slaqtr.f"
		    x[j1] = v[0];
#line 571 "slaqtr.f"
		    x[j2] = v[1];
#line 572 "slaqtr.f"
		    x[*n + j1] = v[2];
#line 573 "slaqtr.f"
		    x[*n + j2] = v[3];

/*                 Scale X(J1), .... to avoid overflow in */
/*                 updating right hand side. */

/* Computing MAX */
#line 578 "slaqtr.f"
		    d__1 = abs(v[0]) + abs(v[2]), d__2 = abs(v[1]) + abs(v[3])
			    ;
#line 578 "slaqtr.f"
		    xj = max(d__1,d__2);
#line 580 "slaqtr.f"
		    if (xj > 1.) {
#line 581 "slaqtr.f"
			rec = 1. / xj;
/* Computing MAX */
#line 582 "slaqtr.f"
			d__1 = work[j1], d__2 = work[j2];
#line 582 "slaqtr.f"
			if (max(d__1,d__2) > (bignum - xmax) * rec) {
#line 584 "slaqtr.f"
			    sscal_(&n2, &rec, &x[1], &c__1);
#line 585 "slaqtr.f"
			    *scale *= rec;
#line 586 "slaqtr.f"
			}
#line 587 "slaqtr.f"
		    }

/*                 Update the right-hand side. */

#line 591 "slaqtr.f"
		    if (j1 > 1) {
#line 592 "slaqtr.f"
			i__1 = j1 - 1;
#line 592 "slaqtr.f"
			d__1 = -x[j1];
#line 592 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);
#line 593 "slaqtr.f"
			i__1 = j1 - 1;
#line 593 "slaqtr.f"
			d__1 = -x[j2];
#line 593 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[1]
				, &c__1);

#line 595 "slaqtr.f"
			i__1 = j1 - 1;
#line 595 "slaqtr.f"
			d__1 = -x[*n + j1];
#line 595 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[*
				n + 1], &c__1);
#line 597 "slaqtr.f"
			i__1 = j1 - 1;
#line 597 "slaqtr.f"
			d__1 = -x[*n + j2];
#line 597 "slaqtr.f"
			saxpy_(&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[*
				n + 1], &c__1);

#line 600 "slaqtr.f"
			x[1] = x[1] + b[j1] * x[*n + j1] + b[j2] * x[*n + j2];
#line 602 "slaqtr.f"
			x[*n + 1] = x[*n + 1] - b[j1] * x[j1] - b[j2] * x[j2];

#line 605 "slaqtr.f"
			xmax = 0.;
#line 606 "slaqtr.f"
			i__1 = j1 - 1;
#line 606 "slaqtr.f"
			for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 607 "slaqtr.f"
			    d__3 = (d__1 = x[k], abs(d__1)) + (d__2 = x[k + *
				    n], abs(d__2));
#line 607 "slaqtr.f"
			    xmax = max(d__3,xmax);
#line 609 "slaqtr.f"
/* L60: */
#line 609 "slaqtr.f"
			}
#line 610 "slaqtr.f"
		    }

#line 612 "slaqtr.f"
		}
#line 613 "slaqtr.f"
L70:
#line 613 "slaqtr.f"
		;
#line 613 "slaqtr.f"
	    }

#line 615 "slaqtr.f"
	} else {

/*           Solve (T + iB)**T*(p+iq) = c+id */

#line 619 "slaqtr.f"
	    jnext = 1;
#line 620 "slaqtr.f"
	    i__1 = *n;
#line 620 "slaqtr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 621 "slaqtr.f"
		if (j < jnext) {
#line 621 "slaqtr.f"
		    goto L80;
#line 621 "slaqtr.f"
		}
#line 623 "slaqtr.f"
		j1 = j;
#line 624 "slaqtr.f"
		j2 = j;
#line 625 "slaqtr.f"
		jnext = j + 1;
#line 626 "slaqtr.f"
		if (j < *n) {
#line 627 "slaqtr.f"
		    if (t[j + 1 + j * t_dim1] != 0.) {
#line 628 "slaqtr.f"
			j2 = j + 1;
#line 629 "slaqtr.f"
			jnext = j + 2;
#line 630 "slaqtr.f"
		    }
#line 631 "slaqtr.f"
		}

#line 633 "slaqtr.f"
		if (j1 == j2) {

/*                 1 by 1 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side element by inner product. */

#line 640 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[j1 + *n], abs(
			    d__2));
#line 641 "slaqtr.f"
		    if (xmax > 1.) {
#line 642 "slaqtr.f"
			rec = 1. / xmax;
#line 643 "slaqtr.f"
			if (work[j1] > (bignum - xj) * rec) {
#line 644 "slaqtr.f"
			    sscal_(&n2, &rec, &x[1], &c__1);
#line 645 "slaqtr.f"
			    *scale *= rec;
#line 646 "slaqtr.f"
			    xmax *= rec;
#line 647 "slaqtr.f"
			}
#line 648 "slaqtr.f"
		    }

#line 650 "slaqtr.f"
		    i__2 = j1 - 1;
#line 650 "slaqtr.f"
		    x[j1] -= sdot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &
			    c__1);
#line 651 "slaqtr.f"
		    i__2 = j1 - 1;
#line 651 "slaqtr.f"
		    x[*n + j1] -= sdot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[
			    *n + 1], &c__1);
#line 653 "slaqtr.f"
		    if (j1 > 1) {
#line 654 "slaqtr.f"
			x[j1] -= b[j1] * x[*n + 1];
#line 655 "slaqtr.f"
			x[*n + j1] += b[j1] * x[1];
#line 656 "slaqtr.f"
		    }
#line 657 "slaqtr.f"
		    xj = (d__1 = x[j1], abs(d__1)) + (d__2 = x[j1 + *n], abs(
			    d__2));

#line 659 "slaqtr.f"
		    z__ = *w;
#line 660 "slaqtr.f"
		    if (j1 == 1) {
#line 660 "slaqtr.f"
			z__ = b[1];
#line 660 "slaqtr.f"
		    }

/*                 Scale if necessary to avoid overflow in */
/*                 complex division */

#line 666 "slaqtr.f"
		    tjj = (d__1 = t[j1 + j1 * t_dim1], abs(d__1)) + abs(z__);
#line 667 "slaqtr.f"
		    tmp = t[j1 + j1 * t_dim1];
#line 668 "slaqtr.f"
		    if (tjj < sminw) {
#line 669 "slaqtr.f"
			tmp = sminw;
#line 670 "slaqtr.f"
			tjj = sminw;
#line 671 "slaqtr.f"
			*info = 1;
#line 672 "slaqtr.f"
		    }

#line 674 "slaqtr.f"
		    if (tjj < 1.) {
#line 675 "slaqtr.f"
			if (xj > bignum * tjj) {
#line 676 "slaqtr.f"
			    rec = 1. / xj;
#line 677 "slaqtr.f"
			    sscal_(&n2, &rec, &x[1], &c__1);
#line 678 "slaqtr.f"
			    *scale *= rec;
#line 679 "slaqtr.f"
			    xmax *= rec;
#line 680 "slaqtr.f"
			}
#line 681 "slaqtr.f"
		    }
#line 682 "slaqtr.f"
		    d__1 = -z__;
#line 682 "slaqtr.f"
		    sladiv_(&x[j1], &x[*n + j1], &tmp, &d__1, &sr, &si);
#line 683 "slaqtr.f"
		    x[j1] = sr;
#line 684 "slaqtr.f"
		    x[j1 + *n] = si;
/* Computing MAX */
#line 685 "slaqtr.f"
		    d__3 = (d__1 = x[j1], abs(d__1)) + (d__2 = x[j1 + *n], 
			    abs(d__2));
#line 685 "slaqtr.f"
		    xmax = max(d__3,xmax);

#line 687 "slaqtr.f"
		} else {

/*                 2 by 2 diagonal block */

/*                 Scale if necessary to avoid overflow in forming the */
/*                 right-hand side element by inner product. */

/* Computing MAX */
#line 694 "slaqtr.f"
		    d__5 = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], 
			    abs(d__2)), d__6 = (d__3 = x[j2], abs(d__3)) + (
			    d__4 = x[*n + j2], abs(d__4));
#line 694 "slaqtr.f"
		    xj = max(d__5,d__6);
#line 696 "slaqtr.f"
		    if (xmax > 1.) {
#line 697 "slaqtr.f"
			rec = 1. / xmax;
/* Computing MAX */
#line 698 "slaqtr.f"
			d__1 = work[j1], d__2 = work[j2];
#line 698 "slaqtr.f"
			if (max(d__1,d__2) > (bignum - xj) / xmax) {
#line 700 "slaqtr.f"
			    sscal_(&n2, &rec, &x[1], &c__1);
#line 701 "slaqtr.f"
			    *scale *= rec;
#line 702 "slaqtr.f"
			    xmax *= rec;
#line 703 "slaqtr.f"
			}
#line 704 "slaqtr.f"
		    }

#line 706 "slaqtr.f"
		    i__2 = j1 - 1;
#line 706 "slaqtr.f"
		    d__[0] = x[j1] - sdot_(&i__2, &t[j1 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);
#line 708 "slaqtr.f"
		    i__2 = j1 - 1;
#line 708 "slaqtr.f"
		    d__[1] = x[j2] - sdot_(&i__2, &t[j2 * t_dim1 + 1], &c__1, 
			    &x[1], &c__1);
#line 710 "slaqtr.f"
		    i__2 = j1 - 1;
#line 710 "slaqtr.f"
		    d__[2] = x[*n + j1] - sdot_(&i__2, &t[j1 * t_dim1 + 1], &
			    c__1, &x[*n + 1], &c__1);
#line 712 "slaqtr.f"
		    i__2 = j1 - 1;
#line 712 "slaqtr.f"
		    d__[3] = x[*n + j2] - sdot_(&i__2, &t[j2 * t_dim1 + 1], &
			    c__1, &x[*n + 1], &c__1);
#line 714 "slaqtr.f"
		    d__[0] -= b[j1] * x[*n + 1];
#line 715 "slaqtr.f"
		    d__[1] -= b[j2] * x[*n + 1];
#line 716 "slaqtr.f"
		    d__[2] += b[j1] * x[1];
#line 717 "slaqtr.f"
		    d__[3] += b[j2] * x[1];

#line 719 "slaqtr.f"
		    slaln2_(&c_true, &c__2, &c__2, &sminw, &c_b21, &t[j1 + j1 
			    * t_dim1], ldt, &c_b21, &c_b21, d__, &c__2, &
			    c_b25, w, v, &c__2, &scaloc, &xnorm, &ierr);
#line 722 "slaqtr.f"
		    if (ierr != 0) {
#line 722 "slaqtr.f"
			*info = 2;
#line 722 "slaqtr.f"
		    }

#line 725 "slaqtr.f"
		    if (scaloc != 1.) {
#line 726 "slaqtr.f"
			sscal_(&n2, &scaloc, &x[1], &c__1);
#line 727 "slaqtr.f"
			*scale = scaloc * *scale;
#line 728 "slaqtr.f"
		    }
#line 729 "slaqtr.f"
		    x[j1] = v[0];
#line 730 "slaqtr.f"
		    x[j2] = v[1];
#line 731 "slaqtr.f"
		    x[*n + j1] = v[2];
#line 732 "slaqtr.f"
		    x[*n + j2] = v[3];
/* Computing MAX */
#line 733 "slaqtr.f"
		    d__5 = (d__1 = x[j1], abs(d__1)) + (d__2 = x[*n + j1], 
			    abs(d__2)), d__6 = (d__3 = x[j2], abs(d__3)) + (
			    d__4 = x[*n + j2], abs(d__4)), d__5 = max(d__5,
			    d__6);
#line 733 "slaqtr.f"
		    xmax = max(d__5,xmax);

#line 736 "slaqtr.f"
		}

#line 738 "slaqtr.f"
L80:
#line 738 "slaqtr.f"
		;
#line 738 "slaqtr.f"
	    }

#line 740 "slaqtr.f"
	}

#line 742 "slaqtr.f"
    }

#line 744 "slaqtr.f"
    return 0;

/*     End of SLAQTR */

} /* slaqtr_ */

