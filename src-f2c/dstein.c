#line 1 "dstein.f"
/* dstein.f -- translated by f2c (version 20100827).
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

#line 1 "dstein.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b DSTEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDZ, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/*      $                   IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEIN computes the eigenvectors of a real symmetric tridiagonal */
/* > matrix T corresponding to specified eigenvalues, using inverse */
/* > iteration. */
/* > */
/* > The maximum number of iterations allowed for each eigenvector is */
/* > specified by an internal parameter MAXITS (currently set to 5). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of the tridiagonal matrix */
/* >          T, in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of eigenvectors to be found.  0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          The first M elements of W contain the eigenvalues for */
/* >          which eigenvectors are to be computed.  The eigenvalues */
/* >          should be grouped by split-off block and ordered from */
/* >          smallest to largest within the block.  ( The output array */
/* >          W from DSTEBZ with ORDER = 'B' is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] IBLOCK */
/* > \verbatim */
/* >          IBLOCK is INTEGER array, dimension (N) */
/* >          The submatrix indices associated with the corresponding */
/* >          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to */
/* >          the first submatrix from the top, =2 if W(i) belongs to */
/* >          the second submatrix, etc.  ( The output array IBLOCK */
/* >          from DSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] ISPLIT */
/* > \verbatim */
/* >          ISPLIT is INTEGER array, dimension (N) */
/* >          The splitting points, at which T breaks up into submatrices. */
/* >          The first submatrix consists of rows/columns 1 to */
/* >          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1 */
/* >          through ISPLIT( 2 ), etc. */
/* >          ( The output array ISPLIT from DSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, M) */
/* >          The computed eigenvectors.  The eigenvector associated */
/* >          with the eigenvalue W(i) is stored in the i-th column of */
/* >          Z.  Any vector which fails to converge is set to its current */
/* >          iterate after MAXITS iterations. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* >          IFAIL is INTEGER array, dimension (M) */
/* >          On normal exit, all elements of IFAIL are zero. */
/* >          If one or more eigenvectors fail to converge after */
/* >          MAXITS iterations, then their indices are stored in */
/* >          array IFAIL. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit. */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, then i eigenvectors failed to converge */
/* >               in MAXITS iterations.  Their indices are stored in */
/* >               array IFAIL. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  MAXITS  INTEGER, default = 5 */
/* >          The maximum number of iterations performed. */
/* > */
/* >  EXTRA   INTEGER, default = 2 */
/* >          The number of iterations performed after norm growth */
/* >          criterion is satisfied, should be at least 1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dstein_(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, b1, j1, bn;
    static doublereal xj, scl, eps, sep, nrm, tol;
    static integer its;
    static doublereal xjm, ztr, eps1;
    static integer jblk, nblk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iseed[4], gpind, iinfo;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal ortol;
    static integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlagtf_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlagts_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer nrmchk;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    static integer blksiz;
    static doublereal onenrm, dtpcrt, pertol;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

/*     Test the input parameters. */

#line 226 "dstein.f"
    /* Parameter adjustments */
#line 226 "dstein.f"
    --d__;
#line 226 "dstein.f"
    --e;
#line 226 "dstein.f"
    --w;
#line 226 "dstein.f"
    --iblock;
#line 226 "dstein.f"
    --isplit;
#line 226 "dstein.f"
    z_dim1 = *ldz;
#line 226 "dstein.f"
    z_offset = 1 + z_dim1;
#line 226 "dstein.f"
    z__ -= z_offset;
#line 226 "dstein.f"
    --work;
#line 226 "dstein.f"
    --iwork;
#line 226 "dstein.f"
    --ifail;
#line 226 "dstein.f"

#line 226 "dstein.f"
    /* Function Body */
#line 226 "dstein.f"
    *info = 0;
#line 227 "dstein.f"
    i__1 = *m;
#line 227 "dstein.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 228 "dstein.f"
	ifail[i__] = 0;
#line 229 "dstein.f"
/* L10: */
#line 229 "dstein.f"
    }

#line 231 "dstein.f"
    if (*n < 0) {
#line 232 "dstein.f"
	*info = -1;
#line 233 "dstein.f"
    } else if (*m < 0 || *m > *n) {
#line 234 "dstein.f"
	*info = -4;
#line 235 "dstein.f"
    } else if (*ldz < max(1,*n)) {
#line 236 "dstein.f"
	*info = -9;
#line 237 "dstein.f"
    } else {
#line 238 "dstein.f"
	i__1 = *m;
#line 238 "dstein.f"
	for (j = 2; j <= i__1; ++j) {
#line 239 "dstein.f"
	    if (iblock[j] < iblock[j - 1]) {
#line 240 "dstein.f"
		*info = -6;
#line 241 "dstein.f"
		goto L30;
#line 242 "dstein.f"
	    }
#line 243 "dstein.f"
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
#line 245 "dstein.f"
		*info = -5;
#line 246 "dstein.f"
		goto L30;
#line 247 "dstein.f"
	    }
#line 248 "dstein.f"
/* L20: */
#line 248 "dstein.f"
	}
#line 249 "dstein.f"
L30:
#line 250 "dstein.f"
	;
#line 250 "dstein.f"
    }

#line 252 "dstein.f"
    if (*info != 0) {
#line 253 "dstein.f"
	i__1 = -(*info);
#line 253 "dstein.f"
	xerbla_("DSTEIN", &i__1, (ftnlen)6);
#line 254 "dstein.f"
	return 0;
#line 255 "dstein.f"
    }

/*     Quick return if possible */

#line 259 "dstein.f"
    if (*n == 0 || *m == 0) {
#line 260 "dstein.f"
	return 0;
#line 261 "dstein.f"
    } else if (*n == 1) {
#line 262 "dstein.f"
	z__[z_dim1 + 1] = 1.;
#line 263 "dstein.f"
	return 0;
#line 264 "dstein.f"
    }

/*     Get machine constants. */

#line 268 "dstein.f"
    eps = dlamch_("Precision", (ftnlen)9);

/*     Initialize seed for random number generator DLARNV. */

#line 272 "dstein.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 273 "dstein.f"
	iseed[i__ - 1] = 1;
#line 274 "dstein.f"
/* L40: */
#line 274 "dstein.f"
    }

/*     Initialize pointers. */

#line 278 "dstein.f"
    indrv1 = 0;
#line 279 "dstein.f"
    indrv2 = indrv1 + *n;
#line 280 "dstein.f"
    indrv3 = indrv2 + *n;
#line 281 "dstein.f"
    indrv4 = indrv3 + *n;
#line 282 "dstein.f"
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

#line 286 "dstein.f"
    j1 = 1;
#line 287 "dstein.f"
    i__1 = iblock[*m];
#line 287 "dstein.f"
    for (nblk = 1; nblk <= i__1; ++nblk) {

/*        Find starting and ending indices of block nblk. */

#line 291 "dstein.f"
	if (nblk == 1) {
#line 292 "dstein.f"
	    b1 = 1;
#line 293 "dstein.f"
	} else {
#line 294 "dstein.f"
	    b1 = isplit[nblk - 1] + 1;
#line 295 "dstein.f"
	}
#line 296 "dstein.f"
	bn = isplit[nblk];
#line 297 "dstein.f"
	blksiz = bn - b1 + 1;
#line 298 "dstein.f"
	if (blksiz == 1) {
#line 298 "dstein.f"
	    goto L60;
#line 298 "dstein.f"
	}
#line 300 "dstein.f"
	gpind = j1;

/*        Compute reorthogonalization criterion and stopping criterion. */

#line 304 "dstein.f"
	onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
/* Computing MAX */
#line 305 "dstein.f"
	d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1],
		 abs(d__2));
#line 305 "dstein.f"
	onenrm = max(d__3,d__4);
#line 306 "dstein.f"
	i__2 = bn - 1;
#line 306 "dstein.f"
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 307 "dstein.f"
	    d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
		    i__ - 1], abs(d__2)) + (d__3 = e[i__], abs(d__3));
#line 307 "dstein.f"
	    onenrm = max(d__4,d__5);
#line 309 "dstein.f"
/* L50: */
#line 309 "dstein.f"
	}
#line 310 "dstein.f"
	ortol = onenrm * .001;

#line 312 "dstein.f"
	dtpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

#line 316 "dstein.f"
L60:
#line 317 "dstein.f"
	jblk = 0;
#line 318 "dstein.f"
	i__2 = *m;
#line 318 "dstein.f"
	for (j = j1; j <= i__2; ++j) {
#line 319 "dstein.f"
	    if (iblock[j] != nblk) {
#line 320 "dstein.f"
		j1 = j;
#line 321 "dstein.f"
		goto L160;
#line 322 "dstein.f"
	    }
#line 323 "dstein.f"
	    ++jblk;
#line 324 "dstein.f"
	    xj = w[j];

/*           Skip all the work if the block size is one. */

#line 328 "dstein.f"
	    if (blksiz == 1) {
#line 329 "dstein.f"
		work[indrv1 + 1] = 1.;
#line 330 "dstein.f"
		goto L120;
#line 331 "dstein.f"
	    }

/*           If eigenvalues j and j-1 are too close, add a relatively */
/*           small perturbation. */

#line 336 "dstein.f"
	    if (jblk > 1) {
#line 337 "dstein.f"
		eps1 = (d__1 = eps * xj, abs(d__1));
#line 338 "dstein.f"
		pertol = eps1 * 10.;
#line 339 "dstein.f"
		sep = xj - xjm;
#line 340 "dstein.f"
		if (sep < pertol) {
#line 340 "dstein.f"
		    xj = xjm + pertol;
#line 340 "dstein.f"
		}
#line 342 "dstein.f"
	    }

#line 344 "dstein.f"
	    its = 0;
#line 345 "dstein.f"
	    nrmchk = 0;

/*           Get random starting vector. */

#line 349 "dstein.f"
	    dlarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

/*           Copy the matrix T so it won't be destroyed in factorization. */

#line 353 "dstein.f"
	    dcopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
#line 354 "dstein.f"
	    i__3 = blksiz - 1;
#line 354 "dstein.f"
	    dcopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
#line 355 "dstein.f"
	    i__3 = blksiz - 1;
#line 355 "dstein.f"
	    dcopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU ) */

#line 359 "dstein.f"
	    tol = 0.;
#line 360 "dstein.f"
	    dlagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

/*           Update iteration count. */

#line 366 "dstein.f"
L70:
#line 367 "dstein.f"
	    ++its;
#line 368 "dstein.f"
	    if (its > 5) {
#line 368 "dstein.f"
		goto L100;
#line 368 "dstein.f"
	    }

/*           Normalize and scale the righthand side vector Pb. */

#line 373 "dstein.f"
	    jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
/* Computing MAX */
#line 374 "dstein.f"
	    d__3 = eps, d__4 = (d__1 = work[indrv4 + blksiz], abs(d__1));
#line 374 "dstein.f"
	    scl = blksiz * onenrm * max(d__3,d__4) / (d__2 = work[indrv1 + 
		    jmax], abs(d__2));
#line 377 "dstein.f"
	    dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);

/*           Solve the system LU = Pb. */

#line 381 "dstein.f"
	    dlagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
/*           close enough. */

#line 388 "dstein.f"
	    if (jblk == 1) {
#line 388 "dstein.f"
		goto L90;
#line 388 "dstein.f"
	    }
#line 390 "dstein.f"
	    if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
#line 390 "dstein.f"
		gpind = j;
#line 390 "dstein.f"
	    }
#line 392 "dstein.f"
	    if (gpind != j) {
#line 393 "dstein.f"
		i__3 = j - 1;
#line 393 "dstein.f"
		for (i__ = gpind; i__ <= i__3; ++i__) {
#line 394 "dstein.f"
		    ztr = -ddot_(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + 
			    i__ * z_dim1], &c__1);
#line 396 "dstein.f"
		    daxpy_(&blksiz, &ztr, &z__[b1 + i__ * z_dim1], &c__1, &
			    work[indrv1 + 1], &c__1);
#line 398 "dstein.f"
/* L80: */
#line 398 "dstein.f"
		}
#line 399 "dstein.f"
	    }

/*           Check the infinity norm of the iterate. */

#line 403 "dstein.f"
L90:
#line 404 "dstein.f"
	    jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 405 "dstein.f"
	    nrm = (d__1 = work[indrv1 + jmax], abs(d__1));

/*           Continue for additional iterations after norm reaches */
/*           stopping criterion. */

#line 410 "dstein.f"
	    if (nrm < dtpcrt) {
#line 410 "dstein.f"
		goto L70;
#line 410 "dstein.f"
	    }
#line 412 "dstein.f"
	    ++nrmchk;
#line 413 "dstein.f"
	    if (nrmchk < 3) {
#line 413 "dstein.f"
		goto L70;
#line 413 "dstein.f"
	    }

#line 416 "dstein.f"
	    goto L110;

/*           If stopping criterion was not satisfied, update info and */
/*           store eigenvector number in array ifail. */

#line 421 "dstein.f"
L100:
#line 422 "dstein.f"
	    ++(*info);
#line 423 "dstein.f"
	    ifail[*info] = j;

/*           Accept iterate as jth eigenvector. */

#line 427 "dstein.f"
L110:
#line 428 "dstein.f"
	    scl = 1. / dnrm2_(&blksiz, &work[indrv1 + 1], &c__1);
#line 429 "dstein.f"
	    jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 430 "dstein.f"
	    if (work[indrv1 + jmax] < 0.) {
#line 430 "dstein.f"
		scl = -scl;
#line 430 "dstein.f"
	    }
#line 432 "dstein.f"
	    dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
#line 433 "dstein.f"
L120:
#line 434 "dstein.f"
	    i__3 = *n;
#line 434 "dstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 435 "dstein.f"
		z__[i__ + j * z_dim1] = 0.;
#line 436 "dstein.f"
/* L130: */
#line 436 "dstein.f"
	    }
#line 437 "dstein.f"
	    i__3 = blksiz;
#line 437 "dstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 438 "dstein.f"
		z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
#line 439 "dstein.f"
/* L140: */
#line 439 "dstein.f"
	    }

/*           Save the shift to check eigenvalue spacing at next */
/*           iteration. */

#line 444 "dstein.f"
	    xjm = xj;

#line 446 "dstein.f"
/* L150: */
#line 446 "dstein.f"
	}
#line 447 "dstein.f"
L160:
#line 447 "dstein.f"
	;
#line 447 "dstein.f"
    }

#line 449 "dstein.f"
    return 0;

/*     End of DSTEIN */

} /* dstein_ */

