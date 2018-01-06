#line 1 "zstein.f"
/* zstein.f -- translated by f2c (version 20100827).
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

#line 1 "zstein.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZSTEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSTEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDZ, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/*      $                   IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ) */
/*       COMPLEX*16         Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSTEIN computes the eigenvectors of a real symmetric tridiagonal */
/* > matrix T corresponding to specified eigenvalues, using inverse */
/* > iteration. */
/* > */
/* > The maximum number of iterations allowed for each eigenvector is */
/* > specified by an internal parameter MAXITS (currently set to 5). */
/* > */
/* > Although the eigenvectors are real, they are stored in a complex */
/* > array, which may be passed to ZUNMTR or ZUPMTR for back */
/* > transformation to the eigenvectors of a complex Hermitian matrix */
/* > which was reduced to tridiagonal form. */
/* > */
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
/* >          T, stored in elements 1 to N-1. */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, M) */
/* >          The computed eigenvectors.  The eigenvector associated */
/* >          with the eigenvalue W(i) is stored in the i-th column of */
/* >          Z.  Any vector which fails to converge is set to its current */
/* >          iterate after MAXITS iterations. */
/* >          The imaginary parts of the eigenvectors are set to zero. */
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
/* >          = 0: successful exit */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zstein_(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, b1, j1, bn, jr;
    static doublereal xj, scl, eps, sep, nrm, tol;
    static integer its;
    static doublereal xjm, ztr, eps1;
    static integer jblk, nblk, jmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iseed[4], gpind, iinfo;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
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


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

/*     Test the input parameters. */

#line 237 "zstein.f"
    /* Parameter adjustments */
#line 237 "zstein.f"
    --d__;
#line 237 "zstein.f"
    --e;
#line 237 "zstein.f"
    --w;
#line 237 "zstein.f"
    --iblock;
#line 237 "zstein.f"
    --isplit;
#line 237 "zstein.f"
    z_dim1 = *ldz;
#line 237 "zstein.f"
    z_offset = 1 + z_dim1;
#line 237 "zstein.f"
    z__ -= z_offset;
#line 237 "zstein.f"
    --work;
#line 237 "zstein.f"
    --iwork;
#line 237 "zstein.f"
    --ifail;
#line 237 "zstein.f"

#line 237 "zstein.f"
    /* Function Body */
#line 237 "zstein.f"
    *info = 0;
#line 238 "zstein.f"
    i__1 = *m;
#line 238 "zstein.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 239 "zstein.f"
	ifail[i__] = 0;
#line 240 "zstein.f"
/* L10: */
#line 240 "zstein.f"
    }

#line 242 "zstein.f"
    if (*n < 0) {
#line 243 "zstein.f"
	*info = -1;
#line 244 "zstein.f"
    } else if (*m < 0 || *m > *n) {
#line 245 "zstein.f"
	*info = -4;
#line 246 "zstein.f"
    } else if (*ldz < max(1,*n)) {
#line 247 "zstein.f"
	*info = -9;
#line 248 "zstein.f"
    } else {
#line 249 "zstein.f"
	i__1 = *m;
#line 249 "zstein.f"
	for (j = 2; j <= i__1; ++j) {
#line 250 "zstein.f"
	    if (iblock[j] < iblock[j - 1]) {
#line 251 "zstein.f"
		*info = -6;
#line 252 "zstein.f"
		goto L30;
#line 253 "zstein.f"
	    }
#line 254 "zstein.f"
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
#line 256 "zstein.f"
		*info = -5;
#line 257 "zstein.f"
		goto L30;
#line 258 "zstein.f"
	    }
#line 259 "zstein.f"
/* L20: */
#line 259 "zstein.f"
	}
#line 260 "zstein.f"
L30:
#line 261 "zstein.f"
	;
#line 261 "zstein.f"
    }

#line 263 "zstein.f"
    if (*info != 0) {
#line 264 "zstein.f"
	i__1 = -(*info);
#line 264 "zstein.f"
	xerbla_("ZSTEIN", &i__1, (ftnlen)6);
#line 265 "zstein.f"
	return 0;
#line 266 "zstein.f"
    }

/*     Quick return if possible */

#line 270 "zstein.f"
    if (*n == 0 || *m == 0) {
#line 271 "zstein.f"
	return 0;
#line 272 "zstein.f"
    } else if (*n == 1) {
#line 273 "zstein.f"
	i__1 = z_dim1 + 1;
#line 273 "zstein.f"
	z__[i__1].r = 1., z__[i__1].i = 0.;
#line 274 "zstein.f"
	return 0;
#line 275 "zstein.f"
    }

/*     Get machine constants. */

#line 279 "zstein.f"
    eps = dlamch_("Precision", (ftnlen)9);

/*     Initialize seed for random number generator DLARNV. */

#line 283 "zstein.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 284 "zstein.f"
	iseed[i__ - 1] = 1;
#line 285 "zstein.f"
/* L40: */
#line 285 "zstein.f"
    }

/*     Initialize pointers. */

#line 289 "zstein.f"
    indrv1 = 0;
#line 290 "zstein.f"
    indrv2 = indrv1 + *n;
#line 291 "zstein.f"
    indrv3 = indrv2 + *n;
#line 292 "zstein.f"
    indrv4 = indrv3 + *n;
#line 293 "zstein.f"
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

#line 297 "zstein.f"
    j1 = 1;
#line 298 "zstein.f"
    i__1 = iblock[*m];
#line 298 "zstein.f"
    for (nblk = 1; nblk <= i__1; ++nblk) {

/*        Find starting and ending indices of block nblk. */

#line 302 "zstein.f"
	if (nblk == 1) {
#line 303 "zstein.f"
	    b1 = 1;
#line 304 "zstein.f"
	} else {
#line 305 "zstein.f"
	    b1 = isplit[nblk - 1] + 1;
#line 306 "zstein.f"
	}
#line 307 "zstein.f"
	bn = isplit[nblk];
#line 308 "zstein.f"
	blksiz = bn - b1 + 1;
#line 309 "zstein.f"
	if (blksiz == 1) {
#line 309 "zstein.f"
	    goto L60;
#line 309 "zstein.f"
	}
#line 311 "zstein.f"
	gpind = b1;

/*        Compute reorthogonalization criterion and stopping criterion. */

#line 315 "zstein.f"
	onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
/* Computing MAX */
#line 316 "zstein.f"
	d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1],
		 abs(d__2));
#line 316 "zstein.f"
	onenrm = max(d__3,d__4);
#line 317 "zstein.f"
	i__2 = bn - 1;
#line 317 "zstein.f"
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 318 "zstein.f"
	    d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
		    i__ - 1], abs(d__2)) + (d__3 = e[i__], abs(d__3));
#line 318 "zstein.f"
	    onenrm = max(d__4,d__5);
#line 320 "zstein.f"
/* L50: */
#line 320 "zstein.f"
	}
#line 321 "zstein.f"
	ortol = onenrm * .001;

#line 323 "zstein.f"
	dtpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

#line 327 "zstein.f"
L60:
#line 328 "zstein.f"
	jblk = 0;
#line 329 "zstein.f"
	i__2 = *m;
#line 329 "zstein.f"
	for (j = j1; j <= i__2; ++j) {
#line 330 "zstein.f"
	    if (iblock[j] != nblk) {
#line 331 "zstein.f"
		j1 = j;
#line 332 "zstein.f"
		goto L180;
#line 333 "zstein.f"
	    }
#line 334 "zstein.f"
	    ++jblk;
#line 335 "zstein.f"
	    xj = w[j];

/*           Skip all the work if the block size is one. */

#line 339 "zstein.f"
	    if (blksiz == 1) {
#line 340 "zstein.f"
		work[indrv1 + 1] = 1.;
#line 341 "zstein.f"
		goto L140;
#line 342 "zstein.f"
	    }

/*           If eigenvalues j and j-1 are too close, add a relatively */
/*           small perturbation. */

#line 347 "zstein.f"
	    if (jblk > 1) {
#line 348 "zstein.f"
		eps1 = (d__1 = eps * xj, abs(d__1));
#line 349 "zstein.f"
		pertol = eps1 * 10.;
#line 350 "zstein.f"
		sep = xj - xjm;
#line 351 "zstein.f"
		if (sep < pertol) {
#line 351 "zstein.f"
		    xj = xjm + pertol;
#line 351 "zstein.f"
		}
#line 353 "zstein.f"
	    }

#line 355 "zstein.f"
	    its = 0;
#line 356 "zstein.f"
	    nrmchk = 0;

/*           Get random starting vector. */

#line 360 "zstein.f"
	    dlarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

/*           Copy the matrix T so it won't be destroyed in factorization. */

#line 364 "zstein.f"
	    dcopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
#line 365 "zstein.f"
	    i__3 = blksiz - 1;
#line 365 "zstein.f"
	    dcopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
#line 366 "zstein.f"
	    i__3 = blksiz - 1;
#line 366 "zstein.f"
	    dcopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU ) */

#line 370 "zstein.f"
	    tol = 0.;
#line 371 "zstein.f"
	    dlagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

/*           Update iteration count. */

#line 377 "zstein.f"
L70:
#line 378 "zstein.f"
	    ++its;
#line 379 "zstein.f"
	    if (its > 5) {
#line 379 "zstein.f"
		goto L120;
#line 379 "zstein.f"
	    }

/*           Normalize and scale the righthand side vector Pb. */

/* Computing MAX */
#line 384 "zstein.f"
	    d__2 = eps, d__3 = (d__1 = work[indrv4 + blksiz], abs(d__1));
#line 384 "zstein.f"
	    scl = blksiz * onenrm * max(d__2,d__3) / dasum_(&blksiz, &work[
		    indrv1 + 1], &c__1);
#line 387 "zstein.f"
	    dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);

/*           Solve the system LU = Pb. */

#line 391 "zstein.f"
	    dlagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
/*           close enough. */

#line 398 "zstein.f"
	    if (jblk == 1) {
#line 398 "zstein.f"
		goto L110;
#line 398 "zstein.f"
	    }
#line 400 "zstein.f"
	    if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
#line 400 "zstein.f"
		gpind = j;
#line 400 "zstein.f"
	    }
#line 402 "zstein.f"
	    if (gpind != j) {
#line 403 "zstein.f"
		i__3 = j - 1;
#line 403 "zstein.f"
		for (i__ = gpind; i__ <= i__3; ++i__) {
#line 404 "zstein.f"
		    ztr = 0.;
#line 405 "zstein.f"
		    i__4 = blksiz;
#line 405 "zstein.f"
		    for (jr = 1; jr <= i__4; ++jr) {
#line 406 "zstein.f"
			i__5 = b1 - 1 + jr + i__ * z_dim1;
#line 406 "zstein.f"
			ztr += work[indrv1 + jr] * z__[i__5].r;
#line 408 "zstein.f"
/* L80: */
#line 408 "zstein.f"
		    }
#line 409 "zstein.f"
		    i__4 = blksiz;
#line 409 "zstein.f"
		    for (jr = 1; jr <= i__4; ++jr) {
#line 410 "zstein.f"
			i__5 = b1 - 1 + jr + i__ * z_dim1;
#line 410 "zstein.f"
			work[indrv1 + jr] -= ztr * z__[i__5].r;
#line 412 "zstein.f"
/* L90: */
#line 412 "zstein.f"
		    }
#line 413 "zstein.f"
/* L100: */
#line 413 "zstein.f"
		}
#line 414 "zstein.f"
	    }

/*           Check the infinity norm of the iterate. */

#line 418 "zstein.f"
L110:
#line 419 "zstein.f"
	    jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 420 "zstein.f"
	    nrm = (d__1 = work[indrv1 + jmax], abs(d__1));

/*           Continue for additional iterations after norm reaches */
/*           stopping criterion. */

#line 425 "zstein.f"
	    if (nrm < dtpcrt) {
#line 425 "zstein.f"
		goto L70;
#line 425 "zstein.f"
	    }
#line 427 "zstein.f"
	    ++nrmchk;
#line 428 "zstein.f"
	    if (nrmchk < 3) {
#line 428 "zstein.f"
		goto L70;
#line 428 "zstein.f"
	    }

#line 431 "zstein.f"
	    goto L130;

/*           If stopping criterion was not satisfied, update info and */
/*           store eigenvector number in array ifail. */

#line 436 "zstein.f"
L120:
#line 437 "zstein.f"
	    ++(*info);
#line 438 "zstein.f"
	    ifail[*info] = j;

/*           Accept iterate as jth eigenvector. */

#line 442 "zstein.f"
L130:
#line 443 "zstein.f"
	    scl = 1. / dnrm2_(&blksiz, &work[indrv1 + 1], &c__1);
#line 444 "zstein.f"
	    jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 445 "zstein.f"
	    if (work[indrv1 + jmax] < 0.) {
#line 445 "zstein.f"
		scl = -scl;
#line 445 "zstein.f"
	    }
#line 447 "zstein.f"
	    dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
#line 448 "zstein.f"
L140:
#line 449 "zstein.f"
	    i__3 = *n;
#line 449 "zstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 450 "zstein.f"
		i__4 = i__ + j * z_dim1;
#line 450 "zstein.f"
		z__[i__4].r = 0., z__[i__4].i = 0.;
#line 451 "zstein.f"
/* L150: */
#line 451 "zstein.f"
	    }
#line 452 "zstein.f"
	    i__3 = blksiz;
#line 452 "zstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 453 "zstein.f"
		i__4 = b1 + i__ - 1 + j * z_dim1;
#line 453 "zstein.f"
		i__5 = indrv1 + i__;
#line 453 "zstein.f"
		z__1.r = work[i__5], z__1.i = 0.;
#line 453 "zstein.f"
		z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 454 "zstein.f"
/* L160: */
#line 454 "zstein.f"
	    }

/*           Save the shift to check eigenvalue spacing at next */
/*           iteration. */

#line 459 "zstein.f"
	    xjm = xj;

#line 461 "zstein.f"
/* L170: */
#line 461 "zstein.f"
	}
#line 462 "zstein.f"
L180:
#line 462 "zstein.f"
	;
#line 462 "zstein.f"
    }

#line 464 "zstein.f"
    return 0;

/*     End of ZSTEIN */

} /* zstein_ */

