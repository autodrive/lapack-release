#line 1 "cstein.f"
/* cstein.f -- translated by f2c (version 20100827).
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

#line 1 "cstein.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CSTEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSTEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cstein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cstein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cstein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDZ, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/*      $                   IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ) */
/*       COMPLEX            Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSTEIN computes the eigenvectors of a real symmetric tridiagonal */
/* > matrix T corresponding to specified eigenvalues, using inverse */
/* > iteration. */
/* > */
/* > The maximum number of iterations allowed for each eigenvector is */
/* > specified by an internal parameter MAXITS (currently set to 5). */
/* > */
/* > Although the eigenvectors are real, they are stored in a complex */
/* > array, which may be passed to CUNMTR or CUPMTR for back */
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
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
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
/* >          W is REAL array, dimension (N) */
/* >          The first M elements of W contain the eigenvalues for */
/* >          which eigenvectors are to be computed.  The eigenvalues */
/* >          should be grouped by split-off block and ordered from */
/* >          smallest to largest within the block.  ( The output array */
/* >          W from SSTEBZ with ORDER = 'B' is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] IBLOCK */
/* > \verbatim */
/* >          IBLOCK is INTEGER array, dimension (N) */
/* >          The submatrix indices associated with the corresponding */
/* >          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to */
/* >          the first submatrix from the top, =2 if W(i) belongs to */
/* >          the second submatrix, etc.  ( The output array IBLOCK */
/* >          from SSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] ISPLIT */
/* > \verbatim */
/* >          ISPLIT is INTEGER array, dimension (N) */
/* >          The splitting points, at which T breaks up into submatrices. */
/* >          The first submatrix consists of rows/columns 1 to */
/* >          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1 */
/* >          through ISPLIT( 2 ), etc. */
/* >          ( The output array ISPLIT from SSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, M) */
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
/* >          WORK is REAL array, dimension (5*N) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cstein_(integer *n, doublereal *d__, doublereal *e, 
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
    static doublereal xj, scl, eps, ctr, sep, nrm, tol;
    static integer its;
    static doublereal xjm, eps1;
    static integer jblk, nblk, jmax;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static integer iseed[4], gpind, iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), scopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static doublereal ortol;
    static integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slagtf_(
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *);
    static integer nrmchk;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slagts_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static integer blksiz;
    static doublereal onenrm, pertol;
    extern /* Subroutine */ int slarnv_(integer *, integer *, integer *, 
	    doublereal *);
    static doublereal stpcrt;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

/*     Test the input parameters. */

#line 237 "cstein.f"
    /* Parameter adjustments */
#line 237 "cstein.f"
    --d__;
#line 237 "cstein.f"
    --e;
#line 237 "cstein.f"
    --w;
#line 237 "cstein.f"
    --iblock;
#line 237 "cstein.f"
    --isplit;
#line 237 "cstein.f"
    z_dim1 = *ldz;
#line 237 "cstein.f"
    z_offset = 1 + z_dim1;
#line 237 "cstein.f"
    z__ -= z_offset;
#line 237 "cstein.f"
    --work;
#line 237 "cstein.f"
    --iwork;
#line 237 "cstein.f"
    --ifail;
#line 237 "cstein.f"

#line 237 "cstein.f"
    /* Function Body */
#line 237 "cstein.f"
    *info = 0;
#line 238 "cstein.f"
    i__1 = *m;
#line 238 "cstein.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 239 "cstein.f"
	ifail[i__] = 0;
#line 240 "cstein.f"
/* L10: */
#line 240 "cstein.f"
    }

#line 242 "cstein.f"
    if (*n < 0) {
#line 243 "cstein.f"
	*info = -1;
#line 244 "cstein.f"
    } else if (*m < 0 || *m > *n) {
#line 245 "cstein.f"
	*info = -4;
#line 246 "cstein.f"
    } else if (*ldz < max(1,*n)) {
#line 247 "cstein.f"
	*info = -9;
#line 248 "cstein.f"
    } else {
#line 249 "cstein.f"
	i__1 = *m;
#line 249 "cstein.f"
	for (j = 2; j <= i__1; ++j) {
#line 250 "cstein.f"
	    if (iblock[j] < iblock[j - 1]) {
#line 251 "cstein.f"
		*info = -6;
#line 252 "cstein.f"
		goto L30;
#line 253 "cstein.f"
	    }
#line 254 "cstein.f"
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
#line 256 "cstein.f"
		*info = -5;
#line 257 "cstein.f"
		goto L30;
#line 258 "cstein.f"
	    }
#line 259 "cstein.f"
/* L20: */
#line 259 "cstein.f"
	}
#line 260 "cstein.f"
L30:
#line 261 "cstein.f"
	;
#line 261 "cstein.f"
    }

#line 263 "cstein.f"
    if (*info != 0) {
#line 264 "cstein.f"
	i__1 = -(*info);
#line 264 "cstein.f"
	xerbla_("CSTEIN", &i__1, (ftnlen)6);
#line 265 "cstein.f"
	return 0;
#line 266 "cstein.f"
    }

/*     Quick return if possible */

#line 270 "cstein.f"
    if (*n == 0 || *m == 0) {
#line 271 "cstein.f"
	return 0;
#line 272 "cstein.f"
    } else if (*n == 1) {
#line 273 "cstein.f"
	i__1 = z_dim1 + 1;
#line 273 "cstein.f"
	z__[i__1].r = 1., z__[i__1].i = 0.;
#line 274 "cstein.f"
	return 0;
#line 275 "cstein.f"
    }

/*     Get machine constants. */

#line 279 "cstein.f"
    eps = slamch_("Precision", (ftnlen)9);

/*     Initialize seed for random number generator SLARNV. */

#line 283 "cstein.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 284 "cstein.f"
	iseed[i__ - 1] = 1;
#line 285 "cstein.f"
/* L40: */
#line 285 "cstein.f"
    }

/*     Initialize pointers. */

#line 289 "cstein.f"
    indrv1 = 0;
#line 290 "cstein.f"
    indrv2 = indrv1 + *n;
#line 291 "cstein.f"
    indrv3 = indrv2 + *n;
#line 292 "cstein.f"
    indrv4 = indrv3 + *n;
#line 293 "cstein.f"
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

#line 297 "cstein.f"
    j1 = 1;
#line 298 "cstein.f"
    i__1 = iblock[*m];
#line 298 "cstein.f"
    for (nblk = 1; nblk <= i__1; ++nblk) {

/*        Find starting and ending indices of block nblk. */

#line 302 "cstein.f"
	if (nblk == 1) {
#line 303 "cstein.f"
	    b1 = 1;
#line 304 "cstein.f"
	} else {
#line 305 "cstein.f"
	    b1 = isplit[nblk - 1] + 1;
#line 306 "cstein.f"
	}
#line 307 "cstein.f"
	bn = isplit[nblk];
#line 308 "cstein.f"
	blksiz = bn - b1 + 1;
#line 309 "cstein.f"
	if (blksiz == 1) {
#line 309 "cstein.f"
	    goto L60;
#line 309 "cstein.f"
	}
#line 311 "cstein.f"
	gpind = j1;

/*        Compute reorthogonalization criterion and stopping criterion. */

#line 315 "cstein.f"
	onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
/* Computing MAX */
#line 316 "cstein.f"
	d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1],
		 abs(d__2));
#line 316 "cstein.f"
	onenrm = max(d__3,d__4);
#line 317 "cstein.f"
	i__2 = bn - 1;
#line 317 "cstein.f"
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 318 "cstein.f"
	    d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
		    i__ - 1], abs(d__2)) + (d__3 = e[i__], abs(d__3));
#line 318 "cstein.f"
	    onenrm = max(d__4,d__5);
#line 320 "cstein.f"
/* L50: */
#line 320 "cstein.f"
	}
#line 321 "cstein.f"
	ortol = onenrm * .001;

#line 323 "cstein.f"
	stpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

#line 327 "cstein.f"
L60:
#line 328 "cstein.f"
	jblk = 0;
#line 329 "cstein.f"
	i__2 = *m;
#line 329 "cstein.f"
	for (j = j1; j <= i__2; ++j) {
#line 330 "cstein.f"
	    if (iblock[j] != nblk) {
#line 331 "cstein.f"
		j1 = j;
#line 332 "cstein.f"
		goto L180;
#line 333 "cstein.f"
	    }
#line 334 "cstein.f"
	    ++jblk;
#line 335 "cstein.f"
	    xj = w[j];

/*           Skip all the work if the block size is one. */

#line 339 "cstein.f"
	    if (blksiz == 1) {
#line 340 "cstein.f"
		work[indrv1 + 1] = 1.;
#line 341 "cstein.f"
		goto L140;
#line 342 "cstein.f"
	    }

/*           If eigenvalues j and j-1 are too close, add a relatively */
/*           small perturbation. */

#line 347 "cstein.f"
	    if (jblk > 1) {
#line 348 "cstein.f"
		eps1 = (d__1 = eps * xj, abs(d__1));
#line 349 "cstein.f"
		pertol = eps1 * 10.;
#line 350 "cstein.f"
		sep = xj - xjm;
#line 351 "cstein.f"
		if (sep < pertol) {
#line 351 "cstein.f"
		    xj = xjm + pertol;
#line 351 "cstein.f"
		}
#line 353 "cstein.f"
	    }

#line 355 "cstein.f"
	    its = 0;
#line 356 "cstein.f"
	    nrmchk = 0;

/*           Get random starting vector. */

#line 360 "cstein.f"
	    slarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

/*           Copy the matrix T so it won't be destroyed in factorization. */

#line 364 "cstein.f"
	    scopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
#line 365 "cstein.f"
	    i__3 = blksiz - 1;
#line 365 "cstein.f"
	    scopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
#line 366 "cstein.f"
	    i__3 = blksiz - 1;
#line 366 "cstein.f"
	    scopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU ) */

#line 370 "cstein.f"
	    tol = 0.;
#line 371 "cstein.f"
	    slagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

/*           Update iteration count. */

#line 377 "cstein.f"
L70:
#line 378 "cstein.f"
	    ++its;
#line 379 "cstein.f"
	    if (its > 5) {
#line 379 "cstein.f"
		goto L120;
#line 379 "cstein.f"
	    }

/*           Normalize and scale the righthand side vector Pb. */

#line 384 "cstein.f"
	    jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
/* Computing MAX */
#line 385 "cstein.f"
	    d__3 = eps, d__4 = (d__1 = work[indrv4 + blksiz], abs(d__1));
#line 385 "cstein.f"
	    scl = blksiz * onenrm * max(d__3,d__4) / (d__2 = work[indrv1 + 
		    jmax], abs(d__2));
#line 388 "cstein.f"
	    sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);

/*           Solve the system LU = Pb. */

#line 392 "cstein.f"
	    slagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
/*           close enough. */

#line 399 "cstein.f"
	    if (jblk == 1) {
#line 399 "cstein.f"
		goto L110;
#line 399 "cstein.f"
	    }
#line 401 "cstein.f"
	    if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
#line 401 "cstein.f"
		gpind = j;
#line 401 "cstein.f"
	    }
#line 403 "cstein.f"
	    if (gpind != j) {
#line 404 "cstein.f"
		i__3 = j - 1;
#line 404 "cstein.f"
		for (i__ = gpind; i__ <= i__3; ++i__) {
#line 405 "cstein.f"
		    ctr = 0.;
#line 406 "cstein.f"
		    i__4 = blksiz;
#line 406 "cstein.f"
		    for (jr = 1; jr <= i__4; ++jr) {
#line 407 "cstein.f"
			i__5 = b1 - 1 + jr + i__ * z_dim1;
#line 407 "cstein.f"
			ctr += work[indrv1 + jr] * z__[i__5].r;
#line 409 "cstein.f"
/* L80: */
#line 409 "cstein.f"
		    }
#line 410 "cstein.f"
		    i__4 = blksiz;
#line 410 "cstein.f"
		    for (jr = 1; jr <= i__4; ++jr) {
#line 411 "cstein.f"
			i__5 = b1 - 1 + jr + i__ * z_dim1;
#line 411 "cstein.f"
			work[indrv1 + jr] -= ctr * z__[i__5].r;
#line 413 "cstein.f"
/* L90: */
#line 413 "cstein.f"
		    }
#line 414 "cstein.f"
/* L100: */
#line 414 "cstein.f"
		}
#line 415 "cstein.f"
	    }

/*           Check the infinity norm of the iterate. */

#line 419 "cstein.f"
L110:
#line 420 "cstein.f"
	    jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 421 "cstein.f"
	    nrm = (d__1 = work[indrv1 + jmax], abs(d__1));

/*           Continue for additional iterations after norm reaches */
/*           stopping criterion. */

#line 426 "cstein.f"
	    if (nrm < stpcrt) {
#line 426 "cstein.f"
		goto L70;
#line 426 "cstein.f"
	    }
#line 428 "cstein.f"
	    ++nrmchk;
#line 429 "cstein.f"
	    if (nrmchk < 3) {
#line 429 "cstein.f"
		goto L70;
#line 429 "cstein.f"
	    }

#line 432 "cstein.f"
	    goto L130;

/*           If stopping criterion was not satisfied, update info and */
/*           store eigenvector number in array ifail. */

#line 437 "cstein.f"
L120:
#line 438 "cstein.f"
	    ++(*info);
#line 439 "cstein.f"
	    ifail[*info] = j;

/*           Accept iterate as jth eigenvector. */

#line 443 "cstein.f"
L130:
#line 444 "cstein.f"
	    scl = 1. / snrm2_(&blksiz, &work[indrv1 + 1], &c__1);
#line 445 "cstein.f"
	    jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 446 "cstein.f"
	    if (work[indrv1 + jmax] < 0.) {
#line 446 "cstein.f"
		scl = -scl;
#line 446 "cstein.f"
	    }
#line 448 "cstein.f"
	    sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
#line 449 "cstein.f"
L140:
#line 450 "cstein.f"
	    i__3 = *n;
#line 450 "cstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 451 "cstein.f"
		i__4 = i__ + j * z_dim1;
#line 451 "cstein.f"
		z__[i__4].r = 0., z__[i__4].i = 0.;
#line 452 "cstein.f"
/* L150: */
#line 452 "cstein.f"
	    }
#line 453 "cstein.f"
	    i__3 = blksiz;
#line 453 "cstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 454 "cstein.f"
		i__4 = b1 + i__ - 1 + j * z_dim1;
#line 454 "cstein.f"
		i__5 = indrv1 + i__;
#line 454 "cstein.f"
		z__1.r = work[i__5], z__1.i = 0.;
#line 454 "cstein.f"
		z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 455 "cstein.f"
/* L160: */
#line 455 "cstein.f"
	    }

/*           Save the shift to check eigenvalue spacing at next */
/*           iteration. */

#line 460 "cstein.f"
	    xjm = xj;

#line 462 "cstein.f"
/* L170: */
#line 462 "cstein.f"
	}
#line 463 "cstein.f"
L180:
#line 463 "cstein.f"
	;
#line 463 "cstein.f"
    }

#line 465 "cstein.f"
    return 0;

/*     End of CSTEIN */

} /* cstein_ */

