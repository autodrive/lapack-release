#line 1 "sstein.f"
/* sstein.f -- translated by f2c (version 20100827).
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

#line 1 "sstein.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b SSTEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/*                          IWORK, IFAIL, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDZ, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/*      $                   IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEIN computes the eigenvectors of a real symmetric tridiagonal */
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
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
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
/* >          Z is REAL array, dimension (LDZ, M) */
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

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sstein_(integer *n, doublereal *d__, doublereal *e, 
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
    static doublereal xj, scl, eps, ctr, sep, nrm, tol;
    static integer its;
    static doublereal xjm, eps1;
    static integer jblk, nblk, jmax;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), snrm2_(integer *, doublereal *, integer *);
    static integer iseed[4], gpind, iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), scopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static doublereal ortol;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
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

#line 226 "sstein.f"
    /* Parameter adjustments */
#line 226 "sstein.f"
    --d__;
#line 226 "sstein.f"
    --e;
#line 226 "sstein.f"
    --w;
#line 226 "sstein.f"
    --iblock;
#line 226 "sstein.f"
    --isplit;
#line 226 "sstein.f"
    z_dim1 = *ldz;
#line 226 "sstein.f"
    z_offset = 1 + z_dim1;
#line 226 "sstein.f"
    z__ -= z_offset;
#line 226 "sstein.f"
    --work;
#line 226 "sstein.f"
    --iwork;
#line 226 "sstein.f"
    --ifail;
#line 226 "sstein.f"

#line 226 "sstein.f"
    /* Function Body */
#line 226 "sstein.f"
    *info = 0;
#line 227 "sstein.f"
    i__1 = *m;
#line 227 "sstein.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 228 "sstein.f"
	ifail[i__] = 0;
#line 229 "sstein.f"
/* L10: */
#line 229 "sstein.f"
    }

#line 231 "sstein.f"
    if (*n < 0) {
#line 232 "sstein.f"
	*info = -1;
#line 233 "sstein.f"
    } else if (*m < 0 || *m > *n) {
#line 234 "sstein.f"
	*info = -4;
#line 235 "sstein.f"
    } else if (*ldz < max(1,*n)) {
#line 236 "sstein.f"
	*info = -9;
#line 237 "sstein.f"
    } else {
#line 238 "sstein.f"
	i__1 = *m;
#line 238 "sstein.f"
	for (j = 2; j <= i__1; ++j) {
#line 239 "sstein.f"
	    if (iblock[j] < iblock[j - 1]) {
#line 240 "sstein.f"
		*info = -6;
#line 241 "sstein.f"
		goto L30;
#line 242 "sstein.f"
	    }
#line 243 "sstein.f"
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
#line 245 "sstein.f"
		*info = -5;
#line 246 "sstein.f"
		goto L30;
#line 247 "sstein.f"
	    }
#line 248 "sstein.f"
/* L20: */
#line 248 "sstein.f"
	}
#line 249 "sstein.f"
L30:
#line 250 "sstein.f"
	;
#line 250 "sstein.f"
    }

#line 252 "sstein.f"
    if (*info != 0) {
#line 253 "sstein.f"
	i__1 = -(*info);
#line 253 "sstein.f"
	xerbla_("SSTEIN", &i__1, (ftnlen)6);
#line 254 "sstein.f"
	return 0;
#line 255 "sstein.f"
    }

/*     Quick return if possible */

#line 259 "sstein.f"
    if (*n == 0 || *m == 0) {
#line 260 "sstein.f"
	return 0;
#line 261 "sstein.f"
    } else if (*n == 1) {
#line 262 "sstein.f"
	z__[z_dim1 + 1] = 1.;
#line 263 "sstein.f"
	return 0;
#line 264 "sstein.f"
    }

/*     Get machine constants. */

#line 268 "sstein.f"
    eps = slamch_("Precision", (ftnlen)9);

/*     Initialize seed for random number generator SLARNV. */

#line 272 "sstein.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 273 "sstein.f"
	iseed[i__ - 1] = 1;
#line 274 "sstein.f"
/* L40: */
#line 274 "sstein.f"
    }

/*     Initialize pointers. */

#line 278 "sstein.f"
    indrv1 = 0;
#line 279 "sstein.f"
    indrv2 = indrv1 + *n;
#line 280 "sstein.f"
    indrv3 = indrv2 + *n;
#line 281 "sstein.f"
    indrv4 = indrv3 + *n;
#line 282 "sstein.f"
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

#line 286 "sstein.f"
    j1 = 1;
#line 287 "sstein.f"
    i__1 = iblock[*m];
#line 287 "sstein.f"
    for (nblk = 1; nblk <= i__1; ++nblk) {

/*        Find starting and ending indices of block nblk. */

#line 291 "sstein.f"
	if (nblk == 1) {
#line 292 "sstein.f"
	    b1 = 1;
#line 293 "sstein.f"
	} else {
#line 294 "sstein.f"
	    b1 = isplit[nblk - 1] + 1;
#line 295 "sstein.f"
	}
#line 296 "sstein.f"
	bn = isplit[nblk];
#line 297 "sstein.f"
	blksiz = bn - b1 + 1;
#line 298 "sstein.f"
	if (blksiz == 1) {
#line 298 "sstein.f"
	    goto L60;
#line 298 "sstein.f"
	}
#line 300 "sstein.f"
	gpind = j1;

/*        Compute reorthogonalization criterion and stopping criterion. */

#line 304 "sstein.f"
	onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
/* Computing MAX */
#line 305 "sstein.f"
	d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1],
		 abs(d__2));
#line 305 "sstein.f"
	onenrm = max(d__3,d__4);
#line 306 "sstein.f"
	i__2 = bn - 1;
#line 306 "sstein.f"
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 307 "sstein.f"
	    d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
		    i__ - 1], abs(d__2)) + (d__3 = e[i__], abs(d__3));
#line 307 "sstein.f"
	    onenrm = max(d__4,d__5);
#line 309 "sstein.f"
/* L50: */
#line 309 "sstein.f"
	}
#line 310 "sstein.f"
	ortol = onenrm * .001;

#line 312 "sstein.f"
	stpcrt = sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

#line 316 "sstein.f"
L60:
#line 317 "sstein.f"
	jblk = 0;
#line 318 "sstein.f"
	i__2 = *m;
#line 318 "sstein.f"
	for (j = j1; j <= i__2; ++j) {
#line 319 "sstein.f"
	    if (iblock[j] != nblk) {
#line 320 "sstein.f"
		j1 = j;
#line 321 "sstein.f"
		goto L160;
#line 322 "sstein.f"
	    }
#line 323 "sstein.f"
	    ++jblk;
#line 324 "sstein.f"
	    xj = w[j];

/*           Skip all the work if the block size is one. */

#line 328 "sstein.f"
	    if (blksiz == 1) {
#line 329 "sstein.f"
		work[indrv1 + 1] = 1.;
#line 330 "sstein.f"
		goto L120;
#line 331 "sstein.f"
	    }

/*           If eigenvalues j and j-1 are too close, add a relatively */
/*           small perturbation. */

#line 336 "sstein.f"
	    if (jblk > 1) {
#line 337 "sstein.f"
		eps1 = (d__1 = eps * xj, abs(d__1));
#line 338 "sstein.f"
		pertol = eps1 * 10.;
#line 339 "sstein.f"
		sep = xj - xjm;
#line 340 "sstein.f"
		if (sep < pertol) {
#line 340 "sstein.f"
		    xj = xjm + pertol;
#line 340 "sstein.f"
		}
#line 342 "sstein.f"
	    }

#line 344 "sstein.f"
	    its = 0;
#line 345 "sstein.f"
	    nrmchk = 0;

/*           Get random starting vector. */

#line 349 "sstein.f"
	    slarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

/*           Copy the matrix T so it won't be destroyed in factorization. */

#line 353 "sstein.f"
	    scopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
#line 354 "sstein.f"
	    i__3 = blksiz - 1;
#line 354 "sstein.f"
	    scopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
#line 355 "sstein.f"
	    i__3 = blksiz - 1;
#line 355 "sstein.f"
	    scopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU ) */

#line 359 "sstein.f"
	    tol = 0.;
#line 360 "sstein.f"
	    slagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

/*           Update iteration count. */

#line 366 "sstein.f"
L70:
#line 367 "sstein.f"
	    ++its;
#line 368 "sstein.f"
	    if (its > 5) {
#line 368 "sstein.f"
		goto L100;
#line 368 "sstein.f"
	    }

/*           Normalize and scale the righthand side vector Pb. */

#line 373 "sstein.f"
	    jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
/* Computing MAX */
#line 374 "sstein.f"
	    d__3 = eps, d__4 = (d__1 = work[indrv4 + blksiz], abs(d__1));
#line 374 "sstein.f"
	    scl = blksiz * onenrm * max(d__3,d__4) / (d__2 = work[indrv1 + 
		    jmax], abs(d__2));
#line 377 "sstein.f"
	    sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);

/*           Solve the system LU = Pb. */

#line 381 "sstein.f"
	    slagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
/*           close enough. */

#line 388 "sstein.f"
	    if (jblk == 1) {
#line 388 "sstein.f"
		goto L90;
#line 388 "sstein.f"
	    }
#line 390 "sstein.f"
	    if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
#line 390 "sstein.f"
		gpind = j;
#line 390 "sstein.f"
	    }
#line 392 "sstein.f"
	    if (gpind != j) {
#line 393 "sstein.f"
		i__3 = j - 1;
#line 393 "sstein.f"
		for (i__ = gpind; i__ <= i__3; ++i__) {
#line 394 "sstein.f"
		    ctr = -sdot_(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + 
			    i__ * z_dim1], &c__1);
#line 396 "sstein.f"
		    saxpy_(&blksiz, &ctr, &z__[b1 + i__ * z_dim1], &c__1, &
			    work[indrv1 + 1], &c__1);
#line 398 "sstein.f"
/* L80: */
#line 398 "sstein.f"
		}
#line 399 "sstein.f"
	    }

/*           Check the infinity norm of the iterate. */

#line 403 "sstein.f"
L90:
#line 404 "sstein.f"
	    jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 405 "sstein.f"
	    nrm = (d__1 = work[indrv1 + jmax], abs(d__1));

/*           Continue for additional iterations after norm reaches */
/*           stopping criterion. */

#line 410 "sstein.f"
	    if (nrm < stpcrt) {
#line 410 "sstein.f"
		goto L70;
#line 410 "sstein.f"
	    }
#line 412 "sstein.f"
	    ++nrmchk;
#line 413 "sstein.f"
	    if (nrmchk < 3) {
#line 413 "sstein.f"
		goto L70;
#line 413 "sstein.f"
	    }

#line 416 "sstein.f"
	    goto L110;

/*           If stopping criterion was not satisfied, update info and */
/*           store eigenvector number in array ifail. */

#line 421 "sstein.f"
L100:
#line 422 "sstein.f"
	    ++(*info);
#line 423 "sstein.f"
	    ifail[*info] = j;

/*           Accept iterate as jth eigenvector. */

#line 427 "sstein.f"
L110:
#line 428 "sstein.f"
	    scl = 1. / snrm2_(&blksiz, &work[indrv1 + 1], &c__1);
#line 429 "sstein.f"
	    jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
#line 430 "sstein.f"
	    if (work[indrv1 + jmax] < 0.) {
#line 430 "sstein.f"
		scl = -scl;
#line 430 "sstein.f"
	    }
#line 432 "sstein.f"
	    sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
#line 433 "sstein.f"
L120:
#line 434 "sstein.f"
	    i__3 = *n;
#line 434 "sstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 435 "sstein.f"
		z__[i__ + j * z_dim1] = 0.;
#line 436 "sstein.f"
/* L130: */
#line 436 "sstein.f"
	    }
#line 437 "sstein.f"
	    i__3 = blksiz;
#line 437 "sstein.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 438 "sstein.f"
		z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
#line 439 "sstein.f"
/* L140: */
#line 439 "sstein.f"
	    }

/*           Save the shift to check eigenvalue spacing at next */
/*           iteration. */

#line 444 "sstein.f"
	    xjm = xj;

#line 446 "sstein.f"
/* L150: */
#line 446 "sstein.f"
	}
#line 447 "sstein.f"
L160:
#line 447 "sstein.f"
	;
#line 447 "sstein.f"
    }

#line 449 "sstein.f"
    return 0;

/*     End of SSTEIN */

} /* sstein_ */

