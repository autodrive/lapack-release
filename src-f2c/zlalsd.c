#line 1 "zlalsd.f"
/* zlalsd.f -- translated by f2c (version 20100827).
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

#line 1 "zlalsd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b10 = 1.;
static doublereal c_b35 = 0.;

/* > \brief \b ZLALSD uses the singular value decomposition of A to solve the least squares problem. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLALSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlalsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlalsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlalsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, */
/*                          RANK, WORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * ) */
/*       COMPLEX*16         B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLALSD uses the singular value decomposition of A to solve the least */
/* > squares problem of finding X to minimize the Euclidean norm of each */
/* > column of A*X-B, where A is N-by-N upper bidiagonal, and X and B */
/* > are N-by-NRHS. The solution X overwrites B. */
/* > */
/* > The singular values of A smaller than RCOND times the largest */
/* > singular value are treated as zero in solving the least squares */
/* > problem; in this case a minimum norm solution is returned. */
/* > The actual singular values are returned in D in ascending order. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >         = 'U': D and E define an upper bidiagonal matrix. */
/* >         = 'L': D and E define a  lower bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLSIZ */
/* > \verbatim */
/* >          SMLSIZ is INTEGER */
/* >         The maximum size of the subproblems at the bottom of the */
/* >         computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The dimension of the  bidiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >         The number of columns of B. NRHS must be at least 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry D contains the main diagonal of the bidiagonal */
/* >         matrix. On exit, if INFO = 0, D contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >         Contains the super-diagonal entries of the bidiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >         On input, B contains the right hand sides of the least */
/* >         squares problem. On output, B contains the solution X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >         The leading dimension of B in the calling subprogram. */
/* >         LDB must be at least max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >         The singular values of A less than or equal to RCOND times */
/* >         the largest singular value are treated as zero in solving */
/* >         the least squares problem. If RCOND is negative, */
/* >         machine precision is used instead. */
/* >         For example, if diag(S)*X=B were the least squares problem, */
/* >         where diag(S) is a diagonal matrix of singular values, the */
/* >         solution would be X(i) = B(i) / S(i) if S(i) is greater than */
/* >         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to */
/* >         RCOND*max(S). */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* >          RANK is INTEGER */
/* >         The number of singular values of A greater than RCOND times */
/* >         the largest singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N * NRHS) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension at least */
/* >         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + */
/* >         MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ), */
/* >         where */
/* >         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension at least */
/* >         (3*N*NLVL + 11*N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >         = 0:  successful exit. */
/* >         < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >         > 0:  The algorithm failed to compute a singular value while */
/* >               working on the submatrix lying in rows and columns */
/* >               INFO/(N+1) through MOD(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int zlalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, doublereal *d__, doublereal *e, doublecomplex *b, integer *ldb,
	 doublereal *rcond, integer *rank, doublecomplex *work, doublereal *
	rwork, integer *iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), log(doublereal), d_sign(doublereal *, 
	    doublereal *);

    /* Local variables */
    static integer c__, i__, j, k;
    static doublereal r__;
    static integer s, u, z__;
    static doublereal cs;
    static integer bx;
    static doublereal sn;
    static integer st, vt, nm1, st1;
    static doublereal eps;
    static integer iwk;
    static doublereal tol;
    static integer difl, difr;
    static doublereal rcnd;
    static integer jcol, irwb, perm, nsub, nlvl, sqre, bxst, jrow, irwu, 
	    jimag;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer jreal, irwib, poles, sizei, irwrb, nsize;
    extern /* Subroutine */ int zdrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *), zcopy_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    ;
    static integer irwvt, icmpq1, icmpq2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlasda_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlasdq_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), dlartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), xerbla_(char *, integer *, ftnlen);
    static integer givcol;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int zlalsa_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), zlascl_(char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen), dlasrt_(char *, 
	    integer *, doublereal *, integer *, ftnlen), zlacpy_(char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, ftnlen), zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static doublereal orgnrm;
    static integer givnum, givptr, nrwork, irwwrk, smlszp;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

/*     Test the input parameters. */

#line 240 "zlalsd.f"
    /* Parameter adjustments */
#line 240 "zlalsd.f"
    --d__;
#line 240 "zlalsd.f"
    --e;
#line 240 "zlalsd.f"
    b_dim1 = *ldb;
#line 240 "zlalsd.f"
    b_offset = 1 + b_dim1;
#line 240 "zlalsd.f"
    b -= b_offset;
#line 240 "zlalsd.f"
    --work;
#line 240 "zlalsd.f"
    --rwork;
#line 240 "zlalsd.f"
    --iwork;
#line 240 "zlalsd.f"

#line 240 "zlalsd.f"
    /* Function Body */
#line 240 "zlalsd.f"
    *info = 0;

#line 242 "zlalsd.f"
    if (*n < 0) {
#line 243 "zlalsd.f"
	*info = -3;
#line 244 "zlalsd.f"
    } else if (*nrhs < 1) {
#line 245 "zlalsd.f"
	*info = -4;
#line 246 "zlalsd.f"
    } else if (*ldb < 1 || *ldb < *n) {
#line 247 "zlalsd.f"
	*info = -8;
#line 248 "zlalsd.f"
    }
#line 249 "zlalsd.f"
    if (*info != 0) {
#line 250 "zlalsd.f"
	i__1 = -(*info);
#line 250 "zlalsd.f"
	xerbla_("ZLALSD", &i__1, (ftnlen)6);
#line 251 "zlalsd.f"
	return 0;
#line 252 "zlalsd.f"
    }

#line 254 "zlalsd.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Set up the tolerance. */

#line 258 "zlalsd.f"
    if (*rcond <= 0. || *rcond >= 1.) {
#line 259 "zlalsd.f"
	rcnd = eps;
#line 260 "zlalsd.f"
    } else {
#line 261 "zlalsd.f"
	rcnd = *rcond;
#line 262 "zlalsd.f"
    }

#line 264 "zlalsd.f"
    *rank = 0;

/*     Quick return if possible. */

#line 268 "zlalsd.f"
    if (*n == 0) {
#line 269 "zlalsd.f"
	return 0;
#line 270 "zlalsd.f"
    } else if (*n == 1) {
#line 271 "zlalsd.f"
	if (d__[1] == 0.) {
#line 272 "zlalsd.f"
	    zlaset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (
		    ftnlen)1);
#line 273 "zlalsd.f"
	} else {
#line 274 "zlalsd.f"
	    *rank = 1;
#line 275 "zlalsd.f"
	    zlascl_("G", &c__0, &c__0, &d__[1], &c_b10, &c__1, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)1);
#line 276 "zlalsd.f"
	    d__[1] = abs(d__[1]);
#line 277 "zlalsd.f"
	}
#line 278 "zlalsd.f"
	return 0;
#line 279 "zlalsd.f"
    }

/*     Rotate the matrix if it is lower bidiagonal. */

#line 283 "zlalsd.f"
    if (*(unsigned char *)uplo == 'L') {
#line 284 "zlalsd.f"
	i__1 = *n - 1;
#line 284 "zlalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "zlalsd.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 286 "zlalsd.f"
	    d__[i__] = r__;
#line 287 "zlalsd.f"
	    e[i__] = sn * d__[i__ + 1];
#line 288 "zlalsd.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 289 "zlalsd.f"
	    if (*nrhs == 1) {
#line 290 "zlalsd.f"
		zdrot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &
			c__1, &cs, &sn);
#line 291 "zlalsd.f"
	    } else {
#line 292 "zlalsd.f"
		rwork[(i__ << 1) - 1] = cs;
#line 293 "zlalsd.f"
		rwork[i__ * 2] = sn;
#line 294 "zlalsd.f"
	    }
#line 295 "zlalsd.f"
/* L10: */
#line 295 "zlalsd.f"
	}
#line 296 "zlalsd.f"
	if (*nrhs > 1) {
#line 297 "zlalsd.f"
	    i__1 = *nrhs;
#line 297 "zlalsd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "zlalsd.f"
		i__2 = *n - 1;
#line 298 "zlalsd.f"
		for (j = 1; j <= i__2; ++j) {
#line 299 "zlalsd.f"
		    cs = rwork[(j << 1) - 1];
#line 300 "zlalsd.f"
		    sn = rwork[j * 2];
#line 301 "zlalsd.f"
		    zdrot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ 
			    * b_dim1], &c__1, &cs, &sn);
#line 302 "zlalsd.f"
/* L20: */
#line 302 "zlalsd.f"
		}
#line 303 "zlalsd.f"
/* L30: */
#line 303 "zlalsd.f"
	    }
#line 304 "zlalsd.f"
	}
#line 305 "zlalsd.f"
    }

/*     Scale. */

#line 309 "zlalsd.f"
    nm1 = *n - 1;
#line 310 "zlalsd.f"
    orgnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 311 "zlalsd.f"
    if (orgnrm == 0.) {
#line 312 "zlalsd.f"
	zlaset_("A", n, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 313 "zlalsd.f"
	return 0;
#line 314 "zlalsd.f"
    }

#line 316 "zlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 317 "zlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, &nm1, &c__1, &e[1], &nm1, 
	    info, (ftnlen)1);

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

#line 322 "zlalsd.f"
    if (*n <= *smlsiz) {
#line 323 "zlalsd.f"
	irwu = 1;
#line 324 "zlalsd.f"
	irwvt = irwu + *n * *n;
#line 325 "zlalsd.f"
	irwwrk = irwvt + *n * *n;
#line 326 "zlalsd.f"
	irwrb = irwwrk;
#line 327 "zlalsd.f"
	irwib = irwrb + *n * *nrhs;
#line 328 "zlalsd.f"
	irwb = irwib + *n * *nrhs;
#line 329 "zlalsd.f"
	dlaset_("A", n, n, &c_b35, &c_b10, &rwork[irwu], n, (ftnlen)1);
#line 330 "zlalsd.f"
	dlaset_("A", n, n, &c_b35, &c_b10, &rwork[irwvt], n, (ftnlen)1);
#line 331 "zlalsd.f"
	dlasdq_("U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &rwork[irwvt], n, 
		&rwork[irwu], n, &rwork[irwwrk], &c__1, &rwork[irwwrk], info, 
		(ftnlen)1);
#line 334 "zlalsd.f"
	if (*info != 0) {
#line 335 "zlalsd.f"
	    return 0;
#line 336 "zlalsd.f"
	}

/*        In the real version, B is passed to DLASDQ and multiplied */
/*        internally by Q**H. Here B is complex and that product is */
/*        computed below in two steps (real and imaginary parts). */

#line 342 "zlalsd.f"
	j = irwb - 1;
#line 343 "zlalsd.f"
	i__1 = *nrhs;
#line 343 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 344 "zlalsd.f"
	    i__2 = *n;
#line 344 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 345 "zlalsd.f"
		++j;
#line 346 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 346 "zlalsd.f"
		rwork[j] = b[i__3].r;
#line 347 "zlalsd.f"
/* L40: */
#line 347 "zlalsd.f"
	    }
#line 348 "zlalsd.f"
/* L50: */
#line 348 "zlalsd.f"
	}
#line 349 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n,
		 &c_b35, &rwork[irwrb], n, (ftnlen)1, (ftnlen)1);
#line 351 "zlalsd.f"
	j = irwb - 1;
#line 352 "zlalsd.f"
	i__1 = *nrhs;
#line 352 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 353 "zlalsd.f"
	    i__2 = *n;
#line 353 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 354 "zlalsd.f"
		++j;
#line 355 "zlalsd.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 356 "zlalsd.f"
/* L60: */
#line 356 "zlalsd.f"
	    }
#line 357 "zlalsd.f"
/* L70: */
#line 357 "zlalsd.f"
	}
#line 358 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n,
		 &c_b35, &rwork[irwib], n, (ftnlen)1, (ftnlen)1);
#line 360 "zlalsd.f"
	jreal = irwrb - 1;
#line 361 "zlalsd.f"
	jimag = irwib - 1;
#line 362 "zlalsd.f"
	i__1 = *nrhs;
#line 362 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 363 "zlalsd.f"
	    i__2 = *n;
#line 363 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 364 "zlalsd.f"
		++jreal;
#line 365 "zlalsd.f"
		++jimag;
#line 366 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 366 "zlalsd.f"
		i__4 = jreal;
#line 366 "zlalsd.f"
		i__5 = jimag;
#line 366 "zlalsd.f"
		z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 366 "zlalsd.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 368 "zlalsd.f"
/* L80: */
#line 368 "zlalsd.f"
	    }
#line 369 "zlalsd.f"
/* L90: */
#line 369 "zlalsd.f"
	}

#line 371 "zlalsd.f"
	tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));
#line 372 "zlalsd.f"
	i__1 = *n;
#line 372 "zlalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 373 "zlalsd.f"
	    if (d__[i__] <= tol) {
#line 374 "zlalsd.f"
		zlaset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 375 "zlalsd.f"
	    } else {
#line 376 "zlalsd.f"
		zlascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &b[
			i__ + b_dim1], ldb, info, (ftnlen)1);
#line 378 "zlalsd.f"
		++(*rank);
#line 379 "zlalsd.f"
	    }
#line 380 "zlalsd.f"
/* L100: */
#line 380 "zlalsd.f"
	}

/*        Since B is complex, the following call to DGEMM is performed */
/*        in two steps (real and imaginary parts). That is for V * B */
/*        (in the real version of the code V**H is stored in WORK). */

/*        CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, */
/*    $               WORK( NWORK ), N ) */

#line 389 "zlalsd.f"
	j = irwb - 1;
#line 390 "zlalsd.f"
	i__1 = *nrhs;
#line 390 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 391 "zlalsd.f"
	    i__2 = *n;
#line 391 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 392 "zlalsd.f"
		++j;
#line 393 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 393 "zlalsd.f"
		rwork[j] = b[i__3].r;
#line 394 "zlalsd.f"
/* L110: */
#line 394 "zlalsd.f"
	    }
#line 395 "zlalsd.f"
/* L120: */
#line 395 "zlalsd.f"
	}
#line 396 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwrb], n, (ftnlen)1, (ftnlen)1);
#line 398 "zlalsd.f"
	j = irwb - 1;
#line 399 "zlalsd.f"
	i__1 = *nrhs;
#line 399 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 400 "zlalsd.f"
	    i__2 = *n;
#line 400 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 401 "zlalsd.f"
		++j;
#line 402 "zlalsd.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 403 "zlalsd.f"
/* L130: */
#line 403 "zlalsd.f"
	    }
#line 404 "zlalsd.f"
/* L140: */
#line 404 "zlalsd.f"
	}
#line 405 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwib], n, (ftnlen)1, (ftnlen)1);
#line 407 "zlalsd.f"
	jreal = irwrb - 1;
#line 408 "zlalsd.f"
	jimag = irwib - 1;
#line 409 "zlalsd.f"
	i__1 = *nrhs;
#line 409 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 410 "zlalsd.f"
	    i__2 = *n;
#line 410 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 411 "zlalsd.f"
		++jreal;
#line 412 "zlalsd.f"
		++jimag;
#line 413 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 413 "zlalsd.f"
		i__4 = jreal;
#line 413 "zlalsd.f"
		i__5 = jimag;
#line 413 "zlalsd.f"
		z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 413 "zlalsd.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 415 "zlalsd.f"
/* L150: */
#line 415 "zlalsd.f"
	    }
#line 416 "zlalsd.f"
/* L160: */
#line 416 "zlalsd.f"
	}

/*        Unscale. */

#line 420 "zlalsd.f"
	dlascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, 
		info, (ftnlen)1);
#line 421 "zlalsd.f"
	dlasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 422 "zlalsd.f"
	zlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);

#line 424 "zlalsd.f"
	return 0;
#line 425 "zlalsd.f"
    }

/*     Book-keeping and setting up some constants. */

#line 429 "zlalsd.f"
    nlvl = (integer) (log((doublereal) (*n) / (doublereal) (*smlsiz + 1)) / 
	    log(2.)) + 1;

#line 431 "zlalsd.f"
    smlszp = *smlsiz + 1;

#line 433 "zlalsd.f"
    u = 1;
#line 434 "zlalsd.f"
    vt = *smlsiz * *n + 1;
#line 435 "zlalsd.f"
    difl = vt + smlszp * *n;
#line 436 "zlalsd.f"
    difr = difl + nlvl * *n;
#line 437 "zlalsd.f"
    z__ = difr + (nlvl * *n << 1);
#line 438 "zlalsd.f"
    c__ = z__ + nlvl * *n;
#line 439 "zlalsd.f"
    s = c__ + *n;
#line 440 "zlalsd.f"
    poles = s + *n;
#line 441 "zlalsd.f"
    givnum = poles + (nlvl << 1) * *n;
#line 442 "zlalsd.f"
    nrwork = givnum + (nlvl << 1) * *n;
#line 443 "zlalsd.f"
    bx = 1;

#line 445 "zlalsd.f"
    irwrb = nrwork;
#line 446 "zlalsd.f"
    irwib = irwrb + *smlsiz * *nrhs;
#line 447 "zlalsd.f"
    irwb = irwib + *smlsiz * *nrhs;

#line 449 "zlalsd.f"
    sizei = *n + 1;
#line 450 "zlalsd.f"
    k = sizei + *n;
#line 451 "zlalsd.f"
    givptr = k + *n;
#line 452 "zlalsd.f"
    perm = givptr + *n;
#line 453 "zlalsd.f"
    givcol = perm + nlvl * *n;
#line 454 "zlalsd.f"
    iwk = givcol + (nlvl * *n << 1);

#line 456 "zlalsd.f"
    st = 1;
#line 457 "zlalsd.f"
    sqre = 0;
#line 458 "zlalsd.f"
    icmpq1 = 1;
#line 459 "zlalsd.f"
    icmpq2 = 0;
#line 460 "zlalsd.f"
    nsub = 0;

#line 462 "zlalsd.f"
    i__1 = *n;
#line 462 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 463 "zlalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) < eps) {
#line 464 "zlalsd.f"
	    d__[i__] = d_sign(&eps, &d__[i__]);
#line 465 "zlalsd.f"
	}
#line 466 "zlalsd.f"
/* L170: */
#line 466 "zlalsd.f"
    }

#line 468 "zlalsd.f"
    i__1 = nm1;
#line 468 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 469 "zlalsd.f"
	if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
#line 470 "zlalsd.f"
	    ++nsub;
#line 471 "zlalsd.f"
	    iwork[nsub] = st;

/*           Subproblem found. First determine its size and then */
/*           apply divide and conquer on it. */

#line 476 "zlalsd.f"
	    if (i__ < nm1) {

/*              A subproblem with E(I) small for I < NM1. */

#line 480 "zlalsd.f"
		nsize = i__ - st + 1;
#line 481 "zlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 482 "zlalsd.f"
	    } else if ((d__1 = e[i__], abs(d__1)) >= eps) {

/*              A subproblem with E(NM1) not too small but I = NM1. */

#line 486 "zlalsd.f"
		nsize = *n - st + 1;
#line 487 "zlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 488 "zlalsd.f"
	    } else {

/*              A subproblem with E(NM1) small. This implies an */
/*              1-by-1 subproblem at D(N), which is not solved */
/*              explicitly. */

#line 494 "zlalsd.f"
		nsize = i__ - st + 1;
#line 495 "zlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 496 "zlalsd.f"
		++nsub;
#line 497 "zlalsd.f"
		iwork[nsub] = *n;
#line 498 "zlalsd.f"
		iwork[sizei + nsub - 1] = 1;
#line 499 "zlalsd.f"
		zcopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
#line 500 "zlalsd.f"
	    }
#line 501 "zlalsd.f"
	    st1 = st - 1;
#line 502 "zlalsd.f"
	    if (nsize == 1) {

/*              This is a 1-by-1 subproblem and is not solved */
/*              explicitly. */

#line 507 "zlalsd.f"
		zcopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
#line 508 "zlalsd.f"
	    } else if (nsize <= *smlsiz) {

/*              This is a small subproblem and is solved by DLASDQ. */

#line 512 "zlalsd.f"
		dlaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[vt + st1],
			 n, (ftnlen)1);
#line 514 "zlalsd.f"
		dlaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[u + st1], 
			n, (ftnlen)1);
#line 516 "zlalsd.f"
		dlasdq_("U", &c__0, &nsize, &nsize, &nsize, &c__0, &d__[st], &
			e[st], &rwork[vt + st1], n, &rwork[u + st1], n, &
			rwork[nrwork], &c__1, &rwork[nrwork], info, (ftnlen)1)
			;
#line 520 "zlalsd.f"
		if (*info != 0) {
#line 521 "zlalsd.f"
		    return 0;
#line 522 "zlalsd.f"
		}

/*              In the real version, B is passed to DLASDQ and multiplied */
/*              internally by Q**H. Here B is complex and that product is */
/*              computed below in two steps (real and imaginary parts). */

#line 528 "zlalsd.f"
		j = irwb - 1;
#line 529 "zlalsd.f"
		i__2 = *nrhs;
#line 529 "zlalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 530 "zlalsd.f"
		    i__3 = st + nsize - 1;
#line 530 "zlalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 531 "zlalsd.f"
			++j;
#line 532 "zlalsd.f"
			i__4 = jrow + jcol * b_dim1;
#line 532 "zlalsd.f"
			rwork[j] = b[i__4].r;
#line 533 "zlalsd.f"
/* L180: */
#line 533 "zlalsd.f"
		    }
#line 534 "zlalsd.f"
/* L190: */
#line 534 "zlalsd.f"
		}
#line 535 "zlalsd.f"
		dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
			, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &
			nsize, (ftnlen)1, (ftnlen)1);
#line 538 "zlalsd.f"
		j = irwb - 1;
#line 539 "zlalsd.f"
		i__2 = *nrhs;
#line 539 "zlalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 540 "zlalsd.f"
		    i__3 = st + nsize - 1;
#line 540 "zlalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 541 "zlalsd.f"
			++j;
#line 542 "zlalsd.f"
			rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 543 "zlalsd.f"
/* L200: */
#line 543 "zlalsd.f"
		    }
#line 544 "zlalsd.f"
/* L210: */
#line 544 "zlalsd.f"
		}
#line 545 "zlalsd.f"
		dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
			, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &
			nsize, (ftnlen)1, (ftnlen)1);
#line 548 "zlalsd.f"
		jreal = irwrb - 1;
#line 549 "zlalsd.f"
		jimag = irwib - 1;
#line 550 "zlalsd.f"
		i__2 = *nrhs;
#line 550 "zlalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 551 "zlalsd.f"
		    i__3 = st + nsize - 1;
#line 551 "zlalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 552 "zlalsd.f"
			++jreal;
#line 553 "zlalsd.f"
			++jimag;
#line 554 "zlalsd.f"
			i__4 = jrow + jcol * b_dim1;
#line 554 "zlalsd.f"
			i__5 = jreal;
#line 554 "zlalsd.f"
			i__6 = jimag;
#line 554 "zlalsd.f"
			z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 554 "zlalsd.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 556 "zlalsd.f"
/* L220: */
#line 556 "zlalsd.f"
		    }
#line 557 "zlalsd.f"
/* L230: */
#line 557 "zlalsd.f"
		}

#line 559 "zlalsd.f"
		zlacpy_("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + 
			st1], n, (ftnlen)1);
#line 561 "zlalsd.f"
	    } else {

/*              A large problem. Solve it using divide and conquer. */

#line 565 "zlalsd.f"
		dlasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &
			rwork[u + st1], n, &rwork[vt + st1], &iwork[k + st1], 
			&rwork[difl + st1], &rwork[difr + st1], &rwork[z__ + 
			st1], &rwork[poles + st1], &iwork[givptr + st1], &
			iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
			givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &
			rwork[nrwork], &iwork[iwk], info);
#line 574 "zlalsd.f"
		if (*info != 0) {
#line 575 "zlalsd.f"
		    return 0;
#line 576 "zlalsd.f"
		}
#line 577 "zlalsd.f"
		bxst = bx + st1;
#line 578 "zlalsd.f"
		zlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &
			work[bxst], n, &rwork[u + st1], n, &rwork[vt + st1], &
			iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1]
			, &rwork[z__ + st1], &rwork[poles + st1], &iwork[
			givptr + st1], &iwork[givcol + st1], n, &iwork[perm + 
			st1], &rwork[givnum + st1], &rwork[c__ + st1], &rwork[
			s + st1], &rwork[nrwork], &iwork[iwk], info);
#line 587 "zlalsd.f"
		if (*info != 0) {
#line 588 "zlalsd.f"
		    return 0;
#line 589 "zlalsd.f"
		}
#line 590 "zlalsd.f"
	    }
#line 591 "zlalsd.f"
	    st = i__ + 1;
#line 592 "zlalsd.f"
	}
#line 593 "zlalsd.f"
/* L240: */
#line 593 "zlalsd.f"
    }

/*     Apply the singular values and treat the tiny ones as zero. */

#line 597 "zlalsd.f"
    tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));

#line 599 "zlalsd.f"
    i__1 = *n;
#line 599 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Some of the elements in D can be negative because 1-by-1 */
/*        subproblems were not solved explicitly. */

#line 604 "zlalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) <= tol) {
#line 605 "zlalsd.f"
	    zlaset_("A", &c__1, nrhs, &c_b1, &c_b1, &work[bx + i__ - 1], n, (
		    ftnlen)1);
#line 606 "zlalsd.f"
	} else {
#line 607 "zlalsd.f"
	    ++(*rank);
#line 608 "zlalsd.f"
	    zlascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &work[
		    bx + i__ - 1], n, info, (ftnlen)1);
#line 610 "zlalsd.f"
	}
#line 611 "zlalsd.f"
	d__[i__] = (d__1 = d__[i__], abs(d__1));
#line 612 "zlalsd.f"
/* L250: */
#line 612 "zlalsd.f"
    }

/*     Now apply back the right singular vectors. */

#line 616 "zlalsd.f"
    icmpq2 = 1;
#line 617 "zlalsd.f"
    i__1 = nsub;
#line 617 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 618 "zlalsd.f"
	st = iwork[i__];
#line 619 "zlalsd.f"
	st1 = st - 1;
#line 620 "zlalsd.f"
	nsize = iwork[sizei + i__ - 1];
#line 621 "zlalsd.f"
	bxst = bx + st1;
#line 622 "zlalsd.f"
	if (nsize == 1) {
#line 623 "zlalsd.f"
	    zcopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
#line 624 "zlalsd.f"
	} else if (nsize <= *smlsiz) {

/*           Since B and BX are complex, the following call to DGEMM */
/*           is performed in two steps (real and imaginary parts). */

/*           CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE, */
/*    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO, */
/*    $                  B( ST, 1 ), LDB ) */

#line 633 "zlalsd.f"
	    j = bxst - *n - 1;
#line 634 "zlalsd.f"
	    jreal = irwb - 1;
#line 635 "zlalsd.f"
	    i__2 = *nrhs;
#line 635 "zlalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 636 "zlalsd.f"
		j += *n;
#line 637 "zlalsd.f"
		i__3 = nsize;
#line 637 "zlalsd.f"
		for (jrow = 1; jrow <= i__3; ++jrow) {
#line 638 "zlalsd.f"
		    ++jreal;
#line 639 "zlalsd.f"
		    i__4 = j + jrow;
#line 639 "zlalsd.f"
		    rwork[jreal] = work[i__4].r;
#line 640 "zlalsd.f"
/* L260: */
#line 640 "zlalsd.f"
		}
#line 641 "zlalsd.f"
/* L270: */
#line 641 "zlalsd.f"
	    }
#line 642 "zlalsd.f"
	    dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &nsize, (
		    ftnlen)1, (ftnlen)1);
#line 645 "zlalsd.f"
	    j = bxst - *n - 1;
#line 646 "zlalsd.f"
	    jimag = irwb - 1;
#line 647 "zlalsd.f"
	    i__2 = *nrhs;
#line 647 "zlalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 648 "zlalsd.f"
		j += *n;
#line 649 "zlalsd.f"
		i__3 = nsize;
#line 649 "zlalsd.f"
		for (jrow = 1; jrow <= i__3; ++jrow) {
#line 650 "zlalsd.f"
		    ++jimag;
#line 651 "zlalsd.f"
		    rwork[jimag] = d_imag(&work[j + jrow]);
#line 652 "zlalsd.f"
/* L280: */
#line 652 "zlalsd.f"
		}
#line 653 "zlalsd.f"
/* L290: */
#line 653 "zlalsd.f"
	    }
#line 654 "zlalsd.f"
	    dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &nsize, (
		    ftnlen)1, (ftnlen)1);
#line 657 "zlalsd.f"
	    jreal = irwrb - 1;
#line 658 "zlalsd.f"
	    jimag = irwib - 1;
#line 659 "zlalsd.f"
	    i__2 = *nrhs;
#line 659 "zlalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 660 "zlalsd.f"
		i__3 = st + nsize - 1;
#line 660 "zlalsd.f"
		for (jrow = st; jrow <= i__3; ++jrow) {
#line 661 "zlalsd.f"
		    ++jreal;
#line 662 "zlalsd.f"
		    ++jimag;
#line 663 "zlalsd.f"
		    i__4 = jrow + jcol * b_dim1;
#line 663 "zlalsd.f"
		    i__5 = jreal;
#line 663 "zlalsd.f"
		    i__6 = jimag;
#line 663 "zlalsd.f"
		    z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 663 "zlalsd.f"
		    b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 665 "zlalsd.f"
/* L300: */
#line 665 "zlalsd.f"
		}
#line 666 "zlalsd.f"
/* L310: */
#line 666 "zlalsd.f"
	    }
#line 667 "zlalsd.f"
	} else {
#line 668 "zlalsd.f"
	    zlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + 
		    b_dim1], ldb, &rwork[u + st1], n, &rwork[vt + st1], &
		    iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1], &
		    rwork[z__ + st1], &rwork[poles + st1], &iwork[givptr + 
		    st1], &iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
		    givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &rwork[
		    nrwork], &iwork[iwk], info);
#line 677 "zlalsd.f"
	    if (*info != 0) {
#line 678 "zlalsd.f"
		return 0;
#line 679 "zlalsd.f"
	    }
#line 680 "zlalsd.f"
	}
#line 681 "zlalsd.f"
/* L320: */
#line 681 "zlalsd.f"
    }

/*     Unscale and sort the singular values. */

#line 685 "zlalsd.f"
    dlascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 686 "zlalsd.f"
    dlasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 687 "zlalsd.f"
    zlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], ldb, 
	    info, (ftnlen)1);

#line 689 "zlalsd.f"
    return 0;

/*     End of ZLALSD */

} /* zlalsd_ */

