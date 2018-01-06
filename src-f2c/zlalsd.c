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
/* >          WORK is COMPLEX*16 array, dimension at least */
/* >         (N * NRHS). */
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

/* > \date December 2016 */

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 241 "zlalsd.f"
    /* Parameter adjustments */
#line 241 "zlalsd.f"
    --d__;
#line 241 "zlalsd.f"
    --e;
#line 241 "zlalsd.f"
    b_dim1 = *ldb;
#line 241 "zlalsd.f"
    b_offset = 1 + b_dim1;
#line 241 "zlalsd.f"
    b -= b_offset;
#line 241 "zlalsd.f"
    --work;
#line 241 "zlalsd.f"
    --rwork;
#line 241 "zlalsd.f"
    --iwork;
#line 241 "zlalsd.f"

#line 241 "zlalsd.f"
    /* Function Body */
#line 241 "zlalsd.f"
    *info = 0;

#line 243 "zlalsd.f"
    if (*n < 0) {
#line 244 "zlalsd.f"
	*info = -3;
#line 245 "zlalsd.f"
    } else if (*nrhs < 1) {
#line 246 "zlalsd.f"
	*info = -4;
#line 247 "zlalsd.f"
    } else if (*ldb < 1 || *ldb < *n) {
#line 248 "zlalsd.f"
	*info = -8;
#line 249 "zlalsd.f"
    }
#line 250 "zlalsd.f"
    if (*info != 0) {
#line 251 "zlalsd.f"
	i__1 = -(*info);
#line 251 "zlalsd.f"
	xerbla_("ZLALSD", &i__1, (ftnlen)6);
#line 252 "zlalsd.f"
	return 0;
#line 253 "zlalsd.f"
    }

#line 255 "zlalsd.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Set up the tolerance. */

#line 259 "zlalsd.f"
    if (*rcond <= 0. || *rcond >= 1.) {
#line 260 "zlalsd.f"
	rcnd = eps;
#line 261 "zlalsd.f"
    } else {
#line 262 "zlalsd.f"
	rcnd = *rcond;
#line 263 "zlalsd.f"
    }

#line 265 "zlalsd.f"
    *rank = 0;

/*     Quick return if possible. */

#line 269 "zlalsd.f"
    if (*n == 0) {
#line 270 "zlalsd.f"
	return 0;
#line 271 "zlalsd.f"
    } else if (*n == 1) {
#line 272 "zlalsd.f"
	if (d__[1] == 0.) {
#line 273 "zlalsd.f"
	    zlaset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (
		    ftnlen)1);
#line 274 "zlalsd.f"
	} else {
#line 275 "zlalsd.f"
	    *rank = 1;
#line 276 "zlalsd.f"
	    zlascl_("G", &c__0, &c__0, &d__[1], &c_b10, &c__1, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)1);
#line 277 "zlalsd.f"
	    d__[1] = abs(d__[1]);
#line 278 "zlalsd.f"
	}
#line 279 "zlalsd.f"
	return 0;
#line 280 "zlalsd.f"
    }

/*     Rotate the matrix if it is lower bidiagonal. */

#line 284 "zlalsd.f"
    if (*(unsigned char *)uplo == 'L') {
#line 285 "zlalsd.f"
	i__1 = *n - 1;
#line 285 "zlalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 286 "zlalsd.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 287 "zlalsd.f"
	    d__[i__] = r__;
#line 288 "zlalsd.f"
	    e[i__] = sn * d__[i__ + 1];
#line 289 "zlalsd.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 290 "zlalsd.f"
	    if (*nrhs == 1) {
#line 291 "zlalsd.f"
		zdrot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &
			c__1, &cs, &sn);
#line 292 "zlalsd.f"
	    } else {
#line 293 "zlalsd.f"
		rwork[(i__ << 1) - 1] = cs;
#line 294 "zlalsd.f"
		rwork[i__ * 2] = sn;
#line 295 "zlalsd.f"
	    }
#line 296 "zlalsd.f"
/* L10: */
#line 296 "zlalsd.f"
	}
#line 297 "zlalsd.f"
	if (*nrhs > 1) {
#line 298 "zlalsd.f"
	    i__1 = *nrhs;
#line 298 "zlalsd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "zlalsd.f"
		i__2 = *n - 1;
#line 299 "zlalsd.f"
		for (j = 1; j <= i__2; ++j) {
#line 300 "zlalsd.f"
		    cs = rwork[(j << 1) - 1];
#line 301 "zlalsd.f"
		    sn = rwork[j * 2];
#line 302 "zlalsd.f"
		    zdrot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ 
			    * b_dim1], &c__1, &cs, &sn);
#line 303 "zlalsd.f"
/* L20: */
#line 303 "zlalsd.f"
		}
#line 304 "zlalsd.f"
/* L30: */
#line 304 "zlalsd.f"
	    }
#line 305 "zlalsd.f"
	}
#line 306 "zlalsd.f"
    }

/*     Scale. */

#line 310 "zlalsd.f"
    nm1 = *n - 1;
#line 311 "zlalsd.f"
    orgnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 312 "zlalsd.f"
    if (orgnrm == 0.) {
#line 313 "zlalsd.f"
	zlaset_("A", n, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 314 "zlalsd.f"
	return 0;
#line 315 "zlalsd.f"
    }

#line 317 "zlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 318 "zlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, &nm1, &c__1, &e[1], &nm1, 
	    info, (ftnlen)1);

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

#line 323 "zlalsd.f"
    if (*n <= *smlsiz) {
#line 324 "zlalsd.f"
	irwu = 1;
#line 325 "zlalsd.f"
	irwvt = irwu + *n * *n;
#line 326 "zlalsd.f"
	irwwrk = irwvt + *n * *n;
#line 327 "zlalsd.f"
	irwrb = irwwrk;
#line 328 "zlalsd.f"
	irwib = irwrb + *n * *nrhs;
#line 329 "zlalsd.f"
	irwb = irwib + *n * *nrhs;
#line 330 "zlalsd.f"
	dlaset_("A", n, n, &c_b35, &c_b10, &rwork[irwu], n, (ftnlen)1);
#line 331 "zlalsd.f"
	dlaset_("A", n, n, &c_b35, &c_b10, &rwork[irwvt], n, (ftnlen)1);
#line 332 "zlalsd.f"
	dlasdq_("U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &rwork[irwvt], n, 
		&rwork[irwu], n, &rwork[irwwrk], &c__1, &rwork[irwwrk], info, 
		(ftnlen)1);
#line 335 "zlalsd.f"
	if (*info != 0) {
#line 336 "zlalsd.f"
	    return 0;
#line 337 "zlalsd.f"
	}

/*        In the real version, B is passed to DLASDQ and multiplied */
/*        internally by Q**H. Here B is complex and that product is */
/*        computed below in two steps (real and imaginary parts). */

#line 343 "zlalsd.f"
	j = irwb - 1;
#line 344 "zlalsd.f"
	i__1 = *nrhs;
#line 344 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 345 "zlalsd.f"
	    i__2 = *n;
#line 345 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 346 "zlalsd.f"
		++j;
#line 347 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 347 "zlalsd.f"
		rwork[j] = b[i__3].r;
#line 348 "zlalsd.f"
/* L40: */
#line 348 "zlalsd.f"
	    }
#line 349 "zlalsd.f"
/* L50: */
#line 349 "zlalsd.f"
	}
#line 350 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n,
		 &c_b35, &rwork[irwrb], n, (ftnlen)1, (ftnlen)1);
#line 352 "zlalsd.f"
	j = irwb - 1;
#line 353 "zlalsd.f"
	i__1 = *nrhs;
#line 353 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 354 "zlalsd.f"
	    i__2 = *n;
#line 354 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 355 "zlalsd.f"
		++j;
#line 356 "zlalsd.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 357 "zlalsd.f"
/* L60: */
#line 357 "zlalsd.f"
	    }
#line 358 "zlalsd.f"
/* L70: */
#line 358 "zlalsd.f"
	}
#line 359 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n,
		 &c_b35, &rwork[irwib], n, (ftnlen)1, (ftnlen)1);
#line 361 "zlalsd.f"
	jreal = irwrb - 1;
#line 362 "zlalsd.f"
	jimag = irwib - 1;
#line 363 "zlalsd.f"
	i__1 = *nrhs;
#line 363 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 364 "zlalsd.f"
	    i__2 = *n;
#line 364 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 365 "zlalsd.f"
		++jreal;
#line 366 "zlalsd.f"
		++jimag;
#line 367 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 367 "zlalsd.f"
		i__4 = jreal;
#line 367 "zlalsd.f"
		i__5 = jimag;
#line 367 "zlalsd.f"
		z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 367 "zlalsd.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 369 "zlalsd.f"
/* L80: */
#line 369 "zlalsd.f"
	    }
#line 370 "zlalsd.f"
/* L90: */
#line 370 "zlalsd.f"
	}

#line 372 "zlalsd.f"
	tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));
#line 373 "zlalsd.f"
	i__1 = *n;
#line 373 "zlalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 374 "zlalsd.f"
	    if (d__[i__] <= tol) {
#line 375 "zlalsd.f"
		zlaset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 376 "zlalsd.f"
	    } else {
#line 377 "zlalsd.f"
		zlascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &b[
			i__ + b_dim1], ldb, info, (ftnlen)1);
#line 379 "zlalsd.f"
		++(*rank);
#line 380 "zlalsd.f"
	    }
#line 381 "zlalsd.f"
/* L100: */
#line 381 "zlalsd.f"
	}

/*        Since B is complex, the following call to DGEMM is performed */
/*        in two steps (real and imaginary parts). That is for V * B */
/*        (in the real version of the code V**H is stored in WORK). */

/*        CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, */
/*    $               WORK( NWORK ), N ) */

#line 390 "zlalsd.f"
	j = irwb - 1;
#line 391 "zlalsd.f"
	i__1 = *nrhs;
#line 391 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 392 "zlalsd.f"
	    i__2 = *n;
#line 392 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 393 "zlalsd.f"
		++j;
#line 394 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 394 "zlalsd.f"
		rwork[j] = b[i__3].r;
#line 395 "zlalsd.f"
/* L110: */
#line 395 "zlalsd.f"
	    }
#line 396 "zlalsd.f"
/* L120: */
#line 396 "zlalsd.f"
	}
#line 397 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwrb], n, (ftnlen)1, (ftnlen)1);
#line 399 "zlalsd.f"
	j = irwb - 1;
#line 400 "zlalsd.f"
	i__1 = *nrhs;
#line 400 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 401 "zlalsd.f"
	    i__2 = *n;
#line 401 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 402 "zlalsd.f"
		++j;
#line 403 "zlalsd.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 404 "zlalsd.f"
/* L130: */
#line 404 "zlalsd.f"
	    }
#line 405 "zlalsd.f"
/* L140: */
#line 405 "zlalsd.f"
	}
#line 406 "zlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwib], n, (ftnlen)1, (ftnlen)1);
#line 408 "zlalsd.f"
	jreal = irwrb - 1;
#line 409 "zlalsd.f"
	jimag = irwib - 1;
#line 410 "zlalsd.f"
	i__1 = *nrhs;
#line 410 "zlalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 411 "zlalsd.f"
	    i__2 = *n;
#line 411 "zlalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 412 "zlalsd.f"
		++jreal;
#line 413 "zlalsd.f"
		++jimag;
#line 414 "zlalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 414 "zlalsd.f"
		i__4 = jreal;
#line 414 "zlalsd.f"
		i__5 = jimag;
#line 414 "zlalsd.f"
		z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 414 "zlalsd.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 416 "zlalsd.f"
/* L150: */
#line 416 "zlalsd.f"
	    }
#line 417 "zlalsd.f"
/* L160: */
#line 417 "zlalsd.f"
	}

/*        Unscale. */

#line 421 "zlalsd.f"
	dlascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, 
		info, (ftnlen)1);
#line 422 "zlalsd.f"
	dlasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 423 "zlalsd.f"
	zlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);

#line 425 "zlalsd.f"
	return 0;
#line 426 "zlalsd.f"
    }

/*     Book-keeping and setting up some constants. */

#line 430 "zlalsd.f"
    nlvl = (integer) (log((doublereal) (*n) / (doublereal) (*smlsiz + 1)) / 
	    log(2.)) + 1;

#line 432 "zlalsd.f"
    smlszp = *smlsiz + 1;

#line 434 "zlalsd.f"
    u = 1;
#line 435 "zlalsd.f"
    vt = *smlsiz * *n + 1;
#line 436 "zlalsd.f"
    difl = vt + smlszp * *n;
#line 437 "zlalsd.f"
    difr = difl + nlvl * *n;
#line 438 "zlalsd.f"
    z__ = difr + (nlvl * *n << 1);
#line 439 "zlalsd.f"
    c__ = z__ + nlvl * *n;
#line 440 "zlalsd.f"
    s = c__ + *n;
#line 441 "zlalsd.f"
    poles = s + *n;
#line 442 "zlalsd.f"
    givnum = poles + (nlvl << 1) * *n;
#line 443 "zlalsd.f"
    nrwork = givnum + (nlvl << 1) * *n;
#line 444 "zlalsd.f"
    bx = 1;

#line 446 "zlalsd.f"
    irwrb = nrwork;
#line 447 "zlalsd.f"
    irwib = irwrb + *smlsiz * *nrhs;
#line 448 "zlalsd.f"
    irwb = irwib + *smlsiz * *nrhs;

#line 450 "zlalsd.f"
    sizei = *n + 1;
#line 451 "zlalsd.f"
    k = sizei + *n;
#line 452 "zlalsd.f"
    givptr = k + *n;
#line 453 "zlalsd.f"
    perm = givptr + *n;
#line 454 "zlalsd.f"
    givcol = perm + nlvl * *n;
#line 455 "zlalsd.f"
    iwk = givcol + (nlvl * *n << 1);

#line 457 "zlalsd.f"
    st = 1;
#line 458 "zlalsd.f"
    sqre = 0;
#line 459 "zlalsd.f"
    icmpq1 = 1;
#line 460 "zlalsd.f"
    icmpq2 = 0;
#line 461 "zlalsd.f"
    nsub = 0;

#line 463 "zlalsd.f"
    i__1 = *n;
#line 463 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 464 "zlalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) < eps) {
#line 465 "zlalsd.f"
	    d__[i__] = d_sign(&eps, &d__[i__]);
#line 466 "zlalsd.f"
	}
#line 467 "zlalsd.f"
/* L170: */
#line 467 "zlalsd.f"
    }

#line 469 "zlalsd.f"
    i__1 = nm1;
#line 469 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 470 "zlalsd.f"
	if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
#line 471 "zlalsd.f"
	    ++nsub;
#line 472 "zlalsd.f"
	    iwork[nsub] = st;

/*           Subproblem found. First determine its size and then */
/*           apply divide and conquer on it. */

#line 477 "zlalsd.f"
	    if (i__ < nm1) {

/*              A subproblem with E(I) small for I < NM1. */

#line 481 "zlalsd.f"
		nsize = i__ - st + 1;
#line 482 "zlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 483 "zlalsd.f"
	    } else if ((d__1 = e[i__], abs(d__1)) >= eps) {

/*              A subproblem with E(NM1) not too small but I = NM1. */

#line 487 "zlalsd.f"
		nsize = *n - st + 1;
#line 488 "zlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 489 "zlalsd.f"
	    } else {

/*              A subproblem with E(NM1) small. This implies an */
/*              1-by-1 subproblem at D(N), which is not solved */
/*              explicitly. */

#line 495 "zlalsd.f"
		nsize = i__ - st + 1;
#line 496 "zlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 497 "zlalsd.f"
		++nsub;
#line 498 "zlalsd.f"
		iwork[nsub] = *n;
#line 499 "zlalsd.f"
		iwork[sizei + nsub - 1] = 1;
#line 500 "zlalsd.f"
		zcopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
#line 501 "zlalsd.f"
	    }
#line 502 "zlalsd.f"
	    st1 = st - 1;
#line 503 "zlalsd.f"
	    if (nsize == 1) {

/*              This is a 1-by-1 subproblem and is not solved */
/*              explicitly. */

#line 508 "zlalsd.f"
		zcopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
#line 509 "zlalsd.f"
	    } else if (nsize <= *smlsiz) {

/*              This is a small subproblem and is solved by DLASDQ. */

#line 513 "zlalsd.f"
		dlaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[vt + st1],
			 n, (ftnlen)1);
#line 515 "zlalsd.f"
		dlaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[u + st1], 
			n, (ftnlen)1);
#line 517 "zlalsd.f"
		dlasdq_("U", &c__0, &nsize, &nsize, &nsize, &c__0, &d__[st], &
			e[st], &rwork[vt + st1], n, &rwork[u + st1], n, &
			rwork[nrwork], &c__1, &rwork[nrwork], info, (ftnlen)1)
			;
#line 521 "zlalsd.f"
		if (*info != 0) {
#line 522 "zlalsd.f"
		    return 0;
#line 523 "zlalsd.f"
		}

/*              In the real version, B is passed to DLASDQ and multiplied */
/*              internally by Q**H. Here B is complex and that product is */
/*              computed below in two steps (real and imaginary parts). */

#line 529 "zlalsd.f"
		j = irwb - 1;
#line 530 "zlalsd.f"
		i__2 = *nrhs;
#line 530 "zlalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 531 "zlalsd.f"
		    i__3 = st + nsize - 1;
#line 531 "zlalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 532 "zlalsd.f"
			++j;
#line 533 "zlalsd.f"
			i__4 = jrow + jcol * b_dim1;
#line 533 "zlalsd.f"
			rwork[j] = b[i__4].r;
#line 534 "zlalsd.f"
/* L180: */
#line 534 "zlalsd.f"
		    }
#line 535 "zlalsd.f"
/* L190: */
#line 535 "zlalsd.f"
		}
#line 536 "zlalsd.f"
		dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
			, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &
			nsize, (ftnlen)1, (ftnlen)1);
#line 539 "zlalsd.f"
		j = irwb - 1;
#line 540 "zlalsd.f"
		i__2 = *nrhs;
#line 540 "zlalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 541 "zlalsd.f"
		    i__3 = st + nsize - 1;
#line 541 "zlalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 542 "zlalsd.f"
			++j;
#line 543 "zlalsd.f"
			rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 544 "zlalsd.f"
/* L200: */
#line 544 "zlalsd.f"
		    }
#line 545 "zlalsd.f"
/* L210: */
#line 545 "zlalsd.f"
		}
#line 546 "zlalsd.f"
		dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
			, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &
			nsize, (ftnlen)1, (ftnlen)1);
#line 549 "zlalsd.f"
		jreal = irwrb - 1;
#line 550 "zlalsd.f"
		jimag = irwib - 1;
#line 551 "zlalsd.f"
		i__2 = *nrhs;
#line 551 "zlalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 552 "zlalsd.f"
		    i__3 = st + nsize - 1;
#line 552 "zlalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 553 "zlalsd.f"
			++jreal;
#line 554 "zlalsd.f"
			++jimag;
#line 555 "zlalsd.f"
			i__4 = jrow + jcol * b_dim1;
#line 555 "zlalsd.f"
			i__5 = jreal;
#line 555 "zlalsd.f"
			i__6 = jimag;
#line 555 "zlalsd.f"
			z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 555 "zlalsd.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 557 "zlalsd.f"
/* L220: */
#line 557 "zlalsd.f"
		    }
#line 558 "zlalsd.f"
/* L230: */
#line 558 "zlalsd.f"
		}

#line 560 "zlalsd.f"
		zlacpy_("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + 
			st1], n, (ftnlen)1);
#line 562 "zlalsd.f"
	    } else {

/*              A large problem. Solve it using divide and conquer. */

#line 566 "zlalsd.f"
		dlasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &
			rwork[u + st1], n, &rwork[vt + st1], &iwork[k + st1], 
			&rwork[difl + st1], &rwork[difr + st1], &rwork[z__ + 
			st1], &rwork[poles + st1], &iwork[givptr + st1], &
			iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
			givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &
			rwork[nrwork], &iwork[iwk], info);
#line 575 "zlalsd.f"
		if (*info != 0) {
#line 576 "zlalsd.f"
		    return 0;
#line 577 "zlalsd.f"
		}
#line 578 "zlalsd.f"
		bxst = bx + st1;
#line 579 "zlalsd.f"
		zlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &
			work[bxst], n, &rwork[u + st1], n, &rwork[vt + st1], &
			iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1]
			, &rwork[z__ + st1], &rwork[poles + st1], &iwork[
			givptr + st1], &iwork[givcol + st1], n, &iwork[perm + 
			st1], &rwork[givnum + st1], &rwork[c__ + st1], &rwork[
			s + st1], &rwork[nrwork], &iwork[iwk], info);
#line 588 "zlalsd.f"
		if (*info != 0) {
#line 589 "zlalsd.f"
		    return 0;
#line 590 "zlalsd.f"
		}
#line 591 "zlalsd.f"
	    }
#line 592 "zlalsd.f"
	    st = i__ + 1;
#line 593 "zlalsd.f"
	}
#line 594 "zlalsd.f"
/* L240: */
#line 594 "zlalsd.f"
    }

/*     Apply the singular values and treat the tiny ones as zero. */

#line 598 "zlalsd.f"
    tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));

#line 600 "zlalsd.f"
    i__1 = *n;
#line 600 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Some of the elements in D can be negative because 1-by-1 */
/*        subproblems were not solved explicitly. */

#line 605 "zlalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) <= tol) {
#line 606 "zlalsd.f"
	    zlaset_("A", &c__1, nrhs, &c_b1, &c_b1, &work[bx + i__ - 1], n, (
		    ftnlen)1);
#line 607 "zlalsd.f"
	} else {
#line 608 "zlalsd.f"
	    ++(*rank);
#line 609 "zlalsd.f"
	    zlascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &work[
		    bx + i__ - 1], n, info, (ftnlen)1);
#line 611 "zlalsd.f"
	}
#line 612 "zlalsd.f"
	d__[i__] = (d__1 = d__[i__], abs(d__1));
#line 613 "zlalsd.f"
/* L250: */
#line 613 "zlalsd.f"
    }

/*     Now apply back the right singular vectors. */

#line 617 "zlalsd.f"
    icmpq2 = 1;
#line 618 "zlalsd.f"
    i__1 = nsub;
#line 618 "zlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 619 "zlalsd.f"
	st = iwork[i__];
#line 620 "zlalsd.f"
	st1 = st - 1;
#line 621 "zlalsd.f"
	nsize = iwork[sizei + i__ - 1];
#line 622 "zlalsd.f"
	bxst = bx + st1;
#line 623 "zlalsd.f"
	if (nsize == 1) {
#line 624 "zlalsd.f"
	    zcopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
#line 625 "zlalsd.f"
	} else if (nsize <= *smlsiz) {

/*           Since B and BX are complex, the following call to DGEMM */
/*           is performed in two steps (real and imaginary parts). */

/*           CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE, */
/*    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO, */
/*    $                  B( ST, 1 ), LDB ) */

#line 634 "zlalsd.f"
	    j = bxst - *n - 1;
#line 635 "zlalsd.f"
	    jreal = irwb - 1;
#line 636 "zlalsd.f"
	    i__2 = *nrhs;
#line 636 "zlalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 637 "zlalsd.f"
		j += *n;
#line 638 "zlalsd.f"
		i__3 = nsize;
#line 638 "zlalsd.f"
		for (jrow = 1; jrow <= i__3; ++jrow) {
#line 639 "zlalsd.f"
		    ++jreal;
#line 640 "zlalsd.f"
		    i__4 = j + jrow;
#line 640 "zlalsd.f"
		    rwork[jreal] = work[i__4].r;
#line 641 "zlalsd.f"
/* L260: */
#line 641 "zlalsd.f"
		}
#line 642 "zlalsd.f"
/* L270: */
#line 642 "zlalsd.f"
	    }
#line 643 "zlalsd.f"
	    dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &nsize, (
		    ftnlen)1, (ftnlen)1);
#line 646 "zlalsd.f"
	    j = bxst - *n - 1;
#line 647 "zlalsd.f"
	    jimag = irwb - 1;
#line 648 "zlalsd.f"
	    i__2 = *nrhs;
#line 648 "zlalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 649 "zlalsd.f"
		j += *n;
#line 650 "zlalsd.f"
		i__3 = nsize;
#line 650 "zlalsd.f"
		for (jrow = 1; jrow <= i__3; ++jrow) {
#line 651 "zlalsd.f"
		    ++jimag;
#line 652 "zlalsd.f"
		    rwork[jimag] = d_imag(&work[j + jrow]);
#line 653 "zlalsd.f"
/* L280: */
#line 653 "zlalsd.f"
		}
#line 654 "zlalsd.f"
/* L290: */
#line 654 "zlalsd.f"
	    }
#line 655 "zlalsd.f"
	    dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &nsize, (
		    ftnlen)1, (ftnlen)1);
#line 658 "zlalsd.f"
	    jreal = irwrb - 1;
#line 659 "zlalsd.f"
	    jimag = irwib - 1;
#line 660 "zlalsd.f"
	    i__2 = *nrhs;
#line 660 "zlalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 661 "zlalsd.f"
		i__3 = st + nsize - 1;
#line 661 "zlalsd.f"
		for (jrow = st; jrow <= i__3; ++jrow) {
#line 662 "zlalsd.f"
		    ++jreal;
#line 663 "zlalsd.f"
		    ++jimag;
#line 664 "zlalsd.f"
		    i__4 = jrow + jcol * b_dim1;
#line 664 "zlalsd.f"
		    i__5 = jreal;
#line 664 "zlalsd.f"
		    i__6 = jimag;
#line 664 "zlalsd.f"
		    z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 664 "zlalsd.f"
		    b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 666 "zlalsd.f"
/* L300: */
#line 666 "zlalsd.f"
		}
#line 667 "zlalsd.f"
/* L310: */
#line 667 "zlalsd.f"
	    }
#line 668 "zlalsd.f"
	} else {
#line 669 "zlalsd.f"
	    zlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + 
		    b_dim1], ldb, &rwork[u + st1], n, &rwork[vt + st1], &
		    iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1], &
		    rwork[z__ + st1], &rwork[poles + st1], &iwork[givptr + 
		    st1], &iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
		    givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &rwork[
		    nrwork], &iwork[iwk], info);
#line 678 "zlalsd.f"
	    if (*info != 0) {
#line 679 "zlalsd.f"
		return 0;
#line 680 "zlalsd.f"
	    }
#line 681 "zlalsd.f"
	}
#line 682 "zlalsd.f"
/* L320: */
#line 682 "zlalsd.f"
    }

/*     Unscale and sort the singular values. */

#line 686 "zlalsd.f"
    dlascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 687 "zlalsd.f"
    dlasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 688 "zlalsd.f"
    zlascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], ldb, 
	    info, (ftnlen)1);

#line 690 "zlalsd.f"
    return 0;

/*     End of ZLALSD */

} /* zlalsd_ */

