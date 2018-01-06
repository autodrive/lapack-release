#line 1 "clalsd.f"
/* clalsd.f -- translated by f2c (version 20100827).
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

#line 1 "clalsd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b10 = 1.;
static doublereal c_b35 = 0.;

/* > \brief \b CLALSD uses the singular value decomposition of A to solve the least squares problem. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLALSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, */
/*                          RANK, WORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), E( * ), RWORK( * ) */
/*       COMPLEX            B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLALSD uses the singular value decomposition of A to solve the least */
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
/* >          D is REAL array, dimension (N) */
/* >         On entry D contains the main diagonal of the bidiagonal */
/* >         matrix. On exit, if INFO = 0, D contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >         Contains the super-diagonal entries of the bidiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          RCOND is REAL */
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
/* >          WORK is COMPLEX array, dimension (N * NRHS). */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension at least */
/* >         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + */
/* >         MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ), */
/* >         where */
/* >         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (3*N*NLVL + 11*N). */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int clalsd_(char *uplo, integer *smlsiz, integer *n, integer 
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
	    jimag, jreal;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer irwib;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer poles, sizei, irwrb, nsize;
    extern /* Subroutine */ int csrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *);
    static integer irwvt, icmpq1, icmpq2;
    extern /* Subroutine */ int clalsa_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), clascl_(char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int slasda_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), clacpy_(char *, integer *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, ftnlen), claset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), slascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static integer givcol;
    extern /* Subroutine */ int slasdq_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), slaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), slartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal orgnrm;
    static integer givnum;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int slasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static integer givptr, nrwork, irwwrk, smlszp;


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

#line 239 "clalsd.f"
    /* Parameter adjustments */
#line 239 "clalsd.f"
    --d__;
#line 239 "clalsd.f"
    --e;
#line 239 "clalsd.f"
    b_dim1 = *ldb;
#line 239 "clalsd.f"
    b_offset = 1 + b_dim1;
#line 239 "clalsd.f"
    b -= b_offset;
#line 239 "clalsd.f"
    --work;
#line 239 "clalsd.f"
    --rwork;
#line 239 "clalsd.f"
    --iwork;
#line 239 "clalsd.f"

#line 239 "clalsd.f"
    /* Function Body */
#line 239 "clalsd.f"
    *info = 0;

#line 241 "clalsd.f"
    if (*n < 0) {
#line 242 "clalsd.f"
	*info = -3;
#line 243 "clalsd.f"
    } else if (*nrhs < 1) {
#line 244 "clalsd.f"
	*info = -4;
#line 245 "clalsd.f"
    } else if (*ldb < 1 || *ldb < *n) {
#line 246 "clalsd.f"
	*info = -8;
#line 247 "clalsd.f"
    }
#line 248 "clalsd.f"
    if (*info != 0) {
#line 249 "clalsd.f"
	i__1 = -(*info);
#line 249 "clalsd.f"
	xerbla_("CLALSD", &i__1, (ftnlen)6);
#line 250 "clalsd.f"
	return 0;
#line 251 "clalsd.f"
    }

#line 253 "clalsd.f"
    eps = slamch_("Epsilon", (ftnlen)7);

/*     Set up the tolerance. */

#line 257 "clalsd.f"
    if (*rcond <= 0. || *rcond >= 1.) {
#line 258 "clalsd.f"
	rcnd = eps;
#line 259 "clalsd.f"
    } else {
#line 260 "clalsd.f"
	rcnd = *rcond;
#line 261 "clalsd.f"
    }

#line 263 "clalsd.f"
    *rank = 0;

/*     Quick return if possible. */

#line 267 "clalsd.f"
    if (*n == 0) {
#line 268 "clalsd.f"
	return 0;
#line 269 "clalsd.f"
    } else if (*n == 1) {
#line 270 "clalsd.f"
	if (d__[1] == 0.) {
#line 271 "clalsd.f"
	    claset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (
		    ftnlen)1);
#line 272 "clalsd.f"
	} else {
#line 273 "clalsd.f"
	    *rank = 1;
#line 274 "clalsd.f"
	    clascl_("G", &c__0, &c__0, &d__[1], &c_b10, &c__1, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)1);
#line 275 "clalsd.f"
	    d__[1] = abs(d__[1]);
#line 276 "clalsd.f"
	}
#line 277 "clalsd.f"
	return 0;
#line 278 "clalsd.f"
    }

/*     Rotate the matrix if it is lower bidiagonal. */

#line 282 "clalsd.f"
    if (*(unsigned char *)uplo == 'L') {
#line 283 "clalsd.f"
	i__1 = *n - 1;
#line 283 "clalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 284 "clalsd.f"
	    slartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 285 "clalsd.f"
	    d__[i__] = r__;
#line 286 "clalsd.f"
	    e[i__] = sn * d__[i__ + 1];
#line 287 "clalsd.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 288 "clalsd.f"
	    if (*nrhs == 1) {
#line 289 "clalsd.f"
		csrot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &
			c__1, &cs, &sn);
#line 290 "clalsd.f"
	    } else {
#line 291 "clalsd.f"
		rwork[(i__ << 1) - 1] = cs;
#line 292 "clalsd.f"
		rwork[i__ * 2] = sn;
#line 293 "clalsd.f"
	    }
#line 294 "clalsd.f"
/* L10: */
#line 294 "clalsd.f"
	}
#line 295 "clalsd.f"
	if (*nrhs > 1) {
#line 296 "clalsd.f"
	    i__1 = *nrhs;
#line 296 "clalsd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "clalsd.f"
		i__2 = *n - 1;
#line 297 "clalsd.f"
		for (j = 1; j <= i__2; ++j) {
#line 298 "clalsd.f"
		    cs = rwork[(j << 1) - 1];
#line 299 "clalsd.f"
		    sn = rwork[j * 2];
#line 300 "clalsd.f"
		    csrot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ 
			    * b_dim1], &c__1, &cs, &sn);
#line 301 "clalsd.f"
/* L20: */
#line 301 "clalsd.f"
		}
#line 302 "clalsd.f"
/* L30: */
#line 302 "clalsd.f"
	    }
#line 303 "clalsd.f"
	}
#line 304 "clalsd.f"
    }

/*     Scale. */

#line 308 "clalsd.f"
    nm1 = *n - 1;
#line 309 "clalsd.f"
    orgnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 310 "clalsd.f"
    if (orgnrm == 0.) {
#line 311 "clalsd.f"
	claset_("A", n, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 312 "clalsd.f"
	return 0;
#line 313 "clalsd.f"
    }

#line 315 "clalsd.f"
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 316 "clalsd.f"
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b10, &nm1, &c__1, &e[1], &nm1, 
	    info, (ftnlen)1);

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

#line 321 "clalsd.f"
    if (*n <= *smlsiz) {
#line 322 "clalsd.f"
	irwu = 1;
#line 323 "clalsd.f"
	irwvt = irwu + *n * *n;
#line 324 "clalsd.f"
	irwwrk = irwvt + *n * *n;
#line 325 "clalsd.f"
	irwrb = irwwrk;
#line 326 "clalsd.f"
	irwib = irwrb + *n * *nrhs;
#line 327 "clalsd.f"
	irwb = irwib + *n * *nrhs;
#line 328 "clalsd.f"
	slaset_("A", n, n, &c_b35, &c_b10, &rwork[irwu], n, (ftnlen)1);
#line 329 "clalsd.f"
	slaset_("A", n, n, &c_b35, &c_b10, &rwork[irwvt], n, (ftnlen)1);
#line 330 "clalsd.f"
	slasdq_("U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &rwork[irwvt], n, 
		&rwork[irwu], n, &rwork[irwwrk], &c__1, &rwork[irwwrk], info, 
		(ftnlen)1);
#line 333 "clalsd.f"
	if (*info != 0) {
#line 334 "clalsd.f"
	    return 0;
#line 335 "clalsd.f"
	}

/*        In the real version, B is passed to SLASDQ and multiplied */
/*        internally by Q**H. Here B is complex and that product is */
/*        computed below in two steps (real and imaginary parts). */

#line 341 "clalsd.f"
	j = irwb - 1;
#line 342 "clalsd.f"
	i__1 = *nrhs;
#line 342 "clalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 343 "clalsd.f"
	    i__2 = *n;
#line 343 "clalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 344 "clalsd.f"
		++j;
#line 345 "clalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 345 "clalsd.f"
		rwork[j] = b[i__3].r;
#line 346 "clalsd.f"
/* L40: */
#line 346 "clalsd.f"
	    }
#line 347 "clalsd.f"
/* L50: */
#line 347 "clalsd.f"
	}
#line 348 "clalsd.f"
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n,
		 &c_b35, &rwork[irwrb], n, (ftnlen)1, (ftnlen)1);
#line 350 "clalsd.f"
	j = irwb - 1;
#line 351 "clalsd.f"
	i__1 = *nrhs;
#line 351 "clalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 352 "clalsd.f"
	    i__2 = *n;
#line 352 "clalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 353 "clalsd.f"
		++j;
#line 354 "clalsd.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 355 "clalsd.f"
/* L60: */
#line 355 "clalsd.f"
	    }
#line 356 "clalsd.f"
/* L70: */
#line 356 "clalsd.f"
	}
#line 357 "clalsd.f"
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n,
		 &c_b35, &rwork[irwib], n, (ftnlen)1, (ftnlen)1);
#line 359 "clalsd.f"
	jreal = irwrb - 1;
#line 360 "clalsd.f"
	jimag = irwib - 1;
#line 361 "clalsd.f"
	i__1 = *nrhs;
#line 361 "clalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 362 "clalsd.f"
	    i__2 = *n;
#line 362 "clalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 363 "clalsd.f"
		++jreal;
#line 364 "clalsd.f"
		++jimag;
#line 365 "clalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 365 "clalsd.f"
		i__4 = jreal;
#line 365 "clalsd.f"
		i__5 = jimag;
#line 365 "clalsd.f"
		z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 365 "clalsd.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 366 "clalsd.f"
/* L80: */
#line 366 "clalsd.f"
	    }
#line 367 "clalsd.f"
/* L90: */
#line 367 "clalsd.f"
	}

#line 369 "clalsd.f"
	tol = rcnd * (d__1 = d__[isamax_(n, &d__[1], &c__1)], abs(d__1));
#line 370 "clalsd.f"
	i__1 = *n;
#line 370 "clalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 371 "clalsd.f"
	    if (d__[i__] <= tol) {
#line 372 "clalsd.f"
		claset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 373 "clalsd.f"
	    } else {
#line 374 "clalsd.f"
		clascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &b[
			i__ + b_dim1], ldb, info, (ftnlen)1);
#line 376 "clalsd.f"
		++(*rank);
#line 377 "clalsd.f"
	    }
#line 378 "clalsd.f"
/* L100: */
#line 378 "clalsd.f"
	}

/*        Since B is complex, the following call to SGEMM is performed */
/*        in two steps (real and imaginary parts). That is for V * B */
/*        (in the real version of the code V**H is stored in WORK). */

/*        CALL SGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, */
/*    $               WORK( NWORK ), N ) */

#line 387 "clalsd.f"
	j = irwb - 1;
#line 388 "clalsd.f"
	i__1 = *nrhs;
#line 388 "clalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 389 "clalsd.f"
	    i__2 = *n;
#line 389 "clalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 390 "clalsd.f"
		++j;
#line 391 "clalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 391 "clalsd.f"
		rwork[j] = b[i__3].r;
#line 392 "clalsd.f"
/* L110: */
#line 392 "clalsd.f"
	    }
#line 393 "clalsd.f"
/* L120: */
#line 393 "clalsd.f"
	}
#line 394 "clalsd.f"
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwrb], n, (ftnlen)1, (ftnlen)1);
#line 396 "clalsd.f"
	j = irwb - 1;
#line 397 "clalsd.f"
	i__1 = *nrhs;
#line 397 "clalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 398 "clalsd.f"
	    i__2 = *n;
#line 398 "clalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 399 "clalsd.f"
		++j;
#line 400 "clalsd.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 401 "clalsd.f"
/* L130: */
#line 401 "clalsd.f"
	    }
#line 402 "clalsd.f"
/* L140: */
#line 402 "clalsd.f"
	}
#line 403 "clalsd.f"
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwib], n, (ftnlen)1, (ftnlen)1);
#line 405 "clalsd.f"
	jreal = irwrb - 1;
#line 406 "clalsd.f"
	jimag = irwib - 1;
#line 407 "clalsd.f"
	i__1 = *nrhs;
#line 407 "clalsd.f"
	for (jcol = 1; jcol <= i__1; ++jcol) {
#line 408 "clalsd.f"
	    i__2 = *n;
#line 408 "clalsd.f"
	    for (jrow = 1; jrow <= i__2; ++jrow) {
#line 409 "clalsd.f"
		++jreal;
#line 410 "clalsd.f"
		++jimag;
#line 411 "clalsd.f"
		i__3 = jrow + jcol * b_dim1;
#line 411 "clalsd.f"
		i__4 = jreal;
#line 411 "clalsd.f"
		i__5 = jimag;
#line 411 "clalsd.f"
		z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 411 "clalsd.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 412 "clalsd.f"
/* L150: */
#line 412 "clalsd.f"
	    }
#line 413 "clalsd.f"
/* L160: */
#line 413 "clalsd.f"
	}

/*        Unscale. */

#line 417 "clalsd.f"
	slascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, 
		info, (ftnlen)1);
#line 418 "clalsd.f"
	slasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 419 "clalsd.f"
	clascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);

#line 421 "clalsd.f"
	return 0;
#line 422 "clalsd.f"
    }

/*     Book-keeping and setting up some constants. */

#line 426 "clalsd.f"
    nlvl = (integer) (log((doublereal) (*n) / (doublereal) (*smlsiz + 1)) / 
	    log(2.)) + 1;

#line 428 "clalsd.f"
    smlszp = *smlsiz + 1;

#line 430 "clalsd.f"
    u = 1;
#line 431 "clalsd.f"
    vt = *smlsiz * *n + 1;
#line 432 "clalsd.f"
    difl = vt + smlszp * *n;
#line 433 "clalsd.f"
    difr = difl + nlvl * *n;
#line 434 "clalsd.f"
    z__ = difr + (nlvl * *n << 1);
#line 435 "clalsd.f"
    c__ = z__ + nlvl * *n;
#line 436 "clalsd.f"
    s = c__ + *n;
#line 437 "clalsd.f"
    poles = s + *n;
#line 438 "clalsd.f"
    givnum = poles + (nlvl << 1) * *n;
#line 439 "clalsd.f"
    nrwork = givnum + (nlvl << 1) * *n;
#line 440 "clalsd.f"
    bx = 1;

#line 442 "clalsd.f"
    irwrb = nrwork;
#line 443 "clalsd.f"
    irwib = irwrb + *smlsiz * *nrhs;
#line 444 "clalsd.f"
    irwb = irwib + *smlsiz * *nrhs;

#line 446 "clalsd.f"
    sizei = *n + 1;
#line 447 "clalsd.f"
    k = sizei + *n;
#line 448 "clalsd.f"
    givptr = k + *n;
#line 449 "clalsd.f"
    perm = givptr + *n;
#line 450 "clalsd.f"
    givcol = perm + nlvl * *n;
#line 451 "clalsd.f"
    iwk = givcol + (nlvl * *n << 1);

#line 453 "clalsd.f"
    st = 1;
#line 454 "clalsd.f"
    sqre = 0;
#line 455 "clalsd.f"
    icmpq1 = 1;
#line 456 "clalsd.f"
    icmpq2 = 0;
#line 457 "clalsd.f"
    nsub = 0;

#line 459 "clalsd.f"
    i__1 = *n;
#line 459 "clalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 460 "clalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) < eps) {
#line 461 "clalsd.f"
	    d__[i__] = d_sign(&eps, &d__[i__]);
#line 462 "clalsd.f"
	}
#line 463 "clalsd.f"
/* L170: */
#line 463 "clalsd.f"
    }

#line 465 "clalsd.f"
    i__1 = nm1;
#line 465 "clalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 466 "clalsd.f"
	if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
#line 467 "clalsd.f"
	    ++nsub;
#line 468 "clalsd.f"
	    iwork[nsub] = st;

/*           Subproblem found. First determine its size and then */
/*           apply divide and conquer on it. */

#line 473 "clalsd.f"
	    if (i__ < nm1) {

/*              A subproblem with E(I) small for I < NM1. */

#line 477 "clalsd.f"
		nsize = i__ - st + 1;
#line 478 "clalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 479 "clalsd.f"
	    } else if ((d__1 = e[i__], abs(d__1)) >= eps) {

/*              A subproblem with E(NM1) not too small but I = NM1. */

#line 483 "clalsd.f"
		nsize = *n - st + 1;
#line 484 "clalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 485 "clalsd.f"
	    } else {

/*              A subproblem with E(NM1) small. This implies an */
/*              1-by-1 subproblem at D(N), which is not solved */
/*              explicitly. */

#line 491 "clalsd.f"
		nsize = i__ - st + 1;
#line 492 "clalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 493 "clalsd.f"
		++nsub;
#line 494 "clalsd.f"
		iwork[nsub] = *n;
#line 495 "clalsd.f"
		iwork[sizei + nsub - 1] = 1;
#line 496 "clalsd.f"
		ccopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
#line 497 "clalsd.f"
	    }
#line 498 "clalsd.f"
	    st1 = st - 1;
#line 499 "clalsd.f"
	    if (nsize == 1) {

/*              This is a 1-by-1 subproblem and is not solved */
/*              explicitly. */

#line 504 "clalsd.f"
		ccopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
#line 505 "clalsd.f"
	    } else if (nsize <= *smlsiz) {

/*              This is a small subproblem and is solved by SLASDQ. */

#line 509 "clalsd.f"
		slaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[vt + st1],
			 n, (ftnlen)1);
#line 511 "clalsd.f"
		slaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[u + st1], 
			n, (ftnlen)1);
#line 513 "clalsd.f"
		slasdq_("U", &c__0, &nsize, &nsize, &nsize, &c__0, &d__[st], &
			e[st], &rwork[vt + st1], n, &rwork[u + st1], n, &
			rwork[nrwork], &c__1, &rwork[nrwork], info, (ftnlen)1)
			;
#line 517 "clalsd.f"
		if (*info != 0) {
#line 518 "clalsd.f"
		    return 0;
#line 519 "clalsd.f"
		}

/*              In the real version, B is passed to SLASDQ and multiplied */
/*              internally by Q**H. Here B is complex and that product is */
/*              computed below in two steps (real and imaginary parts). */

#line 525 "clalsd.f"
		j = irwb - 1;
#line 526 "clalsd.f"
		i__2 = *nrhs;
#line 526 "clalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 527 "clalsd.f"
		    i__3 = st + nsize - 1;
#line 527 "clalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 528 "clalsd.f"
			++j;
#line 529 "clalsd.f"
			i__4 = jrow + jcol * b_dim1;
#line 529 "clalsd.f"
			rwork[j] = b[i__4].r;
#line 530 "clalsd.f"
/* L180: */
#line 530 "clalsd.f"
		    }
#line 531 "clalsd.f"
/* L190: */
#line 531 "clalsd.f"
		}
#line 532 "clalsd.f"
		sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
			, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &
			nsize, (ftnlen)1, (ftnlen)1);
#line 535 "clalsd.f"
		j = irwb - 1;
#line 536 "clalsd.f"
		i__2 = *nrhs;
#line 536 "clalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 537 "clalsd.f"
		    i__3 = st + nsize - 1;
#line 537 "clalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 538 "clalsd.f"
			++j;
#line 539 "clalsd.f"
			rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 540 "clalsd.f"
/* L200: */
#line 540 "clalsd.f"
		    }
#line 541 "clalsd.f"
/* L210: */
#line 541 "clalsd.f"
		}
#line 542 "clalsd.f"
		sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
			, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &
			nsize, (ftnlen)1, (ftnlen)1);
#line 545 "clalsd.f"
		jreal = irwrb - 1;
#line 546 "clalsd.f"
		jimag = irwib - 1;
#line 547 "clalsd.f"
		i__2 = *nrhs;
#line 547 "clalsd.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 548 "clalsd.f"
		    i__3 = st + nsize - 1;
#line 548 "clalsd.f"
		    for (jrow = st; jrow <= i__3; ++jrow) {
#line 549 "clalsd.f"
			++jreal;
#line 550 "clalsd.f"
			++jimag;
#line 551 "clalsd.f"
			i__4 = jrow + jcol * b_dim1;
#line 551 "clalsd.f"
			i__5 = jreal;
#line 551 "clalsd.f"
			i__6 = jimag;
#line 551 "clalsd.f"
			z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 551 "clalsd.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 553 "clalsd.f"
/* L220: */
#line 553 "clalsd.f"
		    }
#line 554 "clalsd.f"
/* L230: */
#line 554 "clalsd.f"
		}

#line 556 "clalsd.f"
		clacpy_("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + 
			st1], n, (ftnlen)1);
#line 558 "clalsd.f"
	    } else {

/*              A large problem. Solve it using divide and conquer. */

#line 562 "clalsd.f"
		slasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &
			rwork[u + st1], n, &rwork[vt + st1], &iwork[k + st1], 
			&rwork[difl + st1], &rwork[difr + st1], &rwork[z__ + 
			st1], &rwork[poles + st1], &iwork[givptr + st1], &
			iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
			givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &
			rwork[nrwork], &iwork[iwk], info);
#line 571 "clalsd.f"
		if (*info != 0) {
#line 572 "clalsd.f"
		    return 0;
#line 573 "clalsd.f"
		}
#line 574 "clalsd.f"
		bxst = bx + st1;
#line 575 "clalsd.f"
		clalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &
			work[bxst], n, &rwork[u + st1], n, &rwork[vt + st1], &
			iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1]
			, &rwork[z__ + st1], &rwork[poles + st1], &iwork[
			givptr + st1], &iwork[givcol + st1], n, &iwork[perm + 
			st1], &rwork[givnum + st1], &rwork[c__ + st1], &rwork[
			s + st1], &rwork[nrwork], &iwork[iwk], info);
#line 584 "clalsd.f"
		if (*info != 0) {
#line 585 "clalsd.f"
		    return 0;
#line 586 "clalsd.f"
		}
#line 587 "clalsd.f"
	    }
#line 588 "clalsd.f"
	    st = i__ + 1;
#line 589 "clalsd.f"
	}
#line 590 "clalsd.f"
/* L240: */
#line 590 "clalsd.f"
    }

/*     Apply the singular values and treat the tiny ones as zero. */

#line 594 "clalsd.f"
    tol = rcnd * (d__1 = d__[isamax_(n, &d__[1], &c__1)], abs(d__1));

#line 596 "clalsd.f"
    i__1 = *n;
#line 596 "clalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Some of the elements in D can be negative because 1-by-1 */
/*        subproblems were not solved explicitly. */

#line 601 "clalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) <= tol) {
#line 602 "clalsd.f"
	    claset_("A", &c__1, nrhs, &c_b1, &c_b1, &work[bx + i__ - 1], n, (
		    ftnlen)1);
#line 603 "clalsd.f"
	} else {
#line 604 "clalsd.f"
	    ++(*rank);
#line 605 "clalsd.f"
	    clascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &work[
		    bx + i__ - 1], n, info, (ftnlen)1);
#line 607 "clalsd.f"
	}
#line 608 "clalsd.f"
	d__[i__] = (d__1 = d__[i__], abs(d__1));
#line 609 "clalsd.f"
/* L250: */
#line 609 "clalsd.f"
    }

/*     Now apply back the right singular vectors. */

#line 613 "clalsd.f"
    icmpq2 = 1;
#line 614 "clalsd.f"
    i__1 = nsub;
#line 614 "clalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 615 "clalsd.f"
	st = iwork[i__];
#line 616 "clalsd.f"
	st1 = st - 1;
#line 617 "clalsd.f"
	nsize = iwork[sizei + i__ - 1];
#line 618 "clalsd.f"
	bxst = bx + st1;
#line 619 "clalsd.f"
	if (nsize == 1) {
#line 620 "clalsd.f"
	    ccopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
#line 621 "clalsd.f"
	} else if (nsize <= *smlsiz) {

/*           Since B and BX are complex, the following call to SGEMM */
/*           is performed in two steps (real and imaginary parts). */

/*           CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE, */
/*    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO, */
/*    $                  B( ST, 1 ), LDB ) */

#line 630 "clalsd.f"
	    j = bxst - *n - 1;
#line 631 "clalsd.f"
	    jreal = irwb - 1;
#line 632 "clalsd.f"
	    i__2 = *nrhs;
#line 632 "clalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 633 "clalsd.f"
		j += *n;
#line 634 "clalsd.f"
		i__3 = nsize;
#line 634 "clalsd.f"
		for (jrow = 1; jrow <= i__3; ++jrow) {
#line 635 "clalsd.f"
		    ++jreal;
#line 636 "clalsd.f"
		    i__4 = j + jrow;
#line 636 "clalsd.f"
		    rwork[jreal] = work[i__4].r;
#line 637 "clalsd.f"
/* L260: */
#line 637 "clalsd.f"
		}
#line 638 "clalsd.f"
/* L270: */
#line 638 "clalsd.f"
	    }
#line 639 "clalsd.f"
	    sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &nsize, (
		    ftnlen)1, (ftnlen)1);
#line 642 "clalsd.f"
	    j = bxst - *n - 1;
#line 643 "clalsd.f"
	    jimag = irwb - 1;
#line 644 "clalsd.f"
	    i__2 = *nrhs;
#line 644 "clalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 645 "clalsd.f"
		j += *n;
#line 646 "clalsd.f"
		i__3 = nsize;
#line 646 "clalsd.f"
		for (jrow = 1; jrow <= i__3; ++jrow) {
#line 647 "clalsd.f"
		    ++jimag;
#line 648 "clalsd.f"
		    rwork[jimag] = d_imag(&work[j + jrow]);
#line 649 "clalsd.f"
/* L280: */
#line 649 "clalsd.f"
		}
#line 650 "clalsd.f"
/* L290: */
#line 650 "clalsd.f"
	    }
#line 651 "clalsd.f"
	    sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &nsize, (
		    ftnlen)1, (ftnlen)1);
#line 654 "clalsd.f"
	    jreal = irwrb - 1;
#line 655 "clalsd.f"
	    jimag = irwib - 1;
#line 656 "clalsd.f"
	    i__2 = *nrhs;
#line 656 "clalsd.f"
	    for (jcol = 1; jcol <= i__2; ++jcol) {
#line 657 "clalsd.f"
		i__3 = st + nsize - 1;
#line 657 "clalsd.f"
		for (jrow = st; jrow <= i__3; ++jrow) {
#line 658 "clalsd.f"
		    ++jreal;
#line 659 "clalsd.f"
		    ++jimag;
#line 660 "clalsd.f"
		    i__4 = jrow + jcol * b_dim1;
#line 660 "clalsd.f"
		    i__5 = jreal;
#line 660 "clalsd.f"
		    i__6 = jimag;
#line 660 "clalsd.f"
		    z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 660 "clalsd.f"
		    b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 662 "clalsd.f"
/* L300: */
#line 662 "clalsd.f"
		}
#line 663 "clalsd.f"
/* L310: */
#line 663 "clalsd.f"
	    }
#line 664 "clalsd.f"
	} else {
#line 665 "clalsd.f"
	    clalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + 
		    b_dim1], ldb, &rwork[u + st1], n, &rwork[vt + st1], &
		    iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1], &
		    rwork[z__ + st1], &rwork[poles + st1], &iwork[givptr + 
		    st1], &iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
		    givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &rwork[
		    nrwork], &iwork[iwk], info);
#line 674 "clalsd.f"
	    if (*info != 0) {
#line 675 "clalsd.f"
		return 0;
#line 676 "clalsd.f"
	    }
#line 677 "clalsd.f"
	}
#line 678 "clalsd.f"
/* L320: */
#line 678 "clalsd.f"
    }

/*     Unscale and sort the singular values. */

#line 682 "clalsd.f"
    slascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 683 "clalsd.f"
    slasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 684 "clalsd.f"
    clascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], ldb, 
	    info, (ftnlen)1);

#line 686 "clalsd.f"
    return 0;

/*     End of CLALSD */

} /* clalsd_ */

