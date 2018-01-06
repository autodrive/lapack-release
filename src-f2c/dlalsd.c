#line 1 "dlalsd.f"
/* dlalsd.f -- translated by f2c (version 20100827).
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

#line 1 "dlalsd.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b6 = 0.;
static integer c__0 = 0;
static doublereal c_b11 = 1.;

/* > \brief \b DLALSD uses the singular value decomposition of A to solve the least squares problem. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLALSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlalsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlalsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlalsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, */
/*                          RANK, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLALSD uses the singular value decomposition of A to solve the least */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* >          WORK is DOUBLE PRECISION array, dimension at least */
/* >         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2), */
/* >         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension at least */
/* >         (3*N*NLVL + 11*N) */
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

/* > \ingroup doubleOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb, 
	doublereal *rcond, integer *rank, doublereal *work, integer *iwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), d_sign(doublereal *, doublereal *);

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
    static integer perm, nsub;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nlvl, sqre, bxst;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dcopy_(integer *, doublereal *, integer *, doublereal *, integer 
	    *);
    static integer poles, sizei, nsize, nwork, icmpq1, icmpq2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlasda_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), dlalsa_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlasdq_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dlaset_(char *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer givcol;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal orgnrm;
    static integer givnum, givptr, smlszp;


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

#line 226 "dlalsd.f"
    /* Parameter adjustments */
#line 226 "dlalsd.f"
    --d__;
#line 226 "dlalsd.f"
    --e;
#line 226 "dlalsd.f"
    b_dim1 = *ldb;
#line 226 "dlalsd.f"
    b_offset = 1 + b_dim1;
#line 226 "dlalsd.f"
    b -= b_offset;
#line 226 "dlalsd.f"
    --work;
#line 226 "dlalsd.f"
    --iwork;
#line 226 "dlalsd.f"

#line 226 "dlalsd.f"
    /* Function Body */
#line 226 "dlalsd.f"
    *info = 0;

#line 228 "dlalsd.f"
    if (*n < 0) {
#line 229 "dlalsd.f"
	*info = -3;
#line 230 "dlalsd.f"
    } else if (*nrhs < 1) {
#line 231 "dlalsd.f"
	*info = -4;
#line 232 "dlalsd.f"
    } else if (*ldb < 1 || *ldb < *n) {
#line 233 "dlalsd.f"
	*info = -8;
#line 234 "dlalsd.f"
    }
#line 235 "dlalsd.f"
    if (*info != 0) {
#line 236 "dlalsd.f"
	i__1 = -(*info);
#line 236 "dlalsd.f"
	xerbla_("DLALSD", &i__1, (ftnlen)6);
#line 237 "dlalsd.f"
	return 0;
#line 238 "dlalsd.f"
    }

#line 240 "dlalsd.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Set up the tolerance. */

#line 244 "dlalsd.f"
    if (*rcond <= 0. || *rcond >= 1.) {
#line 245 "dlalsd.f"
	rcnd = eps;
#line 246 "dlalsd.f"
    } else {
#line 247 "dlalsd.f"
	rcnd = *rcond;
#line 248 "dlalsd.f"
    }

#line 250 "dlalsd.f"
    *rank = 0;

/*     Quick return if possible. */

#line 254 "dlalsd.f"
    if (*n == 0) {
#line 255 "dlalsd.f"
	return 0;
#line 256 "dlalsd.f"
    } else if (*n == 1) {
#line 257 "dlalsd.f"
	if (d__[1] == 0.) {
#line 258 "dlalsd.f"
	    dlaset_("A", &c__1, nrhs, &c_b6, &c_b6, &b[b_offset], ldb, (
		    ftnlen)1);
#line 259 "dlalsd.f"
	} else {
#line 260 "dlalsd.f"
	    *rank = 1;
#line 261 "dlalsd.f"
	    dlascl_("G", &c__0, &c__0, &d__[1], &c_b11, &c__1, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)1);
#line 262 "dlalsd.f"
	    d__[1] = abs(d__[1]);
#line 263 "dlalsd.f"
	}
#line 264 "dlalsd.f"
	return 0;
#line 265 "dlalsd.f"
    }

/*     Rotate the matrix if it is lower bidiagonal. */

#line 269 "dlalsd.f"
    if (*(unsigned char *)uplo == 'L') {
#line 270 "dlalsd.f"
	i__1 = *n - 1;
#line 270 "dlalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "dlalsd.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 272 "dlalsd.f"
	    d__[i__] = r__;
#line 273 "dlalsd.f"
	    e[i__] = sn * d__[i__ + 1];
#line 274 "dlalsd.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 275 "dlalsd.f"
	    if (*nrhs == 1) {
#line 276 "dlalsd.f"
		drot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &
			c__1, &cs, &sn);
#line 277 "dlalsd.f"
	    } else {
#line 278 "dlalsd.f"
		work[(i__ << 1) - 1] = cs;
#line 279 "dlalsd.f"
		work[i__ * 2] = sn;
#line 280 "dlalsd.f"
	    }
#line 281 "dlalsd.f"
/* L10: */
#line 281 "dlalsd.f"
	}
#line 282 "dlalsd.f"
	if (*nrhs > 1) {
#line 283 "dlalsd.f"
	    i__1 = *nrhs;
#line 283 "dlalsd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 284 "dlalsd.f"
		i__2 = *n - 1;
#line 284 "dlalsd.f"
		for (j = 1; j <= i__2; ++j) {
#line 285 "dlalsd.f"
		    cs = work[(j << 1) - 1];
#line 286 "dlalsd.f"
		    sn = work[j * 2];
#line 287 "dlalsd.f"
		    drot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ *
			     b_dim1], &c__1, &cs, &sn);
#line 288 "dlalsd.f"
/* L20: */
#line 288 "dlalsd.f"
		}
#line 289 "dlalsd.f"
/* L30: */
#line 289 "dlalsd.f"
	    }
#line 290 "dlalsd.f"
	}
#line 291 "dlalsd.f"
    }

/*     Scale. */

#line 295 "dlalsd.f"
    nm1 = *n - 1;
#line 296 "dlalsd.f"
    orgnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 297 "dlalsd.f"
    if (orgnrm == 0.) {
#line 298 "dlalsd.f"
	dlaset_("A", n, nrhs, &c_b6, &c_b6, &b[b_offset], ldb, (ftnlen)1);
#line 299 "dlalsd.f"
	return 0;
#line 300 "dlalsd.f"
    }

#line 302 "dlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b11, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 303 "dlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b11, &nm1, &c__1, &e[1], &nm1, 
	    info, (ftnlen)1);

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

#line 308 "dlalsd.f"
    if (*n <= *smlsiz) {
#line 309 "dlalsd.f"
	nwork = *n * *n + 1;
#line 310 "dlalsd.f"
	dlaset_("A", n, n, &c_b6, &c_b11, &work[1], n, (ftnlen)1);
#line 311 "dlalsd.f"
	dlasdq_("U", &c__0, n, n, &c__0, nrhs, &d__[1], &e[1], &work[1], n, &
		work[1], n, &b[b_offset], ldb, &work[nwork], info, (ftnlen)1);
#line 313 "dlalsd.f"
	if (*info != 0) {
#line 314 "dlalsd.f"
	    return 0;
#line 315 "dlalsd.f"
	}
#line 316 "dlalsd.f"
	tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));
#line 317 "dlalsd.f"
	i__1 = *n;
#line 317 "dlalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 318 "dlalsd.f"
	    if (d__[i__] <= tol) {
#line 319 "dlalsd.f"
		dlaset_("A", &c__1, nrhs, &c_b6, &c_b6, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 320 "dlalsd.f"
	    } else {
#line 321 "dlalsd.f"
		dlascl_("G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs, &b[
			i__ + b_dim1], ldb, info, (ftnlen)1);
#line 323 "dlalsd.f"
		++(*rank);
#line 324 "dlalsd.f"
	    }
#line 325 "dlalsd.f"
/* L40: */
#line 325 "dlalsd.f"
	}
#line 326 "dlalsd.f"
	dgemm_("T", "N", n, nrhs, n, &c_b11, &work[1], n, &b[b_offset], ldb, &
		c_b6, &work[nwork], n, (ftnlen)1, (ftnlen)1);
#line 328 "dlalsd.f"
	dlacpy_("A", n, nrhs, &work[nwork], n, &b[b_offset], ldb, (ftnlen)1);

/*        Unscale. */

#line 332 "dlalsd.f"
	dlascl_("G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, 
		info, (ftnlen)1);
#line 333 "dlalsd.f"
	dlasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 334 "dlalsd.f"
	dlascl_("G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);

#line 336 "dlalsd.f"
	return 0;
#line 337 "dlalsd.f"
    }

/*     Book-keeping and setting up some constants. */

#line 341 "dlalsd.f"
    nlvl = (integer) (log((doublereal) (*n) / (doublereal) (*smlsiz + 1)) / 
	    log(2.)) + 1;

#line 343 "dlalsd.f"
    smlszp = *smlsiz + 1;

#line 345 "dlalsd.f"
    u = 1;
#line 346 "dlalsd.f"
    vt = *smlsiz * *n + 1;
#line 347 "dlalsd.f"
    difl = vt + smlszp * *n;
#line 348 "dlalsd.f"
    difr = difl + nlvl * *n;
#line 349 "dlalsd.f"
    z__ = difr + (nlvl * *n << 1);
#line 350 "dlalsd.f"
    c__ = z__ + nlvl * *n;
#line 351 "dlalsd.f"
    s = c__ + *n;
#line 352 "dlalsd.f"
    poles = s + *n;
#line 353 "dlalsd.f"
    givnum = poles + (nlvl << 1) * *n;
#line 354 "dlalsd.f"
    bx = givnum + (nlvl << 1) * *n;
#line 355 "dlalsd.f"
    nwork = bx + *n * *nrhs;

#line 357 "dlalsd.f"
    sizei = *n + 1;
#line 358 "dlalsd.f"
    k = sizei + *n;
#line 359 "dlalsd.f"
    givptr = k + *n;
#line 360 "dlalsd.f"
    perm = givptr + *n;
#line 361 "dlalsd.f"
    givcol = perm + nlvl * *n;
#line 362 "dlalsd.f"
    iwk = givcol + (nlvl * *n << 1);

#line 364 "dlalsd.f"
    st = 1;
#line 365 "dlalsd.f"
    sqre = 0;
#line 366 "dlalsd.f"
    icmpq1 = 1;
#line 367 "dlalsd.f"
    icmpq2 = 0;
#line 368 "dlalsd.f"
    nsub = 0;

#line 370 "dlalsd.f"
    i__1 = *n;
#line 370 "dlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 371 "dlalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) < eps) {
#line 372 "dlalsd.f"
	    d__[i__] = d_sign(&eps, &d__[i__]);
#line 373 "dlalsd.f"
	}
#line 374 "dlalsd.f"
/* L50: */
#line 374 "dlalsd.f"
    }

#line 376 "dlalsd.f"
    i__1 = nm1;
#line 376 "dlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "dlalsd.f"
	if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
#line 378 "dlalsd.f"
	    ++nsub;
#line 379 "dlalsd.f"
	    iwork[nsub] = st;

/*           Subproblem found. First determine its size and then */
/*           apply divide and conquer on it. */

#line 384 "dlalsd.f"
	    if (i__ < nm1) {

/*              A subproblem with E(I) small for I < NM1. */

#line 388 "dlalsd.f"
		nsize = i__ - st + 1;
#line 389 "dlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 390 "dlalsd.f"
	    } else if ((d__1 = e[i__], abs(d__1)) >= eps) {

/*              A subproblem with E(NM1) not too small but I = NM1. */

#line 394 "dlalsd.f"
		nsize = *n - st + 1;
#line 395 "dlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 396 "dlalsd.f"
	    } else {

/*              A subproblem with E(NM1) small. This implies an */
/*              1-by-1 subproblem at D(N), which is not solved */
/*              explicitly. */

#line 402 "dlalsd.f"
		nsize = i__ - st + 1;
#line 403 "dlalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 404 "dlalsd.f"
		++nsub;
#line 405 "dlalsd.f"
		iwork[nsub] = *n;
#line 406 "dlalsd.f"
		iwork[sizei + nsub - 1] = 1;
#line 407 "dlalsd.f"
		dcopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
#line 408 "dlalsd.f"
	    }
#line 409 "dlalsd.f"
	    st1 = st - 1;
#line 410 "dlalsd.f"
	    if (nsize == 1) {

/*              This is a 1-by-1 subproblem and is not solved */
/*              explicitly. */

#line 415 "dlalsd.f"
		dcopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
#line 416 "dlalsd.f"
	    } else if (nsize <= *smlsiz) {

/*              This is a small subproblem and is solved by DLASDQ. */

#line 420 "dlalsd.f"
		dlaset_("A", &nsize, &nsize, &c_b6, &c_b11, &work[vt + st1], 
			n, (ftnlen)1);
#line 422 "dlalsd.f"
		dlasdq_("U", &c__0, &nsize, &nsize, &c__0, nrhs, &d__[st], &e[
			st], &work[vt + st1], n, &work[nwork], n, &b[st + 
			b_dim1], ldb, &work[nwork], info, (ftnlen)1);
#line 425 "dlalsd.f"
		if (*info != 0) {
#line 426 "dlalsd.f"
		    return 0;
#line 427 "dlalsd.f"
		}
#line 428 "dlalsd.f"
		dlacpy_("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + 
			st1], n, (ftnlen)1);
#line 430 "dlalsd.f"
	    } else {

/*              A large problem. Solve it using divide and conquer. */

#line 434 "dlalsd.f"
		dlasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &
			work[u + st1], n, &work[vt + st1], &iwork[k + st1], &
			work[difl + st1], &work[difr + st1], &work[z__ + st1],
			 &work[poles + st1], &iwork[givptr + st1], &iwork[
			givcol + st1], n, &iwork[perm + st1], &work[givnum + 
			st1], &work[c__ + st1], &work[s + st1], &work[nwork], 
			&iwork[iwk], info);
#line 443 "dlalsd.f"
		if (*info != 0) {
#line 444 "dlalsd.f"
		    return 0;
#line 445 "dlalsd.f"
		}
#line 446 "dlalsd.f"
		bxst = bx + st1;
#line 447 "dlalsd.f"
		dlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &
			work[bxst], n, &work[u + st1], n, &work[vt + st1], &
			iwork[k + st1], &work[difl + st1], &work[difr + st1], 
			&work[z__ + st1], &work[poles + st1], &iwork[givptr + 
			st1], &iwork[givcol + st1], n, &iwork[perm + st1], &
			work[givnum + st1], &work[c__ + st1], &work[s + st1], 
			&work[nwork], &iwork[iwk], info);
#line 456 "dlalsd.f"
		if (*info != 0) {
#line 457 "dlalsd.f"
		    return 0;
#line 458 "dlalsd.f"
		}
#line 459 "dlalsd.f"
	    }
#line 460 "dlalsd.f"
	    st = i__ + 1;
#line 461 "dlalsd.f"
	}
#line 462 "dlalsd.f"
/* L60: */
#line 462 "dlalsd.f"
    }

/*     Apply the singular values and treat the tiny ones as zero. */

#line 466 "dlalsd.f"
    tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));

#line 468 "dlalsd.f"
    i__1 = *n;
#line 468 "dlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Some of the elements in D can be negative because 1-by-1 */
/*        subproblems were not solved explicitly. */

#line 473 "dlalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) <= tol) {
#line 474 "dlalsd.f"
	    dlaset_("A", &c__1, nrhs, &c_b6, &c_b6, &work[bx + i__ - 1], n, (
		    ftnlen)1);
#line 475 "dlalsd.f"
	} else {
#line 476 "dlalsd.f"
	    ++(*rank);
#line 477 "dlalsd.f"
	    dlascl_("G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs, &work[
		    bx + i__ - 1], n, info, (ftnlen)1);
#line 479 "dlalsd.f"
	}
#line 480 "dlalsd.f"
	d__[i__] = (d__1 = d__[i__], abs(d__1));
#line 481 "dlalsd.f"
/* L70: */
#line 481 "dlalsd.f"
    }

/*     Now apply back the right singular vectors. */

#line 485 "dlalsd.f"
    icmpq2 = 1;
#line 486 "dlalsd.f"
    i__1 = nsub;
#line 486 "dlalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 487 "dlalsd.f"
	st = iwork[i__];
#line 488 "dlalsd.f"
	st1 = st - 1;
#line 489 "dlalsd.f"
	nsize = iwork[sizei + i__ - 1];
#line 490 "dlalsd.f"
	bxst = bx + st1;
#line 491 "dlalsd.f"
	if (nsize == 1) {
#line 492 "dlalsd.f"
	    dcopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
#line 493 "dlalsd.f"
	} else if (nsize <= *smlsiz) {
#line 494 "dlalsd.f"
	    dgemm_("T", "N", &nsize, nrhs, &nsize, &c_b11, &work[vt + st1], n,
		     &work[bxst], n, &c_b6, &b[st + b_dim1], ldb, (ftnlen)1, (
		    ftnlen)1);
#line 497 "dlalsd.f"
	} else {
#line 498 "dlalsd.f"
	    dlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + 
		    b_dim1], ldb, &work[u + st1], n, &work[vt + st1], &iwork[
		    k + st1], &work[difl + st1], &work[difr + st1], &work[z__ 
		    + st1], &work[poles + st1], &iwork[givptr + st1], &iwork[
		    givcol + st1], n, &iwork[perm + st1], &work[givnum + st1],
		     &work[c__ + st1], &work[s + st1], &work[nwork], &iwork[
		    iwk], info);
#line 507 "dlalsd.f"
	    if (*info != 0) {
#line 508 "dlalsd.f"
		return 0;
#line 509 "dlalsd.f"
	    }
#line 510 "dlalsd.f"
	}
#line 511 "dlalsd.f"
/* L80: */
#line 511 "dlalsd.f"
    }

/*     Unscale and sort the singular values. */

#line 515 "dlalsd.f"
    dlascl_("G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 516 "dlalsd.f"
    dlasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 517 "dlalsd.f"
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], ldb, 
	    info, (ftnlen)1);

#line 519 "dlalsd.f"
    return 0;

/*     End of DLALSD */

} /* dlalsd_ */

