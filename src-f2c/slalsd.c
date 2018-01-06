#line 1 "slalsd.f"
/* slalsd.f -- translated by f2c (version 20100827).
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

#line 1 "slalsd.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b6 = 0.;
static integer c__0 = 0;
static doublereal c_b11 = 1.;

/* > \brief \b SLALSD uses the singular value decomposition of A to solve the least squares problem. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLALSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slalsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slalsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slalsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, */
/*                          RANK, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               B( LDB, * ), D( * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALSD uses the singular value decomposition of A to solve the least */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          WORK is REAL array, dimension at least */
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

/* > \ingroup realOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int slalsd_(char *uplo, integer *smlsiz, integer *n, integer 
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
    static integer perm, nsub, nlvl, sqre, bxst;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), sgemm_(char 
	    *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static integer poles, sizei, nsize;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nwork, icmpq1, icmpq2;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int slasda_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), xerbla_(char *, integer *, ftnlen), slalsa_(integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    slascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static integer givcol;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slasdq_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), slacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), slartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), slaset_(char *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal orgnrm;
    static integer givnum;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int slasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static integer givptr, smlszp;


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

#line 226 "slalsd.f"
    /* Parameter adjustments */
#line 226 "slalsd.f"
    --d__;
#line 226 "slalsd.f"
    --e;
#line 226 "slalsd.f"
    b_dim1 = *ldb;
#line 226 "slalsd.f"
    b_offset = 1 + b_dim1;
#line 226 "slalsd.f"
    b -= b_offset;
#line 226 "slalsd.f"
    --work;
#line 226 "slalsd.f"
    --iwork;
#line 226 "slalsd.f"

#line 226 "slalsd.f"
    /* Function Body */
#line 226 "slalsd.f"
    *info = 0;

#line 228 "slalsd.f"
    if (*n < 0) {
#line 229 "slalsd.f"
	*info = -3;
#line 230 "slalsd.f"
    } else if (*nrhs < 1) {
#line 231 "slalsd.f"
	*info = -4;
#line 232 "slalsd.f"
    } else if (*ldb < 1 || *ldb < *n) {
#line 233 "slalsd.f"
	*info = -8;
#line 234 "slalsd.f"
    }
#line 235 "slalsd.f"
    if (*info != 0) {
#line 236 "slalsd.f"
	i__1 = -(*info);
#line 236 "slalsd.f"
	xerbla_("SLALSD", &i__1, (ftnlen)6);
#line 237 "slalsd.f"
	return 0;
#line 238 "slalsd.f"
    }

#line 240 "slalsd.f"
    eps = slamch_("Epsilon", (ftnlen)7);

/*     Set up the tolerance. */

#line 244 "slalsd.f"
    if (*rcond <= 0. || *rcond >= 1.) {
#line 245 "slalsd.f"
	rcnd = eps;
#line 246 "slalsd.f"
    } else {
#line 247 "slalsd.f"
	rcnd = *rcond;
#line 248 "slalsd.f"
    }

#line 250 "slalsd.f"
    *rank = 0;

/*     Quick return if possible. */

#line 254 "slalsd.f"
    if (*n == 0) {
#line 255 "slalsd.f"
	return 0;
#line 256 "slalsd.f"
    } else if (*n == 1) {
#line 257 "slalsd.f"
	if (d__[1] == 0.) {
#line 258 "slalsd.f"
	    slaset_("A", &c__1, nrhs, &c_b6, &c_b6, &b[b_offset], ldb, (
		    ftnlen)1);
#line 259 "slalsd.f"
	} else {
#line 260 "slalsd.f"
	    *rank = 1;
#line 261 "slalsd.f"
	    slascl_("G", &c__0, &c__0, &d__[1], &c_b11, &c__1, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)1);
#line 262 "slalsd.f"
	    d__[1] = abs(d__[1]);
#line 263 "slalsd.f"
	}
#line 264 "slalsd.f"
	return 0;
#line 265 "slalsd.f"
    }

/*     Rotate the matrix if it is lower bidiagonal. */

#line 269 "slalsd.f"
    if (*(unsigned char *)uplo == 'L') {
#line 270 "slalsd.f"
	i__1 = *n - 1;
#line 270 "slalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "slalsd.f"
	    slartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 272 "slalsd.f"
	    d__[i__] = r__;
#line 273 "slalsd.f"
	    e[i__] = sn * d__[i__ + 1];
#line 274 "slalsd.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 275 "slalsd.f"
	    if (*nrhs == 1) {
#line 276 "slalsd.f"
		srot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &
			c__1, &cs, &sn);
#line 277 "slalsd.f"
	    } else {
#line 278 "slalsd.f"
		work[(i__ << 1) - 1] = cs;
#line 279 "slalsd.f"
		work[i__ * 2] = sn;
#line 280 "slalsd.f"
	    }
#line 281 "slalsd.f"
/* L10: */
#line 281 "slalsd.f"
	}
#line 282 "slalsd.f"
	if (*nrhs > 1) {
#line 283 "slalsd.f"
	    i__1 = *nrhs;
#line 283 "slalsd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 284 "slalsd.f"
		i__2 = *n - 1;
#line 284 "slalsd.f"
		for (j = 1; j <= i__2; ++j) {
#line 285 "slalsd.f"
		    cs = work[(j << 1) - 1];
#line 286 "slalsd.f"
		    sn = work[j * 2];
#line 287 "slalsd.f"
		    srot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ *
			     b_dim1], &c__1, &cs, &sn);
#line 288 "slalsd.f"
/* L20: */
#line 288 "slalsd.f"
		}
#line 289 "slalsd.f"
/* L30: */
#line 289 "slalsd.f"
	    }
#line 290 "slalsd.f"
	}
#line 291 "slalsd.f"
    }

/*     Scale. */

#line 295 "slalsd.f"
    nm1 = *n - 1;
#line 296 "slalsd.f"
    orgnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 297 "slalsd.f"
    if (orgnrm == 0.) {
#line 298 "slalsd.f"
	slaset_("A", n, nrhs, &c_b6, &c_b6, &b[b_offset], ldb, (ftnlen)1);
#line 299 "slalsd.f"
	return 0;
#line 300 "slalsd.f"
    }

#line 302 "slalsd.f"
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b11, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 303 "slalsd.f"
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b11, &nm1, &c__1, &e[1], &nm1, 
	    info, (ftnlen)1);

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

#line 308 "slalsd.f"
    if (*n <= *smlsiz) {
#line 309 "slalsd.f"
	nwork = *n * *n + 1;
#line 310 "slalsd.f"
	slaset_("A", n, n, &c_b6, &c_b11, &work[1], n, (ftnlen)1);
#line 311 "slalsd.f"
	slasdq_("U", &c__0, n, n, &c__0, nrhs, &d__[1], &e[1], &work[1], n, &
		work[1], n, &b[b_offset], ldb, &work[nwork], info, (ftnlen)1);
#line 313 "slalsd.f"
	if (*info != 0) {
#line 314 "slalsd.f"
	    return 0;
#line 315 "slalsd.f"
	}
#line 316 "slalsd.f"
	tol = rcnd * (d__1 = d__[isamax_(n, &d__[1], &c__1)], abs(d__1));
#line 317 "slalsd.f"
	i__1 = *n;
#line 317 "slalsd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 318 "slalsd.f"
	    if (d__[i__] <= tol) {
#line 319 "slalsd.f"
		slaset_("A", &c__1, nrhs, &c_b6, &c_b6, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 320 "slalsd.f"
	    } else {
#line 321 "slalsd.f"
		slascl_("G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs, &b[
			i__ + b_dim1], ldb, info, (ftnlen)1);
#line 323 "slalsd.f"
		++(*rank);
#line 324 "slalsd.f"
	    }
#line 325 "slalsd.f"
/* L40: */
#line 325 "slalsd.f"
	}
#line 326 "slalsd.f"
	sgemm_("T", "N", n, nrhs, n, &c_b11, &work[1], n, &b[b_offset], ldb, &
		c_b6, &work[nwork], n, (ftnlen)1, (ftnlen)1);
#line 328 "slalsd.f"
	slacpy_("A", n, nrhs, &work[nwork], n, &b[b_offset], ldb, (ftnlen)1);

/*        Unscale. */

#line 332 "slalsd.f"
	slascl_("G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, 
		info, (ftnlen)1);
#line 333 "slalsd.f"
	slasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 334 "slalsd.f"
	slascl_("G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], 
		ldb, info, (ftnlen)1);

#line 336 "slalsd.f"
	return 0;
#line 337 "slalsd.f"
    }

/*     Book-keeping and setting up some constants. */

#line 341 "slalsd.f"
    nlvl = (integer) (log((doublereal) (*n) / (doublereal) (*smlsiz + 1)) / 
	    log(2.)) + 1;

#line 343 "slalsd.f"
    smlszp = *smlsiz + 1;

#line 345 "slalsd.f"
    u = 1;
#line 346 "slalsd.f"
    vt = *smlsiz * *n + 1;
#line 347 "slalsd.f"
    difl = vt + smlszp * *n;
#line 348 "slalsd.f"
    difr = difl + nlvl * *n;
#line 349 "slalsd.f"
    z__ = difr + (nlvl * *n << 1);
#line 350 "slalsd.f"
    c__ = z__ + nlvl * *n;
#line 351 "slalsd.f"
    s = c__ + *n;
#line 352 "slalsd.f"
    poles = s + *n;
#line 353 "slalsd.f"
    givnum = poles + (nlvl << 1) * *n;
#line 354 "slalsd.f"
    bx = givnum + (nlvl << 1) * *n;
#line 355 "slalsd.f"
    nwork = bx + *n * *nrhs;

#line 357 "slalsd.f"
    sizei = *n + 1;
#line 358 "slalsd.f"
    k = sizei + *n;
#line 359 "slalsd.f"
    givptr = k + *n;
#line 360 "slalsd.f"
    perm = givptr + *n;
#line 361 "slalsd.f"
    givcol = perm + nlvl * *n;
#line 362 "slalsd.f"
    iwk = givcol + (nlvl * *n << 1);

#line 364 "slalsd.f"
    st = 1;
#line 365 "slalsd.f"
    sqre = 0;
#line 366 "slalsd.f"
    icmpq1 = 1;
#line 367 "slalsd.f"
    icmpq2 = 0;
#line 368 "slalsd.f"
    nsub = 0;

#line 370 "slalsd.f"
    i__1 = *n;
#line 370 "slalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 371 "slalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) < eps) {
#line 372 "slalsd.f"
	    d__[i__] = d_sign(&eps, &d__[i__]);
#line 373 "slalsd.f"
	}
#line 374 "slalsd.f"
/* L50: */
#line 374 "slalsd.f"
    }

#line 376 "slalsd.f"
    i__1 = nm1;
#line 376 "slalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "slalsd.f"
	if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
#line 378 "slalsd.f"
	    ++nsub;
#line 379 "slalsd.f"
	    iwork[nsub] = st;

/*           Subproblem found. First determine its size and then */
/*           apply divide and conquer on it. */

#line 384 "slalsd.f"
	    if (i__ < nm1) {

/*              A subproblem with E(I) small for I < NM1. */

#line 388 "slalsd.f"
		nsize = i__ - st + 1;
#line 389 "slalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 390 "slalsd.f"
	    } else if ((d__1 = e[i__], abs(d__1)) >= eps) {

/*              A subproblem with E(NM1) not too small but I = NM1. */

#line 394 "slalsd.f"
		nsize = *n - st + 1;
#line 395 "slalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 396 "slalsd.f"
	    } else {

/*              A subproblem with E(NM1) small. This implies an */
/*              1-by-1 subproblem at D(N), which is not solved */
/*              explicitly. */

#line 402 "slalsd.f"
		nsize = i__ - st + 1;
#line 403 "slalsd.f"
		iwork[sizei + nsub - 1] = nsize;
#line 404 "slalsd.f"
		++nsub;
#line 405 "slalsd.f"
		iwork[nsub] = *n;
#line 406 "slalsd.f"
		iwork[sizei + nsub - 1] = 1;
#line 407 "slalsd.f"
		scopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
#line 408 "slalsd.f"
	    }
#line 409 "slalsd.f"
	    st1 = st - 1;
#line 410 "slalsd.f"
	    if (nsize == 1) {

/*              This is a 1-by-1 subproblem and is not solved */
/*              explicitly. */

#line 415 "slalsd.f"
		scopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
#line 416 "slalsd.f"
	    } else if (nsize <= *smlsiz) {

/*              This is a small subproblem and is solved by SLASDQ. */

#line 420 "slalsd.f"
		slaset_("A", &nsize, &nsize, &c_b6, &c_b11, &work[vt + st1], 
			n, (ftnlen)1);
#line 422 "slalsd.f"
		slasdq_("U", &c__0, &nsize, &nsize, &c__0, nrhs, &d__[st], &e[
			st], &work[vt + st1], n, &work[nwork], n, &b[st + 
			b_dim1], ldb, &work[nwork], info, (ftnlen)1);
#line 425 "slalsd.f"
		if (*info != 0) {
#line 426 "slalsd.f"
		    return 0;
#line 427 "slalsd.f"
		}
#line 428 "slalsd.f"
		slacpy_("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + 
			st1], n, (ftnlen)1);
#line 430 "slalsd.f"
	    } else {

/*              A large problem. Solve it using divide and conquer. */

#line 434 "slalsd.f"
		slasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &
			work[u + st1], n, &work[vt + st1], &iwork[k + st1], &
			work[difl + st1], &work[difr + st1], &work[z__ + st1],
			 &work[poles + st1], &iwork[givptr + st1], &iwork[
			givcol + st1], n, &iwork[perm + st1], &work[givnum + 
			st1], &work[c__ + st1], &work[s + st1], &work[nwork], 
			&iwork[iwk], info);
#line 443 "slalsd.f"
		if (*info != 0) {
#line 444 "slalsd.f"
		    return 0;
#line 445 "slalsd.f"
		}
#line 446 "slalsd.f"
		bxst = bx + st1;
#line 447 "slalsd.f"
		slalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &
			work[bxst], n, &work[u + st1], n, &work[vt + st1], &
			iwork[k + st1], &work[difl + st1], &work[difr + st1], 
			&work[z__ + st1], &work[poles + st1], &iwork[givptr + 
			st1], &iwork[givcol + st1], n, &iwork[perm + st1], &
			work[givnum + st1], &work[c__ + st1], &work[s + st1], 
			&work[nwork], &iwork[iwk], info);
#line 456 "slalsd.f"
		if (*info != 0) {
#line 457 "slalsd.f"
		    return 0;
#line 458 "slalsd.f"
		}
#line 459 "slalsd.f"
	    }
#line 460 "slalsd.f"
	    st = i__ + 1;
#line 461 "slalsd.f"
	}
#line 462 "slalsd.f"
/* L60: */
#line 462 "slalsd.f"
    }

/*     Apply the singular values and treat the tiny ones as zero. */

#line 466 "slalsd.f"
    tol = rcnd * (d__1 = d__[isamax_(n, &d__[1], &c__1)], abs(d__1));

#line 468 "slalsd.f"
    i__1 = *n;
#line 468 "slalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Some of the elements in D can be negative because 1-by-1 */
/*        subproblems were not solved explicitly. */

#line 473 "slalsd.f"
	if ((d__1 = d__[i__], abs(d__1)) <= tol) {
#line 474 "slalsd.f"
	    slaset_("A", &c__1, nrhs, &c_b6, &c_b6, &work[bx + i__ - 1], n, (
		    ftnlen)1);
#line 475 "slalsd.f"
	} else {
#line 476 "slalsd.f"
	    ++(*rank);
#line 477 "slalsd.f"
	    slascl_("G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs, &work[
		    bx + i__ - 1], n, info, (ftnlen)1);
#line 479 "slalsd.f"
	}
#line 480 "slalsd.f"
	d__[i__] = (d__1 = d__[i__], abs(d__1));
#line 481 "slalsd.f"
/* L70: */
#line 481 "slalsd.f"
    }

/*     Now apply back the right singular vectors. */

#line 485 "slalsd.f"
    icmpq2 = 1;
#line 486 "slalsd.f"
    i__1 = nsub;
#line 486 "slalsd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 487 "slalsd.f"
	st = iwork[i__];
#line 488 "slalsd.f"
	st1 = st - 1;
#line 489 "slalsd.f"
	nsize = iwork[sizei + i__ - 1];
#line 490 "slalsd.f"
	bxst = bx + st1;
#line 491 "slalsd.f"
	if (nsize == 1) {
#line 492 "slalsd.f"
	    scopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
#line 493 "slalsd.f"
	} else if (nsize <= *smlsiz) {
#line 494 "slalsd.f"
	    sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b11, &work[vt + st1], n,
		     &work[bxst], n, &c_b6, &b[st + b_dim1], ldb, (ftnlen)1, (
		    ftnlen)1);
#line 497 "slalsd.f"
	} else {
#line 498 "slalsd.f"
	    slalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + 
		    b_dim1], ldb, &work[u + st1], n, &work[vt + st1], &iwork[
		    k + st1], &work[difl + st1], &work[difr + st1], &work[z__ 
		    + st1], &work[poles + st1], &iwork[givptr + st1], &iwork[
		    givcol + st1], n, &iwork[perm + st1], &work[givnum + st1],
		     &work[c__ + st1], &work[s + st1], &work[nwork], &iwork[
		    iwk], info);
#line 507 "slalsd.f"
	    if (*info != 0) {
#line 508 "slalsd.f"
		return 0;
#line 509 "slalsd.f"
	    }
#line 510 "slalsd.f"
	}
#line 511 "slalsd.f"
/* L80: */
#line 511 "slalsd.f"
    }

/*     Unscale and sort the singular values. */

#line 515 "slalsd.f"
    slascl_("G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, info, (
	    ftnlen)1);
#line 516 "slalsd.f"
    slasrt_("D", n, &d__[1], info, (ftnlen)1);
#line 517 "slalsd.f"
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], ldb, 
	    info, (ftnlen)1);

#line 519 "slalsd.f"
    return 0;

/*     End of SLALSD */

} /* slalsd_ */

