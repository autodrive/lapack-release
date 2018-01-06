#line 1 "csteqr.f"
/* csteqr.f -- translated by f2c (version 20100827).
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

#line 1 "csteqr.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b41 = 1.;

/* > \brief \b CSTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), WORK( * ) */
/*       COMPLEX            Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the implicit QL or QR method. */
/* > The eigenvectors of a full or band complex Hermitian matrix can also */
/* > be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this */
/* > matrix to tridiagonal form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only. */
/* >          = 'V':  Compute eigenvalues and eigenvectors of the original */
/* >                  Hermitian matrix.  On entry, Z must contain the */
/* >                  unitary matrix used to reduce the original matrix */
/* >                  to tridiagonal form. */
/* >          = 'I':  Compute eigenvalues and eigenvectors of the */
/* >                  tridiagonal matrix.  Z is initialized to the identity */
/* >                  matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          On entry, if  COMPZ = 'V', then Z contains the unitary */
/* >          matrix used in the reduction to tridiagonal form. */
/* >          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the */
/* >          orthonormal eigenvectors of the original Hermitian matrix, */
/* >          and if COMPZ = 'I', Z contains the orthonormal eigenvectors */
/* >          of the symmetric tridiagonal matrix. */
/* >          If COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          eigenvectors are desired, then  LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (max(1,2*N-2)) */
/* >          If COMPZ = 'N', then WORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  the algorithm has failed to find all the eigenvalues in */
/* >                a total of 30*N iterations; if INFO = i, then i */
/* >                elements of E have not converged to zero; on exit, D */
/* >                and E contain the elements of a symmetric tridiagonal */
/* >                matrix which is unitarily similar to the original */
/* >                matrix. */
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
/* Subroutine */ int csteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info, ftnlen compz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b, c__, f, g;
    static integer i__, j, k, l, m;
    static doublereal p, r__, s;
    static integer l1, ii, mm, lm1, mm1, nm1;
    static doublereal rt1, rt2, eps;
    static integer lsv;
    static doublereal tst, eps2;
    static integer lend, jtot;
    extern /* Subroutine */ int slae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int clasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer lendm1, lendp1;
    extern /* Subroutine */ int slaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal slapy2_(doublereal *, doublereal *);
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer lendsv;
    extern /* Subroutine */ int slartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal ssfmin;
    static integer nmaxit, icompz;
    static doublereal ssfmax;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int slasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);


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

#line 184 "csteqr.f"
    /* Parameter adjustments */
#line 184 "csteqr.f"
    --d__;
#line 184 "csteqr.f"
    --e;
#line 184 "csteqr.f"
    z_dim1 = *ldz;
#line 184 "csteqr.f"
    z_offset = 1 + z_dim1;
#line 184 "csteqr.f"
    z__ -= z_offset;
#line 184 "csteqr.f"
    --work;
#line 184 "csteqr.f"

#line 184 "csteqr.f"
    /* Function Body */
#line 184 "csteqr.f"
    *info = 0;

#line 186 "csteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 187 "csteqr.f"
	icompz = 0;
#line 188 "csteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 189 "csteqr.f"
	icompz = 1;
#line 190 "csteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 191 "csteqr.f"
	icompz = 2;
#line 192 "csteqr.f"
    } else {
#line 193 "csteqr.f"
	icompz = -1;
#line 194 "csteqr.f"
    }
#line 195 "csteqr.f"
    if (icompz < 0) {
#line 196 "csteqr.f"
	*info = -1;
#line 197 "csteqr.f"
    } else if (*n < 0) {
#line 198 "csteqr.f"
	*info = -2;
#line 199 "csteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 201 "csteqr.f"
	*info = -6;
#line 202 "csteqr.f"
    }
#line 203 "csteqr.f"
    if (*info != 0) {
#line 204 "csteqr.f"
	i__1 = -(*info);
#line 204 "csteqr.f"
	xerbla_("CSTEQR", &i__1, (ftnlen)6);
#line 205 "csteqr.f"
	return 0;
#line 206 "csteqr.f"
    }

/*     Quick return if possible */

#line 210 "csteqr.f"
    if (*n == 0) {
#line 210 "csteqr.f"
	return 0;
#line 210 "csteqr.f"
    }

#line 213 "csteqr.f"
    if (*n == 1) {
#line 214 "csteqr.f"
	if (icompz == 2) {
#line 214 "csteqr.f"
	    i__1 = z_dim1 + 1;
#line 214 "csteqr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 214 "csteqr.f"
	}
#line 216 "csteqr.f"
	return 0;
#line 217 "csteqr.f"
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

#line 221 "csteqr.f"
    eps = slamch_("E", (ftnlen)1);
/* Computing 2nd power */
#line 222 "csteqr.f"
    d__1 = eps;
#line 222 "csteqr.f"
    eps2 = d__1 * d__1;
#line 223 "csteqr.f"
    safmin = slamch_("S", (ftnlen)1);
#line 224 "csteqr.f"
    safmax = 1. / safmin;
#line 225 "csteqr.f"
    ssfmax = sqrt(safmax) / 3.;
#line 226 "csteqr.f"
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal */
/*     matrix. */

#line 231 "csteqr.f"
    if (icompz == 2) {
#line 231 "csteqr.f"
	claset_("Full", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)4);
#line 231 "csteqr.f"
    }

#line 234 "csteqr.f"
    nmaxit = *n * 30;
#line 235 "csteqr.f"
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

#line 241 "csteqr.f"
    l1 = 1;
#line 242 "csteqr.f"
    nm1 = *n - 1;

#line 244 "csteqr.f"
L10:
#line 245 "csteqr.f"
    if (l1 > *n) {
#line 245 "csteqr.f"
	goto L160;
#line 245 "csteqr.f"
    }
#line 247 "csteqr.f"
    if (l1 > 1) {
#line 247 "csteqr.f"
	e[l1 - 1] = 0.;
#line 247 "csteqr.f"
    }
#line 249 "csteqr.f"
    if (l1 <= nm1) {
#line 250 "csteqr.f"
	i__1 = nm1;
#line 250 "csteqr.f"
	for (m = l1; m <= i__1; ++m) {
#line 251 "csteqr.f"
	    tst = (d__1 = e[m], abs(d__1));
#line 252 "csteqr.f"
	    if (tst == 0.) {
#line 252 "csteqr.f"
		goto L30;
#line 252 "csteqr.f"
	    }
#line 254 "csteqr.f"
	    if (tst <= sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m 
		    + 1], abs(d__2))) * eps) {
#line 256 "csteqr.f"
		e[m] = 0.;
#line 257 "csteqr.f"
		goto L30;
#line 258 "csteqr.f"
	    }
#line 259 "csteqr.f"
/* L20: */
#line 259 "csteqr.f"
	}
#line 260 "csteqr.f"
    }
#line 261 "csteqr.f"
    m = *n;

#line 263 "csteqr.f"
L30:
#line 264 "csteqr.f"
    l = l1;
#line 265 "csteqr.f"
    lsv = l;
#line 266 "csteqr.f"
    lend = m;
#line 267 "csteqr.f"
    lendsv = lend;
#line 268 "csteqr.f"
    l1 = m + 1;
#line 269 "csteqr.f"
    if (lend == l) {
#line 269 "csteqr.f"
	goto L10;
#line 269 "csteqr.f"
    }

/*     Scale submatrix in rows and columns L to LEND */

#line 274 "csteqr.f"
    i__1 = lend - l + 1;
#line 274 "csteqr.f"
    anorm = slanst_("I", &i__1, &d__[l], &e[l], (ftnlen)1);
#line 275 "csteqr.f"
    iscale = 0;
#line 276 "csteqr.f"
    if (anorm == 0.) {
#line 276 "csteqr.f"
	goto L10;
#line 276 "csteqr.f"
    }
#line 278 "csteqr.f"
    if (anorm > ssfmax) {
#line 279 "csteqr.f"
	iscale = 1;
#line 280 "csteqr.f"
	i__1 = lend - l + 1;
#line 280 "csteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 282 "csteqr.f"
	i__1 = lend - l;
#line 282 "csteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 284 "csteqr.f"
    } else if (anorm < ssfmin) {
#line 285 "csteqr.f"
	iscale = 2;
#line 286 "csteqr.f"
	i__1 = lend - l + 1;
#line 286 "csteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 288 "csteqr.f"
	i__1 = lend - l;
#line 288 "csteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 290 "csteqr.f"
    }

/*     Choose between QL and QR iteration */

#line 294 "csteqr.f"
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
#line 295 "csteqr.f"
	lend = lsv;
#line 296 "csteqr.f"
	l = lendsv;
#line 297 "csteqr.f"
    }

#line 299 "csteqr.f"
    if (lend > l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

#line 305 "csteqr.f"
L40:
#line 306 "csteqr.f"
	if (l != lend) {
#line 307 "csteqr.f"
	    lendm1 = lend - 1;
#line 308 "csteqr.f"
	    i__1 = lendm1;
#line 308 "csteqr.f"
	    for (m = l; m <= i__1; ++m) {
/* Computing 2nd power */
#line 309 "csteqr.f"
		d__2 = (d__1 = e[m], abs(d__1));
#line 309 "csteqr.f"
		tst = d__2 * d__2;
#line 310 "csteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			+ 1], abs(d__2)) + safmin) {
#line 310 "csteqr.f"
		    goto L60;
#line 310 "csteqr.f"
		}
#line 312 "csteqr.f"
/* L50: */
#line 312 "csteqr.f"
	    }
#line 313 "csteqr.f"
	}

#line 315 "csteqr.f"
	m = lend;

#line 317 "csteqr.f"
L60:
#line 318 "csteqr.f"
	if (m < lend) {
#line 318 "csteqr.f"
	    e[m] = 0.;
#line 318 "csteqr.f"
	}
#line 320 "csteqr.f"
	p = d__[l];
#line 321 "csteqr.f"
	if (m == l) {
#line 321 "csteqr.f"
	    goto L80;
#line 321 "csteqr.f"
	}

/*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 327 "csteqr.f"
	if (m == l + 1) {
#line 328 "csteqr.f"
	    if (icompz > 0) {
#line 329 "csteqr.f"
		slaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
#line 330 "csteqr.f"
		work[l] = c__;
#line 331 "csteqr.f"
		work[*n - 1 + l] = s;
#line 332 "csteqr.f"
		clasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 334 "csteqr.f"
	    } else {
#line 335 "csteqr.f"
		slae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
#line 336 "csteqr.f"
	    }
#line 337 "csteqr.f"
	    d__[l] = rt1;
#line 338 "csteqr.f"
	    d__[l + 1] = rt2;
#line 339 "csteqr.f"
	    e[l] = 0.;
#line 340 "csteqr.f"
	    l += 2;
#line 341 "csteqr.f"
	    if (l <= lend) {
#line 341 "csteqr.f"
		goto L40;
#line 341 "csteqr.f"
	    }
#line 343 "csteqr.f"
	    goto L140;
#line 344 "csteqr.f"
	}

#line 346 "csteqr.f"
	if (jtot == nmaxit) {
#line 346 "csteqr.f"
	    goto L140;
#line 346 "csteqr.f"
	}
#line 348 "csteqr.f"
	++jtot;

/*        Form shift. */

#line 352 "csteqr.f"
	g = (d__[l + 1] - p) / (e[l] * 2.);
#line 353 "csteqr.f"
	r__ = slapy2_(&g, &c_b41);
#line 354 "csteqr.f"
	g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

#line 356 "csteqr.f"
	s = 1.;
#line 357 "csteqr.f"
	c__ = 1.;
#line 358 "csteqr.f"
	p = 0.;

/*        Inner loop */

#line 362 "csteqr.f"
	mm1 = m - 1;
#line 363 "csteqr.f"
	i__1 = l;
#line 363 "csteqr.f"
	for (i__ = mm1; i__ >= i__1; --i__) {
#line 364 "csteqr.f"
	    f = s * e[i__];
#line 365 "csteqr.f"
	    b = c__ * e[i__];
#line 366 "csteqr.f"
	    slartg_(&g, &f, &c__, &s, &r__);
#line 367 "csteqr.f"
	    if (i__ != m - 1) {
#line 367 "csteqr.f"
		e[i__ + 1] = r__;
#line 367 "csteqr.f"
	    }
#line 369 "csteqr.f"
	    g = d__[i__ + 1] - p;
#line 370 "csteqr.f"
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
#line 371 "csteqr.f"
	    p = s * r__;
#line 372 "csteqr.f"
	    d__[i__ + 1] = g + p;
#line 373 "csteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 377 "csteqr.f"
	    if (icompz > 0) {
#line 378 "csteqr.f"
		work[i__] = c__;
#line 379 "csteqr.f"
		work[*n - 1 + i__] = -s;
#line 380 "csteqr.f"
	    }

#line 382 "csteqr.f"
/* L70: */
#line 382 "csteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 386 "csteqr.f"
	if (icompz > 0) {
#line 387 "csteqr.f"
	    mm = m - l + 1;
#line 388 "csteqr.f"
	    clasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 390 "csteqr.f"
	}

#line 392 "csteqr.f"
	d__[l] -= p;
#line 393 "csteqr.f"
	e[l] = g;
#line 394 "csteqr.f"
	goto L40;

/*        Eigenvalue found. */

#line 398 "csteqr.f"
L80:
#line 399 "csteqr.f"
	d__[l] = p;

#line 401 "csteqr.f"
	++l;
#line 402 "csteqr.f"
	if (l <= lend) {
#line 402 "csteqr.f"
	    goto L40;
#line 402 "csteqr.f"
	}
#line 404 "csteqr.f"
	goto L140;

#line 406 "csteqr.f"
    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

#line 412 "csteqr.f"
L90:
#line 413 "csteqr.f"
	if (l != lend) {
#line 414 "csteqr.f"
	    lendp1 = lend + 1;
#line 415 "csteqr.f"
	    i__1 = lendp1;
#line 415 "csteqr.f"
	    for (m = l; m >= i__1; --m) {
/* Computing 2nd power */
#line 416 "csteqr.f"
		d__2 = (d__1 = e[m - 1], abs(d__1));
#line 416 "csteqr.f"
		tst = d__2 * d__2;
#line 417 "csteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			- 1], abs(d__2)) + safmin) {
#line 417 "csteqr.f"
		    goto L110;
#line 417 "csteqr.f"
		}
#line 419 "csteqr.f"
/* L100: */
#line 419 "csteqr.f"
	    }
#line 420 "csteqr.f"
	}

#line 422 "csteqr.f"
	m = lend;

#line 424 "csteqr.f"
L110:
#line 425 "csteqr.f"
	if (m > lend) {
#line 425 "csteqr.f"
	    e[m - 1] = 0.;
#line 425 "csteqr.f"
	}
#line 427 "csteqr.f"
	p = d__[l];
#line 428 "csteqr.f"
	if (m == l) {
#line 428 "csteqr.f"
	    goto L130;
#line 428 "csteqr.f"
	}

/*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 434 "csteqr.f"
	if (m == l - 1) {
#line 435 "csteqr.f"
	    if (icompz > 0) {
#line 436 "csteqr.f"
		slaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
#line 437 "csteqr.f"
		work[m] = c__;
#line 438 "csteqr.f"
		work[*n - 1 + m] = s;
#line 439 "csteqr.f"
		clasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 441 "csteqr.f"
	    } else {
#line 442 "csteqr.f"
		slae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
#line 443 "csteqr.f"
	    }
#line 444 "csteqr.f"
	    d__[l - 1] = rt1;
#line 445 "csteqr.f"
	    d__[l] = rt2;
#line 446 "csteqr.f"
	    e[l - 1] = 0.;
#line 447 "csteqr.f"
	    l += -2;
#line 448 "csteqr.f"
	    if (l >= lend) {
#line 448 "csteqr.f"
		goto L90;
#line 448 "csteqr.f"
	    }
#line 450 "csteqr.f"
	    goto L140;
#line 451 "csteqr.f"
	}

#line 453 "csteqr.f"
	if (jtot == nmaxit) {
#line 453 "csteqr.f"
	    goto L140;
#line 453 "csteqr.f"
	}
#line 455 "csteqr.f"
	++jtot;

/*        Form shift. */

#line 459 "csteqr.f"
	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
#line 460 "csteqr.f"
	r__ = slapy2_(&g, &c_b41);
#line 461 "csteqr.f"
	g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

#line 463 "csteqr.f"
	s = 1.;
#line 464 "csteqr.f"
	c__ = 1.;
#line 465 "csteqr.f"
	p = 0.;

/*        Inner loop */

#line 469 "csteqr.f"
	lm1 = l - 1;
#line 470 "csteqr.f"
	i__1 = lm1;
#line 470 "csteqr.f"
	for (i__ = m; i__ <= i__1; ++i__) {
#line 471 "csteqr.f"
	    f = s * e[i__];
#line 472 "csteqr.f"
	    b = c__ * e[i__];
#line 473 "csteqr.f"
	    slartg_(&g, &f, &c__, &s, &r__);
#line 474 "csteqr.f"
	    if (i__ != m) {
#line 474 "csteqr.f"
		e[i__ - 1] = r__;
#line 474 "csteqr.f"
	    }
#line 476 "csteqr.f"
	    g = d__[i__] - p;
#line 477 "csteqr.f"
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
#line 478 "csteqr.f"
	    p = s * r__;
#line 479 "csteqr.f"
	    d__[i__] = g + p;
#line 480 "csteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 484 "csteqr.f"
	    if (icompz > 0) {
#line 485 "csteqr.f"
		work[i__] = c__;
#line 486 "csteqr.f"
		work[*n - 1 + i__] = s;
#line 487 "csteqr.f"
	    }

#line 489 "csteqr.f"
/* L120: */
#line 489 "csteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 493 "csteqr.f"
	if (icompz > 0) {
#line 494 "csteqr.f"
	    mm = l - m + 1;
#line 495 "csteqr.f"
	    clasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 497 "csteqr.f"
	}

#line 499 "csteqr.f"
	d__[l] -= p;
#line 500 "csteqr.f"
	e[lm1] = g;
#line 501 "csteqr.f"
	goto L90;

/*        Eigenvalue found. */

#line 505 "csteqr.f"
L130:
#line 506 "csteqr.f"
	d__[l] = p;

#line 508 "csteqr.f"
	--l;
#line 509 "csteqr.f"
	if (l >= lend) {
#line 509 "csteqr.f"
	    goto L90;
#line 509 "csteqr.f"
	}
#line 511 "csteqr.f"
	goto L140;

#line 513 "csteqr.f"
    }

/*     Undo scaling if necessary */

#line 517 "csteqr.f"
L140:
#line 518 "csteqr.f"
    if (iscale == 1) {
#line 519 "csteqr.f"
	i__1 = lendsv - lsv + 1;
#line 519 "csteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 521 "csteqr.f"
	i__1 = lendsv - lsv;
#line 521 "csteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 523 "csteqr.f"
    } else if (iscale == 2) {
#line 524 "csteqr.f"
	i__1 = lendsv - lsv + 1;
#line 524 "csteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 526 "csteqr.f"
	i__1 = lendsv - lsv;
#line 526 "csteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 528 "csteqr.f"
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

#line 533 "csteqr.f"
    if (jtot == nmaxit) {
#line 534 "csteqr.f"
	i__1 = *n - 1;
#line 534 "csteqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "csteqr.f"
	    if (e[i__] != 0.) {
#line 535 "csteqr.f"
		++(*info);
#line 535 "csteqr.f"
	    }
#line 537 "csteqr.f"
/* L150: */
#line 537 "csteqr.f"
	}
#line 538 "csteqr.f"
	return 0;
#line 539 "csteqr.f"
    }
#line 540 "csteqr.f"
    goto L10;

/*     Order eigenvalues and eigenvectors. */

#line 544 "csteqr.f"
L160:
#line 545 "csteqr.f"
    if (icompz == 0) {

/*        Use Quick Sort */

#line 549 "csteqr.f"
	slasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 551 "csteqr.f"
    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

#line 555 "csteqr.f"
	i__1 = *n;
#line 555 "csteqr.f"
	for (ii = 2; ii <= i__1; ++ii) {
#line 556 "csteqr.f"
	    i__ = ii - 1;
#line 557 "csteqr.f"
	    k = i__;
#line 558 "csteqr.f"
	    p = d__[i__];
#line 559 "csteqr.f"
	    i__2 = *n;
#line 559 "csteqr.f"
	    for (j = ii; j <= i__2; ++j) {
#line 560 "csteqr.f"
		if (d__[j] < p) {
#line 561 "csteqr.f"
		    k = j;
#line 562 "csteqr.f"
		    p = d__[j];
#line 563 "csteqr.f"
		}
#line 564 "csteqr.f"
/* L170: */
#line 564 "csteqr.f"
	    }
#line 565 "csteqr.f"
	    if (k != i__) {
#line 566 "csteqr.f"
		d__[k] = d__[i__];
#line 567 "csteqr.f"
		d__[i__] = p;
#line 568 "csteqr.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
#line 569 "csteqr.f"
	    }
#line 570 "csteqr.f"
/* L180: */
#line 570 "csteqr.f"
	}
#line 571 "csteqr.f"
    }
#line 572 "csteqr.f"
    return 0;

/*     End of CSTEQR */

} /* csteqr_ */

