#line 1 "zsteqr.f"
/* zsteqr.f -- translated by f2c (version 20100827).
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

#line 1 "zsteqr.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b41 = 1.;

/* > \brief \b ZSTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ) */
/*       COMPLEX*16         Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the implicit QL or QR method. */
/* > The eigenvectors of a full or band complex Hermitian matrix can also */
/* > be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2)) */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsteqr_(char *compz, integer *n, doublereal *d__, 
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
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int zlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), dlaev2_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static integer lendsv;
    static doublereal ssfmin;
    static integer nmaxit, icompz;
    static doublereal ssfmax;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 184 "zsteqr.f"
    /* Parameter adjustments */
#line 184 "zsteqr.f"
    --d__;
#line 184 "zsteqr.f"
    --e;
#line 184 "zsteqr.f"
    z_dim1 = *ldz;
#line 184 "zsteqr.f"
    z_offset = 1 + z_dim1;
#line 184 "zsteqr.f"
    z__ -= z_offset;
#line 184 "zsteqr.f"
    --work;
#line 184 "zsteqr.f"

#line 184 "zsteqr.f"
    /* Function Body */
#line 184 "zsteqr.f"
    *info = 0;

#line 186 "zsteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 187 "zsteqr.f"
	icompz = 0;
#line 188 "zsteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 189 "zsteqr.f"
	icompz = 1;
#line 190 "zsteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 191 "zsteqr.f"
	icompz = 2;
#line 192 "zsteqr.f"
    } else {
#line 193 "zsteqr.f"
	icompz = -1;
#line 194 "zsteqr.f"
    }
#line 195 "zsteqr.f"
    if (icompz < 0) {
#line 196 "zsteqr.f"
	*info = -1;
#line 197 "zsteqr.f"
    } else if (*n < 0) {
#line 198 "zsteqr.f"
	*info = -2;
#line 199 "zsteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 201 "zsteqr.f"
	*info = -6;
#line 202 "zsteqr.f"
    }
#line 203 "zsteqr.f"
    if (*info != 0) {
#line 204 "zsteqr.f"
	i__1 = -(*info);
#line 204 "zsteqr.f"
	xerbla_("ZSTEQR", &i__1, (ftnlen)6);
#line 205 "zsteqr.f"
	return 0;
#line 206 "zsteqr.f"
    }

/*     Quick return if possible */

#line 210 "zsteqr.f"
    if (*n == 0) {
#line 210 "zsteqr.f"
	return 0;
#line 210 "zsteqr.f"
    }

#line 213 "zsteqr.f"
    if (*n == 1) {
#line 214 "zsteqr.f"
	if (icompz == 2) {
#line 214 "zsteqr.f"
	    i__1 = z_dim1 + 1;
#line 214 "zsteqr.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 214 "zsteqr.f"
	}
#line 216 "zsteqr.f"
	return 0;
#line 217 "zsteqr.f"
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

#line 221 "zsteqr.f"
    eps = dlamch_("E", (ftnlen)1);
/* Computing 2nd power */
#line 222 "zsteqr.f"
    d__1 = eps;
#line 222 "zsteqr.f"
    eps2 = d__1 * d__1;
#line 223 "zsteqr.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 224 "zsteqr.f"
    safmax = 1. / safmin;
#line 225 "zsteqr.f"
    ssfmax = sqrt(safmax) / 3.;
#line 226 "zsteqr.f"
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal */
/*     matrix. */

#line 231 "zsteqr.f"
    if (icompz == 2) {
#line 231 "zsteqr.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)4);
#line 231 "zsteqr.f"
    }

#line 234 "zsteqr.f"
    nmaxit = *n * 30;
#line 235 "zsteqr.f"
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

#line 241 "zsteqr.f"
    l1 = 1;
#line 242 "zsteqr.f"
    nm1 = *n - 1;

#line 244 "zsteqr.f"
L10:
#line 245 "zsteqr.f"
    if (l1 > *n) {
#line 245 "zsteqr.f"
	goto L160;
#line 245 "zsteqr.f"
    }
#line 247 "zsteqr.f"
    if (l1 > 1) {
#line 247 "zsteqr.f"
	e[l1 - 1] = 0.;
#line 247 "zsteqr.f"
    }
#line 249 "zsteqr.f"
    if (l1 <= nm1) {
#line 250 "zsteqr.f"
	i__1 = nm1;
#line 250 "zsteqr.f"
	for (m = l1; m <= i__1; ++m) {
#line 251 "zsteqr.f"
	    tst = (d__1 = e[m], abs(d__1));
#line 252 "zsteqr.f"
	    if (tst == 0.) {
#line 252 "zsteqr.f"
		goto L30;
#line 252 "zsteqr.f"
	    }
#line 254 "zsteqr.f"
	    if (tst <= sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m 
		    + 1], abs(d__2))) * eps) {
#line 256 "zsteqr.f"
		e[m] = 0.;
#line 257 "zsteqr.f"
		goto L30;
#line 258 "zsteqr.f"
	    }
#line 259 "zsteqr.f"
/* L20: */
#line 259 "zsteqr.f"
	}
#line 260 "zsteqr.f"
    }
#line 261 "zsteqr.f"
    m = *n;

#line 263 "zsteqr.f"
L30:
#line 264 "zsteqr.f"
    l = l1;
#line 265 "zsteqr.f"
    lsv = l;
#line 266 "zsteqr.f"
    lend = m;
#line 267 "zsteqr.f"
    lendsv = lend;
#line 268 "zsteqr.f"
    l1 = m + 1;
#line 269 "zsteqr.f"
    if (lend == l) {
#line 269 "zsteqr.f"
	goto L10;
#line 269 "zsteqr.f"
    }

/*     Scale submatrix in rows and columns L to LEND */

#line 274 "zsteqr.f"
    i__1 = lend - l + 1;
#line 274 "zsteqr.f"
    anorm = dlanst_("I", &i__1, &d__[l], &e[l], (ftnlen)1);
#line 275 "zsteqr.f"
    iscale = 0;
#line 276 "zsteqr.f"
    if (anorm == 0.) {
#line 276 "zsteqr.f"
	goto L10;
#line 276 "zsteqr.f"
    }
#line 278 "zsteqr.f"
    if (anorm > ssfmax) {
#line 279 "zsteqr.f"
	iscale = 1;
#line 280 "zsteqr.f"
	i__1 = lend - l + 1;
#line 280 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 282 "zsteqr.f"
	i__1 = lend - l;
#line 282 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 284 "zsteqr.f"
    } else if (anorm < ssfmin) {
#line 285 "zsteqr.f"
	iscale = 2;
#line 286 "zsteqr.f"
	i__1 = lend - l + 1;
#line 286 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 288 "zsteqr.f"
	i__1 = lend - l;
#line 288 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 290 "zsteqr.f"
    }

/*     Choose between QL and QR iteration */

#line 294 "zsteqr.f"
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
#line 295 "zsteqr.f"
	lend = lsv;
#line 296 "zsteqr.f"
	l = lendsv;
#line 297 "zsteqr.f"
    }

#line 299 "zsteqr.f"
    if (lend > l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

#line 305 "zsteqr.f"
L40:
#line 306 "zsteqr.f"
	if (l != lend) {
#line 307 "zsteqr.f"
	    lendm1 = lend - 1;
#line 308 "zsteqr.f"
	    i__1 = lendm1;
#line 308 "zsteqr.f"
	    for (m = l; m <= i__1; ++m) {
/* Computing 2nd power */
#line 309 "zsteqr.f"
		d__2 = (d__1 = e[m], abs(d__1));
#line 309 "zsteqr.f"
		tst = d__2 * d__2;
#line 310 "zsteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			+ 1], abs(d__2)) + safmin) {
#line 310 "zsteqr.f"
		    goto L60;
#line 310 "zsteqr.f"
		}
#line 312 "zsteqr.f"
/* L50: */
#line 312 "zsteqr.f"
	    }
#line 313 "zsteqr.f"
	}

#line 315 "zsteqr.f"
	m = lend;

#line 317 "zsteqr.f"
L60:
#line 318 "zsteqr.f"
	if (m < lend) {
#line 318 "zsteqr.f"
	    e[m] = 0.;
#line 318 "zsteqr.f"
	}
#line 320 "zsteqr.f"
	p = d__[l];
#line 321 "zsteqr.f"
	if (m == l) {
#line 321 "zsteqr.f"
	    goto L80;
#line 321 "zsteqr.f"
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 327 "zsteqr.f"
	if (m == l + 1) {
#line 328 "zsteqr.f"
	    if (icompz > 0) {
#line 329 "zsteqr.f"
		dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
#line 330 "zsteqr.f"
		work[l] = c__;
#line 331 "zsteqr.f"
		work[*n - 1 + l] = s;
#line 332 "zsteqr.f"
		zlasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 334 "zsteqr.f"
	    } else {
#line 335 "zsteqr.f"
		dlae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
#line 336 "zsteqr.f"
	    }
#line 337 "zsteqr.f"
	    d__[l] = rt1;
#line 338 "zsteqr.f"
	    d__[l + 1] = rt2;
#line 339 "zsteqr.f"
	    e[l] = 0.;
#line 340 "zsteqr.f"
	    l += 2;
#line 341 "zsteqr.f"
	    if (l <= lend) {
#line 341 "zsteqr.f"
		goto L40;
#line 341 "zsteqr.f"
	    }
#line 343 "zsteqr.f"
	    goto L140;
#line 344 "zsteqr.f"
	}

#line 346 "zsteqr.f"
	if (jtot == nmaxit) {
#line 346 "zsteqr.f"
	    goto L140;
#line 346 "zsteqr.f"
	}
#line 348 "zsteqr.f"
	++jtot;

/*        Form shift. */

#line 352 "zsteqr.f"
	g = (d__[l + 1] - p) / (e[l] * 2.);
#line 353 "zsteqr.f"
	r__ = dlapy2_(&g, &c_b41);
#line 354 "zsteqr.f"
	g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

#line 356 "zsteqr.f"
	s = 1.;
#line 357 "zsteqr.f"
	c__ = 1.;
#line 358 "zsteqr.f"
	p = 0.;

/*        Inner loop */

#line 362 "zsteqr.f"
	mm1 = m - 1;
#line 363 "zsteqr.f"
	i__1 = l;
#line 363 "zsteqr.f"
	for (i__ = mm1; i__ >= i__1; --i__) {
#line 364 "zsteqr.f"
	    f = s * e[i__];
#line 365 "zsteqr.f"
	    b = c__ * e[i__];
#line 366 "zsteqr.f"
	    dlartg_(&g, &f, &c__, &s, &r__);
#line 367 "zsteqr.f"
	    if (i__ != m - 1) {
#line 367 "zsteqr.f"
		e[i__ + 1] = r__;
#line 367 "zsteqr.f"
	    }
#line 369 "zsteqr.f"
	    g = d__[i__ + 1] - p;
#line 370 "zsteqr.f"
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
#line 371 "zsteqr.f"
	    p = s * r__;
#line 372 "zsteqr.f"
	    d__[i__ + 1] = g + p;
#line 373 "zsteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 377 "zsteqr.f"
	    if (icompz > 0) {
#line 378 "zsteqr.f"
		work[i__] = c__;
#line 379 "zsteqr.f"
		work[*n - 1 + i__] = -s;
#line 380 "zsteqr.f"
	    }

#line 382 "zsteqr.f"
/* L70: */
#line 382 "zsteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 386 "zsteqr.f"
	if (icompz > 0) {
#line 387 "zsteqr.f"
	    mm = m - l + 1;
#line 388 "zsteqr.f"
	    zlasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 390 "zsteqr.f"
	}

#line 392 "zsteqr.f"
	d__[l] -= p;
#line 393 "zsteqr.f"
	e[l] = g;
#line 394 "zsteqr.f"
	goto L40;

/*        Eigenvalue found. */

#line 398 "zsteqr.f"
L80:
#line 399 "zsteqr.f"
	d__[l] = p;

#line 401 "zsteqr.f"
	++l;
#line 402 "zsteqr.f"
	if (l <= lend) {
#line 402 "zsteqr.f"
	    goto L40;
#line 402 "zsteqr.f"
	}
#line 404 "zsteqr.f"
	goto L140;

#line 406 "zsteqr.f"
    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

#line 412 "zsteqr.f"
L90:
#line 413 "zsteqr.f"
	if (l != lend) {
#line 414 "zsteqr.f"
	    lendp1 = lend + 1;
#line 415 "zsteqr.f"
	    i__1 = lendp1;
#line 415 "zsteqr.f"
	    for (m = l; m >= i__1; --m) {
/* Computing 2nd power */
#line 416 "zsteqr.f"
		d__2 = (d__1 = e[m - 1], abs(d__1));
#line 416 "zsteqr.f"
		tst = d__2 * d__2;
#line 417 "zsteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			- 1], abs(d__2)) + safmin) {
#line 417 "zsteqr.f"
		    goto L110;
#line 417 "zsteqr.f"
		}
#line 419 "zsteqr.f"
/* L100: */
#line 419 "zsteqr.f"
	    }
#line 420 "zsteqr.f"
	}

#line 422 "zsteqr.f"
	m = lend;

#line 424 "zsteqr.f"
L110:
#line 425 "zsteqr.f"
	if (m > lend) {
#line 425 "zsteqr.f"
	    e[m - 1] = 0.;
#line 425 "zsteqr.f"
	}
#line 427 "zsteqr.f"
	p = d__[l];
#line 428 "zsteqr.f"
	if (m == l) {
#line 428 "zsteqr.f"
	    goto L130;
#line 428 "zsteqr.f"
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 434 "zsteqr.f"
	if (m == l - 1) {
#line 435 "zsteqr.f"
	    if (icompz > 0) {
#line 436 "zsteqr.f"
		dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
#line 437 "zsteqr.f"
		work[m] = c__;
#line 438 "zsteqr.f"
		work[*n - 1 + m] = s;
#line 439 "zsteqr.f"
		zlasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 441 "zsteqr.f"
	    } else {
#line 442 "zsteqr.f"
		dlae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
#line 443 "zsteqr.f"
	    }
#line 444 "zsteqr.f"
	    d__[l - 1] = rt1;
#line 445 "zsteqr.f"
	    d__[l] = rt2;
#line 446 "zsteqr.f"
	    e[l - 1] = 0.;
#line 447 "zsteqr.f"
	    l += -2;
#line 448 "zsteqr.f"
	    if (l >= lend) {
#line 448 "zsteqr.f"
		goto L90;
#line 448 "zsteqr.f"
	    }
#line 450 "zsteqr.f"
	    goto L140;
#line 451 "zsteqr.f"
	}

#line 453 "zsteqr.f"
	if (jtot == nmaxit) {
#line 453 "zsteqr.f"
	    goto L140;
#line 453 "zsteqr.f"
	}
#line 455 "zsteqr.f"
	++jtot;

/*        Form shift. */

#line 459 "zsteqr.f"
	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
#line 460 "zsteqr.f"
	r__ = dlapy2_(&g, &c_b41);
#line 461 "zsteqr.f"
	g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

#line 463 "zsteqr.f"
	s = 1.;
#line 464 "zsteqr.f"
	c__ = 1.;
#line 465 "zsteqr.f"
	p = 0.;

/*        Inner loop */

#line 469 "zsteqr.f"
	lm1 = l - 1;
#line 470 "zsteqr.f"
	i__1 = lm1;
#line 470 "zsteqr.f"
	for (i__ = m; i__ <= i__1; ++i__) {
#line 471 "zsteqr.f"
	    f = s * e[i__];
#line 472 "zsteqr.f"
	    b = c__ * e[i__];
#line 473 "zsteqr.f"
	    dlartg_(&g, &f, &c__, &s, &r__);
#line 474 "zsteqr.f"
	    if (i__ != m) {
#line 474 "zsteqr.f"
		e[i__ - 1] = r__;
#line 474 "zsteqr.f"
	    }
#line 476 "zsteqr.f"
	    g = d__[i__] - p;
#line 477 "zsteqr.f"
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
#line 478 "zsteqr.f"
	    p = s * r__;
#line 479 "zsteqr.f"
	    d__[i__] = g + p;
#line 480 "zsteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 484 "zsteqr.f"
	    if (icompz > 0) {
#line 485 "zsteqr.f"
		work[i__] = c__;
#line 486 "zsteqr.f"
		work[*n - 1 + i__] = s;
#line 487 "zsteqr.f"
	    }

#line 489 "zsteqr.f"
/* L120: */
#line 489 "zsteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 493 "zsteqr.f"
	if (icompz > 0) {
#line 494 "zsteqr.f"
	    mm = l - m + 1;
#line 495 "zsteqr.f"
	    zlasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 497 "zsteqr.f"
	}

#line 499 "zsteqr.f"
	d__[l] -= p;
#line 500 "zsteqr.f"
	e[lm1] = g;
#line 501 "zsteqr.f"
	goto L90;

/*        Eigenvalue found. */

#line 505 "zsteqr.f"
L130:
#line 506 "zsteqr.f"
	d__[l] = p;

#line 508 "zsteqr.f"
	--l;
#line 509 "zsteqr.f"
	if (l >= lend) {
#line 509 "zsteqr.f"
	    goto L90;
#line 509 "zsteqr.f"
	}
#line 511 "zsteqr.f"
	goto L140;

#line 513 "zsteqr.f"
    }

/*     Undo scaling if necessary */

#line 517 "zsteqr.f"
L140:
#line 518 "zsteqr.f"
    if (iscale == 1) {
#line 519 "zsteqr.f"
	i__1 = lendsv - lsv + 1;
#line 519 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 521 "zsteqr.f"
	i__1 = lendsv - lsv;
#line 521 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 523 "zsteqr.f"
    } else if (iscale == 2) {
#line 524 "zsteqr.f"
	i__1 = lendsv - lsv + 1;
#line 524 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 526 "zsteqr.f"
	i__1 = lendsv - lsv;
#line 526 "zsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 528 "zsteqr.f"
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

#line 533 "zsteqr.f"
    if (jtot == nmaxit) {
#line 534 "zsteqr.f"
	i__1 = *n - 1;
#line 534 "zsteqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "zsteqr.f"
	    if (e[i__] != 0.) {
#line 535 "zsteqr.f"
		++(*info);
#line 535 "zsteqr.f"
	    }
#line 537 "zsteqr.f"
/* L150: */
#line 537 "zsteqr.f"
	}
#line 538 "zsteqr.f"
	return 0;
#line 539 "zsteqr.f"
    }
#line 540 "zsteqr.f"
    goto L10;

/*     Order eigenvalues and eigenvectors. */

#line 544 "zsteqr.f"
L160:
#line 545 "zsteqr.f"
    if (icompz == 0) {

/*        Use Quick Sort */

#line 549 "zsteqr.f"
	dlasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 551 "zsteqr.f"
    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

#line 555 "zsteqr.f"
	i__1 = *n;
#line 555 "zsteqr.f"
	for (ii = 2; ii <= i__1; ++ii) {
#line 556 "zsteqr.f"
	    i__ = ii - 1;
#line 557 "zsteqr.f"
	    k = i__;
#line 558 "zsteqr.f"
	    p = d__[i__];
#line 559 "zsteqr.f"
	    i__2 = *n;
#line 559 "zsteqr.f"
	    for (j = ii; j <= i__2; ++j) {
#line 560 "zsteqr.f"
		if (d__[j] < p) {
#line 561 "zsteqr.f"
		    k = j;
#line 562 "zsteqr.f"
		    p = d__[j];
#line 563 "zsteqr.f"
		}
#line 564 "zsteqr.f"
/* L170: */
#line 564 "zsteqr.f"
	    }
#line 565 "zsteqr.f"
	    if (k != i__) {
#line 566 "zsteqr.f"
		d__[k] = d__[i__];
#line 567 "zsteqr.f"
		d__[i__] = p;
#line 568 "zsteqr.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
#line 569 "zsteqr.f"
	    }
#line 570 "zsteqr.f"
/* L180: */
#line 570 "zsteqr.f"
	}
#line 571 "zsteqr.f"
    }
#line 572 "zsteqr.f"
    return 0;

/*     End of ZSTEQR */

} /* zsteqr_ */

