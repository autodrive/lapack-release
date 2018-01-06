#line 1 "dsteqr.f"
/* dsteqr.f -- translated by f2c (version 20100827).
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

#line 1 "dsteqr.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DSTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the implicit QL or QR method. */
/* > The eigenvectors of a full or band symmetric matrix can also be found */
/* > if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to */
/* > tridiagonal form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only. */
/* >          = 'V':  Compute eigenvalues and eigenvectors of the original */
/* >                  symmetric matrix.  On entry, Z must contain the */
/* >                  orthogonal matrix used to reduce the original matrix */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          On entry, if  COMPZ = 'V', then Z contains the orthogonal */
/* >          matrix used in the reduction to tridiagonal form. */
/* >          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the */
/* >          orthonormal eigenvectors of the original symmetric matrix, */
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
/* >                matrix which is orthogonally similar to the original */
/* >                matrix. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
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
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen);
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

#line 179 "dsteqr.f"
    /* Parameter adjustments */
#line 179 "dsteqr.f"
    --d__;
#line 179 "dsteqr.f"
    --e;
#line 179 "dsteqr.f"
    z_dim1 = *ldz;
#line 179 "dsteqr.f"
    z_offset = 1 + z_dim1;
#line 179 "dsteqr.f"
    z__ -= z_offset;
#line 179 "dsteqr.f"
    --work;
#line 179 "dsteqr.f"

#line 179 "dsteqr.f"
    /* Function Body */
#line 179 "dsteqr.f"
    *info = 0;

#line 181 "dsteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 182 "dsteqr.f"
	icompz = 0;
#line 183 "dsteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 184 "dsteqr.f"
	icompz = 1;
#line 185 "dsteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 186 "dsteqr.f"
	icompz = 2;
#line 187 "dsteqr.f"
    } else {
#line 188 "dsteqr.f"
	icompz = -1;
#line 189 "dsteqr.f"
    }
#line 190 "dsteqr.f"
    if (icompz < 0) {
#line 191 "dsteqr.f"
	*info = -1;
#line 192 "dsteqr.f"
    } else if (*n < 0) {
#line 193 "dsteqr.f"
	*info = -2;
#line 194 "dsteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 196 "dsteqr.f"
	*info = -6;
#line 197 "dsteqr.f"
    }
#line 198 "dsteqr.f"
    if (*info != 0) {
#line 199 "dsteqr.f"
	i__1 = -(*info);
#line 199 "dsteqr.f"
	xerbla_("DSTEQR", &i__1, (ftnlen)6);
#line 200 "dsteqr.f"
	return 0;
#line 201 "dsteqr.f"
    }

/*     Quick return if possible */

#line 205 "dsteqr.f"
    if (*n == 0) {
#line 205 "dsteqr.f"
	return 0;
#line 205 "dsteqr.f"
    }

#line 208 "dsteqr.f"
    if (*n == 1) {
#line 209 "dsteqr.f"
	if (icompz == 2) {
#line 209 "dsteqr.f"
	    z__[z_dim1 + 1] = 1.;
#line 209 "dsteqr.f"
	}
#line 211 "dsteqr.f"
	return 0;
#line 212 "dsteqr.f"
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

#line 216 "dsteqr.f"
    eps = dlamch_("E", (ftnlen)1);
/* Computing 2nd power */
#line 217 "dsteqr.f"
    d__1 = eps;
#line 217 "dsteqr.f"
    eps2 = d__1 * d__1;
#line 218 "dsteqr.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 219 "dsteqr.f"
    safmax = 1. / safmin;
#line 220 "dsteqr.f"
    ssfmax = sqrt(safmax) / 3.;
#line 221 "dsteqr.f"
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal */
/*     matrix. */

#line 226 "dsteqr.f"
    if (icompz == 2) {
#line 226 "dsteqr.f"
	dlaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz, (ftnlen)4);
#line 226 "dsteqr.f"
    }

#line 229 "dsteqr.f"
    nmaxit = *n * 30;
#line 230 "dsteqr.f"
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

#line 236 "dsteqr.f"
    l1 = 1;
#line 237 "dsteqr.f"
    nm1 = *n - 1;

#line 239 "dsteqr.f"
L10:
#line 240 "dsteqr.f"
    if (l1 > *n) {
#line 240 "dsteqr.f"
	goto L160;
#line 240 "dsteqr.f"
    }
#line 242 "dsteqr.f"
    if (l1 > 1) {
#line 242 "dsteqr.f"
	e[l1 - 1] = 0.;
#line 242 "dsteqr.f"
    }
#line 244 "dsteqr.f"
    if (l1 <= nm1) {
#line 245 "dsteqr.f"
	i__1 = nm1;
#line 245 "dsteqr.f"
	for (m = l1; m <= i__1; ++m) {
#line 246 "dsteqr.f"
	    tst = (d__1 = e[m], abs(d__1));
#line 247 "dsteqr.f"
	    if (tst == 0.) {
#line 247 "dsteqr.f"
		goto L30;
#line 247 "dsteqr.f"
	    }
#line 249 "dsteqr.f"
	    if (tst <= sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m 
		    + 1], abs(d__2))) * eps) {
#line 251 "dsteqr.f"
		e[m] = 0.;
#line 252 "dsteqr.f"
		goto L30;
#line 253 "dsteqr.f"
	    }
#line 254 "dsteqr.f"
/* L20: */
#line 254 "dsteqr.f"
	}
#line 255 "dsteqr.f"
    }
#line 256 "dsteqr.f"
    m = *n;

#line 258 "dsteqr.f"
L30:
#line 259 "dsteqr.f"
    l = l1;
#line 260 "dsteqr.f"
    lsv = l;
#line 261 "dsteqr.f"
    lend = m;
#line 262 "dsteqr.f"
    lendsv = lend;
#line 263 "dsteqr.f"
    l1 = m + 1;
#line 264 "dsteqr.f"
    if (lend == l) {
#line 264 "dsteqr.f"
	goto L10;
#line 264 "dsteqr.f"
    }

/*     Scale submatrix in rows and columns L to LEND */

#line 269 "dsteqr.f"
    i__1 = lend - l + 1;
#line 269 "dsteqr.f"
    anorm = dlanst_("M", &i__1, &d__[l], &e[l], (ftnlen)1);
#line 270 "dsteqr.f"
    iscale = 0;
#line 271 "dsteqr.f"
    if (anorm == 0.) {
#line 271 "dsteqr.f"
	goto L10;
#line 271 "dsteqr.f"
    }
#line 273 "dsteqr.f"
    if (anorm > ssfmax) {
#line 274 "dsteqr.f"
	iscale = 1;
#line 275 "dsteqr.f"
	i__1 = lend - l + 1;
#line 275 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 277 "dsteqr.f"
	i__1 = lend - l;
#line 277 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 279 "dsteqr.f"
    } else if (anorm < ssfmin) {
#line 280 "dsteqr.f"
	iscale = 2;
#line 281 "dsteqr.f"
	i__1 = lend - l + 1;
#line 281 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 283 "dsteqr.f"
	i__1 = lend - l;
#line 283 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 285 "dsteqr.f"
    }

/*     Choose between QL and QR iteration */

#line 289 "dsteqr.f"
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
#line 290 "dsteqr.f"
	lend = lsv;
#line 291 "dsteqr.f"
	l = lendsv;
#line 292 "dsteqr.f"
    }

#line 294 "dsteqr.f"
    if (lend > l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

#line 300 "dsteqr.f"
L40:
#line 301 "dsteqr.f"
	if (l != lend) {
#line 302 "dsteqr.f"
	    lendm1 = lend - 1;
#line 303 "dsteqr.f"
	    i__1 = lendm1;
#line 303 "dsteqr.f"
	    for (m = l; m <= i__1; ++m) {
/* Computing 2nd power */
#line 304 "dsteqr.f"
		d__2 = (d__1 = e[m], abs(d__1));
#line 304 "dsteqr.f"
		tst = d__2 * d__2;
#line 305 "dsteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			+ 1], abs(d__2)) + safmin) {
#line 305 "dsteqr.f"
		    goto L60;
#line 305 "dsteqr.f"
		}
#line 307 "dsteqr.f"
/* L50: */
#line 307 "dsteqr.f"
	    }
#line 308 "dsteqr.f"
	}

#line 310 "dsteqr.f"
	m = lend;

#line 312 "dsteqr.f"
L60:
#line 313 "dsteqr.f"
	if (m < lend) {
#line 313 "dsteqr.f"
	    e[m] = 0.;
#line 313 "dsteqr.f"
	}
#line 315 "dsteqr.f"
	p = d__[l];
#line 316 "dsteqr.f"
	if (m == l) {
#line 316 "dsteqr.f"
	    goto L80;
#line 316 "dsteqr.f"
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 322 "dsteqr.f"
	if (m == l + 1) {
#line 323 "dsteqr.f"
	    if (icompz > 0) {
#line 324 "dsteqr.f"
		dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
#line 325 "dsteqr.f"
		work[l] = c__;
#line 326 "dsteqr.f"
		work[*n - 1 + l] = s;
#line 327 "dsteqr.f"
		dlasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 329 "dsteqr.f"
	    } else {
#line 330 "dsteqr.f"
		dlae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
#line 331 "dsteqr.f"
	    }
#line 332 "dsteqr.f"
	    d__[l] = rt1;
#line 333 "dsteqr.f"
	    d__[l + 1] = rt2;
#line 334 "dsteqr.f"
	    e[l] = 0.;
#line 335 "dsteqr.f"
	    l += 2;
#line 336 "dsteqr.f"
	    if (l <= lend) {
#line 336 "dsteqr.f"
		goto L40;
#line 336 "dsteqr.f"
	    }
#line 338 "dsteqr.f"
	    goto L140;
#line 339 "dsteqr.f"
	}

#line 341 "dsteqr.f"
	if (jtot == nmaxit) {
#line 341 "dsteqr.f"
	    goto L140;
#line 341 "dsteqr.f"
	}
#line 343 "dsteqr.f"
	++jtot;

/*        Form shift. */

#line 347 "dsteqr.f"
	g = (d__[l + 1] - p) / (e[l] * 2.);
#line 348 "dsteqr.f"
	r__ = dlapy2_(&g, &c_b10);
#line 349 "dsteqr.f"
	g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

#line 351 "dsteqr.f"
	s = 1.;
#line 352 "dsteqr.f"
	c__ = 1.;
#line 353 "dsteqr.f"
	p = 0.;

/*        Inner loop */

#line 357 "dsteqr.f"
	mm1 = m - 1;
#line 358 "dsteqr.f"
	i__1 = l;
#line 358 "dsteqr.f"
	for (i__ = mm1; i__ >= i__1; --i__) {
#line 359 "dsteqr.f"
	    f = s * e[i__];
#line 360 "dsteqr.f"
	    b = c__ * e[i__];
#line 361 "dsteqr.f"
	    dlartg_(&g, &f, &c__, &s, &r__);
#line 362 "dsteqr.f"
	    if (i__ != m - 1) {
#line 362 "dsteqr.f"
		e[i__ + 1] = r__;
#line 362 "dsteqr.f"
	    }
#line 364 "dsteqr.f"
	    g = d__[i__ + 1] - p;
#line 365 "dsteqr.f"
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
#line 366 "dsteqr.f"
	    p = s * r__;
#line 367 "dsteqr.f"
	    d__[i__ + 1] = g + p;
#line 368 "dsteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 372 "dsteqr.f"
	    if (icompz > 0) {
#line 373 "dsteqr.f"
		work[i__] = c__;
#line 374 "dsteqr.f"
		work[*n - 1 + i__] = -s;
#line 375 "dsteqr.f"
	    }

#line 377 "dsteqr.f"
/* L70: */
#line 377 "dsteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 381 "dsteqr.f"
	if (icompz > 0) {
#line 382 "dsteqr.f"
	    mm = m - l + 1;
#line 383 "dsteqr.f"
	    dlasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 385 "dsteqr.f"
	}

#line 387 "dsteqr.f"
	d__[l] -= p;
#line 388 "dsteqr.f"
	e[l] = g;
#line 389 "dsteqr.f"
	goto L40;

/*        Eigenvalue found. */

#line 393 "dsteqr.f"
L80:
#line 394 "dsteqr.f"
	d__[l] = p;

#line 396 "dsteqr.f"
	++l;
#line 397 "dsteqr.f"
	if (l <= lend) {
#line 397 "dsteqr.f"
	    goto L40;
#line 397 "dsteqr.f"
	}
#line 399 "dsteqr.f"
	goto L140;

#line 401 "dsteqr.f"
    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

#line 407 "dsteqr.f"
L90:
#line 408 "dsteqr.f"
	if (l != lend) {
#line 409 "dsteqr.f"
	    lendp1 = lend + 1;
#line 410 "dsteqr.f"
	    i__1 = lendp1;
#line 410 "dsteqr.f"
	    for (m = l; m >= i__1; --m) {
/* Computing 2nd power */
#line 411 "dsteqr.f"
		d__2 = (d__1 = e[m - 1], abs(d__1));
#line 411 "dsteqr.f"
		tst = d__2 * d__2;
#line 412 "dsteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			- 1], abs(d__2)) + safmin) {
#line 412 "dsteqr.f"
		    goto L110;
#line 412 "dsteqr.f"
		}
#line 414 "dsteqr.f"
/* L100: */
#line 414 "dsteqr.f"
	    }
#line 415 "dsteqr.f"
	}

#line 417 "dsteqr.f"
	m = lend;

#line 419 "dsteqr.f"
L110:
#line 420 "dsteqr.f"
	if (m > lend) {
#line 420 "dsteqr.f"
	    e[m - 1] = 0.;
#line 420 "dsteqr.f"
	}
#line 422 "dsteqr.f"
	p = d__[l];
#line 423 "dsteqr.f"
	if (m == l) {
#line 423 "dsteqr.f"
	    goto L130;
#line 423 "dsteqr.f"
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 429 "dsteqr.f"
	if (m == l - 1) {
#line 430 "dsteqr.f"
	    if (icompz > 0) {
#line 431 "dsteqr.f"
		dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
#line 432 "dsteqr.f"
		work[m] = c__;
#line 433 "dsteqr.f"
		work[*n - 1 + m] = s;
#line 434 "dsteqr.f"
		dlasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 436 "dsteqr.f"
	    } else {
#line 437 "dsteqr.f"
		dlae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
#line 438 "dsteqr.f"
	    }
#line 439 "dsteqr.f"
	    d__[l - 1] = rt1;
#line 440 "dsteqr.f"
	    d__[l] = rt2;
#line 441 "dsteqr.f"
	    e[l - 1] = 0.;
#line 442 "dsteqr.f"
	    l += -2;
#line 443 "dsteqr.f"
	    if (l >= lend) {
#line 443 "dsteqr.f"
		goto L90;
#line 443 "dsteqr.f"
	    }
#line 445 "dsteqr.f"
	    goto L140;
#line 446 "dsteqr.f"
	}

#line 448 "dsteqr.f"
	if (jtot == nmaxit) {
#line 448 "dsteqr.f"
	    goto L140;
#line 448 "dsteqr.f"
	}
#line 450 "dsteqr.f"
	++jtot;

/*        Form shift. */

#line 454 "dsteqr.f"
	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
#line 455 "dsteqr.f"
	r__ = dlapy2_(&g, &c_b10);
#line 456 "dsteqr.f"
	g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

#line 458 "dsteqr.f"
	s = 1.;
#line 459 "dsteqr.f"
	c__ = 1.;
#line 460 "dsteqr.f"
	p = 0.;

/*        Inner loop */

#line 464 "dsteqr.f"
	lm1 = l - 1;
#line 465 "dsteqr.f"
	i__1 = lm1;
#line 465 "dsteqr.f"
	for (i__ = m; i__ <= i__1; ++i__) {
#line 466 "dsteqr.f"
	    f = s * e[i__];
#line 467 "dsteqr.f"
	    b = c__ * e[i__];
#line 468 "dsteqr.f"
	    dlartg_(&g, &f, &c__, &s, &r__);
#line 469 "dsteqr.f"
	    if (i__ != m) {
#line 469 "dsteqr.f"
		e[i__ - 1] = r__;
#line 469 "dsteqr.f"
	    }
#line 471 "dsteqr.f"
	    g = d__[i__] - p;
#line 472 "dsteqr.f"
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
#line 473 "dsteqr.f"
	    p = s * r__;
#line 474 "dsteqr.f"
	    d__[i__] = g + p;
#line 475 "dsteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 479 "dsteqr.f"
	    if (icompz > 0) {
#line 480 "dsteqr.f"
		work[i__] = c__;
#line 481 "dsteqr.f"
		work[*n - 1 + i__] = s;
#line 482 "dsteqr.f"
	    }

#line 484 "dsteqr.f"
/* L120: */
#line 484 "dsteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 488 "dsteqr.f"
	if (icompz > 0) {
#line 489 "dsteqr.f"
	    mm = l - m + 1;
#line 490 "dsteqr.f"
	    dlasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 492 "dsteqr.f"
	}

#line 494 "dsteqr.f"
	d__[l] -= p;
#line 495 "dsteqr.f"
	e[lm1] = g;
#line 496 "dsteqr.f"
	goto L90;

/*        Eigenvalue found. */

#line 500 "dsteqr.f"
L130:
#line 501 "dsteqr.f"
	d__[l] = p;

#line 503 "dsteqr.f"
	--l;
#line 504 "dsteqr.f"
	if (l >= lend) {
#line 504 "dsteqr.f"
	    goto L90;
#line 504 "dsteqr.f"
	}
#line 506 "dsteqr.f"
	goto L140;

#line 508 "dsteqr.f"
    }

/*     Undo scaling if necessary */

#line 512 "dsteqr.f"
L140:
#line 513 "dsteqr.f"
    if (iscale == 1) {
#line 514 "dsteqr.f"
	i__1 = lendsv - lsv + 1;
#line 514 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 516 "dsteqr.f"
	i__1 = lendsv - lsv;
#line 516 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 518 "dsteqr.f"
    } else if (iscale == 2) {
#line 519 "dsteqr.f"
	i__1 = lendsv - lsv + 1;
#line 519 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 521 "dsteqr.f"
	i__1 = lendsv - lsv;
#line 521 "dsteqr.f"
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 523 "dsteqr.f"
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

#line 528 "dsteqr.f"
    if (jtot < nmaxit) {
#line 528 "dsteqr.f"
	goto L10;
#line 528 "dsteqr.f"
    }
#line 530 "dsteqr.f"
    i__1 = *n - 1;
#line 530 "dsteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 531 "dsteqr.f"
	if (e[i__] != 0.) {
#line 531 "dsteqr.f"
	    ++(*info);
#line 531 "dsteqr.f"
	}
#line 533 "dsteqr.f"
/* L150: */
#line 533 "dsteqr.f"
    }
#line 534 "dsteqr.f"
    goto L190;

/*     Order eigenvalues and eigenvectors. */

#line 538 "dsteqr.f"
L160:
#line 539 "dsteqr.f"
    if (icompz == 0) {

/*        Use Quick Sort */

#line 543 "dsteqr.f"
	dlasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 545 "dsteqr.f"
    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

#line 549 "dsteqr.f"
	i__1 = *n;
#line 549 "dsteqr.f"
	for (ii = 2; ii <= i__1; ++ii) {
#line 550 "dsteqr.f"
	    i__ = ii - 1;
#line 551 "dsteqr.f"
	    k = i__;
#line 552 "dsteqr.f"
	    p = d__[i__];
#line 553 "dsteqr.f"
	    i__2 = *n;
#line 553 "dsteqr.f"
	    for (j = ii; j <= i__2; ++j) {
#line 554 "dsteqr.f"
		if (d__[j] < p) {
#line 555 "dsteqr.f"
		    k = j;
#line 556 "dsteqr.f"
		    p = d__[j];
#line 557 "dsteqr.f"
		}
#line 558 "dsteqr.f"
/* L170: */
#line 558 "dsteqr.f"
	    }
#line 559 "dsteqr.f"
	    if (k != i__) {
#line 560 "dsteqr.f"
		d__[k] = d__[i__];
#line 561 "dsteqr.f"
		d__[i__] = p;
#line 562 "dsteqr.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
#line 563 "dsteqr.f"
	    }
#line 564 "dsteqr.f"
/* L180: */
#line 564 "dsteqr.f"
	}
#line 565 "dsteqr.f"
    }

#line 567 "dsteqr.f"
L190:
#line 568 "dsteqr.f"
    return 0;

/*     End of DSTEQR */

} /* dsteqr_ */

