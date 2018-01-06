#line 1 "ssteqr.f"
/* ssteqr.f -- translated by f2c (version 20100827).
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

#line 1 "ssteqr.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SSTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the implicit QL or QR method. */
/* > The eigenvectors of a full or band symmetric matrix can also be found */
/* > if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to */
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
/* >          Z is REAL array, dimension (LDZ, N) */
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
/* >                matrix which is orthogonally similar to the original */
/* >                matrix. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssteqr_(char *compz, integer *n, doublereal *d__, 
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
    extern /* Subroutine */ int slae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int slasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), sswap_(integer *, doublereal *, integer *
	    , doublereal *, integer *);
    static integer lendm1, lendp1;
    extern /* Subroutine */ int slaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal slapy2_(doublereal *, doublereal *);
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer lendsv;
    extern /* Subroutine */ int slartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), slaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
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

#line 179 "ssteqr.f"
    /* Parameter adjustments */
#line 179 "ssteqr.f"
    --d__;
#line 179 "ssteqr.f"
    --e;
#line 179 "ssteqr.f"
    z_dim1 = *ldz;
#line 179 "ssteqr.f"
    z_offset = 1 + z_dim1;
#line 179 "ssteqr.f"
    z__ -= z_offset;
#line 179 "ssteqr.f"
    --work;
#line 179 "ssteqr.f"

#line 179 "ssteqr.f"
    /* Function Body */
#line 179 "ssteqr.f"
    *info = 0;

#line 181 "ssteqr.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 182 "ssteqr.f"
	icompz = 0;
#line 183 "ssteqr.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 184 "ssteqr.f"
	icompz = 1;
#line 185 "ssteqr.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 186 "ssteqr.f"
	icompz = 2;
#line 187 "ssteqr.f"
    } else {
#line 188 "ssteqr.f"
	icompz = -1;
#line 189 "ssteqr.f"
    }
#line 190 "ssteqr.f"
    if (icompz < 0) {
#line 191 "ssteqr.f"
	*info = -1;
#line 192 "ssteqr.f"
    } else if (*n < 0) {
#line 193 "ssteqr.f"
	*info = -2;
#line 194 "ssteqr.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 196 "ssteqr.f"
	*info = -6;
#line 197 "ssteqr.f"
    }
#line 198 "ssteqr.f"
    if (*info != 0) {
#line 199 "ssteqr.f"
	i__1 = -(*info);
#line 199 "ssteqr.f"
	xerbla_("SSTEQR", &i__1, (ftnlen)6);
#line 200 "ssteqr.f"
	return 0;
#line 201 "ssteqr.f"
    }

/*     Quick return if possible */

#line 205 "ssteqr.f"
    if (*n == 0) {
#line 205 "ssteqr.f"
	return 0;
#line 205 "ssteqr.f"
    }

#line 208 "ssteqr.f"
    if (*n == 1) {
#line 209 "ssteqr.f"
	if (icompz == 2) {
#line 209 "ssteqr.f"
	    z__[z_dim1 + 1] = 1.;
#line 209 "ssteqr.f"
	}
#line 211 "ssteqr.f"
	return 0;
#line 212 "ssteqr.f"
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

#line 216 "ssteqr.f"
    eps = slamch_("E", (ftnlen)1);
/* Computing 2nd power */
#line 217 "ssteqr.f"
    d__1 = eps;
#line 217 "ssteqr.f"
    eps2 = d__1 * d__1;
#line 218 "ssteqr.f"
    safmin = slamch_("S", (ftnlen)1);
#line 219 "ssteqr.f"
    safmax = 1. / safmin;
#line 220 "ssteqr.f"
    ssfmax = sqrt(safmax) / 3.;
#line 221 "ssteqr.f"
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal */
/*     matrix. */

#line 226 "ssteqr.f"
    if (icompz == 2) {
#line 226 "ssteqr.f"
	slaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz, (ftnlen)4);
#line 226 "ssteqr.f"
    }

#line 229 "ssteqr.f"
    nmaxit = *n * 30;
#line 230 "ssteqr.f"
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

#line 236 "ssteqr.f"
    l1 = 1;
#line 237 "ssteqr.f"
    nm1 = *n - 1;

#line 239 "ssteqr.f"
L10:
#line 240 "ssteqr.f"
    if (l1 > *n) {
#line 240 "ssteqr.f"
	goto L160;
#line 240 "ssteqr.f"
    }
#line 242 "ssteqr.f"
    if (l1 > 1) {
#line 242 "ssteqr.f"
	e[l1 - 1] = 0.;
#line 242 "ssteqr.f"
    }
#line 244 "ssteqr.f"
    if (l1 <= nm1) {
#line 245 "ssteqr.f"
	i__1 = nm1;
#line 245 "ssteqr.f"
	for (m = l1; m <= i__1; ++m) {
#line 246 "ssteqr.f"
	    tst = (d__1 = e[m], abs(d__1));
#line 247 "ssteqr.f"
	    if (tst == 0.) {
#line 247 "ssteqr.f"
		goto L30;
#line 247 "ssteqr.f"
	    }
#line 249 "ssteqr.f"
	    if (tst <= sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m 
		    + 1], abs(d__2))) * eps) {
#line 251 "ssteqr.f"
		e[m] = 0.;
#line 252 "ssteqr.f"
		goto L30;
#line 253 "ssteqr.f"
	    }
#line 254 "ssteqr.f"
/* L20: */
#line 254 "ssteqr.f"
	}
#line 255 "ssteqr.f"
    }
#line 256 "ssteqr.f"
    m = *n;

#line 258 "ssteqr.f"
L30:
#line 259 "ssteqr.f"
    l = l1;
#line 260 "ssteqr.f"
    lsv = l;
#line 261 "ssteqr.f"
    lend = m;
#line 262 "ssteqr.f"
    lendsv = lend;
#line 263 "ssteqr.f"
    l1 = m + 1;
#line 264 "ssteqr.f"
    if (lend == l) {
#line 264 "ssteqr.f"
	goto L10;
#line 264 "ssteqr.f"
    }

/*     Scale submatrix in rows and columns L to LEND */

#line 269 "ssteqr.f"
    i__1 = lend - l + 1;
#line 269 "ssteqr.f"
    anorm = slanst_("M", &i__1, &d__[l], &e[l], (ftnlen)1);
#line 270 "ssteqr.f"
    iscale = 0;
#line 271 "ssteqr.f"
    if (anorm == 0.) {
#line 271 "ssteqr.f"
	goto L10;
#line 271 "ssteqr.f"
    }
#line 273 "ssteqr.f"
    if (anorm > ssfmax) {
#line 274 "ssteqr.f"
	iscale = 1;
#line 275 "ssteqr.f"
	i__1 = lend - l + 1;
#line 275 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 277 "ssteqr.f"
	i__1 = lend - l;
#line 277 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 279 "ssteqr.f"
    } else if (anorm < ssfmin) {
#line 280 "ssteqr.f"
	iscale = 2;
#line 281 "ssteqr.f"
	i__1 = lend - l + 1;
#line 281 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 283 "ssteqr.f"
	i__1 = lend - l;
#line 283 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 285 "ssteqr.f"
    }

/*     Choose between QL and QR iteration */

#line 289 "ssteqr.f"
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
#line 290 "ssteqr.f"
	lend = lsv;
#line 291 "ssteqr.f"
	l = lendsv;
#line 292 "ssteqr.f"
    }

#line 294 "ssteqr.f"
    if (lend > l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

#line 300 "ssteqr.f"
L40:
#line 301 "ssteqr.f"
	if (l != lend) {
#line 302 "ssteqr.f"
	    lendm1 = lend - 1;
#line 303 "ssteqr.f"
	    i__1 = lendm1;
#line 303 "ssteqr.f"
	    for (m = l; m <= i__1; ++m) {
/* Computing 2nd power */
#line 304 "ssteqr.f"
		d__2 = (d__1 = e[m], abs(d__1));
#line 304 "ssteqr.f"
		tst = d__2 * d__2;
#line 305 "ssteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			+ 1], abs(d__2)) + safmin) {
#line 305 "ssteqr.f"
		    goto L60;
#line 305 "ssteqr.f"
		}
#line 307 "ssteqr.f"
/* L50: */
#line 307 "ssteqr.f"
	    }
#line 308 "ssteqr.f"
	}

#line 310 "ssteqr.f"
	m = lend;

#line 312 "ssteqr.f"
L60:
#line 313 "ssteqr.f"
	if (m < lend) {
#line 313 "ssteqr.f"
	    e[m] = 0.;
#line 313 "ssteqr.f"
	}
#line 315 "ssteqr.f"
	p = d__[l];
#line 316 "ssteqr.f"
	if (m == l) {
#line 316 "ssteqr.f"
	    goto L80;
#line 316 "ssteqr.f"
	}

/*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 322 "ssteqr.f"
	if (m == l + 1) {
#line 323 "ssteqr.f"
	    if (icompz > 0) {
#line 324 "ssteqr.f"
		slaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
#line 325 "ssteqr.f"
		work[l] = c__;
#line 326 "ssteqr.f"
		work[*n - 1 + l] = s;
#line 327 "ssteqr.f"
		slasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 329 "ssteqr.f"
	    } else {
#line 330 "ssteqr.f"
		slae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
#line 331 "ssteqr.f"
	    }
#line 332 "ssteqr.f"
	    d__[l] = rt1;
#line 333 "ssteqr.f"
	    d__[l + 1] = rt2;
#line 334 "ssteqr.f"
	    e[l] = 0.;
#line 335 "ssteqr.f"
	    l += 2;
#line 336 "ssteqr.f"
	    if (l <= lend) {
#line 336 "ssteqr.f"
		goto L40;
#line 336 "ssteqr.f"
	    }
#line 338 "ssteqr.f"
	    goto L140;
#line 339 "ssteqr.f"
	}

#line 341 "ssteqr.f"
	if (jtot == nmaxit) {
#line 341 "ssteqr.f"
	    goto L140;
#line 341 "ssteqr.f"
	}
#line 343 "ssteqr.f"
	++jtot;

/*        Form shift. */

#line 347 "ssteqr.f"
	g = (d__[l + 1] - p) / (e[l] * 2.);
#line 348 "ssteqr.f"
	r__ = slapy2_(&g, &c_b10);
#line 349 "ssteqr.f"
	g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

#line 351 "ssteqr.f"
	s = 1.;
#line 352 "ssteqr.f"
	c__ = 1.;
#line 353 "ssteqr.f"
	p = 0.;

/*        Inner loop */

#line 357 "ssteqr.f"
	mm1 = m - 1;
#line 358 "ssteqr.f"
	i__1 = l;
#line 358 "ssteqr.f"
	for (i__ = mm1; i__ >= i__1; --i__) {
#line 359 "ssteqr.f"
	    f = s * e[i__];
#line 360 "ssteqr.f"
	    b = c__ * e[i__];
#line 361 "ssteqr.f"
	    slartg_(&g, &f, &c__, &s, &r__);
#line 362 "ssteqr.f"
	    if (i__ != m - 1) {
#line 362 "ssteqr.f"
		e[i__ + 1] = r__;
#line 362 "ssteqr.f"
	    }
#line 364 "ssteqr.f"
	    g = d__[i__ + 1] - p;
#line 365 "ssteqr.f"
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
#line 366 "ssteqr.f"
	    p = s * r__;
#line 367 "ssteqr.f"
	    d__[i__ + 1] = g + p;
#line 368 "ssteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 372 "ssteqr.f"
	    if (icompz > 0) {
#line 373 "ssteqr.f"
		work[i__] = c__;
#line 374 "ssteqr.f"
		work[*n - 1 + i__] = -s;
#line 375 "ssteqr.f"
	    }

#line 377 "ssteqr.f"
/* L70: */
#line 377 "ssteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 381 "ssteqr.f"
	if (icompz > 0) {
#line 382 "ssteqr.f"
	    mm = m - l + 1;
#line 383 "ssteqr.f"
	    slasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 385 "ssteqr.f"
	}

#line 387 "ssteqr.f"
	d__[l] -= p;
#line 388 "ssteqr.f"
	e[l] = g;
#line 389 "ssteqr.f"
	goto L40;

/*        Eigenvalue found. */

#line 393 "ssteqr.f"
L80:
#line 394 "ssteqr.f"
	d__[l] = p;

#line 396 "ssteqr.f"
	++l;
#line 397 "ssteqr.f"
	if (l <= lend) {
#line 397 "ssteqr.f"
	    goto L40;
#line 397 "ssteqr.f"
	}
#line 399 "ssteqr.f"
	goto L140;

#line 401 "ssteqr.f"
    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

#line 407 "ssteqr.f"
L90:
#line 408 "ssteqr.f"
	if (l != lend) {
#line 409 "ssteqr.f"
	    lendp1 = lend + 1;
#line 410 "ssteqr.f"
	    i__1 = lendp1;
#line 410 "ssteqr.f"
	    for (m = l; m >= i__1; --m) {
/* Computing 2nd power */
#line 411 "ssteqr.f"
		d__2 = (d__1 = e[m - 1], abs(d__1));
#line 411 "ssteqr.f"
		tst = d__2 * d__2;
#line 412 "ssteqr.f"
		if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
			- 1], abs(d__2)) + safmin) {
#line 412 "ssteqr.f"
		    goto L110;
#line 412 "ssteqr.f"
		}
#line 414 "ssteqr.f"
/* L100: */
#line 414 "ssteqr.f"
	    }
#line 415 "ssteqr.f"
	}

#line 417 "ssteqr.f"
	m = lend;

#line 419 "ssteqr.f"
L110:
#line 420 "ssteqr.f"
	if (m > lend) {
#line 420 "ssteqr.f"
	    e[m - 1] = 0.;
#line 420 "ssteqr.f"
	}
#line 422 "ssteqr.f"
	p = d__[l];
#line 423 "ssteqr.f"
	if (m == l) {
#line 423 "ssteqr.f"
	    goto L130;
#line 423 "ssteqr.f"
	}

/*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

#line 429 "ssteqr.f"
	if (m == l - 1) {
#line 430 "ssteqr.f"
	    if (icompz > 0) {
#line 431 "ssteqr.f"
		slaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
#line 432 "ssteqr.f"
		work[m] = c__;
#line 433 "ssteqr.f"
		work[*n - 1 + m] = s;
#line 434 "ssteqr.f"
		slasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 436 "ssteqr.f"
	    } else {
#line 437 "ssteqr.f"
		slae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
#line 438 "ssteqr.f"
	    }
#line 439 "ssteqr.f"
	    d__[l - 1] = rt1;
#line 440 "ssteqr.f"
	    d__[l] = rt2;
#line 441 "ssteqr.f"
	    e[l - 1] = 0.;
#line 442 "ssteqr.f"
	    l += -2;
#line 443 "ssteqr.f"
	    if (l >= lend) {
#line 443 "ssteqr.f"
		goto L90;
#line 443 "ssteqr.f"
	    }
#line 445 "ssteqr.f"
	    goto L140;
#line 446 "ssteqr.f"
	}

#line 448 "ssteqr.f"
	if (jtot == nmaxit) {
#line 448 "ssteqr.f"
	    goto L140;
#line 448 "ssteqr.f"
	}
#line 450 "ssteqr.f"
	++jtot;

/*        Form shift. */

#line 454 "ssteqr.f"
	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
#line 455 "ssteqr.f"
	r__ = slapy2_(&g, &c_b10);
#line 456 "ssteqr.f"
	g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

#line 458 "ssteqr.f"
	s = 1.;
#line 459 "ssteqr.f"
	c__ = 1.;
#line 460 "ssteqr.f"
	p = 0.;

/*        Inner loop */

#line 464 "ssteqr.f"
	lm1 = l - 1;
#line 465 "ssteqr.f"
	i__1 = lm1;
#line 465 "ssteqr.f"
	for (i__ = m; i__ <= i__1; ++i__) {
#line 466 "ssteqr.f"
	    f = s * e[i__];
#line 467 "ssteqr.f"
	    b = c__ * e[i__];
#line 468 "ssteqr.f"
	    slartg_(&g, &f, &c__, &s, &r__);
#line 469 "ssteqr.f"
	    if (i__ != m) {
#line 469 "ssteqr.f"
		e[i__ - 1] = r__;
#line 469 "ssteqr.f"
	    }
#line 471 "ssteqr.f"
	    g = d__[i__] - p;
#line 472 "ssteqr.f"
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
#line 473 "ssteqr.f"
	    p = s * r__;
#line 474 "ssteqr.f"
	    d__[i__] = g + p;
#line 475 "ssteqr.f"
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

#line 479 "ssteqr.f"
	    if (icompz > 0) {
#line 480 "ssteqr.f"
		work[i__] = c__;
#line 481 "ssteqr.f"
		work[*n - 1 + i__] = s;
#line 482 "ssteqr.f"
	    }

#line 484 "ssteqr.f"
/* L120: */
#line 484 "ssteqr.f"
	}

/*        If eigenvectors are desired, then apply saved rotations. */

#line 488 "ssteqr.f"
	if (icompz > 0) {
#line 489 "ssteqr.f"
	    mm = l - m + 1;
#line 490 "ssteqr.f"
	    slasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 492 "ssteqr.f"
	}

#line 494 "ssteqr.f"
	d__[l] -= p;
#line 495 "ssteqr.f"
	e[lm1] = g;
#line 496 "ssteqr.f"
	goto L90;

/*        Eigenvalue found. */

#line 500 "ssteqr.f"
L130:
#line 501 "ssteqr.f"
	d__[l] = p;

#line 503 "ssteqr.f"
	--l;
#line 504 "ssteqr.f"
	if (l >= lend) {
#line 504 "ssteqr.f"
	    goto L90;
#line 504 "ssteqr.f"
	}
#line 506 "ssteqr.f"
	goto L140;

#line 508 "ssteqr.f"
    }

/*     Undo scaling if necessary */

#line 512 "ssteqr.f"
L140:
#line 513 "ssteqr.f"
    if (iscale == 1) {
#line 514 "ssteqr.f"
	i__1 = lendsv - lsv + 1;
#line 514 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 516 "ssteqr.f"
	i__1 = lendsv - lsv;
#line 516 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 518 "ssteqr.f"
    } else if (iscale == 2) {
#line 519 "ssteqr.f"
	i__1 = lendsv - lsv + 1;
#line 519 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 521 "ssteqr.f"
	i__1 = lendsv - lsv;
#line 521 "ssteqr.f"
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info, (ftnlen)1);
#line 523 "ssteqr.f"
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

#line 528 "ssteqr.f"
    if (jtot < nmaxit) {
#line 528 "ssteqr.f"
	goto L10;
#line 528 "ssteqr.f"
    }
#line 530 "ssteqr.f"
    i__1 = *n - 1;
#line 530 "ssteqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 531 "ssteqr.f"
	if (e[i__] != 0.) {
#line 531 "ssteqr.f"
	    ++(*info);
#line 531 "ssteqr.f"
	}
#line 533 "ssteqr.f"
/* L150: */
#line 533 "ssteqr.f"
    }
#line 534 "ssteqr.f"
    goto L190;

/*     Order eigenvalues and eigenvectors. */

#line 538 "ssteqr.f"
L160:
#line 539 "ssteqr.f"
    if (icompz == 0) {

/*        Use Quick Sort */

#line 543 "ssteqr.f"
	slasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 545 "ssteqr.f"
    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

#line 549 "ssteqr.f"
	i__1 = *n;
#line 549 "ssteqr.f"
	for (ii = 2; ii <= i__1; ++ii) {
#line 550 "ssteqr.f"
	    i__ = ii - 1;
#line 551 "ssteqr.f"
	    k = i__;
#line 552 "ssteqr.f"
	    p = d__[i__];
#line 553 "ssteqr.f"
	    i__2 = *n;
#line 553 "ssteqr.f"
	    for (j = ii; j <= i__2; ++j) {
#line 554 "ssteqr.f"
		if (d__[j] < p) {
#line 555 "ssteqr.f"
		    k = j;
#line 556 "ssteqr.f"
		    p = d__[j];
#line 557 "ssteqr.f"
		}
#line 558 "ssteqr.f"
/* L170: */
#line 558 "ssteqr.f"
	    }
#line 559 "ssteqr.f"
	    if (k != i__) {
#line 560 "ssteqr.f"
		d__[k] = d__[i__];
#line 561 "ssteqr.f"
		d__[i__] = p;
#line 562 "ssteqr.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
#line 563 "ssteqr.f"
	    }
#line 564 "ssteqr.f"
/* L180: */
#line 564 "ssteqr.f"
	}
#line 565 "ssteqr.f"
    }

#line 567 "ssteqr.f"
L190:
#line 568 "ssteqr.f"
    return 0;

/*     End of SSTEQR */

} /* ssteqr_ */

