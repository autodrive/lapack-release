#line 1 "dlaein.f"
/* dlaein.f -- translated by f2c (version 20100827).
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

#line 1 "dlaein.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse 
iteration. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B, */
/*                          LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            NOINIT, RIGHTV */
/*       INTEGER            INFO, LDB, LDH, N */
/*       DOUBLE PRECISION   BIGNUM, EPS3, SMLNUM, WI, WR */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), H( LDH, * ), VI( * ), VR( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAEIN uses inverse iteration to find a right or left eigenvector */
/* > corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg */
/* > matrix H. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] RIGHTV */
/* > \verbatim */
/* >          RIGHTV is LOGICAL */
/* >          = .TRUE. : compute right eigenvector; */
/* >          = .FALSE.: compute left eigenvector. */
/* > \endverbatim */
/* > */
/* > \param[in] NOINIT */
/* > \verbatim */
/* >          NOINIT is LOGICAL */
/* >          = .TRUE. : no initial vector supplied in (VR,VI). */
/* >          = .FALSE.: initial vector supplied in (VR,VI). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION array, dimension (LDH,N) */
/* >          The upper Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H.  LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] WR */
/* > \verbatim */
/* >          WR is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION */
/* >          The real and imaginary parts of the eigenvalue of H whose */
/* >          corresponding right or left eigenvector is to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in,out] VI */
/* > \verbatim */
/* >          VI is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain */
/* >          a real starting vector for inverse iteration using the real */
/* >          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI */
/* >          must contain the real and imaginary parts of a complex */
/* >          starting vector for inverse iteration using the complex */
/* >          eigenvalue (WR,WI); otherwise VR and VI need not be set. */
/* >          On exit, if WI = 0.0 (real eigenvalue), VR contains the */
/* >          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue), */
/* >          VR and VI contain the real and imaginary parts of the */
/* >          computed complex eigenvector. The eigenvector is normalized */
/* >          so that the component of largest magnitude has magnitude 1; */
/* >          here the magnitude of a complex number (x,y) is taken to be */
/* >          |x| + |y|. */
/* >          VI is not referenced if WI = 0.0. */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= N+1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] EPS3 */
/* > \verbatim */
/* >          EPS3 is DOUBLE PRECISION */
/* >          A small machine-dependent value which is used to perturb */
/* >          close eigenvalues, and to replace zero pivots. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLNUM */
/* > \verbatim */
/* >          SMLNUM is DOUBLE PRECISION */
/* >          A machine-dependent value close to the underflow threshold. */
/* > \endverbatim */
/* > */
/* > \param[in] BIGNUM */
/* > \verbatim */
/* >          BIGNUM is DOUBLE PRECISION */
/* >          A machine-dependent value close to the overflow threshold. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          = 1:  inverse iteration did not converge; VR is set to the */
/* >                last iterate, and so is VI if WI.ne.0.0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaein_(logical *rightv, logical *noinit, integer *n, 
	doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, 
	doublereal *vr, doublereal *vi, doublereal *b, integer *ldb, 
	doublereal *work, doublereal *eps3, doublereal *smlnum, doublereal *
	bignum, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal w, x, y;
    static integer i1, i2, i3;
    static doublereal w1, ei, ej, xi, xr, rec;
    static integer its, ierr;
    static doublereal temp, norm, vmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static char trans[1];
    static doublereal vcrit, rootn, vnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal absbii, absbjj;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlatrs_(
	    char *, char *, char *, char *, integer *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static char normin[1];
    static doublereal nrmsml, growto;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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

#line 216 "dlaein.f"
    /* Parameter adjustments */
#line 216 "dlaein.f"
    h_dim1 = *ldh;
#line 216 "dlaein.f"
    h_offset = 1 + h_dim1;
#line 216 "dlaein.f"
    h__ -= h_offset;
#line 216 "dlaein.f"
    --vr;
#line 216 "dlaein.f"
    --vi;
#line 216 "dlaein.f"
    b_dim1 = *ldb;
#line 216 "dlaein.f"
    b_offset = 1 + b_dim1;
#line 216 "dlaein.f"
    b -= b_offset;
#line 216 "dlaein.f"
    --work;
#line 216 "dlaein.f"

#line 216 "dlaein.f"
    /* Function Body */
#line 216 "dlaein.f"
    *info = 0;

/*     GROWTO is the threshold used in the acceptance test for an */
/*     eigenvector. */

#line 221 "dlaein.f"
    rootn = sqrt((doublereal) (*n));
#line 222 "dlaein.f"
    growto = .1 / rootn;
/* Computing MAX */
#line 223 "dlaein.f"
    d__1 = 1., d__2 = *eps3 * rootn;
#line 223 "dlaein.f"
    nrmsml = max(d__1,d__2) * *smlnum;

/*     Form B = H - (WR,WI)*I (except that the subdiagonal elements and */
/*     the imaginary parts of the diagonal elements are not stored). */

#line 228 "dlaein.f"
    i__1 = *n;
#line 228 "dlaein.f"
    for (j = 1; j <= i__1; ++j) {
#line 229 "dlaein.f"
	i__2 = j - 1;
#line 229 "dlaein.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 230 "dlaein.f"
	    b[i__ + j * b_dim1] = h__[i__ + j * h_dim1];
#line 231 "dlaein.f"
/* L10: */
#line 231 "dlaein.f"
	}
#line 232 "dlaein.f"
	b[j + j * b_dim1] = h__[j + j * h_dim1] - *wr;
#line 233 "dlaein.f"
/* L20: */
#line 233 "dlaein.f"
    }

#line 235 "dlaein.f"
    if (*wi == 0.) {

/*        Real eigenvalue. */

#line 239 "dlaein.f"
	if (*noinit) {

/*           Set initial vector. */

#line 243 "dlaein.f"
	    i__1 = *n;
#line 243 "dlaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "dlaein.f"
		vr[i__] = *eps3;
#line 245 "dlaein.f"
/* L30: */
#line 245 "dlaein.f"
	    }
#line 246 "dlaein.f"
	} else {

/*           Scale supplied initial vector. */

#line 250 "dlaein.f"
	    vnorm = dnrm2_(n, &vr[1], &c__1);
#line 251 "dlaein.f"
	    d__1 = *eps3 * rootn / max(vnorm,nrmsml);
#line 251 "dlaein.f"
	    dscal_(n, &d__1, &vr[1], &c__1);
#line 253 "dlaein.f"
	}

#line 255 "dlaein.f"
	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacing zero */
/*           pivots by EPS3. */

#line 260 "dlaein.f"
	    i__1 = *n - 1;
#line 260 "dlaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 261 "dlaein.f"
		ei = h__[i__ + 1 + i__ * h_dim1];
#line 262 "dlaein.f"
		if ((d__1 = b[i__ + i__ * b_dim1], abs(d__1)) < abs(ei)) {

/*                 Interchange rows and eliminate. */

#line 266 "dlaein.f"
		    x = b[i__ + i__ * b_dim1] / ei;
#line 267 "dlaein.f"
		    b[i__ + i__ * b_dim1] = ei;
#line 268 "dlaein.f"
		    i__2 = *n;
#line 268 "dlaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 269 "dlaein.f"
			temp = b[i__ + 1 + j * b_dim1];
#line 270 "dlaein.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - x * 
				temp;
#line 271 "dlaein.f"
			b[i__ + j * b_dim1] = temp;
#line 272 "dlaein.f"
/* L40: */
#line 272 "dlaein.f"
		    }
#line 273 "dlaein.f"
		} else {

/*                 Eliminate without interchange. */

#line 277 "dlaein.f"
		    if (b[i__ + i__ * b_dim1] == 0.) {
#line 277 "dlaein.f"
			b[i__ + i__ * b_dim1] = *eps3;
#line 277 "dlaein.f"
		    }
#line 279 "dlaein.f"
		    x = ei / b[i__ + i__ * b_dim1];
#line 280 "dlaein.f"
		    if (x != 0.) {
#line 281 "dlaein.f"
			i__2 = *n;
#line 281 "dlaein.f"
			for (j = i__ + 1; j <= i__2; ++j) {
#line 282 "dlaein.f"
			    b[i__ + 1 + j * b_dim1] -= x * b[i__ + j * b_dim1]
				    ;
#line 283 "dlaein.f"
/* L50: */
#line 283 "dlaein.f"
			}
#line 284 "dlaein.f"
		    }
#line 285 "dlaein.f"
		}
#line 286 "dlaein.f"
/* L60: */
#line 286 "dlaein.f"
	    }
#line 287 "dlaein.f"
	    if (b[*n + *n * b_dim1] == 0.) {
#line 287 "dlaein.f"
		b[*n + *n * b_dim1] = *eps3;
#line 287 "dlaein.f"
	    }

#line 290 "dlaein.f"
	    *(unsigned char *)trans = 'N';

#line 292 "dlaein.f"
	} else {

/*           UL decomposition with partial pivoting of B, replacing zero */
/*           pivots by EPS3. */

#line 297 "dlaein.f"
	    for (j = *n; j >= 2; --j) {
#line 298 "dlaein.f"
		ej = h__[j + (j - 1) * h_dim1];
#line 299 "dlaein.f"
		if ((d__1 = b[j + j * b_dim1], abs(d__1)) < abs(ej)) {

/*                 Interchange columns and eliminate. */

#line 303 "dlaein.f"
		    x = b[j + j * b_dim1] / ej;
#line 304 "dlaein.f"
		    b[j + j * b_dim1] = ej;
#line 305 "dlaein.f"
		    i__1 = j - 1;
#line 305 "dlaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "dlaein.f"
			temp = b[i__ + (j - 1) * b_dim1];
#line 307 "dlaein.f"
			b[i__ + (j - 1) * b_dim1] = b[i__ + j * b_dim1] - x * 
				temp;
#line 308 "dlaein.f"
			b[i__ + j * b_dim1] = temp;
#line 309 "dlaein.f"
/* L70: */
#line 309 "dlaein.f"
		    }
#line 310 "dlaein.f"
		} else {

/*                 Eliminate without interchange. */

#line 314 "dlaein.f"
		    if (b[j + j * b_dim1] == 0.) {
#line 314 "dlaein.f"
			b[j + j * b_dim1] = *eps3;
#line 314 "dlaein.f"
		    }
#line 316 "dlaein.f"
		    x = ej / b[j + j * b_dim1];
#line 317 "dlaein.f"
		    if (x != 0.) {
#line 318 "dlaein.f"
			i__1 = j - 1;
#line 318 "dlaein.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 319 "dlaein.f"
			    b[i__ + (j - 1) * b_dim1] -= x * b[i__ + j * 
				    b_dim1];
#line 320 "dlaein.f"
/* L80: */
#line 320 "dlaein.f"
			}
#line 321 "dlaein.f"
		    }
#line 322 "dlaein.f"
		}
#line 323 "dlaein.f"
/* L90: */
#line 323 "dlaein.f"
	    }
#line 324 "dlaein.f"
	    if (b[b_dim1 + 1] == 0.) {
#line 324 "dlaein.f"
		b[b_dim1 + 1] = *eps3;
#line 324 "dlaein.f"
	    }

#line 327 "dlaein.f"
	    *(unsigned char *)trans = 'T';

#line 329 "dlaein.f"
	}

#line 331 "dlaein.f"
	*(unsigned char *)normin = 'N';
#line 332 "dlaein.f"
	i__1 = *n;
#line 332 "dlaein.f"
	for (its = 1; its <= i__1; ++its) {

/*           Solve U*x = scale*v for a right eigenvector */
/*             or U**T*x = scale*v for a left eigenvector, */
/*           overwriting x on v. */

#line 338 "dlaein.f"
	    dlatrs_("Upper", trans, "Nonunit", normin, n, &b[b_offset], ldb, &
		    vr[1], &scale, &work[1], &ierr, (ftnlen)5, (ftnlen)1, (
		    ftnlen)7, (ftnlen)1);
#line 340 "dlaein.f"
	    *(unsigned char *)normin = 'Y';

/*           Test for sufficient growth in the norm of v. */

#line 344 "dlaein.f"
	    vnorm = dasum_(n, &vr[1], &c__1);
#line 345 "dlaein.f"
	    if (vnorm >= growto * scale) {
#line 345 "dlaein.f"
		goto L120;
#line 345 "dlaein.f"
	    }

/*           Choose new orthogonal starting vector and try again. */

#line 350 "dlaein.f"
	    temp = *eps3 / (rootn + 1.);
#line 351 "dlaein.f"
	    vr[1] = *eps3;
#line 352 "dlaein.f"
	    i__2 = *n;
#line 352 "dlaein.f"
	    for (i__ = 2; i__ <= i__2; ++i__) {
#line 353 "dlaein.f"
		vr[i__] = temp;
#line 354 "dlaein.f"
/* L100: */
#line 354 "dlaein.f"
	    }
#line 355 "dlaein.f"
	    vr[*n - its + 1] -= *eps3 * rootn;
#line 356 "dlaein.f"
/* L110: */
#line 356 "dlaein.f"
	}

/*        Failure to find eigenvector in N iterations. */

#line 360 "dlaein.f"
	*info = 1;

#line 362 "dlaein.f"
L120:

/*        Normalize eigenvector. */

#line 366 "dlaein.f"
	i__ = idamax_(n, &vr[1], &c__1);
#line 367 "dlaein.f"
	d__2 = 1. / (d__1 = vr[i__], abs(d__1));
#line 367 "dlaein.f"
	dscal_(n, &d__2, &vr[1], &c__1);
#line 368 "dlaein.f"
    } else {

/*        Complex eigenvalue. */

#line 372 "dlaein.f"
	if (*noinit) {

/*           Set initial vector. */

#line 376 "dlaein.f"
	    i__1 = *n;
#line 376 "dlaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "dlaein.f"
		vr[i__] = *eps3;
#line 378 "dlaein.f"
		vi[i__] = 0.;
#line 379 "dlaein.f"
/* L130: */
#line 379 "dlaein.f"
	    }
#line 380 "dlaein.f"
	} else {

/*           Scale supplied initial vector. */

#line 384 "dlaein.f"
	    d__1 = dnrm2_(n, &vr[1], &c__1);
#line 384 "dlaein.f"
	    d__2 = dnrm2_(n, &vi[1], &c__1);
#line 384 "dlaein.f"
	    norm = dlapy2_(&d__1, &d__2);
#line 385 "dlaein.f"
	    rec = *eps3 * rootn / max(norm,nrmsml);
#line 386 "dlaein.f"
	    dscal_(n, &rec, &vr[1], &c__1);
#line 387 "dlaein.f"
	    dscal_(n, &rec, &vi[1], &c__1);
#line 388 "dlaein.f"
	}

#line 390 "dlaein.f"
	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacing zero */
/*           pivots by EPS3. */

/*           The imaginary part of the (i,j)-th element of U is stored in */
/*           B(j+1,i). */

#line 398 "dlaein.f"
	    b[b_dim1 + 2] = -(*wi);
#line 399 "dlaein.f"
	    i__1 = *n;
#line 399 "dlaein.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 400 "dlaein.f"
		b[i__ + 1 + b_dim1] = 0.;
#line 401 "dlaein.f"
/* L140: */
#line 401 "dlaein.f"
	    }

#line 403 "dlaein.f"
	    i__1 = *n - 1;
#line 403 "dlaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 404 "dlaein.f"
		absbii = dlapy2_(&b[i__ + i__ * b_dim1], &b[i__ + 1 + i__ * 
			b_dim1]);
#line 405 "dlaein.f"
		ei = h__[i__ + 1 + i__ * h_dim1];
#line 406 "dlaein.f"
		if (absbii < abs(ei)) {

/*                 Interchange rows and eliminate. */

#line 410 "dlaein.f"
		    xr = b[i__ + i__ * b_dim1] / ei;
#line 411 "dlaein.f"
		    xi = b[i__ + 1 + i__ * b_dim1] / ei;
#line 412 "dlaein.f"
		    b[i__ + i__ * b_dim1] = ei;
#line 413 "dlaein.f"
		    b[i__ + 1 + i__ * b_dim1] = 0.;
#line 414 "dlaein.f"
		    i__2 = *n;
#line 414 "dlaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 415 "dlaein.f"
			temp = b[i__ + 1 + j * b_dim1];
#line 416 "dlaein.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - xr * 
				temp;
#line 417 "dlaein.f"
			b[j + 1 + (i__ + 1) * b_dim1] = b[j + 1 + i__ * 
				b_dim1] - xi * temp;
#line 418 "dlaein.f"
			b[i__ + j * b_dim1] = temp;
#line 419 "dlaein.f"
			b[j + 1 + i__ * b_dim1] = 0.;
#line 420 "dlaein.f"
/* L150: */
#line 420 "dlaein.f"
		    }
#line 421 "dlaein.f"
		    b[i__ + 2 + i__ * b_dim1] = -(*wi);
#line 422 "dlaein.f"
		    b[i__ + 1 + (i__ + 1) * b_dim1] -= xi * *wi;
#line 423 "dlaein.f"
		    b[i__ + 2 + (i__ + 1) * b_dim1] += xr * *wi;
#line 424 "dlaein.f"
		} else {

/*                 Eliminate without interchanging rows. */

#line 428 "dlaein.f"
		    if (absbii == 0.) {
#line 429 "dlaein.f"
			b[i__ + i__ * b_dim1] = *eps3;
#line 430 "dlaein.f"
			b[i__ + 1 + i__ * b_dim1] = 0.;
#line 431 "dlaein.f"
			absbii = *eps3;
#line 432 "dlaein.f"
		    }
#line 433 "dlaein.f"
		    ei = ei / absbii / absbii;
#line 434 "dlaein.f"
		    xr = b[i__ + i__ * b_dim1] * ei;
#line 435 "dlaein.f"
		    xi = -b[i__ + 1 + i__ * b_dim1] * ei;
#line 436 "dlaein.f"
		    i__2 = *n;
#line 436 "dlaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 437 "dlaein.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + 1 + j * b_dim1] - 
				xr * b[i__ + j * b_dim1] + xi * b[j + 1 + i__ 
				* b_dim1];
#line 439 "dlaein.f"
			b[j + 1 + (i__ + 1) * b_dim1] = -xr * b[j + 1 + i__ * 
				b_dim1] - xi * b[i__ + j * b_dim1];
#line 440 "dlaein.f"
/* L160: */
#line 440 "dlaein.f"
		    }
#line 441 "dlaein.f"
		    b[i__ + 2 + (i__ + 1) * b_dim1] -= *wi;
#line 442 "dlaein.f"
		}

/*              Compute 1-norm of offdiagonal elements of i-th row. */

#line 446 "dlaein.f"
		i__2 = *n - i__;
#line 446 "dlaein.f"
		i__3 = *n - i__;
#line 446 "dlaein.f"
		work[i__] = dasum_(&i__2, &b[i__ + (i__ + 1) * b_dim1], ldb) 
			+ dasum_(&i__3, &b[i__ + 2 + i__ * b_dim1], &c__1);
#line 448 "dlaein.f"
/* L170: */
#line 448 "dlaein.f"
	    }
#line 449 "dlaein.f"
	    if (b[*n + *n * b_dim1] == 0. && b[*n + 1 + *n * b_dim1] == 0.) {
#line 449 "dlaein.f"
		b[*n + *n * b_dim1] = *eps3;
#line 449 "dlaein.f"
	    }
#line 451 "dlaein.f"
	    work[*n] = 0.;

#line 453 "dlaein.f"
	    i1 = *n;
#line 454 "dlaein.f"
	    i2 = 1;
#line 455 "dlaein.f"
	    i3 = -1;
#line 456 "dlaein.f"
	} else {

/*           UL decomposition with partial pivoting of conjg(B), */
/*           replacing zero pivots by EPS3. */

/*           The imaginary part of the (i,j)-th element of U is stored in */
/*           B(j+1,i). */

#line 464 "dlaein.f"
	    b[*n + 1 + *n * b_dim1] = *wi;
#line 465 "dlaein.f"
	    i__1 = *n - 1;
#line 465 "dlaein.f"
	    for (j = 1; j <= i__1; ++j) {
#line 466 "dlaein.f"
		b[*n + 1 + j * b_dim1] = 0.;
#line 467 "dlaein.f"
/* L180: */
#line 467 "dlaein.f"
	    }

#line 469 "dlaein.f"
	    for (j = *n; j >= 2; --j) {
#line 470 "dlaein.f"
		ej = h__[j + (j - 1) * h_dim1];
#line 471 "dlaein.f"
		absbjj = dlapy2_(&b[j + j * b_dim1], &b[j + 1 + j * b_dim1]);
#line 472 "dlaein.f"
		if (absbjj < abs(ej)) {

/*                 Interchange columns and eliminate */

#line 476 "dlaein.f"
		    xr = b[j + j * b_dim1] / ej;
#line 477 "dlaein.f"
		    xi = b[j + 1 + j * b_dim1] / ej;
#line 478 "dlaein.f"
		    b[j + j * b_dim1] = ej;
#line 479 "dlaein.f"
		    b[j + 1 + j * b_dim1] = 0.;
#line 480 "dlaein.f"
		    i__1 = j - 1;
#line 480 "dlaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 481 "dlaein.f"
			temp = b[i__ + (j - 1) * b_dim1];
#line 482 "dlaein.f"
			b[i__ + (j - 1) * b_dim1] = b[i__ + j * b_dim1] - xr *
				 temp;
#line 483 "dlaein.f"
			b[j + i__ * b_dim1] = b[j + 1 + i__ * b_dim1] - xi * 
				temp;
#line 484 "dlaein.f"
			b[i__ + j * b_dim1] = temp;
#line 485 "dlaein.f"
			b[j + 1 + i__ * b_dim1] = 0.;
#line 486 "dlaein.f"
/* L190: */
#line 486 "dlaein.f"
		    }
#line 487 "dlaein.f"
		    b[j + 1 + (j - 1) * b_dim1] = *wi;
#line 488 "dlaein.f"
		    b[j - 1 + (j - 1) * b_dim1] += xi * *wi;
#line 489 "dlaein.f"
		    b[j + (j - 1) * b_dim1] -= xr * *wi;
#line 490 "dlaein.f"
		} else {

/*                 Eliminate without interchange. */

#line 494 "dlaein.f"
		    if (absbjj == 0.) {
#line 495 "dlaein.f"
			b[j + j * b_dim1] = *eps3;
#line 496 "dlaein.f"
			b[j + 1 + j * b_dim1] = 0.;
#line 497 "dlaein.f"
			absbjj = *eps3;
#line 498 "dlaein.f"
		    }
#line 499 "dlaein.f"
		    ej = ej / absbjj / absbjj;
#line 500 "dlaein.f"
		    xr = b[j + j * b_dim1] * ej;
#line 501 "dlaein.f"
		    xi = -b[j + 1 + j * b_dim1] * ej;
#line 502 "dlaein.f"
		    i__1 = j - 1;
#line 502 "dlaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 503 "dlaein.f"
			b[i__ + (j - 1) * b_dim1] = b[i__ + (j - 1) * b_dim1] 
				- xr * b[i__ + j * b_dim1] + xi * b[j + 1 + 
				i__ * b_dim1];
#line 505 "dlaein.f"
			b[j + i__ * b_dim1] = -xr * b[j + 1 + i__ * b_dim1] - 
				xi * b[i__ + j * b_dim1];
#line 506 "dlaein.f"
/* L200: */
#line 506 "dlaein.f"
		    }
#line 507 "dlaein.f"
		    b[j + (j - 1) * b_dim1] += *wi;
#line 508 "dlaein.f"
		}

/*              Compute 1-norm of offdiagonal elements of j-th column. */

#line 512 "dlaein.f"
		i__1 = j - 1;
#line 512 "dlaein.f"
		i__2 = j - 1;
#line 512 "dlaein.f"
		work[j] = dasum_(&i__1, &b[j * b_dim1 + 1], &c__1) + dasum_(&
			i__2, &b[j + 1 + b_dim1], ldb);
#line 514 "dlaein.f"
/* L210: */
#line 514 "dlaein.f"
	    }
#line 515 "dlaein.f"
	    if (b[b_dim1 + 1] == 0. && b[b_dim1 + 2] == 0.) {
#line 515 "dlaein.f"
		b[b_dim1 + 1] = *eps3;
#line 515 "dlaein.f"
	    }
#line 517 "dlaein.f"
	    work[1] = 0.;

#line 519 "dlaein.f"
	    i1 = 1;
#line 520 "dlaein.f"
	    i2 = *n;
#line 521 "dlaein.f"
	    i3 = 1;
#line 522 "dlaein.f"
	}

#line 524 "dlaein.f"
	i__1 = *n;
#line 524 "dlaein.f"
	for (its = 1; its <= i__1; ++its) {
#line 525 "dlaein.f"
	    scale = 1.;
#line 526 "dlaein.f"
	    vmax = 1.;
#line 527 "dlaein.f"
	    vcrit = *bignum;

/*           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector, */
/*             or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector, */
/*           overwriting (xr,xi) on (vr,vi). */

#line 533 "dlaein.f"
	    i__2 = i2;
#line 533 "dlaein.f"
	    i__3 = i3;
#line 533 "dlaein.f"
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {

#line 535 "dlaein.f"
		if (work[i__] > vcrit) {
#line 536 "dlaein.f"
		    rec = 1. / vmax;
#line 537 "dlaein.f"
		    dscal_(n, &rec, &vr[1], &c__1);
#line 538 "dlaein.f"
		    dscal_(n, &rec, &vi[1], &c__1);
#line 539 "dlaein.f"
		    scale *= rec;
#line 540 "dlaein.f"
		    vmax = 1.;
#line 541 "dlaein.f"
		    vcrit = *bignum;
#line 542 "dlaein.f"
		}

#line 544 "dlaein.f"
		xr = vr[i__];
#line 545 "dlaein.f"
		xi = vi[i__];
#line 546 "dlaein.f"
		if (*rightv) {
#line 547 "dlaein.f"
		    i__4 = *n;
#line 547 "dlaein.f"
		    for (j = i__ + 1; j <= i__4; ++j) {
#line 548 "dlaein.f"
			xr = xr - b[i__ + j * b_dim1] * vr[j] + b[j + 1 + i__ 
				* b_dim1] * vi[j];
#line 549 "dlaein.f"
			xi = xi - b[i__ + j * b_dim1] * vi[j] - b[j + 1 + i__ 
				* b_dim1] * vr[j];
#line 550 "dlaein.f"
/* L220: */
#line 550 "dlaein.f"
		    }
#line 551 "dlaein.f"
		} else {
#line 552 "dlaein.f"
		    i__4 = i__ - 1;
#line 552 "dlaein.f"
		    for (j = 1; j <= i__4; ++j) {
#line 553 "dlaein.f"
			xr = xr - b[j + i__ * b_dim1] * vr[j] + b[i__ + 1 + j 
				* b_dim1] * vi[j];
#line 554 "dlaein.f"
			xi = xi - b[j + i__ * b_dim1] * vi[j] - b[i__ + 1 + j 
				* b_dim1] * vr[j];
#line 555 "dlaein.f"
/* L230: */
#line 555 "dlaein.f"
		    }
#line 556 "dlaein.f"
		}

#line 558 "dlaein.f"
		w = (d__1 = b[i__ + i__ * b_dim1], abs(d__1)) + (d__2 = b[i__ 
			+ 1 + i__ * b_dim1], abs(d__2));
#line 559 "dlaein.f"
		if (w > *smlnum) {
#line 560 "dlaein.f"
		    if (w < 1.) {
#line 561 "dlaein.f"
			w1 = abs(xr) + abs(xi);
#line 562 "dlaein.f"
			if (w1 > w * *bignum) {
#line 563 "dlaein.f"
			    rec = 1. / w1;
#line 564 "dlaein.f"
			    dscal_(n, &rec, &vr[1], &c__1);
#line 565 "dlaein.f"
			    dscal_(n, &rec, &vi[1], &c__1);
#line 566 "dlaein.f"
			    xr = vr[i__];
#line 567 "dlaein.f"
			    xi = vi[i__];
#line 568 "dlaein.f"
			    scale *= rec;
#line 569 "dlaein.f"
			    vmax *= rec;
#line 570 "dlaein.f"
			}
#line 571 "dlaein.f"
		    }

/*                 Divide by diagonal element of B. */

#line 575 "dlaein.f"
		    dladiv_(&xr, &xi, &b[i__ + i__ * b_dim1], &b[i__ + 1 + 
			    i__ * b_dim1], &vr[i__], &vi[i__]);
/* Computing MAX */
#line 577 "dlaein.f"
		    d__3 = (d__1 = vr[i__], abs(d__1)) + (d__2 = vi[i__], abs(
			    d__2));
#line 577 "dlaein.f"
		    vmax = max(d__3,vmax);
#line 578 "dlaein.f"
		    vcrit = *bignum / vmax;
#line 579 "dlaein.f"
		} else {
#line 580 "dlaein.f"
		    i__4 = *n;
#line 580 "dlaein.f"
		    for (j = 1; j <= i__4; ++j) {
#line 581 "dlaein.f"
			vr[j] = 0.;
#line 582 "dlaein.f"
			vi[j] = 0.;
#line 583 "dlaein.f"
/* L240: */
#line 583 "dlaein.f"
		    }
#line 584 "dlaein.f"
		    vr[i__] = 1.;
#line 585 "dlaein.f"
		    vi[i__] = 1.;
#line 586 "dlaein.f"
		    scale = 0.;
#line 587 "dlaein.f"
		    vmax = 1.;
#line 588 "dlaein.f"
		    vcrit = *bignum;
#line 589 "dlaein.f"
		}
#line 590 "dlaein.f"
/* L250: */
#line 590 "dlaein.f"
	    }

/*           Test for sufficient growth in the norm of (VR,VI). */

#line 594 "dlaein.f"
	    vnorm = dasum_(n, &vr[1], &c__1) + dasum_(n, &vi[1], &c__1);
#line 595 "dlaein.f"
	    if (vnorm >= growto * scale) {
#line 595 "dlaein.f"
		goto L280;
#line 595 "dlaein.f"
	    }

/*           Choose a new orthogonal starting vector and try again. */

#line 600 "dlaein.f"
	    y = *eps3 / (rootn + 1.);
#line 601 "dlaein.f"
	    vr[1] = *eps3;
#line 602 "dlaein.f"
	    vi[1] = 0.;

#line 604 "dlaein.f"
	    i__3 = *n;
#line 604 "dlaein.f"
	    for (i__ = 2; i__ <= i__3; ++i__) {
#line 605 "dlaein.f"
		vr[i__] = y;
#line 606 "dlaein.f"
		vi[i__] = 0.;
#line 607 "dlaein.f"
/* L260: */
#line 607 "dlaein.f"
	    }
#line 608 "dlaein.f"
	    vr[*n - its + 1] -= *eps3 * rootn;
#line 609 "dlaein.f"
/* L270: */
#line 609 "dlaein.f"
	}

/*        Failure to find eigenvector in N iterations */

#line 613 "dlaein.f"
	*info = 1;

#line 615 "dlaein.f"
L280:

/*        Normalize eigenvector. */

#line 619 "dlaein.f"
	vnorm = 0.;
#line 620 "dlaein.f"
	i__1 = *n;
#line 620 "dlaein.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 621 "dlaein.f"
	    d__3 = vnorm, d__4 = (d__1 = vr[i__], abs(d__1)) + (d__2 = vi[i__]
		    , abs(d__2));
#line 621 "dlaein.f"
	    vnorm = max(d__3,d__4);
#line 622 "dlaein.f"
/* L290: */
#line 622 "dlaein.f"
	}
#line 623 "dlaein.f"
	d__1 = 1. / vnorm;
#line 623 "dlaein.f"
	dscal_(n, &d__1, &vr[1], &c__1);
#line 624 "dlaein.f"
	d__1 = 1. / vnorm;
#line 624 "dlaein.f"
	dscal_(n, &d__1, &vi[1], &c__1);

#line 626 "dlaein.f"
    }

#line 628 "dlaein.f"
    return 0;

/*     End of DLAEIN */

} /* dlaein_ */

