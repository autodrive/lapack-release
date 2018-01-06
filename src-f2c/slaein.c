#line 1 "slaein.f"
/* slaein.f -- translated by f2c (version 20100827).
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

#line 1 "slaein.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse 
iteration. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B, */
/*                          LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            NOINIT, RIGHTV */
/*       INTEGER            INFO, LDB, LDH, N */
/*       REAL               BIGNUM, EPS3, SMLNUM, WI, WR */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               B( LDB, * ), H( LDH, * ), VI( * ), VR( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAEIN uses inverse iteration to find a right or left eigenvector */
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
/* >          H is REAL array, dimension (LDH,N) */
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
/* >          WR is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* >          WI is REAL */
/* >          The real and imaginary parts of the eigenvalue of H whose */
/* >          corresponding right or left eigenvector is to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in,out] VI */
/* > \verbatim */
/* >          VI is REAL array, dimension (N) */
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
/* >          B is REAL array, dimension (LDB,N) */
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
/* >          WORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] EPS3 */
/* > \verbatim */
/* >          EPS3 is REAL */
/* >          A small machine-dependent value which is used to perturb */
/* >          close eigenvalues, and to replace zero pivots. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLNUM */
/* > \verbatim */
/* >          SMLNUM is REAL */
/* >          A machine-dependent value close to the underflow threshold. */
/* > \endverbatim */
/* > */
/* > \param[in] BIGNUM */
/* > \verbatim */
/* >          BIGNUM is REAL */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaein_(logical *rightv, logical *noinit, integer *n, 
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char trans[1];
    static doublereal vcrit;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal rootn, vnorm;
    extern doublereal slapy2_(doublereal *, doublereal *);
    static doublereal absbii, absbjj;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int sladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static char normin[1];
    static doublereal nrmsml;
    extern /* Subroutine */ int slatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal growto;


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

#line 216 "slaein.f"
    /* Parameter adjustments */
#line 216 "slaein.f"
    h_dim1 = *ldh;
#line 216 "slaein.f"
    h_offset = 1 + h_dim1;
#line 216 "slaein.f"
    h__ -= h_offset;
#line 216 "slaein.f"
    --vr;
#line 216 "slaein.f"
    --vi;
#line 216 "slaein.f"
    b_dim1 = *ldb;
#line 216 "slaein.f"
    b_offset = 1 + b_dim1;
#line 216 "slaein.f"
    b -= b_offset;
#line 216 "slaein.f"
    --work;
#line 216 "slaein.f"

#line 216 "slaein.f"
    /* Function Body */
#line 216 "slaein.f"
    *info = 0;

/*     GROWTO is the threshold used in the acceptance test for an */
/*     eigenvector. */

#line 221 "slaein.f"
    rootn = sqrt((doublereal) (*n));
#line 222 "slaein.f"
    growto = .1 / rootn;
/* Computing MAX */
#line 223 "slaein.f"
    d__1 = 1., d__2 = *eps3 * rootn;
#line 223 "slaein.f"
    nrmsml = max(d__1,d__2) * *smlnum;

/*     Form B = H - (WR,WI)*I (except that the subdiagonal elements and */
/*     the imaginary parts of the diagonal elements are not stored). */

#line 228 "slaein.f"
    i__1 = *n;
#line 228 "slaein.f"
    for (j = 1; j <= i__1; ++j) {
#line 229 "slaein.f"
	i__2 = j - 1;
#line 229 "slaein.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 230 "slaein.f"
	    b[i__ + j * b_dim1] = h__[i__ + j * h_dim1];
#line 231 "slaein.f"
/* L10: */
#line 231 "slaein.f"
	}
#line 232 "slaein.f"
	b[j + j * b_dim1] = h__[j + j * h_dim1] - *wr;
#line 233 "slaein.f"
/* L20: */
#line 233 "slaein.f"
    }

#line 235 "slaein.f"
    if (*wi == 0.) {

/*        Real eigenvalue. */

#line 239 "slaein.f"
	if (*noinit) {

/*           Set initial vector. */

#line 243 "slaein.f"
	    i__1 = *n;
#line 243 "slaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "slaein.f"
		vr[i__] = *eps3;
#line 245 "slaein.f"
/* L30: */
#line 245 "slaein.f"
	    }
#line 246 "slaein.f"
	} else {

/*           Scale supplied initial vector. */

#line 250 "slaein.f"
	    vnorm = snrm2_(n, &vr[1], &c__1);
#line 251 "slaein.f"
	    d__1 = *eps3 * rootn / max(vnorm,nrmsml);
#line 251 "slaein.f"
	    sscal_(n, &d__1, &vr[1], &c__1);
#line 253 "slaein.f"
	}

#line 255 "slaein.f"
	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacing zero */
/*           pivots by EPS3. */

#line 260 "slaein.f"
	    i__1 = *n - 1;
#line 260 "slaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 261 "slaein.f"
		ei = h__[i__ + 1 + i__ * h_dim1];
#line 262 "slaein.f"
		if ((d__1 = b[i__ + i__ * b_dim1], abs(d__1)) < abs(ei)) {

/*                 Interchange rows and eliminate. */

#line 266 "slaein.f"
		    x = b[i__ + i__ * b_dim1] / ei;
#line 267 "slaein.f"
		    b[i__ + i__ * b_dim1] = ei;
#line 268 "slaein.f"
		    i__2 = *n;
#line 268 "slaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 269 "slaein.f"
			temp = b[i__ + 1 + j * b_dim1];
#line 270 "slaein.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - x * 
				temp;
#line 271 "slaein.f"
			b[i__ + j * b_dim1] = temp;
#line 272 "slaein.f"
/* L40: */
#line 272 "slaein.f"
		    }
#line 273 "slaein.f"
		} else {

/*                 Eliminate without interchange. */

#line 277 "slaein.f"
		    if (b[i__ + i__ * b_dim1] == 0.) {
#line 277 "slaein.f"
			b[i__ + i__ * b_dim1] = *eps3;
#line 277 "slaein.f"
		    }
#line 279 "slaein.f"
		    x = ei / b[i__ + i__ * b_dim1];
#line 280 "slaein.f"
		    if (x != 0.) {
#line 281 "slaein.f"
			i__2 = *n;
#line 281 "slaein.f"
			for (j = i__ + 1; j <= i__2; ++j) {
#line 282 "slaein.f"
			    b[i__ + 1 + j * b_dim1] -= x * b[i__ + j * b_dim1]
				    ;
#line 283 "slaein.f"
/* L50: */
#line 283 "slaein.f"
			}
#line 284 "slaein.f"
		    }
#line 285 "slaein.f"
		}
#line 286 "slaein.f"
/* L60: */
#line 286 "slaein.f"
	    }
#line 287 "slaein.f"
	    if (b[*n + *n * b_dim1] == 0.) {
#line 287 "slaein.f"
		b[*n + *n * b_dim1] = *eps3;
#line 287 "slaein.f"
	    }

#line 290 "slaein.f"
	    *(unsigned char *)trans = 'N';

#line 292 "slaein.f"
	} else {

/*           UL decomposition with partial pivoting of B, replacing zero */
/*           pivots by EPS3. */

#line 297 "slaein.f"
	    for (j = *n; j >= 2; --j) {
#line 298 "slaein.f"
		ej = h__[j + (j - 1) * h_dim1];
#line 299 "slaein.f"
		if ((d__1 = b[j + j * b_dim1], abs(d__1)) < abs(ej)) {

/*                 Interchange columns and eliminate. */

#line 303 "slaein.f"
		    x = b[j + j * b_dim1] / ej;
#line 304 "slaein.f"
		    b[j + j * b_dim1] = ej;
#line 305 "slaein.f"
		    i__1 = j - 1;
#line 305 "slaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "slaein.f"
			temp = b[i__ + (j - 1) * b_dim1];
#line 307 "slaein.f"
			b[i__ + (j - 1) * b_dim1] = b[i__ + j * b_dim1] - x * 
				temp;
#line 308 "slaein.f"
			b[i__ + j * b_dim1] = temp;
#line 309 "slaein.f"
/* L70: */
#line 309 "slaein.f"
		    }
#line 310 "slaein.f"
		} else {

/*                 Eliminate without interchange. */

#line 314 "slaein.f"
		    if (b[j + j * b_dim1] == 0.) {
#line 314 "slaein.f"
			b[j + j * b_dim1] = *eps3;
#line 314 "slaein.f"
		    }
#line 316 "slaein.f"
		    x = ej / b[j + j * b_dim1];
#line 317 "slaein.f"
		    if (x != 0.) {
#line 318 "slaein.f"
			i__1 = j - 1;
#line 318 "slaein.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 319 "slaein.f"
			    b[i__ + (j - 1) * b_dim1] -= x * b[i__ + j * 
				    b_dim1];
#line 320 "slaein.f"
/* L80: */
#line 320 "slaein.f"
			}
#line 321 "slaein.f"
		    }
#line 322 "slaein.f"
		}
#line 323 "slaein.f"
/* L90: */
#line 323 "slaein.f"
	    }
#line 324 "slaein.f"
	    if (b[b_dim1 + 1] == 0.) {
#line 324 "slaein.f"
		b[b_dim1 + 1] = *eps3;
#line 324 "slaein.f"
	    }

#line 327 "slaein.f"
	    *(unsigned char *)trans = 'T';

#line 329 "slaein.f"
	}

#line 331 "slaein.f"
	*(unsigned char *)normin = 'N';
#line 332 "slaein.f"
	i__1 = *n;
#line 332 "slaein.f"
	for (its = 1; its <= i__1; ++its) {

/*           Solve U*x = scale*v for a right eigenvector */
/*             or U**T*x = scale*v for a left eigenvector, */
/*           overwriting x on v. */

#line 338 "slaein.f"
	    slatrs_("Upper", trans, "Nonunit", normin, n, &b[b_offset], ldb, &
		    vr[1], &scale, &work[1], &ierr, (ftnlen)5, (ftnlen)1, (
		    ftnlen)7, (ftnlen)1);
#line 340 "slaein.f"
	    *(unsigned char *)normin = 'Y';

/*           Test for sufficient growth in the norm of v. */

#line 344 "slaein.f"
	    vnorm = sasum_(n, &vr[1], &c__1);
#line 345 "slaein.f"
	    if (vnorm >= growto * scale) {
#line 345 "slaein.f"
		goto L120;
#line 345 "slaein.f"
	    }

/*           Choose new orthogonal starting vector and try again. */

#line 350 "slaein.f"
	    temp = *eps3 / (rootn + 1.);
#line 351 "slaein.f"
	    vr[1] = *eps3;
#line 352 "slaein.f"
	    i__2 = *n;
#line 352 "slaein.f"
	    for (i__ = 2; i__ <= i__2; ++i__) {
#line 353 "slaein.f"
		vr[i__] = temp;
#line 354 "slaein.f"
/* L100: */
#line 354 "slaein.f"
	    }
#line 355 "slaein.f"
	    vr[*n - its + 1] -= *eps3 * rootn;
#line 356 "slaein.f"
/* L110: */
#line 356 "slaein.f"
	}

/*        Failure to find eigenvector in N iterations. */

#line 360 "slaein.f"
	*info = 1;

#line 362 "slaein.f"
L120:

/*        Normalize eigenvector. */

#line 366 "slaein.f"
	i__ = isamax_(n, &vr[1], &c__1);
#line 367 "slaein.f"
	d__2 = 1. / (d__1 = vr[i__], abs(d__1));
#line 367 "slaein.f"
	sscal_(n, &d__2, &vr[1], &c__1);
#line 368 "slaein.f"
    } else {

/*        Complex eigenvalue. */

#line 372 "slaein.f"
	if (*noinit) {

/*           Set initial vector. */

#line 376 "slaein.f"
	    i__1 = *n;
#line 376 "slaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "slaein.f"
		vr[i__] = *eps3;
#line 378 "slaein.f"
		vi[i__] = 0.;
#line 379 "slaein.f"
/* L130: */
#line 379 "slaein.f"
	    }
#line 380 "slaein.f"
	} else {

/*           Scale supplied initial vector. */

#line 384 "slaein.f"
	    d__1 = snrm2_(n, &vr[1], &c__1);
#line 384 "slaein.f"
	    d__2 = snrm2_(n, &vi[1], &c__1);
#line 384 "slaein.f"
	    norm = slapy2_(&d__1, &d__2);
#line 385 "slaein.f"
	    rec = *eps3 * rootn / max(norm,nrmsml);
#line 386 "slaein.f"
	    sscal_(n, &rec, &vr[1], &c__1);
#line 387 "slaein.f"
	    sscal_(n, &rec, &vi[1], &c__1);
#line 388 "slaein.f"
	}

#line 390 "slaein.f"
	if (*rightv) {

/*           LU decomposition with partial pivoting of B, replacing zero */
/*           pivots by EPS3. */

/*           The imaginary part of the (i,j)-th element of U is stored in */
/*           B(j+1,i). */

#line 398 "slaein.f"
	    b[b_dim1 + 2] = -(*wi);
#line 399 "slaein.f"
	    i__1 = *n;
#line 399 "slaein.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 400 "slaein.f"
		b[i__ + 1 + b_dim1] = 0.;
#line 401 "slaein.f"
/* L140: */
#line 401 "slaein.f"
	    }

#line 403 "slaein.f"
	    i__1 = *n - 1;
#line 403 "slaein.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 404 "slaein.f"
		absbii = slapy2_(&b[i__ + i__ * b_dim1], &b[i__ + 1 + i__ * 
			b_dim1]);
#line 405 "slaein.f"
		ei = h__[i__ + 1 + i__ * h_dim1];
#line 406 "slaein.f"
		if (absbii < abs(ei)) {

/*                 Interchange rows and eliminate. */

#line 410 "slaein.f"
		    xr = b[i__ + i__ * b_dim1] / ei;
#line 411 "slaein.f"
		    xi = b[i__ + 1 + i__ * b_dim1] / ei;
#line 412 "slaein.f"
		    b[i__ + i__ * b_dim1] = ei;
#line 413 "slaein.f"
		    b[i__ + 1 + i__ * b_dim1] = 0.;
#line 414 "slaein.f"
		    i__2 = *n;
#line 414 "slaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 415 "slaein.f"
			temp = b[i__ + 1 + j * b_dim1];
#line 416 "slaein.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - xr * 
				temp;
#line 417 "slaein.f"
			b[j + 1 + (i__ + 1) * b_dim1] = b[j + 1 + i__ * 
				b_dim1] - xi * temp;
#line 418 "slaein.f"
			b[i__ + j * b_dim1] = temp;
#line 419 "slaein.f"
			b[j + 1 + i__ * b_dim1] = 0.;
#line 420 "slaein.f"
/* L150: */
#line 420 "slaein.f"
		    }
#line 421 "slaein.f"
		    b[i__ + 2 + i__ * b_dim1] = -(*wi);
#line 422 "slaein.f"
		    b[i__ + 1 + (i__ + 1) * b_dim1] -= xi * *wi;
#line 423 "slaein.f"
		    b[i__ + 2 + (i__ + 1) * b_dim1] += xr * *wi;
#line 424 "slaein.f"
		} else {

/*                 Eliminate without interchanging rows. */

#line 428 "slaein.f"
		    if (absbii == 0.) {
#line 429 "slaein.f"
			b[i__ + i__ * b_dim1] = *eps3;
#line 430 "slaein.f"
			b[i__ + 1 + i__ * b_dim1] = 0.;
#line 431 "slaein.f"
			absbii = *eps3;
#line 432 "slaein.f"
		    }
#line 433 "slaein.f"
		    ei = ei / absbii / absbii;
#line 434 "slaein.f"
		    xr = b[i__ + i__ * b_dim1] * ei;
#line 435 "slaein.f"
		    xi = -b[i__ + 1 + i__ * b_dim1] * ei;
#line 436 "slaein.f"
		    i__2 = *n;
#line 436 "slaein.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 437 "slaein.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + 1 + j * b_dim1] - 
				xr * b[i__ + j * b_dim1] + xi * b[j + 1 + i__ 
				* b_dim1];
#line 439 "slaein.f"
			b[j + 1 + (i__ + 1) * b_dim1] = -xr * b[j + 1 + i__ * 
				b_dim1] - xi * b[i__ + j * b_dim1];
#line 440 "slaein.f"
/* L160: */
#line 440 "slaein.f"
		    }
#line 441 "slaein.f"
		    b[i__ + 2 + (i__ + 1) * b_dim1] -= *wi;
#line 442 "slaein.f"
		}

/*              Compute 1-norm of offdiagonal elements of i-th row. */

#line 446 "slaein.f"
		i__2 = *n - i__;
#line 446 "slaein.f"
		i__3 = *n - i__;
#line 446 "slaein.f"
		work[i__] = sasum_(&i__2, &b[i__ + (i__ + 1) * b_dim1], ldb) 
			+ sasum_(&i__3, &b[i__ + 2 + i__ * b_dim1], &c__1);
#line 448 "slaein.f"
/* L170: */
#line 448 "slaein.f"
	    }
#line 449 "slaein.f"
	    if (b[*n + *n * b_dim1] == 0. && b[*n + 1 + *n * b_dim1] == 0.) {
#line 449 "slaein.f"
		b[*n + *n * b_dim1] = *eps3;
#line 449 "slaein.f"
	    }
#line 451 "slaein.f"
	    work[*n] = 0.;

#line 453 "slaein.f"
	    i1 = *n;
#line 454 "slaein.f"
	    i2 = 1;
#line 455 "slaein.f"
	    i3 = -1;
#line 456 "slaein.f"
	} else {

/*           UL decomposition with partial pivoting of conjg(B), */
/*           replacing zero pivots by EPS3. */

/*           The imaginary part of the (i,j)-th element of U is stored in */
/*           B(j+1,i). */

#line 464 "slaein.f"
	    b[*n + 1 + *n * b_dim1] = *wi;
#line 465 "slaein.f"
	    i__1 = *n - 1;
#line 465 "slaein.f"
	    for (j = 1; j <= i__1; ++j) {
#line 466 "slaein.f"
		b[*n + 1 + j * b_dim1] = 0.;
#line 467 "slaein.f"
/* L180: */
#line 467 "slaein.f"
	    }

#line 469 "slaein.f"
	    for (j = *n; j >= 2; --j) {
#line 470 "slaein.f"
		ej = h__[j + (j - 1) * h_dim1];
#line 471 "slaein.f"
		absbjj = slapy2_(&b[j + j * b_dim1], &b[j + 1 + j * b_dim1]);
#line 472 "slaein.f"
		if (absbjj < abs(ej)) {

/*                 Interchange columns and eliminate */

#line 476 "slaein.f"
		    xr = b[j + j * b_dim1] / ej;
#line 477 "slaein.f"
		    xi = b[j + 1 + j * b_dim1] / ej;
#line 478 "slaein.f"
		    b[j + j * b_dim1] = ej;
#line 479 "slaein.f"
		    b[j + 1 + j * b_dim1] = 0.;
#line 480 "slaein.f"
		    i__1 = j - 1;
#line 480 "slaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 481 "slaein.f"
			temp = b[i__ + (j - 1) * b_dim1];
#line 482 "slaein.f"
			b[i__ + (j - 1) * b_dim1] = b[i__ + j * b_dim1] - xr *
				 temp;
#line 483 "slaein.f"
			b[j + i__ * b_dim1] = b[j + 1 + i__ * b_dim1] - xi * 
				temp;
#line 484 "slaein.f"
			b[i__ + j * b_dim1] = temp;
#line 485 "slaein.f"
			b[j + 1 + i__ * b_dim1] = 0.;
#line 486 "slaein.f"
/* L190: */
#line 486 "slaein.f"
		    }
#line 487 "slaein.f"
		    b[j + 1 + (j - 1) * b_dim1] = *wi;
#line 488 "slaein.f"
		    b[j - 1 + (j - 1) * b_dim1] += xi * *wi;
#line 489 "slaein.f"
		    b[j + (j - 1) * b_dim1] -= xr * *wi;
#line 490 "slaein.f"
		} else {

/*                 Eliminate without interchange. */

#line 494 "slaein.f"
		    if (absbjj == 0.) {
#line 495 "slaein.f"
			b[j + j * b_dim1] = *eps3;
#line 496 "slaein.f"
			b[j + 1 + j * b_dim1] = 0.;
#line 497 "slaein.f"
			absbjj = *eps3;
#line 498 "slaein.f"
		    }
#line 499 "slaein.f"
		    ej = ej / absbjj / absbjj;
#line 500 "slaein.f"
		    xr = b[j + j * b_dim1] * ej;
#line 501 "slaein.f"
		    xi = -b[j + 1 + j * b_dim1] * ej;
#line 502 "slaein.f"
		    i__1 = j - 1;
#line 502 "slaein.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 503 "slaein.f"
			b[i__ + (j - 1) * b_dim1] = b[i__ + (j - 1) * b_dim1] 
				- xr * b[i__ + j * b_dim1] + xi * b[j + 1 + 
				i__ * b_dim1];
#line 505 "slaein.f"
			b[j + i__ * b_dim1] = -xr * b[j + 1 + i__ * b_dim1] - 
				xi * b[i__ + j * b_dim1];
#line 506 "slaein.f"
/* L200: */
#line 506 "slaein.f"
		    }
#line 507 "slaein.f"
		    b[j + (j - 1) * b_dim1] += *wi;
#line 508 "slaein.f"
		}

/*              Compute 1-norm of offdiagonal elements of j-th column. */

#line 512 "slaein.f"
		i__1 = j - 1;
#line 512 "slaein.f"
		i__2 = j - 1;
#line 512 "slaein.f"
		work[j] = sasum_(&i__1, &b[j * b_dim1 + 1], &c__1) + sasum_(&
			i__2, &b[j + 1 + b_dim1], ldb);
#line 514 "slaein.f"
/* L210: */
#line 514 "slaein.f"
	    }
#line 515 "slaein.f"
	    if (b[b_dim1 + 1] == 0. && b[b_dim1 + 2] == 0.) {
#line 515 "slaein.f"
		b[b_dim1 + 1] = *eps3;
#line 515 "slaein.f"
	    }
#line 517 "slaein.f"
	    work[1] = 0.;

#line 519 "slaein.f"
	    i1 = 1;
#line 520 "slaein.f"
	    i2 = *n;
#line 521 "slaein.f"
	    i3 = 1;
#line 522 "slaein.f"
	}

#line 524 "slaein.f"
	i__1 = *n;
#line 524 "slaein.f"
	for (its = 1; its <= i__1; ++its) {
#line 525 "slaein.f"
	    scale = 1.;
#line 526 "slaein.f"
	    vmax = 1.;
#line 527 "slaein.f"
	    vcrit = *bignum;

/*           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector, */
/*             or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector, */
/*           overwriting (xr,xi) on (vr,vi). */

#line 533 "slaein.f"
	    i__2 = i2;
#line 533 "slaein.f"
	    i__3 = i3;
#line 533 "slaein.f"
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {

#line 535 "slaein.f"
		if (work[i__] > vcrit) {
#line 536 "slaein.f"
		    rec = 1. / vmax;
#line 537 "slaein.f"
		    sscal_(n, &rec, &vr[1], &c__1);
#line 538 "slaein.f"
		    sscal_(n, &rec, &vi[1], &c__1);
#line 539 "slaein.f"
		    scale *= rec;
#line 540 "slaein.f"
		    vmax = 1.;
#line 541 "slaein.f"
		    vcrit = *bignum;
#line 542 "slaein.f"
		}

#line 544 "slaein.f"
		xr = vr[i__];
#line 545 "slaein.f"
		xi = vi[i__];
#line 546 "slaein.f"
		if (*rightv) {
#line 547 "slaein.f"
		    i__4 = *n;
#line 547 "slaein.f"
		    for (j = i__ + 1; j <= i__4; ++j) {
#line 548 "slaein.f"
			xr = xr - b[i__ + j * b_dim1] * vr[j] + b[j + 1 + i__ 
				* b_dim1] * vi[j];
#line 549 "slaein.f"
			xi = xi - b[i__ + j * b_dim1] * vi[j] - b[j + 1 + i__ 
				* b_dim1] * vr[j];
#line 550 "slaein.f"
/* L220: */
#line 550 "slaein.f"
		    }
#line 551 "slaein.f"
		} else {
#line 552 "slaein.f"
		    i__4 = i__ - 1;
#line 552 "slaein.f"
		    for (j = 1; j <= i__4; ++j) {
#line 553 "slaein.f"
			xr = xr - b[j + i__ * b_dim1] * vr[j] + b[i__ + 1 + j 
				* b_dim1] * vi[j];
#line 554 "slaein.f"
			xi = xi - b[j + i__ * b_dim1] * vi[j] - b[i__ + 1 + j 
				* b_dim1] * vr[j];
#line 555 "slaein.f"
/* L230: */
#line 555 "slaein.f"
		    }
#line 556 "slaein.f"
		}

#line 558 "slaein.f"
		w = (d__1 = b[i__ + i__ * b_dim1], abs(d__1)) + (d__2 = b[i__ 
			+ 1 + i__ * b_dim1], abs(d__2));
#line 559 "slaein.f"
		if (w > *smlnum) {
#line 560 "slaein.f"
		    if (w < 1.) {
#line 561 "slaein.f"
			w1 = abs(xr) + abs(xi);
#line 562 "slaein.f"
			if (w1 > w * *bignum) {
#line 563 "slaein.f"
			    rec = 1. / w1;
#line 564 "slaein.f"
			    sscal_(n, &rec, &vr[1], &c__1);
#line 565 "slaein.f"
			    sscal_(n, &rec, &vi[1], &c__1);
#line 566 "slaein.f"
			    xr = vr[i__];
#line 567 "slaein.f"
			    xi = vi[i__];
#line 568 "slaein.f"
			    scale *= rec;
#line 569 "slaein.f"
			    vmax *= rec;
#line 570 "slaein.f"
			}
#line 571 "slaein.f"
		    }

/*                 Divide by diagonal element of B. */

#line 575 "slaein.f"
		    sladiv_(&xr, &xi, &b[i__ + i__ * b_dim1], &b[i__ + 1 + 
			    i__ * b_dim1], &vr[i__], &vi[i__]);
/* Computing MAX */
#line 577 "slaein.f"
		    d__3 = (d__1 = vr[i__], abs(d__1)) + (d__2 = vi[i__], abs(
			    d__2));
#line 577 "slaein.f"
		    vmax = max(d__3,vmax);
#line 578 "slaein.f"
		    vcrit = *bignum / vmax;
#line 579 "slaein.f"
		} else {
#line 580 "slaein.f"
		    i__4 = *n;
#line 580 "slaein.f"
		    for (j = 1; j <= i__4; ++j) {
#line 581 "slaein.f"
			vr[j] = 0.;
#line 582 "slaein.f"
			vi[j] = 0.;
#line 583 "slaein.f"
/* L240: */
#line 583 "slaein.f"
		    }
#line 584 "slaein.f"
		    vr[i__] = 1.;
#line 585 "slaein.f"
		    vi[i__] = 1.;
#line 586 "slaein.f"
		    scale = 0.;
#line 587 "slaein.f"
		    vmax = 1.;
#line 588 "slaein.f"
		    vcrit = *bignum;
#line 589 "slaein.f"
		}
#line 590 "slaein.f"
/* L250: */
#line 590 "slaein.f"
	    }

/*           Test for sufficient growth in the norm of (VR,VI). */

#line 594 "slaein.f"
	    vnorm = sasum_(n, &vr[1], &c__1) + sasum_(n, &vi[1], &c__1);
#line 595 "slaein.f"
	    if (vnorm >= growto * scale) {
#line 595 "slaein.f"
		goto L280;
#line 595 "slaein.f"
	    }

/*           Choose a new orthogonal starting vector and try again. */

#line 600 "slaein.f"
	    y = *eps3 / (rootn + 1.);
#line 601 "slaein.f"
	    vr[1] = *eps3;
#line 602 "slaein.f"
	    vi[1] = 0.;

#line 604 "slaein.f"
	    i__3 = *n;
#line 604 "slaein.f"
	    for (i__ = 2; i__ <= i__3; ++i__) {
#line 605 "slaein.f"
		vr[i__] = y;
#line 606 "slaein.f"
		vi[i__] = 0.;
#line 607 "slaein.f"
/* L260: */
#line 607 "slaein.f"
	    }
#line 608 "slaein.f"
	    vr[*n - its + 1] -= *eps3 * rootn;
#line 609 "slaein.f"
/* L270: */
#line 609 "slaein.f"
	}

/*        Failure to find eigenvector in N iterations */

#line 613 "slaein.f"
	*info = 1;

#line 615 "slaein.f"
L280:

/*        Normalize eigenvector. */

#line 619 "slaein.f"
	vnorm = 0.;
#line 620 "slaein.f"
	i__1 = *n;
#line 620 "slaein.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 621 "slaein.f"
	    d__3 = vnorm, d__4 = (d__1 = vr[i__], abs(d__1)) + (d__2 = vi[i__]
		    , abs(d__2));
#line 621 "slaein.f"
	    vnorm = max(d__3,d__4);
#line 622 "slaein.f"
/* L290: */
#line 622 "slaein.f"
	}
#line 623 "slaein.f"
	d__1 = 1. / vnorm;
#line 623 "slaein.f"
	sscal_(n, &d__1, &vr[1], &c__1);
#line 624 "slaein.f"
	d__1 = 1. / vnorm;
#line 624 "slaein.f"
	sscal_(n, &d__1, &vi[1], &c__1);

#line 626 "slaein.f"
    }

#line 628 "slaein.f"
    return 0;

/*     End of SLAEIN */

} /* slaein_ */

