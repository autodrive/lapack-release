#line 1 "slahqr.f"
/* slahqr.f -- translated by f2c (version 20100827).
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

#line 1 "slahqr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th
e double-shift/single-shift QR algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, */
/*                          ILOZ, IHIZ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLAHQR is an auxiliary routine called by SHSEQR to update the */
/* >    eigenvalues and Schur decomposition already computed by SHSEQR, by */
/* >    dealing with the Hessenberg submatrix in rows and columns ILO to */
/* >    IHI. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          = .TRUE. : the full Schur form T is required; */
/* >          = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          = .TRUE. : the matrix of Schur vectors Z is required; */
/* >          = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          It is assumed that H is already upper quasi-triangular in */
/* >          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless */
/* >          ILO = 1). SLAHQR works primarily with the Hessenberg */
/* >          submatrix in rows and columns ILO to IHI, but applies */
/* >          transformations to all of H if WANTT is .TRUE.. */
/* >          1 <= ILO <= max(1,IHI); IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is REAL array, dimension (LDH,N) */
/* >          On entry, the upper Hessenberg matrix H. */
/* >          On exit, if INFO is zero and if WANTT is .TRUE., H is upper */
/* >          quasi-triangular in rows and columns ILO:IHI, with any */
/* >          2-by-2 diagonal blocks in standard form. If INFO is zero */
/* >          and WANTT is .FALSE., the contents of H are unspecified on */
/* >          exit.  The output state of H if INFO is nonzero is given */
/* >          below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (N) */
/* >          The real and imaginary parts, respectively, of the computed */
/* >          eigenvalues ILO to IHI are stored in the corresponding */
/* >          elements of WR and WI. If two eigenvalues are computed as a */
/* >          complex conjugate pair, they are stored in consecutive */
/* >          elements of WR and WI, say the i-th and (i+1)th, with */
/* >          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the */
/* >          eigenvalues are stored in the same order as on the diagonal */
/* >          of the Schur form returned in H, with WR(i) = H(i,i), and, if */
/* >          H(i:i+1,i:i+1) is a 2-by-2 diagonal block, */
/* >          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* >          ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* >          IHIZ is INTEGER */
/* >          Specify the rows of Z to which transformations must be */
/* >          applied if WANTZ is .TRUE.. */
/* >          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ,N) */
/* >          If WANTZ is .TRUE., on entry Z must contain the current */
/* >          matrix Z of transformations accumulated by SHSEQR, and on */
/* >          exit Z has been updated; transformations are applied only to */
/* >          the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* >          If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           =   0: successful exit */
/* >          .GT. 0: If INFO = i, SLAHQR failed to compute all the */
/* >                  eigenvalues ILO to IHI in a total of 30 iterations */
/* >                  per eigenvalue; elements i+1:ihi of WR and WI */
/* >                  contain those eigenvalues which have been */
/* >                  successfully computed. */
/* > */
/* >                  If INFO .GT. 0 and WANTT is .FALSE., then on exit, */
/* >                  the remaining unconverged eigenvalues are the */
/* >                  eigenvalues of the upper Hessenberg matrix rows */
/* >                  and columns ILO thorugh INFO of the final, output */
/* >                  value of H. */
/* > */
/* >                  If INFO .GT. 0 and WANTT is .TRUE., then on exit */
/* >          (*)       (initial value of H)*U  = U*(final value of H) */
/* >                  where U is an orthognal matrix.    The final */
/* >                  value of H is upper Hessenberg and triangular in */
/* >                  rows and columns INFO+1 through IHI. */
/* > */
/* >                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit */
/* >                      (final value of Z)  = (initial value of Z)*U */
/* >                  where U is the orthogonal matrix in (*) */
/* >                  (regardless of the value of WANTT.) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     02-96 Based on modifications by */
/* >     David Day, Sandia National Laboratory, USA */
/* > */
/* >     12-04 Further modifications by */
/* >     Ralph Byers, University of Kansas, USA */
/* >     This is a modified version of SLAHQR from LAPACK version 3.0. */
/* >     It is (1) more robust against overflow and underflow and */
/* >     (2) adopts the more conservative Ahues & Tisseur stopping */
/* >     criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	*wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, 
	integer *ldz, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal s, v[3];
    static integer i1, i2;
    static doublereal t1, t2, t3, v2, v3, aa, ab, ba, bb, h11, h12, h21, h22, 
	    cs;
    static integer nh;
    static doublereal sn;
    static integer nr;
    static doublereal tr;
    static integer nz;
    static doublereal det, h21s;
    static integer its;
    static doublereal ulp, sum, tst, rt1i, rt2i, rt1r, rt2r;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), scopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    slanv2_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int slarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal safmax, rtdisc, smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ========================================================= */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 255 "slahqr.f"
    /* Parameter adjustments */
#line 255 "slahqr.f"
    h_dim1 = *ldh;
#line 255 "slahqr.f"
    h_offset = 1 + h_dim1;
#line 255 "slahqr.f"
    h__ -= h_offset;
#line 255 "slahqr.f"
    --wr;
#line 255 "slahqr.f"
    --wi;
#line 255 "slahqr.f"
    z_dim1 = *ldz;
#line 255 "slahqr.f"
    z_offset = 1 + z_dim1;
#line 255 "slahqr.f"
    z__ -= z_offset;
#line 255 "slahqr.f"

#line 255 "slahqr.f"
    /* Function Body */
#line 255 "slahqr.f"
    *info = 0;

/*     Quick return if possible */

#line 259 "slahqr.f"
    if (*n == 0) {
#line 259 "slahqr.f"
	return 0;
#line 259 "slahqr.f"
    }
#line 261 "slahqr.f"
    if (*ilo == *ihi) {
#line 262 "slahqr.f"
	wr[*ilo] = h__[*ilo + *ilo * h_dim1];
#line 263 "slahqr.f"
	wi[*ilo] = 0.;
#line 264 "slahqr.f"
	return 0;
#line 265 "slahqr.f"
    }

/*     ==== clear out the trash ==== */
#line 268 "slahqr.f"
    i__1 = *ihi - 3;
#line 268 "slahqr.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 269 "slahqr.f"
	h__[j + 2 + j * h_dim1] = 0.;
#line 270 "slahqr.f"
	h__[j + 3 + j * h_dim1] = 0.;
#line 271 "slahqr.f"
/* L10: */
#line 271 "slahqr.f"
    }
#line 272 "slahqr.f"
    if (*ilo <= *ihi - 2) {
#line 272 "slahqr.f"
	h__[*ihi + (*ihi - 2) * h_dim1] = 0.;
#line 272 "slahqr.f"
    }

#line 275 "slahqr.f"
    nh = *ihi - *ilo + 1;
#line 276 "slahqr.f"
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

#line 280 "slahqr.f"
    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
#line 281 "slahqr.f"
    safmax = 1. / safmin;
#line 282 "slahqr.f"
    slabad_(&safmin, &safmax);
#line 283 "slahqr.f"
    ulp = slamch_("PRECISION", (ftnlen)9);
#line 284 "slahqr.f"
    smlnum = safmin * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

#line 290 "slahqr.f"
    if (*wantt) {
#line 291 "slahqr.f"
	i1 = 1;
#line 292 "slahqr.f"
	i2 = *n;
#line 293 "slahqr.f"
    }

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

#line 301 "slahqr.f"
    i__ = *ihi;
#line 302 "slahqr.f"
L20:
#line 303 "slahqr.f"
    l = *ilo;
#line 304 "slahqr.f"
    if (i__ < *ilo) {
#line 304 "slahqr.f"
	goto L160;
#line 304 "slahqr.f"
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 or 2 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

#line 311 "slahqr.f"
    for (its = 0; its <= 30; ++its) {

/*        Look for a single small subdiagonal element. */

#line 315 "slahqr.f"
	i__1 = l + 1;
#line 315 "slahqr.f"
	for (k = i__; k >= i__1; --k) {
#line 316 "slahqr.f"
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= smlnum) {
#line 316 "slahqr.f"
		goto L40;
#line 316 "slahqr.f"
	    }
#line 318 "slahqr.f"
	    tst = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs(d__1)) + (d__2 = 
		    h__[k + k * h_dim1], abs(d__2));
#line 319 "slahqr.f"
	    if (tst == 0.) {
#line 320 "slahqr.f"
		if (k - 2 >= *ilo) {
#line 320 "slahqr.f"
		    tst += (d__1 = h__[k - 1 + (k - 2) * h_dim1], abs(d__1));
#line 320 "slahqr.f"
		}
#line 322 "slahqr.f"
		if (k + 1 <= *ihi) {
#line 322 "slahqr.f"
		    tst += (d__1 = h__[k + 1 + k * h_dim1], abs(d__1));
#line 322 "slahqr.f"
		}
#line 324 "slahqr.f"
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some cases.  ==== */
#line 329 "slahqr.f"
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= ulp * tst) {
/* Computing MAX */
#line 330 "slahqr.f"
		d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)), d__4 = (
			d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
#line 330 "slahqr.f"
		ab = max(d__3,d__4);
/* Computing MIN */
#line 331 "slahqr.f"
		d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)), d__4 = (
			d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
#line 331 "slahqr.f"
		ba = min(d__3,d__4);
/* Computing MAX */
#line 332 "slahqr.f"
		d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)), d__4 = (d__2 =
			 h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], 
			abs(d__2));
#line 332 "slahqr.f"
		aa = max(d__3,d__4);
/* Computing MIN */
#line 334 "slahqr.f"
		d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)), d__4 = (d__2 =
			 h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], 
			abs(d__2));
#line 334 "slahqr.f"
		bb = min(d__3,d__4);
#line 336 "slahqr.f"
		s = aa + ab;
/* Computing MAX */
#line 337 "slahqr.f"
		d__1 = smlnum, d__2 = ulp * (bb * (aa / s));
#line 337 "slahqr.f"
		if (ba * (ab / s) <= max(d__1,d__2)) {
#line 337 "slahqr.f"
		    goto L40;
#line 337 "slahqr.f"
		}
#line 339 "slahqr.f"
	    }
#line 340 "slahqr.f"
/* L30: */
#line 340 "slahqr.f"
	}
#line 341 "slahqr.f"
L40:
#line 342 "slahqr.f"
	l = k;
#line 343 "slahqr.f"
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

#line 347 "slahqr.f"
	    h__[l + (l - 1) * h_dim1] = 0.;
#line 348 "slahqr.f"
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

#line 352 "slahqr.f"
	if (l >= i__ - 1) {
#line 352 "slahqr.f"
	    goto L150;
#line 352 "slahqr.f"
	}

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

#line 359 "slahqr.f"
	if (! (*wantt)) {
#line 360 "slahqr.f"
	    i1 = l;
#line 361 "slahqr.f"
	    i2 = i__;
#line 362 "slahqr.f"
	}

#line 364 "slahqr.f"
	if (its == 10) {

/*           Exceptional shift. */

#line 368 "slahqr.f"
	    s = (d__1 = h__[l + 1 + l * h_dim1], abs(d__1)) + (d__2 = h__[l + 
		    2 + (l + 1) * h_dim1], abs(d__2));
#line 369 "slahqr.f"
	    h11 = s * .75 + h__[l + l * h_dim1];
#line 370 "slahqr.f"
	    h12 = s * -.4375;
#line 371 "slahqr.f"
	    h21 = s;
#line 372 "slahqr.f"
	    h22 = h11;
#line 373 "slahqr.f"
	} else if (its == 20) {

/*           Exceptional shift. */

#line 377 "slahqr.f"
	    s = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1)) + (d__2 = 
		    h__[i__ - 1 + (i__ - 2) * h_dim1], abs(d__2));
#line 378 "slahqr.f"
	    h11 = s * .75 + h__[i__ + i__ * h_dim1];
#line 379 "slahqr.f"
	    h12 = s * -.4375;
#line 380 "slahqr.f"
	    h21 = s;
#line 381 "slahqr.f"
	    h22 = h11;
#line 382 "slahqr.f"
	} else {

/*           Prepare to use Francis' double shift */
/*           (i.e. 2nd degree generalized Rayleigh quotient) */

#line 387 "slahqr.f"
	    h11 = h__[i__ - 1 + (i__ - 1) * h_dim1];
#line 388 "slahqr.f"
	    h21 = h__[i__ + (i__ - 1) * h_dim1];
#line 389 "slahqr.f"
	    h12 = h__[i__ - 1 + i__ * h_dim1];
#line 390 "slahqr.f"
	    h22 = h__[i__ + i__ * h_dim1];
#line 391 "slahqr.f"
	}
#line 392 "slahqr.f"
	s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
#line 393 "slahqr.f"
	if (s == 0.) {
#line 394 "slahqr.f"
	    rt1r = 0.;
#line 395 "slahqr.f"
	    rt1i = 0.;
#line 396 "slahqr.f"
	    rt2r = 0.;
#line 397 "slahqr.f"
	    rt2i = 0.;
#line 398 "slahqr.f"
	} else {
#line 399 "slahqr.f"
	    h11 /= s;
#line 400 "slahqr.f"
	    h21 /= s;
#line 401 "slahqr.f"
	    h12 /= s;
#line 402 "slahqr.f"
	    h22 /= s;
#line 403 "slahqr.f"
	    tr = (h11 + h22) / 2.;
#line 404 "slahqr.f"
	    det = (h11 - tr) * (h22 - tr) - h12 * h21;
#line 405 "slahqr.f"
	    rtdisc = sqrt((abs(det)));
#line 406 "slahqr.f"
	    if (det >= 0.) {

/*              ==== complex conjugate shifts ==== */

#line 410 "slahqr.f"
		rt1r = tr * s;
#line 411 "slahqr.f"
		rt2r = rt1r;
#line 412 "slahqr.f"
		rt1i = rtdisc * s;
#line 413 "slahqr.f"
		rt2i = -rt1i;
#line 414 "slahqr.f"
	    } else {

/*              ==== real shifts (use only one of them)  ==== */

#line 418 "slahqr.f"
		rt1r = tr + rtdisc;
#line 419 "slahqr.f"
		rt2r = tr - rtdisc;
#line 420 "slahqr.f"
		if ((d__1 = rt1r - h22, abs(d__1)) <= (d__2 = rt2r - h22, abs(
			d__2))) {
#line 421 "slahqr.f"
		    rt1r *= s;
#line 422 "slahqr.f"
		    rt2r = rt1r;
#line 423 "slahqr.f"
		} else {
#line 424 "slahqr.f"
		    rt2r *= s;
#line 425 "slahqr.f"
		    rt1r = rt2r;
#line 426 "slahqr.f"
		}
#line 427 "slahqr.f"
		rt1i = 0.;
#line 428 "slahqr.f"
		rt2i = 0.;
#line 429 "slahqr.f"
	    }
#line 430 "slahqr.f"
	}

/*        Look for two consecutive small subdiagonal elements. */

#line 434 "slahqr.f"
	i__1 = l;
#line 434 "slahqr.f"
	for (m = i__ - 2; m >= i__1; --m) {
/*           Determine the effect of starting the double-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible.  (The following uses scaling to avoid */
/*           overflows and most underflows.) */

#line 440 "slahqr.f"
	    h21s = h__[m + 1 + m * h_dim1];
#line 441 "slahqr.f"
	    s = (d__1 = h__[m + m * h_dim1] - rt2r, abs(d__1)) + abs(rt2i) + 
		    abs(h21s);
#line 442 "slahqr.f"
	    h21s = h__[m + 1 + m * h_dim1] / s;
#line 443 "slahqr.f"
	    v[0] = h21s * h__[m + (m + 1) * h_dim1] + (h__[m + m * h_dim1] - 
		    rt1r) * ((h__[m + m * h_dim1] - rt2r) / s) - rt1i * (rt2i 
		    / s);
#line 445 "slahqr.f"
	    v[1] = h21s * (h__[m + m * h_dim1] + h__[m + 1 + (m + 1) * h_dim1]
		     - rt1r - rt2r);
#line 446 "slahqr.f"
	    v[2] = h21s * h__[m + 2 + (m + 1) * h_dim1];
#line 447 "slahqr.f"
	    s = abs(v[0]) + abs(v[1]) + abs(v[2]);
#line 448 "slahqr.f"
	    v[0] /= s;
#line 449 "slahqr.f"
	    v[1] /= s;
#line 450 "slahqr.f"
	    v[2] /= s;
#line 451 "slahqr.f"
	    if (m == l) {
#line 451 "slahqr.f"
		goto L60;
#line 451 "slahqr.f"
	    }
#line 453 "slahqr.f"
	    if ((d__1 = h__[m + (m - 1) * h_dim1], abs(d__1)) * (abs(v[1]) + 
		    abs(v[2])) <= ulp * abs(v[0]) * ((d__2 = h__[m - 1 + (m - 
		    1) * h_dim1], abs(d__2)) + (d__3 = h__[m + m * h_dim1], 
		    abs(d__3)) + (d__4 = h__[m + 1 + (m + 1) * h_dim1], abs(
		    d__4)))) {
#line 453 "slahqr.f"
		goto L60;
#line 453 "slahqr.f"
	    }
#line 456 "slahqr.f"
/* L50: */
#line 456 "slahqr.f"
	}
#line 457 "slahqr.f"
L60:

/*        Double-shift QR step */

#line 461 "slahqr.f"
	i__1 = i__ - 1;
#line 461 "slahqr.f"
	for (k = m; k <= i__1; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
#line 472 "slahqr.f"
	    i__2 = 3, i__3 = i__ - k + 1;
#line 472 "slahqr.f"
	    nr = min(i__2,i__3);
#line 473 "slahqr.f"
	    if (k > m) {
#line 473 "slahqr.f"
		scopy_(&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
#line 473 "slahqr.f"
	    }
#line 475 "slahqr.f"
	    slarfg_(&nr, v, &v[1], &c__1, &t1);
#line 476 "slahqr.f"
	    if (k > m) {
#line 477 "slahqr.f"
		h__[k + (k - 1) * h_dim1] = v[0];
#line 478 "slahqr.f"
		h__[k + 1 + (k - 1) * h_dim1] = 0.;
#line 479 "slahqr.f"
		if (k < i__ - 1) {
#line 479 "slahqr.f"
		    h__[k + 2 + (k - 1) * h_dim1] = 0.;
#line 479 "slahqr.f"
		}
#line 481 "slahqr.f"
	    } else if (m > l) {
/*               ==== Use the following instead of */
/*               .    H( K, K-1 ) = -H( K, K-1 ) to */
/*               .    avoid a bug when v(2) and v(3) */
/*               .    underflow. ==== */
#line 486 "slahqr.f"
		h__[k + (k - 1) * h_dim1] *= 1. - t1;
#line 487 "slahqr.f"
	    }
#line 488 "slahqr.f"
	    v2 = v[1];
#line 489 "slahqr.f"
	    t2 = t1 * v2;
#line 490 "slahqr.f"
	    if (nr == 3) {
#line 491 "slahqr.f"
		v3 = v[2];
#line 492 "slahqr.f"
		t3 = t1 * v3;

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

#line 497 "slahqr.f"
		i__2 = i2;
#line 497 "slahqr.f"
		for (j = k; j <= i__2; ++j) {
#line 498 "slahqr.f"
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] 
			    + v3 * h__[k + 2 + j * h_dim1];
#line 499 "slahqr.f"
		    h__[k + j * h_dim1] -= sum * t1;
#line 500 "slahqr.f"
		    h__[k + 1 + j * h_dim1] -= sum * t2;
#line 501 "slahqr.f"
		    h__[k + 2 + j * h_dim1] -= sum * t3;
#line 502 "slahqr.f"
/* L70: */
#line 502 "slahqr.f"
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

/* Computing MIN */
#line 507 "slahqr.f"
		i__3 = k + 3;
#line 507 "slahqr.f"
		i__2 = min(i__3,i__);
#line 507 "slahqr.f"
		for (j = i1; j <= i__2; ++j) {
#line 508 "slahqr.f"
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			     + v3 * h__[j + (k + 2) * h_dim1];
#line 509 "slahqr.f"
		    h__[j + k * h_dim1] -= sum * t1;
#line 510 "slahqr.f"
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
#line 511 "slahqr.f"
		    h__[j + (k + 2) * h_dim1] -= sum * t3;
#line 512 "slahqr.f"
/* L80: */
#line 512 "slahqr.f"
		}

#line 514 "slahqr.f"
		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

#line 518 "slahqr.f"
		    i__2 = *ihiz;
#line 518 "slahqr.f"
		    for (j = *iloz; j <= i__2; ++j) {
#line 519 "slahqr.f"
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1] + v3 * z__[j + (k + 2) * z_dim1];
#line 520 "slahqr.f"
			z__[j + k * z_dim1] -= sum * t1;
#line 521 "slahqr.f"
			z__[j + (k + 1) * z_dim1] -= sum * t2;
#line 522 "slahqr.f"
			z__[j + (k + 2) * z_dim1] -= sum * t3;
#line 523 "slahqr.f"
/* L90: */
#line 523 "slahqr.f"
		    }
#line 524 "slahqr.f"
		}
#line 525 "slahqr.f"
	    } else if (nr == 2) {

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

#line 530 "slahqr.f"
		i__2 = i2;
#line 530 "slahqr.f"
		for (j = k; j <= i__2; ++j) {
#line 531 "slahqr.f"
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
#line 532 "slahqr.f"
		    h__[k + j * h_dim1] -= sum * t1;
#line 533 "slahqr.f"
		    h__[k + 1 + j * h_dim1] -= sum * t2;
#line 534 "slahqr.f"
/* L100: */
#line 534 "slahqr.f"
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

#line 539 "slahqr.f"
		i__2 = i__;
#line 539 "slahqr.f"
		for (j = i1; j <= i__2; ++j) {
#line 540 "slahqr.f"
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			    ;
#line 541 "slahqr.f"
		    h__[j + k * h_dim1] -= sum * t1;
#line 542 "slahqr.f"
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
#line 543 "slahqr.f"
/* L110: */
#line 543 "slahqr.f"
		}

#line 545 "slahqr.f"
		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

#line 549 "slahqr.f"
		    i__2 = *ihiz;
#line 549 "slahqr.f"
		    for (j = *iloz; j <= i__2; ++j) {
#line 550 "slahqr.f"
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1];
#line 551 "slahqr.f"
			z__[j + k * z_dim1] -= sum * t1;
#line 552 "slahqr.f"
			z__[j + (k + 1) * z_dim1] -= sum * t2;
#line 553 "slahqr.f"
/* L120: */
#line 553 "slahqr.f"
		    }
#line 554 "slahqr.f"
		}
#line 555 "slahqr.f"
	    }
#line 556 "slahqr.f"
/* L130: */
#line 556 "slahqr.f"
	}

#line 558 "slahqr.f"
/* L140: */
#line 558 "slahqr.f"
    }

/*     Failure to converge in remaining number of iterations */

#line 562 "slahqr.f"
    *info = i__;
#line 563 "slahqr.f"
    return 0;

#line 565 "slahqr.f"
L150:

#line 567 "slahqr.f"
    if (l == i__) {

/*        H(I,I-1) is negligible: one eigenvalue has converged. */

#line 571 "slahqr.f"
	wr[i__] = h__[i__ + i__ * h_dim1];
#line 572 "slahqr.f"
	wi[i__] = 0.;
#line 573 "slahqr.f"
    } else if (l == i__ - 1) {

/*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged. */

/*        Transform the 2-by-2 submatrix to standard Schur form, */
/*        and compute and store the eigenvalues. */

#line 580 "slahqr.f"
	slanv2_(&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ * 
		h_dim1], &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ * 
		h_dim1], &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__], &cs, 
		&sn);

#line 584 "slahqr.f"
	if (*wantt) {

/*           Apply the transformation to the rest of H. */

#line 588 "slahqr.f"
	    if (i2 > i__) {
#line 588 "slahqr.f"
		i__1 = i2 - i__;
#line 588 "slahqr.f"
		srot_(&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh, &h__[
			i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
#line 588 "slahqr.f"
	    }
#line 591 "slahqr.f"
	    i__1 = i__ - i1 - 1;
#line 591 "slahqr.f"
	    srot_(&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ *
		     h_dim1], &c__1, &cs, &sn);
#line 592 "slahqr.f"
	}
#line 593 "slahqr.f"
	if (*wantz) {

/*           Apply the transformation to Z. */

#line 597 "slahqr.f"
	    srot_(&nz, &z__[*iloz + (i__ - 1) * z_dim1], &c__1, &z__[*iloz + 
		    i__ * z_dim1], &c__1, &cs, &sn);
#line 598 "slahqr.f"
	}
#line 599 "slahqr.f"
    }

/*     return to start of the main loop with new value of I. */

#line 603 "slahqr.f"
    i__ = l - 1;
#line 604 "slahqr.f"
    goto L20;

#line 606 "slahqr.f"
L160:
#line 607 "slahqr.f"
    return 0;

/*     End of SLAHQR */

} /* slahqr_ */

