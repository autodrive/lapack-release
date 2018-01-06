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

/* > \date December 2016 */

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
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
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
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer itmax;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), slabad_(
	    doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int slarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal safmax, rtdisc, smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 253 "slahqr.f"
    /* Parameter adjustments */
#line 253 "slahqr.f"
    h_dim1 = *ldh;
#line 253 "slahqr.f"
    h_offset = 1 + h_dim1;
#line 253 "slahqr.f"
    h__ -= h_offset;
#line 253 "slahqr.f"
    --wr;
#line 253 "slahqr.f"
    --wi;
#line 253 "slahqr.f"
    z_dim1 = *ldz;
#line 253 "slahqr.f"
    z_offset = 1 + z_dim1;
#line 253 "slahqr.f"
    z__ -= z_offset;
#line 253 "slahqr.f"

#line 253 "slahqr.f"
    /* Function Body */
#line 253 "slahqr.f"
    *info = 0;

/*     Quick return if possible */

#line 257 "slahqr.f"
    if (*n == 0) {
#line 257 "slahqr.f"
	return 0;
#line 257 "slahqr.f"
    }
#line 259 "slahqr.f"
    if (*ilo == *ihi) {
#line 260 "slahqr.f"
	wr[*ilo] = h__[*ilo + *ilo * h_dim1];
#line 261 "slahqr.f"
	wi[*ilo] = 0.;
#line 262 "slahqr.f"
	return 0;
#line 263 "slahqr.f"
    }

/*     ==== clear out the trash ==== */
#line 266 "slahqr.f"
    i__1 = *ihi - 3;
#line 266 "slahqr.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 267 "slahqr.f"
	h__[j + 2 + j * h_dim1] = 0.;
#line 268 "slahqr.f"
	h__[j + 3 + j * h_dim1] = 0.;
#line 269 "slahqr.f"
/* L10: */
#line 269 "slahqr.f"
    }
#line 270 "slahqr.f"
    if (*ilo <= *ihi - 2) {
#line 270 "slahqr.f"
	h__[*ihi + (*ihi - 2) * h_dim1] = 0.;
#line 270 "slahqr.f"
    }

#line 273 "slahqr.f"
    nh = *ihi - *ilo + 1;
#line 274 "slahqr.f"
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

#line 278 "slahqr.f"
    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
#line 279 "slahqr.f"
    safmax = 1. / safmin;
#line 280 "slahqr.f"
    slabad_(&safmin, &safmax);
#line 281 "slahqr.f"
    ulp = slamch_("PRECISION", (ftnlen)9);
#line 282 "slahqr.f"
    smlnum = safmin * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

#line 288 "slahqr.f"
    if (*wantt) {
#line 289 "slahqr.f"
	i1 = 1;
#line 290 "slahqr.f"
	i2 = *n;
#line 291 "slahqr.f"
    }

/*     ITMAX is the total number of QR iterations allowed. */

#line 295 "slahqr.f"
    itmax = max(10,nh) * 30;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

#line 303 "slahqr.f"
    i__ = *ihi;
#line 304 "slahqr.f"
L20:
#line 305 "slahqr.f"
    l = *ilo;
#line 306 "slahqr.f"
    if (i__ < *ilo) {
#line 306 "slahqr.f"
	goto L160;
#line 306 "slahqr.f"
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 or 2 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

#line 313 "slahqr.f"
    i__1 = itmax;
#line 313 "slahqr.f"
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

#line 317 "slahqr.f"
	i__2 = l + 1;
#line 317 "slahqr.f"
	for (k = i__; k >= i__2; --k) {
#line 318 "slahqr.f"
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= smlnum) {
#line 318 "slahqr.f"
		goto L40;
#line 318 "slahqr.f"
	    }
#line 320 "slahqr.f"
	    tst = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs(d__1)) + (d__2 = 
		    h__[k + k * h_dim1], abs(d__2));
#line 321 "slahqr.f"
	    if (tst == 0.) {
#line 322 "slahqr.f"
		if (k - 2 >= *ilo) {
#line 322 "slahqr.f"
		    tst += (d__1 = h__[k - 1 + (k - 2) * h_dim1], abs(d__1));
#line 322 "slahqr.f"
		}
#line 324 "slahqr.f"
		if (k + 1 <= *ihi) {
#line 324 "slahqr.f"
		    tst += (d__1 = h__[k + 1 + k * h_dim1], abs(d__1));
#line 324 "slahqr.f"
		}
#line 326 "slahqr.f"
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some cases.  ==== */
#line 331 "slahqr.f"
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= ulp * tst) {
/* Computing MAX */
#line 332 "slahqr.f"
		d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)), d__4 = (
			d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
#line 332 "slahqr.f"
		ab = max(d__3,d__4);
/* Computing MIN */
#line 333 "slahqr.f"
		d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)), d__4 = (
			d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
#line 333 "slahqr.f"
		ba = min(d__3,d__4);
/* Computing MAX */
#line 334 "slahqr.f"
		d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)), d__4 = (d__2 =
			 h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], 
			abs(d__2));
#line 334 "slahqr.f"
		aa = max(d__3,d__4);
/* Computing MIN */
#line 336 "slahqr.f"
		d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)), d__4 = (d__2 =
			 h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], 
			abs(d__2));
#line 336 "slahqr.f"
		bb = min(d__3,d__4);
#line 338 "slahqr.f"
		s = aa + ab;
/* Computing MAX */
#line 339 "slahqr.f"
		d__1 = smlnum, d__2 = ulp * (bb * (aa / s));
#line 339 "slahqr.f"
		if (ba * (ab / s) <= max(d__1,d__2)) {
#line 339 "slahqr.f"
		    goto L40;
#line 339 "slahqr.f"
		}
#line 341 "slahqr.f"
	    }
#line 342 "slahqr.f"
/* L30: */
#line 342 "slahqr.f"
	}
#line 343 "slahqr.f"
L40:
#line 344 "slahqr.f"
	l = k;
#line 345 "slahqr.f"
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

#line 349 "slahqr.f"
	    h__[l + (l - 1) * h_dim1] = 0.;
#line 350 "slahqr.f"
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

#line 354 "slahqr.f"
	if (l >= i__ - 1) {
#line 354 "slahqr.f"
	    goto L150;
#line 354 "slahqr.f"
	}

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

#line 361 "slahqr.f"
	if (! (*wantt)) {
#line 362 "slahqr.f"
	    i1 = l;
#line 363 "slahqr.f"
	    i2 = i__;
#line 364 "slahqr.f"
	}

#line 366 "slahqr.f"
	if (its == 10) {

/*           Exceptional shift. */

#line 370 "slahqr.f"
	    s = (d__1 = h__[l + 1 + l * h_dim1], abs(d__1)) + (d__2 = h__[l + 
		    2 + (l + 1) * h_dim1], abs(d__2));
#line 371 "slahqr.f"
	    h11 = s * .75 + h__[l + l * h_dim1];
#line 372 "slahqr.f"
	    h12 = s * -.4375;
#line 373 "slahqr.f"
	    h21 = s;
#line 374 "slahqr.f"
	    h22 = h11;
#line 375 "slahqr.f"
	} else if (its == 20) {

/*           Exceptional shift. */

#line 379 "slahqr.f"
	    s = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1)) + (d__2 = 
		    h__[i__ - 1 + (i__ - 2) * h_dim1], abs(d__2));
#line 380 "slahqr.f"
	    h11 = s * .75 + h__[i__ + i__ * h_dim1];
#line 381 "slahqr.f"
	    h12 = s * -.4375;
#line 382 "slahqr.f"
	    h21 = s;
#line 383 "slahqr.f"
	    h22 = h11;
#line 384 "slahqr.f"
	} else {

/*           Prepare to use Francis' double shift */
/*           (i.e. 2nd degree generalized Rayleigh quotient) */

#line 389 "slahqr.f"
	    h11 = h__[i__ - 1 + (i__ - 1) * h_dim1];
#line 390 "slahqr.f"
	    h21 = h__[i__ + (i__ - 1) * h_dim1];
#line 391 "slahqr.f"
	    h12 = h__[i__ - 1 + i__ * h_dim1];
#line 392 "slahqr.f"
	    h22 = h__[i__ + i__ * h_dim1];
#line 393 "slahqr.f"
	}
#line 394 "slahqr.f"
	s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
#line 395 "slahqr.f"
	if (s == 0.) {
#line 396 "slahqr.f"
	    rt1r = 0.;
#line 397 "slahqr.f"
	    rt1i = 0.;
#line 398 "slahqr.f"
	    rt2r = 0.;
#line 399 "slahqr.f"
	    rt2i = 0.;
#line 400 "slahqr.f"
	} else {
#line 401 "slahqr.f"
	    h11 /= s;
#line 402 "slahqr.f"
	    h21 /= s;
#line 403 "slahqr.f"
	    h12 /= s;
#line 404 "slahqr.f"
	    h22 /= s;
#line 405 "slahqr.f"
	    tr = (h11 + h22) / 2.;
#line 406 "slahqr.f"
	    det = (h11 - tr) * (h22 - tr) - h12 * h21;
#line 407 "slahqr.f"
	    rtdisc = sqrt((abs(det)));
#line 408 "slahqr.f"
	    if (det >= 0.) {

/*              ==== complex conjugate shifts ==== */

#line 412 "slahqr.f"
		rt1r = tr * s;
#line 413 "slahqr.f"
		rt2r = rt1r;
#line 414 "slahqr.f"
		rt1i = rtdisc * s;
#line 415 "slahqr.f"
		rt2i = -rt1i;
#line 416 "slahqr.f"
	    } else {

/*              ==== real shifts (use only one of them)  ==== */

#line 420 "slahqr.f"
		rt1r = tr + rtdisc;
#line 421 "slahqr.f"
		rt2r = tr - rtdisc;
#line 422 "slahqr.f"
		if ((d__1 = rt1r - h22, abs(d__1)) <= (d__2 = rt2r - h22, abs(
			d__2))) {
#line 423 "slahqr.f"
		    rt1r *= s;
#line 424 "slahqr.f"
		    rt2r = rt1r;
#line 425 "slahqr.f"
		} else {
#line 426 "slahqr.f"
		    rt2r *= s;
#line 427 "slahqr.f"
		    rt1r = rt2r;
#line 428 "slahqr.f"
		}
#line 429 "slahqr.f"
		rt1i = 0.;
#line 430 "slahqr.f"
		rt2i = 0.;
#line 431 "slahqr.f"
	    }
#line 432 "slahqr.f"
	}

/*        Look for two consecutive small subdiagonal elements. */

#line 436 "slahqr.f"
	i__2 = l;
#line 436 "slahqr.f"
	for (m = i__ - 2; m >= i__2; --m) {
/*           Determine the effect of starting the double-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible.  (The following uses scaling to avoid */
/*           overflows and most underflows.) */

#line 442 "slahqr.f"
	    h21s = h__[m + 1 + m * h_dim1];
#line 443 "slahqr.f"
	    s = (d__1 = h__[m + m * h_dim1] - rt2r, abs(d__1)) + abs(rt2i) + 
		    abs(h21s);
#line 444 "slahqr.f"
	    h21s = h__[m + 1 + m * h_dim1] / s;
#line 445 "slahqr.f"
	    v[0] = h21s * h__[m + (m + 1) * h_dim1] + (h__[m + m * h_dim1] - 
		    rt1r) * ((h__[m + m * h_dim1] - rt2r) / s) - rt1i * (rt2i 
		    / s);
#line 447 "slahqr.f"
	    v[1] = h21s * (h__[m + m * h_dim1] + h__[m + 1 + (m + 1) * h_dim1]
		     - rt1r - rt2r);
#line 448 "slahqr.f"
	    v[2] = h21s * h__[m + 2 + (m + 1) * h_dim1];
#line 449 "slahqr.f"
	    s = abs(v[0]) + abs(v[1]) + abs(v[2]);
#line 450 "slahqr.f"
	    v[0] /= s;
#line 451 "slahqr.f"
	    v[1] /= s;
#line 452 "slahqr.f"
	    v[2] /= s;
#line 453 "slahqr.f"
	    if (m == l) {
#line 453 "slahqr.f"
		goto L60;
#line 453 "slahqr.f"
	    }
#line 455 "slahqr.f"
	    if ((d__1 = h__[m + (m - 1) * h_dim1], abs(d__1)) * (abs(v[1]) + 
		    abs(v[2])) <= ulp * abs(v[0]) * ((d__2 = h__[m - 1 + (m - 
		    1) * h_dim1], abs(d__2)) + (d__3 = h__[m + m * h_dim1], 
		    abs(d__3)) + (d__4 = h__[m + 1 + (m + 1) * h_dim1], abs(
		    d__4)))) {
#line 455 "slahqr.f"
		goto L60;
#line 455 "slahqr.f"
	    }
#line 458 "slahqr.f"
/* L50: */
#line 458 "slahqr.f"
	}
#line 459 "slahqr.f"
L60:

/*        Double-shift QR step */

#line 463 "slahqr.f"
	i__2 = i__ - 1;
#line 463 "slahqr.f"
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
#line 474 "slahqr.f"
	    i__3 = 3, i__4 = i__ - k + 1;
#line 474 "slahqr.f"
	    nr = min(i__3,i__4);
#line 475 "slahqr.f"
	    if (k > m) {
#line 475 "slahqr.f"
		scopy_(&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
#line 475 "slahqr.f"
	    }
#line 477 "slahqr.f"
	    slarfg_(&nr, v, &v[1], &c__1, &t1);
#line 478 "slahqr.f"
	    if (k > m) {
#line 479 "slahqr.f"
		h__[k + (k - 1) * h_dim1] = v[0];
#line 480 "slahqr.f"
		h__[k + 1 + (k - 1) * h_dim1] = 0.;
#line 481 "slahqr.f"
		if (k < i__ - 1) {
#line 481 "slahqr.f"
		    h__[k + 2 + (k - 1) * h_dim1] = 0.;
#line 481 "slahqr.f"
		}
#line 483 "slahqr.f"
	    } else if (m > l) {
/*               ==== Use the following instead of */
/*               .    H( K, K-1 ) = -H( K, K-1 ) to */
/*               .    avoid a bug when v(2) and v(3) */
/*               .    underflow. ==== */
#line 488 "slahqr.f"
		h__[k + (k - 1) * h_dim1] *= 1. - t1;
#line 489 "slahqr.f"
	    }
#line 490 "slahqr.f"
	    v2 = v[1];
#line 491 "slahqr.f"
	    t2 = t1 * v2;
#line 492 "slahqr.f"
	    if (nr == 3) {
#line 493 "slahqr.f"
		v3 = v[2];
#line 494 "slahqr.f"
		t3 = t1 * v3;

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

#line 499 "slahqr.f"
		i__3 = i2;
#line 499 "slahqr.f"
		for (j = k; j <= i__3; ++j) {
#line 500 "slahqr.f"
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] 
			    + v3 * h__[k + 2 + j * h_dim1];
#line 501 "slahqr.f"
		    h__[k + j * h_dim1] -= sum * t1;
#line 502 "slahqr.f"
		    h__[k + 1 + j * h_dim1] -= sum * t2;
#line 503 "slahqr.f"
		    h__[k + 2 + j * h_dim1] -= sum * t3;
#line 504 "slahqr.f"
/* L70: */
#line 504 "slahqr.f"
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

/* Computing MIN */
#line 509 "slahqr.f"
		i__4 = k + 3;
#line 509 "slahqr.f"
		i__3 = min(i__4,i__);
#line 509 "slahqr.f"
		for (j = i1; j <= i__3; ++j) {
#line 510 "slahqr.f"
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			     + v3 * h__[j + (k + 2) * h_dim1];
#line 511 "slahqr.f"
		    h__[j + k * h_dim1] -= sum * t1;
#line 512 "slahqr.f"
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
#line 513 "slahqr.f"
		    h__[j + (k + 2) * h_dim1] -= sum * t3;
#line 514 "slahqr.f"
/* L80: */
#line 514 "slahqr.f"
		}

#line 516 "slahqr.f"
		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

#line 520 "slahqr.f"
		    i__3 = *ihiz;
#line 520 "slahqr.f"
		    for (j = *iloz; j <= i__3; ++j) {
#line 521 "slahqr.f"
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1] + v3 * z__[j + (k + 2) * z_dim1];
#line 522 "slahqr.f"
			z__[j + k * z_dim1] -= sum * t1;
#line 523 "slahqr.f"
			z__[j + (k + 1) * z_dim1] -= sum * t2;
#line 524 "slahqr.f"
			z__[j + (k + 2) * z_dim1] -= sum * t3;
#line 525 "slahqr.f"
/* L90: */
#line 525 "slahqr.f"
		    }
#line 526 "slahqr.f"
		}
#line 527 "slahqr.f"
	    } else if (nr == 2) {

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

#line 532 "slahqr.f"
		i__3 = i2;
#line 532 "slahqr.f"
		for (j = k; j <= i__3; ++j) {
#line 533 "slahqr.f"
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
#line 534 "slahqr.f"
		    h__[k + j * h_dim1] -= sum * t1;
#line 535 "slahqr.f"
		    h__[k + 1 + j * h_dim1] -= sum * t2;
#line 536 "slahqr.f"
/* L100: */
#line 536 "slahqr.f"
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

#line 541 "slahqr.f"
		i__3 = i__;
#line 541 "slahqr.f"
		for (j = i1; j <= i__3; ++j) {
#line 542 "slahqr.f"
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			    ;
#line 543 "slahqr.f"
		    h__[j + k * h_dim1] -= sum * t1;
#line 544 "slahqr.f"
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
#line 545 "slahqr.f"
/* L110: */
#line 545 "slahqr.f"
		}

#line 547 "slahqr.f"
		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

#line 551 "slahqr.f"
		    i__3 = *ihiz;
#line 551 "slahqr.f"
		    for (j = *iloz; j <= i__3; ++j) {
#line 552 "slahqr.f"
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1];
#line 553 "slahqr.f"
			z__[j + k * z_dim1] -= sum * t1;
#line 554 "slahqr.f"
			z__[j + (k + 1) * z_dim1] -= sum * t2;
#line 555 "slahqr.f"
/* L120: */
#line 555 "slahqr.f"
		    }
#line 556 "slahqr.f"
		}
#line 557 "slahqr.f"
	    }
#line 558 "slahqr.f"
/* L130: */
#line 558 "slahqr.f"
	}

#line 560 "slahqr.f"
/* L140: */
#line 560 "slahqr.f"
    }

/*     Failure to converge in remaining number of iterations */

#line 564 "slahqr.f"
    *info = i__;
#line 565 "slahqr.f"
    return 0;

#line 567 "slahqr.f"
L150:

#line 569 "slahqr.f"
    if (l == i__) {

/*        H(I,I-1) is negligible: one eigenvalue has converged. */

#line 573 "slahqr.f"
	wr[i__] = h__[i__ + i__ * h_dim1];
#line 574 "slahqr.f"
	wi[i__] = 0.;
#line 575 "slahqr.f"
    } else if (l == i__ - 1) {

/*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged. */

/*        Transform the 2-by-2 submatrix to standard Schur form, */
/*        and compute and store the eigenvalues. */

#line 582 "slahqr.f"
	slanv2_(&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ * 
		h_dim1], &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ * 
		h_dim1], &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__], &cs, 
		&sn);

#line 586 "slahqr.f"
	if (*wantt) {

/*           Apply the transformation to the rest of H. */

#line 590 "slahqr.f"
	    if (i2 > i__) {
#line 590 "slahqr.f"
		i__1 = i2 - i__;
#line 590 "slahqr.f"
		srot_(&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh, &h__[
			i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
#line 590 "slahqr.f"
	    }
#line 593 "slahqr.f"
	    i__1 = i__ - i1 - 1;
#line 593 "slahqr.f"
	    srot_(&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ *
		     h_dim1], &c__1, &cs, &sn);
#line 594 "slahqr.f"
	}
#line 595 "slahqr.f"
	if (*wantz) {

/*           Apply the transformation to Z. */

#line 599 "slahqr.f"
	    srot_(&nz, &z__[*iloz + (i__ - 1) * z_dim1], &c__1, &z__[*iloz + 
		    i__ * z_dim1], &c__1, &cs, &sn);
#line 600 "slahqr.f"
	}
#line 601 "slahqr.f"
    }

/*     return to start of the main loop with new value of I. */

#line 605 "slahqr.f"
    i__ = l - 1;
#line 606 "slahqr.f"
    goto L20;

#line 608 "slahqr.f"
L160:
#line 609 "slahqr.f"
    return 0;

/*     End of SLAHQR */

} /* slahqr_ */

