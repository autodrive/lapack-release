#line 1 "zlahqr.f"
/* zlahqr.f -- translated by f2c (version 20100827).
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

#line 1 "zlahqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b ZLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th
e double-shift/single-shift QR algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, */
/*                          IHIZ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLAHQR is an auxiliary routine called by CHSEQR to update the */
/* >    eigenvalues and Schur decomposition already computed by CHSEQR, by */
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
/* >          It is assumed that H is already upper triangular in rows and */
/* >          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1). */
/* >          ZLAHQR works primarily with the Hessenberg submatrix in rows */
/* >          and columns ILO to IHI, but applies transformations to all of */
/* >          H if WANTT is .TRUE.. */
/* >          1 <= ILO <= max(1,IHI); IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX*16 array, dimension (LDH,N) */
/* >          On entry, the upper Hessenberg matrix H. */
/* >          On exit, if INFO is zero and if WANTT is .TRUE., then H */
/* >          is upper triangular in rows and columns ILO:IHI.  If INFO */
/* >          is zero and if WANTT is .FALSE., then the contents of H */
/* >          are unspecified on exit.  The output state of H in case */
/* >          INF is positive is below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          The computed eigenvalues ILO to IHI are stored in the */
/* >          corresponding elements of W. If WANTT is .TRUE., the */
/* >          eigenvalues are stored in the same order as on the diagonal */
/* >          of the Schur form returned in H, with W(i) = H(i,i). */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ,N) */
/* >          If WANTZ is .TRUE., on entry Z must contain the current */
/* >          matrix Z of transformations accumulated by CHSEQR, and on */
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
/* >          .GT. 0: if INFO = i, ZLAHQR failed to compute all the */
/* >                  eigenvalues ILO to IHI in a total of 30 iterations */
/* >                  per eigenvalue; elements i+1:ihi of W contain */
/* >                  those eigenvalues which have been successfully */
/* >                  computed. */
/* > */
/* >                  If INFO .GT. 0 and WANTT is .FALSE., then on exit, */
/* >                  the remaining unconverged eigenvalues are the */
/* >                  eigenvalues of the upper Hessenberg matrix */
/* >                  rows and columns ILO thorugh INFO of the final, */
/* >                  output value of H. */
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

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >     02-96 Based on modifications by */
/* >     David Day, Sandia National Laboratory, USA */
/* > */
/* >     12-04 Further modifications by */
/* >     Ralph Byers, University of Kansas, USA */
/* >     This is a modified version of ZLAHQR from LAPACK version 3.0. */
/* >     It is (1) more robust against overflow and underflow and */
/* >     (2) adopts the more conservative Ahues & Tisseur stopping */
/* >     criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    void z_sqrt(doublecomplex *, doublecomplex *), pow_zi(doublecomplex *, 
	    doublecomplex *, integer *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal s;
    static doublecomplex t, u, v[2], x, y;
    static integer i1, i2;
    static doublecomplex t1;
    static doublereal t2;
    static doublecomplex v2;
    static doublereal aa, ab, ba, bb, h10;
    static doublecomplex h11;
    static doublereal h21;
    static doublecomplex h22, sc;
    static integer nh, nz;
    static doublereal sx;
    static integer jhi;
    static doublecomplex h11s;
    static integer jlo, its;
    static doublereal ulp;
    static doublecomplex sum;
    static doublereal tst;
    static doublecomplex temp;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer itmax;
    static doublereal rtemp;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin, safmax;
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static doublereal smlnum;


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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 252 "zlahqr.f"
    /* Parameter adjustments */
#line 252 "zlahqr.f"
    h_dim1 = *ldh;
#line 252 "zlahqr.f"
    h_offset = 1 + h_dim1;
#line 252 "zlahqr.f"
    h__ -= h_offset;
#line 252 "zlahqr.f"
    --w;
#line 252 "zlahqr.f"
    z_dim1 = *ldz;
#line 252 "zlahqr.f"
    z_offset = 1 + z_dim1;
#line 252 "zlahqr.f"
    z__ -= z_offset;
#line 252 "zlahqr.f"

#line 252 "zlahqr.f"
    /* Function Body */
#line 252 "zlahqr.f"
    *info = 0;

/*     Quick return if possible */

#line 256 "zlahqr.f"
    if (*n == 0) {
#line 256 "zlahqr.f"
	return 0;
#line 256 "zlahqr.f"
    }
#line 258 "zlahqr.f"
    if (*ilo == *ihi) {
#line 259 "zlahqr.f"
	i__1 = *ilo;
#line 259 "zlahqr.f"
	i__2 = *ilo + *ilo * h_dim1;
#line 259 "zlahqr.f"
	w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
#line 260 "zlahqr.f"
	return 0;
#line 261 "zlahqr.f"
    }

/*     ==== clear out the trash ==== */
#line 264 "zlahqr.f"
    i__1 = *ihi - 3;
#line 264 "zlahqr.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 265 "zlahqr.f"
	i__2 = j + 2 + j * h_dim1;
#line 265 "zlahqr.f"
	h__[i__2].r = 0., h__[i__2].i = 0.;
#line 266 "zlahqr.f"
	i__2 = j + 3 + j * h_dim1;
#line 266 "zlahqr.f"
	h__[i__2].r = 0., h__[i__2].i = 0.;
#line 267 "zlahqr.f"
/* L10: */
#line 267 "zlahqr.f"
    }
#line 268 "zlahqr.f"
    if (*ilo <= *ihi - 2) {
#line 268 "zlahqr.f"
	i__1 = *ihi + (*ihi - 2) * h_dim1;
#line 268 "zlahqr.f"
	h__[i__1].r = 0., h__[i__1].i = 0.;
#line 268 "zlahqr.f"
    }
/*     ==== ensure that subdiagonal entries are real ==== */
#line 271 "zlahqr.f"
    if (*wantt) {
#line 272 "zlahqr.f"
	jlo = 1;
#line 273 "zlahqr.f"
	jhi = *n;
#line 274 "zlahqr.f"
    } else {
#line 275 "zlahqr.f"
	jlo = *ilo;
#line 276 "zlahqr.f"
	jhi = *ihi;
#line 277 "zlahqr.f"
    }
#line 278 "zlahqr.f"
    i__1 = *ihi;
#line 278 "zlahqr.f"
    for (i__ = *ilo + 1; i__ <= i__1; ++i__) {
#line 279 "zlahqr.f"
	if (d_imag(&h__[i__ + (i__ - 1) * h_dim1]) != 0.) {
/*           ==== The following redundant normalization */
/*           .    avoids problems with both gradual and */
/*           .    sudden underflow in ABS(H(I,I-1)) ==== */
#line 283 "zlahqr.f"
	    i__2 = i__ + (i__ - 1) * h_dim1;
#line 283 "zlahqr.f"
	    i__3 = i__ + (i__ - 1) * h_dim1;
#line 283 "zlahqr.f"
	    d__3 = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[i__ 
		    + (i__ - 1) * h_dim1]), abs(d__2));
#line 283 "zlahqr.f"
	    z__1.r = h__[i__2].r / d__3, z__1.i = h__[i__2].i / d__3;
#line 283 "zlahqr.f"
	    sc.r = z__1.r, sc.i = z__1.i;
#line 284 "zlahqr.f"
	    d_cnjg(&z__2, &sc);
#line 284 "zlahqr.f"
	    d__1 = z_abs(&sc);
#line 284 "zlahqr.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 284 "zlahqr.f"
	    sc.r = z__1.r, sc.i = z__1.i;
#line 285 "zlahqr.f"
	    i__2 = i__ + (i__ - 1) * h_dim1;
#line 285 "zlahqr.f"
	    d__1 = z_abs(&h__[i__ + (i__ - 1) * h_dim1]);
#line 285 "zlahqr.f"
	    h__[i__2].r = d__1, h__[i__2].i = 0.;
#line 286 "zlahqr.f"
	    i__2 = jhi - i__ + 1;
#line 286 "zlahqr.f"
	    zscal_(&i__2, &sc, &h__[i__ + i__ * h_dim1], ldh);
/* Computing MIN */
#line 287 "zlahqr.f"
	    i__3 = jhi, i__4 = i__ + 1;
#line 287 "zlahqr.f"
	    i__2 = min(i__3,i__4) - jlo + 1;
#line 287 "zlahqr.f"
	    d_cnjg(&z__1, &sc);
#line 287 "zlahqr.f"
	    zscal_(&i__2, &z__1, &h__[jlo + i__ * h_dim1], &c__1);
#line 289 "zlahqr.f"
	    if (*wantz) {
#line 289 "zlahqr.f"
		i__2 = *ihiz - *iloz + 1;
#line 289 "zlahqr.f"
		d_cnjg(&z__1, &sc);
#line 289 "zlahqr.f"
		zscal_(&i__2, &z__1, &z__[*iloz + i__ * z_dim1], &c__1);
#line 289 "zlahqr.f"
	    }
#line 291 "zlahqr.f"
	}
#line 292 "zlahqr.f"
/* L20: */
#line 292 "zlahqr.f"
    }

#line 294 "zlahqr.f"
    nh = *ihi - *ilo + 1;
#line 295 "zlahqr.f"
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

#line 299 "zlahqr.f"
    safmin = dlamch_("SAFE MINIMUM", (ftnlen)12);
#line 300 "zlahqr.f"
    safmax = 1. / safmin;
#line 301 "zlahqr.f"
    dlabad_(&safmin, &safmax);
#line 302 "zlahqr.f"
    ulp = dlamch_("PRECISION", (ftnlen)9);
#line 303 "zlahqr.f"
    smlnum = safmin * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

#line 309 "zlahqr.f"
    if (*wantt) {
#line 310 "zlahqr.f"
	i1 = 1;
#line 311 "zlahqr.f"
	i2 = *n;
#line 312 "zlahqr.f"
    }

/*     ITMAX is the total number of QR iterations allowed. */

#line 316 "zlahqr.f"
    itmax = max(10,nh) * 30;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or */
/*     H(L,L-1) is negligible so that the matrix splits. */

#line 324 "zlahqr.f"
    i__ = *ihi;
#line 325 "zlahqr.f"
L30:
#line 326 "zlahqr.f"
    if (i__ < *ilo) {
#line 326 "zlahqr.f"
	goto L150;
#line 326 "zlahqr.f"
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

#line 333 "zlahqr.f"
    l = *ilo;
#line 334 "zlahqr.f"
    i__1 = itmax;
#line 334 "zlahqr.f"
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

#line 338 "zlahqr.f"
	i__2 = l + 1;
#line 338 "zlahqr.f"
	for (k = i__; k >= i__2; --k) {
#line 339 "zlahqr.f"
	    i__3 = k + (k - 1) * h_dim1;
#line 339 "zlahqr.f"
	    if ((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[k + (k 
		    - 1) * h_dim1]), abs(d__2)) <= smlnum) {
#line 339 "zlahqr.f"
		goto L50;
#line 339 "zlahqr.f"
	    }
#line 341 "zlahqr.f"
	    i__3 = k - 1 + (k - 1) * h_dim1;
#line 341 "zlahqr.f"
	    i__4 = k + k * h_dim1;
#line 341 "zlahqr.f"
	    tst = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[k - 1 
		    + (k - 1) * h_dim1]), abs(d__2)) + ((d__3 = h__[i__4].r, 
		    abs(d__3)) + (d__4 = d_imag(&h__[k + k * h_dim1]), abs(
		    d__4)));
#line 342 "zlahqr.f"
	    if (tst == 0.) {
#line 343 "zlahqr.f"
		if (k - 2 >= *ilo) {
#line 343 "zlahqr.f"
		    i__3 = k - 1 + (k - 2) * h_dim1;
#line 343 "zlahqr.f"
		    tst += (d__1 = h__[i__3].r, abs(d__1));
#line 343 "zlahqr.f"
		}
#line 345 "zlahqr.f"
		if (k + 1 <= *ihi) {
#line 345 "zlahqr.f"
		    i__3 = k + 1 + k * h_dim1;
#line 345 "zlahqr.f"
		    tst += (d__1 = h__[i__3].r, abs(d__1));
#line 345 "zlahqr.f"
		}
#line 347 "zlahqr.f"
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some examples.  ==== */
#line 352 "zlahqr.f"
	    i__3 = k + (k - 1) * h_dim1;
#line 352 "zlahqr.f"
	    if ((d__1 = h__[i__3].r, abs(d__1)) <= ulp * tst) {
/* Computing MAX */
#line 353 "zlahqr.f"
		i__3 = k + (k - 1) * h_dim1;
#line 353 "zlahqr.f"
		i__4 = k - 1 + k * h_dim1;
#line 353 "zlahqr.f"
		d__5 = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + (k - 1) * h_dim1]), abs(d__2)), d__6 = (d__3 = 
			h__[i__4].r, abs(d__3)) + (d__4 = d_imag(&h__[k - 1 + 
			k * h_dim1]), abs(d__4));
#line 353 "zlahqr.f"
		ab = max(d__5,d__6);
/* Computing MIN */
#line 354 "zlahqr.f"
		i__3 = k + (k - 1) * h_dim1;
#line 354 "zlahqr.f"
		i__4 = k - 1 + k * h_dim1;
#line 354 "zlahqr.f"
		d__5 = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + (k - 1) * h_dim1]), abs(d__2)), d__6 = (d__3 = 
			h__[i__4].r, abs(d__3)) + (d__4 = d_imag(&h__[k - 1 + 
			k * h_dim1]), abs(d__4));
#line 354 "zlahqr.f"
		ba = min(d__5,d__6);
#line 355 "zlahqr.f"
		i__3 = k - 1 + (k - 1) * h_dim1;
#line 355 "zlahqr.f"
		i__4 = k + k * h_dim1;
#line 355 "zlahqr.f"
		z__2.r = h__[i__3].r - h__[i__4].r, z__2.i = h__[i__3].i - 
			h__[i__4].i;
#line 355 "zlahqr.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
/* Computing MAX */
#line 355 "zlahqr.f"
		i__5 = k + k * h_dim1;
#line 355 "zlahqr.f"
		d__5 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + k * h_dim1]), abs(d__2)), d__6 = (d__3 = z__1.r, 
			abs(d__3)) + (d__4 = d_imag(&z__1), abs(d__4));
#line 355 "zlahqr.f"
		aa = max(d__5,d__6);
#line 357 "zlahqr.f"
		i__3 = k - 1 + (k - 1) * h_dim1;
#line 357 "zlahqr.f"
		i__4 = k + k * h_dim1;
#line 357 "zlahqr.f"
		z__2.r = h__[i__3].r - h__[i__4].r, z__2.i = h__[i__3].i - 
			h__[i__4].i;
#line 357 "zlahqr.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
/* Computing MIN */
#line 357 "zlahqr.f"
		i__5 = k + k * h_dim1;
#line 357 "zlahqr.f"
		d__5 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + k * h_dim1]), abs(d__2)), d__6 = (d__3 = z__1.r, 
			abs(d__3)) + (d__4 = d_imag(&z__1), abs(d__4));
#line 357 "zlahqr.f"
		bb = min(d__5,d__6);
#line 359 "zlahqr.f"
		s = aa + ab;
/* Computing MAX */
#line 360 "zlahqr.f"
		d__1 = smlnum, d__2 = ulp * (bb * (aa / s));
#line 360 "zlahqr.f"
		if (ba * (ab / s) <= max(d__1,d__2)) {
#line 360 "zlahqr.f"
		    goto L50;
#line 360 "zlahqr.f"
		}
#line 362 "zlahqr.f"
	    }
#line 363 "zlahqr.f"
/* L40: */
#line 363 "zlahqr.f"
	}
#line 364 "zlahqr.f"
L50:
#line 365 "zlahqr.f"
	l = k;
#line 366 "zlahqr.f"
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

#line 370 "zlahqr.f"
	    i__2 = l + (l - 1) * h_dim1;
#line 370 "zlahqr.f"
	    h__[i__2].r = 0., h__[i__2].i = 0.;
#line 371 "zlahqr.f"
	}

/*        Exit from loop if a submatrix of order 1 has split off. */

#line 375 "zlahqr.f"
	if (l >= i__) {
#line 375 "zlahqr.f"
	    goto L140;
#line 375 "zlahqr.f"
	}

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

#line 382 "zlahqr.f"
	if (! (*wantt)) {
#line 383 "zlahqr.f"
	    i1 = l;
#line 384 "zlahqr.f"
	    i2 = i__;
#line 385 "zlahqr.f"
	}

#line 387 "zlahqr.f"
	if (its == 10) {

/*           Exceptional shift. */

#line 391 "zlahqr.f"
	    i__2 = l + 1 + l * h_dim1;
#line 391 "zlahqr.f"
	    s = (d__1 = h__[i__2].r, abs(d__1)) * .75;
#line 392 "zlahqr.f"
	    i__2 = l + l * h_dim1;
#line 392 "zlahqr.f"
	    z__1.r = s + h__[i__2].r, z__1.i = h__[i__2].i;
#line 392 "zlahqr.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 393 "zlahqr.f"
	} else if (its == 20) {

/*           Exceptional shift. */

#line 397 "zlahqr.f"
	    i__2 = i__ + (i__ - 1) * h_dim1;
#line 397 "zlahqr.f"
	    s = (d__1 = h__[i__2].r, abs(d__1)) * .75;
#line 398 "zlahqr.f"
	    i__2 = i__ + i__ * h_dim1;
#line 398 "zlahqr.f"
	    z__1.r = s + h__[i__2].r, z__1.i = h__[i__2].i;
#line 398 "zlahqr.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 399 "zlahqr.f"
	} else {

/*           Wilkinson's shift. */

#line 403 "zlahqr.f"
	    i__2 = i__ + i__ * h_dim1;
#line 403 "zlahqr.f"
	    t.r = h__[i__2].r, t.i = h__[i__2].i;
#line 404 "zlahqr.f"
	    z_sqrt(&z__2, &h__[i__ - 1 + i__ * h_dim1]);
#line 404 "zlahqr.f"
	    z_sqrt(&z__3, &h__[i__ + (i__ - 1) * h_dim1]);
#line 404 "zlahqr.f"
	    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
		    z__3.i + z__2.i * z__3.r;
#line 404 "zlahqr.f"
	    u.r = z__1.r, u.i = z__1.i;
#line 405 "zlahqr.f"
	    s = (d__1 = u.r, abs(d__1)) + (d__2 = d_imag(&u), abs(d__2));
#line 406 "zlahqr.f"
	    if (s != 0.) {
#line 407 "zlahqr.f"
		i__2 = i__ - 1 + (i__ - 1) * h_dim1;
#line 407 "zlahqr.f"
		z__2.r = h__[i__2].r - t.r, z__2.i = h__[i__2].i - t.i;
#line 407 "zlahqr.f"
		z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
#line 407 "zlahqr.f"
		x.r = z__1.r, x.i = z__1.i;
#line 408 "zlahqr.f"
		sx = (d__1 = x.r, abs(d__1)) + (d__2 = d_imag(&x), abs(d__2));
/* Computing MAX */
#line 409 "zlahqr.f"
		d__3 = s, d__4 = (d__1 = x.r, abs(d__1)) + (d__2 = d_imag(&x),
			 abs(d__2));
#line 409 "zlahqr.f"
		s = max(d__3,d__4);
#line 410 "zlahqr.f"
		z__5.r = x.r / s, z__5.i = x.i / s;
#line 410 "zlahqr.f"
		pow_zi(&z__4, &z__5, &c__2);
#line 410 "zlahqr.f"
		z__7.r = u.r / s, z__7.i = u.i / s;
#line 410 "zlahqr.f"
		pow_zi(&z__6, &z__7, &c__2);
#line 410 "zlahqr.f"
		z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
#line 410 "zlahqr.f"
		z_sqrt(&z__2, &z__3);
#line 410 "zlahqr.f"
		z__1.r = s * z__2.r, z__1.i = s * z__2.i;
#line 410 "zlahqr.f"
		y.r = z__1.r, y.i = z__1.i;
#line 411 "zlahqr.f"
		if (sx > 0.) {
#line 412 "zlahqr.f"
		    z__1.r = x.r / sx, z__1.i = x.i / sx;
#line 412 "zlahqr.f"
		    z__2.r = x.r / sx, z__2.i = x.i / sx;
#line 412 "zlahqr.f"
		    if (z__1.r * y.r + d_imag(&z__2) * d_imag(&y) < 0.) {
#line 412 "zlahqr.f"
			z__3.r = -y.r, z__3.i = -y.i;
#line 412 "zlahqr.f"
			y.r = z__3.r, y.i = z__3.i;
#line 412 "zlahqr.f"
		    }
#line 414 "zlahqr.f"
		}
#line 415 "zlahqr.f"
		z__4.r = x.r + y.r, z__4.i = x.i + y.i;
#line 415 "zlahqr.f"
		zladiv_(&z__3, &u, &z__4);
#line 415 "zlahqr.f"
		z__2.r = u.r * z__3.r - u.i * z__3.i, z__2.i = u.r * z__3.i + 
			u.i * z__3.r;
#line 415 "zlahqr.f"
		z__1.r = t.r - z__2.r, z__1.i = t.i - z__2.i;
#line 415 "zlahqr.f"
		t.r = z__1.r, t.i = z__1.i;
#line 416 "zlahqr.f"
	    }
#line 417 "zlahqr.f"
	}

/*        Look for two consecutive small subdiagonal elements. */

#line 421 "zlahqr.f"
	i__2 = l + 1;
#line 421 "zlahqr.f"
	for (m = i__ - 1; m >= i__2; --m) {

/*           Determine the effect of starting the single-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible. */

#line 427 "zlahqr.f"
	    i__3 = m + m * h_dim1;
#line 427 "zlahqr.f"
	    h11.r = h__[i__3].r, h11.i = h__[i__3].i;
#line 428 "zlahqr.f"
	    i__3 = m + 1 + (m + 1) * h_dim1;
#line 428 "zlahqr.f"
	    h22.r = h__[i__3].r, h22.i = h__[i__3].i;
#line 429 "zlahqr.f"
	    z__1.r = h11.r - t.r, z__1.i = h11.i - t.i;
#line 429 "zlahqr.f"
	    h11s.r = z__1.r, h11s.i = z__1.i;
#line 430 "zlahqr.f"
	    i__3 = m + 1 + m * h_dim1;
#line 430 "zlahqr.f"
	    h21 = h__[i__3].r;
#line 431 "zlahqr.f"
	    s = (d__1 = h11s.r, abs(d__1)) + (d__2 = d_imag(&h11s), abs(d__2))
		     + abs(h21);
#line 432 "zlahqr.f"
	    z__1.r = h11s.r / s, z__1.i = h11s.i / s;
#line 432 "zlahqr.f"
	    h11s.r = z__1.r, h11s.i = z__1.i;
#line 433 "zlahqr.f"
	    h21 /= s;
#line 434 "zlahqr.f"
	    v[0].r = h11s.r, v[0].i = h11s.i;
#line 435 "zlahqr.f"
	    v[1].r = h21, v[1].i = 0.;
#line 436 "zlahqr.f"
	    i__3 = m + (m - 1) * h_dim1;
#line 436 "zlahqr.f"
	    h10 = h__[i__3].r;
#line 437 "zlahqr.f"
	    if (abs(h10) * abs(h21) <= ulp * (((d__1 = h11s.r, abs(d__1)) + (
		    d__2 = d_imag(&h11s), abs(d__2))) * ((d__3 = h11.r, abs(
		    d__3)) + (d__4 = d_imag(&h11), abs(d__4)) + ((d__5 = 
		    h22.r, abs(d__5)) + (d__6 = d_imag(&h22), abs(d__6)))))) {
#line 437 "zlahqr.f"
		goto L70;
#line 437 "zlahqr.f"
	    }
#line 440 "zlahqr.f"
/* L60: */
#line 440 "zlahqr.f"
	}
#line 441 "zlahqr.f"
	i__2 = l + l * h_dim1;
#line 441 "zlahqr.f"
	h11.r = h__[i__2].r, h11.i = h__[i__2].i;
#line 442 "zlahqr.f"
	i__2 = l + 1 + (l + 1) * h_dim1;
#line 442 "zlahqr.f"
	h22.r = h__[i__2].r, h22.i = h__[i__2].i;
#line 443 "zlahqr.f"
	z__1.r = h11.r - t.r, z__1.i = h11.i - t.i;
#line 443 "zlahqr.f"
	h11s.r = z__1.r, h11s.i = z__1.i;
#line 444 "zlahqr.f"
	i__2 = l + 1 + l * h_dim1;
#line 444 "zlahqr.f"
	h21 = h__[i__2].r;
#line 445 "zlahqr.f"
	s = (d__1 = h11s.r, abs(d__1)) + (d__2 = d_imag(&h11s), abs(d__2)) + 
		abs(h21);
#line 446 "zlahqr.f"
	z__1.r = h11s.r / s, z__1.i = h11s.i / s;
#line 446 "zlahqr.f"
	h11s.r = z__1.r, h11s.i = z__1.i;
#line 447 "zlahqr.f"
	h21 /= s;
#line 448 "zlahqr.f"
	v[0].r = h11s.r, v[0].i = h11s.i;
#line 449 "zlahqr.f"
	v[1].r = h21, v[1].i = 0.;
#line 450 "zlahqr.f"
L70:

/*        Single-shift QR step */

#line 454 "zlahqr.f"
	i__2 = i__ - 1;
#line 454 "zlahqr.f"
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. */

/*           V(2) is always real before the call to ZLARFG, and hence */
/*           after the call T2 ( = T1*V(2) ) is also real. */

#line 468 "zlahqr.f"
	    if (k > m) {
#line 468 "zlahqr.f"
		zcopy_(&c__2, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
#line 468 "zlahqr.f"
	    }
#line 470 "zlahqr.f"
	    zlarfg_(&c__2, v, &v[1], &c__1, &t1);
#line 471 "zlahqr.f"
	    if (k > m) {
#line 472 "zlahqr.f"
		i__3 = k + (k - 1) * h_dim1;
#line 472 "zlahqr.f"
		h__[i__3].r = v[0].r, h__[i__3].i = v[0].i;
#line 473 "zlahqr.f"
		i__3 = k + 1 + (k - 1) * h_dim1;
#line 473 "zlahqr.f"
		h__[i__3].r = 0., h__[i__3].i = 0.;
#line 474 "zlahqr.f"
	    }
#line 475 "zlahqr.f"
	    v2.r = v[1].r, v2.i = v[1].i;
#line 476 "zlahqr.f"
	    z__1.r = t1.r * v2.r - t1.i * v2.i, z__1.i = t1.r * v2.i + t1.i * 
		    v2.r;
#line 476 "zlahqr.f"
	    t2 = z__1.r;

/*           Apply G from the left to transform the rows of the matrix */
/*           in columns K to I2. */

#line 481 "zlahqr.f"
	    i__3 = i2;
#line 481 "zlahqr.f"
	    for (j = k; j <= i__3; ++j) {
#line 482 "zlahqr.f"
		d_cnjg(&z__3, &t1);
#line 482 "zlahqr.f"
		i__4 = k + j * h_dim1;
#line 482 "zlahqr.f"
		z__2.r = z__3.r * h__[i__4].r - z__3.i * h__[i__4].i, z__2.i =
			 z__3.r * h__[i__4].i + z__3.i * h__[i__4].r;
#line 482 "zlahqr.f"
		i__5 = k + 1 + j * h_dim1;
#line 482 "zlahqr.f"
		z__4.r = t2 * h__[i__5].r, z__4.i = t2 * h__[i__5].i;
#line 482 "zlahqr.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 482 "zlahqr.f"
		sum.r = z__1.r, sum.i = z__1.i;
#line 483 "zlahqr.f"
		i__4 = k + j * h_dim1;
#line 483 "zlahqr.f"
		i__5 = k + j * h_dim1;
#line 483 "zlahqr.f"
		z__1.r = h__[i__5].r - sum.r, z__1.i = h__[i__5].i - sum.i;
#line 483 "zlahqr.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 484 "zlahqr.f"
		i__4 = k + 1 + j * h_dim1;
#line 484 "zlahqr.f"
		i__5 = k + 1 + j * h_dim1;
#line 484 "zlahqr.f"
		z__2.r = sum.r * v2.r - sum.i * v2.i, z__2.i = sum.r * v2.i + 
			sum.i * v2.r;
#line 484 "zlahqr.f"
		z__1.r = h__[i__5].r - z__2.r, z__1.i = h__[i__5].i - z__2.i;
#line 484 "zlahqr.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 485 "zlahqr.f"
/* L80: */
#line 485 "zlahqr.f"
	    }

/*           Apply G from the right to transform the columns of the */
/*           matrix in rows I1 to min(K+2,I). */

/* Computing MIN */
#line 490 "zlahqr.f"
	    i__4 = k + 2;
#line 490 "zlahqr.f"
	    i__3 = min(i__4,i__);
#line 490 "zlahqr.f"
	    for (j = i1; j <= i__3; ++j) {
#line 491 "zlahqr.f"
		i__4 = j + k * h_dim1;
#line 491 "zlahqr.f"
		z__2.r = t1.r * h__[i__4].r - t1.i * h__[i__4].i, z__2.i = 
			t1.r * h__[i__4].i + t1.i * h__[i__4].r;
#line 491 "zlahqr.f"
		i__5 = j + (k + 1) * h_dim1;
#line 491 "zlahqr.f"
		z__3.r = t2 * h__[i__5].r, z__3.i = t2 * h__[i__5].i;
#line 491 "zlahqr.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 491 "zlahqr.f"
		sum.r = z__1.r, sum.i = z__1.i;
#line 492 "zlahqr.f"
		i__4 = j + k * h_dim1;
#line 492 "zlahqr.f"
		i__5 = j + k * h_dim1;
#line 492 "zlahqr.f"
		z__1.r = h__[i__5].r - sum.r, z__1.i = h__[i__5].i - sum.i;
#line 492 "zlahqr.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 493 "zlahqr.f"
		i__4 = j + (k + 1) * h_dim1;
#line 493 "zlahqr.f"
		i__5 = j + (k + 1) * h_dim1;
#line 493 "zlahqr.f"
		d_cnjg(&z__3, &v2);
#line 493 "zlahqr.f"
		z__2.r = sum.r * z__3.r - sum.i * z__3.i, z__2.i = sum.r * 
			z__3.i + sum.i * z__3.r;
#line 493 "zlahqr.f"
		z__1.r = h__[i__5].r - z__2.r, z__1.i = h__[i__5].i - z__2.i;
#line 493 "zlahqr.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 494 "zlahqr.f"
/* L90: */
#line 494 "zlahqr.f"
	    }

#line 496 "zlahqr.f"
	    if (*wantz) {

/*              Accumulate transformations in the matrix Z */

#line 500 "zlahqr.f"
		i__3 = *ihiz;
#line 500 "zlahqr.f"
		for (j = *iloz; j <= i__3; ++j) {
#line 501 "zlahqr.f"
		    i__4 = j + k * z_dim1;
#line 501 "zlahqr.f"
		    z__2.r = t1.r * z__[i__4].r - t1.i * z__[i__4].i, z__2.i =
			     t1.r * z__[i__4].i + t1.i * z__[i__4].r;
#line 501 "zlahqr.f"
		    i__5 = j + (k + 1) * z_dim1;
#line 501 "zlahqr.f"
		    z__3.r = t2 * z__[i__5].r, z__3.i = t2 * z__[i__5].i;
#line 501 "zlahqr.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 501 "zlahqr.f"
		    sum.r = z__1.r, sum.i = z__1.i;
#line 502 "zlahqr.f"
		    i__4 = j + k * z_dim1;
#line 502 "zlahqr.f"
		    i__5 = j + k * z_dim1;
#line 502 "zlahqr.f"
		    z__1.r = z__[i__5].r - sum.r, z__1.i = z__[i__5].i - 
			    sum.i;
#line 502 "zlahqr.f"
		    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 503 "zlahqr.f"
		    i__4 = j + (k + 1) * z_dim1;
#line 503 "zlahqr.f"
		    i__5 = j + (k + 1) * z_dim1;
#line 503 "zlahqr.f"
		    d_cnjg(&z__3, &v2);
#line 503 "zlahqr.f"
		    z__2.r = sum.r * z__3.r - sum.i * z__3.i, z__2.i = sum.r *
			     z__3.i + sum.i * z__3.r;
#line 503 "zlahqr.f"
		    z__1.r = z__[i__5].r - z__2.r, z__1.i = z__[i__5].i - 
			    z__2.i;
#line 503 "zlahqr.f"
		    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 504 "zlahqr.f"
/* L100: */
#line 504 "zlahqr.f"
		}
#line 505 "zlahqr.f"
	    }

#line 507 "zlahqr.f"
	    if (k == m && m > l) {

/*              If the QR step was started at row M > L because two */
/*              consecutive small subdiagonals were found, then extra */
/*              scaling must be performed to ensure that H(M,M-1) remains */
/*              real. */

#line 514 "zlahqr.f"
		z__1.r = 1. - t1.r, z__1.i = 0. - t1.i;
#line 514 "zlahqr.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 515 "zlahqr.f"
		d__1 = z_abs(&temp);
#line 515 "zlahqr.f"
		z__1.r = temp.r / d__1, z__1.i = temp.i / d__1;
#line 515 "zlahqr.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 516 "zlahqr.f"
		i__3 = m + 1 + m * h_dim1;
#line 516 "zlahqr.f"
		i__4 = m + 1 + m * h_dim1;
#line 516 "zlahqr.f"
		d_cnjg(&z__2, &temp);
#line 516 "zlahqr.f"
		z__1.r = h__[i__4].r * z__2.r - h__[i__4].i * z__2.i, z__1.i =
			 h__[i__4].r * z__2.i + h__[i__4].i * z__2.r;
#line 516 "zlahqr.f"
		h__[i__3].r = z__1.r, h__[i__3].i = z__1.i;
#line 517 "zlahqr.f"
		if (m + 2 <= i__) {
#line 517 "zlahqr.f"
		    i__3 = m + 2 + (m + 1) * h_dim1;
#line 517 "zlahqr.f"
		    i__4 = m + 2 + (m + 1) * h_dim1;
#line 517 "zlahqr.f"
		    z__1.r = h__[i__4].r * temp.r - h__[i__4].i * temp.i, 
			    z__1.i = h__[i__4].r * temp.i + h__[i__4].i * 
			    temp.r;
#line 517 "zlahqr.f"
		    h__[i__3].r = z__1.r, h__[i__3].i = z__1.i;
#line 517 "zlahqr.f"
		}
#line 519 "zlahqr.f"
		i__3 = i__;
#line 519 "zlahqr.f"
		for (j = m; j <= i__3; ++j) {
#line 520 "zlahqr.f"
		    if (j != m + 1) {
#line 521 "zlahqr.f"
			if (i2 > j) {
#line 521 "zlahqr.f"
			    i__4 = i2 - j;
#line 521 "zlahqr.f"
			    zscal_(&i__4, &temp, &h__[j + (j + 1) * h_dim1], 
				    ldh);
#line 521 "zlahqr.f"
			}
#line 523 "zlahqr.f"
			i__4 = j - i1;
#line 523 "zlahqr.f"
			d_cnjg(&z__1, &temp);
#line 523 "zlahqr.f"
			zscal_(&i__4, &z__1, &h__[i1 + j * h_dim1], &c__1);
#line 524 "zlahqr.f"
			if (*wantz) {
#line 525 "zlahqr.f"
			    d_cnjg(&z__1, &temp);
#line 525 "zlahqr.f"
			    zscal_(&nz, &z__1, &z__[*iloz + j * z_dim1], &
				    c__1);
#line 527 "zlahqr.f"
			}
#line 528 "zlahqr.f"
		    }
#line 529 "zlahqr.f"
/* L110: */
#line 529 "zlahqr.f"
		}
#line 530 "zlahqr.f"
	    }
#line 531 "zlahqr.f"
/* L120: */
#line 531 "zlahqr.f"
	}

/*        Ensure that H(I,I-1) is real. */

#line 535 "zlahqr.f"
	i__2 = i__ + (i__ - 1) * h_dim1;
#line 535 "zlahqr.f"
	temp.r = h__[i__2].r, temp.i = h__[i__2].i;
#line 536 "zlahqr.f"
	if (d_imag(&temp) != 0.) {
#line 537 "zlahqr.f"
	    rtemp = z_abs(&temp);
#line 538 "zlahqr.f"
	    i__2 = i__ + (i__ - 1) * h_dim1;
#line 538 "zlahqr.f"
	    h__[i__2].r = rtemp, h__[i__2].i = 0.;
#line 539 "zlahqr.f"
	    z__1.r = temp.r / rtemp, z__1.i = temp.i / rtemp;
#line 539 "zlahqr.f"
	    temp.r = z__1.r, temp.i = z__1.i;
#line 540 "zlahqr.f"
	    if (i2 > i__) {
#line 540 "zlahqr.f"
		i__2 = i2 - i__;
#line 540 "zlahqr.f"
		d_cnjg(&z__1, &temp);
#line 540 "zlahqr.f"
		zscal_(&i__2, &z__1, &h__[i__ + (i__ + 1) * h_dim1], ldh);
#line 540 "zlahqr.f"
	    }
#line 542 "zlahqr.f"
	    i__2 = i__ - i1;
#line 542 "zlahqr.f"
	    zscal_(&i__2, &temp, &h__[i1 + i__ * h_dim1], &c__1);
#line 543 "zlahqr.f"
	    if (*wantz) {
#line 544 "zlahqr.f"
		zscal_(&nz, &temp, &z__[*iloz + i__ * z_dim1], &c__1);
#line 545 "zlahqr.f"
	    }
#line 546 "zlahqr.f"
	}

#line 548 "zlahqr.f"
/* L130: */
#line 548 "zlahqr.f"
    }

/*     Failure to converge in remaining number of iterations */

#line 552 "zlahqr.f"
    *info = i__;
#line 553 "zlahqr.f"
    return 0;

#line 555 "zlahqr.f"
L140:

/*     H(I,I-1) is negligible: one eigenvalue has converged. */

#line 559 "zlahqr.f"
    i__1 = i__;
#line 559 "zlahqr.f"
    i__2 = i__ + i__ * h_dim1;
#line 559 "zlahqr.f"
    w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;

/*     return to start of the main loop with new value of I. */

#line 563 "zlahqr.f"
    i__ = l - 1;
#line 564 "zlahqr.f"
    goto L30;

#line 566 "zlahqr.f"
L150:
#line 567 "zlahqr.f"
    return 0;

/*     End of ZLAHQR */

} /* zlahqr_ */

