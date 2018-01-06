#line 1 "ctgevc.f"
/* ctgevc.f -- translated by f2c (version 20100827).
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

#line 1 "ctgevc.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CTGEVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTGEVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, */
/*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            P( LDP, * ), S( LDS, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */



/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGEVC computes some or all of the right and/or left eigenvectors of */
/* > a pair of complex matrices (S,P), where S and P are upper triangular. */
/* > Matrix pairs of this type are produced by the generalized Schur */
/* > factorization of a complex matrix pair (A,B): */
/* > */
/* >    A = Q*S*Z**H,  B = Q*P*Z**H */
/* > */
/* > as computed by CGGHRD + CHGEQZ. */
/* > */
/* > The right eigenvector x and the left eigenvector y of (S,P) */
/* > corresponding to an eigenvalue w are defined by: */
/* > */
/* >    S*x = w*P*x,  (y**H)*S = w*(y**H)*P, */
/* > */
/* > where y**H denotes the conjugate tranpose of y. */
/* > The eigenvalues are not input to this routine, but are computed */
/* > directly from the diagonal elements of S and P. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of (S,P), or the products Z*X and/or Q*Y, */
/* > where Z and Q are input matrices. */
/* > If Q and Z are the unitary factors from the generalized Schur */
/* > factorization of a matrix pair (A,B), then Z*X and Q*Y */
/* > are the matrices of right and left eigenvectors of (A,B). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'R': compute right eigenvectors only; */
/* >          = 'L': compute left eigenvectors only; */
/* >          = 'B': compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* >          HOWMNY is CHARACTER*1 */
/* >          = 'A': compute all right and/or left eigenvectors; */
/* >          = 'B': compute all right and/or left eigenvectors, */
/* >                 backtransformed by the matrices in VR and/or VL; */
/* >          = 'S': compute selected right and/or left eigenvectors, */
/* >                 specified by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY='S', SELECT specifies the eigenvectors to be */
/* >          computed.  The eigenvector corresponding to the j-th */
/* >          eigenvalue is computed if SELECT(j) = .TRUE.. */
/* >          Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices S and P.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is COMPLEX array, dimension (LDS,N) */
/* >          The upper triangular matrix S from a generalized Schur */
/* >          factorization, as computed by CHGEQZ. */
/* > \endverbatim */
/* > */
/* > \param[in] LDS */
/* > \verbatim */
/* >          LDS is INTEGER */
/* >          The leading dimension of array S.  LDS >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is COMPLEX array, dimension (LDP,N) */
/* >          The upper triangular matrix P from a generalized Schur */
/* >          factorization, as computed by CHGEQZ.  P must have real */
/* >          diagonal elements. */
/* > \endverbatim */
/* > */
/* > \param[in] LDP */
/* > \verbatim */
/* >          LDP is INTEGER */
/* >          The leading dimension of array P.  LDP >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is COMPLEX array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Q */
/* >          of left Schur vectors returned by CHGEQZ). */
/* >          On exit, if SIDE = 'L' or 'B', VL contains: */
/* >          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P); */
/* >          if HOWMNY = 'B', the matrix Q*Y; */
/* >          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by */
/* >                      SELECT, stored consecutively in the columns of */
/* >                      VL, in the same order as their eigenvalues. */
/* >          Not referenced if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of array VL.  LDVL >= 1, and if */
/* >          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is COMPLEX array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Z */
/* >          of right Schur vectors returned by CHGEQZ). */
/* >          On exit, if SIDE = 'R' or 'B', VR contains: */
/* >          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P); */
/* >          if HOWMNY = 'B', the matrix Z*X; */
/* >          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by */
/* >                      SELECT, stored consecutively in the columns of */
/* >                      VR, in the same order as their eigenvalues. */
/* >          Not referenced if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR.  LDVR >= 1, and if */
/* >          SIDE = 'R' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* >          MM is INTEGER */
/* >          The number of columns in the arrays VL and/or VR. MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of columns in the arrays VL and/or VR actually */
/* >          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M */
/* >          is set to N.  Each selected eigenvector occupies one column. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *s, integer *lds, doublecomplex *p, integer 
	*ldp, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork,
	 integer *info, ftnlen side_len, ftnlen howmny_len)
{
    /* System generated locals */
    integer p_dim1, p_offset, s_dim1, s_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex d__;
    static integer i__, j;
    static doublecomplex ca, cb;
    static integer je, im, jr;
    static doublereal big;
    static logical lsa, lsb;
    static doublereal ulp;
    static doublecomplex sum;
    static integer ibeg, ieig, iend;
    static doublereal dmin__;
    static integer isrc;
    static doublereal temp;
    static doublecomplex suma, sumb;
    static doublereal xmax, scale;
    static logical ilall;
    static integer iside;
    static doublereal sbeta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal small;
    static logical compl;
    static doublereal anorm, bnorm;
    static logical compr, ilbbad;
    static doublereal acoefa, bcoefa, acoeff;
    static doublecomplex bcoeff;
    static logical ilback;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    static doublereal ascale, bscale;
    extern /* Double Complex */ VOID cladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    static doublecomplex salpha;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical ilcomp;
    static integer ihwmny;


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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test the input parameters */

#line 280 "ctgevc.f"
    /* Parameter adjustments */
#line 280 "ctgevc.f"
    --select;
#line 280 "ctgevc.f"
    s_dim1 = *lds;
#line 280 "ctgevc.f"
    s_offset = 1 + s_dim1;
#line 280 "ctgevc.f"
    s -= s_offset;
#line 280 "ctgevc.f"
    p_dim1 = *ldp;
#line 280 "ctgevc.f"
    p_offset = 1 + p_dim1;
#line 280 "ctgevc.f"
    p -= p_offset;
#line 280 "ctgevc.f"
    vl_dim1 = *ldvl;
#line 280 "ctgevc.f"
    vl_offset = 1 + vl_dim1;
#line 280 "ctgevc.f"
    vl -= vl_offset;
#line 280 "ctgevc.f"
    vr_dim1 = *ldvr;
#line 280 "ctgevc.f"
    vr_offset = 1 + vr_dim1;
#line 280 "ctgevc.f"
    vr -= vr_offset;
#line 280 "ctgevc.f"
    --work;
#line 280 "ctgevc.f"
    --rwork;
#line 280 "ctgevc.f"

#line 280 "ctgevc.f"
    /* Function Body */
#line 280 "ctgevc.f"
    if (lsame_(howmny, "A", (ftnlen)1, (ftnlen)1)) {
#line 281 "ctgevc.f"
	ihwmny = 1;
#line 282 "ctgevc.f"
	ilall = TRUE_;
#line 283 "ctgevc.f"
	ilback = FALSE_;
#line 284 "ctgevc.f"
    } else if (lsame_(howmny, "S", (ftnlen)1, (ftnlen)1)) {
#line 285 "ctgevc.f"
	ihwmny = 2;
#line 286 "ctgevc.f"
	ilall = FALSE_;
#line 287 "ctgevc.f"
	ilback = FALSE_;
#line 288 "ctgevc.f"
    } else if (lsame_(howmny, "B", (ftnlen)1, (ftnlen)1)) {
#line 289 "ctgevc.f"
	ihwmny = 3;
#line 290 "ctgevc.f"
	ilall = TRUE_;
#line 291 "ctgevc.f"
	ilback = TRUE_;
#line 292 "ctgevc.f"
    } else {
#line 293 "ctgevc.f"
	ihwmny = -1;
#line 294 "ctgevc.f"
    }

#line 296 "ctgevc.f"
    if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 297 "ctgevc.f"
	iside = 1;
#line 298 "ctgevc.f"
	compl = FALSE_;
#line 299 "ctgevc.f"
	compr = TRUE_;
#line 300 "ctgevc.f"
    } else if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 301 "ctgevc.f"
	iside = 2;
#line 302 "ctgevc.f"
	compl = TRUE_;
#line 303 "ctgevc.f"
	compr = FALSE_;
#line 304 "ctgevc.f"
    } else if (lsame_(side, "B", (ftnlen)1, (ftnlen)1)) {
#line 305 "ctgevc.f"
	iside = 3;
#line 306 "ctgevc.f"
	compl = TRUE_;
#line 307 "ctgevc.f"
	compr = TRUE_;
#line 308 "ctgevc.f"
    } else {
#line 309 "ctgevc.f"
	iside = -1;
#line 310 "ctgevc.f"
    }

#line 312 "ctgevc.f"
    *info = 0;
#line 313 "ctgevc.f"
    if (iside < 0) {
#line 314 "ctgevc.f"
	*info = -1;
#line 315 "ctgevc.f"
    } else if (ihwmny < 0) {
#line 316 "ctgevc.f"
	*info = -2;
#line 317 "ctgevc.f"
    } else if (*n < 0) {
#line 318 "ctgevc.f"
	*info = -4;
#line 319 "ctgevc.f"
    } else if (*lds < max(1,*n)) {
#line 320 "ctgevc.f"
	*info = -6;
#line 321 "ctgevc.f"
    } else if (*ldp < max(1,*n)) {
#line 322 "ctgevc.f"
	*info = -8;
#line 323 "ctgevc.f"
    }
#line 324 "ctgevc.f"
    if (*info != 0) {
#line 325 "ctgevc.f"
	i__1 = -(*info);
#line 325 "ctgevc.f"
	xerbla_("CTGEVC", &i__1, (ftnlen)6);
#line 326 "ctgevc.f"
	return 0;
#line 327 "ctgevc.f"
    }

/*     Count the number of eigenvectors */

#line 331 "ctgevc.f"
    if (! ilall) {
#line 332 "ctgevc.f"
	im = 0;
#line 333 "ctgevc.f"
	i__1 = *n;
#line 333 "ctgevc.f"
	for (j = 1; j <= i__1; ++j) {
#line 334 "ctgevc.f"
	    if (select[j]) {
#line 334 "ctgevc.f"
		++im;
#line 334 "ctgevc.f"
	    }
#line 336 "ctgevc.f"
/* L10: */
#line 336 "ctgevc.f"
	}
#line 337 "ctgevc.f"
    } else {
#line 338 "ctgevc.f"
	im = *n;
#line 339 "ctgevc.f"
    }

/*     Check diagonal of B */

#line 343 "ctgevc.f"
    ilbbad = FALSE_;
#line 344 "ctgevc.f"
    i__1 = *n;
#line 344 "ctgevc.f"
    for (j = 1; j <= i__1; ++j) {
#line 345 "ctgevc.f"
	if (d_imag(&p[j + j * p_dim1]) != 0.) {
#line 345 "ctgevc.f"
	    ilbbad = TRUE_;
#line 345 "ctgevc.f"
	}
#line 347 "ctgevc.f"
/* L20: */
#line 347 "ctgevc.f"
    }

#line 349 "ctgevc.f"
    if (ilbbad) {
#line 350 "ctgevc.f"
	*info = -7;
#line 351 "ctgevc.f"
    } else if (compl && *ldvl < *n || *ldvl < 1) {
#line 352 "ctgevc.f"
	*info = -10;
#line 353 "ctgevc.f"
    } else if (compr && *ldvr < *n || *ldvr < 1) {
#line 354 "ctgevc.f"
	*info = -12;
#line 355 "ctgevc.f"
    } else if (*mm < im) {
#line 356 "ctgevc.f"
	*info = -13;
#line 357 "ctgevc.f"
    }
#line 358 "ctgevc.f"
    if (*info != 0) {
#line 359 "ctgevc.f"
	i__1 = -(*info);
#line 359 "ctgevc.f"
	xerbla_("CTGEVC", &i__1, (ftnlen)6);
#line 360 "ctgevc.f"
	return 0;
#line 361 "ctgevc.f"
    }

/*     Quick return if possible */

#line 365 "ctgevc.f"
    *m = im;
#line 366 "ctgevc.f"
    if (*n == 0) {
#line 366 "ctgevc.f"
	return 0;
#line 366 "ctgevc.f"
    }

/*     Machine Constants */

#line 371 "ctgevc.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 372 "ctgevc.f"
    big = 1. / safmin;
#line 373 "ctgevc.f"
    slabad_(&safmin, &big);
#line 374 "ctgevc.f"
    ulp = slamch_("Epsilon", (ftnlen)7) * slamch_("Base", (ftnlen)4);
#line 375 "ctgevc.f"
    small = safmin * *n / ulp;
#line 376 "ctgevc.f"
    big = 1. / small;
#line 377 "ctgevc.f"
    bignum = 1. / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular */
/*     part of A and B to check for possible overflow in the triangular */
/*     solver. */

#line 383 "ctgevc.f"
    i__1 = s_dim1 + 1;
#line 383 "ctgevc.f"
    anorm = (d__1 = s[i__1].r, abs(d__1)) + (d__2 = d_imag(&s[s_dim1 + 1]), 
	    abs(d__2));
#line 384 "ctgevc.f"
    i__1 = p_dim1 + 1;
#line 384 "ctgevc.f"
    bnorm = (d__1 = p[i__1].r, abs(d__1)) + (d__2 = d_imag(&p[p_dim1 + 1]), 
	    abs(d__2));
#line 385 "ctgevc.f"
    rwork[1] = 0.;
#line 386 "ctgevc.f"
    rwork[*n + 1] = 0.;
#line 387 "ctgevc.f"
    i__1 = *n;
#line 387 "ctgevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 388 "ctgevc.f"
	rwork[j] = 0.;
#line 389 "ctgevc.f"
	rwork[*n + j] = 0.;
#line 390 "ctgevc.f"
	i__2 = j - 1;
#line 390 "ctgevc.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 391 "ctgevc.f"
	    i__3 = i__ + j * s_dim1;
#line 391 "ctgevc.f"
	    rwork[j] += (d__1 = s[i__3].r, abs(d__1)) + (d__2 = d_imag(&s[i__ 
		    + j * s_dim1]), abs(d__2));
#line 392 "ctgevc.f"
	    i__3 = i__ + j * p_dim1;
#line 392 "ctgevc.f"
	    rwork[*n + j] += (d__1 = p[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		    p[i__ + j * p_dim1]), abs(d__2));
#line 393 "ctgevc.f"
/* L30: */
#line 393 "ctgevc.f"
	}
/* Computing MAX */
#line 394 "ctgevc.f"
	i__2 = j + j * s_dim1;
#line 394 "ctgevc.f"
	d__3 = anorm, d__4 = rwork[j] + ((d__1 = s[i__2].r, abs(d__1)) + (
		d__2 = d_imag(&s[j + j * s_dim1]), abs(d__2)));
#line 394 "ctgevc.f"
	anorm = max(d__3,d__4);
/* Computing MAX */
#line 395 "ctgevc.f"
	i__2 = j + j * p_dim1;
#line 395 "ctgevc.f"
	d__3 = bnorm, d__4 = rwork[*n + j] + ((d__1 = p[i__2].r, abs(d__1)) + 
		(d__2 = d_imag(&p[j + j * p_dim1]), abs(d__2)));
#line 395 "ctgevc.f"
	bnorm = max(d__3,d__4);
#line 396 "ctgevc.f"
/* L40: */
#line 396 "ctgevc.f"
    }

#line 398 "ctgevc.f"
    ascale = 1. / max(anorm,safmin);
#line 399 "ctgevc.f"
    bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

#line 403 "ctgevc.f"
    if (compl) {
#line 404 "ctgevc.f"
	ieig = 0;

/*        Main loop over eigenvalues */

#line 408 "ctgevc.f"
	i__1 = *n;
#line 408 "ctgevc.f"
	for (je = 1; je <= i__1; ++je) {
#line 409 "ctgevc.f"
	    if (ilall) {
#line 410 "ctgevc.f"
		ilcomp = TRUE_;
#line 411 "ctgevc.f"
	    } else {
#line 412 "ctgevc.f"
		ilcomp = select[je];
#line 413 "ctgevc.f"
	    }
#line 414 "ctgevc.f"
	    if (ilcomp) {
#line 415 "ctgevc.f"
		++ieig;

#line 417 "ctgevc.f"
		i__2 = je + je * s_dim1;
#line 417 "ctgevc.f"
		i__3 = je + je * p_dim1;
#line 417 "ctgevc.f"
		if ((d__2 = s[i__2].r, abs(d__2)) + (d__3 = d_imag(&s[je + je 
			* s_dim1]), abs(d__3)) <= safmin && (d__1 = p[i__3].r,
			 abs(d__1)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

#line 422 "ctgevc.f"
		    i__2 = *n;
#line 422 "ctgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 423 "ctgevc.f"
			i__3 = jr + ieig * vl_dim1;
#line 423 "ctgevc.f"
			vl[i__3].r = 0., vl[i__3].i = 0.;
#line 424 "ctgevc.f"
/* L50: */
#line 424 "ctgevc.f"
		    }
#line 425 "ctgevc.f"
		    i__2 = ieig + ieig * vl_dim1;
#line 425 "ctgevc.f"
		    vl[i__2].r = 1., vl[i__2].i = 0.;
#line 426 "ctgevc.f"
		    goto L140;
#line 427 "ctgevc.f"
		}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */
/*                   H */
/*                 y  ( a A - b B ) = 0 */

/* Computing MAX */
#line 434 "ctgevc.f"
		i__2 = je + je * s_dim1;
#line 434 "ctgevc.f"
		i__3 = je + je * p_dim1;
#line 434 "ctgevc.f"
		d__4 = ((d__2 = s[i__2].r, abs(d__2)) + (d__3 = d_imag(&s[je 
			+ je * s_dim1]), abs(d__3))) * ascale, d__5 = (d__1 = 
			p[i__3].r, abs(d__1)) * bscale, d__4 = max(d__4,d__5);
#line 434 "ctgevc.f"
		temp = 1. / max(d__4,safmin);
#line 436 "ctgevc.f"
		i__2 = je + je * s_dim1;
#line 436 "ctgevc.f"
		z__2.r = temp * s[i__2].r, z__2.i = temp * s[i__2].i;
#line 436 "ctgevc.f"
		z__1.r = ascale * z__2.r, z__1.i = ascale * z__2.i;
#line 436 "ctgevc.f"
		salpha.r = z__1.r, salpha.i = z__1.i;
#line 437 "ctgevc.f"
		i__2 = je + je * p_dim1;
#line 437 "ctgevc.f"
		sbeta = temp * p[i__2].r * bscale;
#line 438 "ctgevc.f"
		acoeff = sbeta * ascale;
#line 439 "ctgevc.f"
		z__1.r = bscale * salpha.r, z__1.i = bscale * salpha.i;
#line 439 "ctgevc.f"
		bcoeff.r = z__1.r, bcoeff.i = z__1.i;

/*              Scale to avoid underflow */

#line 443 "ctgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
#line 444 "ctgevc.f"
		lsb = (d__1 = salpha.r, abs(d__1)) + (d__2 = d_imag(&salpha), 
			abs(d__2)) >= safmin && (d__3 = bcoeff.r, abs(d__3)) 
			+ (d__4 = d_imag(&bcoeff), abs(d__4)) < small;

#line 447 "ctgevc.f"
		scale = 1.;
#line 448 "ctgevc.f"
		if (lsa) {
#line 448 "ctgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 448 "ctgevc.f"
		}
#line 450 "ctgevc.f"
		if (lsb) {
/* Computing MAX */
#line 450 "ctgevc.f"
		    d__3 = scale, d__4 = small / ((d__1 = salpha.r, abs(d__1))
			     + (d__2 = d_imag(&salpha), abs(d__2))) * min(
			    bnorm,big);
#line 450 "ctgevc.f"
		    scale = max(d__3,d__4);
#line 450 "ctgevc.f"
		}
#line 453 "ctgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 454 "ctgevc.f"
		    d__5 = 1., d__6 = abs(acoeff), d__5 = max(d__5,d__6), 
			    d__6 = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = 
			    d_imag(&bcoeff), abs(d__2));
#line 454 "ctgevc.f"
		    d__3 = scale, d__4 = 1. / (safmin * max(d__5,d__6));
#line 454 "ctgevc.f"
		    scale = min(d__3,d__4);
#line 457 "ctgevc.f"
		    if (lsa) {
#line 458 "ctgevc.f"
			acoeff = ascale * (scale * sbeta);
#line 459 "ctgevc.f"
		    } else {
#line 460 "ctgevc.f"
			acoeff = scale * acoeff;
#line 461 "ctgevc.f"
		    }
#line 462 "ctgevc.f"
		    if (lsb) {
#line 463 "ctgevc.f"
			z__2.r = scale * salpha.r, z__2.i = scale * salpha.i;
#line 463 "ctgevc.f"
			z__1.r = bscale * z__2.r, z__1.i = bscale * z__2.i;
#line 463 "ctgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 464 "ctgevc.f"
		    } else {
#line 465 "ctgevc.f"
			z__1.r = scale * bcoeff.r, z__1.i = scale * bcoeff.i;
#line 465 "ctgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 466 "ctgevc.f"
		    }
#line 467 "ctgevc.f"
		}

#line 469 "ctgevc.f"
		acoefa = abs(acoeff);
#line 470 "ctgevc.f"
		bcoefa = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = d_imag(&
			bcoeff), abs(d__2));
#line 471 "ctgevc.f"
		xmax = 1.;
#line 472 "ctgevc.f"
		i__2 = *n;
#line 472 "ctgevc.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 473 "ctgevc.f"
		    i__3 = jr;
#line 473 "ctgevc.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 474 "ctgevc.f"
/* L60: */
#line 474 "ctgevc.f"
		}
#line 475 "ctgevc.f"
		i__2 = je;
#line 475 "ctgevc.f"
		work[i__2].r = 1., work[i__2].i = 0.;
/* Computing MAX */
#line 476 "ctgevc.f"
		d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, 
			d__1 = max(d__1,d__2);
#line 476 "ctgevc.f"
		dmin__ = max(d__1,safmin);

/*                                              H */
/*              Triangular solve of  (a A - b B)  y = 0 */

/*                                      H */
/*              (rowwise in  (a A - b B) , or columnwise in a A - b B) */

#line 484 "ctgevc.f"
		i__2 = *n;
#line 484 "ctgevc.f"
		for (j = je + 1; j <= i__2; ++j) {

/*                 Compute */
/*                       j-1 */
/*                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k) */
/*                       k=je */
/*                 (Scale if necessary) */

#line 492 "ctgevc.f"
		    temp = 1. / xmax;
#line 493 "ctgevc.f"
		    if (acoefa * rwork[j] + bcoefa * rwork[*n + j] > bignum * 
			    temp) {
#line 495 "ctgevc.f"
			i__3 = j - 1;
#line 495 "ctgevc.f"
			for (jr = je; jr <= i__3; ++jr) {
#line 496 "ctgevc.f"
			    i__4 = jr;
#line 496 "ctgevc.f"
			    i__5 = jr;
#line 496 "ctgevc.f"
			    z__1.r = temp * work[i__5].r, z__1.i = temp * 
				    work[i__5].i;
#line 496 "ctgevc.f"
			    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 497 "ctgevc.f"
/* L70: */
#line 497 "ctgevc.f"
			}
#line 498 "ctgevc.f"
			xmax = 1.;
#line 499 "ctgevc.f"
		    }
#line 500 "ctgevc.f"
		    suma.r = 0., suma.i = 0.;
#line 501 "ctgevc.f"
		    sumb.r = 0., sumb.i = 0.;

#line 503 "ctgevc.f"
		    i__3 = j - 1;
#line 503 "ctgevc.f"
		    for (jr = je; jr <= i__3; ++jr) {
#line 504 "ctgevc.f"
			d_cnjg(&z__3, &s[jr + j * s_dim1]);
#line 504 "ctgevc.f"
			i__4 = jr;
#line 504 "ctgevc.f"
			z__2.r = z__3.r * work[i__4].r - z__3.i * work[i__4]
				.i, z__2.i = z__3.r * work[i__4].i + z__3.i * 
				work[i__4].r;
#line 504 "ctgevc.f"
			z__1.r = suma.r + z__2.r, z__1.i = suma.i + z__2.i;
#line 504 "ctgevc.f"
			suma.r = z__1.r, suma.i = z__1.i;
#line 505 "ctgevc.f"
			d_cnjg(&z__3, &p[jr + j * p_dim1]);
#line 505 "ctgevc.f"
			i__4 = jr;
#line 505 "ctgevc.f"
			z__2.r = z__3.r * work[i__4].r - z__3.i * work[i__4]
				.i, z__2.i = z__3.r * work[i__4].i + z__3.i * 
				work[i__4].r;
#line 505 "ctgevc.f"
			z__1.r = sumb.r + z__2.r, z__1.i = sumb.i + z__2.i;
#line 505 "ctgevc.f"
			sumb.r = z__1.r, sumb.i = z__1.i;
#line 506 "ctgevc.f"
/* L80: */
#line 506 "ctgevc.f"
		    }
#line 507 "ctgevc.f"
		    z__2.r = acoeff * suma.r, z__2.i = acoeff * suma.i;
#line 507 "ctgevc.f"
		    d_cnjg(&z__4, &bcoeff);
#line 507 "ctgevc.f"
		    z__3.r = z__4.r * sumb.r - z__4.i * sumb.i, z__3.i = 
			    z__4.r * sumb.i + z__4.i * sumb.r;
#line 507 "ctgevc.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 507 "ctgevc.f"
		    sum.r = z__1.r, sum.i = z__1.i;

/*                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) ) */

/*                 with scaling and perturbation of the denominator */

#line 513 "ctgevc.f"
		    i__3 = j + j * s_dim1;
#line 513 "ctgevc.f"
		    z__3.r = acoeff * s[i__3].r, z__3.i = acoeff * s[i__3].i;
#line 513 "ctgevc.f"
		    i__4 = j + j * p_dim1;
#line 513 "ctgevc.f"
		    z__4.r = bcoeff.r * p[i__4].r - bcoeff.i * p[i__4].i, 
			    z__4.i = bcoeff.r * p[i__4].i + bcoeff.i * p[i__4]
			    .r;
#line 513 "ctgevc.f"
		    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 513 "ctgevc.f"
		    d_cnjg(&z__1, &z__2);
#line 513 "ctgevc.f"
		    d__.r = z__1.r, d__.i = z__1.i;
#line 514 "ctgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) <= dmin__) {
#line 514 "ctgevc.f"
			z__1.r = dmin__, z__1.i = 0.;
#line 514 "ctgevc.f"
			d__.r = z__1.r, d__.i = z__1.i;
#line 514 "ctgevc.f"
		    }

#line 517 "ctgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) < 1.) {
#line 518 "ctgevc.f"
			if ((d__1 = sum.r, abs(d__1)) + (d__2 = d_imag(&sum), 
				abs(d__2)) >= bignum * ((d__3 = d__.r, abs(
				d__3)) + (d__4 = d_imag(&d__), abs(d__4)))) {
#line 519 "ctgevc.f"
			    temp = 1. / ((d__1 = sum.r, abs(d__1)) + (d__2 = 
				    d_imag(&sum), abs(d__2)));
#line 520 "ctgevc.f"
			    i__3 = j - 1;
#line 520 "ctgevc.f"
			    for (jr = je; jr <= i__3; ++jr) {
#line 521 "ctgevc.f"
				i__4 = jr;
#line 521 "ctgevc.f"
				i__5 = jr;
#line 521 "ctgevc.f"
				z__1.r = temp * work[i__5].r, z__1.i = temp * 
					work[i__5].i;
#line 521 "ctgevc.f"
				work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 522 "ctgevc.f"
/* L90: */
#line 522 "ctgevc.f"
			    }
#line 523 "ctgevc.f"
			    xmax = temp * xmax;
#line 524 "ctgevc.f"
			    z__1.r = temp * sum.r, z__1.i = temp * sum.i;
#line 524 "ctgevc.f"
			    sum.r = z__1.r, sum.i = z__1.i;
#line 525 "ctgevc.f"
			}
#line 526 "ctgevc.f"
		    }
#line 527 "ctgevc.f"
		    i__3 = j;
#line 527 "ctgevc.f"
		    z__2.r = -sum.r, z__2.i = -sum.i;
#line 527 "ctgevc.f"
		    cladiv_(&z__1, &z__2, &d__);
#line 527 "ctgevc.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
/* Computing MAX */
#line 528 "ctgevc.f"
		    i__3 = j;
#line 528 "ctgevc.f"
		    d__3 = xmax, d__4 = (d__1 = work[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&work[j]), abs(d__2));
#line 528 "ctgevc.f"
		    xmax = max(d__3,d__4);
#line 529 "ctgevc.f"
/* L100: */
#line 529 "ctgevc.f"
		}

/*              Back transform eigenvector if HOWMNY='B'. */

#line 533 "ctgevc.f"
		if (ilback) {
#line 534 "ctgevc.f"
		    i__2 = *n + 1 - je;
#line 534 "ctgevc.f"
		    cgemv_("N", n, &i__2, &c_b2, &vl[je * vl_dim1 + 1], ldvl, 
			    &work[je], &c__1, &c_b1, &work[*n + 1], &c__1, (
			    ftnlen)1);
#line 536 "ctgevc.f"
		    isrc = 2;
#line 537 "ctgevc.f"
		    ibeg = 1;
#line 538 "ctgevc.f"
		} else {
#line 539 "ctgevc.f"
		    isrc = 1;
#line 540 "ctgevc.f"
		    ibeg = je;
#line 541 "ctgevc.f"
		}

/*              Copy and scale eigenvector into column of VL */

#line 545 "ctgevc.f"
		xmax = 0.;
#line 546 "ctgevc.f"
		i__2 = *n;
#line 546 "ctgevc.f"
		for (jr = ibeg; jr <= i__2; ++jr) {
/* Computing MAX */
#line 547 "ctgevc.f"
		    i__3 = (isrc - 1) * *n + jr;
#line 547 "ctgevc.f"
		    d__3 = xmax, d__4 = (d__1 = work[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&work[(isrc - 1) * *n + jr]), abs(
			    d__2));
#line 547 "ctgevc.f"
		    xmax = max(d__3,d__4);
#line 548 "ctgevc.f"
/* L110: */
#line 548 "ctgevc.f"
		}

#line 550 "ctgevc.f"
		if (xmax > safmin) {
#line 551 "ctgevc.f"
		    temp = 1. / xmax;
#line 552 "ctgevc.f"
		    i__2 = *n;
#line 552 "ctgevc.f"
		    for (jr = ibeg; jr <= i__2; ++jr) {
#line 553 "ctgevc.f"
			i__3 = jr + ieig * vl_dim1;
#line 553 "ctgevc.f"
			i__4 = (isrc - 1) * *n + jr;
#line 553 "ctgevc.f"
			z__1.r = temp * work[i__4].r, z__1.i = temp * work[
				i__4].i;
#line 553 "ctgevc.f"
			vl[i__3].r = z__1.r, vl[i__3].i = z__1.i;
#line 554 "ctgevc.f"
/* L120: */
#line 554 "ctgevc.f"
		    }
#line 555 "ctgevc.f"
		} else {
#line 556 "ctgevc.f"
		    ibeg = *n + 1;
#line 557 "ctgevc.f"
		}

#line 559 "ctgevc.f"
		i__2 = ibeg - 1;
#line 559 "ctgevc.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 560 "ctgevc.f"
		    i__3 = jr + ieig * vl_dim1;
#line 560 "ctgevc.f"
		    vl[i__3].r = 0., vl[i__3].i = 0.;
#line 561 "ctgevc.f"
/* L130: */
#line 561 "ctgevc.f"
		}

#line 563 "ctgevc.f"
	    }
#line 564 "ctgevc.f"
L140:
#line 564 "ctgevc.f"
	    ;
#line 564 "ctgevc.f"
	}
#line 565 "ctgevc.f"
    }

/*     Right eigenvectors */

#line 569 "ctgevc.f"
    if (compr) {
#line 570 "ctgevc.f"
	ieig = im + 1;

/*        Main loop over eigenvalues */

#line 574 "ctgevc.f"
	for (je = *n; je >= 1; --je) {
#line 575 "ctgevc.f"
	    if (ilall) {
#line 576 "ctgevc.f"
		ilcomp = TRUE_;
#line 577 "ctgevc.f"
	    } else {
#line 578 "ctgevc.f"
		ilcomp = select[je];
#line 579 "ctgevc.f"
	    }
#line 580 "ctgevc.f"
	    if (ilcomp) {
#line 581 "ctgevc.f"
		--ieig;

#line 583 "ctgevc.f"
		i__1 = je + je * s_dim1;
#line 583 "ctgevc.f"
		i__2 = je + je * p_dim1;
#line 583 "ctgevc.f"
		if ((d__2 = s[i__1].r, abs(d__2)) + (d__3 = d_imag(&s[je + je 
			* s_dim1]), abs(d__3)) <= safmin && (d__1 = p[i__2].r,
			 abs(d__1)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

#line 588 "ctgevc.f"
		    i__1 = *n;
#line 588 "ctgevc.f"
		    for (jr = 1; jr <= i__1; ++jr) {
#line 589 "ctgevc.f"
			i__2 = jr + ieig * vr_dim1;
#line 589 "ctgevc.f"
			vr[i__2].r = 0., vr[i__2].i = 0.;
#line 590 "ctgevc.f"
/* L150: */
#line 590 "ctgevc.f"
		    }
#line 591 "ctgevc.f"
		    i__1 = ieig + ieig * vr_dim1;
#line 591 "ctgevc.f"
		    vr[i__1].r = 1., vr[i__1].i = 0.;
#line 592 "ctgevc.f"
		    goto L250;
#line 593 "ctgevc.f"
		}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */

/*              ( a A - b B ) x  = 0 */

/* Computing MAX */
#line 600 "ctgevc.f"
		i__1 = je + je * s_dim1;
#line 600 "ctgevc.f"
		i__2 = je + je * p_dim1;
#line 600 "ctgevc.f"
		d__4 = ((d__2 = s[i__1].r, abs(d__2)) + (d__3 = d_imag(&s[je 
			+ je * s_dim1]), abs(d__3))) * ascale, d__5 = (d__1 = 
			p[i__2].r, abs(d__1)) * bscale, d__4 = max(d__4,d__5);
#line 600 "ctgevc.f"
		temp = 1. / max(d__4,safmin);
#line 602 "ctgevc.f"
		i__1 = je + je * s_dim1;
#line 602 "ctgevc.f"
		z__2.r = temp * s[i__1].r, z__2.i = temp * s[i__1].i;
#line 602 "ctgevc.f"
		z__1.r = ascale * z__2.r, z__1.i = ascale * z__2.i;
#line 602 "ctgevc.f"
		salpha.r = z__1.r, salpha.i = z__1.i;
#line 603 "ctgevc.f"
		i__1 = je + je * p_dim1;
#line 603 "ctgevc.f"
		sbeta = temp * p[i__1].r * bscale;
#line 604 "ctgevc.f"
		acoeff = sbeta * ascale;
#line 605 "ctgevc.f"
		z__1.r = bscale * salpha.r, z__1.i = bscale * salpha.i;
#line 605 "ctgevc.f"
		bcoeff.r = z__1.r, bcoeff.i = z__1.i;

/*              Scale to avoid underflow */

#line 609 "ctgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
#line 610 "ctgevc.f"
		lsb = (d__1 = salpha.r, abs(d__1)) + (d__2 = d_imag(&salpha), 
			abs(d__2)) >= safmin && (d__3 = bcoeff.r, abs(d__3)) 
			+ (d__4 = d_imag(&bcoeff), abs(d__4)) < small;

#line 613 "ctgevc.f"
		scale = 1.;
#line 614 "ctgevc.f"
		if (lsa) {
#line 614 "ctgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 614 "ctgevc.f"
		}
#line 616 "ctgevc.f"
		if (lsb) {
/* Computing MAX */
#line 616 "ctgevc.f"
		    d__3 = scale, d__4 = small / ((d__1 = salpha.r, abs(d__1))
			     + (d__2 = d_imag(&salpha), abs(d__2))) * min(
			    bnorm,big);
#line 616 "ctgevc.f"
		    scale = max(d__3,d__4);
#line 616 "ctgevc.f"
		}
#line 619 "ctgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 620 "ctgevc.f"
		    d__5 = 1., d__6 = abs(acoeff), d__5 = max(d__5,d__6), 
			    d__6 = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = 
			    d_imag(&bcoeff), abs(d__2));
#line 620 "ctgevc.f"
		    d__3 = scale, d__4 = 1. / (safmin * max(d__5,d__6));
#line 620 "ctgevc.f"
		    scale = min(d__3,d__4);
#line 623 "ctgevc.f"
		    if (lsa) {
#line 624 "ctgevc.f"
			acoeff = ascale * (scale * sbeta);
#line 625 "ctgevc.f"
		    } else {
#line 626 "ctgevc.f"
			acoeff = scale * acoeff;
#line 627 "ctgevc.f"
		    }
#line 628 "ctgevc.f"
		    if (lsb) {
#line 629 "ctgevc.f"
			z__2.r = scale * salpha.r, z__2.i = scale * salpha.i;
#line 629 "ctgevc.f"
			z__1.r = bscale * z__2.r, z__1.i = bscale * z__2.i;
#line 629 "ctgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 630 "ctgevc.f"
		    } else {
#line 631 "ctgevc.f"
			z__1.r = scale * bcoeff.r, z__1.i = scale * bcoeff.i;
#line 631 "ctgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 632 "ctgevc.f"
		    }
#line 633 "ctgevc.f"
		}

#line 635 "ctgevc.f"
		acoefa = abs(acoeff);
#line 636 "ctgevc.f"
		bcoefa = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = d_imag(&
			bcoeff), abs(d__2));
#line 637 "ctgevc.f"
		xmax = 1.;
#line 638 "ctgevc.f"
		i__1 = *n;
#line 638 "ctgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 639 "ctgevc.f"
		    i__2 = jr;
#line 639 "ctgevc.f"
		    work[i__2].r = 0., work[i__2].i = 0.;
#line 640 "ctgevc.f"
/* L160: */
#line 640 "ctgevc.f"
		}
#line 641 "ctgevc.f"
		i__1 = je;
#line 641 "ctgevc.f"
		work[i__1].r = 1., work[i__1].i = 0.;
/* Computing MAX */
#line 642 "ctgevc.f"
		d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, 
			d__1 = max(d__1,d__2);
#line 642 "ctgevc.f"
		dmin__ = max(d__1,safmin);

/*              Triangular solve of  (a A - b B) x = 0  (columnwise) */

/*              WORK(1:j-1) contains sums w, */
/*              WORK(j+1:JE) contains x */

#line 649 "ctgevc.f"
		i__1 = je - 1;
#line 649 "ctgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 650 "ctgevc.f"
		    i__2 = jr;
#line 650 "ctgevc.f"
		    i__3 = jr + je * s_dim1;
#line 650 "ctgevc.f"
		    z__2.r = acoeff * s[i__3].r, z__2.i = acoeff * s[i__3].i;
#line 650 "ctgevc.f"
		    i__4 = jr + je * p_dim1;
#line 650 "ctgevc.f"
		    z__3.r = bcoeff.r * p[i__4].r - bcoeff.i * p[i__4].i, 
			    z__3.i = bcoeff.r * p[i__4].i + bcoeff.i * p[i__4]
			    .r;
#line 650 "ctgevc.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 650 "ctgevc.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 651 "ctgevc.f"
/* L170: */
#line 651 "ctgevc.f"
		}
#line 652 "ctgevc.f"
		i__1 = je;
#line 652 "ctgevc.f"
		work[i__1].r = 1., work[i__1].i = 0.;

#line 654 "ctgevc.f"
		for (j = je - 1; j >= 1; --j) {

/*                 Form x(j) := - w(j) / d */
/*                 with scaling and perturbation of the denominator */

#line 659 "ctgevc.f"
		    i__1 = j + j * s_dim1;
#line 659 "ctgevc.f"
		    z__2.r = acoeff * s[i__1].r, z__2.i = acoeff * s[i__1].i;
#line 659 "ctgevc.f"
		    i__2 = j + j * p_dim1;
#line 659 "ctgevc.f"
		    z__3.r = bcoeff.r * p[i__2].r - bcoeff.i * p[i__2].i, 
			    z__3.i = bcoeff.r * p[i__2].i + bcoeff.i * p[i__2]
			    .r;
#line 659 "ctgevc.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 659 "ctgevc.f"
		    d__.r = z__1.r, d__.i = z__1.i;
#line 660 "ctgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) <= dmin__) {
#line 660 "ctgevc.f"
			z__1.r = dmin__, z__1.i = 0.;
#line 660 "ctgevc.f"
			d__.r = z__1.r, d__.i = z__1.i;
#line 660 "ctgevc.f"
		    }

#line 663 "ctgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) < 1.) {
#line 664 "ctgevc.f"
			i__1 = j;
#line 664 "ctgevc.f"
			if ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(
				&work[j]), abs(d__2)) >= bignum * ((d__3 = 
				d__.r, abs(d__3)) + (d__4 = d_imag(&d__), abs(
				d__4)))) {
#line 665 "ctgevc.f"
			    i__1 = j;
#line 665 "ctgevc.f"
			    temp = 1. / ((d__1 = work[i__1].r, abs(d__1)) + (
				    d__2 = d_imag(&work[j]), abs(d__2)));
#line 666 "ctgevc.f"
			    i__1 = je;
#line 666 "ctgevc.f"
			    for (jr = 1; jr <= i__1; ++jr) {
#line 667 "ctgevc.f"
				i__2 = jr;
#line 667 "ctgevc.f"
				i__3 = jr;
#line 667 "ctgevc.f"
				z__1.r = temp * work[i__3].r, z__1.i = temp * 
					work[i__3].i;
#line 667 "ctgevc.f"
				work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 668 "ctgevc.f"
/* L180: */
#line 668 "ctgevc.f"
			    }
#line 669 "ctgevc.f"
			}
#line 670 "ctgevc.f"
		    }

#line 672 "ctgevc.f"
		    i__1 = j;
#line 672 "ctgevc.f"
		    i__2 = j;
#line 672 "ctgevc.f"
		    z__2.r = -work[i__2].r, z__2.i = -work[i__2].i;
#line 672 "ctgevc.f"
		    cladiv_(&z__1, &z__2, &d__);
#line 672 "ctgevc.f"
		    work[i__1].r = z__1.r, work[i__1].i = z__1.i;

#line 674 "ctgevc.f"
		    if (j > 1) {

/*                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling */

#line 678 "ctgevc.f"
			i__1 = j;
#line 678 "ctgevc.f"
			if ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(
				&work[j]), abs(d__2)) > 1.) {
#line 679 "ctgevc.f"
			    i__1 = j;
#line 679 "ctgevc.f"
			    temp = 1. / ((d__1 = work[i__1].r, abs(d__1)) + (
				    d__2 = d_imag(&work[j]), abs(d__2)));
#line 680 "ctgevc.f"
			    if (acoefa * rwork[j] + bcoefa * rwork[*n + j] >= 
				    bignum * temp) {
#line 682 "ctgevc.f"
				i__1 = je;
#line 682 "ctgevc.f"
				for (jr = 1; jr <= i__1; ++jr) {
#line 683 "ctgevc.f"
				    i__2 = jr;
#line 683 "ctgevc.f"
				    i__3 = jr;
#line 683 "ctgevc.f"
				    z__1.r = temp * work[i__3].r, z__1.i = 
					    temp * work[i__3].i;
#line 683 "ctgevc.f"
				    work[i__2].r = z__1.r, work[i__2].i = 
					    z__1.i;
#line 684 "ctgevc.f"
/* L190: */
#line 684 "ctgevc.f"
				}
#line 685 "ctgevc.f"
			    }
#line 686 "ctgevc.f"
			}

#line 688 "ctgevc.f"
			i__1 = j;
#line 688 "ctgevc.f"
			z__1.r = acoeff * work[i__1].r, z__1.i = acoeff * 
				work[i__1].i;
#line 688 "ctgevc.f"
			ca.r = z__1.r, ca.i = z__1.i;
#line 689 "ctgevc.f"
			i__1 = j;
#line 689 "ctgevc.f"
			z__1.r = bcoeff.r * work[i__1].r - bcoeff.i * work[
				i__1].i, z__1.i = bcoeff.r * work[i__1].i + 
				bcoeff.i * work[i__1].r;
#line 689 "ctgevc.f"
			cb.r = z__1.r, cb.i = z__1.i;
#line 690 "ctgevc.f"
			i__1 = j - 1;
#line 690 "ctgevc.f"
			for (jr = 1; jr <= i__1; ++jr) {
#line 691 "ctgevc.f"
			    i__2 = jr;
#line 691 "ctgevc.f"
			    i__3 = jr;
#line 691 "ctgevc.f"
			    i__4 = jr + j * s_dim1;
#line 691 "ctgevc.f"
			    z__3.r = ca.r * s[i__4].r - ca.i * s[i__4].i, 
				    z__3.i = ca.r * s[i__4].i + ca.i * s[i__4]
				    .r;
#line 691 "ctgevc.f"
			    z__2.r = work[i__3].r + z__3.r, z__2.i = work[
				    i__3].i + z__3.i;
#line 691 "ctgevc.f"
			    i__5 = jr + j * p_dim1;
#line 691 "ctgevc.f"
			    z__4.r = cb.r * p[i__5].r - cb.i * p[i__5].i, 
				    z__4.i = cb.r * p[i__5].i + cb.i * p[i__5]
				    .r;
#line 691 "ctgevc.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 691 "ctgevc.f"
			    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 693 "ctgevc.f"
/* L200: */
#line 693 "ctgevc.f"
			}
#line 694 "ctgevc.f"
		    }
#line 695 "ctgevc.f"
/* L210: */
#line 695 "ctgevc.f"
		}

/*              Back transform eigenvector if HOWMNY='B'. */

#line 699 "ctgevc.f"
		if (ilback) {
#line 700 "ctgevc.f"
		    cgemv_("N", n, &je, &c_b2, &vr[vr_offset], ldvr, &work[1],
			     &c__1, &c_b1, &work[*n + 1], &c__1, (ftnlen)1);
#line 702 "ctgevc.f"
		    isrc = 2;
#line 703 "ctgevc.f"
		    iend = *n;
#line 704 "ctgevc.f"
		} else {
#line 705 "ctgevc.f"
		    isrc = 1;
#line 706 "ctgevc.f"
		    iend = je;
#line 707 "ctgevc.f"
		}

/*              Copy and scale eigenvector into column of VR */

#line 711 "ctgevc.f"
		xmax = 0.;
#line 712 "ctgevc.f"
		i__1 = iend;
#line 712 "ctgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
/* Computing MAX */
#line 713 "ctgevc.f"
		    i__2 = (isrc - 1) * *n + jr;
#line 713 "ctgevc.f"
		    d__3 = xmax, d__4 = (d__1 = work[i__2].r, abs(d__1)) + (
			    d__2 = d_imag(&work[(isrc - 1) * *n + jr]), abs(
			    d__2));
#line 713 "ctgevc.f"
		    xmax = max(d__3,d__4);
#line 714 "ctgevc.f"
/* L220: */
#line 714 "ctgevc.f"
		}

#line 716 "ctgevc.f"
		if (xmax > safmin) {
#line 717 "ctgevc.f"
		    temp = 1. / xmax;
#line 718 "ctgevc.f"
		    i__1 = iend;
#line 718 "ctgevc.f"
		    for (jr = 1; jr <= i__1; ++jr) {
#line 719 "ctgevc.f"
			i__2 = jr + ieig * vr_dim1;
#line 719 "ctgevc.f"
			i__3 = (isrc - 1) * *n + jr;
#line 719 "ctgevc.f"
			z__1.r = temp * work[i__3].r, z__1.i = temp * work[
				i__3].i;
#line 719 "ctgevc.f"
			vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
#line 720 "ctgevc.f"
/* L230: */
#line 720 "ctgevc.f"
		    }
#line 721 "ctgevc.f"
		} else {
#line 722 "ctgevc.f"
		    iend = 0;
#line 723 "ctgevc.f"
		}

#line 725 "ctgevc.f"
		i__1 = *n;
#line 725 "ctgevc.f"
		for (jr = iend + 1; jr <= i__1; ++jr) {
#line 726 "ctgevc.f"
		    i__2 = jr + ieig * vr_dim1;
#line 726 "ctgevc.f"
		    vr[i__2].r = 0., vr[i__2].i = 0.;
#line 727 "ctgevc.f"
/* L240: */
#line 727 "ctgevc.f"
		}

#line 729 "ctgevc.f"
	    }
#line 730 "ctgevc.f"
L250:
#line 730 "ctgevc.f"
	    ;
#line 730 "ctgevc.f"
	}
#line 731 "ctgevc.f"
    }

#line 733 "ctgevc.f"
    return 0;

/*     End of CTGEVC */

} /* ctgevc_ */

