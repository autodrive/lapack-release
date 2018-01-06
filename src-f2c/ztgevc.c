#line 1 "ztgevc.f"
/* ztgevc.f -- translated by f2c (version 20100827).
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

#line 1 "ztgevc.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZTGEVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTGEVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, */
/*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         P( LDP, * ), S( LDS, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */



/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGEVC computes some or all of the right and/or left eigenvectors of */
/* > a pair of complex matrices (S,P), where S and P are upper triangular. */
/* > Matrix pairs of this type are produced by the generalized Schur */
/* > factorization of a complex matrix pair (A,B): */
/* > */
/* >    A = Q*S*Z**H,  B = Q*P*Z**H */
/* > */
/* > as computed by ZGGHRD + ZHGEQZ. */
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
/* >          S is COMPLEX*16 array, dimension (LDS,N) */
/* >          The upper triangular matrix S from a generalized Schur */
/* >          factorization, as computed by ZHGEQZ. */
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
/* >          P is COMPLEX*16 array, dimension (LDP,N) */
/* >          The upper triangular matrix P from a generalized Schur */
/* >          factorization, as computed by ZHGEQZ.  P must have real */
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
/* >          VL is COMPLEX*16 array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Q */
/* >          of left Schur vectors returned by ZHGEQZ). */
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
/* >          VR is COMPLEX*16 array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Z */
/* >          of right Schur vectors returned by ZHGEQZ). */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztgevc_(char *side, char *howmny, logical *select, 
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
    static doublereal small;
    static logical compl;
    static doublereal anorm, bnorm;
    static logical compr;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    dlabad_(doublereal *, doublereal *);
    static logical ilbbad;
    static doublereal acoefa, bcoefa, acoeff;
    static doublecomplex bcoeff;
    static logical ilback;
    static doublereal ascale, bscale;
    extern doublereal dlamch_(char *, ftnlen);
    static doublecomplex salpha;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical ilcomp;
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
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

#line 280 "ztgevc.f"
    /* Parameter adjustments */
#line 280 "ztgevc.f"
    --select;
#line 280 "ztgevc.f"
    s_dim1 = *lds;
#line 280 "ztgevc.f"
    s_offset = 1 + s_dim1;
#line 280 "ztgevc.f"
    s -= s_offset;
#line 280 "ztgevc.f"
    p_dim1 = *ldp;
#line 280 "ztgevc.f"
    p_offset = 1 + p_dim1;
#line 280 "ztgevc.f"
    p -= p_offset;
#line 280 "ztgevc.f"
    vl_dim1 = *ldvl;
#line 280 "ztgevc.f"
    vl_offset = 1 + vl_dim1;
#line 280 "ztgevc.f"
    vl -= vl_offset;
#line 280 "ztgevc.f"
    vr_dim1 = *ldvr;
#line 280 "ztgevc.f"
    vr_offset = 1 + vr_dim1;
#line 280 "ztgevc.f"
    vr -= vr_offset;
#line 280 "ztgevc.f"
    --work;
#line 280 "ztgevc.f"
    --rwork;
#line 280 "ztgevc.f"

#line 280 "ztgevc.f"
    /* Function Body */
#line 280 "ztgevc.f"
    if (lsame_(howmny, "A", (ftnlen)1, (ftnlen)1)) {
#line 281 "ztgevc.f"
	ihwmny = 1;
#line 282 "ztgevc.f"
	ilall = TRUE_;
#line 283 "ztgevc.f"
	ilback = FALSE_;
#line 284 "ztgevc.f"
    } else if (lsame_(howmny, "S", (ftnlen)1, (ftnlen)1)) {
#line 285 "ztgevc.f"
	ihwmny = 2;
#line 286 "ztgevc.f"
	ilall = FALSE_;
#line 287 "ztgevc.f"
	ilback = FALSE_;
#line 288 "ztgevc.f"
    } else if (lsame_(howmny, "B", (ftnlen)1, (ftnlen)1)) {
#line 289 "ztgevc.f"
	ihwmny = 3;
#line 290 "ztgevc.f"
	ilall = TRUE_;
#line 291 "ztgevc.f"
	ilback = TRUE_;
#line 292 "ztgevc.f"
    } else {
#line 293 "ztgevc.f"
	ihwmny = -1;
#line 294 "ztgevc.f"
    }

#line 296 "ztgevc.f"
    if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 297 "ztgevc.f"
	iside = 1;
#line 298 "ztgevc.f"
	compl = FALSE_;
#line 299 "ztgevc.f"
	compr = TRUE_;
#line 300 "ztgevc.f"
    } else if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 301 "ztgevc.f"
	iside = 2;
#line 302 "ztgevc.f"
	compl = TRUE_;
#line 303 "ztgevc.f"
	compr = FALSE_;
#line 304 "ztgevc.f"
    } else if (lsame_(side, "B", (ftnlen)1, (ftnlen)1)) {
#line 305 "ztgevc.f"
	iside = 3;
#line 306 "ztgevc.f"
	compl = TRUE_;
#line 307 "ztgevc.f"
	compr = TRUE_;
#line 308 "ztgevc.f"
    } else {
#line 309 "ztgevc.f"
	iside = -1;
#line 310 "ztgevc.f"
    }

#line 312 "ztgevc.f"
    *info = 0;
#line 313 "ztgevc.f"
    if (iside < 0) {
#line 314 "ztgevc.f"
	*info = -1;
#line 315 "ztgevc.f"
    } else if (ihwmny < 0) {
#line 316 "ztgevc.f"
	*info = -2;
#line 317 "ztgevc.f"
    } else if (*n < 0) {
#line 318 "ztgevc.f"
	*info = -4;
#line 319 "ztgevc.f"
    } else if (*lds < max(1,*n)) {
#line 320 "ztgevc.f"
	*info = -6;
#line 321 "ztgevc.f"
    } else if (*ldp < max(1,*n)) {
#line 322 "ztgevc.f"
	*info = -8;
#line 323 "ztgevc.f"
    }
#line 324 "ztgevc.f"
    if (*info != 0) {
#line 325 "ztgevc.f"
	i__1 = -(*info);
#line 325 "ztgevc.f"
	xerbla_("ZTGEVC", &i__1, (ftnlen)6);
#line 326 "ztgevc.f"
	return 0;
#line 327 "ztgevc.f"
    }

/*     Count the number of eigenvectors */

#line 331 "ztgevc.f"
    if (! ilall) {
#line 332 "ztgevc.f"
	im = 0;
#line 333 "ztgevc.f"
	i__1 = *n;
#line 333 "ztgevc.f"
	for (j = 1; j <= i__1; ++j) {
#line 334 "ztgevc.f"
	    if (select[j]) {
#line 334 "ztgevc.f"
		++im;
#line 334 "ztgevc.f"
	    }
#line 336 "ztgevc.f"
/* L10: */
#line 336 "ztgevc.f"
	}
#line 337 "ztgevc.f"
    } else {
#line 338 "ztgevc.f"
	im = *n;
#line 339 "ztgevc.f"
    }

/*     Check diagonal of B */

#line 343 "ztgevc.f"
    ilbbad = FALSE_;
#line 344 "ztgevc.f"
    i__1 = *n;
#line 344 "ztgevc.f"
    for (j = 1; j <= i__1; ++j) {
#line 345 "ztgevc.f"
	if (d_imag(&p[j + j * p_dim1]) != 0.) {
#line 345 "ztgevc.f"
	    ilbbad = TRUE_;
#line 345 "ztgevc.f"
	}
#line 347 "ztgevc.f"
/* L20: */
#line 347 "ztgevc.f"
    }

#line 349 "ztgevc.f"
    if (ilbbad) {
#line 350 "ztgevc.f"
	*info = -7;
#line 351 "ztgevc.f"
    } else if (compl && *ldvl < *n || *ldvl < 1) {
#line 352 "ztgevc.f"
	*info = -10;
#line 353 "ztgevc.f"
    } else if (compr && *ldvr < *n || *ldvr < 1) {
#line 354 "ztgevc.f"
	*info = -12;
#line 355 "ztgevc.f"
    } else if (*mm < im) {
#line 356 "ztgevc.f"
	*info = -13;
#line 357 "ztgevc.f"
    }
#line 358 "ztgevc.f"
    if (*info != 0) {
#line 359 "ztgevc.f"
	i__1 = -(*info);
#line 359 "ztgevc.f"
	xerbla_("ZTGEVC", &i__1, (ftnlen)6);
#line 360 "ztgevc.f"
	return 0;
#line 361 "ztgevc.f"
    }

/*     Quick return if possible */

#line 365 "ztgevc.f"
    *m = im;
#line 366 "ztgevc.f"
    if (*n == 0) {
#line 366 "ztgevc.f"
	return 0;
#line 366 "ztgevc.f"
    }

/*     Machine Constants */

#line 371 "ztgevc.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 372 "ztgevc.f"
    big = 1. / safmin;
#line 373 "ztgevc.f"
    dlabad_(&safmin, &big);
#line 374 "ztgevc.f"
    ulp = dlamch_("Epsilon", (ftnlen)7) * dlamch_("Base", (ftnlen)4);
#line 375 "ztgevc.f"
    small = safmin * *n / ulp;
#line 376 "ztgevc.f"
    big = 1. / small;
#line 377 "ztgevc.f"
    bignum = 1. / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular */
/*     part of A and B to check for possible overflow in the triangular */
/*     solver. */

#line 383 "ztgevc.f"
    i__1 = s_dim1 + 1;
#line 383 "ztgevc.f"
    anorm = (d__1 = s[i__1].r, abs(d__1)) + (d__2 = d_imag(&s[s_dim1 + 1]), 
	    abs(d__2));
#line 384 "ztgevc.f"
    i__1 = p_dim1 + 1;
#line 384 "ztgevc.f"
    bnorm = (d__1 = p[i__1].r, abs(d__1)) + (d__2 = d_imag(&p[p_dim1 + 1]), 
	    abs(d__2));
#line 385 "ztgevc.f"
    rwork[1] = 0.;
#line 386 "ztgevc.f"
    rwork[*n + 1] = 0.;
#line 387 "ztgevc.f"
    i__1 = *n;
#line 387 "ztgevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 388 "ztgevc.f"
	rwork[j] = 0.;
#line 389 "ztgevc.f"
	rwork[*n + j] = 0.;
#line 390 "ztgevc.f"
	i__2 = j - 1;
#line 390 "ztgevc.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 391 "ztgevc.f"
	    i__3 = i__ + j * s_dim1;
#line 391 "ztgevc.f"
	    rwork[j] += (d__1 = s[i__3].r, abs(d__1)) + (d__2 = d_imag(&s[i__ 
		    + j * s_dim1]), abs(d__2));
#line 392 "ztgevc.f"
	    i__3 = i__ + j * p_dim1;
#line 392 "ztgevc.f"
	    rwork[*n + j] += (d__1 = p[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		    p[i__ + j * p_dim1]), abs(d__2));
#line 393 "ztgevc.f"
/* L30: */
#line 393 "ztgevc.f"
	}
/* Computing MAX */
#line 394 "ztgevc.f"
	i__2 = j + j * s_dim1;
#line 394 "ztgevc.f"
	d__3 = anorm, d__4 = rwork[j] + ((d__1 = s[i__2].r, abs(d__1)) + (
		d__2 = d_imag(&s[j + j * s_dim1]), abs(d__2)));
#line 394 "ztgevc.f"
	anorm = max(d__3,d__4);
/* Computing MAX */
#line 395 "ztgevc.f"
	i__2 = j + j * p_dim1;
#line 395 "ztgevc.f"
	d__3 = bnorm, d__4 = rwork[*n + j] + ((d__1 = p[i__2].r, abs(d__1)) + 
		(d__2 = d_imag(&p[j + j * p_dim1]), abs(d__2)));
#line 395 "ztgevc.f"
	bnorm = max(d__3,d__4);
#line 396 "ztgevc.f"
/* L40: */
#line 396 "ztgevc.f"
    }

#line 398 "ztgevc.f"
    ascale = 1. / max(anorm,safmin);
#line 399 "ztgevc.f"
    bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

#line 403 "ztgevc.f"
    if (compl) {
#line 404 "ztgevc.f"
	ieig = 0;

/*        Main loop over eigenvalues */

#line 408 "ztgevc.f"
	i__1 = *n;
#line 408 "ztgevc.f"
	for (je = 1; je <= i__1; ++je) {
#line 409 "ztgevc.f"
	    if (ilall) {
#line 410 "ztgevc.f"
		ilcomp = TRUE_;
#line 411 "ztgevc.f"
	    } else {
#line 412 "ztgevc.f"
		ilcomp = select[je];
#line 413 "ztgevc.f"
	    }
#line 414 "ztgevc.f"
	    if (ilcomp) {
#line 415 "ztgevc.f"
		++ieig;

#line 417 "ztgevc.f"
		i__2 = je + je * s_dim1;
#line 417 "ztgevc.f"
		i__3 = je + je * p_dim1;
#line 417 "ztgevc.f"
		if ((d__2 = s[i__2].r, abs(d__2)) + (d__3 = d_imag(&s[je + je 
			* s_dim1]), abs(d__3)) <= safmin && (d__1 = p[i__3].r,
			 abs(d__1)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

#line 422 "ztgevc.f"
		    i__2 = *n;
#line 422 "ztgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 423 "ztgevc.f"
			i__3 = jr + ieig * vl_dim1;
#line 423 "ztgevc.f"
			vl[i__3].r = 0., vl[i__3].i = 0.;
#line 424 "ztgevc.f"
/* L50: */
#line 424 "ztgevc.f"
		    }
#line 425 "ztgevc.f"
		    i__2 = ieig + ieig * vl_dim1;
#line 425 "ztgevc.f"
		    vl[i__2].r = 1., vl[i__2].i = 0.;
#line 426 "ztgevc.f"
		    goto L140;
#line 427 "ztgevc.f"
		}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */
/*                   H */
/*                 y  ( a A - b B ) = 0 */

/* Computing MAX */
#line 434 "ztgevc.f"
		i__2 = je + je * s_dim1;
#line 434 "ztgevc.f"
		i__3 = je + je * p_dim1;
#line 434 "ztgevc.f"
		d__4 = ((d__2 = s[i__2].r, abs(d__2)) + (d__3 = d_imag(&s[je 
			+ je * s_dim1]), abs(d__3))) * ascale, d__5 = (d__1 = 
			p[i__3].r, abs(d__1)) * bscale, d__4 = max(d__4,d__5);
#line 434 "ztgevc.f"
		temp = 1. / max(d__4,safmin);
#line 436 "ztgevc.f"
		i__2 = je + je * s_dim1;
#line 436 "ztgevc.f"
		z__2.r = temp * s[i__2].r, z__2.i = temp * s[i__2].i;
#line 436 "ztgevc.f"
		z__1.r = ascale * z__2.r, z__1.i = ascale * z__2.i;
#line 436 "ztgevc.f"
		salpha.r = z__1.r, salpha.i = z__1.i;
#line 437 "ztgevc.f"
		i__2 = je + je * p_dim1;
#line 437 "ztgevc.f"
		sbeta = temp * p[i__2].r * bscale;
#line 438 "ztgevc.f"
		acoeff = sbeta * ascale;
#line 439 "ztgevc.f"
		z__1.r = bscale * salpha.r, z__1.i = bscale * salpha.i;
#line 439 "ztgevc.f"
		bcoeff.r = z__1.r, bcoeff.i = z__1.i;

/*              Scale to avoid underflow */

#line 443 "ztgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
#line 444 "ztgevc.f"
		lsb = (d__1 = salpha.r, abs(d__1)) + (d__2 = d_imag(&salpha), 
			abs(d__2)) >= safmin && (d__3 = bcoeff.r, abs(d__3)) 
			+ (d__4 = d_imag(&bcoeff), abs(d__4)) < small;

#line 447 "ztgevc.f"
		scale = 1.;
#line 448 "ztgevc.f"
		if (lsa) {
#line 448 "ztgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 448 "ztgevc.f"
		}
#line 450 "ztgevc.f"
		if (lsb) {
/* Computing MAX */
#line 450 "ztgevc.f"
		    d__3 = scale, d__4 = small / ((d__1 = salpha.r, abs(d__1))
			     + (d__2 = d_imag(&salpha), abs(d__2))) * min(
			    bnorm,big);
#line 450 "ztgevc.f"
		    scale = max(d__3,d__4);
#line 450 "ztgevc.f"
		}
#line 453 "ztgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 454 "ztgevc.f"
		    d__5 = 1., d__6 = abs(acoeff), d__5 = max(d__5,d__6), 
			    d__6 = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = 
			    d_imag(&bcoeff), abs(d__2));
#line 454 "ztgevc.f"
		    d__3 = scale, d__4 = 1. / (safmin * max(d__5,d__6));
#line 454 "ztgevc.f"
		    scale = min(d__3,d__4);
#line 457 "ztgevc.f"
		    if (lsa) {
#line 458 "ztgevc.f"
			acoeff = ascale * (scale * sbeta);
#line 459 "ztgevc.f"
		    } else {
#line 460 "ztgevc.f"
			acoeff = scale * acoeff;
#line 461 "ztgevc.f"
		    }
#line 462 "ztgevc.f"
		    if (lsb) {
#line 463 "ztgevc.f"
			z__2.r = scale * salpha.r, z__2.i = scale * salpha.i;
#line 463 "ztgevc.f"
			z__1.r = bscale * z__2.r, z__1.i = bscale * z__2.i;
#line 463 "ztgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 464 "ztgevc.f"
		    } else {
#line 465 "ztgevc.f"
			z__1.r = scale * bcoeff.r, z__1.i = scale * bcoeff.i;
#line 465 "ztgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 466 "ztgevc.f"
		    }
#line 467 "ztgevc.f"
		}

#line 469 "ztgevc.f"
		acoefa = abs(acoeff);
#line 470 "ztgevc.f"
		bcoefa = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = d_imag(&
			bcoeff), abs(d__2));
#line 471 "ztgevc.f"
		xmax = 1.;
#line 472 "ztgevc.f"
		i__2 = *n;
#line 472 "ztgevc.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 473 "ztgevc.f"
		    i__3 = jr;
#line 473 "ztgevc.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 474 "ztgevc.f"
/* L60: */
#line 474 "ztgevc.f"
		}
#line 475 "ztgevc.f"
		i__2 = je;
#line 475 "ztgevc.f"
		work[i__2].r = 1., work[i__2].i = 0.;
/* Computing MAX */
#line 476 "ztgevc.f"
		d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, 
			d__1 = max(d__1,d__2);
#line 476 "ztgevc.f"
		dmin__ = max(d__1,safmin);

/*                                              H */
/*              Triangular solve of  (a A - b B)  y = 0 */

/*                                      H */
/*              (rowwise in  (a A - b B) , or columnwise in a A - b B) */

#line 484 "ztgevc.f"
		i__2 = *n;
#line 484 "ztgevc.f"
		for (j = je + 1; j <= i__2; ++j) {

/*                 Compute */
/*                       j-1 */
/*                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k) */
/*                       k=je */
/*                 (Scale if necessary) */

#line 492 "ztgevc.f"
		    temp = 1. / xmax;
#line 493 "ztgevc.f"
		    if (acoefa * rwork[j] + bcoefa * rwork[*n + j] > bignum * 
			    temp) {
#line 495 "ztgevc.f"
			i__3 = j - 1;
#line 495 "ztgevc.f"
			for (jr = je; jr <= i__3; ++jr) {
#line 496 "ztgevc.f"
			    i__4 = jr;
#line 496 "ztgevc.f"
			    i__5 = jr;
#line 496 "ztgevc.f"
			    z__1.r = temp * work[i__5].r, z__1.i = temp * 
				    work[i__5].i;
#line 496 "ztgevc.f"
			    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 497 "ztgevc.f"
/* L70: */
#line 497 "ztgevc.f"
			}
#line 498 "ztgevc.f"
			xmax = 1.;
#line 499 "ztgevc.f"
		    }
#line 500 "ztgevc.f"
		    suma.r = 0., suma.i = 0.;
#line 501 "ztgevc.f"
		    sumb.r = 0., sumb.i = 0.;

#line 503 "ztgevc.f"
		    i__3 = j - 1;
#line 503 "ztgevc.f"
		    for (jr = je; jr <= i__3; ++jr) {
#line 504 "ztgevc.f"
			d_cnjg(&z__3, &s[jr + j * s_dim1]);
#line 504 "ztgevc.f"
			i__4 = jr;
#line 504 "ztgevc.f"
			z__2.r = z__3.r * work[i__4].r - z__3.i * work[i__4]
				.i, z__2.i = z__3.r * work[i__4].i + z__3.i * 
				work[i__4].r;
#line 504 "ztgevc.f"
			z__1.r = suma.r + z__2.r, z__1.i = suma.i + z__2.i;
#line 504 "ztgevc.f"
			suma.r = z__1.r, suma.i = z__1.i;
#line 505 "ztgevc.f"
			d_cnjg(&z__3, &p[jr + j * p_dim1]);
#line 505 "ztgevc.f"
			i__4 = jr;
#line 505 "ztgevc.f"
			z__2.r = z__3.r * work[i__4].r - z__3.i * work[i__4]
				.i, z__2.i = z__3.r * work[i__4].i + z__3.i * 
				work[i__4].r;
#line 505 "ztgevc.f"
			z__1.r = sumb.r + z__2.r, z__1.i = sumb.i + z__2.i;
#line 505 "ztgevc.f"
			sumb.r = z__1.r, sumb.i = z__1.i;
#line 506 "ztgevc.f"
/* L80: */
#line 506 "ztgevc.f"
		    }
#line 507 "ztgevc.f"
		    z__2.r = acoeff * suma.r, z__2.i = acoeff * suma.i;
#line 507 "ztgevc.f"
		    d_cnjg(&z__4, &bcoeff);
#line 507 "ztgevc.f"
		    z__3.r = z__4.r * sumb.r - z__4.i * sumb.i, z__3.i = 
			    z__4.r * sumb.i + z__4.i * sumb.r;
#line 507 "ztgevc.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 507 "ztgevc.f"
		    sum.r = z__1.r, sum.i = z__1.i;

/*                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) ) */

/*                 with scaling and perturbation of the denominator */

#line 513 "ztgevc.f"
		    i__3 = j + j * s_dim1;
#line 513 "ztgevc.f"
		    z__3.r = acoeff * s[i__3].r, z__3.i = acoeff * s[i__3].i;
#line 513 "ztgevc.f"
		    i__4 = j + j * p_dim1;
#line 513 "ztgevc.f"
		    z__4.r = bcoeff.r * p[i__4].r - bcoeff.i * p[i__4].i, 
			    z__4.i = bcoeff.r * p[i__4].i + bcoeff.i * p[i__4]
			    .r;
#line 513 "ztgevc.f"
		    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 513 "ztgevc.f"
		    d_cnjg(&z__1, &z__2);
#line 513 "ztgevc.f"
		    d__.r = z__1.r, d__.i = z__1.i;
#line 514 "ztgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) <= dmin__) {
#line 514 "ztgevc.f"
			z__1.r = dmin__, z__1.i = 0.;
#line 514 "ztgevc.f"
			d__.r = z__1.r, d__.i = z__1.i;
#line 514 "ztgevc.f"
		    }

#line 517 "ztgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) < 1.) {
#line 518 "ztgevc.f"
			if ((d__1 = sum.r, abs(d__1)) + (d__2 = d_imag(&sum), 
				abs(d__2)) >= bignum * ((d__3 = d__.r, abs(
				d__3)) + (d__4 = d_imag(&d__), abs(d__4)))) {
#line 519 "ztgevc.f"
			    temp = 1. / ((d__1 = sum.r, abs(d__1)) + (d__2 = 
				    d_imag(&sum), abs(d__2)));
#line 520 "ztgevc.f"
			    i__3 = j - 1;
#line 520 "ztgevc.f"
			    for (jr = je; jr <= i__3; ++jr) {
#line 521 "ztgevc.f"
				i__4 = jr;
#line 521 "ztgevc.f"
				i__5 = jr;
#line 521 "ztgevc.f"
				z__1.r = temp * work[i__5].r, z__1.i = temp * 
					work[i__5].i;
#line 521 "ztgevc.f"
				work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 522 "ztgevc.f"
/* L90: */
#line 522 "ztgevc.f"
			    }
#line 523 "ztgevc.f"
			    xmax = temp * xmax;
#line 524 "ztgevc.f"
			    z__1.r = temp * sum.r, z__1.i = temp * sum.i;
#line 524 "ztgevc.f"
			    sum.r = z__1.r, sum.i = z__1.i;
#line 525 "ztgevc.f"
			}
#line 526 "ztgevc.f"
		    }
#line 527 "ztgevc.f"
		    i__3 = j;
#line 527 "ztgevc.f"
		    z__2.r = -sum.r, z__2.i = -sum.i;
#line 527 "ztgevc.f"
		    zladiv_(&z__1, &z__2, &d__);
#line 527 "ztgevc.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
/* Computing MAX */
#line 528 "ztgevc.f"
		    i__3 = j;
#line 528 "ztgevc.f"
		    d__3 = xmax, d__4 = (d__1 = work[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&work[j]), abs(d__2));
#line 528 "ztgevc.f"
		    xmax = max(d__3,d__4);
#line 529 "ztgevc.f"
/* L100: */
#line 529 "ztgevc.f"
		}

/*              Back transform eigenvector if HOWMNY='B'. */

#line 533 "ztgevc.f"
		if (ilback) {
#line 534 "ztgevc.f"
		    i__2 = *n + 1 - je;
#line 534 "ztgevc.f"
		    zgemv_("N", n, &i__2, &c_b2, &vl[je * vl_dim1 + 1], ldvl, 
			    &work[je], &c__1, &c_b1, &work[*n + 1], &c__1, (
			    ftnlen)1);
#line 536 "ztgevc.f"
		    isrc = 2;
#line 537 "ztgevc.f"
		    ibeg = 1;
#line 538 "ztgevc.f"
		} else {
#line 539 "ztgevc.f"
		    isrc = 1;
#line 540 "ztgevc.f"
		    ibeg = je;
#line 541 "ztgevc.f"
		}

/*              Copy and scale eigenvector into column of VL */

#line 545 "ztgevc.f"
		xmax = 0.;
#line 546 "ztgevc.f"
		i__2 = *n;
#line 546 "ztgevc.f"
		for (jr = ibeg; jr <= i__2; ++jr) {
/* Computing MAX */
#line 547 "ztgevc.f"
		    i__3 = (isrc - 1) * *n + jr;
#line 547 "ztgevc.f"
		    d__3 = xmax, d__4 = (d__1 = work[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&work[(isrc - 1) * *n + jr]), abs(
			    d__2));
#line 547 "ztgevc.f"
		    xmax = max(d__3,d__4);
#line 548 "ztgevc.f"
/* L110: */
#line 548 "ztgevc.f"
		}

#line 550 "ztgevc.f"
		if (xmax > safmin) {
#line 551 "ztgevc.f"
		    temp = 1. / xmax;
#line 552 "ztgevc.f"
		    i__2 = *n;
#line 552 "ztgevc.f"
		    for (jr = ibeg; jr <= i__2; ++jr) {
#line 553 "ztgevc.f"
			i__3 = jr + ieig * vl_dim1;
#line 553 "ztgevc.f"
			i__4 = (isrc - 1) * *n + jr;
#line 553 "ztgevc.f"
			z__1.r = temp * work[i__4].r, z__1.i = temp * work[
				i__4].i;
#line 553 "ztgevc.f"
			vl[i__3].r = z__1.r, vl[i__3].i = z__1.i;
#line 554 "ztgevc.f"
/* L120: */
#line 554 "ztgevc.f"
		    }
#line 555 "ztgevc.f"
		} else {
#line 556 "ztgevc.f"
		    ibeg = *n + 1;
#line 557 "ztgevc.f"
		}

#line 559 "ztgevc.f"
		i__2 = ibeg - 1;
#line 559 "ztgevc.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 560 "ztgevc.f"
		    i__3 = jr + ieig * vl_dim1;
#line 560 "ztgevc.f"
		    vl[i__3].r = 0., vl[i__3].i = 0.;
#line 561 "ztgevc.f"
/* L130: */
#line 561 "ztgevc.f"
		}

#line 563 "ztgevc.f"
	    }
#line 564 "ztgevc.f"
L140:
#line 564 "ztgevc.f"
	    ;
#line 564 "ztgevc.f"
	}
#line 565 "ztgevc.f"
    }

/*     Right eigenvectors */

#line 569 "ztgevc.f"
    if (compr) {
#line 570 "ztgevc.f"
	ieig = im + 1;

/*        Main loop over eigenvalues */

#line 574 "ztgevc.f"
	for (je = *n; je >= 1; --je) {
#line 575 "ztgevc.f"
	    if (ilall) {
#line 576 "ztgevc.f"
		ilcomp = TRUE_;
#line 577 "ztgevc.f"
	    } else {
#line 578 "ztgevc.f"
		ilcomp = select[je];
#line 579 "ztgevc.f"
	    }
#line 580 "ztgevc.f"
	    if (ilcomp) {
#line 581 "ztgevc.f"
		--ieig;

#line 583 "ztgevc.f"
		i__1 = je + je * s_dim1;
#line 583 "ztgevc.f"
		i__2 = je + je * p_dim1;
#line 583 "ztgevc.f"
		if ((d__2 = s[i__1].r, abs(d__2)) + (d__3 = d_imag(&s[je + je 
			* s_dim1]), abs(d__3)) <= safmin && (d__1 = p[i__2].r,
			 abs(d__1)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

#line 588 "ztgevc.f"
		    i__1 = *n;
#line 588 "ztgevc.f"
		    for (jr = 1; jr <= i__1; ++jr) {
#line 589 "ztgevc.f"
			i__2 = jr + ieig * vr_dim1;
#line 589 "ztgevc.f"
			vr[i__2].r = 0., vr[i__2].i = 0.;
#line 590 "ztgevc.f"
/* L150: */
#line 590 "ztgevc.f"
		    }
#line 591 "ztgevc.f"
		    i__1 = ieig + ieig * vr_dim1;
#line 591 "ztgevc.f"
		    vr[i__1].r = 1., vr[i__1].i = 0.;
#line 592 "ztgevc.f"
		    goto L250;
#line 593 "ztgevc.f"
		}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */

/*              ( a A - b B ) x  = 0 */

/* Computing MAX */
#line 600 "ztgevc.f"
		i__1 = je + je * s_dim1;
#line 600 "ztgevc.f"
		i__2 = je + je * p_dim1;
#line 600 "ztgevc.f"
		d__4 = ((d__2 = s[i__1].r, abs(d__2)) + (d__3 = d_imag(&s[je 
			+ je * s_dim1]), abs(d__3))) * ascale, d__5 = (d__1 = 
			p[i__2].r, abs(d__1)) * bscale, d__4 = max(d__4,d__5);
#line 600 "ztgevc.f"
		temp = 1. / max(d__4,safmin);
#line 602 "ztgevc.f"
		i__1 = je + je * s_dim1;
#line 602 "ztgevc.f"
		z__2.r = temp * s[i__1].r, z__2.i = temp * s[i__1].i;
#line 602 "ztgevc.f"
		z__1.r = ascale * z__2.r, z__1.i = ascale * z__2.i;
#line 602 "ztgevc.f"
		salpha.r = z__1.r, salpha.i = z__1.i;
#line 603 "ztgevc.f"
		i__1 = je + je * p_dim1;
#line 603 "ztgevc.f"
		sbeta = temp * p[i__1].r * bscale;
#line 604 "ztgevc.f"
		acoeff = sbeta * ascale;
#line 605 "ztgevc.f"
		z__1.r = bscale * salpha.r, z__1.i = bscale * salpha.i;
#line 605 "ztgevc.f"
		bcoeff.r = z__1.r, bcoeff.i = z__1.i;

/*              Scale to avoid underflow */

#line 609 "ztgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
#line 610 "ztgevc.f"
		lsb = (d__1 = salpha.r, abs(d__1)) + (d__2 = d_imag(&salpha), 
			abs(d__2)) >= safmin && (d__3 = bcoeff.r, abs(d__3)) 
			+ (d__4 = d_imag(&bcoeff), abs(d__4)) < small;

#line 613 "ztgevc.f"
		scale = 1.;
#line 614 "ztgevc.f"
		if (lsa) {
#line 614 "ztgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 614 "ztgevc.f"
		}
#line 616 "ztgevc.f"
		if (lsb) {
/* Computing MAX */
#line 616 "ztgevc.f"
		    d__3 = scale, d__4 = small / ((d__1 = salpha.r, abs(d__1))
			     + (d__2 = d_imag(&salpha), abs(d__2))) * min(
			    bnorm,big);
#line 616 "ztgevc.f"
		    scale = max(d__3,d__4);
#line 616 "ztgevc.f"
		}
#line 619 "ztgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 620 "ztgevc.f"
		    d__5 = 1., d__6 = abs(acoeff), d__5 = max(d__5,d__6), 
			    d__6 = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = 
			    d_imag(&bcoeff), abs(d__2));
#line 620 "ztgevc.f"
		    d__3 = scale, d__4 = 1. / (safmin * max(d__5,d__6));
#line 620 "ztgevc.f"
		    scale = min(d__3,d__4);
#line 623 "ztgevc.f"
		    if (lsa) {
#line 624 "ztgevc.f"
			acoeff = ascale * (scale * sbeta);
#line 625 "ztgevc.f"
		    } else {
#line 626 "ztgevc.f"
			acoeff = scale * acoeff;
#line 627 "ztgevc.f"
		    }
#line 628 "ztgevc.f"
		    if (lsb) {
#line 629 "ztgevc.f"
			z__2.r = scale * salpha.r, z__2.i = scale * salpha.i;
#line 629 "ztgevc.f"
			z__1.r = bscale * z__2.r, z__1.i = bscale * z__2.i;
#line 629 "ztgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 630 "ztgevc.f"
		    } else {
#line 631 "ztgevc.f"
			z__1.r = scale * bcoeff.r, z__1.i = scale * bcoeff.i;
#line 631 "ztgevc.f"
			bcoeff.r = z__1.r, bcoeff.i = z__1.i;
#line 632 "ztgevc.f"
		    }
#line 633 "ztgevc.f"
		}

#line 635 "ztgevc.f"
		acoefa = abs(acoeff);
#line 636 "ztgevc.f"
		bcoefa = (d__1 = bcoeff.r, abs(d__1)) + (d__2 = d_imag(&
			bcoeff), abs(d__2));
#line 637 "ztgevc.f"
		xmax = 1.;
#line 638 "ztgevc.f"
		i__1 = *n;
#line 638 "ztgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 639 "ztgevc.f"
		    i__2 = jr;
#line 639 "ztgevc.f"
		    work[i__2].r = 0., work[i__2].i = 0.;
#line 640 "ztgevc.f"
/* L160: */
#line 640 "ztgevc.f"
		}
#line 641 "ztgevc.f"
		i__1 = je;
#line 641 "ztgevc.f"
		work[i__1].r = 1., work[i__1].i = 0.;
/* Computing MAX */
#line 642 "ztgevc.f"
		d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, 
			d__1 = max(d__1,d__2);
#line 642 "ztgevc.f"
		dmin__ = max(d__1,safmin);

/*              Triangular solve of  (a A - b B) x = 0  (columnwise) */

/*              WORK(1:j-1) contains sums w, */
/*              WORK(j+1:JE) contains x */

#line 649 "ztgevc.f"
		i__1 = je - 1;
#line 649 "ztgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 650 "ztgevc.f"
		    i__2 = jr;
#line 650 "ztgevc.f"
		    i__3 = jr + je * s_dim1;
#line 650 "ztgevc.f"
		    z__2.r = acoeff * s[i__3].r, z__2.i = acoeff * s[i__3].i;
#line 650 "ztgevc.f"
		    i__4 = jr + je * p_dim1;
#line 650 "ztgevc.f"
		    z__3.r = bcoeff.r * p[i__4].r - bcoeff.i * p[i__4].i, 
			    z__3.i = bcoeff.r * p[i__4].i + bcoeff.i * p[i__4]
			    .r;
#line 650 "ztgevc.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 650 "ztgevc.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 651 "ztgevc.f"
/* L170: */
#line 651 "ztgevc.f"
		}
#line 652 "ztgevc.f"
		i__1 = je;
#line 652 "ztgevc.f"
		work[i__1].r = 1., work[i__1].i = 0.;

#line 654 "ztgevc.f"
		for (j = je - 1; j >= 1; --j) {

/*                 Form x(j) := - w(j) / d */
/*                 with scaling and perturbation of the denominator */

#line 659 "ztgevc.f"
		    i__1 = j + j * s_dim1;
#line 659 "ztgevc.f"
		    z__2.r = acoeff * s[i__1].r, z__2.i = acoeff * s[i__1].i;
#line 659 "ztgevc.f"
		    i__2 = j + j * p_dim1;
#line 659 "ztgevc.f"
		    z__3.r = bcoeff.r * p[i__2].r - bcoeff.i * p[i__2].i, 
			    z__3.i = bcoeff.r * p[i__2].i + bcoeff.i * p[i__2]
			    .r;
#line 659 "ztgevc.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 659 "ztgevc.f"
		    d__.r = z__1.r, d__.i = z__1.i;
#line 660 "ztgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) <= dmin__) {
#line 660 "ztgevc.f"
			z__1.r = dmin__, z__1.i = 0.;
#line 660 "ztgevc.f"
			d__.r = z__1.r, d__.i = z__1.i;
#line 660 "ztgevc.f"
		    }

#line 663 "ztgevc.f"
		    if ((d__1 = d__.r, abs(d__1)) + (d__2 = d_imag(&d__), abs(
			    d__2)) < 1.) {
#line 664 "ztgevc.f"
			i__1 = j;
#line 664 "ztgevc.f"
			if ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(
				&work[j]), abs(d__2)) >= bignum * ((d__3 = 
				d__.r, abs(d__3)) + (d__4 = d_imag(&d__), abs(
				d__4)))) {
#line 665 "ztgevc.f"
			    i__1 = j;
#line 665 "ztgevc.f"
			    temp = 1. / ((d__1 = work[i__1].r, abs(d__1)) + (
				    d__2 = d_imag(&work[j]), abs(d__2)));
#line 666 "ztgevc.f"
			    i__1 = je;
#line 666 "ztgevc.f"
			    for (jr = 1; jr <= i__1; ++jr) {
#line 667 "ztgevc.f"
				i__2 = jr;
#line 667 "ztgevc.f"
				i__3 = jr;
#line 667 "ztgevc.f"
				z__1.r = temp * work[i__3].r, z__1.i = temp * 
					work[i__3].i;
#line 667 "ztgevc.f"
				work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 668 "ztgevc.f"
/* L180: */
#line 668 "ztgevc.f"
			    }
#line 669 "ztgevc.f"
			}
#line 670 "ztgevc.f"
		    }

#line 672 "ztgevc.f"
		    i__1 = j;
#line 672 "ztgevc.f"
		    i__2 = j;
#line 672 "ztgevc.f"
		    z__2.r = -work[i__2].r, z__2.i = -work[i__2].i;
#line 672 "ztgevc.f"
		    zladiv_(&z__1, &z__2, &d__);
#line 672 "ztgevc.f"
		    work[i__1].r = z__1.r, work[i__1].i = z__1.i;

#line 674 "ztgevc.f"
		    if (j > 1) {

/*                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling */

#line 678 "ztgevc.f"
			i__1 = j;
#line 678 "ztgevc.f"
			if ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(
				&work[j]), abs(d__2)) > 1.) {
#line 679 "ztgevc.f"
			    i__1 = j;
#line 679 "ztgevc.f"
			    temp = 1. / ((d__1 = work[i__1].r, abs(d__1)) + (
				    d__2 = d_imag(&work[j]), abs(d__2)));
#line 680 "ztgevc.f"
			    if (acoefa * rwork[j] + bcoefa * rwork[*n + j] >= 
				    bignum * temp) {
#line 682 "ztgevc.f"
				i__1 = je;
#line 682 "ztgevc.f"
				for (jr = 1; jr <= i__1; ++jr) {
#line 683 "ztgevc.f"
				    i__2 = jr;
#line 683 "ztgevc.f"
				    i__3 = jr;
#line 683 "ztgevc.f"
				    z__1.r = temp * work[i__3].r, z__1.i = 
					    temp * work[i__3].i;
#line 683 "ztgevc.f"
				    work[i__2].r = z__1.r, work[i__2].i = 
					    z__1.i;
#line 684 "ztgevc.f"
/* L190: */
#line 684 "ztgevc.f"
				}
#line 685 "ztgevc.f"
			    }
#line 686 "ztgevc.f"
			}

#line 688 "ztgevc.f"
			i__1 = j;
#line 688 "ztgevc.f"
			z__1.r = acoeff * work[i__1].r, z__1.i = acoeff * 
				work[i__1].i;
#line 688 "ztgevc.f"
			ca.r = z__1.r, ca.i = z__1.i;
#line 689 "ztgevc.f"
			i__1 = j;
#line 689 "ztgevc.f"
			z__1.r = bcoeff.r * work[i__1].r - bcoeff.i * work[
				i__1].i, z__1.i = bcoeff.r * work[i__1].i + 
				bcoeff.i * work[i__1].r;
#line 689 "ztgevc.f"
			cb.r = z__1.r, cb.i = z__1.i;
#line 690 "ztgevc.f"
			i__1 = j - 1;
#line 690 "ztgevc.f"
			for (jr = 1; jr <= i__1; ++jr) {
#line 691 "ztgevc.f"
			    i__2 = jr;
#line 691 "ztgevc.f"
			    i__3 = jr;
#line 691 "ztgevc.f"
			    i__4 = jr + j * s_dim1;
#line 691 "ztgevc.f"
			    z__3.r = ca.r * s[i__4].r - ca.i * s[i__4].i, 
				    z__3.i = ca.r * s[i__4].i + ca.i * s[i__4]
				    .r;
#line 691 "ztgevc.f"
			    z__2.r = work[i__3].r + z__3.r, z__2.i = work[
				    i__3].i + z__3.i;
#line 691 "ztgevc.f"
			    i__5 = jr + j * p_dim1;
#line 691 "ztgevc.f"
			    z__4.r = cb.r * p[i__5].r - cb.i * p[i__5].i, 
				    z__4.i = cb.r * p[i__5].i + cb.i * p[i__5]
				    .r;
#line 691 "ztgevc.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 691 "ztgevc.f"
			    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 693 "ztgevc.f"
/* L200: */
#line 693 "ztgevc.f"
			}
#line 694 "ztgevc.f"
		    }
#line 695 "ztgevc.f"
/* L210: */
#line 695 "ztgevc.f"
		}

/*              Back transform eigenvector if HOWMNY='B'. */

#line 699 "ztgevc.f"
		if (ilback) {
#line 700 "ztgevc.f"
		    zgemv_("N", n, &je, &c_b2, &vr[vr_offset], ldvr, &work[1],
			     &c__1, &c_b1, &work[*n + 1], &c__1, (ftnlen)1);
#line 702 "ztgevc.f"
		    isrc = 2;
#line 703 "ztgevc.f"
		    iend = *n;
#line 704 "ztgevc.f"
		} else {
#line 705 "ztgevc.f"
		    isrc = 1;
#line 706 "ztgevc.f"
		    iend = je;
#line 707 "ztgevc.f"
		}

/*              Copy and scale eigenvector into column of VR */

#line 711 "ztgevc.f"
		xmax = 0.;
#line 712 "ztgevc.f"
		i__1 = iend;
#line 712 "ztgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
/* Computing MAX */
#line 713 "ztgevc.f"
		    i__2 = (isrc - 1) * *n + jr;
#line 713 "ztgevc.f"
		    d__3 = xmax, d__4 = (d__1 = work[i__2].r, abs(d__1)) + (
			    d__2 = d_imag(&work[(isrc - 1) * *n + jr]), abs(
			    d__2));
#line 713 "ztgevc.f"
		    xmax = max(d__3,d__4);
#line 714 "ztgevc.f"
/* L220: */
#line 714 "ztgevc.f"
		}

#line 716 "ztgevc.f"
		if (xmax > safmin) {
#line 717 "ztgevc.f"
		    temp = 1. / xmax;
#line 718 "ztgevc.f"
		    i__1 = iend;
#line 718 "ztgevc.f"
		    for (jr = 1; jr <= i__1; ++jr) {
#line 719 "ztgevc.f"
			i__2 = jr + ieig * vr_dim1;
#line 719 "ztgevc.f"
			i__3 = (isrc - 1) * *n + jr;
#line 719 "ztgevc.f"
			z__1.r = temp * work[i__3].r, z__1.i = temp * work[
				i__3].i;
#line 719 "ztgevc.f"
			vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
#line 720 "ztgevc.f"
/* L230: */
#line 720 "ztgevc.f"
		    }
#line 721 "ztgevc.f"
		} else {
#line 722 "ztgevc.f"
		    iend = 0;
#line 723 "ztgevc.f"
		}

#line 725 "ztgevc.f"
		i__1 = *n;
#line 725 "ztgevc.f"
		for (jr = iend + 1; jr <= i__1; ++jr) {
#line 726 "ztgevc.f"
		    i__2 = jr + ieig * vr_dim1;
#line 726 "ztgevc.f"
		    vr[i__2].r = 0., vr[i__2].i = 0.;
#line 727 "ztgevc.f"
/* L240: */
#line 727 "ztgevc.f"
		}

#line 729 "ztgevc.f"
	    }
#line 730 "ztgevc.f"
L250:
#line 730 "ztgevc.f"
	    ;
#line 730 "ztgevc.f"
	}
#line 731 "ztgevc.f"
    }

#line 733 "ztgevc.f"
    return 0;

/*     End of ZTGEVC */

} /* ztgevc_ */

