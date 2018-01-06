#line 1 "ctrevc3.f"
/* ctrevc3.f -- translated by f2c (version 20100827).
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

#line 1 "ctrevc3.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CTREVC3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTREVC3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                           LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       REAL   RWORK( * ) */
/*       COMPLEX         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTREVC3 computes some or all of the right and/or left eigenvectors of */
/* > a complex upper triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR. */
/* > */
/* > The right eigenvector x and the left eigenvector y of T corresponding */
/* > to an eigenvalue w are defined by: */
/* > */
/* >              T*x = w*x,     (y**H)*T = w*(y**H) */
/* > */
/* > where y**H denotes the conjugate transpose of the vector y. */
/* > The eigenvalues are not input to this routine, but are read directly */
/* > from the diagonal of T. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an */
/* > input matrix. If Q is the unitary factor that reduces a matrix A to */
/* > Schur form T, then Q*X and Q*Y are the matrices of right and left */
/* > eigenvectors of A. */
/* > */
/* > This uses a Level 3 BLAS version of the back transformation. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'R':  compute right eigenvectors only; */
/* >          = 'L':  compute left eigenvectors only; */
/* >          = 'B':  compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* >          HOWMNY is CHARACTER*1 */
/* >          = 'A':  compute all right and/or left eigenvectors; */
/* >          = 'B':  compute all right and/or left eigenvectors, */
/* >                  backtransformed using the matrices supplied in */
/* >                  VR and/or VL; */
/* >          = 'S':  compute selected right and/or left eigenvectors, */
/* >                  as indicated by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
/* >          computed. */
/* >          The eigenvector corresponding to the j-th eigenvalue is */
/* >          computed if SELECT(j) = .TRUE.. */
/* >          Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,N) */
/* >          The upper triangular matrix T.  T is modified, but restored */
/* >          on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is COMPLEX array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Q of */
/* >          Schur vectors returned by CHSEQR). */
/* >          On exit, if SIDE = 'L' or 'B', VL contains: */
/* >          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; */
/* >          if HOWMNY = 'B', the matrix Q*Y; */
/* >          if HOWMNY = 'S', the left eigenvectors of T specified by */
/* >                           SELECT, stored consecutively in the columns */
/* >                           of VL, in the same order as their */
/* >                           eigenvalues. */
/* >          Not referenced if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL. */
/* >          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is COMPLEX array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Q of */
/* >          Schur vectors returned by CHSEQR). */
/* >          On exit, if SIDE = 'R' or 'B', VR contains: */
/* >          if HOWMNY = 'A', the matrix X of right eigenvectors of T; */
/* >          if HOWMNY = 'B', the matrix Q*X; */
/* >          if HOWMNY = 'S', the right eigenvectors of T specified by */
/* >                           SELECT, stored consecutively in the columns */
/* >                           of VR, in the same order as their */
/* >                           eigenvalues. */
/* >          Not referenced if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR. */
/* >          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N. */
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
/* >          used to store the eigenvectors. */
/* >          If HOWMNY = 'A' or 'B', M is set to N. */
/* >          Each selected eigenvector occupies one column. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of array WORK. LWORK >= max(1,2*N). */
/* >          For optimum performance, LWORK >= N + 2*N*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (LRWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of array RWORK. LRWORK >= max(1,N). */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the RWORK array, returns */
/* >          this value as the first entry of the RWORK array, and no error */
/* >          message related to LRWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/*  @generated from ztrevc3.f, fortran z -> c, Tue Apr 19 01:47:44 2016 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The algorithm used in this program is basically backward (forward) */
/* >  substitution, with scaling to make the the code robust against */
/* >  possible overflow. */
/* > */
/* >  Each eigenvector is normalized so that the element of largest */
/* >  magnitude has magnitude 1; here the magnitude of a complex number */
/* >  (x,y) is taken to be |x| + |y|. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctrevc3_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer 
	*m, doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *info, ftnlen side_len, ftnlen howmny_len)
{
    /* System generated locals */
    address a__1[2];
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2[2], i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, nb, ii, ki, is, iv;
    static doublereal ulp;
    static logical allv;
    static doublereal unfl, ovfl, smin;
    static logical over;
    static doublereal scale;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal remax;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical leftv, bothv, somev;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), claset_(char *, integer *, integer *,
	     doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int clatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);
    static logical rightv;
    static integer maxwrk;
    static doublereal smlnum;
    static logical lquery;


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

/*     Decode and test the input parameters */

#line 306 "ctrevc3.f"
    /* Parameter adjustments */
#line 306 "ctrevc3.f"
    --select;
#line 306 "ctrevc3.f"
    t_dim1 = *ldt;
#line 306 "ctrevc3.f"
    t_offset = 1 + t_dim1;
#line 306 "ctrevc3.f"
    t -= t_offset;
#line 306 "ctrevc3.f"
    vl_dim1 = *ldvl;
#line 306 "ctrevc3.f"
    vl_offset = 1 + vl_dim1;
#line 306 "ctrevc3.f"
    vl -= vl_offset;
#line 306 "ctrevc3.f"
    vr_dim1 = *ldvr;
#line 306 "ctrevc3.f"
    vr_offset = 1 + vr_dim1;
#line 306 "ctrevc3.f"
    vr -= vr_offset;
#line 306 "ctrevc3.f"
    --work;
#line 306 "ctrevc3.f"
    --rwork;
#line 306 "ctrevc3.f"

#line 306 "ctrevc3.f"
    /* Function Body */
#line 306 "ctrevc3.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 307 "ctrevc3.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 308 "ctrevc3.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 310 "ctrevc3.f"
    allv = lsame_(howmny, "A", (ftnlen)1, (ftnlen)1);
#line 311 "ctrevc3.f"
    over = lsame_(howmny, "B", (ftnlen)1, (ftnlen)1);
#line 312 "ctrevc3.f"
    somev = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of columns required to store the selected */
/*     eigenvectors. */

#line 317 "ctrevc3.f"
    if (somev) {
#line 318 "ctrevc3.f"
	*m = 0;
#line 319 "ctrevc3.f"
	i__1 = *n;
#line 319 "ctrevc3.f"
	for (j = 1; j <= i__1; ++j) {
#line 320 "ctrevc3.f"
	    if (select[j]) {
#line 320 "ctrevc3.f"
		++(*m);
#line 320 "ctrevc3.f"
	    }
#line 322 "ctrevc3.f"
/* L10: */
#line 322 "ctrevc3.f"
	}
#line 323 "ctrevc3.f"
    } else {
#line 324 "ctrevc3.f"
	*m = *n;
#line 325 "ctrevc3.f"
    }

#line 327 "ctrevc3.f"
    *info = 0;
/* Writing concatenation */
#line 328 "ctrevc3.f"
    i__2[0] = 1, a__1[0] = side;
#line 328 "ctrevc3.f"
    i__2[1] = 1, a__1[1] = howmny;
#line 328 "ctrevc3.f"
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
#line 328 "ctrevc3.f"
    nb = ilaenv_(&c__1, "CTREVC", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
#line 329 "ctrevc3.f"
    maxwrk = *n + (*n << 1) * nb;
#line 330 "ctrevc3.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 331 "ctrevc3.f"
    rwork[1] = (doublereal) (*n);
#line 332 "ctrevc3.f"
    lquery = *lwork == -1 || *lrwork == -1;
#line 333 "ctrevc3.f"
    if (! rightv && ! leftv) {
#line 334 "ctrevc3.f"
	*info = -1;
#line 335 "ctrevc3.f"
    } else if (! allv && ! over && ! somev) {
#line 336 "ctrevc3.f"
	*info = -2;
#line 337 "ctrevc3.f"
    } else if (*n < 0) {
#line 338 "ctrevc3.f"
	*info = -4;
#line 339 "ctrevc3.f"
    } else if (*ldt < max(1,*n)) {
#line 340 "ctrevc3.f"
	*info = -6;
#line 341 "ctrevc3.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 342 "ctrevc3.f"
	*info = -8;
#line 343 "ctrevc3.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 344 "ctrevc3.f"
	*info = -10;
#line 345 "ctrevc3.f"
    } else if (*mm < *m) {
#line 346 "ctrevc3.f"
	*info = -11;
#line 347 "ctrevc3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 347 "ctrevc3.f"
	i__1 = 1, i__3 = *n << 1;
#line 347 "ctrevc3.f"
	if (*lwork < max(i__1,i__3) && ! lquery) {
#line 348 "ctrevc3.f"
	    *info = -14;
#line 349 "ctrevc3.f"
	} else if (*lrwork < max(1,*n) && ! lquery) {
#line 350 "ctrevc3.f"
	    *info = -16;
#line 351 "ctrevc3.f"
	}
#line 351 "ctrevc3.f"
    }
#line 352 "ctrevc3.f"
    if (*info != 0) {
#line 353 "ctrevc3.f"
	i__1 = -(*info);
#line 353 "ctrevc3.f"
	xerbla_("CTREVC3", &i__1, (ftnlen)7);
#line 354 "ctrevc3.f"
	return 0;
#line 355 "ctrevc3.f"
    } else if (lquery) {
#line 356 "ctrevc3.f"
	return 0;
#line 357 "ctrevc3.f"
    }

/*     Quick return if possible. */

#line 361 "ctrevc3.f"
    if (*n == 0) {
#line 361 "ctrevc3.f"
	return 0;
#line 361 "ctrevc3.f"
    }

/*     Use blocked version of back-transformation if sufficient workspace. */
/*     Zero-out the workspace to avoid potential NaN propagation. */

#line 367 "ctrevc3.f"
    if (over && *lwork >= *n + (*n << 4)) {
#line 368 "ctrevc3.f"
	nb = (*lwork - *n) / (*n << 1);
#line 369 "ctrevc3.f"
	nb = min(nb,128);
#line 370 "ctrevc3.f"
	i__1 = (nb << 1) + 1;
#line 370 "ctrevc3.f"
	claset_("F", n, &i__1, &c_b1, &c_b1, &work[1], n, (ftnlen)1);
#line 371 "ctrevc3.f"
    } else {
#line 372 "ctrevc3.f"
	nb = 1;
#line 373 "ctrevc3.f"
    }

/*     Set the constants to control overflow. */

#line 377 "ctrevc3.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);
#line 378 "ctrevc3.f"
    ovfl = 1. / unfl;
#line 379 "ctrevc3.f"
    slabad_(&unfl, &ovfl);
#line 380 "ctrevc3.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 381 "ctrevc3.f"
    smlnum = unfl * (*n / ulp);

/*     Store the diagonal elements of T in working array WORK. */

#line 385 "ctrevc3.f"
    i__1 = *n;
#line 385 "ctrevc3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "ctrevc3.f"
	i__3 = i__;
#line 386 "ctrevc3.f"
	i__4 = i__ + i__ * t_dim1;
#line 386 "ctrevc3.f"
	work[i__3].r = t[i__4].r, work[i__3].i = t[i__4].i;
#line 387 "ctrevc3.f"
/* L20: */
#line 387 "ctrevc3.f"
    }

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 392 "ctrevc3.f"
    rwork[1] = 0.;
#line 393 "ctrevc3.f"
    i__1 = *n;
#line 393 "ctrevc3.f"
    for (j = 2; j <= i__1; ++j) {
#line 394 "ctrevc3.f"
	i__3 = j - 1;
#line 394 "ctrevc3.f"
	rwork[j] = scasum_(&i__3, &t[j * t_dim1 + 1], &c__1);
#line 395 "ctrevc3.f"
/* L30: */
#line 395 "ctrevc3.f"
    }

#line 397 "ctrevc3.f"
    if (rightv) {

/*        ============================================================ */
/*        Compute right eigenvectors. */

/*        IV is index of column in current block. */
/*        Non-blocked version always uses IV=NB=1; */
/*        blocked     version starts with IV=NB, goes down to 1. */
/*        (Note the "0-th" column is used to store the original diagonal.) */
#line 406 "ctrevc3.f"
	iv = nb;
#line 407 "ctrevc3.f"
	is = *m;
#line 408 "ctrevc3.f"
	for (ki = *n; ki >= 1; --ki) {
#line 409 "ctrevc3.f"
	    if (somev) {
#line 410 "ctrevc3.f"
		if (! select[ki]) {
#line 410 "ctrevc3.f"
		    goto L80;
#line 410 "ctrevc3.f"
		}
#line 412 "ctrevc3.f"
	    }
/* Computing MAX */
#line 413 "ctrevc3.f"
	    i__1 = ki + ki * t_dim1;
#line 413 "ctrevc3.f"
	    d__3 = ulp * ((d__1 = t[i__1].r, abs(d__1)) + (d__2 = d_imag(&t[
		    ki + ki * t_dim1]), abs(d__2)));
#line 413 "ctrevc3.f"
	    smin = max(d__3,smlnum);

/*           -------------------------------------------------------- */
/*           Complex right eigenvector */

#line 418 "ctrevc3.f"
	    i__1 = ki + iv * *n;
#line 418 "ctrevc3.f"
	    work[i__1].r = 1., work[i__1].i = 0.;

/*           Form right-hand side. */

#line 422 "ctrevc3.f"
	    i__1 = ki - 1;
#line 422 "ctrevc3.f"
	    for (k = 1; k <= i__1; ++k) {
#line 423 "ctrevc3.f"
		i__3 = k + iv * *n;
#line 423 "ctrevc3.f"
		i__4 = k + ki * t_dim1;
#line 423 "ctrevc3.f"
		z__1.r = -t[i__4].r, z__1.i = -t[i__4].i;
#line 423 "ctrevc3.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 424 "ctrevc3.f"
/* L40: */
#line 424 "ctrevc3.f"
	    }

/*           Solve upper triangular system: */
/*           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK. */

#line 429 "ctrevc3.f"
	    i__1 = ki - 1;
#line 429 "ctrevc3.f"
	    for (k = 1; k <= i__1; ++k) {
#line 430 "ctrevc3.f"
		i__3 = k + k * t_dim1;
#line 430 "ctrevc3.f"
		i__4 = k + k * t_dim1;
#line 430 "ctrevc3.f"
		i__5 = ki + ki * t_dim1;
#line 430 "ctrevc3.f"
		z__1.r = t[i__4].r - t[i__5].r, z__1.i = t[i__4].i - t[i__5]
			.i;
#line 430 "ctrevc3.f"
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
#line 431 "ctrevc3.f"
		i__3 = k + k * t_dim1;
#line 431 "ctrevc3.f"
		if ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * 
			t_dim1]), abs(d__2)) < smin) {
#line 431 "ctrevc3.f"
		    i__4 = k + k * t_dim1;
#line 431 "ctrevc3.f"
		    t[i__4].r = smin, t[i__4].i = 0.;
#line 431 "ctrevc3.f"
		}
#line 433 "ctrevc3.f"
/* L50: */
#line 433 "ctrevc3.f"
	    }

#line 435 "ctrevc3.f"
	    if (ki > 1) {
#line 436 "ctrevc3.f"
		i__1 = ki - 1;
#line 436 "ctrevc3.f"
		clatrs_("Upper", "No transpose", "Non-unit", "Y", &i__1, &t[
			t_offset], ldt, &work[iv * *n + 1], &scale, &rwork[1],
			 info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 439 "ctrevc3.f"
		i__1 = ki + iv * *n;
#line 439 "ctrevc3.f"
		work[i__1].r = scale, work[i__1].i = 0.;
#line 440 "ctrevc3.f"
	    }

/*           Copy the vector x or Q*x to VR and normalize. */

#line 444 "ctrevc3.f"
	    if (! over) {
/*              ------------------------------ */
/*              no back-transform: copy x to VR and normalize. */
#line 447 "ctrevc3.f"
		ccopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 1], 
			&c__1);

#line 449 "ctrevc3.f"
		ii = icamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
#line 450 "ctrevc3.f"
		i__1 = ii + is * vr_dim1;
#line 450 "ctrevc3.f"
		remax = 1. / ((d__1 = vr[i__1].r, abs(d__1)) + (d__2 = d_imag(
			&vr[ii + is * vr_dim1]), abs(d__2)));
#line 451 "ctrevc3.f"
		csscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 453 "ctrevc3.f"
		i__1 = *n;
#line 453 "ctrevc3.f"
		for (k = ki + 1; k <= i__1; ++k) {
#line 454 "ctrevc3.f"
		    i__3 = k + is * vr_dim1;
#line 454 "ctrevc3.f"
		    vr[i__3].r = 0., vr[i__3].i = 0.;
#line 455 "ctrevc3.f"
/* L60: */
#line 455 "ctrevc3.f"
		}

#line 457 "ctrevc3.f"
	    } else if (nb == 1) {
/*              ------------------------------ */
/*              version 1: back-transform each vector with GEMV, Q*x. */
#line 460 "ctrevc3.f"
		if (ki > 1) {
#line 460 "ctrevc3.f"
		    i__1 = ki - 1;
#line 460 "ctrevc3.f"
		    z__1.r = scale, z__1.i = 0.;
#line 460 "ctrevc3.f"
		    cgemv_("N", n, &i__1, &c_b2, &vr[vr_offset], ldvr, &work[
			    iv * *n + 1], &c__1, &z__1, &vr[ki * vr_dim1 + 1],
			     &c__1, (ftnlen)1);
#line 460 "ctrevc3.f"
		}

#line 465 "ctrevc3.f"
		ii = icamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
#line 466 "ctrevc3.f"
		i__1 = ii + ki * vr_dim1;
#line 466 "ctrevc3.f"
		remax = 1. / ((d__1 = vr[i__1].r, abs(d__1)) + (d__2 = d_imag(
			&vr[ii + ki * vr_dim1]), abs(d__2)));
#line 467 "ctrevc3.f"
		csscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);

#line 469 "ctrevc3.f"
	    } else {
/*              ------------------------------ */
/*              version 2: back-transform block of vectors with GEMM */
/*              zero out below vector */
#line 473 "ctrevc3.f"
		i__1 = *n;
#line 473 "ctrevc3.f"
		for (k = ki + 1; k <= i__1; ++k) {
#line 474 "ctrevc3.f"
		    i__3 = k + iv * *n;
#line 474 "ctrevc3.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 475 "ctrevc3.f"
		}

/*              Columns IV:NB of work are valid vectors. */
/*              When the number of vectors stored reaches NB, */
/*              or if this was last vector, do the GEMM */
#line 480 "ctrevc3.f"
		if (iv == 1 || ki == 1) {
#line 481 "ctrevc3.f"
		    i__1 = nb - iv + 1;
#line 481 "ctrevc3.f"
		    i__3 = ki + nb - iv;
#line 481 "ctrevc3.f"
		    cgemm_("N", "N", n, &i__1, &i__3, &c_b2, &vr[vr_offset], 
			    ldvr, &work[iv * *n + 1], n, &c_b1, &work[(nb + 
			    iv) * *n + 1], n, (ftnlen)1, (ftnlen)1);
/*                 normalize vectors */
#line 487 "ctrevc3.f"
		    i__1 = nb;
#line 487 "ctrevc3.f"
		    for (k = iv; k <= i__1; ++k) {
#line 488 "ctrevc3.f"
			ii = icamax_(n, &work[(nb + k) * *n + 1], &c__1);
#line 489 "ctrevc3.f"
			i__3 = ii + (nb + k) * *n;
#line 489 "ctrevc3.f"
			remax = 1. / ((d__1 = work[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&work[ii + (nb + k) * *n]), abs(
				d__2)));
#line 490 "ctrevc3.f"
			csscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
#line 491 "ctrevc3.f"
		    }
#line 492 "ctrevc3.f"
		    i__1 = nb - iv + 1;
#line 492 "ctrevc3.f"
		    clacpy_("F", n, &i__1, &work[(nb + iv) * *n + 1], n, &vr[
			    ki * vr_dim1 + 1], ldvr, (ftnlen)1);
#line 495 "ctrevc3.f"
		    iv = nb;
#line 496 "ctrevc3.f"
		} else {
#line 497 "ctrevc3.f"
		    --iv;
#line 498 "ctrevc3.f"
		}
#line 499 "ctrevc3.f"
	    }

/*           Restore the original diagonal elements of T. */

#line 503 "ctrevc3.f"
	    i__1 = ki - 1;
#line 503 "ctrevc3.f"
	    for (k = 1; k <= i__1; ++k) {
#line 504 "ctrevc3.f"
		i__3 = k + k * t_dim1;
#line 504 "ctrevc3.f"
		i__4 = k;
#line 504 "ctrevc3.f"
		t[i__3].r = work[i__4].r, t[i__3].i = work[i__4].i;
#line 505 "ctrevc3.f"
/* L70: */
#line 505 "ctrevc3.f"
	    }

#line 507 "ctrevc3.f"
	    --is;
#line 508 "ctrevc3.f"
L80:
#line 508 "ctrevc3.f"
	    ;
#line 508 "ctrevc3.f"
	}
#line 509 "ctrevc3.f"
    }

#line 511 "ctrevc3.f"
    if (leftv) {

/*        ============================================================ */
/*        Compute left eigenvectors. */

/*        IV is index of column in current block. */
/*        Non-blocked version always uses IV=1; */
/*        blocked     version starts with IV=1, goes up to NB. */
/*        (Note the "0-th" column is used to store the original diagonal.) */
#line 520 "ctrevc3.f"
	iv = 1;
#line 521 "ctrevc3.f"
	is = 1;
#line 522 "ctrevc3.f"
	i__1 = *n;
#line 522 "ctrevc3.f"
	for (ki = 1; ki <= i__1; ++ki) {

#line 524 "ctrevc3.f"
	    if (somev) {
#line 525 "ctrevc3.f"
		if (! select[ki]) {
#line 525 "ctrevc3.f"
		    goto L130;
#line 525 "ctrevc3.f"
		}
#line 527 "ctrevc3.f"
	    }
/* Computing MAX */
#line 528 "ctrevc3.f"
	    i__3 = ki + ki * t_dim1;
#line 528 "ctrevc3.f"
	    d__3 = ulp * ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = d_imag(&t[
		    ki + ki * t_dim1]), abs(d__2)));
#line 528 "ctrevc3.f"
	    smin = max(d__3,smlnum);

/*           -------------------------------------------------------- */
/*           Complex left eigenvector */

#line 533 "ctrevc3.f"
	    i__3 = ki + iv * *n;
#line 533 "ctrevc3.f"
	    work[i__3].r = 1., work[i__3].i = 0.;

/*           Form right-hand side. */

#line 537 "ctrevc3.f"
	    i__3 = *n;
#line 537 "ctrevc3.f"
	    for (k = ki + 1; k <= i__3; ++k) {
#line 538 "ctrevc3.f"
		i__4 = k + iv * *n;
#line 538 "ctrevc3.f"
		d_cnjg(&z__2, &t[ki + k * t_dim1]);
#line 538 "ctrevc3.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 538 "ctrevc3.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 539 "ctrevc3.f"
/* L90: */
#line 539 "ctrevc3.f"
	    }

/*           Solve conjugate-transposed triangular system: */
/*           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK. */

#line 544 "ctrevc3.f"
	    i__3 = *n;
#line 544 "ctrevc3.f"
	    for (k = ki + 1; k <= i__3; ++k) {
#line 545 "ctrevc3.f"
		i__4 = k + k * t_dim1;
#line 545 "ctrevc3.f"
		i__5 = k + k * t_dim1;
#line 545 "ctrevc3.f"
		i__6 = ki + ki * t_dim1;
#line 545 "ctrevc3.f"
		z__1.r = t[i__5].r - t[i__6].r, z__1.i = t[i__5].i - t[i__6]
			.i;
#line 545 "ctrevc3.f"
		t[i__4].r = z__1.r, t[i__4].i = z__1.i;
#line 546 "ctrevc3.f"
		i__4 = k + k * t_dim1;
#line 546 "ctrevc3.f"
		if ((d__1 = t[i__4].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * 
			t_dim1]), abs(d__2)) < smin) {
#line 546 "ctrevc3.f"
		    i__5 = k + k * t_dim1;
#line 546 "ctrevc3.f"
		    t[i__5].r = smin, t[i__5].i = 0.;
#line 546 "ctrevc3.f"
		}
#line 548 "ctrevc3.f"
/* L100: */
#line 548 "ctrevc3.f"
	    }

#line 550 "ctrevc3.f"
	    if (ki < *n) {
#line 551 "ctrevc3.f"
		i__3 = *n - ki;
#line 551 "ctrevc3.f"
		clatrs_("Upper", "Conjugate transpose", "Non-unit", "Y", &
			i__3, &t[ki + 1 + (ki + 1) * t_dim1], ldt, &work[ki + 
			1 + iv * *n], &scale, &rwork[1], info, (ftnlen)5, (
			ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 554 "ctrevc3.f"
		i__3 = ki + iv * *n;
#line 554 "ctrevc3.f"
		work[i__3].r = scale, work[i__3].i = 0.;
#line 555 "ctrevc3.f"
	    }

/*           Copy the vector x or Q*x to VL and normalize. */

#line 559 "ctrevc3.f"
	    if (! over) {
/*              ------------------------------ */
/*              no back-transform: copy x to VL and normalize. */
#line 562 "ctrevc3.f"
		i__3 = *n - ki + 1;
#line 562 "ctrevc3.f"
		ccopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * 
			vl_dim1], &c__1);

#line 564 "ctrevc3.f"
		i__3 = *n - ki + 1;
#line 564 "ctrevc3.f"
		ii = icamax_(&i__3, &vl[ki + is * vl_dim1], &c__1) + ki - 1;
#line 565 "ctrevc3.f"
		i__3 = ii + is * vl_dim1;
#line 565 "ctrevc3.f"
		remax = 1. / ((d__1 = vl[i__3].r, abs(d__1)) + (d__2 = d_imag(
			&vl[ii + is * vl_dim1]), abs(d__2)));
#line 566 "ctrevc3.f"
		i__3 = *n - ki + 1;
#line 566 "ctrevc3.f"
		csscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);

#line 568 "ctrevc3.f"
		i__3 = ki - 1;
#line 568 "ctrevc3.f"
		for (k = 1; k <= i__3; ++k) {
#line 569 "ctrevc3.f"
		    i__4 = k + is * vl_dim1;
#line 569 "ctrevc3.f"
		    vl[i__4].r = 0., vl[i__4].i = 0.;
#line 570 "ctrevc3.f"
/* L110: */
#line 570 "ctrevc3.f"
		}

#line 572 "ctrevc3.f"
	    } else if (nb == 1) {
/*              ------------------------------ */
/*              version 1: back-transform each vector with GEMV, Q*x. */
#line 575 "ctrevc3.f"
		if (ki < *n) {
#line 575 "ctrevc3.f"
		    i__3 = *n - ki;
#line 575 "ctrevc3.f"
		    z__1.r = scale, z__1.i = 0.;
#line 575 "ctrevc3.f"
		    cgemv_("N", n, &i__3, &c_b2, &vl[(ki + 1) * vl_dim1 + 1], 
			    ldvl, &work[ki + 1 + iv * *n], &c__1, &z__1, &vl[
			    ki * vl_dim1 + 1], &c__1, (ftnlen)1);
#line 575 "ctrevc3.f"
		}

#line 580 "ctrevc3.f"
		ii = icamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
#line 581 "ctrevc3.f"
		i__3 = ii + ki * vl_dim1;
#line 581 "ctrevc3.f"
		remax = 1. / ((d__1 = vl[i__3].r, abs(d__1)) + (d__2 = d_imag(
			&vl[ii + ki * vl_dim1]), abs(d__2)));
#line 582 "ctrevc3.f"
		csscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);

#line 584 "ctrevc3.f"
	    } else {
/*              ------------------------------ */
/*              version 2: back-transform block of vectors with GEMM */
/*              zero out above vector */
/*              could go from KI-NV+1 to KI-1 */
#line 589 "ctrevc3.f"
		i__3 = ki - 1;
#line 589 "ctrevc3.f"
		for (k = 1; k <= i__3; ++k) {
#line 590 "ctrevc3.f"
		    i__4 = k + iv * *n;
#line 590 "ctrevc3.f"
		    work[i__4].r = 0., work[i__4].i = 0.;
#line 591 "ctrevc3.f"
		}

/*              Columns 1:IV of work are valid vectors. */
/*              When the number of vectors stored reaches NB, */
/*              or if this was last vector, do the GEMM */
#line 596 "ctrevc3.f"
		if (iv == nb || ki == *n) {
#line 597 "ctrevc3.f"
		    i__3 = *n - ki + iv;
#line 597 "ctrevc3.f"
		    cgemm_("N", "N", n, &iv, &i__3, &c_b2, &vl[(ki - iv + 1) *
			     vl_dim1 + 1], ldvl, &work[ki - iv + 1 + *n], n, &
			    c_b1, &work[(nb + 1) * *n + 1], n, (ftnlen)1, (
			    ftnlen)1);
/*                 normalize vectors */
#line 603 "ctrevc3.f"
		    i__3 = iv;
#line 603 "ctrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 604 "ctrevc3.f"
			ii = icamax_(n, &work[(nb + k) * *n + 1], &c__1);
#line 605 "ctrevc3.f"
			i__4 = ii + (nb + k) * *n;
#line 605 "ctrevc3.f"
			remax = 1. / ((d__1 = work[i__4].r, abs(d__1)) + (
				d__2 = d_imag(&work[ii + (nb + k) * *n]), abs(
				d__2)));
#line 606 "ctrevc3.f"
			csscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
#line 607 "ctrevc3.f"
		    }
#line 608 "ctrevc3.f"
		    clacpy_("F", n, &iv, &work[(nb + 1) * *n + 1], n, &vl[(ki 
			    - iv + 1) * vl_dim1 + 1], ldvl, (ftnlen)1);
#line 611 "ctrevc3.f"
		    iv = 1;
#line 612 "ctrevc3.f"
		} else {
#line 613 "ctrevc3.f"
		    ++iv;
#line 614 "ctrevc3.f"
		}
#line 615 "ctrevc3.f"
	    }

/*           Restore the original diagonal elements of T. */

#line 619 "ctrevc3.f"
	    i__3 = *n;
#line 619 "ctrevc3.f"
	    for (k = ki + 1; k <= i__3; ++k) {
#line 620 "ctrevc3.f"
		i__4 = k + k * t_dim1;
#line 620 "ctrevc3.f"
		i__5 = k;
#line 620 "ctrevc3.f"
		t[i__4].r = work[i__5].r, t[i__4].i = work[i__5].i;
#line 621 "ctrevc3.f"
/* L120: */
#line 621 "ctrevc3.f"
	    }

#line 623 "ctrevc3.f"
	    ++is;
#line 624 "ctrevc3.f"
L130:
#line 624 "ctrevc3.f"
	    ;
#line 624 "ctrevc3.f"
	}
#line 625 "ctrevc3.f"
    }

#line 627 "ctrevc3.f"
    return 0;

/*     End of CTREVC3 */

} /* ctrevc3_ */

