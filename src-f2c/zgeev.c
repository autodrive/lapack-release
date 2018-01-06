#line 1 "zgeev.f"
/* zgeev.f -- translated by f2c (version 20100827).
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

#line 1 "zgeev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, */
/*                         WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the */
/* > eigenvalues and, optionally, the left and/or right eigenvectors. */
/* > */
/* > The right eigenvector v(j) of A satisfies */
/* >                  A * v(j) = lambda(j) * v(j) */
/* > where lambda(j) is its eigenvalue. */
/* > The left eigenvector u(j) of A satisfies */
/* >               u(j)**H * A = lambda(j) * u(j)**H */
/* > where u(j)**H denotes the conjugate transpose of u(j). */
/* > */
/* > The computed eigenvectors are normalized to have Euclidean norm */
/* > equal to 1 and largest component real. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBVL */
/* > \verbatim */
/* >          JOBVL is CHARACTER*1 */
/* >          = 'N': left eigenvectors of A are not computed; */
/* >          = 'V': left eigenvectors of are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* >          JOBVR is CHARACTER*1 */
/* >          = 'N': right eigenvectors of A are not computed; */
/* >          = 'V': right eigenvectors of A are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A has been overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          W contains the computed eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* >          after another in the columns of VL, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVL = 'N', VL is not referenced. */
/* >          u(j) = VL(:,j), the j-th column of VL. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL.  LDVL >= 1; if */
/* >          JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* >          VR is COMPLEX*16 array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* >          after another in the columns of VR, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVR = 'N', VR is not referenced. */
/* >          v(j) = VR(:,j), the j-th column of VR. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR.  LDVR >= 1; if */
/* >          JOBVR = 'V', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the QR algorithm failed to compute all the */
/* >                eigenvalues, and no eigenvectors have been computed; */
/* >                elements and i+1:N of W contain eigenvalues which have */
/* >                converged. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/*  @precisions fortran z -> c */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgeev_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info, ftnlen jobvl_len, 
	ftnlen jobvr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k, ihi;
    static doublereal scl;
    static integer ilo;
    static doublereal dum[1], eps;
    static doublecomplex tmp;
    static integer lwork_trevc__, ibal;
    static char side[1];
    static doublereal anrm;
    static integer ierr, itau, iwrk, nout;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    static logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int zgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen, ftnlen), zgebal_(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int zgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), zlacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);
    static integer minwrk, maxwrk;
    static logical wantvl;
    static doublereal smlnum;
    static integer hswork, irwork;
    extern /* Subroutine */ int zhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunghr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static logical lquery, wantvr;
    extern /* Subroutine */ int ztrevc3_(char *, char *, logical *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 233 "zgeev.f"
    /* Parameter adjustments */
#line 233 "zgeev.f"
    a_dim1 = *lda;
#line 233 "zgeev.f"
    a_offset = 1 + a_dim1;
#line 233 "zgeev.f"
    a -= a_offset;
#line 233 "zgeev.f"
    --w;
#line 233 "zgeev.f"
    vl_dim1 = *ldvl;
#line 233 "zgeev.f"
    vl_offset = 1 + vl_dim1;
#line 233 "zgeev.f"
    vl -= vl_offset;
#line 233 "zgeev.f"
    vr_dim1 = *ldvr;
#line 233 "zgeev.f"
    vr_offset = 1 + vr_dim1;
#line 233 "zgeev.f"
    vr -= vr_offset;
#line 233 "zgeev.f"
    --work;
#line 233 "zgeev.f"
    --rwork;
#line 233 "zgeev.f"

#line 233 "zgeev.f"
    /* Function Body */
#line 233 "zgeev.f"
    *info = 0;
#line 234 "zgeev.f"
    lquery = *lwork == -1;
#line 235 "zgeev.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 236 "zgeev.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 237 "zgeev.f"
    if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "zgeev.f"
	*info = -1;
#line 239 "zgeev.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 240 "zgeev.f"
	*info = -2;
#line 241 "zgeev.f"
    } else if (*n < 0) {
#line 242 "zgeev.f"
	*info = -3;
#line 243 "zgeev.f"
    } else if (*lda < max(1,*n)) {
#line 244 "zgeev.f"
	*info = -5;
#line 245 "zgeev.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 246 "zgeev.f"
	*info = -8;
#line 247 "zgeev.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 248 "zgeev.f"
	*info = -10;
#line 249 "zgeev.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by ZHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 262 "zgeev.f"
    if (*info == 0) {
#line 263 "zgeev.f"
	if (*n == 0) {
#line 264 "zgeev.f"
	    minwrk = 1;
#line 265 "zgeev.f"
	    maxwrk = 1;
#line 266 "zgeev.f"
	} else {
#line 267 "zgeev.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);
#line 268 "zgeev.f"
	    minwrk = *n << 1;
#line 269 "zgeev.f"
	    if (wantvl) {
/* Computing MAX */
#line 270 "zgeev.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 270 "zgeev.f"
		maxwrk = max(i__1,i__2);
#line 272 "zgeev.f"
		ztrevc3_("L", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &rwork[1], &c_n1, &ierr, (ftnlen)1, (
			ftnlen)1);
#line 275 "zgeev.f"
		lwork_trevc__ = (integer) work[1].r;
/* Computing MAX */
#line 276 "zgeev.f"
		i__1 = maxwrk, i__2 = *n + lwork_trevc__;
#line 276 "zgeev.f"
		maxwrk = max(i__1,i__2);
#line 277 "zgeev.f"
		zhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[
			vl_offset], ldvl, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 279 "zgeev.f"
	    } else if (wantvr) {
/* Computing MAX */
#line 280 "zgeev.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 280 "zgeev.f"
		maxwrk = max(i__1,i__2);
#line 282 "zgeev.f"
		ztrevc3_("R", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &rwork[1], &c_n1, &ierr, (ftnlen)1, (
			ftnlen)1);
#line 285 "zgeev.f"
		lwork_trevc__ = (integer) work[1].r;
/* Computing MAX */
#line 286 "zgeev.f"
		i__1 = maxwrk, i__2 = *n + lwork_trevc__;
#line 286 "zgeev.f"
		maxwrk = max(i__1,i__2);
#line 287 "zgeev.f"
		zhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 289 "zgeev.f"
	    } else {
#line 290 "zgeev.f"
		zhseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 292 "zgeev.f"
	    }
#line 293 "zgeev.f"
	    hswork = (integer) work[1].r;
/* Computing MAX */
#line 294 "zgeev.f"
	    i__1 = max(maxwrk,hswork);
#line 294 "zgeev.f"
	    maxwrk = max(i__1,minwrk);
#line 295 "zgeev.f"
	}
#line 296 "zgeev.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 298 "zgeev.f"
	if (*lwork < minwrk && ! lquery) {
#line 299 "zgeev.f"
	    *info = -12;
#line 300 "zgeev.f"
	}
#line 301 "zgeev.f"
    }

#line 303 "zgeev.f"
    if (*info != 0) {
#line 304 "zgeev.f"
	i__1 = -(*info);
#line 304 "zgeev.f"
	xerbla_("ZGEEV ", &i__1, (ftnlen)6);
#line 305 "zgeev.f"
	return 0;
#line 306 "zgeev.f"
    } else if (lquery) {
#line 307 "zgeev.f"
	return 0;
#line 308 "zgeev.f"
    }

/*     Quick return if possible */

#line 312 "zgeev.f"
    if (*n == 0) {
#line 312 "zgeev.f"
	return 0;
#line 312 "zgeev.f"
    }

/*     Get machine constants */

#line 317 "zgeev.f"
    eps = dlamch_("P", (ftnlen)1);
#line 318 "zgeev.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 319 "zgeev.f"
    bignum = 1. / smlnum;
#line 320 "zgeev.f"
    dlabad_(&smlnum, &bignum);
#line 321 "zgeev.f"
    smlnum = sqrt(smlnum) / eps;
#line 322 "zgeev.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 326 "zgeev.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 327 "zgeev.f"
    scalea = FALSE_;
#line 328 "zgeev.f"
    if (anrm > 0. && anrm < smlnum) {
#line 329 "zgeev.f"
	scalea = TRUE_;
#line 330 "zgeev.f"
	cscale = smlnum;
#line 331 "zgeev.f"
    } else if (anrm > bignum) {
#line 332 "zgeev.f"
	scalea = TRUE_;
#line 333 "zgeev.f"
	cscale = bignum;
#line 334 "zgeev.f"
    }
#line 335 "zgeev.f"
    if (scalea) {
#line 335 "zgeev.f"
	zlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 335 "zgeev.f"
    }

/*     Balance the matrix */
/*     (CWorkspace: none) */
/*     (RWorkspace: need N) */

#line 342 "zgeev.f"
    ibal = 1;
#line 343 "zgeev.f"
    zgebal_("B", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 349 "zgeev.f"
    itau = 1;
#line 350 "zgeev.f"
    iwrk = itau + *n;
#line 351 "zgeev.f"
    i__1 = *lwork - iwrk + 1;
#line 351 "zgeev.f"
    zgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 354 "zgeev.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 359 "zgeev.f"
	*(unsigned char *)side = 'L';
#line 360 "zgeev.f"
	zlacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate unitary matrix in VL */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 366 "zgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 366 "zgeev.f"
	zunghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 373 "zgeev.f"
	iwrk = itau;
#line 374 "zgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 374 "zgeev.f"
	zhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 377 "zgeev.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 382 "zgeev.f"
	    *(unsigned char *)side = 'B';
#line 383 "zgeev.f"
	    zlacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 384 "zgeev.f"
	}

#line 386 "zgeev.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 391 "zgeev.f"
	*(unsigned char *)side = 'R';
#line 392 "zgeev.f"
	zlacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate unitary matrix in VR */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 398 "zgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 398 "zgeev.f"
	zunghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 405 "zgeev.f"
	iwrk = itau;
#line 406 "zgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 406 "zgeev.f"
	zhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 409 "zgeev.f"
    } else {

/*        Compute eigenvalues only */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 415 "zgeev.f"
	iwrk = itau;
#line 416 "zgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 416 "zgeev.f"
	zhseqr_("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 418 "zgeev.f"
    }

/*     If INFO .NE. 0 from ZHSEQR, then quit */

#line 422 "zgeev.f"
    if (*info != 0) {
#line 422 "zgeev.f"
	goto L50;
#line 422 "zgeev.f"
    }

#line 425 "zgeev.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (CWorkspace: need 2*N, prefer N + 2*N*NB) */
/*        (RWorkspace: need 2*N) */

#line 431 "zgeev.f"
	irwork = ibal + *n;
#line 432 "zgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 432 "zgeev.f"
	ztrevc3_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &i__1, &
		rwork[irwork], n, &ierr, (ftnlen)1, (ftnlen)1);
#line 435 "zgeev.f"
    }

#line 437 "zgeev.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 443 "zgeev.f"
	zgebak_("B", "L", n, &ilo, &ihi, &rwork[ibal], n, &vl[vl_offset], 
		ldvl, &ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 448 "zgeev.f"
	i__1 = *n;
#line 448 "zgeev.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 449 "zgeev.f"
	    scl = 1. / dznrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 450 "zgeev.f"
	    zdscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 451 "zgeev.f"
	    i__2 = *n;
#line 451 "zgeev.f"
	    for (k = 1; k <= i__2; ++k) {
#line 452 "zgeev.f"
		i__3 = k + i__ * vl_dim1;
/* Computing 2nd power */
#line 452 "zgeev.f"
		d__1 = vl[i__3].r;
/* Computing 2nd power */
#line 452 "zgeev.f"
		d__2 = d_imag(&vl[k + i__ * vl_dim1]);
#line 452 "zgeev.f"
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
#line 454 "zgeev.f"
/* L10: */
#line 454 "zgeev.f"
	    }
#line 455 "zgeev.f"
	    k = idamax_(n, &rwork[irwork], &c__1);
#line 456 "zgeev.f"
	    d_cnjg(&z__2, &vl[k + i__ * vl_dim1]);
#line 456 "zgeev.f"
	    d__1 = sqrt(rwork[irwork + k - 1]);
#line 456 "zgeev.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 456 "zgeev.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 457 "zgeev.f"
	    zscal_(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
#line 458 "zgeev.f"
	    i__2 = k + i__ * vl_dim1;
#line 458 "zgeev.f"
	    i__3 = k + i__ * vl_dim1;
#line 458 "zgeev.f"
	    d__1 = vl[i__3].r;
#line 458 "zgeev.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 458 "zgeev.f"
	    vl[i__2].r = z__1.r, vl[i__2].i = z__1.i;
#line 459 "zgeev.f"
/* L20: */
#line 459 "zgeev.f"
	}
#line 460 "zgeev.f"
    }

#line 462 "zgeev.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 468 "zgeev.f"
	zgebak_("B", "R", n, &ilo, &ihi, &rwork[ibal], n, &vr[vr_offset], 
		ldvr, &ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 473 "zgeev.f"
	i__1 = *n;
#line 473 "zgeev.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 474 "zgeev.f"
	    scl = 1. / dznrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 475 "zgeev.f"
	    zdscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 476 "zgeev.f"
	    i__2 = *n;
#line 476 "zgeev.f"
	    for (k = 1; k <= i__2; ++k) {
#line 477 "zgeev.f"
		i__3 = k + i__ * vr_dim1;
/* Computing 2nd power */
#line 477 "zgeev.f"
		d__1 = vr[i__3].r;
/* Computing 2nd power */
#line 477 "zgeev.f"
		d__2 = d_imag(&vr[k + i__ * vr_dim1]);
#line 477 "zgeev.f"
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
#line 479 "zgeev.f"
/* L30: */
#line 479 "zgeev.f"
	    }
#line 480 "zgeev.f"
	    k = idamax_(n, &rwork[irwork], &c__1);
#line 481 "zgeev.f"
	    d_cnjg(&z__2, &vr[k + i__ * vr_dim1]);
#line 481 "zgeev.f"
	    d__1 = sqrt(rwork[irwork + k - 1]);
#line 481 "zgeev.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 481 "zgeev.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 482 "zgeev.f"
	    zscal_(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
#line 483 "zgeev.f"
	    i__2 = k + i__ * vr_dim1;
#line 483 "zgeev.f"
	    i__3 = k + i__ * vr_dim1;
#line 483 "zgeev.f"
	    d__1 = vr[i__3].r;
#line 483 "zgeev.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 483 "zgeev.f"
	    vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
#line 484 "zgeev.f"
/* L40: */
#line 484 "zgeev.f"
	}
#line 485 "zgeev.f"
    }

/*     Undo scaling if necessary */

#line 489 "zgeev.f"
L50:
#line 490 "zgeev.f"
    if (scalea) {
#line 491 "zgeev.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 491 "zgeev.f"
	i__3 = *n - *info;
#line 491 "zgeev.f"
	i__2 = max(i__3,1);
#line 491 "zgeev.f"
	zlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1]
		, &i__2, &ierr, (ftnlen)1);
#line 493 "zgeev.f"
	if (*info > 0) {
#line 494 "zgeev.f"
	    i__1 = ilo - 1;
#line 494 "zgeev.f"
	    zlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n,
		     &ierr, (ftnlen)1);
#line 495 "zgeev.f"
	}
#line 496 "zgeev.f"
    }

#line 498 "zgeev.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 499 "zgeev.f"
    return 0;

/*     End of ZGEEV */

} /* zgeev_ */

