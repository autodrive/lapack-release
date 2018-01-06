#line 1 "cgeev.f"
/* cgeev.f -- translated by f2c (version 20100827).
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

#line 1 "cgeev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> CGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, */
/*                         WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL   RWORK( * ) */
/*       COMPLEX         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEEV computes for an N-by-N complex nonsymmetric matrix A, the */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          W is COMPLEX array, dimension (N) */
/* >          W contains the computed eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is COMPLEX array, dimension (LDVL,N) */
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
/* >          VR is COMPLEX array, dimension (LDVR,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is REAL array, dimension (2*N) */
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

/*  @generated from zgeev.f, fortran z -> c, Tue Apr 19 01:47:44 2016 */

/* > \ingroup complexGEeigen */

/*  ===================================================================== */
/* Subroutine */ int cgeev_(char *jobvl, char *jobvr, integer *n, 
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen, ftnlen), cgebal_(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen), slabad_(doublereal *, doublereal *);
    static logical scalea;
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int cgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), clacpy_(char *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int chseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), cunghr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer minwrk, maxwrk;
    static logical wantvl;
    static doublereal smlnum;
    static integer hswork, irwork;
    static logical lquery, wantvr;
    extern /* Subroutine */ int ctrevc3_(char *, char *, logical *, integer *,
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

#line 233 "cgeev.f"
    /* Parameter adjustments */
#line 233 "cgeev.f"
    a_dim1 = *lda;
#line 233 "cgeev.f"
    a_offset = 1 + a_dim1;
#line 233 "cgeev.f"
    a -= a_offset;
#line 233 "cgeev.f"
    --w;
#line 233 "cgeev.f"
    vl_dim1 = *ldvl;
#line 233 "cgeev.f"
    vl_offset = 1 + vl_dim1;
#line 233 "cgeev.f"
    vl -= vl_offset;
#line 233 "cgeev.f"
    vr_dim1 = *ldvr;
#line 233 "cgeev.f"
    vr_offset = 1 + vr_dim1;
#line 233 "cgeev.f"
    vr -= vr_offset;
#line 233 "cgeev.f"
    --work;
#line 233 "cgeev.f"
    --rwork;
#line 233 "cgeev.f"

#line 233 "cgeev.f"
    /* Function Body */
#line 233 "cgeev.f"
    *info = 0;
#line 234 "cgeev.f"
    lquery = *lwork == -1;
#line 235 "cgeev.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 236 "cgeev.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 237 "cgeev.f"
    if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "cgeev.f"
	*info = -1;
#line 239 "cgeev.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 240 "cgeev.f"
	*info = -2;
#line 241 "cgeev.f"
    } else if (*n < 0) {
#line 242 "cgeev.f"
	*info = -3;
#line 243 "cgeev.f"
    } else if (*lda < max(1,*n)) {
#line 244 "cgeev.f"
	*info = -5;
#line 245 "cgeev.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 246 "cgeev.f"
	*info = -8;
#line 247 "cgeev.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 248 "cgeev.f"
	*info = -10;
#line 249 "cgeev.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by CHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 262 "cgeev.f"
    if (*info == 0) {
#line 263 "cgeev.f"
	if (*n == 0) {
#line 264 "cgeev.f"
	    minwrk = 1;
#line 265 "cgeev.f"
	    maxwrk = 1;
#line 266 "cgeev.f"
	} else {
#line 267 "cgeev.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "CGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);
#line 268 "cgeev.f"
	    minwrk = *n << 1;
#line 269 "cgeev.f"
	    if (wantvl) {
/* Computing MAX */
#line 270 "cgeev.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 270 "cgeev.f"
		maxwrk = max(i__1,i__2);
#line 272 "cgeev.f"
		ctrevc3_("L", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &rwork[1], &c_n1, &ierr, (ftnlen)1, (
			ftnlen)1);
#line 275 "cgeev.f"
		lwork_trevc__ = (integer) work[1].r;
/* Computing MAX */
#line 276 "cgeev.f"
		i__1 = maxwrk, i__2 = *n + lwork_trevc__;
#line 276 "cgeev.f"
		maxwrk = max(i__1,i__2);
#line 277 "cgeev.f"
		chseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[
			vl_offset], ldvl, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 279 "cgeev.f"
	    } else if (wantvr) {
/* Computing MAX */
#line 280 "cgeev.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 280 "cgeev.f"
		maxwrk = max(i__1,i__2);
#line 282 "cgeev.f"
		ctrevc3_("R", "B", select, n, &a[a_offset], lda, &vl[
			vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &
			work[1], &c_n1, &rwork[1], &c_n1, &ierr, (ftnlen)1, (
			ftnlen)1);
#line 285 "cgeev.f"
		lwork_trevc__ = (integer) work[1].r;
/* Computing MAX */
#line 286 "cgeev.f"
		i__1 = maxwrk, i__2 = *n + lwork_trevc__;
#line 286 "cgeev.f"
		maxwrk = max(i__1,i__2);
#line 287 "cgeev.f"
		chseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 289 "cgeev.f"
	    } else {
#line 290 "cgeev.f"
		chseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info, (ftnlen)1, (
			ftnlen)1);
#line 292 "cgeev.f"
	    }
#line 293 "cgeev.f"
	    hswork = (integer) work[1].r;
/* Computing MAX */
#line 294 "cgeev.f"
	    i__1 = max(maxwrk,hswork);
#line 294 "cgeev.f"
	    maxwrk = max(i__1,minwrk);
#line 295 "cgeev.f"
	}
#line 296 "cgeev.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 298 "cgeev.f"
	if (*lwork < minwrk && ! lquery) {
#line 299 "cgeev.f"
	    *info = -12;
#line 300 "cgeev.f"
	}
#line 301 "cgeev.f"
    }

#line 303 "cgeev.f"
    if (*info != 0) {
#line 304 "cgeev.f"
	i__1 = -(*info);
#line 304 "cgeev.f"
	xerbla_("CGEEV ", &i__1, (ftnlen)6);
#line 305 "cgeev.f"
	return 0;
#line 306 "cgeev.f"
    } else if (lquery) {
#line 307 "cgeev.f"
	return 0;
#line 308 "cgeev.f"
    }

/*     Quick return if possible */

#line 312 "cgeev.f"
    if (*n == 0) {
#line 312 "cgeev.f"
	return 0;
#line 312 "cgeev.f"
    }

/*     Get machine constants */

#line 317 "cgeev.f"
    eps = slamch_("P", (ftnlen)1);
#line 318 "cgeev.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 319 "cgeev.f"
    bignum = 1. / smlnum;
#line 320 "cgeev.f"
    slabad_(&smlnum, &bignum);
#line 321 "cgeev.f"
    smlnum = sqrt(smlnum) / eps;
#line 322 "cgeev.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 326 "cgeev.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 327 "cgeev.f"
    scalea = FALSE_;
#line 328 "cgeev.f"
    if (anrm > 0. && anrm < smlnum) {
#line 329 "cgeev.f"
	scalea = TRUE_;
#line 330 "cgeev.f"
	cscale = smlnum;
#line 331 "cgeev.f"
    } else if (anrm > bignum) {
#line 332 "cgeev.f"
	scalea = TRUE_;
#line 333 "cgeev.f"
	cscale = bignum;
#line 334 "cgeev.f"
    }
#line 335 "cgeev.f"
    if (scalea) {
#line 335 "cgeev.f"
	clascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 335 "cgeev.f"
    }

/*     Balance the matrix */
/*     (CWorkspace: none) */
/*     (RWorkspace: need N) */

#line 342 "cgeev.f"
    ibal = 1;
#line 343 "cgeev.f"
    cgebal_("B", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 349 "cgeev.f"
    itau = 1;
#line 350 "cgeev.f"
    iwrk = itau + *n;
#line 351 "cgeev.f"
    i__1 = *lwork - iwrk + 1;
#line 351 "cgeev.f"
    cgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 354 "cgeev.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 359 "cgeev.f"
	*(unsigned char *)side = 'L';
#line 360 "cgeev.f"
	clacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate unitary matrix in VL */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 366 "cgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 366 "cgeev.f"
	cunghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 373 "cgeev.f"
	iwrk = itau;
#line 374 "cgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 374 "cgeev.f"
	chseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 377 "cgeev.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 382 "cgeev.f"
	    *(unsigned char *)side = 'B';
#line 383 "cgeev.f"
	    clacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 384 "cgeev.f"
	}

#line 386 "cgeev.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 391 "cgeev.f"
	*(unsigned char *)side = 'R';
#line 392 "cgeev.f"
	clacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate unitary matrix in VR */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 398 "cgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 398 "cgeev.f"
	cunghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 405 "cgeev.f"
	iwrk = itau;
#line 406 "cgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 406 "cgeev.f"
	chseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 409 "cgeev.f"
    } else {

/*        Compute eigenvalues only */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

#line 415 "cgeev.f"
	iwrk = itau;
#line 416 "cgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 416 "cgeev.f"
	chseqr_("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 418 "cgeev.f"
    }

/*     If INFO .NE. 0 from CHSEQR, then quit */

#line 422 "cgeev.f"
    if (*info != 0) {
#line 422 "cgeev.f"
	goto L50;
#line 422 "cgeev.f"
    }

#line 425 "cgeev.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (CWorkspace: need 2*N, prefer N + 2*N*NB) */
/*        (RWorkspace: need 2*N) */

#line 431 "cgeev.f"
	irwork = ibal + *n;
#line 432 "cgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 432 "cgeev.f"
	ctrevc3_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], 
		ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &i__1, &
		rwork[irwork], n, &ierr, (ftnlen)1, (ftnlen)1);
#line 435 "cgeev.f"
    }

#line 437 "cgeev.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 443 "cgeev.f"
	cgebak_("B", "L", n, &ilo, &ihi, &rwork[ibal], n, &vl[vl_offset], 
		ldvl, &ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 448 "cgeev.f"
	i__1 = *n;
#line 448 "cgeev.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 449 "cgeev.f"
	    scl = 1. / scnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 450 "cgeev.f"
	    csscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 451 "cgeev.f"
	    i__2 = *n;
#line 451 "cgeev.f"
	    for (k = 1; k <= i__2; ++k) {
#line 452 "cgeev.f"
		i__3 = k + i__ * vl_dim1;
/* Computing 2nd power */
#line 452 "cgeev.f"
		d__1 = vl[i__3].r;
/* Computing 2nd power */
#line 452 "cgeev.f"
		d__2 = d_imag(&vl[k + i__ * vl_dim1]);
#line 452 "cgeev.f"
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
#line 454 "cgeev.f"
/* L10: */
#line 454 "cgeev.f"
	    }
#line 455 "cgeev.f"
	    k = isamax_(n, &rwork[irwork], &c__1);
#line 456 "cgeev.f"
	    d_cnjg(&z__2, &vl[k + i__ * vl_dim1]);
#line 456 "cgeev.f"
	    d__1 = sqrt(rwork[irwork + k - 1]);
#line 456 "cgeev.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 456 "cgeev.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 457 "cgeev.f"
	    cscal_(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
#line 458 "cgeev.f"
	    i__2 = k + i__ * vl_dim1;
#line 458 "cgeev.f"
	    i__3 = k + i__ * vl_dim1;
#line 458 "cgeev.f"
	    d__1 = vl[i__3].r;
#line 458 "cgeev.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 458 "cgeev.f"
	    vl[i__2].r = z__1.r, vl[i__2].i = z__1.i;
#line 459 "cgeev.f"
/* L20: */
#line 459 "cgeev.f"
	}
#line 460 "cgeev.f"
    }

#line 462 "cgeev.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 468 "cgeev.f"
	cgebak_("B", "R", n, &ilo, &ihi, &rwork[ibal], n, &vr[vr_offset], 
		ldvr, &ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 473 "cgeev.f"
	i__1 = *n;
#line 473 "cgeev.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 474 "cgeev.f"
	    scl = 1. / scnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 475 "cgeev.f"
	    csscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 476 "cgeev.f"
	    i__2 = *n;
#line 476 "cgeev.f"
	    for (k = 1; k <= i__2; ++k) {
#line 477 "cgeev.f"
		i__3 = k + i__ * vr_dim1;
/* Computing 2nd power */
#line 477 "cgeev.f"
		d__1 = vr[i__3].r;
/* Computing 2nd power */
#line 477 "cgeev.f"
		d__2 = d_imag(&vr[k + i__ * vr_dim1]);
#line 477 "cgeev.f"
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
#line 479 "cgeev.f"
/* L30: */
#line 479 "cgeev.f"
	    }
#line 480 "cgeev.f"
	    k = isamax_(n, &rwork[irwork], &c__1);
#line 481 "cgeev.f"
	    d_cnjg(&z__2, &vr[k + i__ * vr_dim1]);
#line 481 "cgeev.f"
	    d__1 = sqrt(rwork[irwork + k - 1]);
#line 481 "cgeev.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 481 "cgeev.f"
	    tmp.r = z__1.r, tmp.i = z__1.i;
#line 482 "cgeev.f"
	    cscal_(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
#line 483 "cgeev.f"
	    i__2 = k + i__ * vr_dim1;
#line 483 "cgeev.f"
	    i__3 = k + i__ * vr_dim1;
#line 483 "cgeev.f"
	    d__1 = vr[i__3].r;
#line 483 "cgeev.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 483 "cgeev.f"
	    vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
#line 484 "cgeev.f"
/* L40: */
#line 484 "cgeev.f"
	}
#line 485 "cgeev.f"
    }

/*     Undo scaling if necessary */

#line 489 "cgeev.f"
L50:
#line 490 "cgeev.f"
    if (scalea) {
#line 491 "cgeev.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 491 "cgeev.f"
	i__3 = *n - *info;
#line 491 "cgeev.f"
	i__2 = max(i__3,1);
#line 491 "cgeev.f"
	clascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1]
		, &i__2, &ierr, (ftnlen)1);
#line 493 "cgeev.f"
	if (*info > 0) {
#line 494 "cgeev.f"
	    i__1 = ilo - 1;
#line 494 "cgeev.f"
	    clascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n,
		     &ierr, (ftnlen)1);
#line 495 "cgeev.f"
	}
#line 496 "cgeev.f"
    }

#line 498 "cgeev.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 499 "cgeev.f"
    return 0;

/*     End of CGEEV */

} /* cgeev_ */

