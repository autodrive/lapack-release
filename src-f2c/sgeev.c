#line 1 "sgeev.f"
/* sgeev.f -- translated by f2c (version 20100827).
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

#line 1 "sgeev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> SGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, */
/*                         LDVR, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WI( * ), WORK( * ), WR( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEEV computes for an N-by-N real nonsymmetric matrix A, the */
/* > eigenvalues and, optionally, the left and/or right eigenvectors. */
/* > */
/* > The right eigenvector v(j) of A satisfies */
/* >                  A * v(j) = lambda(j) * v(j) */
/* > where lambda(j) is its eigenvalue. */
/* > The left eigenvector u(j) of A satisfies */
/* >               u(j)**H * A = lambda(j) * u(j)**H */
/* > where u(j)**H denotes the conjugate-transpose of u(j). */
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
/* >          = 'V': left eigenvectors of A are computed. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (N) */
/* >          WR and WI contain the real and imaginary parts, */
/* >          respectively, of the computed eigenvalues.  Complex */
/* >          conjugate pairs of eigenvalues appear consecutively */
/* >          with the eigenvalue having the positive imaginary part */
/* >          first. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is REAL array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* >          after another in the columns of VL, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVL = 'N', VL is not referenced. */
/* >          If the j-th eigenvalue is real, then u(j) = VL(:,j), */
/* >          the j-th column of VL. */
/* >          If the j-th and (j+1)-st eigenvalues form a complex */
/* >          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and */
/* >          u(j+1) = VL(:,j) - i*VL(:,j+1). */
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
/* >          VR is REAL array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* >          after another in the columns of VR, in the same order */
/* >          as their eigenvalues. */
/* >          If JOBVR = 'N', VR is not referenced. */
/* >          If the j-th eigenvalue is real, then v(j) = VR(:,j), */
/* >          the j-th column of VR. */
/* >          If the j-th and (j+1)-st eigenvalues form a complex */
/* >          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and */
/* >          v(j+1) = VR(:,j) - i*VR(:,j+1). */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,3*N), and */
/* >          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good */
/* >          performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the QR algorithm failed to compute all the */
/* >                eigenvalues, and no eigenvectors have been computed; */
/* >                elements i+1:N of WR and WI contain eigenvalues which */
/* >                have converged. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGEeigen */

/*  ===================================================================== */
/* Subroutine */ int sgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, 
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, 
	integer *lwork, integer *info, ftnlen jobvl_len, ftnlen jobvr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal r__, cs, sn;
    static integer ihi;
    static doublereal scl;
    static integer ilo;
    static doublereal dum[1], eps;
    static integer ibal;
    static char side[1];
    static doublereal anrm;
    static integer ierr, itau, iwrk, nout;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal slapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    static logical scalea;
    static doublereal cscale;
    extern /* Subroutine */ int sgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sgebal_(char *, integer *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int sgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), sorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), shseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), strevc_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static integer minwrk, maxwrk;
    static logical wantvl;
    static doublereal smlnum;
    static integer hswork;
    static logical lquery, wantvr;


/*  -- LAPACK driver routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 243 "sgeev.f"
    /* Parameter adjustments */
#line 243 "sgeev.f"
    a_dim1 = *lda;
#line 243 "sgeev.f"
    a_offset = 1 + a_dim1;
#line 243 "sgeev.f"
    a -= a_offset;
#line 243 "sgeev.f"
    --wr;
#line 243 "sgeev.f"
    --wi;
#line 243 "sgeev.f"
    vl_dim1 = *ldvl;
#line 243 "sgeev.f"
    vl_offset = 1 + vl_dim1;
#line 243 "sgeev.f"
    vl -= vl_offset;
#line 243 "sgeev.f"
    vr_dim1 = *ldvr;
#line 243 "sgeev.f"
    vr_offset = 1 + vr_dim1;
#line 243 "sgeev.f"
    vr -= vr_offset;
#line 243 "sgeev.f"
    --work;
#line 243 "sgeev.f"

#line 243 "sgeev.f"
    /* Function Body */
#line 243 "sgeev.f"
    *info = 0;
#line 244 "sgeev.f"
    lquery = *lwork == -1;
#line 245 "sgeev.f"
    wantvl = lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1);
#line 246 "sgeev.f"
    wantvr = lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1);
#line 247 "sgeev.f"
    if (! wantvl && ! lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 248 "sgeev.f"
	*info = -1;
#line 249 "sgeev.f"
    } else if (! wantvr && ! lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 250 "sgeev.f"
	*info = -2;
#line 251 "sgeev.f"
    } else if (*n < 0) {
#line 252 "sgeev.f"
	*info = -3;
#line 253 "sgeev.f"
    } else if (*lda < max(1,*n)) {
#line 254 "sgeev.f"
	*info = -5;
#line 255 "sgeev.f"
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
#line 256 "sgeev.f"
	*info = -9;
#line 257 "sgeev.f"
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
#line 258 "sgeev.f"
	*info = -11;
#line 259 "sgeev.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by SHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 271 "sgeev.f"
    if (*info == 0) {
#line 272 "sgeev.f"
	if (*n == 0) {
#line 273 "sgeev.f"
	    minwrk = 1;
#line 274 "sgeev.f"
	    maxwrk = 1;
#line 275 "sgeev.f"
	} else {
#line 276 "sgeev.f"
	    maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "SGEHRD", " ", n, &c__1, 
		    n, &c__0, (ftnlen)6, (ftnlen)1);
#line 277 "sgeev.f"
	    if (wantvl) {
#line 278 "sgeev.f"
		minwrk = *n << 2;
/* Computing MAX */
#line 279 "sgeev.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"SORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
#line 279 "sgeev.f"
		maxwrk = max(i__1,i__2);
#line 281 "sgeev.f"
		shseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vl[vl_offset], ldvl, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 283 "sgeev.f"
		hswork = (integer) work[1];
/* Computing MAX */
#line 284 "sgeev.f"
		i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *
			n + hswork;
#line 284 "sgeev.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 285 "sgeev.f"
		i__1 = maxwrk, i__2 = *n << 2;
#line 285 "sgeev.f"
		maxwrk = max(i__1,i__2);
#line 286 "sgeev.f"
	    } else if (wantvr) {
#line 287 "sgeev.f"
		minwrk = *n << 2;
/* Computing MAX */
#line 288 "sgeev.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"SORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
#line 288 "sgeev.f"
		maxwrk = max(i__1,i__2);
#line 290 "sgeev.f"
		shseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 292 "sgeev.f"
		hswork = (integer) work[1];
/* Computing MAX */
#line 293 "sgeev.f"
		i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *
			n + hswork;
#line 293 "sgeev.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 294 "sgeev.f"
		i__1 = maxwrk, i__2 = *n << 2;
#line 294 "sgeev.f"
		maxwrk = max(i__1,i__2);
#line 295 "sgeev.f"
	    } else {
#line 296 "sgeev.f"
		minwrk = *n * 3;
#line 297 "sgeev.f"
		shseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[
			1], &vr[vr_offset], ldvr, &work[1], &c_n1, info, (
			ftnlen)1, (ftnlen)1);
#line 299 "sgeev.f"
		hswork = (integer) work[1];
/* Computing MAX */
#line 300 "sgeev.f"
		i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *
			n + hswork;
#line 300 "sgeev.f"
		maxwrk = max(i__1,i__2);
#line 301 "sgeev.f"
	    }
#line 302 "sgeev.f"
	    maxwrk = max(maxwrk,minwrk);
#line 303 "sgeev.f"
	}
#line 304 "sgeev.f"
	work[1] = (doublereal) maxwrk;

#line 306 "sgeev.f"
	if (*lwork < minwrk && ! lquery) {
#line 307 "sgeev.f"
	    *info = -13;
#line 308 "sgeev.f"
	}
#line 309 "sgeev.f"
    }

#line 311 "sgeev.f"
    if (*info != 0) {
#line 312 "sgeev.f"
	i__1 = -(*info);
#line 312 "sgeev.f"
	xerbla_("SGEEV ", &i__1, (ftnlen)6);
#line 313 "sgeev.f"
	return 0;
#line 314 "sgeev.f"
    } else if (lquery) {
#line 315 "sgeev.f"
	return 0;
#line 316 "sgeev.f"
    }

/*     Quick return if possible */

#line 320 "sgeev.f"
    if (*n == 0) {
#line 320 "sgeev.f"
	return 0;
#line 320 "sgeev.f"
    }

/*     Get machine constants */

#line 325 "sgeev.f"
    eps = slamch_("P", (ftnlen)1);
#line 326 "sgeev.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 327 "sgeev.f"
    bignum = 1. / smlnum;
#line 328 "sgeev.f"
    slabad_(&smlnum, &bignum);
#line 329 "sgeev.f"
    smlnum = sqrt(smlnum) / eps;
#line 330 "sgeev.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 334 "sgeev.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 335 "sgeev.f"
    scalea = FALSE_;
#line 336 "sgeev.f"
    if (anrm > 0. && anrm < smlnum) {
#line 337 "sgeev.f"
	scalea = TRUE_;
#line 338 "sgeev.f"
	cscale = smlnum;
#line 339 "sgeev.f"
    } else if (anrm > bignum) {
#line 340 "sgeev.f"
	scalea = TRUE_;
#line 341 "sgeev.f"
	cscale = bignum;
#line 342 "sgeev.f"
    }
#line 343 "sgeev.f"
    if (scalea) {
#line 343 "sgeev.f"
	slascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 343 "sgeev.f"
    }

/*     Balance the matrix */
/*     (Workspace: need N) */

#line 349 "sgeev.f"
    ibal = 1;
#line 350 "sgeev.f"
    sgebal_("B", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (Workspace: need 3*N, prefer 2*N+N*NB) */

#line 355 "sgeev.f"
    itau = ibal + *n;
#line 356 "sgeev.f"
    iwrk = itau + *n;
#line 357 "sgeev.f"
    i__1 = *lwork - iwrk + 1;
#line 357 "sgeev.f"
    sgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 360 "sgeev.f"
    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

#line 365 "sgeev.f"
	*(unsigned char *)side = 'L';
#line 366 "sgeev.f"
	slacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VL */
/*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

#line 371 "sgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 371 "sgeev.f"
	sorghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

#line 377 "sgeev.f"
	iwrk = itau;
#line 378 "sgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 378 "sgeev.f"
	shseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &
		vl[vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 381 "sgeev.f"
	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

#line 386 "sgeev.f"
	    *(unsigned char *)side = 'B';
#line 387 "sgeev.f"
	    slacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (
		    ftnlen)1);
#line 388 "sgeev.f"
	}

#line 390 "sgeev.f"
    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

#line 395 "sgeev.f"
	*(unsigned char *)side = 'R';
#line 396 "sgeev.f"
	slacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VR */
/*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

#line 401 "sgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 401 "sgeev.f"
	sorghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk],
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

#line 407 "sgeev.f"
	iwrk = itau;
#line 408 "sgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 408 "sgeev.f"
	shseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &
		vr[vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);

#line 411 "sgeev.f"
    } else {

/*        Compute eigenvalues only */
/*        (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

#line 416 "sgeev.f"
	iwrk = itau;
#line 417 "sgeev.f"
	i__1 = *lwork - iwrk + 1;
#line 417 "sgeev.f"
	shseqr_("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &
		vr[vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen)1, (
		ftnlen)1);
#line 419 "sgeev.f"
    }

/*     If INFO > 0 from SHSEQR, then quit */

#line 423 "sgeev.f"
    if (*info > 0) {
#line 423 "sgeev.f"
	goto L50;
#line 423 "sgeev.f"
    }

#line 426 "sgeev.f"
    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (Workspace: need 4*N) */

#line 431 "sgeev.f"
	strevc_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
		 &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &ierr, (ftnlen)
		1, (ftnlen)1);
#line 433 "sgeev.f"
    }

#line 435 "sgeev.f"
    if (wantvl) {

/*        Undo balancing of left eigenvectors */
/*        (Workspace: need N) */

#line 440 "sgeev.f"
	sgebak_("B", "L", n, &ilo, &ihi, &work[ibal], n, &vl[vl_offset], ldvl,
		 &ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize left eigenvectors and make largest component real */

#line 445 "sgeev.f"
	i__1 = *n;
#line 445 "sgeev.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 446 "sgeev.f"
	    if (wi[i__] == 0.) {
#line 447 "sgeev.f"
		scl = 1. / snrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 448 "sgeev.f"
		sscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 449 "sgeev.f"
	    } else if (wi[i__] > 0.) {
#line 450 "sgeev.f"
		d__1 = snrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
#line 450 "sgeev.f"
		d__2 = snrm2_(n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
#line 450 "sgeev.f"
		scl = 1. / slapy2_(&d__1, &d__2);
#line 452 "sgeev.f"
		sscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
#line 453 "sgeev.f"
		sscal_(n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
#line 454 "sgeev.f"
		i__2 = *n;
#line 454 "sgeev.f"
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
#line 455 "sgeev.f"
		    d__1 = vl[k + i__ * vl_dim1];
/* Computing 2nd power */
#line 455 "sgeev.f"
		    d__2 = vl[k + (i__ + 1) * vl_dim1];
#line 455 "sgeev.f"
		    work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
#line 456 "sgeev.f"
/* L10: */
#line 456 "sgeev.f"
		}
#line 457 "sgeev.f"
		k = isamax_(n, &work[iwrk], &c__1);
#line 458 "sgeev.f"
		slartg_(&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], 
			&cs, &sn, &r__);
#line 459 "sgeev.f"
		srot_(n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * 
			vl_dim1 + 1], &c__1, &cs, &sn);
#line 460 "sgeev.f"
		vl[k + (i__ + 1) * vl_dim1] = 0.;
#line 461 "sgeev.f"
	    }
#line 462 "sgeev.f"
/* L20: */
#line 462 "sgeev.f"
	}
#line 463 "sgeev.f"
    }

#line 465 "sgeev.f"
    if (wantvr) {

/*        Undo balancing of right eigenvectors */
/*        (Workspace: need N) */

#line 470 "sgeev.f"
	sgebak_("B", "R", n, &ilo, &ihi, &work[ibal], n, &vr[vr_offset], ldvr,
		 &ierr, (ftnlen)1, (ftnlen)1);

/*        Normalize right eigenvectors and make largest component real */

#line 475 "sgeev.f"
	i__1 = *n;
#line 475 "sgeev.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 476 "sgeev.f"
	    if (wi[i__] == 0.) {
#line 477 "sgeev.f"
		scl = 1. / snrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 478 "sgeev.f"
		sscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 479 "sgeev.f"
	    } else if (wi[i__] > 0.) {
#line 480 "sgeev.f"
		d__1 = snrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
#line 480 "sgeev.f"
		d__2 = snrm2_(n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
#line 480 "sgeev.f"
		scl = 1. / slapy2_(&d__1, &d__2);
#line 482 "sgeev.f"
		sscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
#line 483 "sgeev.f"
		sscal_(n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
#line 484 "sgeev.f"
		i__2 = *n;
#line 484 "sgeev.f"
		for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
#line 485 "sgeev.f"
		    d__1 = vr[k + i__ * vr_dim1];
/* Computing 2nd power */
#line 485 "sgeev.f"
		    d__2 = vr[k + (i__ + 1) * vr_dim1];
#line 485 "sgeev.f"
		    work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
#line 486 "sgeev.f"
/* L30: */
#line 486 "sgeev.f"
		}
#line 487 "sgeev.f"
		k = isamax_(n, &work[iwrk], &c__1);
#line 488 "sgeev.f"
		slartg_(&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], 
			&cs, &sn, &r__);
#line 489 "sgeev.f"
		srot_(n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * 
			vr_dim1 + 1], &c__1, &cs, &sn);
#line 490 "sgeev.f"
		vr[k + (i__ + 1) * vr_dim1] = 0.;
#line 491 "sgeev.f"
	    }
#line 492 "sgeev.f"
/* L40: */
#line 492 "sgeev.f"
	}
#line 493 "sgeev.f"
    }

/*     Undo scaling if necessary */

#line 497 "sgeev.f"
L50:
#line 498 "sgeev.f"
    if (scalea) {
#line 499 "sgeev.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 499 "sgeev.f"
	i__3 = *n - *info;
#line 499 "sgeev.f"
	i__2 = max(i__3,1);
#line 499 "sgeev.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 
		1], &i__2, &ierr, (ftnlen)1);
#line 501 "sgeev.f"
	i__1 = *n - *info;
/* Computing MAX */
#line 501 "sgeev.f"
	i__3 = *n - *info;
#line 501 "sgeev.f"
	i__2 = max(i__3,1);
#line 501 "sgeev.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 
		1], &i__2, &ierr, (ftnlen)1);
#line 503 "sgeev.f"
	if (*info > 0) {
#line 504 "sgeev.f"
	    i__1 = ilo - 1;
#line 504 "sgeev.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], 
		    n, &ierr, (ftnlen)1);
#line 506 "sgeev.f"
	    i__1 = ilo - 1;
#line 506 "sgeev.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], 
		    n, &ierr, (ftnlen)1);
#line 508 "sgeev.f"
	}
#line 509 "sgeev.f"
    }

#line 511 "sgeev.f"
    work[1] = (doublereal) maxwrk;
#line 512 "sgeev.f"
    return 0;

/*     End of SGEEV */

} /* sgeev_ */

