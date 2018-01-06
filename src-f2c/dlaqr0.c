#line 1 "dlaqr0.f"
/* dlaqr0.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqr0.f"
/* Table of constant values */

static integer c__13 = 13;
static integer c__15 = 15;
static integer c_n1 = -1;
static integer c__12 = 12;
static integer c__14 = 14;
static integer c__16 = 16;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__3 = 3;

/* > \brief \b DLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Sc
hur decomposition. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQR0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, */
/*                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLAQR0 computes the eigenvalues of a Hessenberg matrix H */
/* >    and, optionally, the matrices T and Z from the Schur decomposition */
/* >    H = Z T Z**T, where T is an upper quasi-triangular matrix (the */
/* >    Schur form), and Z is the orthogonal matrix of Schur vectors. */
/* > */
/* >    Optionally Z may be postmultiplied into an input orthogonal */
/* >    matrix Q so that this routine can give the Schur factorization */
/* >    of a matrix A which has been reduced to the Hessenberg form H */
/* >    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T. */
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
/* >           The order of the matrix H.  N .GE. 0. */
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
/* >           It is assumed that H is already upper triangular in rows */
/* >           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1, */
/* >           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a */
/* >           previous call to DGEBAL, and then passed to DGEHRD when the */
/* >           matrix output by DGEBAL is reduced to Hessenberg form. */
/* >           Otherwise, ILO and IHI should be set to 1 and N, */
/* >           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N. */
/* >           If N = 0, then ILO = 1 and IHI = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION array, dimension (LDH,N) */
/* >           On entry, the upper Hessenberg matrix H. */
/* >           On exit, if INFO = 0 and WANTT is .TRUE., then H contains */
/* >           the upper quasi-triangular matrix T from the Schur */
/* >           decomposition (the Schur form); 2-by-2 diagonal blocks */
/* >           (corresponding to complex conjugate pairs of eigenvalues) */
/* >           are returned in standard form, with H(i,i) = H(i+1,i+1) */
/* >           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is */
/* >           .FALSE., then the contents of H are unspecified on exit. */
/* >           (The output value of H when INFO.GT.0 is given under the */
/* >           description of INFO below.) */
/* > */
/* >           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and */
/* >           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >           The leading dimension of the array H. LDH .GE. max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is DOUBLE PRECISION array, dimension (IHI) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION array, dimension (IHI) */
/* >           The real and imaginary parts, respectively, of the computed */
/* >           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI) */
/* >           and WI(ILO:IHI). If two eigenvalues are computed as a */
/* >           complex conjugate pair, they are stored in consecutive */
/* >           elements of WR and WI, say the i-th and (i+1)th, with */
/* >           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then */
/* >           the eigenvalues are stored in the same order as on the */
/* >           diagonal of the Schur form returned in H, with */
/* >           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal */
/* >           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and */
/* >           WI(i+1) = -WI(i). */
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
/* >           Specify the rows of Z to which transformations must be */
/* >           applied if WANTZ is .TRUE.. */
/* >           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ,IHI) */
/* >           If WANTZ is .FALSE., then Z is not referenced. */
/* >           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is */
/* >           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the */
/* >           orthogonal Schur factor of H(ILO:IHI,ILO:IHI). */
/* >           (The output value of Z when INFO.GT.0 is given under */
/* >           the description of INFO below.) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >           The leading dimension of the array Z.  if WANTZ is .TRUE. */
/* >           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension LWORK */
/* >           On exit, if LWORK = -1, WORK(1) returns an estimate of */
/* >           the optimal value for LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >           The dimension of the array WORK.  LWORK .GE. max(1,N) */
/* >           is sufficient, but LWORK typically as large as 6*N may */
/* >           be required for optimal performance.  A workspace query */
/* >           to determine the optimal workspace size is recommended. */
/* > */
/* >           If LWORK = -1, then DLAQR0 does a workspace query. */
/* >           In this case, DLAQR0 checks the input parameters and */
/* >           estimates the optimal workspace size for the given */
/* >           values of N, ILO and IHI.  The estimate is returned */
/* >           in WORK(1).  No error message related to LWORK is */
/* >           issued by XERBLA.  Neither H nor Z are accessed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >             =  0:  successful exit */
/* >           .GT. 0:  if INFO = i, DLAQR0 failed to compute all of */
/* >                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR */
/* >                and WI contain those eigenvalues which have been */
/* >                successfully computed.  (Failures are rare.) */
/* > */
/* >                If INFO .GT. 0 and WANT is .FALSE., then on exit, */
/* >                the remaining unconverged eigenvalues are the eigen- */
/* >                values of the upper Hessenberg matrix rows and */
/* >                columns ILO through INFO of the final, output */
/* >                value of H. */
/* > */
/* >                If INFO .GT. 0 and WANTT is .TRUE., then on exit */
/* > */
/* >           (*)  (initial value of H)*U  = U*(final value of H) */
/* > */
/* >                where U is an orthogonal matrix.  The final */
/* >                value of H is upper Hessenberg and quasi-triangular */
/* >                in rows and columns INFO+1 through IHI. */
/* > */
/* >                If INFO .GT. 0 and WANTZ is .TRUE., then on exit */
/* > */
/* >                  (final value of Z(ILO:IHI,ILOZ:IHIZ) */
/* >                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U */
/* > */
/* >                where U is the orthogonal matrix in (*) (regard- */
/* >                less of the value of WANTT.) */
/* > */
/* >                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not */
/* >                accessed. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */

/* > \par References: */
/*  ================ */
/* > */
/* >       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/* >       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3 */
/* >       Performance, SIAM Journal of Matrix Analysis, volume 23, pages */
/* >       929--947, 2002. */
/* > \n */
/* >       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/* >       Algorithm Part II: Aggressive Early Deflation, SIAM Journal */
/* >       of Matrix Analysis, volume 23, pages 948--973, 2002. */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaqr0_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	*wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, k;
    static doublereal aa, bb, cc, dd;
    static integer ld;
    static doublereal cs;
    static integer nh, it, ks, kt;
    static doublereal sn;
    static integer ku, kv, ls, ns;
    static doublereal ss;
    static integer nw, inf, kdu, nho, nve, kwh, nsr, nwr, kwv, ndec, ndfl, 
	    kbot, nmin;
    static doublereal swap;
    static integer ktop;
    static doublereal zdum[1]	/* was [1][1] */;
    static integer kacc22, itmax, nsmax, nwmax, kwtop;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlaqr3_(
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    dlaqr4_(logical *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), dlaqr5_(logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    static integer nibble;
    extern /* Subroutine */ int dlahqr_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char jbcmpz[2];
    static integer nwupbd;
    static logical sorted;
    static integer lwkopt;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ================================================================ */

/*     .. Parameters .. */

/*     ==== Matrices of order NTINY or smaller must be processed by */
/*     .    DLAHQR because of insufficient subdiagonal scratch space. */
/*     .    (This is a hard limit.) ==== */

/*     ==== Exceptional deflation windows:  try to cure rare */
/*     .    slow convergence by varying the size of the */
/*     .    deflation window after KEXNW iterations. ==== */

/*     ==== Exceptional shifts: try to cure rare slow convergence */
/*     .    with ad-hoc exceptional shifts every KEXSH iterations. */
/*     .    ==== */

/*     ==== The constants WILK1 and WILK2 are used to form the */
/*     .    exceptional shifts. ==== */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 325 "dlaqr0.f"
    /* Parameter adjustments */
#line 325 "dlaqr0.f"
    h_dim1 = *ldh;
#line 325 "dlaqr0.f"
    h_offset = 1 + h_dim1;
#line 325 "dlaqr0.f"
    h__ -= h_offset;
#line 325 "dlaqr0.f"
    --wr;
#line 325 "dlaqr0.f"
    --wi;
#line 325 "dlaqr0.f"
    z_dim1 = *ldz;
#line 325 "dlaqr0.f"
    z_offset = 1 + z_dim1;
#line 325 "dlaqr0.f"
    z__ -= z_offset;
#line 325 "dlaqr0.f"
    --work;
#line 325 "dlaqr0.f"

#line 325 "dlaqr0.f"
    /* Function Body */
#line 325 "dlaqr0.f"
    *info = 0;

/*     ==== Quick return for N = 0: nothing to do. ==== */

#line 329 "dlaqr0.f"
    if (*n == 0) {
#line 330 "dlaqr0.f"
	work[1] = 1.;
#line 331 "dlaqr0.f"
	return 0;
#line 332 "dlaqr0.f"
    }

#line 334 "dlaqr0.f"
    if (*n <= 11) {

/*        ==== Tiny matrices must use DLAHQR. ==== */

#line 338 "dlaqr0.f"
	lwkopt = 1;
#line 339 "dlaqr0.f"
	if (*lwork != -1) {
#line 339 "dlaqr0.f"
	    dlahqr_(wantt, wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &
		    wi[1], iloz, ihiz, &z__[z_offset], ldz, info);
#line 339 "dlaqr0.f"
	}
#line 342 "dlaqr0.f"
    } else {

/*        ==== Use small bulge multi-shift QR with aggressive early */
/*        .    deflation on larger-than-tiny matrices. ==== */

/*        ==== Hope for the best. ==== */

#line 349 "dlaqr0.f"
	*info = 0;

/*        ==== Set up job flags for ILAENV. ==== */

#line 353 "dlaqr0.f"
	if (*wantt) {
#line 354 "dlaqr0.f"
	    *(unsigned char *)jbcmpz = 'S';
#line 355 "dlaqr0.f"
	} else {
#line 356 "dlaqr0.f"
	    *(unsigned char *)jbcmpz = 'E';
#line 357 "dlaqr0.f"
	}
#line 358 "dlaqr0.f"
	if (*wantz) {
#line 359 "dlaqr0.f"
	    *(unsigned char *)&jbcmpz[1] = 'V';
#line 360 "dlaqr0.f"
	} else {
#line 361 "dlaqr0.f"
	    *(unsigned char *)&jbcmpz[1] = 'N';
#line 362 "dlaqr0.f"
	}

/*        ==== NWR = recommended deflation window size.  At this */
/*        .    point,  N .GT. NTINY = 11, so there is enough */
/*        .    subdiagonal workspace for NWR.GE.2 as required. */
/*        .    (In fact, there is enough subdiagonal space for */
/*        .    NWR.GE.3.) ==== */

#line 370 "dlaqr0.f"
	nwr = ilaenv_(&c__13, "DLAQR0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
#line 371 "dlaqr0.f"
	nwr = max(2,nwr);
/* Computing MIN */
#line 372 "dlaqr0.f"
	i__1 = *ihi - *ilo + 1, i__2 = (*n - 1) / 3, i__1 = min(i__1,i__2);
#line 372 "dlaqr0.f"
	nwr = min(i__1,nwr);

/*        ==== NSR = recommended number of simultaneous shifts. */
/*        .    At this point N .GT. NTINY = 11, so there is at */
/*        .    enough subdiagonal workspace for NSR to be even */
/*        .    and greater than or equal to two as required. ==== */

#line 379 "dlaqr0.f"
	nsr = ilaenv_(&c__15, "DLAQR0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
/* Computing MIN */
#line 380 "dlaqr0.f"
	i__1 = nsr, i__2 = (*n + 6) / 9, i__1 = min(i__1,i__2), i__2 = *ihi - 
		*ilo;
#line 380 "dlaqr0.f"
	nsr = min(i__1,i__2);
/* Computing MAX */
#line 381 "dlaqr0.f"
	i__1 = 2, i__2 = nsr - nsr % 2;
#line 381 "dlaqr0.f"
	nsr = max(i__1,i__2);

/*        ==== Estimate optimal workspace ==== */

/*        ==== Workspace query call to DLAQR3 ==== */

#line 387 "dlaqr0.f"
	i__1 = nwr + 1;
#line 387 "dlaqr0.f"
	dlaqr3_(wantt, wantz, n, ilo, ihi, &i__1, &h__[h_offset], ldh, iloz, 
		ihiz, &z__[z_offset], ldz, &ls, &ld, &wr[1], &wi[1], &h__[
		h_offset], ldh, n, &h__[h_offset], ldh, n, &h__[h_offset], 
		ldh, &work[1], &c_n1);

/*        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ==== */

/* Computing MAX */
#line 393 "dlaqr0.f"
	i__1 = nsr * 3 / 2, i__2 = (integer) work[1];
#line 393 "dlaqr0.f"
	lwkopt = max(i__1,i__2);

/*        ==== Quick return in case of workspace query. ==== */

#line 397 "dlaqr0.f"
	if (*lwork == -1) {
#line 398 "dlaqr0.f"
	    work[1] = (doublereal) lwkopt;
#line 399 "dlaqr0.f"
	    return 0;
#line 400 "dlaqr0.f"
	}

/*        ==== DLAHQR/DLAQR0 crossover point ==== */

#line 404 "dlaqr0.f"
	nmin = ilaenv_(&c__12, "DLAQR0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)
		6, (ftnlen)2);
#line 405 "dlaqr0.f"
	nmin = max(11,nmin);

/*        ==== Nibble crossover point ==== */

#line 409 "dlaqr0.f"
	nibble = ilaenv_(&c__14, "DLAQR0", jbcmpz, n, ilo, ihi, lwork, (
		ftnlen)6, (ftnlen)2);
#line 410 "dlaqr0.f"
	nibble = max(0,nibble);

/*        ==== Accumulate reflections during ttswp?  Use block */
/*        .    2-by-2 structure during matrix-matrix multiply? ==== */

#line 415 "dlaqr0.f"
	kacc22 = ilaenv_(&c__16, "DLAQR0", jbcmpz, n, ilo, ihi, lwork, (
		ftnlen)6, (ftnlen)2);
#line 416 "dlaqr0.f"
	kacc22 = max(0,kacc22);
#line 417 "dlaqr0.f"
	kacc22 = min(2,kacc22);

/*        ==== NWMAX = the largest possible deflation window for */
/*        .    which there is sufficient workspace. ==== */

/* Computing MIN */
#line 422 "dlaqr0.f"
	i__1 = (*n - 1) / 3, i__2 = *lwork / 2;
#line 422 "dlaqr0.f"
	nwmax = min(i__1,i__2);
#line 423 "dlaqr0.f"
	nw = nwmax;

/*        ==== NSMAX = the Largest number of simultaneous shifts */
/*        .    for which there is sufficient workspace. ==== */

/* Computing MIN */
#line 428 "dlaqr0.f"
	i__1 = (*n + 6) / 9, i__2 = (*lwork << 1) / 3;
#line 428 "dlaqr0.f"
	nsmax = min(i__1,i__2);
#line 429 "dlaqr0.f"
	nsmax -= nsmax % 2;

/*        ==== NDFL: an iteration count restarted at deflation. ==== */

#line 433 "dlaqr0.f"
	ndfl = 1;

/*        ==== ITMAX = iteration limit ==== */

/* Computing MAX */
#line 437 "dlaqr0.f"
	i__1 = 10, i__2 = *ihi - *ilo + 1;
#line 437 "dlaqr0.f"
	itmax = max(i__1,i__2) * 30;

/*        ==== Last row and column in the active block ==== */

#line 441 "dlaqr0.f"
	kbot = *ihi;

/*        ==== Main Loop ==== */

#line 445 "dlaqr0.f"
	i__1 = itmax;
#line 445 "dlaqr0.f"
	for (it = 1; it <= i__1; ++it) {

/*           ==== Done when KBOT falls below ILO ==== */

#line 449 "dlaqr0.f"
	    if (kbot < *ilo) {
#line 449 "dlaqr0.f"
		goto L90;
#line 449 "dlaqr0.f"
	    }

/*           ==== Locate active block ==== */

#line 454 "dlaqr0.f"
	    i__2 = *ilo + 1;
#line 454 "dlaqr0.f"
	    for (k = kbot; k >= i__2; --k) {
#line 455 "dlaqr0.f"
		if (h__[k + (k - 1) * h_dim1] == 0.) {
#line 455 "dlaqr0.f"
		    goto L20;
#line 455 "dlaqr0.f"
		}
#line 457 "dlaqr0.f"
/* L10: */
#line 457 "dlaqr0.f"
	    }
#line 458 "dlaqr0.f"
	    k = *ilo;
#line 459 "dlaqr0.f"
L20:
#line 460 "dlaqr0.f"
	    ktop = k;

/*           ==== Select deflation window size: */
/*           .    Typical Case: */
/*           .      If possible and advisable, nibble the entire */
/*           .      active block.  If not, use size MIN(NWR,NWMAX) */
/*           .      or MIN(NWR+1,NWMAX) depending upon which has */
/*           .      the smaller corresponding subdiagonal entry */
/*           .      (a heuristic). */
/*           . */
/*           .    Exceptional Case: */
/*           .      If there have been no deflations in KEXNW or */
/*           .      more iterations, then vary the deflation window */
/*           .      size.   At first, because, larger windows are, */
/*           .      in general, more powerful than smaller ones, */
/*           .      rapidly increase the window to the maximum possible. */
/*           .      Then, gradually reduce the window size. ==== */

#line 478 "dlaqr0.f"
	    nh = kbot - ktop + 1;
#line 479 "dlaqr0.f"
	    nwupbd = min(nh,nwmax);
#line 480 "dlaqr0.f"
	    if (ndfl < 5) {
#line 481 "dlaqr0.f"
		nw = min(nwupbd,nwr);
#line 482 "dlaqr0.f"
	    } else {
/* Computing MIN */
#line 483 "dlaqr0.f"
		i__2 = nwupbd, i__3 = nw << 1;
#line 483 "dlaqr0.f"
		nw = min(i__2,i__3);
#line 484 "dlaqr0.f"
	    }
#line 485 "dlaqr0.f"
	    if (nw < nwmax) {
#line 486 "dlaqr0.f"
		if (nw >= nh - 1) {
#line 487 "dlaqr0.f"
		    nw = nh;
#line 488 "dlaqr0.f"
		} else {
#line 489 "dlaqr0.f"
		    kwtop = kbot - nw + 1;
#line 490 "dlaqr0.f"
		    if ((d__1 = h__[kwtop + (kwtop - 1) * h_dim1], abs(d__1)) 
			    > (d__2 = h__[kwtop - 1 + (kwtop - 2) * h_dim1], 
			    abs(d__2))) {
#line 490 "dlaqr0.f"
			++nw;
#line 490 "dlaqr0.f"
		    }
#line 492 "dlaqr0.f"
		}
#line 493 "dlaqr0.f"
	    }
#line 494 "dlaqr0.f"
	    if (ndfl < 5) {
#line 495 "dlaqr0.f"
		ndec = -1;
#line 496 "dlaqr0.f"
	    } else if (ndec >= 0 || nw >= nwupbd) {
#line 497 "dlaqr0.f"
		++ndec;
#line 498 "dlaqr0.f"
		if (nw - ndec < 2) {
#line 498 "dlaqr0.f"
		    ndec = 0;
#line 498 "dlaqr0.f"
		}
#line 500 "dlaqr0.f"
		nw -= ndec;
#line 501 "dlaqr0.f"
	    }

/*           ==== Aggressive early deflation: */
/*           .    split workspace under the subdiagonal into */
/*           .      - an nw-by-nw work array V in the lower */
/*           .        left-hand-corner, */
/*           .      - an NW-by-at-least-NW-but-more-is-better */
/*           .        (NW-by-NHO) horizontal work array along */
/*           .        the bottom edge, */
/*           .      - an at-least-NW-but-more-is-better (NHV-by-NW) */
/*           .        vertical work array along the left-hand-edge. */
/*           .        ==== */

#line 514 "dlaqr0.f"
	    kv = *n - nw + 1;
#line 515 "dlaqr0.f"
	    kt = nw + 1;
#line 516 "dlaqr0.f"
	    nho = *n - nw - 1 - kt + 1;
#line 517 "dlaqr0.f"
	    kwv = nw + 2;
#line 518 "dlaqr0.f"
	    nve = *n - nw - kwv + 1;

/*           ==== Aggressive early deflation ==== */

#line 522 "dlaqr0.f"
	    dlaqr3_(wantt, wantz, n, &ktop, &kbot, &nw, &h__[h_offset], ldh, 
		    iloz, ihiz, &z__[z_offset], ldz, &ls, &ld, &wr[1], &wi[1],
		     &h__[kv + h_dim1], ldh, &nho, &h__[kv + kt * h_dim1], 
		    ldh, &nve, &h__[kwv + h_dim1], ldh, &work[1], lwork);

/*           ==== Adjust KBOT accounting for new deflations. ==== */

#line 529 "dlaqr0.f"
	    kbot -= ld;

/*           ==== KS points to the shifts. ==== */

#line 533 "dlaqr0.f"
	    ks = kbot - ls + 1;

/*           ==== Skip an expensive QR sweep if there is a (partly */
/*           .    heuristic) reason to expect that many eigenvalues */
/*           .    will deflate without it.  Here, the QR sweep is */
/*           .    skipped if many eigenvalues have just been deflated */
/*           .    or if the remaining active block is small. */

#line 541 "dlaqr0.f"
	    if (ld == 0 || ld * 100 <= nw * nibble && kbot - ktop + 1 > min(
		    nmin,nwmax)) {

/*              ==== NS = nominal number of simultaneous shifts. */
/*              .    This may be lowered (slightly) if DLAQR3 */
/*              .    did not provide that many shifts. ==== */

/* Computing MIN */
/* Computing MAX */
#line 548 "dlaqr0.f"
		i__4 = 2, i__5 = kbot - ktop;
#line 548 "dlaqr0.f"
		i__2 = min(nsmax,nsr), i__3 = max(i__4,i__5);
#line 548 "dlaqr0.f"
		ns = min(i__2,i__3);
#line 549 "dlaqr0.f"
		ns -= ns % 2;

/*              ==== If there have been no deflations */
/*              .    in a multiple of KEXSH iterations, */
/*              .    then try exceptional shifts. */
/*              .    Otherwise use shifts provided by */
/*              .    DLAQR3 above or from the eigenvalues */
/*              .    of a trailing principal submatrix. ==== */

#line 558 "dlaqr0.f"
		if (ndfl % 6 == 0) {
#line 559 "dlaqr0.f"
		    ks = kbot - ns + 1;
/* Computing MAX */
#line 560 "dlaqr0.f"
		    i__3 = ks + 1, i__4 = ktop + 2;
#line 560 "dlaqr0.f"
		    i__2 = max(i__3,i__4);
#line 560 "dlaqr0.f"
		    for (i__ = kbot; i__ >= i__2; i__ += -2) {
#line 561 "dlaqr0.f"
			ss = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1))
				 + (d__2 = h__[i__ - 1 + (i__ - 2) * h_dim1], 
				abs(d__2));
#line 562 "dlaqr0.f"
			aa = ss * .75 + h__[i__ + i__ * h_dim1];
#line 563 "dlaqr0.f"
			bb = ss;
#line 564 "dlaqr0.f"
			cc = ss * -.4375;
#line 565 "dlaqr0.f"
			dd = aa;
#line 566 "dlaqr0.f"
			dlanv2_(&aa, &bb, &cc, &dd, &wr[i__ - 1], &wi[i__ - 1]
				, &wr[i__], &wi[i__], &cs, &sn);
#line 568 "dlaqr0.f"
/* L30: */
#line 568 "dlaqr0.f"
		    }
#line 569 "dlaqr0.f"
		    if (ks == ktop) {
#line 570 "dlaqr0.f"
			wr[ks + 1] = h__[ks + 1 + (ks + 1) * h_dim1];
#line 571 "dlaqr0.f"
			wi[ks + 1] = 0.;
#line 572 "dlaqr0.f"
			wr[ks] = wr[ks + 1];
#line 573 "dlaqr0.f"
			wi[ks] = wi[ks + 1];
#line 574 "dlaqr0.f"
		    }
#line 575 "dlaqr0.f"
		} else {

/*                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or */
/*                 .    DLAHQR on a trailing principal submatrix to */
/*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9, */
/*                 .    there is enough space below the subdiagonal */
/*                 .    to fit an NS-by-NS scratch array.) ==== */

#line 583 "dlaqr0.f"
		    if (kbot - ks + 1 <= ns / 2) {
#line 584 "dlaqr0.f"
			ks = kbot - ns + 1;
#line 585 "dlaqr0.f"
			kt = *n - ns + 1;
#line 586 "dlaqr0.f"
			dlacpy_("A", &ns, &ns, &h__[ks + ks * h_dim1], ldh, &
				h__[kt + h_dim1], ldh, (ftnlen)1);
#line 588 "dlaqr0.f"
			if (ns > nmin) {
#line 589 "dlaqr0.f"
			    dlaqr4_(&c_false, &c_false, &ns, &c__1, &ns, &h__[
				    kt + h_dim1], ldh, &wr[ks], &wi[ks], &
				    c__1, &c__1, zdum, &c__1, &work[1], lwork,
				     &inf);
#line 593 "dlaqr0.f"
			} else {
#line 594 "dlaqr0.f"
			    dlahqr_(&c_false, &c_false, &ns, &c__1, &ns, &h__[
				    kt + h_dim1], ldh, &wr[ks], &wi[ks], &
				    c__1, &c__1, zdum, &c__1, &inf);
#line 597 "dlaqr0.f"
			}
#line 598 "dlaqr0.f"
			ks += inf;

/*                    ==== In case of a rare QR failure use */
/*                    .    eigenvalues of the trailing 2-by-2 */
/*                    .    principal submatrix.  ==== */

#line 604 "dlaqr0.f"
			if (ks >= kbot) {
#line 605 "dlaqr0.f"
			    aa = h__[kbot - 1 + (kbot - 1) * h_dim1];
#line 606 "dlaqr0.f"
			    cc = h__[kbot + (kbot - 1) * h_dim1];
#line 607 "dlaqr0.f"
			    bb = h__[kbot - 1 + kbot * h_dim1];
#line 608 "dlaqr0.f"
			    dd = h__[kbot + kbot * h_dim1];
#line 609 "dlaqr0.f"
			    dlanv2_(&aa, &bb, &cc, &dd, &wr[kbot - 1], &wi[
				    kbot - 1], &wr[kbot], &wi[kbot], &cs, &sn)
				    ;
#line 612 "dlaqr0.f"
			    ks = kbot - 1;
#line 613 "dlaqr0.f"
			}
#line 614 "dlaqr0.f"
		    }

#line 616 "dlaqr0.f"
		    if (kbot - ks + 1 > ns) {

/*                    ==== Sort the shifts (Helps a little) */
/*                    .    Bubble sort keeps complex conjugate */
/*                    .    pairs together. ==== */

#line 622 "dlaqr0.f"
			sorted = FALSE_;
#line 623 "dlaqr0.f"
			i__2 = ks + 1;
#line 623 "dlaqr0.f"
			for (k = kbot; k >= i__2; --k) {
#line 624 "dlaqr0.f"
			    if (sorted) {
#line 624 "dlaqr0.f"
				goto L60;
#line 624 "dlaqr0.f"
			    }
#line 626 "dlaqr0.f"
			    sorted = TRUE_;
#line 627 "dlaqr0.f"
			    i__3 = k - 1;
#line 627 "dlaqr0.f"
			    for (i__ = ks; i__ <= i__3; ++i__) {
#line 628 "dlaqr0.f"
				if ((d__1 = wr[i__], abs(d__1)) + (d__2 = wi[
					i__], abs(d__2)) < (d__3 = wr[i__ + 1]
					, abs(d__3)) + (d__4 = wi[i__ + 1], 
					abs(d__4))) {
#line 630 "dlaqr0.f"
				    sorted = FALSE_;

#line 632 "dlaqr0.f"
				    swap = wr[i__];
#line 633 "dlaqr0.f"
				    wr[i__] = wr[i__ + 1];
#line 634 "dlaqr0.f"
				    wr[i__ + 1] = swap;

#line 636 "dlaqr0.f"
				    swap = wi[i__];
#line 637 "dlaqr0.f"
				    wi[i__] = wi[i__ + 1];
#line 638 "dlaqr0.f"
				    wi[i__ + 1] = swap;
#line 639 "dlaqr0.f"
				}
#line 640 "dlaqr0.f"
/* L40: */
#line 640 "dlaqr0.f"
			    }
#line 641 "dlaqr0.f"
/* L50: */
#line 641 "dlaqr0.f"
			}
#line 642 "dlaqr0.f"
L60:
#line 643 "dlaqr0.f"
			;
#line 643 "dlaqr0.f"
		    }

/*                 ==== Shuffle shifts into pairs of real shifts */
/*                 .    and pairs of complex conjugate shifts */
/*                 .    assuming complex conjugate shifts are */
/*                 .    already adjacent to one another. (Yes, */
/*                 .    they are.)  ==== */

#line 651 "dlaqr0.f"
		    i__2 = ks + 2;
#line 651 "dlaqr0.f"
		    for (i__ = kbot; i__ >= i__2; i__ += -2) {
#line 652 "dlaqr0.f"
			if (wi[i__] != -wi[i__ - 1]) {

#line 654 "dlaqr0.f"
			    swap = wr[i__];
#line 655 "dlaqr0.f"
			    wr[i__] = wr[i__ - 1];
#line 656 "dlaqr0.f"
			    wr[i__ - 1] = wr[i__ - 2];
#line 657 "dlaqr0.f"
			    wr[i__ - 2] = swap;

#line 659 "dlaqr0.f"
			    swap = wi[i__];
#line 660 "dlaqr0.f"
			    wi[i__] = wi[i__ - 1];
#line 661 "dlaqr0.f"
			    wi[i__ - 1] = wi[i__ - 2];
#line 662 "dlaqr0.f"
			    wi[i__ - 2] = swap;
#line 663 "dlaqr0.f"
			}
#line 664 "dlaqr0.f"
/* L70: */
#line 664 "dlaqr0.f"
		    }
#line 665 "dlaqr0.f"
		}

/*              ==== If there are only two shifts and both are */
/*              .    real, then use only one.  ==== */

#line 670 "dlaqr0.f"
		if (kbot - ks + 1 == 2) {
#line 671 "dlaqr0.f"
		    if (wi[kbot] == 0.) {
#line 672 "dlaqr0.f"
			if ((d__1 = wr[kbot] - h__[kbot + kbot * h_dim1], abs(
				d__1)) < (d__2 = wr[kbot - 1] - h__[kbot + 
				kbot * h_dim1], abs(d__2))) {
#line 674 "dlaqr0.f"
			    wr[kbot - 1] = wr[kbot];
#line 675 "dlaqr0.f"
			} else {
#line 676 "dlaqr0.f"
			    wr[kbot] = wr[kbot - 1];
#line 677 "dlaqr0.f"
			}
#line 678 "dlaqr0.f"
		    }
#line 679 "dlaqr0.f"
		}

/*              ==== Use up to NS of the the smallest magnatiude */
/*              .    shifts.  If there aren't NS shifts available, */
/*              .    then use them all, possibly dropping one to */
/*              .    make the number of shifts even. ==== */

/* Computing MIN */
#line 686 "dlaqr0.f"
		i__2 = ns, i__3 = kbot - ks + 1;
#line 686 "dlaqr0.f"
		ns = min(i__2,i__3);
#line 687 "dlaqr0.f"
		ns -= ns % 2;
#line 688 "dlaqr0.f"
		ks = kbot - ns + 1;

/*              ==== Small-bulge multi-shift QR sweep: */
/*              .    split workspace under the subdiagonal into */
/*              .    - a KDU-by-KDU work array U in the lower */
/*              .      left-hand-corner, */
/*              .    - a KDU-by-at-least-KDU-but-more-is-better */
/*              .      (KDU-by-NHo) horizontal work array WH along */
/*              .      the bottom edge, */
/*              .    - and an at-least-KDU-but-more-is-better-by-KDU */
/*              .      (NVE-by-KDU) vertical work WV arrow along */
/*              .      the left-hand-edge. ==== */

#line 701 "dlaqr0.f"
		kdu = ns * 3 - 3;
#line 702 "dlaqr0.f"
		ku = *n - kdu + 1;
#line 703 "dlaqr0.f"
		kwh = kdu + 1;
#line 704 "dlaqr0.f"
		nho = *n - kdu - 3 - (kdu + 1) + 1;
#line 705 "dlaqr0.f"
		kwv = kdu + 4;
#line 706 "dlaqr0.f"
		nve = *n - kdu - kwv + 1;

/*              ==== Small-bulge multi-shift QR sweep ==== */

#line 710 "dlaqr0.f"
		dlaqr5_(wantt, wantz, &kacc22, n, &ktop, &kbot, &ns, &wr[ks], 
			&wi[ks], &h__[h_offset], ldh, iloz, ihiz, &z__[
			z_offset], ldz, &work[1], &c__3, &h__[ku + h_dim1], 
			ldh, &nve, &h__[kwv + h_dim1], ldh, &nho, &h__[ku + 
			kwh * h_dim1], ldh);
#line 714 "dlaqr0.f"
	    }

/*           ==== Note progress (or the lack of it). ==== */

#line 718 "dlaqr0.f"
	    if (ld > 0) {
#line 719 "dlaqr0.f"
		ndfl = 1;
#line 720 "dlaqr0.f"
	    } else {
#line 721 "dlaqr0.f"
		++ndfl;
#line 722 "dlaqr0.f"
	    }

/*           ==== End of main loop ==== */
#line 725 "dlaqr0.f"
/* L80: */
#line 725 "dlaqr0.f"
	}

/*        ==== Iteration limit exceeded.  Set INFO to show where */
/*        .    the problem occurred and exit. ==== */

#line 730 "dlaqr0.f"
	*info = kbot;
#line 731 "dlaqr0.f"
L90:
#line 732 "dlaqr0.f"
	;
#line 732 "dlaqr0.f"
    }

/*     ==== Return the optimal value of LWORK. ==== */

#line 736 "dlaqr0.f"
    work[1] = (doublereal) lwkopt;

/*     ==== End of DLAQR0 ==== */

#line 740 "dlaqr0.f"
    return 0;
} /* dlaqr0_ */

