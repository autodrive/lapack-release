#line 1 "slaqr0.f"
/* slaqr0.f -- translated by f2c (version 20100827).
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

#line 1 "slaqr0.f"
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

/* > \brief \b SLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Sc
hur decomposition. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQR0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, */
/*                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               H( LDH, * ), WI( * ), WORK( * ), WR( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLAQR0 computes the eigenvalues of a Hessenberg matrix H */
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
/* >           previous call to SGEBAL, and then passed to SGEHRD when the */
/* >           matrix output by SGEBAL is reduced to Hessenberg form. */
/* >           Otherwise, ILO and IHI should be set to 1 and N, */
/* >           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N. */
/* >           If N = 0, then ILO = 1 and IHI = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is REAL array, dimension (LDH,N) */
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
/* >          WR is REAL array, dimension (IHI) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (IHI) */
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
/* >          Z is REAL array, dimension (LDZ,IHI) */
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
/* >          WORK is REAL array, dimension LWORK */
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
/* >           If LWORK = -1, then SLAQR0 does a workspace query. */
/* >           In this case, SLAQR0 checks the input parameters and */
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
/* >           .GT. 0:  if INFO = i, SLAQR0 failed to compute all of */
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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

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
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaqr0_(logical *wantt, logical *wantz, integer *n, 
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
    extern /* Subroutine */ int slanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), slaqr3_(
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    slaqr4_(logical *, logical *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), slaqr5_(logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *);
    static integer nibble;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char jbcmpz[2];
    extern /* Subroutine */ int slahqr_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), slacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    static integer nwupbd;
    static logical sorted;
    static integer lwkopt;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ================================================================ */
/*     .. Parameters .. */

/*     ==== Matrices of order NTINY or smaller must be processed by */
/*     .    SLAHQR because of insufficient subdiagonal scratch space. */
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
#line 324 "slaqr0.f"
    /* Parameter adjustments */
#line 324 "slaqr0.f"
    h_dim1 = *ldh;
#line 324 "slaqr0.f"
    h_offset = 1 + h_dim1;
#line 324 "slaqr0.f"
    h__ -= h_offset;
#line 324 "slaqr0.f"
    --wr;
#line 324 "slaqr0.f"
    --wi;
#line 324 "slaqr0.f"
    z_dim1 = *ldz;
#line 324 "slaqr0.f"
    z_offset = 1 + z_dim1;
#line 324 "slaqr0.f"
    z__ -= z_offset;
#line 324 "slaqr0.f"
    --work;
#line 324 "slaqr0.f"

#line 324 "slaqr0.f"
    /* Function Body */
#line 324 "slaqr0.f"
    *info = 0;

/*     ==== Quick return for N = 0: nothing to do. ==== */

#line 328 "slaqr0.f"
    if (*n == 0) {
#line 329 "slaqr0.f"
	work[1] = 1.;
#line 330 "slaqr0.f"
	return 0;
#line 331 "slaqr0.f"
    }

#line 333 "slaqr0.f"
    if (*n <= 11) {

/*        ==== Tiny matrices must use SLAHQR. ==== */

#line 337 "slaqr0.f"
	lwkopt = 1;
#line 338 "slaqr0.f"
	if (*lwork != -1) {
#line 338 "slaqr0.f"
	    slahqr_(wantt, wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &
		    wi[1], iloz, ihiz, &z__[z_offset], ldz, info);
#line 338 "slaqr0.f"
	}
#line 341 "slaqr0.f"
    } else {

/*        ==== Use small bulge multi-shift QR with aggressive early */
/*        .    deflation on larger-than-tiny matrices. ==== */

/*        ==== Hope for the best. ==== */

#line 348 "slaqr0.f"
	*info = 0;

/*        ==== Set up job flags for ILAENV. ==== */

#line 352 "slaqr0.f"
	if (*wantt) {
#line 353 "slaqr0.f"
	    *(unsigned char *)jbcmpz = 'S';
#line 354 "slaqr0.f"
	} else {
#line 355 "slaqr0.f"
	    *(unsigned char *)jbcmpz = 'E';
#line 356 "slaqr0.f"
	}
#line 357 "slaqr0.f"
	if (*wantz) {
#line 358 "slaqr0.f"
	    *(unsigned char *)&jbcmpz[1] = 'V';
#line 359 "slaqr0.f"
	} else {
#line 360 "slaqr0.f"
	    *(unsigned char *)&jbcmpz[1] = 'N';
#line 361 "slaqr0.f"
	}

/*        ==== NWR = recommended deflation window size.  At this */
/*        .    point,  N .GT. NTINY = 11, so there is enough */
/*        .    subdiagonal workspace for NWR.GE.2 as required. */
/*        .    (In fact, there is enough subdiagonal space for */
/*        .    NWR.GE.3.) ==== */

#line 369 "slaqr0.f"
	nwr = ilaenv_(&c__13, "SLAQR0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
#line 370 "slaqr0.f"
	nwr = max(2,nwr);
/* Computing MIN */
#line 371 "slaqr0.f"
	i__1 = *ihi - *ilo + 1, i__2 = (*n - 1) / 3, i__1 = min(i__1,i__2);
#line 371 "slaqr0.f"
	nwr = min(i__1,nwr);

/*        ==== NSR = recommended number of simultaneous shifts. */
/*        .    At this point N .GT. NTINY = 11, so there is at */
/*        .    enough subdiagonal workspace for NSR to be even */
/*        .    and greater than or equal to two as required. ==== */

#line 378 "slaqr0.f"
	nsr = ilaenv_(&c__15, "SLAQR0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
/* Computing MIN */
#line 379 "slaqr0.f"
	i__1 = nsr, i__2 = (*n + 6) / 9, i__1 = min(i__1,i__2), i__2 = *ihi - 
		*ilo;
#line 379 "slaqr0.f"
	nsr = min(i__1,i__2);
/* Computing MAX */
#line 380 "slaqr0.f"
	i__1 = 2, i__2 = nsr - nsr % 2;
#line 380 "slaqr0.f"
	nsr = max(i__1,i__2);

/*        ==== Estimate optimal workspace ==== */

/*        ==== Workspace query call to SLAQR3 ==== */

#line 386 "slaqr0.f"
	i__1 = nwr + 1;
#line 386 "slaqr0.f"
	slaqr3_(wantt, wantz, n, ilo, ihi, &i__1, &h__[h_offset], ldh, iloz, 
		ihiz, &z__[z_offset], ldz, &ls, &ld, &wr[1], &wi[1], &h__[
		h_offset], ldh, n, &h__[h_offset], ldh, n, &h__[h_offset], 
		ldh, &work[1], &c_n1);

/*        ==== Optimal workspace = MAX(SLAQR5, SLAQR3) ==== */

/* Computing MAX */
#line 392 "slaqr0.f"
	i__1 = nsr * 3 / 2, i__2 = (integer) work[1];
#line 392 "slaqr0.f"
	lwkopt = max(i__1,i__2);

/*        ==== Quick return in case of workspace query. ==== */

#line 396 "slaqr0.f"
	if (*lwork == -1) {
#line 397 "slaqr0.f"
	    work[1] = (doublereal) lwkopt;
#line 398 "slaqr0.f"
	    return 0;
#line 399 "slaqr0.f"
	}

/*        ==== SLAHQR/SLAQR0 crossover point ==== */

#line 403 "slaqr0.f"
	nmin = ilaenv_(&c__12, "SLAQR0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)
		6, (ftnlen)2);
#line 404 "slaqr0.f"
	nmin = max(11,nmin);

/*        ==== Nibble crossover point ==== */

#line 408 "slaqr0.f"
	nibble = ilaenv_(&c__14, "SLAQR0", jbcmpz, n, ilo, ihi, lwork, (
		ftnlen)6, (ftnlen)2);
#line 409 "slaqr0.f"
	nibble = max(0,nibble);

/*        ==== Accumulate reflections during ttswp?  Use block */
/*        .    2-by-2 structure during matrix-matrix multiply? ==== */

#line 414 "slaqr0.f"
	kacc22 = ilaenv_(&c__16, "SLAQR0", jbcmpz, n, ilo, ihi, lwork, (
		ftnlen)6, (ftnlen)2);
#line 415 "slaqr0.f"
	kacc22 = max(0,kacc22);
#line 416 "slaqr0.f"
	kacc22 = min(2,kacc22);

/*        ==== NWMAX = the largest possible deflation window for */
/*        .    which there is sufficient workspace. ==== */

/* Computing MIN */
#line 421 "slaqr0.f"
	i__1 = (*n - 1) / 3, i__2 = *lwork / 2;
#line 421 "slaqr0.f"
	nwmax = min(i__1,i__2);
#line 422 "slaqr0.f"
	nw = nwmax;

/*        ==== NSMAX = the Largest number of simultaneous shifts */
/*        .    for which there is sufficient workspace. ==== */

/* Computing MIN */
#line 427 "slaqr0.f"
	i__1 = (*n + 6) / 9, i__2 = (*lwork << 1) / 3;
#line 427 "slaqr0.f"
	nsmax = min(i__1,i__2);
#line 428 "slaqr0.f"
	nsmax -= nsmax % 2;

/*        ==== NDFL: an iteration count restarted at deflation. ==== */

#line 432 "slaqr0.f"
	ndfl = 1;

/*        ==== ITMAX = iteration limit ==== */

/* Computing MAX */
#line 436 "slaqr0.f"
	i__1 = 10, i__2 = *ihi - *ilo + 1;
#line 436 "slaqr0.f"
	itmax = max(i__1,i__2) * 30;

/*        ==== Last row and column in the active block ==== */

#line 440 "slaqr0.f"
	kbot = *ihi;

/*        ==== Main Loop ==== */

#line 444 "slaqr0.f"
	i__1 = itmax;
#line 444 "slaqr0.f"
	for (it = 1; it <= i__1; ++it) {

/*           ==== Done when KBOT falls below ILO ==== */

#line 448 "slaqr0.f"
	    if (kbot < *ilo) {
#line 448 "slaqr0.f"
		goto L90;
#line 448 "slaqr0.f"
	    }

/*           ==== Locate active block ==== */

#line 453 "slaqr0.f"
	    i__2 = *ilo + 1;
#line 453 "slaqr0.f"
	    for (k = kbot; k >= i__2; --k) {
#line 454 "slaqr0.f"
		if (h__[k + (k - 1) * h_dim1] == 0.) {
#line 454 "slaqr0.f"
		    goto L20;
#line 454 "slaqr0.f"
		}
#line 456 "slaqr0.f"
/* L10: */
#line 456 "slaqr0.f"
	    }
#line 457 "slaqr0.f"
	    k = *ilo;
#line 458 "slaqr0.f"
L20:
#line 459 "slaqr0.f"
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

#line 477 "slaqr0.f"
	    nh = kbot - ktop + 1;
#line 478 "slaqr0.f"
	    nwupbd = min(nh,nwmax);
#line 479 "slaqr0.f"
	    if (ndfl < 5) {
#line 480 "slaqr0.f"
		nw = min(nwupbd,nwr);
#line 481 "slaqr0.f"
	    } else {
/* Computing MIN */
#line 482 "slaqr0.f"
		i__2 = nwupbd, i__3 = nw << 1;
#line 482 "slaqr0.f"
		nw = min(i__2,i__3);
#line 483 "slaqr0.f"
	    }
#line 484 "slaqr0.f"
	    if (nw < nwmax) {
#line 485 "slaqr0.f"
		if (nw >= nh - 1) {
#line 486 "slaqr0.f"
		    nw = nh;
#line 487 "slaqr0.f"
		} else {
#line 488 "slaqr0.f"
		    kwtop = kbot - nw + 1;
#line 489 "slaqr0.f"
		    if ((d__1 = h__[kwtop + (kwtop - 1) * h_dim1], abs(d__1)) 
			    > (d__2 = h__[kwtop - 1 + (kwtop - 2) * h_dim1], 
			    abs(d__2))) {
#line 489 "slaqr0.f"
			++nw;
#line 489 "slaqr0.f"
		    }
#line 491 "slaqr0.f"
		}
#line 492 "slaqr0.f"
	    }
#line 493 "slaqr0.f"
	    if (ndfl < 5) {
#line 494 "slaqr0.f"
		ndec = -1;
#line 495 "slaqr0.f"
	    } else if (ndec >= 0 || nw >= nwupbd) {
#line 496 "slaqr0.f"
		++ndec;
#line 497 "slaqr0.f"
		if (nw - ndec < 2) {
#line 497 "slaqr0.f"
		    ndec = 0;
#line 497 "slaqr0.f"
		}
#line 499 "slaqr0.f"
		nw -= ndec;
#line 500 "slaqr0.f"
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

#line 513 "slaqr0.f"
	    kv = *n - nw + 1;
#line 514 "slaqr0.f"
	    kt = nw + 1;
#line 515 "slaqr0.f"
	    nho = *n - nw - 1 - kt + 1;
#line 516 "slaqr0.f"
	    kwv = nw + 2;
#line 517 "slaqr0.f"
	    nve = *n - nw - kwv + 1;

/*           ==== Aggressive early deflation ==== */

#line 521 "slaqr0.f"
	    slaqr3_(wantt, wantz, n, &ktop, &kbot, &nw, &h__[h_offset], ldh, 
		    iloz, ihiz, &z__[z_offset], ldz, &ls, &ld, &wr[1], &wi[1],
		     &h__[kv + h_dim1], ldh, &nho, &h__[kv + kt * h_dim1], 
		    ldh, &nve, &h__[kwv + h_dim1], ldh, &work[1], lwork);

/*           ==== Adjust KBOT accounting for new deflations. ==== */

#line 528 "slaqr0.f"
	    kbot -= ld;

/*           ==== KS points to the shifts. ==== */

#line 532 "slaqr0.f"
	    ks = kbot - ls + 1;

/*           ==== Skip an expensive QR sweep if there is a (partly */
/*           .    heuristic) reason to expect that many eigenvalues */
/*           .    will deflate without it.  Here, the QR sweep is */
/*           .    skipped if many eigenvalues have just been deflated */
/*           .    or if the remaining active block is small. */

#line 540 "slaqr0.f"
	    if (ld == 0 || ld * 100 <= nw * nibble && kbot - ktop + 1 > min(
		    nmin,nwmax)) {

/*              ==== NS = nominal number of simultaneous shifts. */
/*              .    This may be lowered (slightly) if SLAQR3 */
/*              .    did not provide that many shifts. ==== */

/* Computing MIN */
/* Computing MAX */
#line 547 "slaqr0.f"
		i__4 = 2, i__5 = kbot - ktop;
#line 547 "slaqr0.f"
		i__2 = min(nsmax,nsr), i__3 = max(i__4,i__5);
#line 547 "slaqr0.f"
		ns = min(i__2,i__3);
#line 548 "slaqr0.f"
		ns -= ns % 2;

/*              ==== If there have been no deflations */
/*              .    in a multiple of KEXSH iterations, */
/*              .    then try exceptional shifts. */
/*              .    Otherwise use shifts provided by */
/*              .    SLAQR3 above or from the eigenvalues */
/*              .    of a trailing principal submatrix. ==== */

#line 557 "slaqr0.f"
		if (ndfl % 6 == 0) {
#line 558 "slaqr0.f"
		    ks = kbot - ns + 1;
/* Computing MAX */
#line 559 "slaqr0.f"
		    i__3 = ks + 1, i__4 = ktop + 2;
#line 559 "slaqr0.f"
		    i__2 = max(i__3,i__4);
#line 559 "slaqr0.f"
		    for (i__ = kbot; i__ >= i__2; i__ += -2) {
#line 560 "slaqr0.f"
			ss = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1))
				 + (d__2 = h__[i__ - 1 + (i__ - 2) * h_dim1], 
				abs(d__2));
#line 561 "slaqr0.f"
			aa = ss * .75 + h__[i__ + i__ * h_dim1];
#line 562 "slaqr0.f"
			bb = ss;
#line 563 "slaqr0.f"
			cc = ss * -.4375;
#line 564 "slaqr0.f"
			dd = aa;
#line 565 "slaqr0.f"
			slanv2_(&aa, &bb, &cc, &dd, &wr[i__ - 1], &wi[i__ - 1]
				, &wr[i__], &wi[i__], &cs, &sn);
#line 567 "slaqr0.f"
/* L30: */
#line 567 "slaqr0.f"
		    }
#line 568 "slaqr0.f"
		    if (ks == ktop) {
#line 569 "slaqr0.f"
			wr[ks + 1] = h__[ks + 1 + (ks + 1) * h_dim1];
#line 570 "slaqr0.f"
			wi[ks + 1] = 0.;
#line 571 "slaqr0.f"
			wr[ks] = wr[ks + 1];
#line 572 "slaqr0.f"
			wi[ks] = wi[ks + 1];
#line 573 "slaqr0.f"
		    }
#line 574 "slaqr0.f"
		} else {

/*                 ==== Got NS/2 or fewer shifts? Use SLAQR4 or */
/*                 .    SLAHQR on a trailing principal submatrix to */
/*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9, */
/*                 .    there is enough space below the subdiagonal */
/*                 .    to fit an NS-by-NS scratch array.) ==== */

#line 582 "slaqr0.f"
		    if (kbot - ks + 1 <= ns / 2) {
#line 583 "slaqr0.f"
			ks = kbot - ns + 1;
#line 584 "slaqr0.f"
			kt = *n - ns + 1;
#line 585 "slaqr0.f"
			slacpy_("A", &ns, &ns, &h__[ks + ks * h_dim1], ldh, &
				h__[kt + h_dim1], ldh, (ftnlen)1);
#line 587 "slaqr0.f"
			if (ns > nmin) {
#line 588 "slaqr0.f"
			    slaqr4_(&c_false, &c_false, &ns, &c__1, &ns, &h__[
				    kt + h_dim1], ldh, &wr[ks], &wi[ks], &
				    c__1, &c__1, zdum, &c__1, &work[1], lwork,
				     &inf);
#line 592 "slaqr0.f"
			} else {
#line 593 "slaqr0.f"
			    slahqr_(&c_false, &c_false, &ns, &c__1, &ns, &h__[
				    kt + h_dim1], ldh, &wr[ks], &wi[ks], &
				    c__1, &c__1, zdum, &c__1, &inf);
#line 596 "slaqr0.f"
			}
#line 597 "slaqr0.f"
			ks += inf;

/*                    ==== In case of a rare QR failure use */
/*                    .    eigenvalues of the trailing 2-by-2 */
/*                    .    principal submatrix.  ==== */

#line 603 "slaqr0.f"
			if (ks >= kbot) {
#line 604 "slaqr0.f"
			    aa = h__[kbot - 1 + (kbot - 1) * h_dim1];
#line 605 "slaqr0.f"
			    cc = h__[kbot + (kbot - 1) * h_dim1];
#line 606 "slaqr0.f"
			    bb = h__[kbot - 1 + kbot * h_dim1];
#line 607 "slaqr0.f"
			    dd = h__[kbot + kbot * h_dim1];
#line 608 "slaqr0.f"
			    slanv2_(&aa, &bb, &cc, &dd, &wr[kbot - 1], &wi[
				    kbot - 1], &wr[kbot], &wi[kbot], &cs, &sn)
				    ;
#line 611 "slaqr0.f"
			    ks = kbot - 1;
#line 612 "slaqr0.f"
			}
#line 613 "slaqr0.f"
		    }

#line 615 "slaqr0.f"
		    if (kbot - ks + 1 > ns) {

/*                    ==== Sort the shifts (Helps a little) */
/*                    .    Bubble sort keeps complex conjugate */
/*                    .    pairs together. ==== */

#line 621 "slaqr0.f"
			sorted = FALSE_;
#line 622 "slaqr0.f"
			i__2 = ks + 1;
#line 622 "slaqr0.f"
			for (k = kbot; k >= i__2; --k) {
#line 623 "slaqr0.f"
			    if (sorted) {
#line 623 "slaqr0.f"
				goto L60;
#line 623 "slaqr0.f"
			    }
#line 625 "slaqr0.f"
			    sorted = TRUE_;
#line 626 "slaqr0.f"
			    i__3 = k - 1;
#line 626 "slaqr0.f"
			    for (i__ = ks; i__ <= i__3; ++i__) {
#line 627 "slaqr0.f"
				if ((d__1 = wr[i__], abs(d__1)) + (d__2 = wi[
					i__], abs(d__2)) < (d__3 = wr[i__ + 1]
					, abs(d__3)) + (d__4 = wi[i__ + 1], 
					abs(d__4))) {
#line 629 "slaqr0.f"
				    sorted = FALSE_;

#line 631 "slaqr0.f"
				    swap = wr[i__];
#line 632 "slaqr0.f"
				    wr[i__] = wr[i__ + 1];
#line 633 "slaqr0.f"
				    wr[i__ + 1] = swap;

#line 635 "slaqr0.f"
				    swap = wi[i__];
#line 636 "slaqr0.f"
				    wi[i__] = wi[i__ + 1];
#line 637 "slaqr0.f"
				    wi[i__ + 1] = swap;
#line 638 "slaqr0.f"
				}
#line 639 "slaqr0.f"
/* L40: */
#line 639 "slaqr0.f"
			    }
#line 640 "slaqr0.f"
/* L50: */
#line 640 "slaqr0.f"
			}
#line 641 "slaqr0.f"
L60:
#line 642 "slaqr0.f"
			;
#line 642 "slaqr0.f"
		    }

/*                 ==== Shuffle shifts into pairs of real shifts */
/*                 .    and pairs of complex conjugate shifts */
/*                 .    assuming complex conjugate shifts are */
/*                 .    already adjacent to one another. (Yes, */
/*                 .    they are.)  ==== */

#line 650 "slaqr0.f"
		    i__2 = ks + 2;
#line 650 "slaqr0.f"
		    for (i__ = kbot; i__ >= i__2; i__ += -2) {
#line 651 "slaqr0.f"
			if (wi[i__] != -wi[i__ - 1]) {

#line 653 "slaqr0.f"
			    swap = wr[i__];
#line 654 "slaqr0.f"
			    wr[i__] = wr[i__ - 1];
#line 655 "slaqr0.f"
			    wr[i__ - 1] = wr[i__ - 2];
#line 656 "slaqr0.f"
			    wr[i__ - 2] = swap;

#line 658 "slaqr0.f"
			    swap = wi[i__];
#line 659 "slaqr0.f"
			    wi[i__] = wi[i__ - 1];
#line 660 "slaqr0.f"
			    wi[i__ - 1] = wi[i__ - 2];
#line 661 "slaqr0.f"
			    wi[i__ - 2] = swap;
#line 662 "slaqr0.f"
			}
#line 663 "slaqr0.f"
/* L70: */
#line 663 "slaqr0.f"
		    }
#line 664 "slaqr0.f"
		}

/*              ==== If there are only two shifts and both are */
/*              .    real, then use only one.  ==== */

#line 669 "slaqr0.f"
		if (kbot - ks + 1 == 2) {
#line 670 "slaqr0.f"
		    if (wi[kbot] == 0.) {
#line 671 "slaqr0.f"
			if ((d__1 = wr[kbot] - h__[kbot + kbot * h_dim1], abs(
				d__1)) < (d__2 = wr[kbot - 1] - h__[kbot + 
				kbot * h_dim1], abs(d__2))) {
#line 673 "slaqr0.f"
			    wr[kbot - 1] = wr[kbot];
#line 674 "slaqr0.f"
			} else {
#line 675 "slaqr0.f"
			    wr[kbot] = wr[kbot - 1];
#line 676 "slaqr0.f"
			}
#line 677 "slaqr0.f"
		    }
#line 678 "slaqr0.f"
		}

/*              ==== Use up to NS of the the smallest magnatiude */
/*              .    shifts.  If there aren't NS shifts available, */
/*              .    then use them all, possibly dropping one to */
/*              .    make the number of shifts even. ==== */

/* Computing MIN */
#line 685 "slaqr0.f"
		i__2 = ns, i__3 = kbot - ks + 1;
#line 685 "slaqr0.f"
		ns = min(i__2,i__3);
#line 686 "slaqr0.f"
		ns -= ns % 2;
#line 687 "slaqr0.f"
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

#line 700 "slaqr0.f"
		kdu = ns * 3 - 3;
#line 701 "slaqr0.f"
		ku = *n - kdu + 1;
#line 702 "slaqr0.f"
		kwh = kdu + 1;
#line 703 "slaqr0.f"
		nho = *n - kdu - 3 - (kdu + 1) + 1;
#line 704 "slaqr0.f"
		kwv = kdu + 4;
#line 705 "slaqr0.f"
		nve = *n - kdu - kwv + 1;

/*              ==== Small-bulge multi-shift QR sweep ==== */

#line 709 "slaqr0.f"
		slaqr5_(wantt, wantz, &kacc22, n, &ktop, &kbot, &ns, &wr[ks], 
			&wi[ks], &h__[h_offset], ldh, iloz, ihiz, &z__[
			z_offset], ldz, &work[1], &c__3, &h__[ku + h_dim1], 
			ldh, &nve, &h__[kwv + h_dim1], ldh, &nho, &h__[ku + 
			kwh * h_dim1], ldh);
#line 713 "slaqr0.f"
	    }

/*           ==== Note progress (or the lack of it). ==== */

#line 717 "slaqr0.f"
	    if (ld > 0) {
#line 718 "slaqr0.f"
		ndfl = 1;
#line 719 "slaqr0.f"
	    } else {
#line 720 "slaqr0.f"
		++ndfl;
#line 721 "slaqr0.f"
	    }

/*           ==== End of main loop ==== */
#line 724 "slaqr0.f"
/* L80: */
#line 724 "slaqr0.f"
	}

/*        ==== Iteration limit exceeded.  Set INFO to show where */
/*        .    the problem occurred and exit. ==== */

#line 729 "slaqr0.f"
	*info = kbot;
#line 730 "slaqr0.f"
L90:
#line 731 "slaqr0.f"
	;
#line 731 "slaqr0.f"
    }

/*     ==== Return the optimal value of LWORK. ==== */

#line 735 "slaqr0.f"
    work[1] = (doublereal) lwkopt;

/*     ==== End of SLAQR0 ==== */

#line 739 "slaqr0.f"
    return 0;
} /* slaqr0_ */

