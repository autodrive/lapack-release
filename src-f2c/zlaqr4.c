#line 1 "zlaqr4.f"
/* zlaqr4.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqr4.f"
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

/* > \brief \b ZLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Sc
hur decomposition. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQR4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, */
/*                          IHIZ, Z, LDZ, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLAQR4 implements one level of recursion for ZLAQR0. */
/* >    It is a complete implementation of the small bulge multi-shift */
/* >    QR algorithm.  It may be called by ZLAQR0 and, for large enough */
/* >    deflation window size, it may be called by ZLAQR3.  This */
/* >    subroutine is identical to ZLAQR0 except that it calls ZLAQR2 */
/* >    instead of ZLAQR3. */
/* > */
/* >    ZLAQR4 computes the eigenvalues of a Hessenberg matrix H */
/* >    and, optionally, the matrices T and Z from the Schur decomposition */
/* >    H = Z T Z**H, where T is an upper triangular matrix (the */
/* >    Schur form), and Z is the unitary matrix of Schur vectors. */
/* > */
/* >    Optionally Z may be postmultiplied into an input unitary */
/* >    matrix Q so that this routine can give the Schur factorization */
/* >    of a matrix A which has been reduced to the Hessenberg form H */
/* >    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H. */
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
/* >           previous call to ZGEBAL, and then passed to ZGEHRD when the */
/* >           matrix output by ZGEBAL is reduced to Hessenberg form. */
/* >           Otherwise, ILO and IHI should be set to 1 and N, */
/* >           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N. */
/* >           If N = 0, then ILO = 1 and IHI = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX*16 array, dimension (LDH,N) */
/* >           On entry, the upper Hessenberg matrix H. */
/* >           On exit, if INFO = 0 and WANTT is .TRUE., then H */
/* >           contains the upper triangular matrix T from the Schur */
/* >           decomposition (the Schur form). If INFO = 0 and WANT is */
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
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX*16 array, dimension (N) */
/* >           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored */
/* >           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are */
/* >           stored in the same order as on the diagonal of the Schur */
/* >           form returned in H, with W(i) = H(i,i). */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ,IHI) */
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
/* >          WORK is COMPLEX*16 array, dimension LWORK */
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
/* >           If LWORK = -1, then ZLAQR4 does a workspace query. */
/* >           In this case, ZLAQR4 checks the input parameters and */
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
/* >           .GT. 0:  if INFO = i, ZLAQR4 failed to compute all of */
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
/* >                where U is a unitary matrix.  The final */
/* >                value of  H is upper Hessenberg and triangular in */
/* >                rows and columns INFO+1 through IHI. */
/* > */
/* >                If INFO .GT. 0 and WANTZ is .TRUE., then on exit */
/* > */
/* >                  (final value of Z(ILO:IHI,ILOZ:IHIZ) */
/* >                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U */
/* > */
/* >                where U is the unitary matrix in (*) (regard- */
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

/* > \ingroup complex16OTHERauxiliary */

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
/* Subroutine */ int zlaqr4_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_sqrt(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k;
    static doublereal s;
    static doublecomplex aa, bb, cc, dd;
    static integer ld, nh, it, ks, kt, ku, kv, ls, ns, nw;
    static doublecomplex tr2, det;
    static integer inf, kdu, nho, nve, kwh, nsr, nwr, kwv, ndec, ndfl, kbot, 
	    nmin;
    static doublecomplex swap;
    static integer ktop;
    static doublecomplex zdum[1]	/* was [1][1] */;
    static integer kacc22, itmax, nsmax, nwmax, kwtop;
    extern /* Subroutine */ int zlaqr2_(logical *, logical *, integer *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *,
	     doublecomplex *, integer *, integer *, doublecomplex *, integer *
	    , doublecomplex *, integer *), zlaqr5_(logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, doublecomplex *, integer *);
    static integer nibble;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char jbcmpz[2];
    static doublecomplex rtdisc;
    static integer nwupbd;
    static logical sorted;
    extern /* Subroutine */ int zlahqr_(logical *, logical *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, integer *, integer *), 
	    zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen);
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
/*     .    ZLAHQR because of insufficient subdiagonal scratch space. */
/*     .    (This is a hard limit.) ==== */

/*     ==== Exceptional deflation windows:  try to cure rare */
/*     .    slow convergence by varying the size of the */
/*     .    deflation window after KEXNW iterations. ==== */

/*     ==== Exceptional shifts: try to cure rare slow convergence */
/*     .    with ad-hoc exceptional shifts every KEXSH iterations. */
/*     .    ==== */

/*     ==== The constant WILK1 is used to form the exceptional */
/*     .    shifts. ==== */
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 326 "zlaqr4.f"
    /* Parameter adjustments */
#line 326 "zlaqr4.f"
    h_dim1 = *ldh;
#line 326 "zlaqr4.f"
    h_offset = 1 + h_dim1;
#line 326 "zlaqr4.f"
    h__ -= h_offset;
#line 326 "zlaqr4.f"
    --w;
#line 326 "zlaqr4.f"
    z_dim1 = *ldz;
#line 326 "zlaqr4.f"
    z_offset = 1 + z_dim1;
#line 326 "zlaqr4.f"
    z__ -= z_offset;
#line 326 "zlaqr4.f"
    --work;
#line 326 "zlaqr4.f"

#line 326 "zlaqr4.f"
    /* Function Body */
#line 326 "zlaqr4.f"
    *info = 0;

/*     ==== Quick return for N = 0: nothing to do. ==== */

#line 330 "zlaqr4.f"
    if (*n == 0) {
#line 331 "zlaqr4.f"
	work[1].r = 1., work[1].i = 0.;
#line 332 "zlaqr4.f"
	return 0;
#line 333 "zlaqr4.f"
    }

#line 335 "zlaqr4.f"
    if (*n <= 11) {

/*        ==== Tiny matrices must use ZLAHQR. ==== */

#line 339 "zlaqr4.f"
	lwkopt = 1;
#line 340 "zlaqr4.f"
	if (*lwork != -1) {
#line 340 "zlaqr4.f"
	    zlahqr_(wantt, wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], 
		    iloz, ihiz, &z__[z_offset], ldz, info);
#line 340 "zlaqr4.f"
	}
#line 343 "zlaqr4.f"
    } else {

/*        ==== Use small bulge multi-shift QR with aggressive early */
/*        .    deflation on larger-than-tiny matrices. ==== */

/*        ==== Hope for the best. ==== */

#line 350 "zlaqr4.f"
	*info = 0;

/*        ==== Set up job flags for ILAENV. ==== */

#line 354 "zlaqr4.f"
	if (*wantt) {
#line 355 "zlaqr4.f"
	    *(unsigned char *)jbcmpz = 'S';
#line 356 "zlaqr4.f"
	} else {
#line 357 "zlaqr4.f"
	    *(unsigned char *)jbcmpz = 'E';
#line 358 "zlaqr4.f"
	}
#line 359 "zlaqr4.f"
	if (*wantz) {
#line 360 "zlaqr4.f"
	    *(unsigned char *)&jbcmpz[1] = 'V';
#line 361 "zlaqr4.f"
	} else {
#line 362 "zlaqr4.f"
	    *(unsigned char *)&jbcmpz[1] = 'N';
#line 363 "zlaqr4.f"
	}

/*        ==== NWR = recommended deflation window size.  At this */
/*        .    point,  N .GT. NTINY = 11, so there is enough */
/*        .    subdiagonal workspace for NWR.GE.2 as required. */
/*        .    (In fact, there is enough subdiagonal space for */
/*        .    NWR.GE.3.) ==== */

#line 371 "zlaqr4.f"
	nwr = ilaenv_(&c__13, "ZLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
#line 372 "zlaqr4.f"
	nwr = max(2,nwr);
/* Computing MIN */
#line 373 "zlaqr4.f"
	i__1 = *ihi - *ilo + 1, i__2 = (*n - 1) / 3, i__1 = min(i__1,i__2);
#line 373 "zlaqr4.f"
	nwr = min(i__1,nwr);

/*        ==== NSR = recommended number of simultaneous shifts. */
/*        .    At this point N .GT. NTINY = 11, so there is at */
/*        .    enough subdiagonal workspace for NSR to be even */
/*        .    and greater than or equal to two as required. ==== */

#line 380 "zlaqr4.f"
	nsr = ilaenv_(&c__15, "ZLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
/* Computing MIN */
#line 381 "zlaqr4.f"
	i__1 = nsr, i__2 = (*n + 6) / 9, i__1 = min(i__1,i__2), i__2 = *ihi - 
		*ilo;
#line 381 "zlaqr4.f"
	nsr = min(i__1,i__2);
/* Computing MAX */
#line 382 "zlaqr4.f"
	i__1 = 2, i__2 = nsr - nsr % 2;
#line 382 "zlaqr4.f"
	nsr = max(i__1,i__2);

/*        ==== Estimate optimal workspace ==== */

/*        ==== Workspace query call to ZLAQR2 ==== */

#line 388 "zlaqr4.f"
	i__1 = nwr + 1;
#line 388 "zlaqr4.f"
	zlaqr2_(wantt, wantz, n, ilo, ihi, &i__1, &h__[h_offset], ldh, iloz, 
		ihiz, &z__[z_offset], ldz, &ls, &ld, &w[1], &h__[h_offset], 
		ldh, n, &h__[h_offset], ldh, n, &h__[h_offset], ldh, &work[1],
		 &c_n1);

/*        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR2) ==== */

/* Computing MAX */
#line 394 "zlaqr4.f"
	i__1 = nsr * 3 / 2, i__2 = (integer) work[1].r;
#line 394 "zlaqr4.f"
	lwkopt = max(i__1,i__2);

/*        ==== Quick return in case of workspace query. ==== */

#line 398 "zlaqr4.f"
	if (*lwork == -1) {
#line 399 "zlaqr4.f"
	    d__1 = (doublereal) lwkopt;
#line 399 "zlaqr4.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 399 "zlaqr4.f"
	    work[1].r = z__1.r, work[1].i = z__1.i;
#line 400 "zlaqr4.f"
	    return 0;
#line 401 "zlaqr4.f"
	}

/*        ==== ZLAHQR/ZLAQR0 crossover point ==== */

#line 405 "zlaqr4.f"
	nmin = ilaenv_(&c__12, "ZLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)
		6, (ftnlen)2);
#line 406 "zlaqr4.f"
	nmin = max(11,nmin);

/*        ==== Nibble crossover point ==== */

#line 410 "zlaqr4.f"
	nibble = ilaenv_(&c__14, "ZLAQR4", jbcmpz, n, ilo, ihi, lwork, (
		ftnlen)6, (ftnlen)2);
#line 411 "zlaqr4.f"
	nibble = max(0,nibble);

/*        ==== Accumulate reflections during ttswp?  Use block */
/*        .    2-by-2 structure during matrix-matrix multiply? ==== */

#line 416 "zlaqr4.f"
	kacc22 = ilaenv_(&c__16, "ZLAQR4", jbcmpz, n, ilo, ihi, lwork, (
		ftnlen)6, (ftnlen)2);
#line 417 "zlaqr4.f"
	kacc22 = max(0,kacc22);
#line 418 "zlaqr4.f"
	kacc22 = min(2,kacc22);

/*        ==== NWMAX = the largest possible deflation window for */
/*        .    which there is sufficient workspace. ==== */

/* Computing MIN */
#line 423 "zlaqr4.f"
	i__1 = (*n - 1) / 3, i__2 = *lwork / 2;
#line 423 "zlaqr4.f"
	nwmax = min(i__1,i__2);
#line 424 "zlaqr4.f"
	nw = nwmax;

/*        ==== NSMAX = the Largest number of simultaneous shifts */
/*        .    for which there is sufficient workspace. ==== */

/* Computing MIN */
#line 429 "zlaqr4.f"
	i__1 = (*n + 6) / 9, i__2 = (*lwork << 1) / 3;
#line 429 "zlaqr4.f"
	nsmax = min(i__1,i__2);
#line 430 "zlaqr4.f"
	nsmax -= nsmax % 2;

/*        ==== NDFL: an iteration count restarted at deflation. ==== */

#line 434 "zlaqr4.f"
	ndfl = 1;

/*        ==== ITMAX = iteration limit ==== */

/* Computing MAX */
#line 438 "zlaqr4.f"
	i__1 = 10, i__2 = *ihi - *ilo + 1;
#line 438 "zlaqr4.f"
	itmax = max(i__1,i__2) * 30;

/*        ==== Last row and column in the active block ==== */

#line 442 "zlaqr4.f"
	kbot = *ihi;

/*        ==== Main Loop ==== */

#line 446 "zlaqr4.f"
	i__1 = itmax;
#line 446 "zlaqr4.f"
	for (it = 1; it <= i__1; ++it) {

/*           ==== Done when KBOT falls below ILO ==== */

#line 450 "zlaqr4.f"
	    if (kbot < *ilo) {
#line 450 "zlaqr4.f"
		goto L80;
#line 450 "zlaqr4.f"
	    }

/*           ==== Locate active block ==== */

#line 455 "zlaqr4.f"
	    i__2 = *ilo + 1;
#line 455 "zlaqr4.f"
	    for (k = kbot; k >= i__2; --k) {
#line 456 "zlaqr4.f"
		i__3 = k + (k - 1) * h_dim1;
#line 456 "zlaqr4.f"
		if (h__[i__3].r == 0. && h__[i__3].i == 0.) {
#line 456 "zlaqr4.f"
		    goto L20;
#line 456 "zlaqr4.f"
		}
#line 458 "zlaqr4.f"
/* L10: */
#line 458 "zlaqr4.f"
	    }
#line 459 "zlaqr4.f"
	    k = *ilo;
#line 460 "zlaqr4.f"
L20:
#line 461 "zlaqr4.f"
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

#line 479 "zlaqr4.f"
	    nh = kbot - ktop + 1;
#line 480 "zlaqr4.f"
	    nwupbd = min(nh,nwmax);
#line 481 "zlaqr4.f"
	    if (ndfl < 5) {
#line 482 "zlaqr4.f"
		nw = min(nwupbd,nwr);
#line 483 "zlaqr4.f"
	    } else {
/* Computing MIN */
#line 484 "zlaqr4.f"
		i__2 = nwupbd, i__3 = nw << 1;
#line 484 "zlaqr4.f"
		nw = min(i__2,i__3);
#line 485 "zlaqr4.f"
	    }
#line 486 "zlaqr4.f"
	    if (nw < nwmax) {
#line 487 "zlaqr4.f"
		if (nw >= nh - 1) {
#line 488 "zlaqr4.f"
		    nw = nh;
#line 489 "zlaqr4.f"
		} else {
#line 490 "zlaqr4.f"
		    kwtop = kbot - nw + 1;
#line 491 "zlaqr4.f"
		    i__2 = kwtop + (kwtop - 1) * h_dim1;
#line 491 "zlaqr4.f"
		    i__3 = kwtop - 1 + (kwtop - 2) * h_dim1;
#line 491 "zlaqr4.f"
		    if ((d__1 = h__[i__2].r, abs(d__1)) + (d__2 = d_imag(&h__[
			    kwtop + (kwtop - 1) * h_dim1]), abs(d__2)) > (
			    d__3 = h__[i__3].r, abs(d__3)) + (d__4 = d_imag(&
			    h__[kwtop - 1 + (kwtop - 2) * h_dim1]), abs(d__4))
			    ) {
#line 491 "zlaqr4.f"
			++nw;
#line 491 "zlaqr4.f"
		    }
#line 493 "zlaqr4.f"
		}
#line 494 "zlaqr4.f"
	    }
#line 495 "zlaqr4.f"
	    if (ndfl < 5) {
#line 496 "zlaqr4.f"
		ndec = -1;
#line 497 "zlaqr4.f"
	    } else if (ndec >= 0 || nw >= nwupbd) {
#line 498 "zlaqr4.f"
		++ndec;
#line 499 "zlaqr4.f"
		if (nw - ndec < 2) {
#line 499 "zlaqr4.f"
		    ndec = 0;
#line 499 "zlaqr4.f"
		}
#line 501 "zlaqr4.f"
		nw -= ndec;
#line 502 "zlaqr4.f"
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

#line 515 "zlaqr4.f"
	    kv = *n - nw + 1;
#line 516 "zlaqr4.f"
	    kt = nw + 1;
#line 517 "zlaqr4.f"
	    nho = *n - nw - 1 - kt + 1;
#line 518 "zlaqr4.f"
	    kwv = nw + 2;
#line 519 "zlaqr4.f"
	    nve = *n - nw - kwv + 1;

/*           ==== Aggressive early deflation ==== */

#line 523 "zlaqr4.f"
	    zlaqr2_(wantt, wantz, n, &ktop, &kbot, &nw, &h__[h_offset], ldh, 
		    iloz, ihiz, &z__[z_offset], ldz, &ls, &ld, &w[1], &h__[kv 
		    + h_dim1], ldh, &nho, &h__[kv + kt * h_dim1], ldh, &nve, &
		    h__[kwv + h_dim1], ldh, &work[1], lwork);

/*           ==== Adjust KBOT accounting for new deflations. ==== */

#line 530 "zlaqr4.f"
	    kbot -= ld;

/*           ==== KS points to the shifts. ==== */

#line 534 "zlaqr4.f"
	    ks = kbot - ls + 1;

/*           ==== Skip an expensive QR sweep if there is a (partly */
/*           .    heuristic) reason to expect that many eigenvalues */
/*           .    will deflate without it.  Here, the QR sweep is */
/*           .    skipped if many eigenvalues have just been deflated */
/*           .    or if the remaining active block is small. */

#line 542 "zlaqr4.f"
	    if (ld == 0 || ld * 100 <= nw * nibble && kbot - ktop + 1 > min(
		    nmin,nwmax)) {

/*              ==== NS = nominal number of simultaneous shifts. */
/*              .    This may be lowered (slightly) if ZLAQR2 */
/*              .    did not provide that many shifts. ==== */

/* Computing MIN */
/* Computing MAX */
#line 549 "zlaqr4.f"
		i__4 = 2, i__5 = kbot - ktop;
#line 549 "zlaqr4.f"
		i__2 = min(nsmax,nsr), i__3 = max(i__4,i__5);
#line 549 "zlaqr4.f"
		ns = min(i__2,i__3);
#line 550 "zlaqr4.f"
		ns -= ns % 2;

/*              ==== If there have been no deflations */
/*              .    in a multiple of KEXSH iterations, */
/*              .    then try exceptional shifts. */
/*              .    Otherwise use shifts provided by */
/*              .    ZLAQR2 above or from the eigenvalues */
/*              .    of a trailing principal submatrix. ==== */

#line 559 "zlaqr4.f"
		if (ndfl % 6 == 0) {
#line 560 "zlaqr4.f"
		    ks = kbot - ns + 1;
#line 561 "zlaqr4.f"
		    i__2 = ks + 1;
#line 561 "zlaqr4.f"
		    for (i__ = kbot; i__ >= i__2; i__ += -2) {
#line 562 "zlaqr4.f"
			i__3 = i__;
#line 562 "zlaqr4.f"
			i__4 = i__ + i__ * h_dim1;
#line 562 "zlaqr4.f"
			i__5 = i__ + (i__ - 1) * h_dim1;
#line 562 "zlaqr4.f"
			d__3 = ((d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				d_imag(&h__[i__ + (i__ - 1) * h_dim1]), abs(
				d__2))) * .75;
#line 562 "zlaqr4.f"
			z__1.r = h__[i__4].r + d__3, z__1.i = h__[i__4].i;
#line 562 "zlaqr4.f"
			w[i__3].r = z__1.r, w[i__3].i = z__1.i;
#line 563 "zlaqr4.f"
			i__3 = i__ - 1;
#line 563 "zlaqr4.f"
			i__4 = i__;
#line 563 "zlaqr4.f"
			w[i__3].r = w[i__4].r, w[i__3].i = w[i__4].i;
#line 564 "zlaqr4.f"
/* L30: */
#line 564 "zlaqr4.f"
		    }
#line 565 "zlaqr4.f"
		} else {

/*                 ==== Got NS/2 or fewer shifts? Use ZLAHQR */
/*                 .    on a trailing principal submatrix to */
/*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9, */
/*                 .    there is enough space below the subdiagonal */
/*                 .    to fit an NS-by-NS scratch array.) ==== */

#line 573 "zlaqr4.f"
		    if (kbot - ks + 1 <= ns / 2) {
#line 574 "zlaqr4.f"
			ks = kbot - ns + 1;
#line 575 "zlaqr4.f"
			kt = *n - ns + 1;
#line 576 "zlaqr4.f"
			zlacpy_("A", &ns, &ns, &h__[ks + ks * h_dim1], ldh, &
				h__[kt + h_dim1], ldh, (ftnlen)1);
#line 578 "zlaqr4.f"
			zlahqr_(&c_false, &c_false, &ns, &c__1, &ns, &h__[kt 
				+ h_dim1], ldh, &w[ks], &c__1, &c__1, zdum, &
				c__1, &inf);
#line 581 "zlaqr4.f"
			ks += inf;

/*                    ==== In case of a rare QR failure use */
/*                    .    eigenvalues of the trailing 2-by-2 */
/*                    .    principal submatrix.  Scale to avoid */
/*                    .    overflows, underflows and subnormals. */
/*                    .    (The scale factor S can not be zero, */
/*                    .    because H(KBOT,KBOT-1) is nonzero.) ==== */

#line 590 "zlaqr4.f"
			if (ks >= kbot) {
#line 591 "zlaqr4.f"
			    i__2 = kbot - 1 + (kbot - 1) * h_dim1;
#line 591 "zlaqr4.f"
			    i__3 = kbot + (kbot - 1) * h_dim1;
#line 591 "zlaqr4.f"
			    i__4 = kbot - 1 + kbot * h_dim1;
#line 591 "zlaqr4.f"
			    i__5 = kbot + kbot * h_dim1;
#line 591 "zlaqr4.f"
			    s = (d__1 = h__[i__2].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[kbot - 1 + (kbot - 1) * 
				    h_dim1]), abs(d__2)) + ((d__3 = h__[i__3]
				    .r, abs(d__3)) + (d__4 = d_imag(&h__[kbot 
				    + (kbot - 1) * h_dim1]), abs(d__4))) + ((
				    d__5 = h__[i__4].r, abs(d__5)) + (d__6 = 
				    d_imag(&h__[kbot - 1 + kbot * h_dim1]), 
				    abs(d__6))) + ((d__7 = h__[i__5].r, abs(
				    d__7)) + (d__8 = d_imag(&h__[kbot + kbot *
				     h_dim1]), abs(d__8)));
#line 595 "zlaqr4.f"
			    i__2 = kbot - 1 + (kbot - 1) * h_dim1;
#line 595 "zlaqr4.f"
			    z__1.r = h__[i__2].r / s, z__1.i = h__[i__2].i / 
				    s;
#line 595 "zlaqr4.f"
			    aa.r = z__1.r, aa.i = z__1.i;
#line 596 "zlaqr4.f"
			    i__2 = kbot + (kbot - 1) * h_dim1;
#line 596 "zlaqr4.f"
			    z__1.r = h__[i__2].r / s, z__1.i = h__[i__2].i / 
				    s;
#line 596 "zlaqr4.f"
			    cc.r = z__1.r, cc.i = z__1.i;
#line 597 "zlaqr4.f"
			    i__2 = kbot - 1 + kbot * h_dim1;
#line 597 "zlaqr4.f"
			    z__1.r = h__[i__2].r / s, z__1.i = h__[i__2].i / 
				    s;
#line 597 "zlaqr4.f"
			    bb.r = z__1.r, bb.i = z__1.i;
#line 598 "zlaqr4.f"
			    i__2 = kbot + kbot * h_dim1;
#line 598 "zlaqr4.f"
			    z__1.r = h__[i__2].r / s, z__1.i = h__[i__2].i / 
				    s;
#line 598 "zlaqr4.f"
			    dd.r = z__1.r, dd.i = z__1.i;
#line 599 "zlaqr4.f"
			    z__2.r = aa.r + dd.r, z__2.i = aa.i + dd.i;
#line 599 "zlaqr4.f"
			    z__1.r = z__2.r / 2., z__1.i = z__2.i / 2.;
#line 599 "zlaqr4.f"
			    tr2.r = z__1.r, tr2.i = z__1.i;
#line 600 "zlaqr4.f"
			    z__3.r = aa.r - tr2.r, z__3.i = aa.i - tr2.i;
#line 600 "zlaqr4.f"
			    z__4.r = dd.r - tr2.r, z__4.i = dd.i - tr2.i;
#line 600 "zlaqr4.f"
			    z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, 
				    z__2.i = z__3.r * z__4.i + z__3.i * 
				    z__4.r;
#line 600 "zlaqr4.f"
			    z__5.r = bb.r * cc.r - bb.i * cc.i, z__5.i = bb.r 
				    * cc.i + bb.i * cc.r;
#line 600 "zlaqr4.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 600 "zlaqr4.f"
			    det.r = z__1.r, det.i = z__1.i;
#line 601 "zlaqr4.f"
			    z__2.r = -det.r, z__2.i = -det.i;
#line 601 "zlaqr4.f"
			    z_sqrt(&z__1, &z__2);
#line 601 "zlaqr4.f"
			    rtdisc.r = z__1.r, rtdisc.i = z__1.i;
#line 602 "zlaqr4.f"
			    i__2 = kbot - 1;
#line 602 "zlaqr4.f"
			    z__2.r = tr2.r + rtdisc.r, z__2.i = tr2.i + 
				    rtdisc.i;
#line 602 "zlaqr4.f"
			    z__1.r = s * z__2.r, z__1.i = s * z__2.i;
#line 602 "zlaqr4.f"
			    w[i__2].r = z__1.r, w[i__2].i = z__1.i;
#line 603 "zlaqr4.f"
			    i__2 = kbot;
#line 603 "zlaqr4.f"
			    z__2.r = tr2.r - rtdisc.r, z__2.i = tr2.i - 
				    rtdisc.i;
#line 603 "zlaqr4.f"
			    z__1.r = s * z__2.r, z__1.i = s * z__2.i;
#line 603 "zlaqr4.f"
			    w[i__2].r = z__1.r, w[i__2].i = z__1.i;

#line 605 "zlaqr4.f"
			    ks = kbot - 1;
#line 606 "zlaqr4.f"
			}
#line 607 "zlaqr4.f"
		    }

#line 609 "zlaqr4.f"
		    if (kbot - ks + 1 > ns) {

/*                    ==== Sort the shifts (Helps a little) ==== */

#line 613 "zlaqr4.f"
			sorted = FALSE_;
#line 614 "zlaqr4.f"
			i__2 = ks + 1;
#line 614 "zlaqr4.f"
			for (k = kbot; k >= i__2; --k) {
#line 615 "zlaqr4.f"
			    if (sorted) {
#line 615 "zlaqr4.f"
				goto L60;
#line 615 "zlaqr4.f"
			    }
#line 617 "zlaqr4.f"
			    sorted = TRUE_;
#line 618 "zlaqr4.f"
			    i__3 = k - 1;
#line 618 "zlaqr4.f"
			    for (i__ = ks; i__ <= i__3; ++i__) {
#line 619 "zlaqr4.f"
				i__4 = i__;
#line 619 "zlaqr4.f"
				i__5 = i__ + 1;
#line 619 "zlaqr4.f"
				if ((d__1 = w[i__4].r, abs(d__1)) + (d__2 = 
					d_imag(&w[i__]), abs(d__2)) < (d__3 = 
					w[i__5].r, abs(d__3)) + (d__4 = 
					d_imag(&w[i__ + 1]), abs(d__4))) {
#line 621 "zlaqr4.f"
				    sorted = FALSE_;
#line 622 "zlaqr4.f"
				    i__4 = i__;
#line 622 "zlaqr4.f"
				    swap.r = w[i__4].r, swap.i = w[i__4].i;
#line 623 "zlaqr4.f"
				    i__4 = i__;
#line 623 "zlaqr4.f"
				    i__5 = i__ + 1;
#line 623 "zlaqr4.f"
				    w[i__4].r = w[i__5].r, w[i__4].i = w[i__5]
					    .i;
#line 624 "zlaqr4.f"
				    i__4 = i__ + 1;
#line 624 "zlaqr4.f"
				    w[i__4].r = swap.r, w[i__4].i = swap.i;
#line 625 "zlaqr4.f"
				}
#line 626 "zlaqr4.f"
/* L40: */
#line 626 "zlaqr4.f"
			    }
#line 627 "zlaqr4.f"
/* L50: */
#line 627 "zlaqr4.f"
			}
#line 628 "zlaqr4.f"
L60:
#line 629 "zlaqr4.f"
			;
#line 629 "zlaqr4.f"
		    }
#line 630 "zlaqr4.f"
		}

/*              ==== If there are only two shifts, then use */
/*              .    only one.  ==== */

#line 635 "zlaqr4.f"
		if (kbot - ks + 1 == 2) {
#line 636 "zlaqr4.f"
		    i__2 = kbot;
#line 636 "zlaqr4.f"
		    i__3 = kbot + kbot * h_dim1;
#line 636 "zlaqr4.f"
		    z__2.r = w[i__2].r - h__[i__3].r, z__2.i = w[i__2].i - 
			    h__[i__3].i;
#line 636 "zlaqr4.f"
		    z__1.r = z__2.r, z__1.i = z__2.i;
#line 636 "zlaqr4.f"
		    i__4 = kbot - 1;
#line 636 "zlaqr4.f"
		    i__5 = kbot + kbot * h_dim1;
#line 636 "zlaqr4.f"
		    z__4.r = w[i__4].r - h__[i__5].r, z__4.i = w[i__4].i - 
			    h__[i__5].i;
#line 636 "zlaqr4.f"
		    z__3.r = z__4.r, z__3.i = z__4.i;
#line 636 "zlaqr4.f"
		    if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			    abs(d__2)) < (d__3 = z__3.r, abs(d__3)) + (d__4 = 
			    d_imag(&z__3), abs(d__4))) {
#line 638 "zlaqr4.f"
			i__2 = kbot - 1;
#line 638 "zlaqr4.f"
			i__3 = kbot;
#line 638 "zlaqr4.f"
			w[i__2].r = w[i__3].r, w[i__2].i = w[i__3].i;
#line 639 "zlaqr4.f"
		    } else {
#line 640 "zlaqr4.f"
			i__2 = kbot;
#line 640 "zlaqr4.f"
			i__3 = kbot - 1;
#line 640 "zlaqr4.f"
			w[i__2].r = w[i__3].r, w[i__2].i = w[i__3].i;
#line 641 "zlaqr4.f"
		    }
#line 642 "zlaqr4.f"
		}

/*              ==== Use up to NS of the the smallest magnatiude */
/*              .    shifts.  If there aren't NS shifts available, */
/*              .    then use them all, possibly dropping one to */
/*              .    make the number of shifts even. ==== */

/* Computing MIN */
#line 649 "zlaqr4.f"
		i__2 = ns, i__3 = kbot - ks + 1;
#line 649 "zlaqr4.f"
		ns = min(i__2,i__3);
#line 650 "zlaqr4.f"
		ns -= ns % 2;
#line 651 "zlaqr4.f"
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

#line 664 "zlaqr4.f"
		kdu = ns * 3 - 3;
#line 665 "zlaqr4.f"
		ku = *n - kdu + 1;
#line 666 "zlaqr4.f"
		kwh = kdu + 1;
#line 667 "zlaqr4.f"
		nho = *n - kdu - 3 - (kdu + 1) + 1;
#line 668 "zlaqr4.f"
		kwv = kdu + 4;
#line 669 "zlaqr4.f"
		nve = *n - kdu - kwv + 1;

/*              ==== Small-bulge multi-shift QR sweep ==== */

#line 673 "zlaqr4.f"
		zlaqr5_(wantt, wantz, &kacc22, n, &ktop, &kbot, &ns, &w[ks], &
			h__[h_offset], ldh, iloz, ihiz, &z__[z_offset], ldz, &
			work[1], &c__3, &h__[ku + h_dim1], ldh, &nve, &h__[
			kwv + h_dim1], ldh, &nho, &h__[ku + kwh * h_dim1], 
			ldh);
#line 677 "zlaqr4.f"
	    }

/*           ==== Note progress (or the lack of it). ==== */

#line 681 "zlaqr4.f"
	    if (ld > 0) {
#line 682 "zlaqr4.f"
		ndfl = 1;
#line 683 "zlaqr4.f"
	    } else {
#line 684 "zlaqr4.f"
		++ndfl;
#line 685 "zlaqr4.f"
	    }

/*           ==== End of main loop ==== */
#line 688 "zlaqr4.f"
/* L70: */
#line 688 "zlaqr4.f"
	}

/*        ==== Iteration limit exceeded.  Set INFO to show where */
/*        .    the problem occurred and exit. ==== */

#line 693 "zlaqr4.f"
	*info = kbot;
#line 694 "zlaqr4.f"
L80:
#line 695 "zlaqr4.f"
	;
#line 695 "zlaqr4.f"
    }

/*     ==== Return the optimal value of LWORK. ==== */

#line 699 "zlaqr4.f"
    d__1 = (doublereal) lwkopt;
#line 699 "zlaqr4.f"
    z__1.r = d__1, z__1.i = 0.;
#line 699 "zlaqr4.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

/*     ==== End of ZLAQR4 ==== */

#line 703 "zlaqr4.f"
    return 0;
} /* zlaqr4_ */

