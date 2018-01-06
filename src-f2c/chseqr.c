#line 1 "chseqr.f"
/* chseqr.f -- translated by f2c (version 20100827).
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

#line 1 "chseqr.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__12 = 12;
static integer c__2 = 2;
static integer c__49 = 49;

/* > \brief \b CHSEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHSEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chseqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chseqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chseqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N */
/*       CHARACTER          COMPZ, JOB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CHSEQR computes the eigenvalues of a Hessenberg matrix H */
/* >    and, optionally, the matrices T and Z from the Schur decomposition */
/* >    H = Z T Z**H, where T is an upper triangular matrix (the */
/* >    Schur form), and Z is the unitary matrix of Schur vectors. */
/* > */
/* >    Optionally Z may be postmultiplied into an input unitary */
/* >    matrix Q so that this routine can give the Schur factorization */
/* >    of a matrix A which has been reduced to the Hessenberg form H */
/* >    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >           = 'E':  compute eigenvalues only; */
/* >           = 'S':  compute eigenvalues and the Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >           = 'N':  no Schur vectors are computed; */
/* >           = 'I':  Z is initialized to the unit matrix and the matrix Z */
/* >                   of Schur vectors of H is returned; */
/* >           = 'V':  Z must contain an unitary matrix Q on entry, and */
/* >                   the product Q*Z is returned. */
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
/* > */
/* >           It is assumed that H is already upper triangular in rows */
/* >           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/* >           set by a previous call to CGEBAL, and then passed to ZGEHRD */
/* >           when the matrix output by CGEBAL is reduced to Hessenberg */
/* >           form. Otherwise ILO and IHI should be set to 1 and N */
/* >           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N. */
/* >           If N = 0, then ILO = 1 and IHI = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX array, dimension (LDH,N) */
/* >           On entry, the upper Hessenberg matrix H. */
/* >           On exit, if INFO = 0 and JOB = 'S', H contains the upper */
/* >           triangular matrix T from the Schur decomposition (the */
/* >           Schur form). If INFO = 0 and JOB = 'E', the contents of */
/* >           H are unspecified on exit.  (The output value of H when */
/* >           INFO.GT.0 is given under the description of INFO below.) */
/* > */
/* >           Unlike earlier versions of CHSEQR, this subroutine may */
/* >           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1 */
/* >           or j = IHI+1, IHI+2, ... N. */
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
/* >          W is COMPLEX array, dimension (N) */
/* >           The computed eigenvalues. If JOB = 'S', the eigenvalues are */
/* >           stored in the same order as on the diagonal of the Schur */
/* >           form returned in H, with W(i) = H(i,i). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ,N) */
/* >           If COMPZ = 'N', Z is not referenced. */
/* >           If COMPZ = 'I', on entry Z need not be set and on exit, */
/* >           if INFO = 0, Z contains the unitary matrix Z of the Schur */
/* >           vectors of H.  If COMPZ = 'V', on entry Z must contain an */
/* >           N-by-N matrix Q, which is assumed to be equal to the unit */
/* >           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit, */
/* >           if INFO = 0, Z contains Q*Z. */
/* >           Normally Q is the unitary matrix generated by CUNGHR */
/* >           after the call to CGEHRD which formed the Hessenberg matrix */
/* >           H. (The output value of Z when INFO.GT.0 is given under */
/* >           the description of INFO below.) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >           The leading dimension of the array Z.  if COMPZ = 'I' or */
/* >           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (LWORK) */
/* >           On exit, if INFO = 0, WORK(1) returns an estimate of */
/* >           the optimal value for LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >           The dimension of the array WORK.  LWORK .GE. max(1,N) */
/* >           is sufficient and delivers very good and sometimes */
/* >           optimal performance.  However, LWORK as large as 11*N */
/* >           may be required for optimal performance.  A workspace */
/* >           query is recommended to determine the optimal workspace */
/* >           size. */
/* > */
/* >           If LWORK = -1, then CHSEQR does a workspace query. */
/* >           In this case, CHSEQR checks the input parameters and */
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
/* >           .LT. 0:  if INFO = -i, the i-th argument had an illegal */
/* >                    value */
/* >           .GT. 0:  if INFO = i, CHSEQR failed to compute all of */
/* >                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR */
/* >                and WI contain those eigenvalues which have been */
/* >                successfully computed.  (Failures are rare.) */
/* > */
/* >                If INFO .GT. 0 and JOB = 'E', then on exit, the */
/* >                remaining unconverged eigenvalues are the eigen- */
/* >                values of the upper Hessenberg matrix rows and */
/* >                columns ILO through INFO of the final, output */
/* >                value of H. */
/* > */
/* >                If INFO .GT. 0 and JOB   = 'S', then on exit */
/* > */
/* >           (*)  (initial value of H)*U  = U*(final value of H) */
/* > */
/* >                where U is a unitary matrix.  The final */
/* >                value of  H is upper Hessenberg and triangular in */
/* >                rows and columns INFO+1 through IHI. */
/* > */
/* >                If INFO .GT. 0 and COMPZ = 'V', then on exit */
/* > */
/* >                  (final value of Z)  =  (initial value of Z)*U */
/* > */
/* >                where U is the unitary matrix in (*) (regard- */
/* >                less of the value of JOB.) */
/* > */
/* >                If INFO .GT. 0 and COMPZ = 'I', then on exit */
/* >                      (final value of Z)  = U */
/* >                where U is the unitary matrix in (*) (regard- */
/* >                less of the value of JOB.) */
/* > */
/* >                If INFO .GT. 0 and COMPZ = 'N', then Z is not */
/* >                accessed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >             Default values supplied by */
/* >             ILAENV(ISPEC,'CHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK). */
/* >             It is suggested that these defaults be adjusted in order */
/* >             to attain best performance in each particular */
/* >             computational environment. */
/* > */
/* >            ISPEC=12: The CLAHQR vs CLAQR0 crossover point. */
/* >                      Default: 75. (Must be at least 11.) */
/* > */
/* >            ISPEC=13: Recommended deflation window size. */
/* >                      This depends on ILO, IHI and NS.  NS is the */
/* >                      number of simultaneous shifts returned */
/* >                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.) */
/* >                      The default for (IHI-ILO+1).LE.500 is NS. */
/* >                      The default for (IHI-ILO+1).GT.500 is 3*NS/2. */
/* > */
/* >            ISPEC=14: Nibble crossover point. (See IPARMQ for */
/* >                      details.)  Default: 14% of deflation window */
/* >                      size. */
/* > */
/* >            ISPEC=15: Number of simultaneous shifts in a multishift */
/* >                      QR iteration. */
/* > */
/* >                      If IHI-ILO+1 is ... */
/* > */
/* >                      greater than      ...but less    ... the */
/* >                      or equal to ...      than        default is */
/* > */
/* >                           1               30          NS =   2(+) */
/* >                          30               60          NS =   4(+) */
/* >                          60              150          NS =  10(+) */
/* >                         150              590          NS =  ** */
/* >                         590             3000          NS =  64 */
/* >                        3000             6000          NS = 128 */
/* >                        6000             infinity      NS = 256 */
/* > */
/* >                  (+)  By default some or all matrices of this order */
/* >                       are passed to the implicit double shift routine */
/* >                       CLAHQR and this parameter is ignored.  See */
/* >                       ISPEC=12 above and comments in IPARMQ for */
/* >                       details. */
/* > */
/* >                 (**)  The asterisks (**) indicate an ad-hoc */
/* >                       function of N increasing from 10 to 64. */
/* > */
/* >            ISPEC=16: Select structured matrix multiply. */
/* >                      If the number of simultaneous shifts (specified */
/* >                      by ISPEC=15) is less than 14, then the default */
/* >                      for ISPEC=16 is 0.  Otherwise the default for */
/* >                      ISPEC=16 is 2. */
/* > \endverbatim */

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

/*  ===================================================================== */
/* Subroutine */ int chseqr_(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, doublecomplex *h__, integer *ldh, doublecomplex *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork,
	 integer *info, ftnlen job_len, ftnlen compz_len)
{
    /* System generated locals */
    address a__1[2];
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3[2];
    doublereal d__1, d__2, d__3;
    doublecomplex z__1;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static doublecomplex hl[2401]	/* was [49][49] */;
    static integer kbot, nmin;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical initz;
    static doublecomplex workl[49];
    static logical wantt, wantz;
    extern /* Subroutine */ int claqr0_(logical *, logical *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, integer *), clahqr_(logical *, logical *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, integer *, integer *), 
	    clacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), claset_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

/*     ==== Matrices of order NTINY or smaller must be processed by */
/*     .    CLAHQR because of insufficient subdiagonal scratch space. */
/*     .    (This is a hard limit.) ==== */

/*     ==== NL allocates some local workspace to help small matrices */
/*     .    through a rare CLAHQR failure.  NL .GT. NTINY = 11 is */
/*     .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom- */
/*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49 */
/*     .    allows up to six simultaneous shifts and a 16-by-16 */
/*     .    deflation window.  ==== */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== Decode and check the input parameters. ==== */

#line 361 "chseqr.f"
    /* Parameter adjustments */
#line 361 "chseqr.f"
    h_dim1 = *ldh;
#line 361 "chseqr.f"
    h_offset = 1 + h_dim1;
#line 361 "chseqr.f"
    h__ -= h_offset;
#line 361 "chseqr.f"
    --w;
#line 361 "chseqr.f"
    z_dim1 = *ldz;
#line 361 "chseqr.f"
    z_offset = 1 + z_dim1;
#line 361 "chseqr.f"
    z__ -= z_offset;
#line 361 "chseqr.f"
    --work;
#line 361 "chseqr.f"

#line 361 "chseqr.f"
    /* Function Body */
#line 361 "chseqr.f"
    wantt = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 362 "chseqr.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 363 "chseqr.f"
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);
#line 364 "chseqr.f"
    d__1 = (doublereal) max(1,*n);
#line 364 "chseqr.f"
    z__1.r = d__1, z__1.i = 0.;
#line 364 "chseqr.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 365 "chseqr.f"
    lquery = *lwork == -1;

#line 367 "chseqr.f"
    *info = 0;
#line 368 "chseqr.f"
    if (! lsame_(job, "E", (ftnlen)1, (ftnlen)1) && ! wantt) {
#line 369 "chseqr.f"
	*info = -1;
#line 370 "chseqr.f"
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
#line 371 "chseqr.f"
	*info = -2;
#line 372 "chseqr.f"
    } else if (*n < 0) {
#line 373 "chseqr.f"
	*info = -3;
#line 374 "chseqr.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 375 "chseqr.f"
	*info = -4;
#line 376 "chseqr.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 377 "chseqr.f"
	*info = -5;
#line 378 "chseqr.f"
    } else if (*ldh < max(1,*n)) {
#line 379 "chseqr.f"
	*info = -7;
#line 380 "chseqr.f"
    } else if (*ldz < 1 || wantz && *ldz < max(1,*n)) {
#line 381 "chseqr.f"
	*info = -10;
#line 382 "chseqr.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 383 "chseqr.f"
	*info = -12;
#line 384 "chseqr.f"
    }

#line 386 "chseqr.f"
    if (*info != 0) {

/*        ==== Quick return in case of invalid argument. ==== */

#line 390 "chseqr.f"
	i__1 = -(*info);
#line 390 "chseqr.f"
	xerbla_("CHSEQR", &i__1, (ftnlen)6);
#line 391 "chseqr.f"
	return 0;

#line 393 "chseqr.f"
    } else if (*n == 0) {

/*        ==== Quick return in case N = 0; nothing to do. ==== */

#line 397 "chseqr.f"
	return 0;

#line 399 "chseqr.f"
    } else if (lquery) {

/*        ==== Quick return in case of a workspace query ==== */

#line 403 "chseqr.f"
	claqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], ilo, 
		ihi, &z__[z_offset], ldz, &work[1], lwork, info);
/*        ==== Ensure reported workspace size is backward-compatible with */
/*        .    previous LAPACK versions. ==== */
/* Computing MAX */
#line 407 "chseqr.f"
	d__2 = work[1].r, d__3 = (doublereal) max(1,*n);
#line 407 "chseqr.f"
	d__1 = max(d__2,d__3);
#line 407 "chseqr.f"
	z__1.r = d__1, z__1.i = 0.;
#line 407 "chseqr.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 409 "chseqr.f"
	return 0;

#line 411 "chseqr.f"
    } else {

/*        ==== copy eigenvalues isolated by CGEBAL ==== */

#line 415 "chseqr.f"
	if (*ilo > 1) {
#line 415 "chseqr.f"
	    i__1 = *ilo - 1;
#line 415 "chseqr.f"
	    i__2 = *ldh + 1;
#line 415 "chseqr.f"
	    ccopy_(&i__1, &h__[h_offset], &i__2, &w[1], &c__1);
#line 415 "chseqr.f"
	}
#line 417 "chseqr.f"
	if (*ihi < *n) {
#line 417 "chseqr.f"
	    i__1 = *n - *ihi;
#line 417 "chseqr.f"
	    i__2 = *ldh + 1;
#line 417 "chseqr.f"
	    ccopy_(&i__1, &h__[*ihi + 1 + (*ihi + 1) * h_dim1], &i__2, &w[*
		    ihi + 1], &c__1);
#line 417 "chseqr.f"
	}

/*        ==== Initialize Z, if requested ==== */

#line 422 "chseqr.f"
	if (initz) {
#line 422 "chseqr.f"
	    claset_("A", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)1);
#line 422 "chseqr.f"
	}

/*        ==== Quick return if possible ==== */

#line 427 "chseqr.f"
	if (*ilo == *ihi) {
#line 428 "chseqr.f"
	    i__1 = *ilo;
#line 428 "chseqr.f"
	    i__2 = *ilo + *ilo * h_dim1;
#line 428 "chseqr.f"
	    w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
#line 429 "chseqr.f"
	    return 0;
#line 430 "chseqr.f"
	}

/*        ==== CLAHQR/CLAQR0 crossover point ==== */

/* Writing concatenation */
#line 434 "chseqr.f"
	i__3[0] = 1, a__1[0] = job;
#line 434 "chseqr.f"
	i__3[1] = 1, a__1[1] = compz;
#line 434 "chseqr.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 434 "chseqr.f"
	nmin = ilaenv_(&c__12, "CHSEQR", ch__1, n, ilo, ihi, lwork, (ftnlen)6,
		 (ftnlen)2);
#line 436 "chseqr.f"
	nmin = max(11,nmin);

/*        ==== CLAQR0 for big matrices; CLAHQR for small ones ==== */

#line 440 "chseqr.f"
	if (*n > nmin) {
#line 441 "chseqr.f"
	    claqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], 
		    ilo, ihi, &z__[z_offset], ldz, &work[1], lwork, info);
#line 443 "chseqr.f"
	} else {

/*           ==== Small matrix ==== */

#line 447 "chseqr.f"
	    clahqr_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], 
		    ilo, ihi, &z__[z_offset], ldz, info);

#line 450 "chseqr.f"
	    if (*info > 0) {

/*              ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds */
/*              .    when CLAHQR fails. ==== */

#line 455 "chseqr.f"
		kbot = *info;

#line 457 "chseqr.f"
		if (*n >= 49) {

/*                 ==== Larger matrices have enough subdiagonal scratch */
/*                 .    space to call CLAQR0 directly. ==== */

#line 462 "chseqr.f"
		    claqr0_(&wantt, &wantz, n, ilo, &kbot, &h__[h_offset], 
			    ldh, &w[1], ilo, ihi, &z__[z_offset], ldz, &work[
			    1], lwork, info);

#line 465 "chseqr.f"
		} else {

/*                 ==== Tiny matrices don't have enough subdiagonal */
/*                 .    scratch space to benefit from CLAQR0.  Hence, */
/*                 .    tiny matrices must be copied into a larger */
/*                 .    array before calling CLAQR0. ==== */

#line 472 "chseqr.f"
		    clacpy_("A", n, n, &h__[h_offset], ldh, hl, &c__49, (
			    ftnlen)1);
#line 473 "chseqr.f"
		    i__1 = *n + 1 + *n * 49 - 50;
#line 473 "chseqr.f"
		    hl[i__1].r = 0., hl[i__1].i = 0.;
#line 474 "chseqr.f"
		    i__1 = 49 - *n;
#line 474 "chseqr.f"
		    claset_("A", &c__49, &i__1, &c_b1, &c_b1, &hl[(*n + 1) * 
			    49 - 49], &c__49, (ftnlen)1);
#line 476 "chseqr.f"
		    claqr0_(&wantt, &wantz, &c__49, ilo, &kbot, hl, &c__49, &
			    w[1], ilo, ihi, &z__[z_offset], ldz, workl, &
			    c__49, info);
#line 478 "chseqr.f"
		    if (wantt || *info != 0) {
#line 478 "chseqr.f"
			clacpy_("A", n, n, hl, &c__49, &h__[h_offset], ldh, (
				ftnlen)1);
#line 478 "chseqr.f"
		    }
#line 480 "chseqr.f"
		}
#line 481 "chseqr.f"
	    }
#line 482 "chseqr.f"
	}

/*        ==== Clear out the trash, if necessary. ==== */

#line 486 "chseqr.f"
	if ((wantt || *info != 0) && *n > 2) {
#line 486 "chseqr.f"
	    i__1 = *n - 2;
#line 486 "chseqr.f"
	    i__2 = *n - 2;
#line 486 "chseqr.f"
	    claset_("L", &i__1, &i__2, &c_b1, &c_b1, &h__[h_dim1 + 3], ldh, (
		    ftnlen)1);
#line 486 "chseqr.f"
	}

/*        ==== Ensure reported workspace size is backward-compatible with */
/*        .    previous LAPACK versions. ==== */

/* Computing MAX */
#line 492 "chseqr.f"
	d__2 = (doublereal) max(1,*n), d__3 = work[1].r;
#line 492 "chseqr.f"
	d__1 = max(d__2,d__3);
#line 492 "chseqr.f"
	z__1.r = d__1, z__1.i = 0.;
#line 492 "chseqr.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 494 "chseqr.f"
    }

/*     ==== End of CHSEQR ==== */

#line 498 "chseqr.f"
    return 0;
} /* chseqr_ */

