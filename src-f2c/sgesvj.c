#line 1 "sgesvj.f"
/* sgesvj.f -- translated by f2c (version 20100827).
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

#line 1 "sgesvj.f"
/* Table of constant values */

static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/* > \brief \b SGESVJ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESVJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvj.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvj.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvj.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, */
/*                          LDV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N */
/*       CHARACTER*1        JOBA, JOBU, JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), SVA( N ), V( LDV, * ), */
/*      $                   WORK( LWORK ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGESVJ computes the singular value decomposition (SVD) of a real */
/* > M-by-N matrix A, where M >= N. The SVD of A is written as */
/* >                                    [++]   [xx]   [x0]   [xx] */
/* >              A = U * SIGMA * V^t,  [++] = [xx] * [ox] * [xx] */
/* >                                    [++]   [xx] */
/* > where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal */
/* > matrix, and V is an N-by-N orthogonal matrix. The diagonal elements */
/* > of SIGMA are the singular values of A. The columns of U and V are the */
/* > left and the right singular vectors of A, respectively. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBA */
/* > \verbatim */
/* >          JOBA is CHARACTER* 1 */
/* >          Specifies the structure of A. */
/* >          = 'L': The input matrix A is lower triangular; */
/* >          = 'U': The input matrix A is upper triangular; */
/* >          = 'G': The input matrix A is general M-by-N matrix, M >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          Specifies whether to compute the left singular vectors */
/* >          (columns of U): */
/* >          = 'U': The left singular vectors corresponding to the nonzero */
/* >                 singular values are computed and returned in the leading */
/* >                 columns of A. See more details in the description of A. */
/* >                 The default numerical orthogonality threshold is set to */
/* >                 approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=SLAMCH('E'). */
/* >          = 'C': Analogous to JOBU='U', except that user can control the */
/* >                 level of numerical orthogonality of the computed left */
/* >                 singular vectors. TOL can be set to TOL = CTOL*EPS, where */
/* >                 CTOL is given on input in the array WORK. */
/* >                 No CTOL smaller than ONE is allowed. CTOL greater */
/* >                 than 1 / EPS is meaningless. The option 'C' */
/* >                 can be used if M*EPS is satisfactory orthogonality */
/* >                 of the computed left singular vectors, so CTOL=M could */
/* >                 save few sweeps of Jacobi rotations. */
/* >                 See the descriptions of A and WORK(1). */
/* >          = 'N': The matrix U is not computed. However, see the */
/* >                 description of A. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >          Specifies whether to compute the right singular vectors, that */
/* >          is, the matrix V: */
/* >          = 'V' : the matrix V is computed and returned in the array V */
/* >          = 'A' : the Jacobi rotations are applied to the MV-by-N */
/* >                  array V. In other words, the right singular vector */
/* >                  matrix V is not computed explicitly; instead it is */
/* >                  applied to an MV-by-N matrix initially stored in the */
/* >                  first MV rows of V. */
/* >          = 'N' : the matrix V is not computed and the array V is not */
/* >                  referenced */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the input matrix A. 1/SLAMCH('E') > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the input matrix A. */
/* >          M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          If JOBU .EQ. 'U' .OR. JOBU .EQ. 'C': */
/* >                 If INFO .EQ. 0 : */
/* >                 RANKA orthonormal columns of U are returned in the */
/* >                 leading RANKA columns of the array A. Here RANKA <= N */
/* >                 is the number of computed singular values of A that are */
/* >                 above the underflow threshold SLAMCH('S'). The singular */
/* >                 vectors corresponding to underflowed or zero singular */
/* >                 values are not computed. The value of RANKA is returned */
/* >                 in the array WORK as RANKA=NINT(WORK(2)). Also see the */
/* >                 descriptions of SVA and WORK. The computed columns of U */
/* >                 are mutually numerically orthogonal up to approximately */
/* >                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU.EQ.'C'), */
/* >                 see the description of JOBU. */
/* >                 If INFO .GT. 0, */
/* >                 the procedure SGESVJ did not converge in the given number */
/* >                 of iterations (sweeps). In that case, the computed */
/* >                 columns of U may not be orthogonal up to TOL. The output */
/* >                 U (stored in A), SIGMA (given by the computed singular */
/* >                 values in SVA(1:N)) and V is still a decomposition of the */
/* >                 input matrix A in the sense that the residual */
/* >                 ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small. */
/* >          If JOBU .EQ. 'N': */
/* >                 If INFO .EQ. 0 : */
/* >                 Note that the left singular vectors are 'for free' in the */
/* >                 one-sided Jacobi SVD algorithm. However, if only the */
/* >                 singular values are needed, the level of numerical */
/* >                 orthogonality of U is not an issue and iterations are */
/* >                 stopped when the columns of the iterated matrix are */
/* >                 numerically orthogonal up to approximately M*EPS. Thus, */
/* >                 on exit, A contains the columns of U scaled with the */
/* >                 corresponding singular values. */
/* >                 If INFO .GT. 0 : */
/* >                 the procedure SGESVJ did not converge in the given number */
/* >                 of iterations (sweeps). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] SVA */
/* > \verbatim */
/* >          SVA is REAL array, dimension (N) */
/* >          On exit, */
/* >          If INFO .EQ. 0 : */
/* >          depending on the value SCALE = WORK(1), we have: */
/* >                 If SCALE .EQ. ONE: */
/* >                 SVA(1:N) contains the computed singular values of A. */
/* >                 During the computation SVA contains the Euclidean column */
/* >                 norms of the iterated matrices in the array A. */
/* >                 If SCALE .NE. ONE: */
/* >                 The singular values of A are SCALE*SVA(1:N), and this */
/* >                 factored representation is due to the fact that some of the */
/* >                 singular values of A might underflow or overflow. */
/* > */
/* >          If INFO .GT. 0 : */
/* >          the procedure SGESVJ did not converge in the given number of */
/* >          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate. */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* >          MV is INTEGER */
/* >          If JOBV .EQ. 'A', then the product of Jacobi rotations in SGESVJ */
/* >          is applied to the first MV rows of V. See the description of JOBV. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is REAL array, dimension (LDV,N) */
/* >          If JOBV = 'V', then V contains on exit the N-by-N matrix of */
/* >                         the right singular vectors; */
/* >          If JOBV = 'A', then V contains the product of the computed right */
/* >                         singular vector matrix and the initial matrix in */
/* >                         the array V. */
/* >          If JOBV = 'N', then V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V, LDV .GE. 1. */
/* >          If JOBV .EQ. 'V', then LDV .GE. max(1,N). */
/* >          If JOBV .EQ. 'A', then LDV .GE. max(1,MV) . */
/* > \endverbatim */
/* > */
/* > \param[in,out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension max(4,M+N). */
/* >          On entry, */
/* >          If JOBU .EQ. 'C' : */
/* >          WORK(1) = CTOL, where CTOL defines the threshold for convergence. */
/* >                    The process stops if all columns of A are mutually */
/* >                    orthogonal up to CTOL*EPS, EPS=SLAMCH('E'). */
/* >                    It is required that CTOL >= ONE, i.e. it is not */
/* >                    allowed to force the routine to obtain orthogonality */
/* >                    below EPSILON. */
/* >          On exit, */
/* >          WORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N) */
/* >                    are the computed singular vcalues of A. */
/* >                    (See description of SVA().) */
/* >          WORK(2) = NINT(WORK(2)) is the number of the computed nonzero */
/* >                    singular values. */
/* >          WORK(3) = NINT(WORK(3)) is the number of the computed singular */
/* >                    values that are larger than the underflow threshold. */
/* >          WORK(4) = NINT(WORK(4)) is the number of sweeps of Jacobi */
/* >                    rotations needed for numerical convergence. */
/* >          WORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep. */
/* >                    This is useful information in cases when SGESVJ did */
/* >                    not converge, as it can be used to estimate whether */
/* >                    the output is stil useful and for post festum analysis. */
/* >          WORK(6) = the largest absolute value over all sines of the */
/* >                    Jacobi rotation angles in the last sweep. It can be */
/* >                    useful for a post festum analysis. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >         length of WORK, WORK >= MAX(6,M+N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0 : successful exit. */
/* >          < 0 : if INFO = -i, then the i-th argument had an illegal value */
/* >          > 0 : SGESVJ did not converge in the maximal allowed number (30) */
/* >                of sweeps. The output may still be useful. See the */
/* >                description of WORK. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane */
/* > rotations. The rotations are implemented as fast scaled rotations of */
/* > Anda and Park [1]. In the case of underflow of the Jacobi angle, a */
/* > modified Jacobi transformation of Drmac [4] is used. Pivot strategy uses */
/* > column interchanges of de Rijk [2]. The relative accuracy of the computed */
/* > singular values and the accuracy of the computed singular vectors (in */
/* > angle metric) is as guaranteed by the theory of Demmel and Veselic [3]. */
/* > The condition number that determines the accuracy in the full rank case */
/* > is essentially min_{D=diag} kappa(A*D), where kappa(.) is the */
/* > spectral condition number. The best performance of this Jacobi SVD */
/* > procedure is achieved if used in an  accelerated version of Drmac and */
/* > Veselic [5,6], and it is the kernel routine in the SIGMA library [7]. */
/* > Some tunning parameters (marked with [TP]) are available for the */
/* > implementer. \n */
/* > The computational range for the nonzero singular values is the  machine */
/* > number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even */
/* > denormalized singular values can be computed with the corresponding */
/* > gradual loss of accurate digits. */
/* > */
/* > \par Contributors: */
/*  ================== */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */
/* > */
/* > \par References: */
/*  ================ */
/* > */
/* > [1] A. A. Anda and H. Park: Fast plane rotations with dynamic scaling. \n */
/* >    SIAM J. matrix Anal. Appl., Vol. 15 (1994), pp. 162-174. \n\n */
/* > [2] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the */
/* >    singular value decomposition on a vector computer. \n */
/* >    SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371. \n\n */
/* > [3] J. Demmel and K. Veselic: Jacobi method is more accurate than QR. \n */
/* > [4] Z. Drmac: Implementation of Jacobi rotations for accurate singular */
/* >    value computation in floating point arithmetic. \n */
/* >    SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222. \n\n */
/* > [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. \n */
/* >    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. \n */
/* >    LAPACK Working note 169. \n\n */
/* > [6] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. \n */
/* >    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. \n */
/* >    LAPACK Working note 170. \n\n */
/* > [7] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV, */
/* >    QSVD, (H,K)-SVD computations.\n */
/* >    Department of Mathematics, University of Zagreb, 2008. */
/* > */
/* > \par Bugs, Examples and Comments: */
/*  ================================= */
/* > */
/* > Please report all bugs and send interesting test examples and comments to */
/* > drmac@math.hr. Thank you. */

/*  ===================================================================== */
/* Subroutine */ int sgesvj_(char *joba, char *jobu, char *jobv, integer *m, 
	integer *n, doublereal *a, integer *lda, doublereal *sva, integer *mv,
	 doublereal *v, integer *ldv, doublereal *work, integer *lwork, 
	integer *info, ftnlen joba_len, ftnlen jobu_len, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal bigtheta;
    static integer pskipped, i__, p, q;
    static doublereal t;
    static integer n2, n4;
    static doublereal rootsfmin;
    static integer n34;
    static doublereal cs, sn;
    static integer ir1, jbc;
    static doublereal big;
    static integer kbl, igl, ibr, jgl, nbl;
    static doublereal skl, tol;
    static integer mvl;
    static doublereal aapp, aapq, aaqq, ctol;
    static integer ierr;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal aapp0, temp1;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal large, apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal small, sfmin;
    static logical lsvec;
    static doublereal fastr[5], epsln;
    static logical applv, rsvec, uctol, lower, upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical rotok;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), srotm_(integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *)
	    , sgsvj0_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen), sgsvj1_(char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static integer blskip;
    static doublereal mxaapq;
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal thsign;
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal mxsinj;
    static integer emptsw, notrot, iswrot, lkahead;
    static logical goscale, noscale;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     from BLAS */
/*     from LAPACK */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     from BLAS */
/*     from LAPACK */

/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 391 "sgesvj.f"
    /* Parameter adjustments */
#line 391 "sgesvj.f"
    --sva;
#line 391 "sgesvj.f"
    a_dim1 = *lda;
#line 391 "sgesvj.f"
    a_offset = 1 + a_dim1;
#line 391 "sgesvj.f"
    a -= a_offset;
#line 391 "sgesvj.f"
    v_dim1 = *ldv;
#line 391 "sgesvj.f"
    v_offset = 1 + v_dim1;
#line 391 "sgesvj.f"
    v -= v_offset;
#line 391 "sgesvj.f"
    --work;
#line 391 "sgesvj.f"

#line 391 "sgesvj.f"
    /* Function Body */
#line 391 "sgesvj.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 392 "sgesvj.f"
    uctol = lsame_(jobu, "C", (ftnlen)1, (ftnlen)1);
#line 393 "sgesvj.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 394 "sgesvj.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 395 "sgesvj.f"
    upper = lsame_(joba, "U", (ftnlen)1, (ftnlen)1);
#line 396 "sgesvj.f"
    lower = lsame_(joba, "L", (ftnlen)1, (ftnlen)1);

#line 398 "sgesvj.f"
    if (! (upper || lower || lsame_(joba, "G", (ftnlen)1, (ftnlen)1))) {
#line 399 "sgesvj.f"
	*info = -1;
#line 400 "sgesvj.f"
    } else if (! (lsvec || uctol || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 401 "sgesvj.f"
	*info = -2;
#line 402 "sgesvj.f"
    } else if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 403 "sgesvj.f"
	*info = -3;
#line 404 "sgesvj.f"
    } else if (*m < 0) {
#line 405 "sgesvj.f"
	*info = -4;
#line 406 "sgesvj.f"
    } else if (*n < 0 || *n > *m) {
#line 407 "sgesvj.f"
	*info = -5;
#line 408 "sgesvj.f"
    } else if (*lda < *m) {
#line 409 "sgesvj.f"
	*info = -7;
#line 410 "sgesvj.f"
    } else if (*mv < 0) {
#line 411 "sgesvj.f"
	*info = -9;
#line 412 "sgesvj.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 414 "sgesvj.f"
	*info = -11;
#line 415 "sgesvj.f"
    } else if (uctol && work[1] <= 1.) {
#line 416 "sgesvj.f"
	*info = -12;
#line 417 "sgesvj.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 417 "sgesvj.f"
	i__1 = *m + *n;
#line 417 "sgesvj.f"
	if (*lwork < max(i__1,6)) {
#line 418 "sgesvj.f"
	    *info = -13;
#line 419 "sgesvj.f"
	} else {
#line 420 "sgesvj.f"
	    *info = 0;
#line 421 "sgesvj.f"
	}
#line 421 "sgesvj.f"
    }

/*     #:( */
#line 424 "sgesvj.f"
    if (*info != 0) {
#line 425 "sgesvj.f"
	i__1 = -(*info);
#line 425 "sgesvj.f"
	xerbla_("SGESVJ", &i__1, (ftnlen)6);
#line 426 "sgesvj.f"
	return 0;
#line 427 "sgesvj.f"
    }

/* #:) Quick return for void matrix */

#line 431 "sgesvj.f"
    if (*m == 0 || *n == 0) {
#line 431 "sgesvj.f"
	return 0;
#line 431 "sgesvj.f"
    }

/*     Set numerical parameters */
/*     The stopping criterion for Jacobi rotations is */

/*     max_{i<>j}|A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS */

/*     where EPS is the round-off and CTOL is defined as follows: */

#line 440 "sgesvj.f"
    if (uctol) {
/*        ... user controlled */
#line 442 "sgesvj.f"
	ctol = work[1];
#line 443 "sgesvj.f"
    } else {
/*        ... default */
#line 445 "sgesvj.f"
	if (lsvec || rsvec || applv) {
#line 446 "sgesvj.f"
	    ctol = sqrt((doublereal) (*m));
#line 447 "sgesvj.f"
	} else {
#line 448 "sgesvj.f"
	    ctol = (doublereal) (*m);
#line 449 "sgesvj.f"
	}
#line 450 "sgesvj.f"
    }
/*     ... and the machine dependent parameters are */
/* [!]  (Make sure that SLAMCH() works properly on the target machine.) */

#line 454 "sgesvj.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 455 "sgesvj.f"
    rooteps = sqrt(epsln);
#line 456 "sgesvj.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 457 "sgesvj.f"
    rootsfmin = sqrt(sfmin);
#line 458 "sgesvj.f"
    small = sfmin / epsln;
#line 459 "sgesvj.f"
    big = slamch_("Overflow", (ftnlen)8);
/*     BIG         = ONE    / SFMIN */
#line 461 "sgesvj.f"
    rootbig = 1. / rootsfmin;
#line 462 "sgesvj.f"
    large = big / sqrt((doublereal) (*m * *n));
#line 463 "sgesvj.f"
    bigtheta = 1. / rooteps;

#line 465 "sgesvj.f"
    tol = ctol * epsln;
#line 466 "sgesvj.f"
    roottol = sqrt(tol);

#line 468 "sgesvj.f"
    if ((doublereal) (*m) * epsln >= 1.) {
#line 469 "sgesvj.f"
	*info = -4;
#line 470 "sgesvj.f"
	i__1 = -(*info);
#line 470 "sgesvj.f"
	xerbla_("SGESVJ", &i__1, (ftnlen)6);
#line 471 "sgesvj.f"
	return 0;
#line 472 "sgesvj.f"
    }

/*     Initialize the right singular vector matrix. */

#line 476 "sgesvj.f"
    if (rsvec) {
#line 477 "sgesvj.f"
	mvl = *n;
#line 478 "sgesvj.f"
	slaset_("A", &mvl, n, &c_b17, &c_b18, &v[v_offset], ldv, (ftnlen)1);
#line 479 "sgesvj.f"
    } else if (applv) {
#line 480 "sgesvj.f"
	mvl = *mv;
#line 481 "sgesvj.f"
    }
#line 482 "sgesvj.f"
    rsvec = rsvec || applv;

/*     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
/* (!)  If necessary, scale A to protect the largest singular value */
/*     from overflow. It is possible that saving the largest singular */
/*     value destroys the information about the small ones. */
/*     This initial scaling is almost minimal in the sense that the */
/*     goal is to make sure that no column norm overflows, and that */
/*     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
/*     in A are detected, the procedure returns with INFO=-6. */

#line 493 "sgesvj.f"
    skl = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 494 "sgesvj.f"
    noscale = TRUE_;
#line 495 "sgesvj.f"
    goscale = TRUE_;

#line 497 "sgesvj.f"
    if (lower) {
/*        the input matrix is M-by-N lower triangular (trapezoidal) */
#line 499 "sgesvj.f"
	i__1 = *n;
#line 499 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 500 "sgesvj.f"
	    aapp = 0.;
#line 501 "sgesvj.f"
	    aaqq = 1.;
#line 502 "sgesvj.f"
	    i__2 = *m - p + 1;
#line 502 "sgesvj.f"
	    slassq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
#line 503 "sgesvj.f"
	    if (aapp > big) {
#line 504 "sgesvj.f"
		*info = -6;
#line 505 "sgesvj.f"
		i__2 = -(*info);
#line 505 "sgesvj.f"
		xerbla_("SGESVJ", &i__2, (ftnlen)6);
#line 506 "sgesvj.f"
		return 0;
#line 507 "sgesvj.f"
	    }
#line 508 "sgesvj.f"
	    aaqq = sqrt(aaqq);
#line 509 "sgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 510 "sgesvj.f"
		sva[p] = aapp * aaqq;
#line 511 "sgesvj.f"
	    } else {
#line 512 "sgesvj.f"
		noscale = FALSE_;
#line 513 "sgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 514 "sgesvj.f"
		if (goscale) {
#line 515 "sgesvj.f"
		    goscale = FALSE_;
#line 516 "sgesvj.f"
		    i__2 = p - 1;
#line 516 "sgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 517 "sgesvj.f"
			sva[q] *= skl;
#line 518 "sgesvj.f"
/* L1873: */
#line 518 "sgesvj.f"
		    }
#line 519 "sgesvj.f"
		}
#line 520 "sgesvj.f"
	    }
#line 521 "sgesvj.f"
/* L1874: */
#line 521 "sgesvj.f"
	}
#line 522 "sgesvj.f"
    } else if (upper) {
/*        the input matrix is M-by-N upper triangular (trapezoidal) */
#line 524 "sgesvj.f"
	i__1 = *n;
#line 524 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 525 "sgesvj.f"
	    aapp = 0.;
#line 526 "sgesvj.f"
	    aaqq = 1.;
#line 527 "sgesvj.f"
	    slassq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 528 "sgesvj.f"
	    if (aapp > big) {
#line 529 "sgesvj.f"
		*info = -6;
#line 530 "sgesvj.f"
		i__2 = -(*info);
#line 530 "sgesvj.f"
		xerbla_("SGESVJ", &i__2, (ftnlen)6);
#line 531 "sgesvj.f"
		return 0;
#line 532 "sgesvj.f"
	    }
#line 533 "sgesvj.f"
	    aaqq = sqrt(aaqq);
#line 534 "sgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 535 "sgesvj.f"
		sva[p] = aapp * aaqq;
#line 536 "sgesvj.f"
	    } else {
#line 537 "sgesvj.f"
		noscale = FALSE_;
#line 538 "sgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 539 "sgesvj.f"
		if (goscale) {
#line 540 "sgesvj.f"
		    goscale = FALSE_;
#line 541 "sgesvj.f"
		    i__2 = p - 1;
#line 541 "sgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 542 "sgesvj.f"
			sva[q] *= skl;
#line 543 "sgesvj.f"
/* L2873: */
#line 543 "sgesvj.f"
		    }
#line 544 "sgesvj.f"
		}
#line 545 "sgesvj.f"
	    }
#line 546 "sgesvj.f"
/* L2874: */
#line 546 "sgesvj.f"
	}
#line 547 "sgesvj.f"
    } else {
/*        the input matrix is M-by-N general dense */
#line 549 "sgesvj.f"
	i__1 = *n;
#line 549 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 550 "sgesvj.f"
	    aapp = 0.;
#line 551 "sgesvj.f"
	    aaqq = 1.;
#line 552 "sgesvj.f"
	    slassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 553 "sgesvj.f"
	    if (aapp > big) {
#line 554 "sgesvj.f"
		*info = -6;
#line 555 "sgesvj.f"
		i__2 = -(*info);
#line 555 "sgesvj.f"
		xerbla_("SGESVJ", &i__2, (ftnlen)6);
#line 556 "sgesvj.f"
		return 0;
#line 557 "sgesvj.f"
	    }
#line 558 "sgesvj.f"
	    aaqq = sqrt(aaqq);
#line 559 "sgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 560 "sgesvj.f"
		sva[p] = aapp * aaqq;
#line 561 "sgesvj.f"
	    } else {
#line 562 "sgesvj.f"
		noscale = FALSE_;
#line 563 "sgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 564 "sgesvj.f"
		if (goscale) {
#line 565 "sgesvj.f"
		    goscale = FALSE_;
#line 566 "sgesvj.f"
		    i__2 = p - 1;
#line 566 "sgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 567 "sgesvj.f"
			sva[q] *= skl;
#line 568 "sgesvj.f"
/* L3873: */
#line 568 "sgesvj.f"
		    }
#line 569 "sgesvj.f"
		}
#line 570 "sgesvj.f"
	    }
#line 571 "sgesvj.f"
/* L3874: */
#line 571 "sgesvj.f"
	}
#line 572 "sgesvj.f"
    }

#line 574 "sgesvj.f"
    if (noscale) {
#line 574 "sgesvj.f"
	skl = 1.;
#line 574 "sgesvj.f"
    }

/*     Move the smaller part of the spectrum from the underflow threshold */
/* (!)  Start by determining the position of the nonzero entries of the */
/*     array SVA() relative to ( SFMIN, BIG ). */

#line 580 "sgesvj.f"
    aapp = 0.;
#line 581 "sgesvj.f"
    aaqq = big;
#line 582 "sgesvj.f"
    i__1 = *n;
#line 582 "sgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 583 "sgesvj.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 583 "sgesvj.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 583 "sgesvj.f"
	    aaqq = min(d__1,d__2);
#line 583 "sgesvj.f"
	}
/* Computing MAX */
#line 584 "sgesvj.f"
	d__1 = aapp, d__2 = sva[p];
#line 584 "sgesvj.f"
	aapp = max(d__1,d__2);
#line 585 "sgesvj.f"
/* L4781: */
#line 585 "sgesvj.f"
    }

/* #:) Quick return for zero matrix */

#line 589 "sgesvj.f"
    if (aapp == 0.) {
#line 590 "sgesvj.f"
	if (lsvec) {
#line 590 "sgesvj.f"
	    slaset_("G", m, n, &c_b17, &c_b18, &a[a_offset], lda, (ftnlen)1);
#line 590 "sgesvj.f"
	}
#line 591 "sgesvj.f"
	work[1] = 1.;
#line 592 "sgesvj.f"
	work[2] = 0.;
#line 593 "sgesvj.f"
	work[3] = 0.;
#line 594 "sgesvj.f"
	work[4] = 0.;
#line 595 "sgesvj.f"
	work[5] = 0.;
#line 596 "sgesvj.f"
	work[6] = 0.;
#line 597 "sgesvj.f"
	return 0;
#line 598 "sgesvj.f"
    }

/* #:) Quick return for one-column matrix */

#line 602 "sgesvj.f"
    if (*n == 1) {
#line 603 "sgesvj.f"
	if (lsvec) {
#line 603 "sgesvj.f"
	    slascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1]
		    , lda, &ierr, (ftnlen)1);
#line 603 "sgesvj.f"
	}
#line 605 "sgesvj.f"
	work[1] = 1. / skl;
#line 606 "sgesvj.f"
	if (sva[1] >= sfmin) {
#line 607 "sgesvj.f"
	    work[2] = 1.;
#line 608 "sgesvj.f"
	} else {
#line 609 "sgesvj.f"
	    work[2] = 0.;
#line 610 "sgesvj.f"
	}
#line 611 "sgesvj.f"
	work[3] = 0.;
#line 612 "sgesvj.f"
	work[4] = 0.;
#line 613 "sgesvj.f"
	work[5] = 0.;
#line 614 "sgesvj.f"
	work[6] = 0.;
#line 615 "sgesvj.f"
	return 0;
#line 616 "sgesvj.f"
    }

/*     Protect small singular values from underflow, and try to */
/*     avoid underflows/overflows in computing Jacobi rotations. */

#line 621 "sgesvj.f"
    sn = sqrt(sfmin / epsln);
#line 622 "sgesvj.f"
    temp1 = sqrt(big / (doublereal) (*n));
#line 623 "sgesvj.f"
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1) {
/* Computing MIN */
#line 625 "sgesvj.f"
	d__1 = big, d__2 = temp1 / aapp;
#line 625 "sgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 628 "sgesvj.f"
    } else if (aaqq <= sn && aapp <= temp1) {
/* Computing MIN */
#line 629 "sgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (aapp * sqrt((doublereal) (*n)));
#line 629 "sgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 632 "sgesvj.f"
    } else if (aaqq >= sn && aapp >= temp1) {
/* Computing MAX */
#line 633 "sgesvj.f"
	d__1 = sn / aaqq, d__2 = temp1 / aapp;
#line 633 "sgesvj.f"
	temp1 = max(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 636 "sgesvj.f"
    } else if (aaqq <= sn && aapp >= temp1) {
/* Computing MIN */
#line 637 "sgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (sqrt((doublereal) (*n)) * aapp);
#line 637 "sgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 640 "sgesvj.f"
    } else {
#line 641 "sgesvj.f"
	temp1 = 1.;
#line 642 "sgesvj.f"
    }

/*     Scale, if necessary */

#line 646 "sgesvj.f"
    if (temp1 != 1.) {
#line 647 "sgesvj.f"
	slascl_("G", &c__0, &c__0, &c_b18, &temp1, n, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 648 "sgesvj.f"
    }
#line 649 "sgesvj.f"
    skl = temp1 * skl;
#line 650 "sgesvj.f"
    if (skl != 1.) {
#line 651 "sgesvj.f"
	slascl_(joba, &c__0, &c__0, &c_b18, &skl, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 652 "sgesvj.f"
	skl = 1. / skl;
#line 653 "sgesvj.f"
    }

/*     Row-cyclic Jacobi SVD algorithm with column pivoting */

#line 657 "sgesvj.f"
    emptsw = *n * (*n - 1) / 2;
#line 658 "sgesvj.f"
    notrot = 0;
#line 659 "sgesvj.f"
    fastr[0] = 0.;

/*     A is represented in factored form A = A * diag(WORK), where diag(WORK) */
/*     is initialized to identity. WORK is updated during fast scaled */
/*     rotations. */

#line 665 "sgesvj.f"
    i__1 = *n;
#line 665 "sgesvj.f"
    for (q = 1; q <= i__1; ++q) {
#line 666 "sgesvj.f"
	work[q] = 1.;
#line 667 "sgesvj.f"
/* L1868: */
#line 667 "sgesvj.f"
    }


#line 670 "sgesvj.f"
    swband = 3;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if SGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 678 "sgesvj.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 684 "sgesvj.f"
    nbl = *n / kbl;
#line 685 "sgesvj.f"
    if (nbl * kbl != *n) {
#line 685 "sgesvj.f"
	++nbl;
#line 685 "sgesvj.f"
    }

/* Computing 2nd power */
#line 687 "sgesvj.f"
    i__1 = kbl;
#line 687 "sgesvj.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 690 "sgesvj.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 693 "sgesvj.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */

/* Computing MAX */
#line 701 "sgesvj.f"
    i__1 = 64, i__2 = kbl << 2;
#line 701 "sgesvj.f"
    if ((lower || upper) && *n > max(i__1,i__2)) {
/* [TP] The number of partition levels and the actual partition are */
/*     tuning parameters. */
#line 704 "sgesvj.f"
	n4 = *n / 4;
#line 705 "sgesvj.f"
	n2 = *n / 2;
#line 706 "sgesvj.f"
	n34 = n4 * 3;
#line 707 "sgesvj.f"
	if (applv) {
#line 708 "sgesvj.f"
	    q = 0;
#line 709 "sgesvj.f"
	} else {
#line 710 "sgesvj.f"
	    q = 1;
#line 711 "sgesvj.f"
	}

#line 713 "sgesvj.f"
	if (lower) {

/*     This works very well on lower triangular matrices, in particular */
/*     in the framework of the preconditioned Jacobi SVD (xGEJSV). */
/*     The idea is simple: */
/*     [+ 0 0 0]   Note that Jacobi transformations of [0 0] */
/*     [+ + 0 0]                                       [0 0] */
/*     [+ + x 0]   actually work on [x 0]              [x 0] */
/*     [+ + x x]                    [x x].             [x x] */

#line 723 "sgesvj.f"
	    i__1 = *m - n34;
#line 723 "sgesvj.f"
	    i__2 = *n - n34;
#line 723 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 723 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda,
		     &work[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + (
		    n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &
		    work[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 728 "sgesvj.f"
	    i__1 = *m - n2;
#line 728 "sgesvj.f"
	    i__2 = n34 - n2;
#line 728 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 728 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &
		    work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1)
		     * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &work[*n + 
		    1], &i__3, &ierr, (ftnlen)1);

#line 733 "sgesvj.f"
	    i__1 = *m - n2;
#line 733 "sgesvj.f"
	    i__2 = *n - n2;
#line 733 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 733 "sgesvj.f"
	    sgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], 
		    lda, &work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (
		    n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 738 "sgesvj.f"
	    i__1 = *m - n4;
#line 738 "sgesvj.f"
	    i__2 = n2 - n4;
#line 738 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 738 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, &
		    work[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1)
		     * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 
		    1], &i__3, &ierr, (ftnlen)1);

#line 743 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 743 "sgesvj.f"
	    sgsvj0_(jobv, m, &n4, &a[a_offset], lda, &work[1], &sva[1], &mvl, 
		    &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n 
		    + 1], &i__1, &ierr, (ftnlen)1);

#line 747 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 747 "sgesvj.f"
	    sgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);


#line 752 "sgesvj.f"
	} else if (upper) {


#line 755 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 755 "sgesvj.f"
	    sgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &work[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 759 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 759 "sgesvj.f"
	    sgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &work[n4 
		    + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], 
		    &i__1, &ierr, (ftnlen)1);

#line 764 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 764 "sgesvj.f"
	    sgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1],
		     &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 768 "sgesvj.f"
	    i__1 = n2 + n4;
#line 768 "sgesvj.f"
	    i__2 = *lwork - *n;
#line 768 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &work[
		    n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], 
		    &i__2, &ierr, (ftnlen)1);
#line 773 "sgesvj.f"
	}

#line 775 "sgesvj.f"
    }

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 779 "sgesvj.f"
    for (i__ = 1; i__ <= 30; ++i__) {

/*     .. go go go ... */

#line 783 "sgesvj.f"
	mxaapq = 0.;
#line 784 "sgesvj.f"
	mxsinj = 0.;
#line 785 "sgesvj.f"
	iswrot = 0;

#line 787 "sgesvj.f"
	notrot = 0;
#line 788 "sgesvj.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 795 "sgesvj.f"
	i__1 = nbl;
#line 795 "sgesvj.f"
	for (ibr = 1; ibr <= i__1; ++ibr) {

#line 797 "sgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 799 "sgesvj.f"
	    i__3 = lkahead, i__4 = nbl - ibr;
#line 799 "sgesvj.f"
	    i__2 = min(i__3,i__4);
#line 799 "sgesvj.f"
	    for (ir1 = 0; ir1 <= i__2; ++ir1) {

#line 801 "sgesvj.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 803 "sgesvj.f"
		i__4 = igl + kbl - 1, i__5 = *n - 1;
#line 803 "sgesvj.f"
		i__3 = min(i__4,i__5);
#line 803 "sgesvj.f"
		for (p = igl; p <= i__3; ++p) {

/*     .. de Rijk's pivoting */

#line 807 "sgesvj.f"
		    i__4 = *n - p + 1;
#line 807 "sgesvj.f"
		    q = isamax_(&i__4, &sva[p], &c__1) + p - 1;
#line 808 "sgesvj.f"
		    if (p != q) {
#line 809 "sgesvj.f"
			sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 810 "sgesvj.f"
			if (rsvec) {
#line 810 "sgesvj.f"
			    sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 810 "sgesvj.f"
			}
#line 812 "sgesvj.f"
			temp1 = sva[p];
#line 813 "sgesvj.f"
			sva[p] = sva[q];
#line 814 "sgesvj.f"
			sva[q] = temp1;
#line 815 "sgesvj.f"
			temp1 = work[p];
#line 816 "sgesvj.f"
			work[p] = work[q];
#line 817 "sgesvj.f"
			work[q] = temp1;
#line 818 "sgesvj.f"
		    }

#line 820 "sgesvj.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/*        Caveat: */
/*        Unfortunately, some BLAS implementations compute SNRM2(M,A(1,p),1) */
/*        as SQRT(SDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to */
/*        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
/*        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
/*        Hence, SNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented SNRM2 is available, the IF-THEN-ELSE */
/*        below should read "AAPP = SNRM2( M, A(1,p), 1 ) * WORK(p)". */

#line 834 "sgesvj.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 836 "sgesvj.f"
			    sva[p] = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * 
				    work[p];
#line 837 "sgesvj.f"
			} else {
#line 838 "sgesvj.f"
			    temp1 = 0.;
#line 839 "sgesvj.f"
			    aapp = 1.;
#line 840 "sgesvj.f"
			    slassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 841 "sgesvj.f"
			    sva[p] = temp1 * sqrt(aapp) * work[p];
#line 842 "sgesvj.f"
			}
#line 843 "sgesvj.f"
			aapp = sva[p];
#line 844 "sgesvj.f"
		    } else {
#line 845 "sgesvj.f"
			aapp = sva[p];
#line 846 "sgesvj.f"
		    }

#line 848 "sgesvj.f"
		    if (aapp > 0.) {

#line 850 "sgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 852 "sgesvj.f"
			i__5 = igl + kbl - 1;
#line 852 "sgesvj.f"
			i__4 = min(i__5,*n);
#line 852 "sgesvj.f"
			for (q = p + 1; q <= i__4; ++q) {

#line 854 "sgesvj.f"
			    aaqq = sva[q];

#line 856 "sgesvj.f"
			    if (aaqq > 0.) {

#line 858 "sgesvj.f"
				aapp0 = aapp;
#line 859 "sgesvj.f"
				if (aaqq >= 1.) {
#line 860 "sgesvj.f"
				    rotok = small * aapp <= aaqq;
#line 861 "sgesvj.f"
				    if (aapp < big / aaqq) {
#line 862 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 865 "sgesvj.f"
				    } else {
#line 866 "sgesvj.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 868 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						work[p], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 871 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[q * a_dim1 + 1], &c__1) * 
						work[q] / aaqq;
#line 873 "sgesvj.f"
				    }
#line 874 "sgesvj.f"
				} else {
#line 875 "sgesvj.f"
				    rotok = aapp <= aaqq / small;
#line 876 "sgesvj.f"
				    if (aapp > small / aaqq) {
#line 877 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 880 "sgesvj.f"
				    } else {
#line 881 "sgesvj.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 883 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						work[q], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 886 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[p * a_dim1 + 1], &c__1) * 
						work[p] / aapp;
#line 888 "sgesvj.f"
				    }
#line 889 "sgesvj.f"
				}

/* Computing MAX */
#line 891 "sgesvj.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 891 "sgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 895 "sgesvj.f"
				if (abs(aapq) > tol) {

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 900 "sgesvj.f"
				    if (ir1 == 0) {
#line 901 "sgesvj.f"
					notrot = 0;
#line 902 "sgesvj.f"
					pskipped = 0;
#line 903 "sgesvj.f"
					++iswrot;
#line 904 "sgesvj.f"
				    }

#line 906 "sgesvj.f"
				    if (rotok) {

#line 908 "sgesvj.f"
					aqoap = aaqq / aapp;
#line 909 "sgesvj.f"
					apoaq = aapp / aaqq;
#line 910 "sgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;

#line 912 "sgesvj.f"
					if (abs(theta) > bigtheta) {

#line 914 "sgesvj.f"
					    t = .5 / theta;
#line 915 "sgesvj.f"
					    fastr[2] = t * work[p] / work[q];
#line 916 "sgesvj.f"
					    fastr[3] = -t * work[q] / work[p];
#line 918 "sgesvj.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 920 "sgesvj.f"
					    if (rsvec) {
#line 920 "sgesvj.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 920 "sgesvj.f"
					    }
/* Computing MAX */
#line 924 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 924 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 926 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 926 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 928 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 928 "sgesvj.f"
					    mxsinj = max(d__1,d__2);

#line 930 "sgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 934 "sgesvj.f"
					    thsign = -d_sign(&c_b18, &aapq);
#line 935 "sgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 937 "sgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 938 "sgesvj.f"
					    sn = t * cs;

/* Computing MAX */
#line 940 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 940 "sgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 941 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 941 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 943 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 943 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 946 "sgesvj.f"
					    apoaq = work[p] / work[q];
#line 947 "sgesvj.f"
					    aqoap = work[q] / work[p];
#line 948 "sgesvj.f"
					    if (work[p] >= 1.) {
#line 949 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 950 "sgesvj.f"
			      fastr[2] = t * apoaq;
#line 951 "sgesvj.f"
			      fastr[3] = -t * aqoap;
#line 952 "sgesvj.f"
			      work[p] *= cs;
#line 953 "sgesvj.f"
			      work[q] *= cs;
#line 954 "sgesvj.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 957 "sgesvj.f"
			      if (rsvec) {
#line 957 "sgesvj.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 957 "sgesvj.f"
			      }
#line 960 "sgesvj.f"
			  } else {
#line 961 "sgesvj.f"
			      d__1 = -t * aqoap;
#line 961 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 964 "sgesvj.f"
			      d__1 = cs * sn * apoaq;
#line 964 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 967 "sgesvj.f"
			      work[p] *= cs;
#line 968 "sgesvj.f"
			      work[q] /= cs;
#line 969 "sgesvj.f"
			      if (rsvec) {
#line 970 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 970 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 973 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 973 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 977 "sgesvj.f"
			      }
#line 978 "sgesvj.f"
			  }
#line 979 "sgesvj.f"
					    } else {
#line 980 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 981 "sgesvj.f"
			      d__1 = t * apoaq;
#line 981 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 984 "sgesvj.f"
			      d__1 = -cs * sn * aqoap;
#line 984 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 987 "sgesvj.f"
			      work[p] /= cs;
#line 988 "sgesvj.f"
			      work[q] *= cs;
#line 989 "sgesvj.f"
			      if (rsvec) {
#line 990 "sgesvj.f"
				  d__1 = t * apoaq;
#line 990 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 993 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 993 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 997 "sgesvj.f"
			      }
#line 998 "sgesvj.f"
			  } else {
#line 999 "sgesvj.f"
			      if (work[p] >= work[q]) {
#line 1001 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 1001 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1004 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1004 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1007 "sgesvj.f"
				  work[p] *= cs;
#line 1008 "sgesvj.f"
				  work[q] /= cs;
#line 1009 "sgesvj.f"
				  if (rsvec) {
#line 1010 "sgesvj.f"
				      d__1 = -t * aqoap;
#line 1010 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1014 "sgesvj.f"
				      d__1 = cs * sn * apoaq;
#line 1014 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1018 "sgesvj.f"
				  }
#line 1019 "sgesvj.f"
			      } else {
#line 1020 "sgesvj.f"
				  d__1 = t * apoaq;
#line 1020 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1023 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1023 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1027 "sgesvj.f"
				  work[p] /= cs;
#line 1028 "sgesvj.f"
				  work[q] *= cs;
#line 1029 "sgesvj.f"
				  if (rsvec) {
#line 1030 "sgesvj.f"
				      d__1 = t * apoaq;
#line 1030 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1033 "sgesvj.f"
				      d__1 = -cs * sn * aqoap;
#line 1033 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1037 "sgesvj.f"
				  }
#line 1038 "sgesvj.f"
			      }
#line 1039 "sgesvj.f"
			  }
#line 1040 "sgesvj.f"
					    }
#line 1041 "sgesvj.f"
					}

#line 1043 "sgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 1045 "sgesvj.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1047 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						c_b18, m, &c__1, &work[*n + 1]
						, lda, &ierr, (ftnlen)1);
#line 1050 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						c_b18, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 1052 "sgesvj.f"
					temp1 = -aapq * work[p] / work[q];
#line 1053 "sgesvj.f"
					saxpy_(m, &temp1, &work[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1055 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &c_b18, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 1057 "sgesvj.f"
					d__1 = 0., d__2 = 1. - aapq * aapq;
#line 1057 "sgesvj.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 1059 "sgesvj.f"
					mxsinj = max(mxsinj,sfmin);
#line 1060 "sgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 1066 "sgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1066 "sgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1068 "sgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1070 "sgesvj.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * work[q];
#line 1072 "sgesvj.f"
					} else {
#line 1073 "sgesvj.f"
					    t = 0.;
#line 1074 "sgesvj.f"
					    aaqq = 1.;
#line 1075 "sgesvj.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1077 "sgesvj.f"
					    sva[q] = t * sqrt(aaqq) * work[q];
#line 1078 "sgesvj.f"
					}
#line 1079 "sgesvj.f"
				    }
#line 1080 "sgesvj.f"
				    if (aapp / aapp0 <= rooteps) {
#line 1081 "sgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1083 "sgesvj.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * work[p];
#line 1085 "sgesvj.f"
					} else {
#line 1086 "sgesvj.f"
					    t = 0.;
#line 1087 "sgesvj.f"
					    aapp = 1.;
#line 1088 "sgesvj.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1090 "sgesvj.f"
					    aapp = t * sqrt(aapp) * work[p];
#line 1091 "sgesvj.f"
					}
#line 1092 "sgesvj.f"
					sva[p] = aapp;
#line 1093 "sgesvj.f"
				    }

#line 1095 "sgesvj.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 1097 "sgesvj.f"
				    if (ir1 == 0) {
#line 1097 "sgesvj.f"
					++notrot;
#line 1097 "sgesvj.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1099 "sgesvj.f"
				    ++pskipped;
#line 1100 "sgesvj.f"
				}
#line 1101 "sgesvj.f"
			    } else {
/*        A(:,q) is zero column */
#line 1103 "sgesvj.f"
				if (ir1 == 0) {
#line 1103 "sgesvj.f"
				    ++notrot;
#line 1103 "sgesvj.f"
				}
#line 1104 "sgesvj.f"
				++pskipped;
#line 1105 "sgesvj.f"
			    }

#line 1107 "sgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1109 "sgesvj.f"
				if (ir1 == 0) {
#line 1109 "sgesvj.f"
				    aapp = -aapp;
#line 1109 "sgesvj.f"
				}
#line 1110 "sgesvj.f"
				notrot = 0;
#line 1111 "sgesvj.f"
				goto L2103;
#line 1112 "sgesvj.f"
			    }

#line 1114 "sgesvj.f"
/* L2002: */
#line 1114 "sgesvj.f"
			}
/*     END q-LOOP */

#line 1117 "sgesvj.f"
L2103:
/*     bailed out of q-loop */

#line 1120 "sgesvj.f"
			sva[p] = aapp;

#line 1122 "sgesvj.f"
		    } else {
#line 1123 "sgesvj.f"
			sva[p] = aapp;
#line 1124 "sgesvj.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 1124 "sgesvj.f"
			    i__4 = igl + kbl - 1;
#line 1124 "sgesvj.f"
			    notrot = notrot + min(i__4,*n) - p;
#line 1124 "sgesvj.f"
			}
#line 1126 "sgesvj.f"
		    }

#line 1128 "sgesvj.f"
/* L2001: */
#line 1128 "sgesvj.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 1131 "sgesvj.f"
/* L1002: */
#line 1131 "sgesvj.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 1136 "sgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

#line 1138 "sgesvj.f"
	    i__2 = nbl;
#line 1138 "sgesvj.f"
	    for (jbc = ibr + 1; jbc <= i__2; ++jbc) {

#line 1140 "sgesvj.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 1144 "sgesvj.f"
		ijblsk = 0;
/* Computing MIN */
#line 1145 "sgesvj.f"
		i__4 = igl + kbl - 1;
#line 1145 "sgesvj.f"
		i__3 = min(i__4,*n);
#line 1145 "sgesvj.f"
		for (p = igl; p <= i__3; ++p) {

#line 1147 "sgesvj.f"
		    aapp = sva[p];
#line 1148 "sgesvj.f"
		    if (aapp > 0.) {

#line 1150 "sgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 1152 "sgesvj.f"
			i__5 = jgl + kbl - 1;
#line 1152 "sgesvj.f"
			i__4 = min(i__5,*n);
#line 1152 "sgesvj.f"
			for (q = jgl; q <= i__4; ++q) {

#line 1154 "sgesvj.f"
			    aaqq = sva[q];
#line 1155 "sgesvj.f"
			    if (aaqq > 0.) {
#line 1156 "sgesvj.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 1162 "sgesvj.f"
				if (aaqq >= 1.) {
#line 1163 "sgesvj.f"
				    if (aapp >= aaqq) {
#line 1164 "sgesvj.f"
					rotok = small * aapp <= aaqq;
#line 1165 "sgesvj.f"
				    } else {
#line 1166 "sgesvj.f"
					rotok = small * aaqq <= aapp;
#line 1167 "sgesvj.f"
				    }
#line 1168 "sgesvj.f"
				    if (aapp < big / aaqq) {
#line 1169 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 1172 "sgesvj.f"
				    } else {
#line 1173 "sgesvj.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1175 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						work[p], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1178 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[q * a_dim1 + 1], &c__1) * 
						work[q] / aaqq;
#line 1180 "sgesvj.f"
				    }
#line 1181 "sgesvj.f"
				} else {
#line 1182 "sgesvj.f"
				    if (aapp >= aaqq) {
#line 1183 "sgesvj.f"
					rotok = aapp <= aaqq / small;
#line 1184 "sgesvj.f"
				    } else {
#line 1185 "sgesvj.f"
					rotok = aaqq <= aapp / small;
#line 1186 "sgesvj.f"
				    }
#line 1187 "sgesvj.f"
				    if (aapp > small / aaqq) {
#line 1188 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 1191 "sgesvj.f"
				    } else {
#line 1192 "sgesvj.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1194 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						work[q], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1197 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[p * a_dim1 + 1], &c__1) * 
						work[p] / aapp;
#line 1199 "sgesvj.f"
				    }
#line 1200 "sgesvj.f"
				}

/* Computing MAX */
#line 1202 "sgesvj.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 1202 "sgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 1206 "sgesvj.f"
				if (abs(aapq) > tol) {
#line 1207 "sgesvj.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 1209 "sgesvj.f"
				    pskipped = 0;
#line 1210 "sgesvj.f"
				    ++iswrot;

#line 1212 "sgesvj.f"
				    if (rotok) {

#line 1214 "sgesvj.f"
					aqoap = aaqq / aapp;
#line 1215 "sgesvj.f"
					apoaq = aapp / aaqq;
#line 1216 "sgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 1217 "sgesvj.f"
					if (aaqq > aapp0) {
#line 1217 "sgesvj.f"
					    theta = -theta;
#line 1217 "sgesvj.f"
					}

#line 1219 "sgesvj.f"
					if (abs(theta) > bigtheta) {
#line 1220 "sgesvj.f"
					    t = .5 / theta;
#line 1221 "sgesvj.f"
					    fastr[2] = t * work[p] / work[q];
#line 1222 "sgesvj.f"
					    fastr[3] = -t * work[q] / work[p];
#line 1224 "sgesvj.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 1226 "sgesvj.f"
					    if (rsvec) {
#line 1226 "sgesvj.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 1226 "sgesvj.f"
					    }
/* Computing MAX */
#line 1230 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 1230 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1232 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 1232 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 1234 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 1234 "sgesvj.f"
					    mxsinj = max(d__1,d__2);
#line 1235 "sgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 1239 "sgesvj.f"
					    thsign = -d_sign(&c_b18, &aapq);
#line 1240 "sgesvj.f"
					    if (aaqq > aapp0) {
#line 1240 "sgesvj.f"
			  thsign = -thsign;
#line 1240 "sgesvj.f"
					    }
#line 1241 "sgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 1243 "sgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 1244 "sgesvj.f"
					    sn = t * cs;
/* Computing MAX */
#line 1245 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 1245 "sgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 1246 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 1246 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1248 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 1248 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 1251 "sgesvj.f"
					    apoaq = work[p] / work[q];
#line 1252 "sgesvj.f"
					    aqoap = work[q] / work[p];
#line 1253 "sgesvj.f"
					    if (work[p] >= 1.) {

#line 1255 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 1256 "sgesvj.f"
			      fastr[2] = t * apoaq;
#line 1257 "sgesvj.f"
			      fastr[3] = -t * aqoap;
#line 1258 "sgesvj.f"
			      work[p] *= cs;
#line 1259 "sgesvj.f"
			      work[q] *= cs;
#line 1260 "sgesvj.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 1263 "sgesvj.f"
			      if (rsvec) {
#line 1263 "sgesvj.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 1263 "sgesvj.f"
			      }
#line 1266 "sgesvj.f"
			  } else {
#line 1267 "sgesvj.f"
			      d__1 = -t * aqoap;
#line 1267 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1270 "sgesvj.f"
			      d__1 = cs * sn * apoaq;
#line 1270 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 1273 "sgesvj.f"
			      if (rsvec) {
#line 1274 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 1274 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1277 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1277 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1281 "sgesvj.f"
			      }
#line 1282 "sgesvj.f"
			      work[p] *= cs;
#line 1283 "sgesvj.f"
			      work[q] /= cs;
#line 1284 "sgesvj.f"
			  }
#line 1285 "sgesvj.f"
					    } else {
#line 1286 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 1287 "sgesvj.f"
			      d__1 = t * apoaq;
#line 1287 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 1290 "sgesvj.f"
			      d__1 = -cs * sn * aqoap;
#line 1290 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1293 "sgesvj.f"
			      if (rsvec) {
#line 1294 "sgesvj.f"
				  d__1 = t * apoaq;
#line 1294 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1297 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1297 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1301 "sgesvj.f"
			      }
#line 1302 "sgesvj.f"
			      work[p] /= cs;
#line 1303 "sgesvj.f"
			      work[q] *= cs;
#line 1304 "sgesvj.f"
			  } else {
#line 1305 "sgesvj.f"
			      if (work[p] >= work[q]) {
#line 1307 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 1307 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1310 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1310 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1313 "sgesvj.f"
				  work[p] *= cs;
#line 1314 "sgesvj.f"
				  work[q] /= cs;
#line 1315 "sgesvj.f"
				  if (rsvec) {
#line 1316 "sgesvj.f"
				      d__1 = -t * aqoap;
#line 1316 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1320 "sgesvj.f"
				      d__1 = cs * sn * apoaq;
#line 1320 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1324 "sgesvj.f"
				  }
#line 1325 "sgesvj.f"
			      } else {
#line 1326 "sgesvj.f"
				  d__1 = t * apoaq;
#line 1326 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1329 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1329 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1333 "sgesvj.f"
				  work[p] /= cs;
#line 1334 "sgesvj.f"
				  work[q] *= cs;
#line 1335 "sgesvj.f"
				  if (rsvec) {
#line 1336 "sgesvj.f"
				      d__1 = t * apoaq;
#line 1336 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1339 "sgesvj.f"
				      d__1 = -cs * sn * aqoap;
#line 1339 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1343 "sgesvj.f"
				  }
#line 1344 "sgesvj.f"
			      }
#line 1345 "sgesvj.f"
			  }
#line 1346 "sgesvj.f"
					    }
#line 1347 "sgesvj.f"
					}

#line 1349 "sgesvj.f"
				    } else {
#line 1350 "sgesvj.f"
					if (aapp > aaqq) {
#line 1351 "sgesvj.f"
					    scopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[*n + 1], &
						    c__1);
#line 1353 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &work[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1356 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1359 "sgesvj.f"
					    temp1 = -aapq * work[p] / work[q];
#line 1360 "sgesvj.f"
					    saxpy_(m, &temp1, &work[*n + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1);
#line 1362 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &c_b18,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1365 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 1365 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 1367 "sgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1368 "sgesvj.f"
					} else {
#line 1369 "sgesvj.f"
					    scopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[*n + 1], &
						    c__1);
#line 1371 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &work[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1374 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1377 "sgesvj.f"
					    temp1 = -aapq * work[q] / work[p];
#line 1378 "sgesvj.f"
					    saxpy_(m, &temp1, &work[*n + 1], &
						    c__1, &a[p * a_dim1 + 1], 
						    &c__1);
#line 1380 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &c_b18,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1383 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 1383 "sgesvj.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 1385 "sgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1386 "sgesvj.f"
					}
#line 1387 "sgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 1392 "sgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1392 "sgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1394 "sgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1396 "sgesvj.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * work[q];
#line 1398 "sgesvj.f"
					} else {
#line 1399 "sgesvj.f"
					    t = 0.;
#line 1400 "sgesvj.f"
					    aaqq = 1.;
#line 1401 "sgesvj.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1403 "sgesvj.f"
					    sva[q] = t * sqrt(aaqq) * work[q];
#line 1404 "sgesvj.f"
					}
#line 1405 "sgesvj.f"
				    }
/* Computing 2nd power */
#line 1406 "sgesvj.f"
				    d__1 = aapp / aapp0;
#line 1406 "sgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1407 "sgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1409 "sgesvj.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * work[p];
#line 1411 "sgesvj.f"
					} else {
#line 1412 "sgesvj.f"
					    t = 0.;
#line 1413 "sgesvj.f"
					    aapp = 1.;
#line 1414 "sgesvj.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1416 "sgesvj.f"
					    aapp = t * sqrt(aapp) * work[p];
#line 1417 "sgesvj.f"
					}
#line 1418 "sgesvj.f"
					sva[p] = aapp;
#line 1419 "sgesvj.f"
				    }
/*              end of OK rotation */
#line 1421 "sgesvj.f"
				} else {
#line 1422 "sgesvj.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1424 "sgesvj.f"
				    ++pskipped;
#line 1425 "sgesvj.f"
				    ++ijblsk;
#line 1426 "sgesvj.f"
				}
#line 1427 "sgesvj.f"
			    } else {
#line 1428 "sgesvj.f"
				++notrot;
#line 1429 "sgesvj.f"
				++pskipped;
#line 1430 "sgesvj.f"
				++ijblsk;
#line 1431 "sgesvj.f"
			    }

#line 1433 "sgesvj.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 1435 "sgesvj.f"
				sva[p] = aapp;
#line 1436 "sgesvj.f"
				notrot = 0;
#line 1437 "sgesvj.f"
				goto L2011;
#line 1438 "sgesvj.f"
			    }
#line 1439 "sgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1441 "sgesvj.f"
				aapp = -aapp;
#line 1442 "sgesvj.f"
				notrot = 0;
#line 1443 "sgesvj.f"
				goto L2203;
#line 1444 "sgesvj.f"
			    }

#line 1446 "sgesvj.f"
/* L2200: */
#line 1446 "sgesvj.f"
			}
/*        end of the q-loop */
#line 1448 "sgesvj.f"
L2203:

#line 1450 "sgesvj.f"
			sva[p] = aapp;

#line 1452 "sgesvj.f"
		    } else {

#line 1454 "sgesvj.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1454 "sgesvj.f"
			    i__4 = jgl + kbl - 1;
#line 1454 "sgesvj.f"
			    notrot = notrot + min(i__4,*n) - jgl + 1;
#line 1454 "sgesvj.f"
			}
#line 1456 "sgesvj.f"
			if (aapp < 0.) {
#line 1456 "sgesvj.f"
			    notrot = 0;
#line 1456 "sgesvj.f"
			}

#line 1458 "sgesvj.f"
		    }

#line 1460 "sgesvj.f"
/* L2100: */
#line 1460 "sgesvj.f"
		}
/*     end of the p-loop */
#line 1462 "sgesvj.f"
/* L2010: */
#line 1462 "sgesvj.f"
	    }
/*     end of the jbc-loop */
#line 1464 "sgesvj.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1466 "sgesvj.f"
	    i__3 = igl + kbl - 1;
#line 1466 "sgesvj.f"
	    i__2 = min(i__3,*n);
#line 1466 "sgesvj.f"
	    for (p = igl; p <= i__2; ++p) {
#line 1467 "sgesvj.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1468 "sgesvj.f"
/* L2012: */
#line 1468 "sgesvj.f"
	    }
/* ** */
#line 1470 "sgesvj.f"
/* L2000: */
#line 1470 "sgesvj.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1474 "sgesvj.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1476 "sgesvj.f"
	    sva[*n] = snrm2_(m, &a[*n * a_dim1 + 1], &c__1) * work[*n];
#line 1477 "sgesvj.f"
	} else {
#line 1478 "sgesvj.f"
	    t = 0.;
#line 1479 "sgesvj.f"
	    aapp = 1.;
#line 1480 "sgesvj.f"
	    slassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1481 "sgesvj.f"
	    sva[*n] = t * sqrt(aapp) * work[*n];
#line 1482 "sgesvj.f"
	}

/*     Additional steering devices */

#line 1486 "sgesvj.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1486 "sgesvj.f"
	    swband = i__;
#line 1486 "sgesvj.f"
	}

#line 1489 "sgesvj.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * tol && (
		doublereal) (*n) * mxaapq * mxsinj < tol) {
#line 1491 "sgesvj.f"
	    goto L1994;
#line 1492 "sgesvj.f"
	}

#line 1494 "sgesvj.f"
	if (notrot >= emptsw) {
#line 1494 "sgesvj.f"
	    goto L1994;
#line 1494 "sgesvj.f"
	}

#line 1496 "sgesvj.f"
/* L1993: */
#line 1496 "sgesvj.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 1500 "sgesvj.f"
    *info = 29;
#line 1501 "sgesvj.f"
    goto L1995;

#line 1503 "sgesvj.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 1507 "sgesvj.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1509 "sgesvj.f"
L1995:

/*     Sort the singular values and find how many are above */
/*     the underflow threshold. */

#line 1514 "sgesvj.f"
    n2 = 0;
#line 1515 "sgesvj.f"
    n4 = 0;
#line 1516 "sgesvj.f"
    i__1 = *n - 1;
#line 1516 "sgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 1517 "sgesvj.f"
	i__2 = *n - p + 1;
#line 1517 "sgesvj.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1518 "sgesvj.f"
	if (p != q) {
#line 1519 "sgesvj.f"
	    temp1 = sva[p];
#line 1520 "sgesvj.f"
	    sva[p] = sva[q];
#line 1521 "sgesvj.f"
	    sva[q] = temp1;
#line 1522 "sgesvj.f"
	    temp1 = work[p];
#line 1523 "sgesvj.f"
	    work[p] = work[q];
#line 1524 "sgesvj.f"
	    work[q] = temp1;
#line 1525 "sgesvj.f"
	    sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1526 "sgesvj.f"
	    if (rsvec) {
#line 1526 "sgesvj.f"
		sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1526 "sgesvj.f"
	    }
#line 1527 "sgesvj.f"
	}
#line 1528 "sgesvj.f"
	if (sva[p] != 0.) {
#line 1529 "sgesvj.f"
	    ++n4;
#line 1530 "sgesvj.f"
	    if (sva[p] * skl > sfmin) {
#line 1530 "sgesvj.f"
		++n2;
#line 1530 "sgesvj.f"
	    }
#line 1531 "sgesvj.f"
	}
#line 1532 "sgesvj.f"
/* L5991: */
#line 1532 "sgesvj.f"
    }
#line 1533 "sgesvj.f"
    if (sva[*n] != 0.) {
#line 1534 "sgesvj.f"
	++n4;
#line 1535 "sgesvj.f"
	if (sva[*n] * skl > sfmin) {
#line 1535 "sgesvj.f"
	    ++n2;
#line 1535 "sgesvj.f"
	}
#line 1536 "sgesvj.f"
    }

/*     Normalize the left singular vectors. */

#line 1540 "sgesvj.f"
    if (lsvec || uctol) {
#line 1541 "sgesvj.f"
	i__1 = n2;
#line 1541 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1542 "sgesvj.f"
	    d__1 = work[p] / sva[p];
#line 1542 "sgesvj.f"
	    sscal_(m, &d__1, &a[p * a_dim1 + 1], &c__1);
#line 1543 "sgesvj.f"
/* L1998: */
#line 1543 "sgesvj.f"
	}
#line 1544 "sgesvj.f"
    }

/*     Scale the product of Jacobi rotations (assemble the fast rotations). */

#line 1548 "sgesvj.f"
    if (rsvec) {
#line 1549 "sgesvj.f"
	if (applv) {
#line 1550 "sgesvj.f"
	    i__1 = *n;
#line 1550 "sgesvj.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1551 "sgesvj.f"
		sscal_(&mvl, &work[p], &v[p * v_dim1 + 1], &c__1);
#line 1552 "sgesvj.f"
/* L2398: */
#line 1552 "sgesvj.f"
	    }
#line 1553 "sgesvj.f"
	} else {
#line 1554 "sgesvj.f"
	    i__1 = *n;
#line 1554 "sgesvj.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1555 "sgesvj.f"
		temp1 = 1. / snrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
#line 1556 "sgesvj.f"
		sscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
#line 1557 "sgesvj.f"
/* L2399: */
#line 1557 "sgesvj.f"
	    }
#line 1558 "sgesvj.f"
	}
#line 1559 "sgesvj.f"
    }

/*     Undo scaling, if necessary (and possible). */
#line 1562 "sgesvj.f"
    if (skl > 1. && sva[1] < big / skl || skl < 1. && sva[max(n2,1)] > sfmin /
	     skl) {
#line 1565 "sgesvj.f"
	i__1 = *n;
#line 1565 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1566 "sgesvj.f"
	    sva[p] = skl * sva[p];
#line 1567 "sgesvj.f"
/* L2400: */
#line 1567 "sgesvj.f"
	}
#line 1568 "sgesvj.f"
	skl = 1.;
#line 1569 "sgesvj.f"
    }

#line 1571 "sgesvj.f"
    work[1] = skl;
/*     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
/*     then some of the singular values may overflow or underflow and */
/*     the spectrum is given in this factored representation. */

#line 1576 "sgesvj.f"
    work[2] = (doublereal) n4;
/*     N4 is the number of computed nonzero singular values of A. */

#line 1579 "sgesvj.f"
    work[3] = (doublereal) n2;
/*     N2 is the number of singular values of A greater than SFMIN. */
/*     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
/*     that may carry some information. */

#line 1584 "sgesvj.f"
    work[4] = (doublereal) i__;
/*     i is the index of the last sweep before declaring convergence. */

#line 1587 "sgesvj.f"
    work[5] = mxaapq;
/*     MXAAPQ is the largest absolute value of scaled pivots in the */
/*     last sweep */

#line 1591 "sgesvj.f"
    work[6] = mxsinj;
/*     MXSINJ is the largest absolute value of the sines of Jacobi angles */
/*     in the last sweep */

#line 1595 "sgesvj.f"
    return 0;
/*     .. */
/*     .. END OF SGESVJ */
/*     .. */
} /* sgesvj_ */

