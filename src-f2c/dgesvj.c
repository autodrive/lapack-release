#line 1 "dgesvj.f"
/* dgesvj.f -- translated by f2c (version 20100827).
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

#line 1 "dgesvj.f"
/* Table of constant values */

static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/* > \brief \b DGESVJ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGESVJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvj.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvj.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvj.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, */
/*                          LDV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N */
/*       CHARACTER*1        JOBA, JOBU, JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), SVA( N ), V( LDV, * ), */
/*      $                   WORK( LWORK ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGESVJ computes the singular value decomposition (SVD) of a real */
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
/* >                 approximately TOL=CTOL*EPS, CTOL=DSQRT(M), EPS=DLAMCH('E'). */
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
/* >                  matrix V is not computed explicitly, instead it is */
/* >                  applied to an MV-by-N matrix initially stored in the */
/* >                  first MV rows of V. */
/* >          = 'N' : the matrix V is not computed and the array V is not */
/* >                  referenced */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the input matrix A. 1/DLAMCH('E') > M >= 0. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit : */
/* >          If JOBU .EQ. 'U' .OR. JOBU .EQ. 'C' : */
/* >                 If INFO .EQ. 0 : */
/* >                 RANKA orthonormal columns of U are returned in the */
/* >                 leading RANKA columns of the array A. Here RANKA <= N */
/* >                 is the number of computed singular values of A that are */
/* >                 above the underflow threshold DLAMCH('S'). The singular */
/* >                 vectors corresponding to underflowed or zero singular */
/* >                 values are not computed. The value of RANKA is returned */
/* >                 in the array WORK as RANKA=NINT(WORK(2)). Also see the */
/* >                 descriptions of SVA and WORK. The computed columns of U */
/* >                 are mutually numerically orthogonal up to approximately */
/* >                 TOL=DSQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU.EQ.'C'), */
/* >                 see the description of JOBU. */
/* >                 If INFO .GT. 0 : */
/* >                 the procedure DGESVJ did not converge in the given number */
/* >                 of iterations (sweeps). In that case, the computed */
/* >                 columns of U may not be orthogonal up to TOL. The output */
/* >                 U (stored in A), SIGMA (given by the computed singular */
/* >                 values in SVA(1:N)) and V is still a decomposition of the */
/* >                 input matrix A in the sense that the residual */
/* >                 ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small. */
/* > */
/* >          If JOBU .EQ. 'N' : */
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
/* >                 the procedure DGESVJ did not converge in the given number */
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
/* >          SVA is DOUBLE PRECISION array, dimension (N) */
/* >          On exit : */
/* >          If INFO .EQ. 0 : */
/* >          depending on the value SCALE = WORK(1), we have: */
/* >                 If SCALE .EQ. ONE : */
/* >                 SVA(1:N) contains the computed singular values of A. */
/* >                 During the computation SVA contains the Euclidean column */
/* >                 norms of the iterated matrices in the array A. */
/* >                 If SCALE .NE. ONE : */
/* >                 The singular values of A are SCALE*SVA(1:N), and this */
/* >                 factored representation is due to the fact that some of the */
/* >                 singular values of A might underflow or overflow. */
/* >          If INFO .GT. 0 : */
/* >          the procedure DGESVJ did not converge in the given number of */
/* >          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate. */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* >          MV is INTEGER */
/* >          If JOBV .EQ. 'A', then the product of Jacobi rotations in DGESVJ */
/* >          is applied to the first MV rows of V. See the description of JOBV. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension (LDV,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension max(4,M+N). */
/* >          On entry : */
/* >          If JOBU .EQ. 'C' : */
/* >          WORK(1) = CTOL, where CTOL defines the threshold for convergence. */
/* >                    The process stops if all columns of A are mutually */
/* >                    orthogonal up to CTOL*EPS, EPS=DLAMCH('E'). */
/* >                    It is required that CTOL >= ONE, i.e. it is not */
/* >                    allowed to force the routine to obtain orthogonality */
/* >                    below EPS. */
/* >          On exit : */
/* >          WORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N) */
/* >                    are the computed singular values of A. */
/* >                    (See description of SVA().) */
/* >          WORK(2) = NINT(WORK(2)) is the number of the computed nonzero */
/* >                    singular values. */
/* >          WORK(3) = NINT(WORK(3)) is the number of the computed singular */
/* >                    values that are larger than the underflow threshold. */
/* >          WORK(4) = NINT(WORK(4)) is the number of sweeps of Jacobi */
/* >                    rotations needed for numerical convergence. */
/* >          WORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep. */
/* >                    This is useful information in cases when DGESVJ did */
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
/* >          length of WORK, WORK >= MAX(6,M+N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0 : successful exit. */
/* >          < 0 : if INFO = -i, then the i-th argument had an illegal value */
/* >          > 0 : DGESVJ did not converge in the maximal allowed number (30) */
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

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane */
/* >  rotations. The rotations are implemented as fast scaled rotations of */
/* >  Anda and Park [1]. In the case of underflow of the Jacobi angle, a */
/* >  modified Jacobi transformation of Drmac [4] is used. Pivot strategy uses */
/* >  column interchanges of de Rijk [2]. The relative accuracy of the computed */
/* >  singular values and the accuracy of the computed singular vectors (in */
/* >  angle metric) is as guaranteed by the theory of Demmel and Veselic [3]. */
/* >  The condition number that determines the accuracy in the full rank case */
/* >  is essentially min_{D=diag} kappa(A*D), where kappa(.) is the */
/* >  spectral condition number. The best performance of this Jacobi SVD */
/* >  procedure is achieved if used in an  accelerated version of Drmac and */
/* >  Veselic [5,6], and it is the kernel routine in the SIGMA library [7]. */
/* >  Some tunning parameters (marked with [TP]) are available for the */
/* >  implementer. */
/* >  The computational range for the nonzero singular values is the  machine */
/* >  number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even */
/* >  denormalized singular values can be computed with the corresponding */
/* >  gradual loss of accurate digits. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  ============ */
/* > */
/* >  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* > [1] A. A. Anda and H. Park: Fast plane rotations with dynamic scaling. */
/* >     SIAM J. matrix Anal. Appl., Vol. 15 (1994), pp. 162-174. */
/* > [2] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the */
/* >     singular value decomposition on a vector computer. */
/* >     SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371. */
/* > [3] J. Demmel and K. Veselic: Jacobi method is more accurate than QR. */
/* > [4] Z. Drmac: Implementation of Jacobi rotations for accurate singular */
/* >     value computation in floating point arithmetic. */
/* >     SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222. */
/* > [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. */
/* >     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. */
/* >     LAPACK Working note 169. */
/* > [6] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. */
/* >     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. */
/* >     LAPACK Working note 170. */
/* > [7] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV, */
/* >     QSVD, (H,K)-SVD computations. */
/* >     Department of Mathematics, University of Zagreb, 2008. */
/* > \endverbatim */

/* >  \par Bugs, examples and comments: */
/*   ================================= */
/* > */
/* > \verbatim */
/* >  =========================== */
/* >  Please report all bugs and send interesting test examples and comments to */
/* >  drmac@math.hr. Thank you. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgesvj_(char *joba, char *jobu, char *jobv, integer *m, 
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
    static doublereal aapp, aapq, aaqq;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal ctol;
    static integer ierr;
    static doublereal aapp0;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal large, apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small, sfmin;
    static logical lsvec;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal fastr[5];
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal epsln;
    static logical applv, rsvec;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical uctol;
    extern /* Subroutine */ int drotm_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *);
    static logical lower, upper, rotok;
    extern /* Subroutine */ int dgsvj0_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), dgsvj1_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband, blskip;
    static doublereal mxaapq;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal thsign, mxsinj;
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

#line 405 "dgesvj.f"
    /* Parameter adjustments */
#line 405 "dgesvj.f"
    --sva;
#line 405 "dgesvj.f"
    a_dim1 = *lda;
#line 405 "dgesvj.f"
    a_offset = 1 + a_dim1;
#line 405 "dgesvj.f"
    a -= a_offset;
#line 405 "dgesvj.f"
    v_dim1 = *ldv;
#line 405 "dgesvj.f"
    v_offset = 1 + v_dim1;
#line 405 "dgesvj.f"
    v -= v_offset;
#line 405 "dgesvj.f"
    --work;
#line 405 "dgesvj.f"

#line 405 "dgesvj.f"
    /* Function Body */
#line 405 "dgesvj.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 406 "dgesvj.f"
    uctol = lsame_(jobu, "C", (ftnlen)1, (ftnlen)1);
#line 407 "dgesvj.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 408 "dgesvj.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 409 "dgesvj.f"
    upper = lsame_(joba, "U", (ftnlen)1, (ftnlen)1);
#line 410 "dgesvj.f"
    lower = lsame_(joba, "L", (ftnlen)1, (ftnlen)1);

#line 412 "dgesvj.f"
    if (! (upper || lower || lsame_(joba, "G", (ftnlen)1, (ftnlen)1))) {
#line 413 "dgesvj.f"
	*info = -1;
#line 414 "dgesvj.f"
    } else if (! (lsvec || uctol || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 415 "dgesvj.f"
	*info = -2;
#line 416 "dgesvj.f"
    } else if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 417 "dgesvj.f"
	*info = -3;
#line 418 "dgesvj.f"
    } else if (*m < 0) {
#line 419 "dgesvj.f"
	*info = -4;
#line 420 "dgesvj.f"
    } else if (*n < 0 || *n > *m) {
#line 421 "dgesvj.f"
	*info = -5;
#line 422 "dgesvj.f"
    } else if (*lda < *m) {
#line 423 "dgesvj.f"
	*info = -7;
#line 424 "dgesvj.f"
    } else if (*mv < 0) {
#line 425 "dgesvj.f"
	*info = -9;
#line 426 "dgesvj.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 428 "dgesvj.f"
	*info = -11;
#line 429 "dgesvj.f"
    } else if (uctol && work[1] <= 1.) {
#line 430 "dgesvj.f"
	*info = -12;
#line 431 "dgesvj.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 431 "dgesvj.f"
	i__1 = *m + *n;
#line 431 "dgesvj.f"
	if (*lwork < max(i__1,6)) {
#line 432 "dgesvj.f"
	    *info = -13;
#line 433 "dgesvj.f"
	} else {
#line 434 "dgesvj.f"
	    *info = 0;
#line 435 "dgesvj.f"
	}
#line 435 "dgesvj.f"
    }

/*     #:( */
#line 438 "dgesvj.f"
    if (*info != 0) {
#line 439 "dgesvj.f"
	i__1 = -(*info);
#line 439 "dgesvj.f"
	xerbla_("DGESVJ", &i__1, (ftnlen)6);
#line 440 "dgesvj.f"
	return 0;
#line 441 "dgesvj.f"
    }

/* #:) Quick return for void matrix */

#line 445 "dgesvj.f"
    if (*m == 0 || *n == 0) {
#line 445 "dgesvj.f"
	return 0;
#line 445 "dgesvj.f"
    }

/*     Set numerical parameters */
/*     The stopping criterion for Jacobi rotations is */

/*     max_{i<>j}|A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS */

/*     where EPS is the round-off and CTOL is defined as follows: */

#line 454 "dgesvj.f"
    if (uctol) {
/*        ... user controlled */
#line 456 "dgesvj.f"
	ctol = work[1];
#line 457 "dgesvj.f"
    } else {
/*        ... default */
#line 459 "dgesvj.f"
	if (lsvec || rsvec || applv) {
#line 460 "dgesvj.f"
	    ctol = sqrt((doublereal) (*m));
#line 461 "dgesvj.f"
	} else {
#line 462 "dgesvj.f"
	    ctol = (doublereal) (*m);
#line 463 "dgesvj.f"
	}
#line 464 "dgesvj.f"
    }
/*     ... and the machine dependent parameters are */
/* [!]  (Make sure that DLAMCH() works properly on the target machine.) */

#line 468 "dgesvj.f"
    epsln = dlamch_("Epsilon", (ftnlen)7);
#line 469 "dgesvj.f"
    rooteps = sqrt(epsln);
#line 470 "dgesvj.f"
    sfmin = dlamch_("SafeMinimum", (ftnlen)11);
#line 471 "dgesvj.f"
    rootsfmin = sqrt(sfmin);
#line 472 "dgesvj.f"
    small = sfmin / epsln;
#line 473 "dgesvj.f"
    big = dlamch_("Overflow", (ftnlen)8);
/*     BIG         = ONE    / SFMIN */
#line 475 "dgesvj.f"
    rootbig = 1. / rootsfmin;
#line 476 "dgesvj.f"
    large = big / sqrt((doublereal) (*m * *n));
#line 477 "dgesvj.f"
    bigtheta = 1. / rooteps;

#line 479 "dgesvj.f"
    tol = ctol * epsln;
#line 480 "dgesvj.f"
    roottol = sqrt(tol);

#line 482 "dgesvj.f"
    if ((doublereal) (*m) * epsln >= 1.) {
#line 483 "dgesvj.f"
	*info = -4;
#line 484 "dgesvj.f"
	i__1 = -(*info);
#line 484 "dgesvj.f"
	xerbla_("DGESVJ", &i__1, (ftnlen)6);
#line 485 "dgesvj.f"
	return 0;
#line 486 "dgesvj.f"
    }

/*     Initialize the right singular vector matrix. */

#line 490 "dgesvj.f"
    if (rsvec) {
#line 491 "dgesvj.f"
	mvl = *n;
#line 492 "dgesvj.f"
	dlaset_("A", &mvl, n, &c_b17, &c_b18, &v[v_offset], ldv, (ftnlen)1);
#line 493 "dgesvj.f"
    } else if (applv) {
#line 494 "dgesvj.f"
	mvl = *mv;
#line 495 "dgesvj.f"
    }
#line 496 "dgesvj.f"
    rsvec = rsvec || applv;

/*     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
/* (!)  If necessary, scale A to protect the largest singular value */
/*     from overflow. It is possible that saving the largest singular */
/*     value destroys the information about the small ones. */
/*     This initial scaling is almost minimal in the sense that the */
/*     goal is to make sure that no column norm overflows, and that */
/*     DSQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
/*     in A are detected, the procedure returns with INFO=-6. */

#line 507 "dgesvj.f"
    skl = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 508 "dgesvj.f"
    noscale = TRUE_;
#line 509 "dgesvj.f"
    goscale = TRUE_;

#line 511 "dgesvj.f"
    if (lower) {
/*        the input matrix is M-by-N lower triangular (trapezoidal) */
#line 513 "dgesvj.f"
	i__1 = *n;
#line 513 "dgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 514 "dgesvj.f"
	    aapp = 0.;
#line 515 "dgesvj.f"
	    aaqq = 1.;
#line 516 "dgesvj.f"
	    i__2 = *m - p + 1;
#line 516 "dgesvj.f"
	    dlassq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
#line 517 "dgesvj.f"
	    if (aapp > big) {
#line 518 "dgesvj.f"
		*info = -6;
#line 519 "dgesvj.f"
		i__2 = -(*info);
#line 519 "dgesvj.f"
		xerbla_("DGESVJ", &i__2, (ftnlen)6);
#line 520 "dgesvj.f"
		return 0;
#line 521 "dgesvj.f"
	    }
#line 522 "dgesvj.f"
	    aaqq = sqrt(aaqq);
#line 523 "dgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 524 "dgesvj.f"
		sva[p] = aapp * aaqq;
#line 525 "dgesvj.f"
	    } else {
#line 526 "dgesvj.f"
		noscale = FALSE_;
#line 527 "dgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 528 "dgesvj.f"
		if (goscale) {
#line 529 "dgesvj.f"
		    goscale = FALSE_;
#line 530 "dgesvj.f"
		    i__2 = p - 1;
#line 530 "dgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 531 "dgesvj.f"
			sva[q] *= skl;
#line 532 "dgesvj.f"
/* L1873: */
#line 532 "dgesvj.f"
		    }
#line 533 "dgesvj.f"
		}
#line 534 "dgesvj.f"
	    }
#line 535 "dgesvj.f"
/* L1874: */
#line 535 "dgesvj.f"
	}
#line 536 "dgesvj.f"
    } else if (upper) {
/*        the input matrix is M-by-N upper triangular (trapezoidal) */
#line 538 "dgesvj.f"
	i__1 = *n;
#line 538 "dgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 539 "dgesvj.f"
	    aapp = 0.;
#line 540 "dgesvj.f"
	    aaqq = 1.;
#line 541 "dgesvj.f"
	    dlassq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 542 "dgesvj.f"
	    if (aapp > big) {
#line 543 "dgesvj.f"
		*info = -6;
#line 544 "dgesvj.f"
		i__2 = -(*info);
#line 544 "dgesvj.f"
		xerbla_("DGESVJ", &i__2, (ftnlen)6);
#line 545 "dgesvj.f"
		return 0;
#line 546 "dgesvj.f"
	    }
#line 547 "dgesvj.f"
	    aaqq = sqrt(aaqq);
#line 548 "dgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 549 "dgesvj.f"
		sva[p] = aapp * aaqq;
#line 550 "dgesvj.f"
	    } else {
#line 551 "dgesvj.f"
		noscale = FALSE_;
#line 552 "dgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 553 "dgesvj.f"
		if (goscale) {
#line 554 "dgesvj.f"
		    goscale = FALSE_;
#line 555 "dgesvj.f"
		    i__2 = p - 1;
#line 555 "dgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 556 "dgesvj.f"
			sva[q] *= skl;
#line 557 "dgesvj.f"
/* L2873: */
#line 557 "dgesvj.f"
		    }
#line 558 "dgesvj.f"
		}
#line 559 "dgesvj.f"
	    }
#line 560 "dgesvj.f"
/* L2874: */
#line 560 "dgesvj.f"
	}
#line 561 "dgesvj.f"
    } else {
/*        the input matrix is M-by-N general dense */
#line 563 "dgesvj.f"
	i__1 = *n;
#line 563 "dgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 564 "dgesvj.f"
	    aapp = 0.;
#line 565 "dgesvj.f"
	    aaqq = 1.;
#line 566 "dgesvj.f"
	    dlassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 567 "dgesvj.f"
	    if (aapp > big) {
#line 568 "dgesvj.f"
		*info = -6;
#line 569 "dgesvj.f"
		i__2 = -(*info);
#line 569 "dgesvj.f"
		xerbla_("DGESVJ", &i__2, (ftnlen)6);
#line 570 "dgesvj.f"
		return 0;
#line 571 "dgesvj.f"
	    }
#line 572 "dgesvj.f"
	    aaqq = sqrt(aaqq);
#line 573 "dgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 574 "dgesvj.f"
		sva[p] = aapp * aaqq;
#line 575 "dgesvj.f"
	    } else {
#line 576 "dgesvj.f"
		noscale = FALSE_;
#line 577 "dgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 578 "dgesvj.f"
		if (goscale) {
#line 579 "dgesvj.f"
		    goscale = FALSE_;
#line 580 "dgesvj.f"
		    i__2 = p - 1;
#line 580 "dgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 581 "dgesvj.f"
			sva[q] *= skl;
#line 582 "dgesvj.f"
/* L3873: */
#line 582 "dgesvj.f"
		    }
#line 583 "dgesvj.f"
		}
#line 584 "dgesvj.f"
	    }
#line 585 "dgesvj.f"
/* L3874: */
#line 585 "dgesvj.f"
	}
#line 586 "dgesvj.f"
    }

#line 588 "dgesvj.f"
    if (noscale) {
#line 588 "dgesvj.f"
	skl = 1.;
#line 588 "dgesvj.f"
    }

/*     Move the smaller part of the spectrum from the underflow threshold */
/* (!)  Start by determining the position of the nonzero entries of the */
/*     array SVA() relative to ( SFMIN, BIG ). */

#line 594 "dgesvj.f"
    aapp = 0.;
#line 595 "dgesvj.f"
    aaqq = big;
#line 596 "dgesvj.f"
    i__1 = *n;
#line 596 "dgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 597 "dgesvj.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 597 "dgesvj.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 597 "dgesvj.f"
	    aaqq = min(d__1,d__2);
#line 597 "dgesvj.f"
	}
/* Computing MAX */
#line 598 "dgesvj.f"
	d__1 = aapp, d__2 = sva[p];
#line 598 "dgesvj.f"
	aapp = max(d__1,d__2);
#line 599 "dgesvj.f"
/* L4781: */
#line 599 "dgesvj.f"
    }

/* #:) Quick return for zero matrix */

#line 603 "dgesvj.f"
    if (aapp == 0.) {
#line 604 "dgesvj.f"
	if (lsvec) {
#line 604 "dgesvj.f"
	    dlaset_("G", m, n, &c_b17, &c_b18, &a[a_offset], lda, (ftnlen)1);
#line 604 "dgesvj.f"
	}
#line 605 "dgesvj.f"
	work[1] = 1.;
#line 606 "dgesvj.f"
	work[2] = 0.;
#line 607 "dgesvj.f"
	work[3] = 0.;
#line 608 "dgesvj.f"
	work[4] = 0.;
#line 609 "dgesvj.f"
	work[5] = 0.;
#line 610 "dgesvj.f"
	work[6] = 0.;
#line 611 "dgesvj.f"
	return 0;
#line 612 "dgesvj.f"
    }

/* #:) Quick return for one-column matrix */

#line 616 "dgesvj.f"
    if (*n == 1) {
#line 617 "dgesvj.f"
	if (lsvec) {
#line 617 "dgesvj.f"
	    dlascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1]
		    , lda, &ierr, (ftnlen)1);
#line 617 "dgesvj.f"
	}
#line 619 "dgesvj.f"
	work[1] = 1. / skl;
#line 620 "dgesvj.f"
	if (sva[1] >= sfmin) {
#line 621 "dgesvj.f"
	    work[2] = 1.;
#line 622 "dgesvj.f"
	} else {
#line 623 "dgesvj.f"
	    work[2] = 0.;
#line 624 "dgesvj.f"
	}
#line 625 "dgesvj.f"
	work[3] = 0.;
#line 626 "dgesvj.f"
	work[4] = 0.;
#line 627 "dgesvj.f"
	work[5] = 0.;
#line 628 "dgesvj.f"
	work[6] = 0.;
#line 629 "dgesvj.f"
	return 0;
#line 630 "dgesvj.f"
    }

/*     Protect small singular values from underflow, and try to */
/*     avoid underflows/overflows in computing Jacobi rotations. */

#line 635 "dgesvj.f"
    sn = sqrt(sfmin / epsln);
#line 636 "dgesvj.f"
    temp1 = sqrt(big / (doublereal) (*n));
#line 637 "dgesvj.f"
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1) {
/* Computing MIN */
#line 639 "dgesvj.f"
	d__1 = big, d__2 = temp1 / aapp;
#line 639 "dgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 642 "dgesvj.f"
    } else if (aaqq <= sn && aapp <= temp1) {
/* Computing MIN */
#line 643 "dgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (aapp * sqrt((doublereal) (*n)));
#line 643 "dgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 646 "dgesvj.f"
    } else if (aaqq >= sn && aapp >= temp1) {
/* Computing MAX */
#line 647 "dgesvj.f"
	d__1 = sn / aaqq, d__2 = temp1 / aapp;
#line 647 "dgesvj.f"
	temp1 = max(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 650 "dgesvj.f"
    } else if (aaqq <= sn && aapp >= temp1) {
/* Computing MIN */
#line 651 "dgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (sqrt((doublereal) (*n)) * aapp);
#line 651 "dgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 654 "dgesvj.f"
    } else {
#line 655 "dgesvj.f"
	temp1 = 1.;
#line 656 "dgesvj.f"
    }

/*     Scale, if necessary */

#line 660 "dgesvj.f"
    if (temp1 != 1.) {
#line 661 "dgesvj.f"
	dlascl_("G", &c__0, &c__0, &c_b18, &temp1, n, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 662 "dgesvj.f"
    }
#line 663 "dgesvj.f"
    skl = temp1 * skl;
#line 664 "dgesvj.f"
    if (skl != 1.) {
#line 665 "dgesvj.f"
	dlascl_(joba, &c__0, &c__0, &c_b18, &skl, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 666 "dgesvj.f"
	skl = 1. / skl;
#line 667 "dgesvj.f"
    }

/*     Row-cyclic Jacobi SVD algorithm with column pivoting */

#line 671 "dgesvj.f"
    emptsw = *n * (*n - 1) / 2;
#line 672 "dgesvj.f"
    notrot = 0;
#line 673 "dgesvj.f"
    fastr[0] = 0.;

/*     A is represented in factored form A = A * diag(WORK), where diag(WORK) */
/*     is initialized to identity. WORK is updated during fast scaled */
/*     rotations. */

#line 679 "dgesvj.f"
    i__1 = *n;
#line 679 "dgesvj.f"
    for (q = 1; q <= i__1; ++q) {
#line 680 "dgesvj.f"
	work[q] = 1.;
#line 681 "dgesvj.f"
/* L1868: */
#line 681 "dgesvj.f"
    }


#line 684 "dgesvj.f"
    swband = 3;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if DGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm DGESVJ. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 692 "dgesvj.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 698 "dgesvj.f"
    nbl = *n / kbl;
#line 699 "dgesvj.f"
    if (nbl * kbl != *n) {
#line 699 "dgesvj.f"
	++nbl;
#line 699 "dgesvj.f"
    }

/* Computing 2nd power */
#line 701 "dgesvj.f"
    i__1 = kbl;
#line 701 "dgesvj.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 704 "dgesvj.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 707 "dgesvj.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */

/* Computing MAX */
#line 715 "dgesvj.f"
    i__1 = 64, i__2 = kbl << 2;
#line 715 "dgesvj.f"
    if ((lower || upper) && *n > max(i__1,i__2)) {
/* [TP] The number of partition levels and the actual partition are */
/*     tuning parameters. */
#line 718 "dgesvj.f"
	n4 = *n / 4;
#line 719 "dgesvj.f"
	n2 = *n / 2;
#line 720 "dgesvj.f"
	n34 = n4 * 3;
#line 721 "dgesvj.f"
	if (applv) {
#line 722 "dgesvj.f"
	    q = 0;
#line 723 "dgesvj.f"
	} else {
#line 724 "dgesvj.f"
	    q = 1;
#line 725 "dgesvj.f"
	}

#line 727 "dgesvj.f"
	if (lower) {

/*     This works very well on lower triangular matrices, in particular */
/*     in the framework of the preconditioned Jacobi SVD (xGEJSV). */
/*     The idea is simple: */
/*     [+ 0 0 0]   Note that Jacobi transformations of [0 0] */
/*     [+ + 0 0]                                       [0 0] */
/*     [+ + x 0]   actually work on [x 0]              [x 0] */
/*     [+ + x x]                    [x x].             [x x] */

#line 737 "dgesvj.f"
	    i__1 = *m - n34;
#line 737 "dgesvj.f"
	    i__2 = *n - n34;
#line 737 "dgesvj.f"
	    i__3 = *lwork - *n;
#line 737 "dgesvj.f"
	    dgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda,
		     &work[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + (
		    n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &
		    work[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 742 "dgesvj.f"
	    i__1 = *m - n2;
#line 742 "dgesvj.f"
	    i__2 = n34 - n2;
#line 742 "dgesvj.f"
	    i__3 = *lwork - *n;
#line 742 "dgesvj.f"
	    dgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &
		    work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1)
		     * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &work[*n + 
		    1], &i__3, &ierr, (ftnlen)1);

#line 747 "dgesvj.f"
	    i__1 = *m - n2;
#line 747 "dgesvj.f"
	    i__2 = *n - n2;
#line 747 "dgesvj.f"
	    i__3 = *lwork - *n;
#line 747 "dgesvj.f"
	    dgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], 
		    lda, &work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (
		    n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 752 "dgesvj.f"
	    i__1 = *m - n4;
#line 752 "dgesvj.f"
	    i__2 = n2 - n4;
#line 752 "dgesvj.f"
	    i__3 = *lwork - *n;
#line 752 "dgesvj.f"
	    dgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, &
		    work[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1)
		     * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 
		    1], &i__3, &ierr, (ftnlen)1);

#line 757 "dgesvj.f"
	    i__1 = *lwork - *n;
#line 757 "dgesvj.f"
	    dgsvj0_(jobv, m, &n4, &a[a_offset], lda, &work[1], &sva[1], &mvl, 
		    &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n 
		    + 1], &i__1, &ierr, (ftnlen)1);

#line 761 "dgesvj.f"
	    i__1 = *lwork - *n;
#line 761 "dgesvj.f"
	    dgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);


#line 766 "dgesvj.f"
	} else if (upper) {


#line 769 "dgesvj.f"
	    i__1 = *lwork - *n;
#line 769 "dgesvj.f"
	    dgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &work[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 773 "dgesvj.f"
	    i__1 = *lwork - *n;
#line 773 "dgesvj.f"
	    dgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &work[n4 
		    + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], 
		    &i__1, &ierr, (ftnlen)1);

#line 778 "dgesvj.f"
	    i__1 = *lwork - *n;
#line 778 "dgesvj.f"
	    dgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1],
		     &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 782 "dgesvj.f"
	    i__1 = n2 + n4;
#line 782 "dgesvj.f"
	    i__2 = *lwork - *n;
#line 782 "dgesvj.f"
	    dgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &work[
		    n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], 
		    &i__2, &ierr, (ftnlen)1);
#line 787 "dgesvj.f"
	}

#line 789 "dgesvj.f"
    }

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 793 "dgesvj.f"
    for (i__ = 1; i__ <= 30; ++i__) {

/*     .. go go go ... */

#line 797 "dgesvj.f"
	mxaapq = 0.;
#line 798 "dgesvj.f"
	mxsinj = 0.;
#line 799 "dgesvj.f"
	iswrot = 0;

#line 801 "dgesvj.f"
	notrot = 0;
#line 802 "dgesvj.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 809 "dgesvj.f"
	i__1 = nbl;
#line 809 "dgesvj.f"
	for (ibr = 1; ibr <= i__1; ++ibr) {

#line 811 "dgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 813 "dgesvj.f"
	    i__3 = lkahead, i__4 = nbl - ibr;
#line 813 "dgesvj.f"
	    i__2 = min(i__3,i__4);
#line 813 "dgesvj.f"
	    for (ir1 = 0; ir1 <= i__2; ++ir1) {

#line 815 "dgesvj.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 817 "dgesvj.f"
		i__4 = igl + kbl - 1, i__5 = *n - 1;
#line 817 "dgesvj.f"
		i__3 = min(i__4,i__5);
#line 817 "dgesvj.f"
		for (p = igl; p <= i__3; ++p) {

/*     .. de Rijk's pivoting */

#line 821 "dgesvj.f"
		    i__4 = *n - p + 1;
#line 821 "dgesvj.f"
		    q = idamax_(&i__4, &sva[p], &c__1) + p - 1;
#line 822 "dgesvj.f"
		    if (p != q) {
#line 823 "dgesvj.f"
			dswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 824 "dgesvj.f"
			if (rsvec) {
#line 824 "dgesvj.f"
			    dswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 824 "dgesvj.f"
			}
#line 826 "dgesvj.f"
			temp1 = sva[p];
#line 827 "dgesvj.f"
			sva[p] = sva[q];
#line 828 "dgesvj.f"
			sva[q] = temp1;
#line 829 "dgesvj.f"
			temp1 = work[p];
#line 830 "dgesvj.f"
			work[p] = work[q];
#line 831 "dgesvj.f"
			work[q] = temp1;
#line 832 "dgesvj.f"
		    }

#line 834 "dgesvj.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/*        Caveat: */
/*        Unfortunately, some BLAS implementations compute DNRM2(M,A(1,p),1) */
/*        as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to */
/*        overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and to */
/*        underflow for ||A(:,p)||_2 < DSQRT(underflow_threshold). */
/*        Hence, DNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented DNRM2 is available, the IF-THEN-ELSE */
/*        below should read "AAPP = DNRM2( M, A(1,p), 1 ) * WORK(p)". */

#line 848 "dgesvj.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 850 "dgesvj.f"
			    sva[p] = dnrm2_(m, &a[p * a_dim1 + 1], &c__1) * 
				    work[p];
#line 851 "dgesvj.f"
			} else {
#line 852 "dgesvj.f"
			    temp1 = 0.;
#line 853 "dgesvj.f"
			    aapp = 1.;
#line 854 "dgesvj.f"
			    dlassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 855 "dgesvj.f"
			    sva[p] = temp1 * sqrt(aapp) * work[p];
#line 856 "dgesvj.f"
			}
#line 857 "dgesvj.f"
			aapp = sva[p];
#line 858 "dgesvj.f"
		    } else {
#line 859 "dgesvj.f"
			aapp = sva[p];
#line 860 "dgesvj.f"
		    }

#line 862 "dgesvj.f"
		    if (aapp > 0.) {

#line 864 "dgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 866 "dgesvj.f"
			i__5 = igl + kbl - 1;
#line 866 "dgesvj.f"
			i__4 = min(i__5,*n);
#line 866 "dgesvj.f"
			for (q = p + 1; q <= i__4; ++q) {

#line 868 "dgesvj.f"
			    aaqq = sva[q];

#line 870 "dgesvj.f"
			    if (aaqq > 0.) {

#line 872 "dgesvj.f"
				aapp0 = aapp;
#line 873 "dgesvj.f"
				if (aaqq >= 1.) {
#line 874 "dgesvj.f"
				    rotok = small * aapp <= aaqq;
#line 875 "dgesvj.f"
				    if (aapp < big / aaqq) {
#line 876 "dgesvj.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 879 "dgesvj.f"
				    } else {
#line 880 "dgesvj.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 882 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						work[p], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 885 "dgesvj.f"
					aapq = ddot_(m, &work[*n + 1], &c__1, 
						&a[q * a_dim1 + 1], &c__1) * 
						work[q] / aaqq;
#line 887 "dgesvj.f"
				    }
#line 888 "dgesvj.f"
				} else {
#line 889 "dgesvj.f"
				    rotok = aapp <= aaqq / small;
#line 890 "dgesvj.f"
				    if (aapp > small / aaqq) {
#line 891 "dgesvj.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 894 "dgesvj.f"
				    } else {
#line 895 "dgesvj.f"
					dcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 897 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						work[q], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 900 "dgesvj.f"
					aapq = ddot_(m, &work[*n + 1], &c__1, 
						&a[p * a_dim1 + 1], &c__1) * 
						work[p] / aapp;
#line 902 "dgesvj.f"
				    }
#line 903 "dgesvj.f"
				}

/* Computing MAX */
#line 905 "dgesvj.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 905 "dgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 909 "dgesvj.f"
				if (abs(aapq) > tol) {

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 914 "dgesvj.f"
				    if (ir1 == 0) {
#line 915 "dgesvj.f"
					notrot = 0;
#line 916 "dgesvj.f"
					pskipped = 0;
#line 917 "dgesvj.f"
					++iswrot;
#line 918 "dgesvj.f"
				    }

#line 920 "dgesvj.f"
				    if (rotok) {

#line 922 "dgesvj.f"
					aqoap = aaqq / aapp;
#line 923 "dgesvj.f"
					apoaq = aapp / aaqq;
#line 924 "dgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;

#line 926 "dgesvj.f"
					if (abs(theta) > bigtheta) {

#line 928 "dgesvj.f"
					    t = .5 / theta;
#line 929 "dgesvj.f"
					    fastr[2] = t * work[p] / work[q];
#line 930 "dgesvj.f"
					    fastr[3] = -t * work[q] / work[p];
#line 932 "dgesvj.f"
					    drotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 934 "dgesvj.f"
					    if (rsvec) {
#line 934 "dgesvj.f"
			  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 934 "dgesvj.f"
					    }
/* Computing MAX */
#line 938 "dgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 938 "dgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 940 "dgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 940 "dgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 942 "dgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 942 "dgesvj.f"
					    mxsinj = max(d__1,d__2);

#line 944 "dgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 948 "dgesvj.f"
					    thsign = -d_sign(&c_b18, &aapq);
#line 949 "dgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 951 "dgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 952 "dgesvj.f"
					    sn = t * cs;

/* Computing MAX */
#line 954 "dgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 954 "dgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 955 "dgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 955 "dgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 957 "dgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 957 "dgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 960 "dgesvj.f"
					    apoaq = work[p] / work[q];
#line 961 "dgesvj.f"
					    aqoap = work[q] / work[p];
#line 962 "dgesvj.f"
					    if (work[p] >= 1.) {
#line 963 "dgesvj.f"
			  if (work[q] >= 1.) {
#line 964 "dgesvj.f"
			      fastr[2] = t * apoaq;
#line 965 "dgesvj.f"
			      fastr[3] = -t * aqoap;
#line 966 "dgesvj.f"
			      work[p] *= cs;
#line 967 "dgesvj.f"
			      work[q] *= cs;
#line 968 "dgesvj.f"
			      drotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 971 "dgesvj.f"
			      if (rsvec) {
#line 971 "dgesvj.f"
				  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 971 "dgesvj.f"
			      }
#line 974 "dgesvj.f"
			  } else {
#line 975 "dgesvj.f"
			      d__1 = -t * aqoap;
#line 975 "dgesvj.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 978 "dgesvj.f"
			      d__1 = cs * sn * apoaq;
#line 978 "dgesvj.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 981 "dgesvj.f"
			      work[p] *= cs;
#line 982 "dgesvj.f"
			      work[q] /= cs;
#line 983 "dgesvj.f"
			      if (rsvec) {
#line 984 "dgesvj.f"
				  d__1 = -t * aqoap;
#line 984 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 987 "dgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 987 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 991 "dgesvj.f"
			      }
#line 992 "dgesvj.f"
			  }
#line 993 "dgesvj.f"
					    } else {
#line 994 "dgesvj.f"
			  if (work[q] >= 1.) {
#line 995 "dgesvj.f"
			      d__1 = t * apoaq;
#line 995 "dgesvj.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 998 "dgesvj.f"
			      d__1 = -cs * sn * aqoap;
#line 998 "dgesvj.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1001 "dgesvj.f"
			      work[p] /= cs;
#line 1002 "dgesvj.f"
			      work[q] *= cs;
#line 1003 "dgesvj.f"
			      if (rsvec) {
#line 1004 "dgesvj.f"
				  d__1 = t * apoaq;
#line 1004 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1007 "dgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1007 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1011 "dgesvj.f"
			      }
#line 1012 "dgesvj.f"
			  } else {
#line 1013 "dgesvj.f"
			      if (work[p] >= work[q]) {
#line 1015 "dgesvj.f"
				  d__1 = -t * aqoap;
#line 1015 "dgesvj.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1018 "dgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1018 "dgesvj.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1021 "dgesvj.f"
				  work[p] *= cs;
#line 1022 "dgesvj.f"
				  work[q] /= cs;
#line 1023 "dgesvj.f"
				  if (rsvec) {
#line 1024 "dgesvj.f"
				      d__1 = -t * aqoap;
#line 1024 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1028 "dgesvj.f"
				      d__1 = cs * sn * apoaq;
#line 1028 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1032 "dgesvj.f"
				  }
#line 1033 "dgesvj.f"
			      } else {
#line 1034 "dgesvj.f"
				  d__1 = t * apoaq;
#line 1034 "dgesvj.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1037 "dgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1037 "dgesvj.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1041 "dgesvj.f"
				  work[p] /= cs;
#line 1042 "dgesvj.f"
				  work[q] *= cs;
#line 1043 "dgesvj.f"
				  if (rsvec) {
#line 1044 "dgesvj.f"
				      d__1 = t * apoaq;
#line 1044 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1047 "dgesvj.f"
				      d__1 = -cs * sn * aqoap;
#line 1047 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1051 "dgesvj.f"
				  }
#line 1052 "dgesvj.f"
			      }
#line 1053 "dgesvj.f"
			  }
#line 1054 "dgesvj.f"
					    }
#line 1055 "dgesvj.f"
					}

#line 1057 "dgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 1059 "dgesvj.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1061 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						c_b18, m, &c__1, &work[*n + 1]
						, lda, &ierr, (ftnlen)1);
#line 1064 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						c_b18, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 1066 "dgesvj.f"
					temp1 = -aapq * work[p] / work[q];
#line 1067 "dgesvj.f"
					daxpy_(m, &temp1, &work[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1069 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &c_b18, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 1071 "dgesvj.f"
					d__1 = 0., d__2 = 1. - aapq * aapq;
#line 1071 "dgesvj.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 1073 "dgesvj.f"
					mxsinj = max(mxsinj,sfmin);
#line 1074 "dgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 1080 "dgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1080 "dgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1082 "dgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1084 "dgesvj.f"
					    sva[q] = dnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * work[q];
#line 1086 "dgesvj.f"
					} else {
#line 1087 "dgesvj.f"
					    t = 0.;
#line 1088 "dgesvj.f"
					    aaqq = 1.;
#line 1089 "dgesvj.f"
					    dlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1091 "dgesvj.f"
					    sva[q] = t * sqrt(aaqq) * work[q];
#line 1092 "dgesvj.f"
					}
#line 1093 "dgesvj.f"
				    }
#line 1094 "dgesvj.f"
				    if (aapp / aapp0 <= rooteps) {
#line 1095 "dgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1097 "dgesvj.f"
					    aapp = dnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * work[p];
#line 1099 "dgesvj.f"
					} else {
#line 1100 "dgesvj.f"
					    t = 0.;
#line 1101 "dgesvj.f"
					    aapp = 1.;
#line 1102 "dgesvj.f"
					    dlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1104 "dgesvj.f"
					    aapp = t * sqrt(aapp) * work[p];
#line 1105 "dgesvj.f"
					}
#line 1106 "dgesvj.f"
					sva[p] = aapp;
#line 1107 "dgesvj.f"
				    }

#line 1109 "dgesvj.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 1111 "dgesvj.f"
				    if (ir1 == 0) {
#line 1111 "dgesvj.f"
					++notrot;
#line 1111 "dgesvj.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1113 "dgesvj.f"
				    ++pskipped;
#line 1114 "dgesvj.f"
				}
#line 1115 "dgesvj.f"
			    } else {
/*        A(:,q) is zero column */
#line 1117 "dgesvj.f"
				if (ir1 == 0) {
#line 1117 "dgesvj.f"
				    ++notrot;
#line 1117 "dgesvj.f"
				}
#line 1118 "dgesvj.f"
				++pskipped;
#line 1119 "dgesvj.f"
			    }

#line 1121 "dgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1123 "dgesvj.f"
				if (ir1 == 0) {
#line 1123 "dgesvj.f"
				    aapp = -aapp;
#line 1123 "dgesvj.f"
				}
#line 1124 "dgesvj.f"
				notrot = 0;
#line 1125 "dgesvj.f"
				goto L2103;
#line 1126 "dgesvj.f"
			    }

#line 1128 "dgesvj.f"
/* L2002: */
#line 1128 "dgesvj.f"
			}
/*     END q-LOOP */

#line 1131 "dgesvj.f"
L2103:
/*     bailed out of q-loop */

#line 1134 "dgesvj.f"
			sva[p] = aapp;

#line 1136 "dgesvj.f"
		    } else {
#line 1137 "dgesvj.f"
			sva[p] = aapp;
#line 1138 "dgesvj.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 1138 "dgesvj.f"
			    i__4 = igl + kbl - 1;
#line 1138 "dgesvj.f"
			    notrot = notrot + min(i__4,*n) - p;
#line 1138 "dgesvj.f"
			}
#line 1140 "dgesvj.f"
		    }

#line 1142 "dgesvj.f"
/* L2001: */
#line 1142 "dgesvj.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 1145 "dgesvj.f"
/* L1002: */
#line 1145 "dgesvj.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 1150 "dgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

#line 1152 "dgesvj.f"
	    i__2 = nbl;
#line 1152 "dgesvj.f"
	    for (jbc = ibr + 1; jbc <= i__2; ++jbc) {

#line 1154 "dgesvj.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 1158 "dgesvj.f"
		ijblsk = 0;
/* Computing MIN */
#line 1159 "dgesvj.f"
		i__4 = igl + kbl - 1;
#line 1159 "dgesvj.f"
		i__3 = min(i__4,*n);
#line 1159 "dgesvj.f"
		for (p = igl; p <= i__3; ++p) {

#line 1161 "dgesvj.f"
		    aapp = sva[p];
#line 1162 "dgesvj.f"
		    if (aapp > 0.) {

#line 1164 "dgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 1166 "dgesvj.f"
			i__5 = jgl + kbl - 1;
#line 1166 "dgesvj.f"
			i__4 = min(i__5,*n);
#line 1166 "dgesvj.f"
			for (q = jgl; q <= i__4; ++q) {

#line 1168 "dgesvj.f"
			    aaqq = sva[q];
#line 1169 "dgesvj.f"
			    if (aaqq > 0.) {
#line 1170 "dgesvj.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 1176 "dgesvj.f"
				if (aaqq >= 1.) {
#line 1177 "dgesvj.f"
				    if (aapp >= aaqq) {
#line 1178 "dgesvj.f"
					rotok = small * aapp <= aaqq;
#line 1179 "dgesvj.f"
				    } else {
#line 1180 "dgesvj.f"
					rotok = small * aaqq <= aapp;
#line 1181 "dgesvj.f"
				    }
#line 1182 "dgesvj.f"
				    if (aapp < big / aaqq) {
#line 1183 "dgesvj.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 1186 "dgesvj.f"
				    } else {
#line 1187 "dgesvj.f"
					dcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1189 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &aapp, &
						work[p], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1192 "dgesvj.f"
					aapq = ddot_(m, &work[*n + 1], &c__1, 
						&a[q * a_dim1 + 1], &c__1) * 
						work[q] / aaqq;
#line 1194 "dgesvj.f"
				    }
#line 1195 "dgesvj.f"
				} else {
#line 1196 "dgesvj.f"
				    if (aapp >= aaqq) {
#line 1197 "dgesvj.f"
					rotok = aapp <= aaqq / small;
#line 1198 "dgesvj.f"
				    } else {
#line 1199 "dgesvj.f"
					rotok = aaqq <= aapp / small;
#line 1200 "dgesvj.f"
				    }
#line 1201 "dgesvj.f"
				    if (aapp > small / aaqq) {
#line 1202 "dgesvj.f"
					aapq = ddot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 1205 "dgesvj.f"
				    } else {
#line 1206 "dgesvj.f"
					dcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1208 "dgesvj.f"
					dlascl_("G", &c__0, &c__0, &aaqq, &
						work[q], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1211 "dgesvj.f"
					aapq = ddot_(m, &work[*n + 1], &c__1, 
						&a[p * a_dim1 + 1], &c__1) * 
						work[p] / aapp;
#line 1213 "dgesvj.f"
				    }
#line 1214 "dgesvj.f"
				}

/* Computing MAX */
#line 1216 "dgesvj.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 1216 "dgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 1220 "dgesvj.f"
				if (abs(aapq) > tol) {
#line 1221 "dgesvj.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 1223 "dgesvj.f"
				    pskipped = 0;
#line 1224 "dgesvj.f"
				    ++iswrot;

#line 1226 "dgesvj.f"
				    if (rotok) {

#line 1228 "dgesvj.f"
					aqoap = aaqq / aapp;
#line 1229 "dgesvj.f"
					apoaq = aapp / aaqq;
#line 1230 "dgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 1231 "dgesvj.f"
					if (aaqq > aapp0) {
#line 1231 "dgesvj.f"
					    theta = -theta;
#line 1231 "dgesvj.f"
					}

#line 1233 "dgesvj.f"
					if (abs(theta) > bigtheta) {
#line 1234 "dgesvj.f"
					    t = .5 / theta;
#line 1235 "dgesvj.f"
					    fastr[2] = t * work[p] / work[q];
#line 1236 "dgesvj.f"
					    fastr[3] = -t * work[q] / work[p];
#line 1238 "dgesvj.f"
					    drotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 1240 "dgesvj.f"
					    if (rsvec) {
#line 1240 "dgesvj.f"
			  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 1240 "dgesvj.f"
					    }
/* Computing MAX */
#line 1244 "dgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 1244 "dgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1246 "dgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 1246 "dgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 1248 "dgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 1248 "dgesvj.f"
					    mxsinj = max(d__1,d__2);
#line 1249 "dgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 1253 "dgesvj.f"
					    thsign = -d_sign(&c_b18, &aapq);
#line 1254 "dgesvj.f"
					    if (aaqq > aapp0) {
#line 1254 "dgesvj.f"
			  thsign = -thsign;
#line 1254 "dgesvj.f"
					    }
#line 1255 "dgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 1257 "dgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 1258 "dgesvj.f"
					    sn = t * cs;
/* Computing MAX */
#line 1259 "dgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 1259 "dgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 1260 "dgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 1260 "dgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1262 "dgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 1262 "dgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 1265 "dgesvj.f"
					    apoaq = work[p] / work[q];
#line 1266 "dgesvj.f"
					    aqoap = work[q] / work[p];
#line 1267 "dgesvj.f"
					    if (work[p] >= 1.) {

#line 1269 "dgesvj.f"
			  if (work[q] >= 1.) {
#line 1270 "dgesvj.f"
			      fastr[2] = t * apoaq;
#line 1271 "dgesvj.f"
			      fastr[3] = -t * aqoap;
#line 1272 "dgesvj.f"
			      work[p] *= cs;
#line 1273 "dgesvj.f"
			      work[q] *= cs;
#line 1274 "dgesvj.f"
			      drotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 1277 "dgesvj.f"
			      if (rsvec) {
#line 1277 "dgesvj.f"
				  drotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 1277 "dgesvj.f"
			      }
#line 1280 "dgesvj.f"
			  } else {
#line 1281 "dgesvj.f"
			      d__1 = -t * aqoap;
#line 1281 "dgesvj.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1284 "dgesvj.f"
			      d__1 = cs * sn * apoaq;
#line 1284 "dgesvj.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 1287 "dgesvj.f"
			      if (rsvec) {
#line 1288 "dgesvj.f"
				  d__1 = -t * aqoap;
#line 1288 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1291 "dgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1291 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1295 "dgesvj.f"
			      }
#line 1296 "dgesvj.f"
			      work[p] *= cs;
#line 1297 "dgesvj.f"
			      work[q] /= cs;
#line 1298 "dgesvj.f"
			  }
#line 1299 "dgesvj.f"
					    } else {
#line 1300 "dgesvj.f"
			  if (work[q] >= 1.) {
#line 1301 "dgesvj.f"
			      d__1 = t * apoaq;
#line 1301 "dgesvj.f"
			      daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 1304 "dgesvj.f"
			      d__1 = -cs * sn * aqoap;
#line 1304 "dgesvj.f"
			      daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1307 "dgesvj.f"
			      if (rsvec) {
#line 1308 "dgesvj.f"
				  d__1 = t * apoaq;
#line 1308 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1311 "dgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1311 "dgesvj.f"
				  daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1315 "dgesvj.f"
			      }
#line 1316 "dgesvj.f"
			      work[p] /= cs;
#line 1317 "dgesvj.f"
			      work[q] *= cs;
#line 1318 "dgesvj.f"
			  } else {
#line 1319 "dgesvj.f"
			      if (work[p] >= work[q]) {
#line 1321 "dgesvj.f"
				  d__1 = -t * aqoap;
#line 1321 "dgesvj.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1324 "dgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1324 "dgesvj.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1327 "dgesvj.f"
				  work[p] *= cs;
#line 1328 "dgesvj.f"
				  work[q] /= cs;
#line 1329 "dgesvj.f"
				  if (rsvec) {
#line 1330 "dgesvj.f"
				      d__1 = -t * aqoap;
#line 1330 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1334 "dgesvj.f"
				      d__1 = cs * sn * apoaq;
#line 1334 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1338 "dgesvj.f"
				  }
#line 1339 "dgesvj.f"
			      } else {
#line 1340 "dgesvj.f"
				  d__1 = t * apoaq;
#line 1340 "dgesvj.f"
				  daxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1343 "dgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1343 "dgesvj.f"
				  daxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1347 "dgesvj.f"
				  work[p] /= cs;
#line 1348 "dgesvj.f"
				  work[q] *= cs;
#line 1349 "dgesvj.f"
				  if (rsvec) {
#line 1350 "dgesvj.f"
				      d__1 = t * apoaq;
#line 1350 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1353 "dgesvj.f"
				      d__1 = -cs * sn * aqoap;
#line 1353 "dgesvj.f"
				      daxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1357 "dgesvj.f"
				  }
#line 1358 "dgesvj.f"
			      }
#line 1359 "dgesvj.f"
			  }
#line 1360 "dgesvj.f"
					    }
#line 1361 "dgesvj.f"
					}

#line 1363 "dgesvj.f"
				    } else {
#line 1364 "dgesvj.f"
					if (aapp > aaqq) {
#line 1365 "dgesvj.f"
					    dcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[*n + 1], &
						    c__1);
#line 1367 "dgesvj.f"
					    dlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &work[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1370 "dgesvj.f"
					    dlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1373 "dgesvj.f"
					    temp1 = -aapq * work[p] / work[q];
#line 1374 "dgesvj.f"
					    daxpy_(m, &temp1, &work[*n + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1);
#line 1376 "dgesvj.f"
					    dlascl_("G", &c__0, &c__0, &c_b18,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1379 "dgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 1379 "dgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 1381 "dgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1382 "dgesvj.f"
					} else {
#line 1383 "dgesvj.f"
					    dcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[*n + 1], &
						    c__1);
#line 1385 "dgesvj.f"
					    dlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &work[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1388 "dgesvj.f"
					    dlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1391 "dgesvj.f"
					    temp1 = -aapq * work[q] / work[p];
#line 1392 "dgesvj.f"
					    daxpy_(m, &temp1, &work[*n + 1], &
						    c__1, &a[p * a_dim1 + 1], 
						    &c__1);
#line 1394 "dgesvj.f"
					    dlascl_("G", &c__0, &c__0, &c_b18,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1397 "dgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 1397 "dgesvj.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 1399 "dgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1400 "dgesvj.f"
					}
#line 1401 "dgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 1406 "dgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1406 "dgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1408 "dgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1410 "dgesvj.f"
					    sva[q] = dnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * work[q];
#line 1412 "dgesvj.f"
					} else {
#line 1413 "dgesvj.f"
					    t = 0.;
#line 1414 "dgesvj.f"
					    aaqq = 1.;
#line 1415 "dgesvj.f"
					    dlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1417 "dgesvj.f"
					    sva[q] = t * sqrt(aaqq) * work[q];
#line 1418 "dgesvj.f"
					}
#line 1419 "dgesvj.f"
				    }
/* Computing 2nd power */
#line 1420 "dgesvj.f"
				    d__1 = aapp / aapp0;
#line 1420 "dgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1421 "dgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1423 "dgesvj.f"
					    aapp = dnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * work[p];
#line 1425 "dgesvj.f"
					} else {
#line 1426 "dgesvj.f"
					    t = 0.;
#line 1427 "dgesvj.f"
					    aapp = 1.;
#line 1428 "dgesvj.f"
					    dlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1430 "dgesvj.f"
					    aapp = t * sqrt(aapp) * work[p];
#line 1431 "dgesvj.f"
					}
#line 1432 "dgesvj.f"
					sva[p] = aapp;
#line 1433 "dgesvj.f"
				    }
/*              end of OK rotation */
#line 1435 "dgesvj.f"
				} else {
#line 1436 "dgesvj.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1438 "dgesvj.f"
				    ++pskipped;
#line 1439 "dgesvj.f"
				    ++ijblsk;
#line 1440 "dgesvj.f"
				}
#line 1441 "dgesvj.f"
			    } else {
#line 1442 "dgesvj.f"
				++notrot;
#line 1443 "dgesvj.f"
				++pskipped;
#line 1444 "dgesvj.f"
				++ijblsk;
#line 1445 "dgesvj.f"
			    }

#line 1447 "dgesvj.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 1449 "dgesvj.f"
				sva[p] = aapp;
#line 1450 "dgesvj.f"
				notrot = 0;
#line 1451 "dgesvj.f"
				goto L2011;
#line 1452 "dgesvj.f"
			    }
#line 1453 "dgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1455 "dgesvj.f"
				aapp = -aapp;
#line 1456 "dgesvj.f"
				notrot = 0;
#line 1457 "dgesvj.f"
				goto L2203;
#line 1458 "dgesvj.f"
			    }

#line 1460 "dgesvj.f"
/* L2200: */
#line 1460 "dgesvj.f"
			}
/*        end of the q-loop */
#line 1462 "dgesvj.f"
L2203:

#line 1464 "dgesvj.f"
			sva[p] = aapp;

#line 1466 "dgesvj.f"
		    } else {

#line 1468 "dgesvj.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1468 "dgesvj.f"
			    i__4 = jgl + kbl - 1;
#line 1468 "dgesvj.f"
			    notrot = notrot + min(i__4,*n) - jgl + 1;
#line 1468 "dgesvj.f"
			}
#line 1470 "dgesvj.f"
			if (aapp < 0.) {
#line 1470 "dgesvj.f"
			    notrot = 0;
#line 1470 "dgesvj.f"
			}

#line 1472 "dgesvj.f"
		    }

#line 1474 "dgesvj.f"
/* L2100: */
#line 1474 "dgesvj.f"
		}
/*     end of the p-loop */
#line 1476 "dgesvj.f"
/* L2010: */
#line 1476 "dgesvj.f"
	    }
/*     end of the jbc-loop */
#line 1478 "dgesvj.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1480 "dgesvj.f"
	    i__3 = igl + kbl - 1;
#line 1480 "dgesvj.f"
	    i__2 = min(i__3,*n);
#line 1480 "dgesvj.f"
	    for (p = igl; p <= i__2; ++p) {
#line 1481 "dgesvj.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1482 "dgesvj.f"
/* L2012: */
#line 1482 "dgesvj.f"
	    }
/* ** */
#line 1484 "dgesvj.f"
/* L2000: */
#line 1484 "dgesvj.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1488 "dgesvj.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1490 "dgesvj.f"
	    sva[*n] = dnrm2_(m, &a[*n * a_dim1 + 1], &c__1) * work[*n];
#line 1491 "dgesvj.f"
	} else {
#line 1492 "dgesvj.f"
	    t = 0.;
#line 1493 "dgesvj.f"
	    aapp = 1.;
#line 1494 "dgesvj.f"
	    dlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1495 "dgesvj.f"
	    sva[*n] = t * sqrt(aapp) * work[*n];
#line 1496 "dgesvj.f"
	}

/*     Additional steering devices */

#line 1500 "dgesvj.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1500 "dgesvj.f"
	    swband = i__;
#line 1500 "dgesvj.f"
	}

#line 1503 "dgesvj.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * tol && (
		doublereal) (*n) * mxaapq * mxsinj < tol) {
#line 1505 "dgesvj.f"
	    goto L1994;
#line 1506 "dgesvj.f"
	}

#line 1508 "dgesvj.f"
	if (notrot >= emptsw) {
#line 1508 "dgesvj.f"
	    goto L1994;
#line 1508 "dgesvj.f"
	}

#line 1510 "dgesvj.f"
/* L1993: */
#line 1510 "dgesvj.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 1514 "dgesvj.f"
    *info = 29;
#line 1515 "dgesvj.f"
    goto L1995;

#line 1517 "dgesvj.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 1521 "dgesvj.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1523 "dgesvj.f"
L1995:

/*     Sort the singular values and find how many are above */
/*     the underflow threshold. */

#line 1528 "dgesvj.f"
    n2 = 0;
#line 1529 "dgesvj.f"
    n4 = 0;
#line 1530 "dgesvj.f"
    i__1 = *n - 1;
#line 1530 "dgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 1531 "dgesvj.f"
	i__2 = *n - p + 1;
#line 1531 "dgesvj.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1532 "dgesvj.f"
	if (p != q) {
#line 1533 "dgesvj.f"
	    temp1 = sva[p];
#line 1534 "dgesvj.f"
	    sva[p] = sva[q];
#line 1535 "dgesvj.f"
	    sva[q] = temp1;
#line 1536 "dgesvj.f"
	    temp1 = work[p];
#line 1537 "dgesvj.f"
	    work[p] = work[q];
#line 1538 "dgesvj.f"
	    work[q] = temp1;
#line 1539 "dgesvj.f"
	    dswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1540 "dgesvj.f"
	    if (rsvec) {
#line 1540 "dgesvj.f"
		dswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1540 "dgesvj.f"
	    }
#line 1541 "dgesvj.f"
	}
#line 1542 "dgesvj.f"
	if (sva[p] != 0.) {
#line 1543 "dgesvj.f"
	    ++n4;
#line 1544 "dgesvj.f"
	    if (sva[p] * skl > sfmin) {
#line 1544 "dgesvj.f"
		++n2;
#line 1544 "dgesvj.f"
	    }
#line 1545 "dgesvj.f"
	}
#line 1546 "dgesvj.f"
/* L5991: */
#line 1546 "dgesvj.f"
    }
#line 1547 "dgesvj.f"
    if (sva[*n] != 0.) {
#line 1548 "dgesvj.f"
	++n4;
#line 1549 "dgesvj.f"
	if (sva[*n] * skl > sfmin) {
#line 1549 "dgesvj.f"
	    ++n2;
#line 1549 "dgesvj.f"
	}
#line 1550 "dgesvj.f"
    }

/*     Normalize the left singular vectors. */

#line 1554 "dgesvj.f"
    if (lsvec || uctol) {
#line 1555 "dgesvj.f"
	i__1 = n2;
#line 1555 "dgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1556 "dgesvj.f"
	    d__1 = work[p] / sva[p];
#line 1556 "dgesvj.f"
	    dscal_(m, &d__1, &a[p * a_dim1 + 1], &c__1);
#line 1557 "dgesvj.f"
/* L1998: */
#line 1557 "dgesvj.f"
	}
#line 1558 "dgesvj.f"
    }

/*     Scale the product of Jacobi rotations (assemble the fast rotations). */

#line 1562 "dgesvj.f"
    if (rsvec) {
#line 1563 "dgesvj.f"
	if (applv) {
#line 1564 "dgesvj.f"
	    i__1 = *n;
#line 1564 "dgesvj.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1565 "dgesvj.f"
		dscal_(&mvl, &work[p], &v[p * v_dim1 + 1], &c__1);
#line 1566 "dgesvj.f"
/* L2398: */
#line 1566 "dgesvj.f"
	    }
#line 1567 "dgesvj.f"
	} else {
#line 1568 "dgesvj.f"
	    i__1 = *n;
#line 1568 "dgesvj.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1569 "dgesvj.f"
		temp1 = 1. / dnrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
#line 1570 "dgesvj.f"
		dscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
#line 1571 "dgesvj.f"
/* L2399: */
#line 1571 "dgesvj.f"
	    }
#line 1572 "dgesvj.f"
	}
#line 1573 "dgesvj.f"
    }

/*     Undo scaling, if necessary (and possible). */
#line 1576 "dgesvj.f"
    if (skl > 1. && sva[1] < big / skl || skl < 1. && sva[max(n2,1)] > sfmin /
	     skl) {
#line 1579 "dgesvj.f"
	i__1 = *n;
#line 1579 "dgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1580 "dgesvj.f"
	    sva[p] = skl * sva[p];
#line 1581 "dgesvj.f"
/* L2400: */
#line 1581 "dgesvj.f"
	}
#line 1582 "dgesvj.f"
	skl = 1.;
#line 1583 "dgesvj.f"
    }

#line 1585 "dgesvj.f"
    work[1] = skl;
/*     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
/*     then some of the singular values may overflow or underflow and */
/*     the spectrum is given in this factored representation. */

#line 1590 "dgesvj.f"
    work[2] = (doublereal) n4;
/*     N4 is the number of computed nonzero singular values of A. */

#line 1593 "dgesvj.f"
    work[3] = (doublereal) n2;
/*     N2 is the number of singular values of A greater than SFMIN. */
/*     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
/*     that may carry some information. */

#line 1598 "dgesvj.f"
    work[4] = (doublereal) i__;
/*     i is the index of the last sweep before declaring convergence. */

#line 1601 "dgesvj.f"
    work[5] = mxaapq;
/*     MXAAPQ is the largest absolute value of scaled pivots in the */
/*     last sweep */

#line 1605 "dgesvj.f"
    work[6] = mxsinj;
/*     MXSINJ is the largest absolute value of the sines of Jacobi angles */
/*     in the last sweep */

#line 1609 "dgesvj.f"
    return 0;
/*     .. */
/*     .. END OF DGESVJ */
/*     .. */
} /* dgesvj_ */

