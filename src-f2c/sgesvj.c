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
/* > SGESVJ can sometimes compute tiny singular values and their singular vectors much */
/* > more accurately than other SVD routines, see below under Further Details. */
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
/* >          WORK is REAL array, dimension MAX(6,M+N). */
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

/* > \date December 2016 */

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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 393 "sgesvj.f"
    /* Parameter adjustments */
#line 393 "sgesvj.f"
    --sva;
#line 393 "sgesvj.f"
    a_dim1 = *lda;
#line 393 "sgesvj.f"
    a_offset = 1 + a_dim1;
#line 393 "sgesvj.f"
    a -= a_offset;
#line 393 "sgesvj.f"
    v_dim1 = *ldv;
#line 393 "sgesvj.f"
    v_offset = 1 + v_dim1;
#line 393 "sgesvj.f"
    v -= v_offset;
#line 393 "sgesvj.f"
    --work;
#line 393 "sgesvj.f"

#line 393 "sgesvj.f"
    /* Function Body */
#line 393 "sgesvj.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 394 "sgesvj.f"
    uctol = lsame_(jobu, "C", (ftnlen)1, (ftnlen)1);
#line 395 "sgesvj.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 396 "sgesvj.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 397 "sgesvj.f"
    upper = lsame_(joba, "U", (ftnlen)1, (ftnlen)1);
#line 398 "sgesvj.f"
    lower = lsame_(joba, "L", (ftnlen)1, (ftnlen)1);

#line 400 "sgesvj.f"
    if (! (upper || lower || lsame_(joba, "G", (ftnlen)1, (ftnlen)1))) {
#line 401 "sgesvj.f"
	*info = -1;
#line 402 "sgesvj.f"
    } else if (! (lsvec || uctol || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 403 "sgesvj.f"
	*info = -2;
#line 404 "sgesvj.f"
    } else if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 405 "sgesvj.f"
	*info = -3;
#line 406 "sgesvj.f"
    } else if (*m < 0) {
#line 407 "sgesvj.f"
	*info = -4;
#line 408 "sgesvj.f"
    } else if (*n < 0 || *n > *m) {
#line 409 "sgesvj.f"
	*info = -5;
#line 410 "sgesvj.f"
    } else if (*lda < *m) {
#line 411 "sgesvj.f"
	*info = -7;
#line 412 "sgesvj.f"
    } else if (*mv < 0) {
#line 413 "sgesvj.f"
	*info = -9;
#line 414 "sgesvj.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 416 "sgesvj.f"
	*info = -11;
#line 417 "sgesvj.f"
    } else if (uctol && work[1] <= 1.) {
#line 418 "sgesvj.f"
	*info = -12;
#line 419 "sgesvj.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 419 "sgesvj.f"
	i__1 = *m + *n;
#line 419 "sgesvj.f"
	if (*lwork < max(i__1,6)) {
#line 420 "sgesvj.f"
	    *info = -13;
#line 421 "sgesvj.f"
	} else {
#line 422 "sgesvj.f"
	    *info = 0;
#line 423 "sgesvj.f"
	}
#line 423 "sgesvj.f"
    }

/*     #:( */
#line 426 "sgesvj.f"
    if (*info != 0) {
#line 427 "sgesvj.f"
	i__1 = -(*info);
#line 427 "sgesvj.f"
	xerbla_("SGESVJ", &i__1, (ftnlen)6);
#line 428 "sgesvj.f"
	return 0;
#line 429 "sgesvj.f"
    }

/* #:) Quick return for void matrix */

#line 433 "sgesvj.f"
    if (*m == 0 || *n == 0) {
#line 433 "sgesvj.f"
	return 0;
#line 433 "sgesvj.f"
    }

/*     Set numerical parameters */
/*     The stopping criterion for Jacobi rotations is */

/*     max_{i<>j}|A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS */

/*     where EPS is the round-off and CTOL is defined as follows: */

#line 442 "sgesvj.f"
    if (uctol) {
/*        ... user controlled */
#line 444 "sgesvj.f"
	ctol = work[1];
#line 445 "sgesvj.f"
    } else {
/*        ... default */
#line 447 "sgesvj.f"
	if (lsvec || rsvec || applv) {
#line 448 "sgesvj.f"
	    ctol = sqrt((doublereal) (*m));
#line 449 "sgesvj.f"
	} else {
#line 450 "sgesvj.f"
	    ctol = (doublereal) (*m);
#line 451 "sgesvj.f"
	}
#line 452 "sgesvj.f"
    }
/*     ... and the machine dependent parameters are */
/* [!]  (Make sure that SLAMCH() works properly on the target machine.) */

#line 456 "sgesvj.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 457 "sgesvj.f"
    rooteps = sqrt(epsln);
#line 458 "sgesvj.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 459 "sgesvj.f"
    rootsfmin = sqrt(sfmin);
#line 460 "sgesvj.f"
    small = sfmin / epsln;
#line 461 "sgesvj.f"
    big = slamch_("Overflow", (ftnlen)8);
/*     BIG         = ONE    / SFMIN */
#line 463 "sgesvj.f"
    rootbig = 1. / rootsfmin;
#line 464 "sgesvj.f"
    large = big / sqrt((doublereal) (*m * *n));
#line 465 "sgesvj.f"
    bigtheta = 1. / rooteps;

#line 467 "sgesvj.f"
    tol = ctol * epsln;
#line 468 "sgesvj.f"
    roottol = sqrt(tol);

#line 470 "sgesvj.f"
    if ((doublereal) (*m) * epsln >= 1.) {
#line 471 "sgesvj.f"
	*info = -4;
#line 472 "sgesvj.f"
	i__1 = -(*info);
#line 472 "sgesvj.f"
	xerbla_("SGESVJ", &i__1, (ftnlen)6);
#line 473 "sgesvj.f"
	return 0;
#line 474 "sgesvj.f"
    }

/*     Initialize the right singular vector matrix. */

#line 478 "sgesvj.f"
    if (rsvec) {
#line 479 "sgesvj.f"
	mvl = *n;
#line 480 "sgesvj.f"
	slaset_("A", &mvl, n, &c_b17, &c_b18, &v[v_offset], ldv, (ftnlen)1);
#line 481 "sgesvj.f"
    } else if (applv) {
#line 482 "sgesvj.f"
	mvl = *mv;
#line 483 "sgesvj.f"
    }
#line 484 "sgesvj.f"
    rsvec = rsvec || applv;

/*     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
/* (!)  If necessary, scale A to protect the largest singular value */
/*     from overflow. It is possible that saving the largest singular */
/*     value destroys the information about the small ones. */
/*     This initial scaling is almost minimal in the sense that the */
/*     goal is to make sure that no column norm overflows, and that */
/*     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
/*     in A are detected, the procedure returns with INFO=-6. */

#line 495 "sgesvj.f"
    skl = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 496 "sgesvj.f"
    noscale = TRUE_;
#line 497 "sgesvj.f"
    goscale = TRUE_;

#line 499 "sgesvj.f"
    if (lower) {
/*        the input matrix is M-by-N lower triangular (trapezoidal) */
#line 501 "sgesvj.f"
	i__1 = *n;
#line 501 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 502 "sgesvj.f"
	    aapp = 0.;
#line 503 "sgesvj.f"
	    aaqq = 1.;
#line 504 "sgesvj.f"
	    i__2 = *m - p + 1;
#line 504 "sgesvj.f"
	    slassq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
#line 505 "sgesvj.f"
	    if (aapp > big) {
#line 506 "sgesvj.f"
		*info = -6;
#line 507 "sgesvj.f"
		i__2 = -(*info);
#line 507 "sgesvj.f"
		xerbla_("SGESVJ", &i__2, (ftnlen)6);
#line 508 "sgesvj.f"
		return 0;
#line 509 "sgesvj.f"
	    }
#line 510 "sgesvj.f"
	    aaqq = sqrt(aaqq);
#line 511 "sgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 512 "sgesvj.f"
		sva[p] = aapp * aaqq;
#line 513 "sgesvj.f"
	    } else {
#line 514 "sgesvj.f"
		noscale = FALSE_;
#line 515 "sgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 516 "sgesvj.f"
		if (goscale) {
#line 517 "sgesvj.f"
		    goscale = FALSE_;
#line 518 "sgesvj.f"
		    i__2 = p - 1;
#line 518 "sgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 519 "sgesvj.f"
			sva[q] *= skl;
#line 520 "sgesvj.f"
/* L1873: */
#line 520 "sgesvj.f"
		    }
#line 521 "sgesvj.f"
		}
#line 522 "sgesvj.f"
	    }
#line 523 "sgesvj.f"
/* L1874: */
#line 523 "sgesvj.f"
	}
#line 524 "sgesvj.f"
    } else if (upper) {
/*        the input matrix is M-by-N upper triangular (trapezoidal) */
#line 526 "sgesvj.f"
	i__1 = *n;
#line 526 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 527 "sgesvj.f"
	    aapp = 0.;
#line 528 "sgesvj.f"
	    aaqq = 1.;
#line 529 "sgesvj.f"
	    slassq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 530 "sgesvj.f"
	    if (aapp > big) {
#line 531 "sgesvj.f"
		*info = -6;
#line 532 "sgesvj.f"
		i__2 = -(*info);
#line 532 "sgesvj.f"
		xerbla_("SGESVJ", &i__2, (ftnlen)6);
#line 533 "sgesvj.f"
		return 0;
#line 534 "sgesvj.f"
	    }
#line 535 "sgesvj.f"
	    aaqq = sqrt(aaqq);
#line 536 "sgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 537 "sgesvj.f"
		sva[p] = aapp * aaqq;
#line 538 "sgesvj.f"
	    } else {
#line 539 "sgesvj.f"
		noscale = FALSE_;
#line 540 "sgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 541 "sgesvj.f"
		if (goscale) {
#line 542 "sgesvj.f"
		    goscale = FALSE_;
#line 543 "sgesvj.f"
		    i__2 = p - 1;
#line 543 "sgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 544 "sgesvj.f"
			sva[q] *= skl;
#line 545 "sgesvj.f"
/* L2873: */
#line 545 "sgesvj.f"
		    }
#line 546 "sgesvj.f"
		}
#line 547 "sgesvj.f"
	    }
#line 548 "sgesvj.f"
/* L2874: */
#line 548 "sgesvj.f"
	}
#line 549 "sgesvj.f"
    } else {
/*        the input matrix is M-by-N general dense */
#line 551 "sgesvj.f"
	i__1 = *n;
#line 551 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 552 "sgesvj.f"
	    aapp = 0.;
#line 553 "sgesvj.f"
	    aaqq = 1.;
#line 554 "sgesvj.f"
	    slassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 555 "sgesvj.f"
	    if (aapp > big) {
#line 556 "sgesvj.f"
		*info = -6;
#line 557 "sgesvj.f"
		i__2 = -(*info);
#line 557 "sgesvj.f"
		xerbla_("SGESVJ", &i__2, (ftnlen)6);
#line 558 "sgesvj.f"
		return 0;
#line 559 "sgesvj.f"
	    }
#line 560 "sgesvj.f"
	    aaqq = sqrt(aaqq);
#line 561 "sgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 562 "sgesvj.f"
		sva[p] = aapp * aaqq;
#line 563 "sgesvj.f"
	    } else {
#line 564 "sgesvj.f"
		noscale = FALSE_;
#line 565 "sgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 566 "sgesvj.f"
		if (goscale) {
#line 567 "sgesvj.f"
		    goscale = FALSE_;
#line 568 "sgesvj.f"
		    i__2 = p - 1;
#line 568 "sgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 569 "sgesvj.f"
			sva[q] *= skl;
#line 570 "sgesvj.f"
/* L3873: */
#line 570 "sgesvj.f"
		    }
#line 571 "sgesvj.f"
		}
#line 572 "sgesvj.f"
	    }
#line 573 "sgesvj.f"
/* L3874: */
#line 573 "sgesvj.f"
	}
#line 574 "sgesvj.f"
    }

#line 576 "sgesvj.f"
    if (noscale) {
#line 576 "sgesvj.f"
	skl = 1.;
#line 576 "sgesvj.f"
    }

/*     Move the smaller part of the spectrum from the underflow threshold */
/* (!)  Start by determining the position of the nonzero entries of the */
/*     array SVA() relative to ( SFMIN, BIG ). */

#line 582 "sgesvj.f"
    aapp = 0.;
#line 583 "sgesvj.f"
    aaqq = big;
#line 584 "sgesvj.f"
    i__1 = *n;
#line 584 "sgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 585 "sgesvj.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 585 "sgesvj.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 585 "sgesvj.f"
	    aaqq = min(d__1,d__2);
#line 585 "sgesvj.f"
	}
/* Computing MAX */
#line 586 "sgesvj.f"
	d__1 = aapp, d__2 = sva[p];
#line 586 "sgesvj.f"
	aapp = max(d__1,d__2);
#line 587 "sgesvj.f"
/* L4781: */
#line 587 "sgesvj.f"
    }

/* #:) Quick return for zero matrix */

#line 591 "sgesvj.f"
    if (aapp == 0.) {
#line 592 "sgesvj.f"
	if (lsvec) {
#line 592 "sgesvj.f"
	    slaset_("G", m, n, &c_b17, &c_b18, &a[a_offset], lda, (ftnlen)1);
#line 592 "sgesvj.f"
	}
#line 593 "sgesvj.f"
	work[1] = 1.;
#line 594 "sgesvj.f"
	work[2] = 0.;
#line 595 "sgesvj.f"
	work[3] = 0.;
#line 596 "sgesvj.f"
	work[4] = 0.;
#line 597 "sgesvj.f"
	work[5] = 0.;
#line 598 "sgesvj.f"
	work[6] = 0.;
#line 599 "sgesvj.f"
	return 0;
#line 600 "sgesvj.f"
    }

/* #:) Quick return for one-column matrix */

#line 604 "sgesvj.f"
    if (*n == 1) {
#line 605 "sgesvj.f"
	if (lsvec) {
#line 605 "sgesvj.f"
	    slascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1]
		    , lda, &ierr, (ftnlen)1);
#line 605 "sgesvj.f"
	}
#line 607 "sgesvj.f"
	work[1] = 1. / skl;
#line 608 "sgesvj.f"
	if (sva[1] >= sfmin) {
#line 609 "sgesvj.f"
	    work[2] = 1.;
#line 610 "sgesvj.f"
	} else {
#line 611 "sgesvj.f"
	    work[2] = 0.;
#line 612 "sgesvj.f"
	}
#line 613 "sgesvj.f"
	work[3] = 0.;
#line 614 "sgesvj.f"
	work[4] = 0.;
#line 615 "sgesvj.f"
	work[5] = 0.;
#line 616 "sgesvj.f"
	work[6] = 0.;
#line 617 "sgesvj.f"
	return 0;
#line 618 "sgesvj.f"
    }

/*     Protect small singular values from underflow, and try to */
/*     avoid underflows/overflows in computing Jacobi rotations. */

#line 623 "sgesvj.f"
    sn = sqrt(sfmin / epsln);
#line 624 "sgesvj.f"
    temp1 = sqrt(big / (doublereal) (*n));
#line 625 "sgesvj.f"
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1) {
/* Computing MIN */
#line 627 "sgesvj.f"
	d__1 = big, d__2 = temp1 / aapp;
#line 627 "sgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 630 "sgesvj.f"
    } else if (aaqq <= sn && aapp <= temp1) {
/* Computing MIN */
#line 631 "sgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (aapp * sqrt((doublereal) (*n)));
#line 631 "sgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 634 "sgesvj.f"
    } else if (aaqq >= sn && aapp >= temp1) {
/* Computing MAX */
#line 635 "sgesvj.f"
	d__1 = sn / aaqq, d__2 = temp1 / aapp;
#line 635 "sgesvj.f"
	temp1 = max(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 638 "sgesvj.f"
    } else if (aaqq <= sn && aapp >= temp1) {
/* Computing MIN */
#line 639 "sgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (sqrt((doublereal) (*n)) * aapp);
#line 639 "sgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 642 "sgesvj.f"
    } else {
#line 643 "sgesvj.f"
	temp1 = 1.;
#line 644 "sgesvj.f"
    }

/*     Scale, if necessary */

#line 648 "sgesvj.f"
    if (temp1 != 1.) {
#line 649 "sgesvj.f"
	slascl_("G", &c__0, &c__0, &c_b18, &temp1, n, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 650 "sgesvj.f"
    }
#line 651 "sgesvj.f"
    skl = temp1 * skl;
#line 652 "sgesvj.f"
    if (skl != 1.) {
#line 653 "sgesvj.f"
	slascl_(joba, &c__0, &c__0, &c_b18, &skl, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 654 "sgesvj.f"
	skl = 1. / skl;
#line 655 "sgesvj.f"
    }

/*     Row-cyclic Jacobi SVD algorithm with column pivoting */

#line 659 "sgesvj.f"
    emptsw = *n * (*n - 1) / 2;
#line 660 "sgesvj.f"
    notrot = 0;
#line 661 "sgesvj.f"
    fastr[0] = 0.;

/*     A is represented in factored form A = A * diag(WORK), where diag(WORK) */
/*     is initialized to identity. WORK is updated during fast scaled */
/*     rotations. */

#line 667 "sgesvj.f"
    i__1 = *n;
#line 667 "sgesvj.f"
    for (q = 1; q <= i__1; ++q) {
#line 668 "sgesvj.f"
	work[q] = 1.;
#line 669 "sgesvj.f"
/* L1868: */
#line 669 "sgesvj.f"
    }


#line 672 "sgesvj.f"
    swband = 3;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if SGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 680 "sgesvj.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 686 "sgesvj.f"
    nbl = *n / kbl;
#line 687 "sgesvj.f"
    if (nbl * kbl != *n) {
#line 687 "sgesvj.f"
	++nbl;
#line 687 "sgesvj.f"
    }

/* Computing 2nd power */
#line 689 "sgesvj.f"
    i__1 = kbl;
#line 689 "sgesvj.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 692 "sgesvj.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 695 "sgesvj.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */

/* Computing MAX */
#line 703 "sgesvj.f"
    i__1 = 64, i__2 = kbl << 2;
#line 703 "sgesvj.f"
    if ((lower || upper) && *n > max(i__1,i__2)) {
/* [TP] The number of partition levels and the actual partition are */
/*     tuning parameters. */
#line 706 "sgesvj.f"
	n4 = *n / 4;
#line 707 "sgesvj.f"
	n2 = *n / 2;
#line 708 "sgesvj.f"
	n34 = n4 * 3;
#line 709 "sgesvj.f"
	if (applv) {
#line 710 "sgesvj.f"
	    q = 0;
#line 711 "sgesvj.f"
	} else {
#line 712 "sgesvj.f"
	    q = 1;
#line 713 "sgesvj.f"
	}

#line 715 "sgesvj.f"
	if (lower) {

/*     This works very well on lower triangular matrices, in particular */
/*     in the framework of the preconditioned Jacobi SVD (xGEJSV). */
/*     The idea is simple: */
/*     [+ 0 0 0]   Note that Jacobi transformations of [0 0] */
/*     [+ + 0 0]                                       [0 0] */
/*     [+ + x 0]   actually work on [x 0]              [x 0] */
/*     [+ + x x]                    [x x].             [x x] */

#line 725 "sgesvj.f"
	    i__1 = *m - n34;
#line 725 "sgesvj.f"
	    i__2 = *n - n34;
#line 725 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 725 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda,
		     &work[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + (
		    n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &
		    work[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 730 "sgesvj.f"
	    i__1 = *m - n2;
#line 730 "sgesvj.f"
	    i__2 = n34 - n2;
#line 730 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 730 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &
		    work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1)
		     * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &work[*n + 
		    1], &i__3, &ierr, (ftnlen)1);

#line 735 "sgesvj.f"
	    i__1 = *m - n2;
#line 735 "sgesvj.f"
	    i__2 = *n - n2;
#line 735 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 735 "sgesvj.f"
	    sgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], 
		    lda, &work[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (
		    n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 740 "sgesvj.f"
	    i__1 = *m - n4;
#line 740 "sgesvj.f"
	    i__2 = n2 - n4;
#line 740 "sgesvj.f"
	    i__3 = *lwork - *n;
#line 740 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, &
		    work[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1)
		     * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 
		    1], &i__3, &ierr, (ftnlen)1);

#line 745 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 745 "sgesvj.f"
	    sgsvj0_(jobv, m, &n4, &a[a_offset], lda, &work[1], &sva[1], &mvl, 
		    &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n 
		    + 1], &i__1, &ierr, (ftnlen)1);

#line 749 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 749 "sgesvj.f"
	    sgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);


#line 754 "sgesvj.f"
	} else if (upper) {


#line 757 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 757 "sgesvj.f"
	    sgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &work[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 761 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 761 "sgesvj.f"
	    sgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &work[n4 
		    + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], 
		    &i__1, &ierr, (ftnlen)1);

#line 766 "sgesvj.f"
	    i__1 = *lwork - *n;
#line 766 "sgesvj.f"
	    sgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &work[1], &sva[1],
		     &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    work[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 770 "sgesvj.f"
	    i__1 = n2 + n4;
#line 770 "sgesvj.f"
	    i__2 = *lwork - *n;
#line 770 "sgesvj.f"
	    sgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &work[
		    n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &work[*n + 1], 
		    &i__2, &ierr, (ftnlen)1);
#line 775 "sgesvj.f"
	}

#line 777 "sgesvj.f"
    }

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 781 "sgesvj.f"
    for (i__ = 1; i__ <= 30; ++i__) {

/*     .. go go go ... */

#line 785 "sgesvj.f"
	mxaapq = 0.;
#line 786 "sgesvj.f"
	mxsinj = 0.;
#line 787 "sgesvj.f"
	iswrot = 0;

#line 789 "sgesvj.f"
	notrot = 0;
#line 790 "sgesvj.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 797 "sgesvj.f"
	i__1 = nbl;
#line 797 "sgesvj.f"
	for (ibr = 1; ibr <= i__1; ++ibr) {

#line 799 "sgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 801 "sgesvj.f"
	    i__3 = lkahead, i__4 = nbl - ibr;
#line 801 "sgesvj.f"
	    i__2 = min(i__3,i__4);
#line 801 "sgesvj.f"
	    for (ir1 = 0; ir1 <= i__2; ++ir1) {

#line 803 "sgesvj.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 805 "sgesvj.f"
		i__4 = igl + kbl - 1, i__5 = *n - 1;
#line 805 "sgesvj.f"
		i__3 = min(i__4,i__5);
#line 805 "sgesvj.f"
		for (p = igl; p <= i__3; ++p) {

/*     .. de Rijk's pivoting */

#line 809 "sgesvj.f"
		    i__4 = *n - p + 1;
#line 809 "sgesvj.f"
		    q = isamax_(&i__4, &sva[p], &c__1) + p - 1;
#line 810 "sgesvj.f"
		    if (p != q) {
#line 811 "sgesvj.f"
			sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 812 "sgesvj.f"
			if (rsvec) {
#line 812 "sgesvj.f"
			    sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 812 "sgesvj.f"
			}
#line 814 "sgesvj.f"
			temp1 = sva[p];
#line 815 "sgesvj.f"
			sva[p] = sva[q];
#line 816 "sgesvj.f"
			sva[q] = temp1;
#line 817 "sgesvj.f"
			temp1 = work[p];
#line 818 "sgesvj.f"
			work[p] = work[q];
#line 819 "sgesvj.f"
			work[q] = temp1;
#line 820 "sgesvj.f"
		    }

#line 822 "sgesvj.f"
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

#line 836 "sgesvj.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 838 "sgesvj.f"
			    sva[p] = snrm2_(m, &a[p * a_dim1 + 1], &c__1) * 
				    work[p];
#line 839 "sgesvj.f"
			} else {
#line 840 "sgesvj.f"
			    temp1 = 0.;
#line 841 "sgesvj.f"
			    aapp = 1.;
#line 842 "sgesvj.f"
			    slassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 843 "sgesvj.f"
			    sva[p] = temp1 * sqrt(aapp) * work[p];
#line 844 "sgesvj.f"
			}
#line 845 "sgesvj.f"
			aapp = sva[p];
#line 846 "sgesvj.f"
		    } else {
#line 847 "sgesvj.f"
			aapp = sva[p];
#line 848 "sgesvj.f"
		    }

#line 850 "sgesvj.f"
		    if (aapp > 0.) {

#line 852 "sgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 854 "sgesvj.f"
			i__5 = igl + kbl - 1;
#line 854 "sgesvj.f"
			i__4 = min(i__5,*n);
#line 854 "sgesvj.f"
			for (q = p + 1; q <= i__4; ++q) {

#line 856 "sgesvj.f"
			    aaqq = sva[q];

#line 858 "sgesvj.f"
			    if (aaqq > 0.) {

#line 860 "sgesvj.f"
				aapp0 = aapp;
#line 861 "sgesvj.f"
				if (aaqq >= 1.) {
#line 862 "sgesvj.f"
				    rotok = small * aapp <= aaqq;
#line 863 "sgesvj.f"
				    if (aapp < big / aaqq) {
#line 864 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 867 "sgesvj.f"
				    } else {
#line 868 "sgesvj.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 870 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						work[p], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 873 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[q * a_dim1 + 1], &c__1) * 
						work[q] / aaqq;
#line 875 "sgesvj.f"
				    }
#line 876 "sgesvj.f"
				} else {
#line 877 "sgesvj.f"
				    rotok = aapp <= aaqq / small;
#line 878 "sgesvj.f"
				    if (aapp > small / aaqq) {
#line 879 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 882 "sgesvj.f"
				    } else {
#line 883 "sgesvj.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 885 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						work[q], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 888 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[p * a_dim1 + 1], &c__1) * 
						work[p] / aapp;
#line 890 "sgesvj.f"
				    }
#line 891 "sgesvj.f"
				}

/* Computing MAX */
#line 893 "sgesvj.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 893 "sgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 897 "sgesvj.f"
				if (abs(aapq) > tol) {

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 902 "sgesvj.f"
				    if (ir1 == 0) {
#line 903 "sgesvj.f"
					notrot = 0;
#line 904 "sgesvj.f"
					pskipped = 0;
#line 905 "sgesvj.f"
					++iswrot;
#line 906 "sgesvj.f"
				    }

#line 908 "sgesvj.f"
				    if (rotok) {

#line 910 "sgesvj.f"
					aqoap = aaqq / aapp;
#line 911 "sgesvj.f"
					apoaq = aapp / aaqq;
#line 912 "sgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;

#line 914 "sgesvj.f"
					if (abs(theta) > bigtheta) {

#line 916 "sgesvj.f"
					    t = .5 / theta;
#line 917 "sgesvj.f"
					    fastr[2] = t * work[p] / work[q];
#line 918 "sgesvj.f"
					    fastr[3] = -t * work[q] / work[p];
#line 920 "sgesvj.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 922 "sgesvj.f"
					    if (rsvec) {
#line 922 "sgesvj.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 922 "sgesvj.f"
					    }
/* Computing MAX */
#line 926 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 926 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 928 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 928 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 930 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 930 "sgesvj.f"
					    mxsinj = max(d__1,d__2);

#line 932 "sgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 936 "sgesvj.f"
					    thsign = -d_sign(&c_b18, &aapq);
#line 937 "sgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 939 "sgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 940 "sgesvj.f"
					    sn = t * cs;

/* Computing MAX */
#line 942 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 942 "sgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 943 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 943 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 945 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 945 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 948 "sgesvj.f"
					    apoaq = work[p] / work[q];
#line 949 "sgesvj.f"
					    aqoap = work[q] / work[p];
#line 950 "sgesvj.f"
					    if (work[p] >= 1.) {
#line 951 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 952 "sgesvj.f"
			      fastr[2] = t * apoaq;
#line 953 "sgesvj.f"
			      fastr[3] = -t * aqoap;
#line 954 "sgesvj.f"
			      work[p] *= cs;
#line 955 "sgesvj.f"
			      work[q] *= cs;
#line 956 "sgesvj.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 959 "sgesvj.f"
			      if (rsvec) {
#line 959 "sgesvj.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 959 "sgesvj.f"
			      }
#line 962 "sgesvj.f"
			  } else {
#line 963 "sgesvj.f"
			      d__1 = -t * aqoap;
#line 963 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 966 "sgesvj.f"
			      d__1 = cs * sn * apoaq;
#line 966 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 969 "sgesvj.f"
			      work[p] *= cs;
#line 970 "sgesvj.f"
			      work[q] /= cs;
#line 971 "sgesvj.f"
			      if (rsvec) {
#line 972 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 972 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 975 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 975 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 979 "sgesvj.f"
			      }
#line 980 "sgesvj.f"
			  }
#line 981 "sgesvj.f"
					    } else {
#line 982 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 983 "sgesvj.f"
			      d__1 = t * apoaq;
#line 983 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 986 "sgesvj.f"
			      d__1 = -cs * sn * aqoap;
#line 986 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 989 "sgesvj.f"
			      work[p] /= cs;
#line 990 "sgesvj.f"
			      work[q] *= cs;
#line 991 "sgesvj.f"
			      if (rsvec) {
#line 992 "sgesvj.f"
				  d__1 = t * apoaq;
#line 992 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 995 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 995 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 999 "sgesvj.f"
			      }
#line 1000 "sgesvj.f"
			  } else {
#line 1001 "sgesvj.f"
			      if (work[p] >= work[q]) {
#line 1003 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 1003 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1006 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1006 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1009 "sgesvj.f"
				  work[p] *= cs;
#line 1010 "sgesvj.f"
				  work[q] /= cs;
#line 1011 "sgesvj.f"
				  if (rsvec) {
#line 1012 "sgesvj.f"
				      d__1 = -t * aqoap;
#line 1012 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1016 "sgesvj.f"
				      d__1 = cs * sn * apoaq;
#line 1016 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1020 "sgesvj.f"
				  }
#line 1021 "sgesvj.f"
			      } else {
#line 1022 "sgesvj.f"
				  d__1 = t * apoaq;
#line 1022 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1025 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1025 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1029 "sgesvj.f"
				  work[p] /= cs;
#line 1030 "sgesvj.f"
				  work[q] *= cs;
#line 1031 "sgesvj.f"
				  if (rsvec) {
#line 1032 "sgesvj.f"
				      d__1 = t * apoaq;
#line 1032 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1035 "sgesvj.f"
				      d__1 = -cs * sn * aqoap;
#line 1035 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1039 "sgesvj.f"
				  }
#line 1040 "sgesvj.f"
			      }
#line 1041 "sgesvj.f"
			  }
#line 1042 "sgesvj.f"
					    }
#line 1043 "sgesvj.f"
					}

#line 1045 "sgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 1047 "sgesvj.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1049 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						c_b18, m, &c__1, &work[*n + 1]
						, lda, &ierr, (ftnlen)1);
#line 1052 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						c_b18, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 1054 "sgesvj.f"
					temp1 = -aapq * work[p] / work[q];
#line 1055 "sgesvj.f"
					saxpy_(m, &temp1, &work[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1057 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &c_b18, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 1059 "sgesvj.f"
					d__1 = 0., d__2 = 1. - aapq * aapq;
#line 1059 "sgesvj.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 1061 "sgesvj.f"
					mxsinj = max(mxsinj,sfmin);
#line 1062 "sgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 1068 "sgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1068 "sgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1070 "sgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1072 "sgesvj.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * work[q];
#line 1074 "sgesvj.f"
					} else {
#line 1075 "sgesvj.f"
					    t = 0.;
#line 1076 "sgesvj.f"
					    aaqq = 1.;
#line 1077 "sgesvj.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1079 "sgesvj.f"
					    sva[q] = t * sqrt(aaqq) * work[q];
#line 1080 "sgesvj.f"
					}
#line 1081 "sgesvj.f"
				    }
#line 1082 "sgesvj.f"
				    if (aapp / aapp0 <= rooteps) {
#line 1083 "sgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1085 "sgesvj.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * work[p];
#line 1087 "sgesvj.f"
					} else {
#line 1088 "sgesvj.f"
					    t = 0.;
#line 1089 "sgesvj.f"
					    aapp = 1.;
#line 1090 "sgesvj.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1092 "sgesvj.f"
					    aapp = t * sqrt(aapp) * work[p];
#line 1093 "sgesvj.f"
					}
#line 1094 "sgesvj.f"
					sva[p] = aapp;
#line 1095 "sgesvj.f"
				    }

#line 1097 "sgesvj.f"
				} else {
/*        A(:,p) and A(:,q) already numerically orthogonal */
#line 1099 "sgesvj.f"
				    if (ir1 == 0) {
#line 1099 "sgesvj.f"
					++notrot;
#line 1099 "sgesvj.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1101 "sgesvj.f"
				    ++pskipped;
#line 1102 "sgesvj.f"
				}
#line 1103 "sgesvj.f"
			    } else {
/*        A(:,q) is zero column */
#line 1105 "sgesvj.f"
				if (ir1 == 0) {
#line 1105 "sgesvj.f"
				    ++notrot;
#line 1105 "sgesvj.f"
				}
#line 1106 "sgesvj.f"
				++pskipped;
#line 1107 "sgesvj.f"
			    }

#line 1109 "sgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1111 "sgesvj.f"
				if (ir1 == 0) {
#line 1111 "sgesvj.f"
				    aapp = -aapp;
#line 1111 "sgesvj.f"
				}
#line 1112 "sgesvj.f"
				notrot = 0;
#line 1113 "sgesvj.f"
				goto L2103;
#line 1114 "sgesvj.f"
			    }

#line 1116 "sgesvj.f"
/* L2002: */
#line 1116 "sgesvj.f"
			}
/*     END q-LOOP */

#line 1119 "sgesvj.f"
L2103:
/*     bailed out of q-loop */

#line 1122 "sgesvj.f"
			sva[p] = aapp;

#line 1124 "sgesvj.f"
		    } else {
#line 1125 "sgesvj.f"
			sva[p] = aapp;
#line 1126 "sgesvj.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 1126 "sgesvj.f"
			    i__4 = igl + kbl - 1;
#line 1126 "sgesvj.f"
			    notrot = notrot + min(i__4,*n) - p;
#line 1126 "sgesvj.f"
			}
#line 1128 "sgesvj.f"
		    }

#line 1130 "sgesvj.f"
/* L2001: */
#line 1130 "sgesvj.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 1133 "sgesvj.f"
/* L1002: */
#line 1133 "sgesvj.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 1138 "sgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

#line 1140 "sgesvj.f"
	    i__2 = nbl;
#line 1140 "sgesvj.f"
	    for (jbc = ibr + 1; jbc <= i__2; ++jbc) {

#line 1142 "sgesvj.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 1146 "sgesvj.f"
		ijblsk = 0;
/* Computing MIN */
#line 1147 "sgesvj.f"
		i__4 = igl + kbl - 1;
#line 1147 "sgesvj.f"
		i__3 = min(i__4,*n);
#line 1147 "sgesvj.f"
		for (p = igl; p <= i__3; ++p) {

#line 1149 "sgesvj.f"
		    aapp = sva[p];
#line 1150 "sgesvj.f"
		    if (aapp > 0.) {

#line 1152 "sgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 1154 "sgesvj.f"
			i__5 = jgl + kbl - 1;
#line 1154 "sgesvj.f"
			i__4 = min(i__5,*n);
#line 1154 "sgesvj.f"
			for (q = jgl; q <= i__4; ++q) {

#line 1156 "sgesvj.f"
			    aaqq = sva[q];
#line 1157 "sgesvj.f"
			    if (aaqq > 0.) {
#line 1158 "sgesvj.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 1164 "sgesvj.f"
				if (aaqq >= 1.) {
#line 1165 "sgesvj.f"
				    if (aapp >= aaqq) {
#line 1166 "sgesvj.f"
					rotok = small * aapp <= aaqq;
#line 1167 "sgesvj.f"
				    } else {
#line 1168 "sgesvj.f"
					rotok = small * aaqq <= aapp;
#line 1169 "sgesvj.f"
				    }
#line 1170 "sgesvj.f"
				    if (aapp < big / aaqq) {
#line 1171 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 1174 "sgesvj.f"
				    } else {
#line 1175 "sgesvj.f"
					scopy_(m, &a[p * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1177 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aapp, &
						work[p], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1180 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[q * a_dim1 + 1], &c__1) * 
						work[q] / aaqq;
#line 1182 "sgesvj.f"
				    }
#line 1183 "sgesvj.f"
				} else {
#line 1184 "sgesvj.f"
				    if (aapp >= aaqq) {
#line 1185 "sgesvj.f"
					rotok = aapp <= aaqq / small;
#line 1186 "sgesvj.f"
				    } else {
#line 1187 "sgesvj.f"
					rotok = aaqq <= aapp / small;
#line 1188 "sgesvj.f"
				    }
#line 1189 "sgesvj.f"
				    if (aapp > small / aaqq) {
#line 1190 "sgesvj.f"
					aapq = sdot_(m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1) * work[p] * work[q] / 
						aaqq / aapp;
#line 1193 "sgesvj.f"
				    } else {
#line 1194 "sgesvj.f"
					scopy_(m, &a[q * a_dim1 + 1], &c__1, &
						work[*n + 1], &c__1);
#line 1196 "sgesvj.f"
					slascl_("G", &c__0, &c__0, &aaqq, &
						work[q], m, &c__1, &work[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1199 "sgesvj.f"
					aapq = sdot_(m, &work[*n + 1], &c__1, 
						&a[p * a_dim1 + 1], &c__1) * 
						work[p] / aapp;
#line 1201 "sgesvj.f"
				    }
#line 1202 "sgesvj.f"
				}

/* Computing MAX */
#line 1204 "sgesvj.f"
				d__1 = mxaapq, d__2 = abs(aapq);
#line 1204 "sgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 1208 "sgesvj.f"
				if (abs(aapq) > tol) {
#line 1209 "sgesvj.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 1211 "sgesvj.f"
				    pskipped = 0;
#line 1212 "sgesvj.f"
				    ++iswrot;

#line 1214 "sgesvj.f"
				    if (rotok) {

#line 1216 "sgesvj.f"
					aqoap = aaqq / aapp;
#line 1217 "sgesvj.f"
					apoaq = aapp / aaqq;
#line 1218 "sgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq;
#line 1219 "sgesvj.f"
					if (aaqq > aapp0) {
#line 1219 "sgesvj.f"
					    theta = -theta;
#line 1219 "sgesvj.f"
					}

#line 1221 "sgesvj.f"
					if (abs(theta) > bigtheta) {
#line 1222 "sgesvj.f"
					    t = .5 / theta;
#line 1223 "sgesvj.f"
					    fastr[2] = t * work[p] / work[q];
#line 1224 "sgesvj.f"
					    fastr[3] = -t * work[q] / work[p];
#line 1226 "sgesvj.f"
					    srotm_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, fastr);
#line 1228 "sgesvj.f"
					    if (rsvec) {
#line 1228 "sgesvj.f"
			  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, fastr);
#line 1228 "sgesvj.f"
					    }
/* Computing MAX */
#line 1232 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 1232 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1234 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 1234 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 1236 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 1236 "sgesvj.f"
					    mxsinj = max(d__1,d__2);
#line 1237 "sgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 1241 "sgesvj.f"
					    thsign = -d_sign(&c_b18, &aapq);
#line 1242 "sgesvj.f"
					    if (aaqq > aapp0) {
#line 1242 "sgesvj.f"
			  thsign = -thsign;
#line 1242 "sgesvj.f"
					    }
#line 1243 "sgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 1245 "sgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 1246 "sgesvj.f"
					    sn = t * cs;
/* Computing MAX */
#line 1247 "sgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 1247 "sgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 1248 "sgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq + 1.;
#line 1248 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1250 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq;
#line 1250 "sgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 1253 "sgesvj.f"
					    apoaq = work[p] / work[q];
#line 1254 "sgesvj.f"
					    aqoap = work[q] / work[p];
#line 1255 "sgesvj.f"
					    if (work[p] >= 1.) {

#line 1257 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 1258 "sgesvj.f"
			      fastr[2] = t * apoaq;
#line 1259 "sgesvj.f"
			      fastr[3] = -t * aqoap;
#line 1260 "sgesvj.f"
			      work[p] *= cs;
#line 1261 "sgesvj.f"
			      work[q] *= cs;
#line 1262 "sgesvj.f"
			      srotm_(m, &a[p * a_dim1 + 1], &c__1, &a[q * 
				      a_dim1 + 1], &c__1, fastr);
#line 1265 "sgesvj.f"
			      if (rsvec) {
#line 1265 "sgesvj.f"
				  srotm_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[
					  q * v_dim1 + 1], &c__1, fastr);
#line 1265 "sgesvj.f"
			      }
#line 1268 "sgesvj.f"
			  } else {
#line 1269 "sgesvj.f"
			      d__1 = -t * aqoap;
#line 1269 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1272 "sgesvj.f"
			      d__1 = cs * sn * apoaq;
#line 1272 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 1275 "sgesvj.f"
			      if (rsvec) {
#line 1276 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 1276 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1279 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1279 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1283 "sgesvj.f"
			      }
#line 1284 "sgesvj.f"
			      work[p] *= cs;
#line 1285 "sgesvj.f"
			      work[q] /= cs;
#line 1286 "sgesvj.f"
			  }
#line 1287 "sgesvj.f"
					    } else {
#line 1288 "sgesvj.f"
			  if (work[q] >= 1.) {
#line 1289 "sgesvj.f"
			      d__1 = t * apoaq;
#line 1289 "sgesvj.f"
			      saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, &a[
				      q * a_dim1 + 1], &c__1);
#line 1292 "sgesvj.f"
			      d__1 = -cs * sn * aqoap;
#line 1292 "sgesvj.f"
			      saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, &a[
				      p * a_dim1 + 1], &c__1);
#line 1295 "sgesvj.f"
			      if (rsvec) {
#line 1296 "sgesvj.f"
				  d__1 = t * apoaq;
#line 1296 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], &
					  c__1, &v[q * v_dim1 + 1], &c__1);
#line 1299 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1299 "sgesvj.f"
				  saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], &
					  c__1, &v[p * v_dim1 + 1], &c__1);
#line 1303 "sgesvj.f"
			      }
#line 1304 "sgesvj.f"
			      work[p] /= cs;
#line 1305 "sgesvj.f"
			      work[q] *= cs;
#line 1306 "sgesvj.f"
			  } else {
#line 1307 "sgesvj.f"
			      if (work[p] >= work[q]) {
#line 1309 "sgesvj.f"
				  d__1 = -t * aqoap;
#line 1309 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1312 "sgesvj.f"
				  d__1 = cs * sn * apoaq;
#line 1312 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1315 "sgesvj.f"
				  work[p] *= cs;
#line 1316 "sgesvj.f"
				  work[q] /= cs;
#line 1317 "sgesvj.f"
				  if (rsvec) {
#line 1318 "sgesvj.f"
				      d__1 = -t * aqoap;
#line 1318 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1322 "sgesvj.f"
				      d__1 = cs * sn * apoaq;
#line 1322 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1326 "sgesvj.f"
				  }
#line 1327 "sgesvj.f"
			      } else {
#line 1328 "sgesvj.f"
				  d__1 = t * apoaq;
#line 1328 "sgesvj.f"
				  saxpy_(m, &d__1, &a[p * a_dim1 + 1], &c__1, 
					  &a[q * a_dim1 + 1], &c__1);
#line 1331 "sgesvj.f"
				  d__1 = -cs * sn * aqoap;
#line 1331 "sgesvj.f"
				  saxpy_(m, &d__1, &a[q * a_dim1 + 1], &c__1, 
					  &a[p * a_dim1 + 1], &c__1);
#line 1335 "sgesvj.f"
				  work[p] /= cs;
#line 1336 "sgesvj.f"
				  work[q] *= cs;
#line 1337 "sgesvj.f"
				  if (rsvec) {
#line 1338 "sgesvj.f"
				      d__1 = t * apoaq;
#line 1338 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[p * v_dim1 + 1], 
					      &c__1, &v[q * v_dim1 + 1], &
					      c__1);
#line 1341 "sgesvj.f"
				      d__1 = -cs * sn * aqoap;
#line 1341 "sgesvj.f"
				      saxpy_(&mvl, &d__1, &v[q * v_dim1 + 1], 
					      &c__1, &v[p * v_dim1 + 1], &
					      c__1);
#line 1345 "sgesvj.f"
				  }
#line 1346 "sgesvj.f"
			      }
#line 1347 "sgesvj.f"
			  }
#line 1348 "sgesvj.f"
					    }
#line 1349 "sgesvj.f"
					}

#line 1351 "sgesvj.f"
				    } else {
#line 1352 "sgesvj.f"
					if (aapp > aaqq) {
#line 1353 "sgesvj.f"
					    scopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &work[*n + 1], &
						    c__1);
#line 1355 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &work[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1358 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1361 "sgesvj.f"
					    temp1 = -aapq * work[p] / work[q];
#line 1362 "sgesvj.f"
					    saxpy_(m, &temp1, &work[*n + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1);
#line 1364 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &c_b18,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1367 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 1367 "sgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 1369 "sgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1370 "sgesvj.f"
					} else {
#line 1371 "sgesvj.f"
					    scopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &work[*n + 1], &
						    c__1);
#line 1373 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b18, m, &c__1, &work[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1376 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &aapp, 
						    &c_b18, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1379 "sgesvj.f"
					    temp1 = -aapq * work[q] / work[p];
#line 1380 "sgesvj.f"
					    saxpy_(m, &temp1, &work[*n + 1], &
						    c__1, &a[p * a_dim1 + 1], 
						    &c__1);
#line 1382 "sgesvj.f"
					    slascl_("G", &c__0, &c__0, &c_b18,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1385 "sgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq * 
						    aapq;
#line 1385 "sgesvj.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 1387 "sgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1388 "sgesvj.f"
					}
#line 1389 "sgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q) */
/*           .. recompute SVA(q) */
/* Computing 2nd power */
#line 1394 "sgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1394 "sgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1396 "sgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1398 "sgesvj.f"
					    sva[q] = snrm2_(m, &a[q * a_dim1 
						    + 1], &c__1) * work[q];
#line 1400 "sgesvj.f"
					} else {
#line 1401 "sgesvj.f"
					    t = 0.;
#line 1402 "sgesvj.f"
					    aaqq = 1.;
#line 1403 "sgesvj.f"
					    slassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1405 "sgesvj.f"
					    sva[q] = t * sqrt(aaqq) * work[q];
#line 1406 "sgesvj.f"
					}
#line 1407 "sgesvj.f"
				    }
/* Computing 2nd power */
#line 1408 "sgesvj.f"
				    d__1 = aapp / aapp0;
#line 1408 "sgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1409 "sgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1411 "sgesvj.f"
					    aapp = snrm2_(m, &a[p * a_dim1 + 
						    1], &c__1) * work[p];
#line 1413 "sgesvj.f"
					} else {
#line 1414 "sgesvj.f"
					    t = 0.;
#line 1415 "sgesvj.f"
					    aapp = 1.;
#line 1416 "sgesvj.f"
					    slassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1418 "sgesvj.f"
					    aapp = t * sqrt(aapp) * work[p];
#line 1419 "sgesvj.f"
					}
#line 1420 "sgesvj.f"
					sva[p] = aapp;
#line 1421 "sgesvj.f"
				    }
/*              end of OK rotation */
#line 1423 "sgesvj.f"
				} else {
#line 1424 "sgesvj.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1426 "sgesvj.f"
				    ++pskipped;
#line 1427 "sgesvj.f"
				    ++ijblsk;
#line 1428 "sgesvj.f"
				}
#line 1429 "sgesvj.f"
			    } else {
#line 1430 "sgesvj.f"
				++notrot;
#line 1431 "sgesvj.f"
				++pskipped;
#line 1432 "sgesvj.f"
				++ijblsk;
#line 1433 "sgesvj.f"
			    }

#line 1435 "sgesvj.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 1437 "sgesvj.f"
				sva[p] = aapp;
#line 1438 "sgesvj.f"
				notrot = 0;
#line 1439 "sgesvj.f"
				goto L2011;
#line 1440 "sgesvj.f"
			    }
#line 1441 "sgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1443 "sgesvj.f"
				aapp = -aapp;
#line 1444 "sgesvj.f"
				notrot = 0;
#line 1445 "sgesvj.f"
				goto L2203;
#line 1446 "sgesvj.f"
			    }

#line 1448 "sgesvj.f"
/* L2200: */
#line 1448 "sgesvj.f"
			}
/*        end of the q-loop */
#line 1450 "sgesvj.f"
L2203:

#line 1452 "sgesvj.f"
			sva[p] = aapp;

#line 1454 "sgesvj.f"
		    } else {

#line 1456 "sgesvj.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1456 "sgesvj.f"
			    i__4 = jgl + kbl - 1;
#line 1456 "sgesvj.f"
			    notrot = notrot + min(i__4,*n) - jgl + 1;
#line 1456 "sgesvj.f"
			}
#line 1458 "sgesvj.f"
			if (aapp < 0.) {
#line 1458 "sgesvj.f"
			    notrot = 0;
#line 1458 "sgesvj.f"
			}

#line 1460 "sgesvj.f"
		    }

#line 1462 "sgesvj.f"
/* L2100: */
#line 1462 "sgesvj.f"
		}
/*     end of the p-loop */
#line 1464 "sgesvj.f"
/* L2010: */
#line 1464 "sgesvj.f"
	    }
/*     end of the jbc-loop */
#line 1466 "sgesvj.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1468 "sgesvj.f"
	    i__3 = igl + kbl - 1;
#line 1468 "sgesvj.f"
	    i__2 = min(i__3,*n);
#line 1468 "sgesvj.f"
	    for (p = igl; p <= i__2; ++p) {
#line 1469 "sgesvj.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1470 "sgesvj.f"
/* L2012: */
#line 1470 "sgesvj.f"
	    }
/* ** */
#line 1472 "sgesvj.f"
/* L2000: */
#line 1472 "sgesvj.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1476 "sgesvj.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1478 "sgesvj.f"
	    sva[*n] = snrm2_(m, &a[*n * a_dim1 + 1], &c__1) * work[*n];
#line 1479 "sgesvj.f"
	} else {
#line 1480 "sgesvj.f"
	    t = 0.;
#line 1481 "sgesvj.f"
	    aapp = 1.;
#line 1482 "sgesvj.f"
	    slassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1483 "sgesvj.f"
	    sva[*n] = t * sqrt(aapp) * work[*n];
#line 1484 "sgesvj.f"
	}

/*     Additional steering devices */

#line 1488 "sgesvj.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1488 "sgesvj.f"
	    swband = i__;
#line 1488 "sgesvj.f"
	}

#line 1491 "sgesvj.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * tol && (
		doublereal) (*n) * mxaapq * mxsinj < tol) {
#line 1493 "sgesvj.f"
	    goto L1994;
#line 1494 "sgesvj.f"
	}

#line 1496 "sgesvj.f"
	if (notrot >= emptsw) {
#line 1496 "sgesvj.f"
	    goto L1994;
#line 1496 "sgesvj.f"
	}

#line 1498 "sgesvj.f"
/* L1993: */
#line 1498 "sgesvj.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 1502 "sgesvj.f"
    *info = 29;
#line 1503 "sgesvj.f"
    goto L1995;

#line 1505 "sgesvj.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 1509 "sgesvj.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1511 "sgesvj.f"
L1995:

/*     Sort the singular values and find how many are above */
/*     the underflow threshold. */

#line 1516 "sgesvj.f"
    n2 = 0;
#line 1517 "sgesvj.f"
    n4 = 0;
#line 1518 "sgesvj.f"
    i__1 = *n - 1;
#line 1518 "sgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 1519 "sgesvj.f"
	i__2 = *n - p + 1;
#line 1519 "sgesvj.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1520 "sgesvj.f"
	if (p != q) {
#line 1521 "sgesvj.f"
	    temp1 = sva[p];
#line 1522 "sgesvj.f"
	    sva[p] = sva[q];
#line 1523 "sgesvj.f"
	    sva[q] = temp1;
#line 1524 "sgesvj.f"
	    temp1 = work[p];
#line 1525 "sgesvj.f"
	    work[p] = work[q];
#line 1526 "sgesvj.f"
	    work[q] = temp1;
#line 1527 "sgesvj.f"
	    sswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1528 "sgesvj.f"
	    if (rsvec) {
#line 1528 "sgesvj.f"
		sswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1528 "sgesvj.f"
	    }
#line 1529 "sgesvj.f"
	}
#line 1530 "sgesvj.f"
	if (sva[p] != 0.) {
#line 1531 "sgesvj.f"
	    ++n4;
#line 1532 "sgesvj.f"
	    if (sva[p] * skl > sfmin) {
#line 1532 "sgesvj.f"
		++n2;
#line 1532 "sgesvj.f"
	    }
#line 1533 "sgesvj.f"
	}
#line 1534 "sgesvj.f"
/* L5991: */
#line 1534 "sgesvj.f"
    }
#line 1535 "sgesvj.f"
    if (sva[*n] != 0.) {
#line 1536 "sgesvj.f"
	++n4;
#line 1537 "sgesvj.f"
	if (sva[*n] * skl > sfmin) {
#line 1537 "sgesvj.f"
	    ++n2;
#line 1537 "sgesvj.f"
	}
#line 1538 "sgesvj.f"
    }

/*     Normalize the left singular vectors. */

#line 1542 "sgesvj.f"
    if (lsvec || uctol) {
#line 1543 "sgesvj.f"
	i__1 = n2;
#line 1543 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1544 "sgesvj.f"
	    d__1 = work[p] / sva[p];
#line 1544 "sgesvj.f"
	    sscal_(m, &d__1, &a[p * a_dim1 + 1], &c__1);
#line 1545 "sgesvj.f"
/* L1998: */
#line 1545 "sgesvj.f"
	}
#line 1546 "sgesvj.f"
    }

/*     Scale the product of Jacobi rotations (assemble the fast rotations). */

#line 1550 "sgesvj.f"
    if (rsvec) {
#line 1551 "sgesvj.f"
	if (applv) {
#line 1552 "sgesvj.f"
	    i__1 = *n;
#line 1552 "sgesvj.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1553 "sgesvj.f"
		sscal_(&mvl, &work[p], &v[p * v_dim1 + 1], &c__1);
#line 1554 "sgesvj.f"
/* L2398: */
#line 1554 "sgesvj.f"
	    }
#line 1555 "sgesvj.f"
	} else {
#line 1556 "sgesvj.f"
	    i__1 = *n;
#line 1556 "sgesvj.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1557 "sgesvj.f"
		temp1 = 1. / snrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
#line 1558 "sgesvj.f"
		sscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
#line 1559 "sgesvj.f"
/* L2399: */
#line 1559 "sgesvj.f"
	    }
#line 1560 "sgesvj.f"
	}
#line 1561 "sgesvj.f"
    }

/*     Undo scaling, if necessary (and possible). */
#line 1564 "sgesvj.f"
    if (skl > 1. && sva[1] < big / skl || skl < 1. && sva[max(n2,1)] > sfmin /
	     skl) {
#line 1567 "sgesvj.f"
	i__1 = *n;
#line 1567 "sgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1568 "sgesvj.f"
	    sva[p] = skl * sva[p];
#line 1569 "sgesvj.f"
/* L2400: */
#line 1569 "sgesvj.f"
	}
#line 1570 "sgesvj.f"
	skl = 1.;
#line 1571 "sgesvj.f"
    }

#line 1573 "sgesvj.f"
    work[1] = skl;
/*     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
/*     then some of the singular values may overflow or underflow and */
/*     the spectrum is given in this factored representation. */

#line 1578 "sgesvj.f"
    work[2] = (doublereal) n4;
/*     N4 is the number of computed nonzero singular values of A. */

#line 1581 "sgesvj.f"
    work[3] = (doublereal) n2;
/*     N2 is the number of singular values of A greater than SFMIN. */
/*     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
/*     that may carry some information. */

#line 1586 "sgesvj.f"
    work[4] = (doublereal) i__;
/*     i is the index of the last sweep before declaring convergence. */

#line 1589 "sgesvj.f"
    work[5] = mxaapq;
/*     MXAAPQ is the largest absolute value of scaled pivots in the */
/*     last sweep */

#line 1593 "sgesvj.f"
    work[6] = mxsinj;
/*     MXSINJ is the largest absolute value of the sines of Jacobi angles */
/*     in the last sweep */

#line 1597 "sgesvj.f"
    return 0;
/*     .. */
/*     .. END OF SGESVJ */
/*     .. */
} /* sgesvj_ */

