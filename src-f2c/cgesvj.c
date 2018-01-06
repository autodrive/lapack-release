#line 1 "cgesvj.f"
/* cgesvj.f -- translated by f2c (version 20100827).
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

#line 1 "cgesvj.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b40 = 1.;
static integer c__2 = 2;

/* > \brief \b CGESVJ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESVJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvj.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvj.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvj.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, */
/*                          LDV, CWORK, LWORK, RWORK, LRWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, LRWORK, M, MV, N */
/*       CHARACTER*1        JOBA, JOBU, JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ),  V( LDV, * ), CWORK( LWORK ) */
/*       REAL               RWORK( LRWORK ),  SVA( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* CGESVJ computes the singular value decomposition (SVD) of a complex */
/* M-by-N matrix A, where M >= N. The SVD of A is written as */
/*                                    [++]   [xx]   [x0]   [xx] */
/*              A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx] */
/*                                    [++]   [xx] */
/* where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal */
/* matrix, and V is an N-by-N unitary matrix. The diagonal elements */
/* of SIGMA are the singular values of A. The columns of U and V are the */
/* left and the right singular vectors of A, respectively. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >                 in the array RWORK as RANKA=NINT(RWORK(2)). Also see the */
/* >                 descriptions of SVA and RWORK. The computed columns of U */
/* >                 are mutually numerically orthogonal up to approximately */
/* >                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU.EQ.'C'), */
/* >                 see the description of JOBU. */
/* >                 If INFO .GT. 0, */
/* >                 the procedure CGESVJ did not converge in the given number */
/* >                 of iterations (sweeps). In that case, the computed */
/* >                 columns of U may not be orthogonal up to TOL. The output */
/* >                 U (stored in A), SIGMA (given by the computed singular */
/* >                 values in SVA(1:N)) and V is still a decomposition of the */
/* >                 input matrix A in the sense that the residual */
/* >                 || A - SCALE * U * SIGMA * V^* ||_2 / ||A||_2 is small. */
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
/* >                 the procedure CGESVJ did not converge in the given number */
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
/* >          depending on the value SCALE = RWORK(1), we have: */
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
/* >          the procedure CGESVJ did not converge in the given number of */
/* >          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate. */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* >          MV is INTEGER */
/* >          If JOBV .EQ. 'A', then the product of Jacobi rotations in CGESVJ */
/* >          is applied to the first MV rows of V. See the description of JOBV. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,N) */
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
/* > \param[in,out] CWORK */
/* >          CWORK is COMPLEX array, dimension M+N. */
/* >          Used as work space. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          Length of CWORK, LWORK >= M+N. */
/* > */
/* > \param[in,out] RWORK */
/* >          RWORK is REAL array, dimension max(6,M+N). */
/* >          On entry, */
/* >          If JOBU .EQ. 'C' : */
/* >          RWORK(1) = CTOL, where CTOL defines the threshold for convergence. */
/* >                    The process stops if all columns of A are mutually */
/* >                    orthogonal up to CTOL*EPS, EPS=SLAMCH('E'). */
/* >                    It is required that CTOL >= ONE, i.e. it is not */
/* >                    allowed to force the routine to obtain orthogonality */
/* >                    below EPSILON. */
/* >          On exit, */
/* >          RWORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N) */
/* >                    are the computed singular values of A. */
/* >                    (See description of SVA().) */
/* >          RWORK(2) = NINT(RWORK(2)) is the number of the computed nonzero */
/* >                    singular values. */
/* >          RWORK(3) = NINT(RWORK(3)) is the number of the computed singular */
/* >                    values that are larger than the underflow threshold. */
/* >          RWORK(4) = NINT(RWORK(4)) is the number of sweeps of Jacobi */
/* >                    rotations needed for numerical convergence. */
/* >          RWORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep. */
/* >                    This is useful information in cases when CGESVJ did */
/* >                    not converge, as it can be used to estimate whether */
/* >                    the output is stil useful and for post festum analysis. */
/* >          RWORK(6) = the largest absolute value over all sines of the */
/* >                    Jacobi rotation angles in the last sweep. It can be */
/* >                    useful for a post festum analysis. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* >          LRWORK is INTEGER */
/* >         Length of RWORK, LRWORK >= MAX(6,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0 : successful exit. */
/* >          < 0 : if INFO = -i, then the i-th argument had an illegal value */
/* >          > 0 : CGESVJ did not converge in the maximal allowed number */
/* >                (NSWEEP=30) of sweeps. The output may still be useful. */
/* >                See the description of RWORK. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane */
/* > rotations. In the case of underflow of the tangent of the Jacobi angle, a */
/* > modified Jacobi transformation of Drmac [3] is used. Pivot strategy uses */
/* > column interchanges of de Rijk [1]. The relative accuracy of the computed */
/* > singular values and the accuracy of the computed singular vectors (in */
/* > angle metric) is as guaranteed by the theory of Demmel and Veselic [2]. */
/* > The condition number that determines the accuracy in the full rank case */
/* > is essentially min_{D=diag} kappa(A*D), where kappa(.) is the */
/* > spectral condition number. The best performance of this Jacobi SVD */
/* > procedure is achieved if used in an  accelerated version of Drmac and */
/* > Veselic [4,5], and it is the kernel routine in the SIGMA library [6]. */
/* > Some tunning parameters (marked with [TP]) are available for the */
/* > implementer. */
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
/* > [1] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the */
/* >    singular value decomposition on a vector computer. */
/* >    SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371. */
/* > [2] J. Demmel and K. Veselic: Jacobi method is more accurate than QR. */
/* > [3] Z. Drmac: Implementation of Jacobi rotations for accurate singular */
/* >    value computation in floating point arithmetic. */
/* >    SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222. */
/* > [4] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. */
/* >    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. */
/* >    LAPACK Working note 169. */
/* > [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. */
/* >    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. */
/* >    LAPACK Working note 170. */
/* > [6] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV, */
/* >    QSVD, (H,K)-SVD computations. */
/* >    Department of Mathematics, University of Zagreb, 2008, 2015. */
/* > */
/* > \par Bugs, Examples and Comments: */
/*  ================================= */
/* > */
/* > Please report all bugs and send interesting test examples and comments to */
/* > drmac@math.hr. Thank you. */

/*  ===================================================================== */
/* Subroutine */ int cgesvj_(char *joba, char *jobu, char *jobv, integer *m, 
	integer *n, doublecomplex *a, integer *lda, doublereal *sva, integer *
	mv, doublecomplex *v, integer *ldv, doublecomplex *cwork, integer *
	lwork, doublereal *rwork, integer *lrwork, integer *info, ftnlen 
	joba_len, ftnlen jobu_len, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_sign(doublereal *, doublereal *);

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
    static doublereal aapp;
    static doublecomplex aapq;
    static doublereal aaqq, ctol;
    static integer ierr;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublecomplex ompq;
    static doublereal aapp0, aapq1, temp1, large;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small, sfmin;
    static logical lsvec;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static doublereal epsln;
    static logical applv, rsvec, uctol;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical lower, upper, rotok;
    extern /* Subroutine */ int cgsvj0_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), cgsvj1_(char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), claset_(char *, integer *, integer *,
	     doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer blskip;
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static doublereal mxaapq, thsign, mxsinj;
    static integer emptsw, notrot, iswrot, lkahead;
    static logical goscale, noscale;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
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

#line 397 "cgesvj.f"
    /* Parameter adjustments */
#line 397 "cgesvj.f"
    --sva;
#line 397 "cgesvj.f"
    a_dim1 = *lda;
#line 397 "cgesvj.f"
    a_offset = 1 + a_dim1;
#line 397 "cgesvj.f"
    a -= a_offset;
#line 397 "cgesvj.f"
    v_dim1 = *ldv;
#line 397 "cgesvj.f"
    v_offset = 1 + v_dim1;
#line 397 "cgesvj.f"
    v -= v_offset;
#line 397 "cgesvj.f"
    --cwork;
#line 397 "cgesvj.f"
    --rwork;
#line 397 "cgesvj.f"

#line 397 "cgesvj.f"
    /* Function Body */
#line 397 "cgesvj.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 398 "cgesvj.f"
    uctol = lsame_(jobu, "C", (ftnlen)1, (ftnlen)1);
#line 399 "cgesvj.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 400 "cgesvj.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 401 "cgesvj.f"
    upper = lsame_(joba, "U", (ftnlen)1, (ftnlen)1);
#line 402 "cgesvj.f"
    lower = lsame_(joba, "L", (ftnlen)1, (ftnlen)1);

#line 404 "cgesvj.f"
    if (! (upper || lower || lsame_(joba, "G", (ftnlen)1, (ftnlen)1))) {
#line 405 "cgesvj.f"
	*info = -1;
#line 406 "cgesvj.f"
    } else if (! (lsvec || uctol || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 407 "cgesvj.f"
	*info = -2;
#line 408 "cgesvj.f"
    } else if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 409 "cgesvj.f"
	*info = -3;
#line 410 "cgesvj.f"
    } else if (*m < 0) {
#line 411 "cgesvj.f"
	*info = -4;
#line 412 "cgesvj.f"
    } else if (*n < 0 || *n > *m) {
#line 413 "cgesvj.f"
	*info = -5;
#line 414 "cgesvj.f"
    } else if (*lda < *m) {
#line 415 "cgesvj.f"
	*info = -7;
#line 416 "cgesvj.f"
    } else if (*mv < 0) {
#line 417 "cgesvj.f"
	*info = -9;
#line 418 "cgesvj.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 420 "cgesvj.f"
	*info = -11;
#line 421 "cgesvj.f"
    } else if (uctol && rwork[1] <= 1.) {
#line 422 "cgesvj.f"
	*info = -12;
#line 423 "cgesvj.f"
    } else if (*lwork < *m + *n) {
#line 424 "cgesvj.f"
	*info = -13;
#line 425 "cgesvj.f"
    } else if (*lrwork < max(*n,6)) {
#line 426 "cgesvj.f"
	*info = -15;
#line 427 "cgesvj.f"
    } else {
#line 428 "cgesvj.f"
	*info = 0;
#line 429 "cgesvj.f"
    }

/*     #:( */
#line 432 "cgesvj.f"
    if (*info != 0) {
#line 433 "cgesvj.f"
	i__1 = -(*info);
#line 433 "cgesvj.f"
	xerbla_("CGESVJ", &i__1, (ftnlen)6);
#line 434 "cgesvj.f"
	return 0;
#line 435 "cgesvj.f"
    }

/* #:) Quick return for void matrix */

#line 439 "cgesvj.f"
    if (*m == 0 || *n == 0) {
#line 439 "cgesvj.f"
	return 0;
#line 439 "cgesvj.f"
    }

/*     Set numerical parameters */
/*     The stopping criterion for Jacobi rotations is */

/*     max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS */

/*     where EPS is the round-off and CTOL is defined as follows: */

#line 448 "cgesvj.f"
    if (uctol) {
/*        ... user controlled */
#line 450 "cgesvj.f"
	ctol = rwork[1];
#line 451 "cgesvj.f"
    } else {
/*        ... default */
#line 453 "cgesvj.f"
	if (lsvec || rsvec || applv) {
#line 454 "cgesvj.f"
	    ctol = sqrt((doublereal) (*m));
#line 455 "cgesvj.f"
	} else {
#line 456 "cgesvj.f"
	    ctol = (doublereal) (*m);
#line 457 "cgesvj.f"
	}
#line 458 "cgesvj.f"
    }
/*     ... and the machine dependent parameters are */
/* [!]  (Make sure that SLAMCH() works properly on the target machine.) */

#line 462 "cgesvj.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 463 "cgesvj.f"
    rooteps = sqrt(epsln);
#line 464 "cgesvj.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 465 "cgesvj.f"
    rootsfmin = sqrt(sfmin);
#line 466 "cgesvj.f"
    small = sfmin / epsln;
#line 467 "cgesvj.f"
    big = slamch_("Overflow", (ftnlen)8);
/*     BIG         = ONE    / SFMIN */
#line 469 "cgesvj.f"
    rootbig = 1. / rootsfmin;
#line 470 "cgesvj.f"
    large = big / sqrt((doublereal) (*m * *n));
#line 471 "cgesvj.f"
    bigtheta = 1. / rooteps;

#line 473 "cgesvj.f"
    tol = ctol * epsln;
#line 474 "cgesvj.f"
    roottol = sqrt(tol);

#line 476 "cgesvj.f"
    if ((doublereal) (*m) * epsln >= 1.) {
#line 477 "cgesvj.f"
	*info = -4;
#line 478 "cgesvj.f"
	i__1 = -(*info);
#line 478 "cgesvj.f"
	xerbla_("CGESVJ", &i__1, (ftnlen)6);
#line 479 "cgesvj.f"
	return 0;
#line 480 "cgesvj.f"
    }

/*     Initialize the right singular vector matrix. */

#line 484 "cgesvj.f"
    if (rsvec) {
#line 485 "cgesvj.f"
	mvl = *n;
#line 486 "cgesvj.f"
	claset_("A", &mvl, n, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 487 "cgesvj.f"
    } else if (applv) {
#line 488 "cgesvj.f"
	mvl = *mv;
#line 489 "cgesvj.f"
    }
#line 490 "cgesvj.f"
    rsvec = rsvec || applv;

/*     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
/* (!)  If necessary, scale A to protect the largest singular value */
/*     from overflow. It is possible that saving the largest singular */
/*     value destroys the information about the small ones. */
/*     This initial scaling is almost minimal in the sense that the */
/*     goal is to make sure that no column norm overflows, and that */
/*     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
/*     in A are detected, the procedure returns with INFO=-6. */

#line 501 "cgesvj.f"
    skl = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 502 "cgesvj.f"
    noscale = TRUE_;
#line 503 "cgesvj.f"
    goscale = TRUE_;

#line 505 "cgesvj.f"
    if (lower) {
/*        the input matrix is M-by-N lower triangular (trapezoidal) */
#line 507 "cgesvj.f"
	i__1 = *n;
#line 507 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 508 "cgesvj.f"
	    aapp = 0.;
#line 509 "cgesvj.f"
	    aaqq = 1.;
#line 510 "cgesvj.f"
	    i__2 = *m - p + 1;
#line 510 "cgesvj.f"
	    classq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
#line 511 "cgesvj.f"
	    if (aapp > big) {
#line 512 "cgesvj.f"
		*info = -6;
#line 513 "cgesvj.f"
		i__2 = -(*info);
#line 513 "cgesvj.f"
		xerbla_("CGESVJ", &i__2, (ftnlen)6);
#line 514 "cgesvj.f"
		return 0;
#line 515 "cgesvj.f"
	    }
#line 516 "cgesvj.f"
	    aaqq = sqrt(aaqq);
#line 517 "cgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 518 "cgesvj.f"
		sva[p] = aapp * aaqq;
#line 519 "cgesvj.f"
	    } else {
#line 520 "cgesvj.f"
		noscale = FALSE_;
#line 521 "cgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 522 "cgesvj.f"
		if (goscale) {
#line 523 "cgesvj.f"
		    goscale = FALSE_;
#line 524 "cgesvj.f"
		    i__2 = p - 1;
#line 524 "cgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 525 "cgesvj.f"
			sva[q] *= skl;
#line 526 "cgesvj.f"
/* L1873: */
#line 526 "cgesvj.f"
		    }
#line 527 "cgesvj.f"
		}
#line 528 "cgesvj.f"
	    }
#line 529 "cgesvj.f"
/* L1874: */
#line 529 "cgesvj.f"
	}
#line 530 "cgesvj.f"
    } else if (upper) {
/*        the input matrix is M-by-N upper triangular (trapezoidal) */
#line 532 "cgesvj.f"
	i__1 = *n;
#line 532 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 533 "cgesvj.f"
	    aapp = 0.;
#line 534 "cgesvj.f"
	    aaqq = 1.;
#line 535 "cgesvj.f"
	    classq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 536 "cgesvj.f"
	    if (aapp > big) {
#line 537 "cgesvj.f"
		*info = -6;
#line 538 "cgesvj.f"
		i__2 = -(*info);
#line 538 "cgesvj.f"
		xerbla_("CGESVJ", &i__2, (ftnlen)6);
#line 539 "cgesvj.f"
		return 0;
#line 540 "cgesvj.f"
	    }
#line 541 "cgesvj.f"
	    aaqq = sqrt(aaqq);
#line 542 "cgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 543 "cgesvj.f"
		sva[p] = aapp * aaqq;
#line 544 "cgesvj.f"
	    } else {
#line 545 "cgesvj.f"
		noscale = FALSE_;
#line 546 "cgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 547 "cgesvj.f"
		if (goscale) {
#line 548 "cgesvj.f"
		    goscale = FALSE_;
#line 549 "cgesvj.f"
		    i__2 = p - 1;
#line 549 "cgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 550 "cgesvj.f"
			sva[q] *= skl;
#line 551 "cgesvj.f"
/* L2873: */
#line 551 "cgesvj.f"
		    }
#line 552 "cgesvj.f"
		}
#line 553 "cgesvj.f"
	    }
#line 554 "cgesvj.f"
/* L2874: */
#line 554 "cgesvj.f"
	}
#line 555 "cgesvj.f"
    } else {
/*        the input matrix is M-by-N general dense */
#line 557 "cgesvj.f"
	i__1 = *n;
#line 557 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 558 "cgesvj.f"
	    aapp = 0.;
#line 559 "cgesvj.f"
	    aaqq = 1.;
#line 560 "cgesvj.f"
	    classq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 561 "cgesvj.f"
	    if (aapp > big) {
#line 562 "cgesvj.f"
		*info = -6;
#line 563 "cgesvj.f"
		i__2 = -(*info);
#line 563 "cgesvj.f"
		xerbla_("CGESVJ", &i__2, (ftnlen)6);
#line 564 "cgesvj.f"
		return 0;
#line 565 "cgesvj.f"
	    }
#line 566 "cgesvj.f"
	    aaqq = sqrt(aaqq);
#line 567 "cgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 568 "cgesvj.f"
		sva[p] = aapp * aaqq;
#line 569 "cgesvj.f"
	    } else {
#line 570 "cgesvj.f"
		noscale = FALSE_;
#line 571 "cgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 572 "cgesvj.f"
		if (goscale) {
#line 573 "cgesvj.f"
		    goscale = FALSE_;
#line 574 "cgesvj.f"
		    i__2 = p - 1;
#line 574 "cgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 575 "cgesvj.f"
			sva[q] *= skl;
#line 576 "cgesvj.f"
/* L3873: */
#line 576 "cgesvj.f"
		    }
#line 577 "cgesvj.f"
		}
#line 578 "cgesvj.f"
	    }
#line 579 "cgesvj.f"
/* L3874: */
#line 579 "cgesvj.f"
	}
#line 580 "cgesvj.f"
    }

#line 582 "cgesvj.f"
    if (noscale) {
#line 582 "cgesvj.f"
	skl = 1.;
#line 582 "cgesvj.f"
    }

/*     Move the smaller part of the spectrum from the underflow threshold */
/* (!)  Start by determining the position of the nonzero entries of the */
/*     array SVA() relative to ( SFMIN, BIG ). */

#line 588 "cgesvj.f"
    aapp = 0.;
#line 589 "cgesvj.f"
    aaqq = big;
#line 590 "cgesvj.f"
    i__1 = *n;
#line 590 "cgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 591 "cgesvj.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 591 "cgesvj.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 591 "cgesvj.f"
	    aaqq = min(d__1,d__2);
#line 591 "cgesvj.f"
	}
/* Computing MAX */
#line 592 "cgesvj.f"
	d__1 = aapp, d__2 = sva[p];
#line 592 "cgesvj.f"
	aapp = max(d__1,d__2);
#line 593 "cgesvj.f"
/* L4781: */
#line 593 "cgesvj.f"
    }

/* #:) Quick return for zero matrix */

#line 597 "cgesvj.f"
    if (aapp == 0.) {
#line 598 "cgesvj.f"
	if (lsvec) {
#line 598 "cgesvj.f"
	    claset_("G", m, n, &c_b1, &c_b2, &a[a_offset], lda, (ftnlen)1);
#line 598 "cgesvj.f"
	}
#line 599 "cgesvj.f"
	rwork[1] = 1.;
#line 600 "cgesvj.f"
	rwork[2] = 0.;
#line 601 "cgesvj.f"
	rwork[3] = 0.;
#line 602 "cgesvj.f"
	rwork[4] = 0.;
#line 603 "cgesvj.f"
	rwork[5] = 0.;
#line 604 "cgesvj.f"
	rwork[6] = 0.;
#line 605 "cgesvj.f"
	return 0;
#line 606 "cgesvj.f"
    }

/* #:) Quick return for one-column matrix */

#line 610 "cgesvj.f"
    if (*n == 1) {
#line 611 "cgesvj.f"
	if (lsvec) {
#line 611 "cgesvj.f"
	    clascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1]
		    , lda, &ierr, (ftnlen)1);
#line 611 "cgesvj.f"
	}
#line 613 "cgesvj.f"
	rwork[1] = 1. / skl;
#line 614 "cgesvj.f"
	if (sva[1] >= sfmin) {
#line 615 "cgesvj.f"
	    rwork[2] = 1.;
#line 616 "cgesvj.f"
	} else {
#line 617 "cgesvj.f"
	    rwork[2] = 0.;
#line 618 "cgesvj.f"
	}
#line 619 "cgesvj.f"
	rwork[3] = 0.;
#line 620 "cgesvj.f"
	rwork[4] = 0.;
#line 621 "cgesvj.f"
	rwork[5] = 0.;
#line 622 "cgesvj.f"
	rwork[6] = 0.;
#line 623 "cgesvj.f"
	return 0;
#line 624 "cgesvj.f"
    }

/*     Protect small singular values from underflow, and try to */
/*     avoid underflows/overflows in computing Jacobi rotations. */

#line 629 "cgesvj.f"
    sn = sqrt(sfmin / epsln);
#line 630 "cgesvj.f"
    temp1 = sqrt(big / (doublereal) (*n));
#line 631 "cgesvj.f"
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1) {
/* Computing MIN */
#line 633 "cgesvj.f"
	d__1 = big, d__2 = temp1 / aapp;
#line 633 "cgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 636 "cgesvj.f"
    } else if (aaqq <= sn && aapp <= temp1) {
/* Computing MIN */
#line 637 "cgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (aapp * sqrt((doublereal) (*n)));
#line 637 "cgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 640 "cgesvj.f"
    } else if (aaqq >= sn && aapp >= temp1) {
/* Computing MAX */
#line 641 "cgesvj.f"
	d__1 = sn / aaqq, d__2 = temp1 / aapp;
#line 641 "cgesvj.f"
	temp1 = max(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 644 "cgesvj.f"
    } else if (aaqq <= sn && aapp >= temp1) {
/* Computing MIN */
#line 645 "cgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (sqrt((doublereal) (*n)) * aapp);
#line 645 "cgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 648 "cgesvj.f"
    } else {
#line 649 "cgesvj.f"
	temp1 = 1.;
#line 650 "cgesvj.f"
    }

/*     Scale, if necessary */

#line 654 "cgesvj.f"
    if (temp1 != 1.) {
#line 655 "cgesvj.f"
	slascl_("G", &c__0, &c__0, &c_b40, &temp1, n, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 656 "cgesvj.f"
    }
#line 657 "cgesvj.f"
    skl = temp1 * skl;
#line 658 "cgesvj.f"
    if (skl != 1.) {
#line 659 "cgesvj.f"
	clascl_(joba, &c__0, &c__0, &c_b40, &skl, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 660 "cgesvj.f"
	skl = 1. / skl;
#line 661 "cgesvj.f"
    }

/*     Row-cyclic Jacobi SVD algorithm with column pivoting */

#line 665 "cgesvj.f"
    emptsw = *n * (*n - 1) / 2;
#line 666 "cgesvj.f"
    notrot = 0;
#line 668 "cgesvj.f"
    i__1 = *n;
#line 668 "cgesvj.f"
    for (q = 1; q <= i__1; ++q) {
#line 669 "cgesvj.f"
	i__2 = q;
#line 669 "cgesvj.f"
	cwork[i__2].r = 1., cwork[i__2].i = 0.;
#line 670 "cgesvj.f"
/* L1868: */
#line 670 "cgesvj.f"
    }



#line 674 "cgesvj.f"
    swband = 3;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if CGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm CGEJSV. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 682 "cgesvj.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 688 "cgesvj.f"
    nbl = *n / kbl;
#line 689 "cgesvj.f"
    if (nbl * kbl != *n) {
#line 689 "cgesvj.f"
	++nbl;
#line 689 "cgesvj.f"
    }

/* Computing 2nd power */
#line 691 "cgesvj.f"
    i__1 = kbl;
#line 691 "cgesvj.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 694 "cgesvj.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 697 "cgesvj.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */

/* Computing MAX */
#line 705 "cgesvj.f"
    i__1 = 64, i__2 = kbl << 2;
#line 705 "cgesvj.f"
    if ((lower || upper) && *n > max(i__1,i__2)) {
/* [TP] The number of partition levels and the actual partition are */
/*     tuning parameters. */
#line 708 "cgesvj.f"
	n4 = *n / 4;
#line 709 "cgesvj.f"
	n2 = *n / 2;
#line 710 "cgesvj.f"
	n34 = n4 * 3;
#line 711 "cgesvj.f"
	if (applv) {
#line 712 "cgesvj.f"
	    q = 0;
#line 713 "cgesvj.f"
	} else {
#line 714 "cgesvj.f"
	    q = 1;
#line 715 "cgesvj.f"
	}

#line 717 "cgesvj.f"
	if (lower) {

/*     This works very well on lower triangular matrices, in particular */
/*     in the framework of the preconditioned Jacobi SVD (xGEJSV). */
/*     The idea is simple: */
/*     [+ 0 0 0]   Note that Jacobi transformations of [0 0] */
/*     [+ + 0 0]                                       [0 0] */
/*     [+ + x 0]   actually work on [x 0]              [x 0] */
/*     [+ + x x]                    [x x].             [x x] */

#line 727 "cgesvj.f"
	    i__1 = *m - n34;
#line 727 "cgesvj.f"
	    i__2 = *n - n34;
#line 727 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 727 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda,
		     &cwork[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + (
		    n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &
		    cwork[*n + 1], &i__3, &ierr, (ftnlen)1);
#line 732 "cgesvj.f"
	    i__1 = *m - n2;
#line 732 "cgesvj.f"
	    i__2 = n34 - n2;
#line 732 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 732 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &
		    cwork[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 
		    1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &cwork[*n 
		    + 1], &i__3, &ierr, (ftnlen)1);
#line 737 "cgesvj.f"
	    i__1 = *m - n2;
#line 737 "cgesvj.f"
	    i__2 = *n - n2;
#line 737 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 737 "cgesvj.f"
	    cgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], 
		    lda, &cwork[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (
		    n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 742 "cgesvj.f"
	    i__1 = *m - n4;
#line 742 "cgesvj.f"
	    i__2 = n2 - n4;
#line 742 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 742 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, &
		    cwork[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 
		    1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n 
		    + 1], &i__3, &ierr, (ftnlen)1);

#line 747 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 747 "cgesvj.f"
	    cgsvj0_(jobv, m, &n4, &a[a_offset], lda, &cwork[1], &sva[1], &mvl,
		     &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*
		    n + 1], &i__1, &ierr, (ftnlen)1);

#line 751 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 751 "cgesvj.f"
	    cgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &cwork[1], &sva[1], 
		    &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);


#line 756 "cgesvj.f"
	} else if (upper) {


#line 759 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 759 "cgesvj.f"
	    cgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &cwork[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 763 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 763 "cgesvj.f"
	    cgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &cwork[n4 
		    + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n + 1],
		     &i__1, &ierr, (ftnlen)1);

#line 768 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 768 "cgesvj.f"
	    cgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &cwork[1], &sva[1]
		    , &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 772 "cgesvj.f"
	    i__1 = n2 + n4;
#line 772 "cgesvj.f"
	    i__2 = *lwork - *n;
#line 772 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &cwork[
		    n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n + 1],
		     &i__2, &ierr, (ftnlen)1);
#line 777 "cgesvj.f"
	}

#line 779 "cgesvj.f"
    }

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 783 "cgesvj.f"
    for (i__ = 1; i__ <= 30; ++i__) {

/*     .. go go go ... */

#line 787 "cgesvj.f"
	mxaapq = 0.;
#line 788 "cgesvj.f"
	mxsinj = 0.;
#line 789 "cgesvj.f"
	iswrot = 0;

#line 791 "cgesvj.f"
	notrot = 0;
#line 792 "cgesvj.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 799 "cgesvj.f"
	i__1 = nbl;
#line 799 "cgesvj.f"
	for (ibr = 1; ibr <= i__1; ++ibr) {

#line 801 "cgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 803 "cgesvj.f"
	    i__3 = lkahead, i__4 = nbl - ibr;
#line 803 "cgesvj.f"
	    i__2 = min(i__3,i__4);
#line 803 "cgesvj.f"
	    for (ir1 = 0; ir1 <= i__2; ++ir1) {

#line 805 "cgesvj.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 807 "cgesvj.f"
		i__4 = igl + kbl - 1, i__5 = *n - 1;
#line 807 "cgesvj.f"
		i__3 = min(i__4,i__5);
#line 807 "cgesvj.f"
		for (p = igl; p <= i__3; ++p) {

/*     .. de Rijk's pivoting */

#line 811 "cgesvj.f"
		    i__4 = *n - p + 1;
#line 811 "cgesvj.f"
		    q = isamax_(&i__4, &sva[p], &c__1) + p - 1;
#line 812 "cgesvj.f"
		    if (p != q) {
#line 813 "cgesvj.f"
			cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 814 "cgesvj.f"
			if (rsvec) {
#line 814 "cgesvj.f"
			    cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 814 "cgesvj.f"
			}
#line 816 "cgesvj.f"
			temp1 = sva[p];
#line 817 "cgesvj.f"
			sva[p] = sva[q];
#line 818 "cgesvj.f"
			sva[q] = temp1;
#line 819 "cgesvj.f"
			i__4 = p;
#line 819 "cgesvj.f"
			aapq.r = cwork[i__4].r, aapq.i = cwork[i__4].i;
#line 820 "cgesvj.f"
			i__4 = p;
#line 820 "cgesvj.f"
			i__5 = q;
#line 820 "cgesvj.f"
			cwork[i__4].r = cwork[i__5].r, cwork[i__4].i = cwork[
				i__5].i;
#line 821 "cgesvj.f"
			i__4 = q;
#line 821 "cgesvj.f"
			cwork[i__4].r = aapq.r, cwork[i__4].i = aapq.i;
#line 822 "cgesvj.f"
		    }

#line 824 "cgesvj.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/* [!]     Caveat: */
/*        Unfortunately, some BLAS implementations compute SCNRM2(M,A(1,p),1) */
/*        as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to */
/*        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
/*        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
/*        Hence, SCNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF */
/*        below should be replaced with "AAPP = SCNRM2( M, A(1,p), 1 )". */

#line 838 "cgesvj.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 840 "cgesvj.f"
			    sva[p] = scnrm2_(m, &a[p * a_dim1 + 1], &c__1);
#line 841 "cgesvj.f"
			} else {
#line 842 "cgesvj.f"
			    temp1 = 0.;
#line 843 "cgesvj.f"
			    aapp = 1.;
#line 844 "cgesvj.f"
			    classq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 845 "cgesvj.f"
			    sva[p] = temp1 * sqrt(aapp);
#line 846 "cgesvj.f"
			}
#line 847 "cgesvj.f"
			aapp = sva[p];
#line 848 "cgesvj.f"
		    } else {
#line 849 "cgesvj.f"
			aapp = sva[p];
#line 850 "cgesvj.f"
		    }

#line 852 "cgesvj.f"
		    if (aapp > 0.) {

#line 854 "cgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 856 "cgesvj.f"
			i__5 = igl + kbl - 1;
#line 856 "cgesvj.f"
			i__4 = min(i__5,*n);
#line 856 "cgesvj.f"
			for (q = p + 1; q <= i__4; ++q) {

#line 858 "cgesvj.f"
			    aaqq = sva[q];

#line 860 "cgesvj.f"
			    if (aaqq > 0.) {

#line 862 "cgesvj.f"
				aapp0 = aapp;
#line 863 "cgesvj.f"
				if (aaqq >= 1.) {
#line 864 "cgesvj.f"
				    rotok = small * aapp <= aaqq;
#line 865 "cgesvj.f"
				    if (aapp < big / aaqq) {
#line 866 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 866 "cgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 866 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 866 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 868 "cgesvj.f"
				    } else {
#line 869 "cgesvj.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 871 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b40, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 873 "cgesvj.f"
					cdotc_(&z__2, m, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 873 "cgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 873 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 875 "cgesvj.f"
				    }
#line 876 "cgesvj.f"
				} else {
#line 877 "cgesvj.f"
				    rotok = aapp <= aaqq / small;
#line 878 "cgesvj.f"
				    if (aapp > small / aaqq) {
#line 879 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 879 "cgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 879 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 879 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 881 "cgesvj.f"
				    } else {
#line 882 "cgesvj.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 884 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b40, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 887 "cgesvj.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &cwork[*n + 1], &c__1);
#line 887 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 887 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 889 "cgesvj.f"
				    }
#line 890 "cgesvj.f"
				}

#line 892 "cgesvj.f"
				d__1 = z_abs(&aapq);
#line 892 "cgesvj.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 892 "cgesvj.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q) */
#line 894 "cgesvj.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 895 "cgesvj.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 895 "cgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 899 "cgesvj.f"
				if (abs(aapq1) > tol) {

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 904 "cgesvj.f"
				    if (ir1 == 0) {
#line 905 "cgesvj.f"
					notrot = 0;
#line 906 "cgesvj.f"
					pskipped = 0;
#line 907 "cgesvj.f"
					++iswrot;
#line 908 "cgesvj.f"
				    }

#line 910 "cgesvj.f"
				    if (rotok) {

#line 912 "cgesvj.f"
					aqoap = aaqq / aapp;
#line 913 "cgesvj.f"
					apoaq = aapp / aaqq;
#line 914 "cgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;

#line 916 "cgesvj.f"
					if (abs(theta) > bigtheta) {

#line 918 "cgesvj.f"
					    t = .5 / theta;
#line 919 "cgesvj.f"
					    cs = 1.;
#line 921 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 921 "cgesvj.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 921 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 923 "cgesvj.f"
					    if (rsvec) {
#line 924 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 924 "cgesvj.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 924 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 926 "cgesvj.f"
					    }
/* Computing MAX */
#line 928 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 928 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 930 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 930 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 932 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 932 "cgesvj.f"
					    mxsinj = max(d__1,d__2);

#line 934 "cgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 938 "cgesvj.f"
					    thsign = -d_sign(&c_b40, &aapq1);
#line 939 "cgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 941 "cgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 942 "cgesvj.f"
					    sn = t * cs;

/* Computing MAX */
#line 944 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 944 "cgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 945 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 945 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 947 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 947 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 950 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 950 "cgesvj.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 950 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 952 "cgesvj.f"
					    if (rsvec) {
#line 953 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 953 "cgesvj.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 953 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 955 "cgesvj.f"
					    }
#line 956 "cgesvj.f"
					}
#line 957 "cgesvj.f"
					i__5 = p;
#line 957 "cgesvj.f"
					i__6 = q;
#line 957 "cgesvj.f"
					z__2.r = -cwork[i__6].r, z__2.i = 
						-cwork[i__6].i;
#line 957 "cgesvj.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 957 "cgesvj.f"
					cwork[i__5].r = z__1.r, cwork[i__5].i 
						= z__1.i;

#line 959 "cgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 961 "cgesvj.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 963 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b40, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 966 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b40, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 968 "cgesvj.f"
					z__1.r = -aapq.r, z__1.i = -aapq.i;
#line 968 "cgesvj.f"
					caxpy_(m, &z__1, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 970 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &c_b40, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 972 "cgesvj.f"
					d__1 = 0., d__2 = 1. - aapq1 * aapq1;
#line 972 "cgesvj.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 974 "cgesvj.f"
					mxsinj = max(mxsinj,sfmin);
#line 975 "cgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 981 "cgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 981 "cgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 983 "cgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 985 "cgesvj.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 986 "cgesvj.f"
					} else {
#line 987 "cgesvj.f"
					    t = 0.;
#line 988 "cgesvj.f"
					    aaqq = 1.;
#line 989 "cgesvj.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 991 "cgesvj.f"
					    sva[q] = t * sqrt(aaqq);
#line 992 "cgesvj.f"
					}
#line 993 "cgesvj.f"
				    }
#line 994 "cgesvj.f"
				    if (aapp / aapp0 <= rooteps) {
#line 995 "cgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 997 "cgesvj.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 998 "cgesvj.f"
					} else {
#line 999 "cgesvj.f"
					    t = 0.;
#line 1000 "cgesvj.f"
					    aapp = 1.;
#line 1001 "cgesvj.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1003 "cgesvj.f"
					    aapp = t * sqrt(aapp);
#line 1004 "cgesvj.f"
					}
#line 1005 "cgesvj.f"
					sva[p] = aapp;
#line 1006 "cgesvj.f"
				    }

#line 1008 "cgesvj.f"
				} else {
/*                             A(:,p) and A(:,q) already numerically orthogonal */
#line 1010 "cgesvj.f"
				    if (ir1 == 0) {
#line 1010 "cgesvj.f"
					++notrot;
#line 1010 "cgesvj.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED + 1 */
#line 1012 "cgesvj.f"
				    ++pskipped;
#line 1013 "cgesvj.f"
				}
#line 1014 "cgesvj.f"
			    } else {
/*                          A(:,q) is zero column */
#line 1016 "cgesvj.f"
				if (ir1 == 0) {
#line 1016 "cgesvj.f"
				    ++notrot;
#line 1016 "cgesvj.f"
				}
#line 1017 "cgesvj.f"
				++pskipped;
#line 1018 "cgesvj.f"
			    }

#line 1020 "cgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1022 "cgesvj.f"
				if (ir1 == 0) {
#line 1022 "cgesvj.f"
				    aapp = -aapp;
#line 1022 "cgesvj.f"
				}
#line 1023 "cgesvj.f"
				notrot = 0;
#line 1024 "cgesvj.f"
				goto L2103;
#line 1025 "cgesvj.f"
			    }

#line 1027 "cgesvj.f"
/* L2002: */
#line 1027 "cgesvj.f"
			}
/*     END q-LOOP */

#line 1030 "cgesvj.f"
L2103:
/*     bailed out of q-loop */

#line 1033 "cgesvj.f"
			sva[p] = aapp;

#line 1035 "cgesvj.f"
		    } else {
#line 1036 "cgesvj.f"
			sva[p] = aapp;
#line 1037 "cgesvj.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 1037 "cgesvj.f"
			    i__4 = igl + kbl - 1;
#line 1037 "cgesvj.f"
			    notrot = notrot + min(i__4,*n) - p;
#line 1037 "cgesvj.f"
			}
#line 1039 "cgesvj.f"
		    }

#line 1041 "cgesvj.f"
/* L2001: */
#line 1041 "cgesvj.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 1044 "cgesvj.f"
/* L1002: */
#line 1044 "cgesvj.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 1049 "cgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

#line 1051 "cgesvj.f"
	    i__2 = nbl;
#line 1051 "cgesvj.f"
	    for (jbc = ibr + 1; jbc <= i__2; ++jbc) {

#line 1053 "cgesvj.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 1057 "cgesvj.f"
		ijblsk = 0;
/* Computing MIN */
#line 1058 "cgesvj.f"
		i__4 = igl + kbl - 1;
#line 1058 "cgesvj.f"
		i__3 = min(i__4,*n);
#line 1058 "cgesvj.f"
		for (p = igl; p <= i__3; ++p) {

#line 1060 "cgesvj.f"
		    aapp = sva[p];
#line 1061 "cgesvj.f"
		    if (aapp > 0.) {

#line 1063 "cgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 1065 "cgesvj.f"
			i__5 = jgl + kbl - 1;
#line 1065 "cgesvj.f"
			i__4 = min(i__5,*n);
#line 1065 "cgesvj.f"
			for (q = jgl; q <= i__4; ++q) {

#line 1067 "cgesvj.f"
			    aaqq = sva[q];
#line 1068 "cgesvj.f"
			    if (aaqq > 0.) {
#line 1069 "cgesvj.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 1075 "cgesvj.f"
				if (aaqq >= 1.) {
#line 1076 "cgesvj.f"
				    if (aapp >= aaqq) {
#line 1077 "cgesvj.f"
					rotok = small * aapp <= aaqq;
#line 1078 "cgesvj.f"
				    } else {
#line 1079 "cgesvj.f"
					rotok = small * aaqq <= aapp;
#line 1080 "cgesvj.f"
				    }
#line 1081 "cgesvj.f"
				    if (aapp < big / aaqq) {
#line 1082 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1082 "cgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 1082 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1082 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1084 "cgesvj.f"
				    } else {
#line 1085 "cgesvj.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 1087 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b40, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1090 "cgesvj.f"
					cdotc_(&z__2, m, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1090 "cgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 1090 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1092 "cgesvj.f"
				    }
#line 1093 "cgesvj.f"
				} else {
#line 1094 "cgesvj.f"
				    if (aapp >= aaqq) {
#line 1095 "cgesvj.f"
					rotok = aapp <= aaqq / small;
#line 1096 "cgesvj.f"
				    } else {
#line 1097 "cgesvj.f"
					rotok = aaqq <= aapp / small;
#line 1098 "cgesvj.f"
				    }
#line 1099 "cgesvj.f"
				    if (aapp > small / aaqq) {
#line 1100 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1100 "cgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 1100 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1100 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1102 "cgesvj.f"
				    } else {
#line 1103 "cgesvj.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 1105 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b40, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1108 "cgesvj.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &cwork[*n + 1], &c__1);
#line 1108 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1108 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1110 "cgesvj.f"
				    }
#line 1111 "cgesvj.f"
				}

#line 1113 "cgesvj.f"
				d__1 = z_abs(&aapq);
#line 1113 "cgesvj.f"
				z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					d__1;
#line 1113 "cgesvj.f"
				ompq.r = z__1.r, ompq.i = z__1.i;
/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 1115 "cgesvj.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 1116 "cgesvj.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 1116 "cgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 1120 "cgesvj.f"
				if (abs(aapq1) > tol) {
#line 1121 "cgesvj.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 1123 "cgesvj.f"
				    pskipped = 0;
#line 1124 "cgesvj.f"
				    ++iswrot;

#line 1126 "cgesvj.f"
				    if (rotok) {

#line 1128 "cgesvj.f"
					aqoap = aaqq / aapp;
#line 1129 "cgesvj.f"
					apoaq = aapp / aaqq;
#line 1130 "cgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 1131 "cgesvj.f"
					if (aaqq > aapp0) {
#line 1131 "cgesvj.f"
					    theta = -theta;
#line 1131 "cgesvj.f"
					}

#line 1133 "cgesvj.f"
					if (abs(theta) > bigtheta) {
#line 1134 "cgesvj.f"
					    t = .5 / theta;
#line 1135 "cgesvj.f"
					    cs = 1.;
#line 1136 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 1136 "cgesvj.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 1136 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 1138 "cgesvj.f"
					    if (rsvec) {
#line 1139 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 1139 "cgesvj.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1139 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 1141 "cgesvj.f"
					    }
/* Computing MAX */
#line 1142 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 1142 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1144 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 1144 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 1146 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 1146 "cgesvj.f"
					    mxsinj = max(d__1,d__2);
#line 1147 "cgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 1151 "cgesvj.f"
					    thsign = -d_sign(&c_b40, &aapq1);
#line 1152 "cgesvj.f"
					    if (aaqq > aapp0) {
#line 1152 "cgesvj.f"
			  thsign = -thsign;
#line 1152 "cgesvj.f"
					    }
#line 1153 "cgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 1155 "cgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 1156 "cgesvj.f"
					    sn = t * cs;
/* Computing MAX */
#line 1157 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 1157 "cgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 1158 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 1158 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1160 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 1160 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 1163 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 1163 "cgesvj.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 1163 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 1165 "cgesvj.f"
					    if (rsvec) {
#line 1166 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 1166 "cgesvj.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 1166 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 1168 "cgesvj.f"
					    }
#line 1169 "cgesvj.f"
					}
#line 1170 "cgesvj.f"
					i__5 = p;
#line 1170 "cgesvj.f"
					i__6 = q;
#line 1170 "cgesvj.f"
					z__2.r = -cwork[i__6].r, z__2.i = 
						-cwork[i__6].i;
#line 1170 "cgesvj.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 1170 "cgesvj.f"
					cwork[i__5].r = z__1.r, cwork[i__5].i 
						= z__1.i;

#line 1172 "cgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 1174 "cgesvj.f"
					if (aapp > aaqq) {
#line 1175 "cgesvj.f"
					    ccopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &cwork[*n + 1], &
						    c__1);
#line 1177 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b40, m, &c__1, &cwork[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1180 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b40, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1183 "cgesvj.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 1183 "cgesvj.f"
					    caxpy_(m, &z__1, &cwork[*n + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1);
#line 1185 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &c_b40,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1188 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 1188 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 1190 "cgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1191 "cgesvj.f"
					} else {
#line 1192 "cgesvj.f"
					    ccopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &cwork[*n + 1], &
						    c__1);
#line 1194 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b40, m, &c__1, &cwork[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1197 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b40, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1200 "cgesvj.f"
					    d_cnjg(&z__2, &aapq);
#line 1200 "cgesvj.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 1200 "cgesvj.f"
					    caxpy_(m, &z__1, &cwork[*n + 1], &
						    c__1, &a[p * a_dim1 + 1], 
						    &c__1);
#line 1202 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &c_b40,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1205 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 1205 "cgesvj.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 1207 "cgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1208 "cgesvj.f"
					}
#line 1209 "cgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 1214 "cgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1214 "cgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1216 "cgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1218 "cgesvj.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 1219 "cgesvj.f"
					} else {
#line 1220 "cgesvj.f"
					    t = 0.;
#line 1221 "cgesvj.f"
					    aaqq = 1.;
#line 1222 "cgesvj.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1224 "cgesvj.f"
					    sva[q] = t * sqrt(aaqq);
#line 1225 "cgesvj.f"
					}
#line 1226 "cgesvj.f"
				    }
/* Computing 2nd power */
#line 1227 "cgesvj.f"
				    d__1 = aapp / aapp0;
#line 1227 "cgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1228 "cgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1230 "cgesvj.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 1231 "cgesvj.f"
					} else {
#line 1232 "cgesvj.f"
					    t = 0.;
#line 1233 "cgesvj.f"
					    aapp = 1.;
#line 1234 "cgesvj.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1236 "cgesvj.f"
					    aapp = t * sqrt(aapp);
#line 1237 "cgesvj.f"
					}
#line 1238 "cgesvj.f"
					sva[p] = aapp;
#line 1239 "cgesvj.f"
				    }
/*              end of OK rotation */
#line 1241 "cgesvj.f"
				} else {
#line 1242 "cgesvj.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1244 "cgesvj.f"
				    ++pskipped;
#line 1245 "cgesvj.f"
				    ++ijblsk;
#line 1246 "cgesvj.f"
				}
#line 1247 "cgesvj.f"
			    } else {
#line 1248 "cgesvj.f"
				++notrot;
#line 1249 "cgesvj.f"
				++pskipped;
#line 1250 "cgesvj.f"
				++ijblsk;
#line 1251 "cgesvj.f"
			    }

#line 1253 "cgesvj.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 1255 "cgesvj.f"
				sva[p] = aapp;
#line 1256 "cgesvj.f"
				notrot = 0;
#line 1257 "cgesvj.f"
				goto L2011;
#line 1258 "cgesvj.f"
			    }
#line 1259 "cgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1261 "cgesvj.f"
				aapp = -aapp;
#line 1262 "cgesvj.f"
				notrot = 0;
#line 1263 "cgesvj.f"
				goto L2203;
#line 1264 "cgesvj.f"
			    }

#line 1266 "cgesvj.f"
/* L2200: */
#line 1266 "cgesvj.f"
			}
/*        end of the q-loop */
#line 1268 "cgesvj.f"
L2203:

#line 1270 "cgesvj.f"
			sva[p] = aapp;

#line 1272 "cgesvj.f"
		    } else {

#line 1274 "cgesvj.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1274 "cgesvj.f"
			    i__4 = jgl + kbl - 1;
#line 1274 "cgesvj.f"
			    notrot = notrot + min(i__4,*n) - jgl + 1;
#line 1274 "cgesvj.f"
			}
#line 1276 "cgesvj.f"
			if (aapp < 0.) {
#line 1276 "cgesvj.f"
			    notrot = 0;
#line 1276 "cgesvj.f"
			}

#line 1278 "cgesvj.f"
		    }

#line 1280 "cgesvj.f"
/* L2100: */
#line 1280 "cgesvj.f"
		}
/*     end of the p-loop */
#line 1282 "cgesvj.f"
/* L2010: */
#line 1282 "cgesvj.f"
	    }
/*     end of the jbc-loop */
#line 1284 "cgesvj.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1286 "cgesvj.f"
	    i__3 = igl + kbl - 1;
#line 1286 "cgesvj.f"
	    i__2 = min(i__3,*n);
#line 1286 "cgesvj.f"
	    for (p = igl; p <= i__2; ++p) {
#line 1287 "cgesvj.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1288 "cgesvj.f"
/* L2012: */
#line 1288 "cgesvj.f"
	    }
/* ** */
#line 1290 "cgesvj.f"
/* L2000: */
#line 1290 "cgesvj.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1294 "cgesvj.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1296 "cgesvj.f"
	    sva[*n] = scnrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 1297 "cgesvj.f"
	} else {
#line 1298 "cgesvj.f"
	    t = 0.;
#line 1299 "cgesvj.f"
	    aapp = 1.;
#line 1300 "cgesvj.f"
	    classq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1301 "cgesvj.f"
	    sva[*n] = t * sqrt(aapp);
#line 1302 "cgesvj.f"
	}

/*     Additional steering devices */

#line 1306 "cgesvj.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1306 "cgesvj.f"
	    swband = i__;
#line 1306 "cgesvj.f"
	}

#line 1309 "cgesvj.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * tol && (
		doublereal) (*n) * mxaapq * mxsinj < tol) {
#line 1311 "cgesvj.f"
	    goto L1994;
#line 1312 "cgesvj.f"
	}

#line 1314 "cgesvj.f"
	if (notrot >= emptsw) {
#line 1314 "cgesvj.f"
	    goto L1994;
#line 1314 "cgesvj.f"
	}

#line 1316 "cgesvj.f"
/* L1993: */
#line 1316 "cgesvj.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 1320 "cgesvj.f"
    *info = 29;
#line 1321 "cgesvj.f"
    goto L1995;

#line 1323 "cgesvj.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 1327 "cgesvj.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1329 "cgesvj.f"
L1995:

/*     Sort the singular values and find how many are above */
/*     the underflow threshold. */

#line 1334 "cgesvj.f"
    n2 = 0;
#line 1335 "cgesvj.f"
    n4 = 0;
#line 1336 "cgesvj.f"
    i__1 = *n - 1;
#line 1336 "cgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 1337 "cgesvj.f"
	i__2 = *n - p + 1;
#line 1337 "cgesvj.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1338 "cgesvj.f"
	if (p != q) {
#line 1339 "cgesvj.f"
	    temp1 = sva[p];
#line 1340 "cgesvj.f"
	    sva[p] = sva[q];
#line 1341 "cgesvj.f"
	    sva[q] = temp1;
#line 1342 "cgesvj.f"
	    cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1343 "cgesvj.f"
	    if (rsvec) {
#line 1343 "cgesvj.f"
		cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1343 "cgesvj.f"
	    }
#line 1344 "cgesvj.f"
	}
#line 1345 "cgesvj.f"
	if (sva[p] != 0.) {
#line 1346 "cgesvj.f"
	    ++n4;
#line 1347 "cgesvj.f"
	    if (sva[p] * skl > sfmin) {
#line 1347 "cgesvj.f"
		++n2;
#line 1347 "cgesvj.f"
	    }
#line 1348 "cgesvj.f"
	}
#line 1349 "cgesvj.f"
/* L5991: */
#line 1349 "cgesvj.f"
    }
#line 1350 "cgesvj.f"
    if (sva[*n] != 0.) {
#line 1351 "cgesvj.f"
	++n4;
#line 1352 "cgesvj.f"
	if (sva[*n] * skl > sfmin) {
#line 1352 "cgesvj.f"
	    ++n2;
#line 1352 "cgesvj.f"
	}
#line 1353 "cgesvj.f"
    }

/*     Normalize the left singular vectors. */

#line 1357 "cgesvj.f"
    if (lsvec || uctol) {
#line 1358 "cgesvj.f"
	i__1 = n2;
#line 1358 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1359 "cgesvj.f"
	    d__1 = 1. / sva[p];
#line 1359 "cgesvj.f"
	    csscal_(m, &d__1, &a[p * a_dim1 + 1], &c__1);
#line 1360 "cgesvj.f"
/* L1998: */
#line 1360 "cgesvj.f"
	}
#line 1361 "cgesvj.f"
    }

/*     Scale the product of Jacobi rotations. */

#line 1365 "cgesvj.f"
    if (rsvec) {
#line 1366 "cgesvj.f"
	i__1 = *n;
#line 1366 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1367 "cgesvj.f"
	    temp1 = 1. / scnrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
#line 1368 "cgesvj.f"
	    csscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
#line 1369 "cgesvj.f"
/* L2399: */
#line 1369 "cgesvj.f"
	}
#line 1370 "cgesvj.f"
    }

/*     Undo scaling, if necessary (and possible). */
#line 1373 "cgesvj.f"
    if (skl > 1. && sva[1] < big / skl || skl < 1. && sva[max(n2,1)] > sfmin /
	     skl) {
#line 1376 "cgesvj.f"
	i__1 = *n;
#line 1376 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1377 "cgesvj.f"
	    sva[p] = skl * sva[p];
#line 1378 "cgesvj.f"
/* L2400: */
#line 1378 "cgesvj.f"
	}
#line 1379 "cgesvj.f"
	skl = 1.;
#line 1380 "cgesvj.f"
    }

#line 1382 "cgesvj.f"
    rwork[1] = skl;
/*     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
/*     then some of the singular values may overflow or underflow and */
/*     the spectrum is given in this factored representation. */

#line 1387 "cgesvj.f"
    rwork[2] = (doublereal) n4;
/*     N4 is the number of computed nonzero singular values of A. */

#line 1390 "cgesvj.f"
    rwork[3] = (doublereal) n2;
/*     N2 is the number of singular values of A greater than SFMIN. */
/*     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
/*     that may carry some information. */

#line 1395 "cgesvj.f"
    rwork[4] = (doublereal) i__;
/*     i is the index of the last sweep before declaring convergence. */

#line 1398 "cgesvj.f"
    rwork[5] = mxaapq;
/*     MXAAPQ is the largest absolute value of scaled pivots in the */
/*     last sweep */

#line 1402 "cgesvj.f"
    rwork[6] = mxsinj;
/*     MXSINJ is the largest absolute value of the sines of Jacobi angles */
/*     in the last sweep */

#line 1406 "cgesvj.f"
    return 0;
/*     .. */
/*     .. END OF CGESVJ */
/*     .. */
} /* cgesvj_ */

