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
static doublereal c_b41 = 1.;
static integer c__2 = 2;

/* > \brief <b> CGESVJ </b> */

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
/* > CGESVJ computes the singular value decomposition (SVD) of a complex */
/* > M-by-N matrix A, where M >= N. The SVD of A is written as */
/* >                                    [++]   [xx]   [x0]   [xx] */
/* >              A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx] */
/* >                                    [++]   [xx] */
/* > where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal */
/* > matrix, and V is an N-by-N unitary matrix. The diagonal elements */
/* > of SIGMA are the singular values of A. The columns of U and V are the */
/* > left and the right singular vectors of A, respectively. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBA */
/* > \verbatim */
/* >          JOBA is CHARACTER*1 */
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
/* >          = 'U' or 'F': The left singular vectors corresponding to the nonzero */
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
/* >          = 'V' or 'J': the matrix V is computed and returned in the array V */
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
/* > \verbatim */
/* >          CWORK is COMPLEX array, dimension (max(1,LWORK)) */
/* >          Used as workspace. */
/* >          If on entry LWORK .EQ. -1, then a workspace query is assumed and */
/* >          no computation is done; CWORK(1) is set to the minial (and optimal) */
/* >          length of CWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER. */
/* >          Length of CWORK, LWORK >= M+N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (max(6,LRWORK)) */
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
/* >         If on entry LRWORK .EQ. -1, then a workspace query is assumed and */
/* >         no computation is done; RWORK(1) is set to the minial (and optimal) */
/* >         length of RWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >         LRWORK is INTEGER */
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
/* > */
/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
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
/* > \endverbatim */

/* > \par Contributor: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  ============ */
/* > */
/* >  Zlatko Drmac (Zagreb, Croatia) */
/* > */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
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
/* > \endverbatim */

/* > \par Bugs, examples and comments: */
/*  ================================= */
/* > */
/* > \verbatim */
/* >  =========================== */
/* >  Please report all bugs and send interesting test examples and comments to */
/* >  drmac@math.hr. Thank you. */
/* > \endverbatim */
/* > */
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
    static doublereal aapp0, aapq1, temp1;
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
    static integer emptsw;
    static logical lquery;
    static integer notrot, iswrot, lkahead;
    static logical goscale, noscale;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 421 "cgesvj.f"
    /* Parameter adjustments */
#line 421 "cgesvj.f"
    --sva;
#line 421 "cgesvj.f"
    a_dim1 = *lda;
#line 421 "cgesvj.f"
    a_offset = 1 + a_dim1;
#line 421 "cgesvj.f"
    a -= a_offset;
#line 421 "cgesvj.f"
    v_dim1 = *ldv;
#line 421 "cgesvj.f"
    v_offset = 1 + v_dim1;
#line 421 "cgesvj.f"
    v -= v_offset;
#line 421 "cgesvj.f"
    --cwork;
#line 421 "cgesvj.f"
    --rwork;
#line 421 "cgesvj.f"

#line 421 "cgesvj.f"
    /* Function Body */
#line 421 "cgesvj.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 422 "cgesvj.f"
    uctol = lsame_(jobu, "C", (ftnlen)1, (ftnlen)1);
#line 423 "cgesvj.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || lsame_(jobv, "J", (
	    ftnlen)1, (ftnlen)1);
#line 424 "cgesvj.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 425 "cgesvj.f"
    upper = lsame_(joba, "U", (ftnlen)1, (ftnlen)1);
#line 426 "cgesvj.f"
    lower = lsame_(joba, "L", (ftnlen)1, (ftnlen)1);

#line 428 "cgesvj.f"
    lquery = *lwork == -1 || *lrwork == -1;
#line 429 "cgesvj.f"
    if (! (upper || lower || lsame_(joba, "G", (ftnlen)1, (ftnlen)1))) {
#line 430 "cgesvj.f"
	*info = -1;
#line 431 "cgesvj.f"
    } else if (! (lsvec || uctol || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 432 "cgesvj.f"
	*info = -2;
#line 433 "cgesvj.f"
    } else if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 434 "cgesvj.f"
	*info = -3;
#line 435 "cgesvj.f"
    } else if (*m < 0) {
#line 436 "cgesvj.f"
	*info = -4;
#line 437 "cgesvj.f"
    } else if (*n < 0 || *n > *m) {
#line 438 "cgesvj.f"
	*info = -5;
#line 439 "cgesvj.f"
    } else if (*lda < *m) {
#line 440 "cgesvj.f"
	*info = -7;
#line 441 "cgesvj.f"
    } else if (*mv < 0) {
#line 442 "cgesvj.f"
	*info = -9;
#line 443 "cgesvj.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 445 "cgesvj.f"
	*info = -11;
#line 446 "cgesvj.f"
    } else if (uctol && rwork[1] <= 1.) {
#line 447 "cgesvj.f"
	*info = -12;
#line 448 "cgesvj.f"
    } else if (*lwork < *m + *n && ! lquery) {
#line 449 "cgesvj.f"
	*info = -13;
#line 450 "cgesvj.f"
    } else if (*lrwork < max(*n,6) && ! lquery) {
#line 451 "cgesvj.f"
	*info = -15;
#line 452 "cgesvj.f"
    } else {
#line 453 "cgesvj.f"
	*info = 0;
#line 454 "cgesvj.f"
    }

/*     #:( */
#line 457 "cgesvj.f"
    if (*info != 0) {
#line 458 "cgesvj.f"
	i__1 = -(*info);
#line 458 "cgesvj.f"
	xerbla_("CGESVJ", &i__1, (ftnlen)6);
#line 459 "cgesvj.f"
	return 0;
#line 460 "cgesvj.f"
    } else if (lquery) {
#line 461 "cgesvj.f"
	i__1 = *m + *n;
#line 461 "cgesvj.f"
	cwork[1].r = (doublereal) i__1, cwork[1].i = 0.;
#line 462 "cgesvj.f"
	rwork[1] = (doublereal) max(*n,6);
#line 463 "cgesvj.f"
	return 0;
#line 464 "cgesvj.f"
    }

/* #:) Quick return for void matrix */

#line 468 "cgesvj.f"
    if (*m == 0 || *n == 0) {
#line 468 "cgesvj.f"
	return 0;
#line 468 "cgesvj.f"
    }

/*     Set numerical parameters */
/*     The stopping criterion for Jacobi rotations is */

/*     max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS */

/*     where EPS is the round-off and CTOL is defined as follows: */

#line 477 "cgesvj.f"
    if (uctol) {
/*        ... user controlled */
#line 479 "cgesvj.f"
	ctol = rwork[1];
#line 480 "cgesvj.f"
    } else {
/*        ... default */
#line 482 "cgesvj.f"
	if (lsvec || rsvec || applv) {
#line 483 "cgesvj.f"
	    ctol = sqrt((doublereal) (*m));
#line 484 "cgesvj.f"
	} else {
#line 485 "cgesvj.f"
	    ctol = (doublereal) (*m);
#line 486 "cgesvj.f"
	}
#line 487 "cgesvj.f"
    }
/*     ... and the machine dependent parameters are */
/* [!]  (Make sure that SLAMCH() works properly on the target machine.) */

#line 491 "cgesvj.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 492 "cgesvj.f"
    rooteps = sqrt(epsln);
#line 493 "cgesvj.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 494 "cgesvj.f"
    rootsfmin = sqrt(sfmin);
#line 495 "cgesvj.f"
    small = sfmin / epsln;
/*      BIG = SLAMCH( 'Overflow' ) */
#line 497 "cgesvj.f"
    big = 1. / sfmin;
#line 498 "cgesvj.f"
    rootbig = 1. / rootsfmin;
/*     LARGE = BIG / SQRT( REAL( M*N ) ) */
#line 500 "cgesvj.f"
    bigtheta = 1. / rooteps;

#line 502 "cgesvj.f"
    tol = ctol * epsln;
#line 503 "cgesvj.f"
    roottol = sqrt(tol);

#line 505 "cgesvj.f"
    if ((doublereal) (*m) * epsln >= 1.) {
#line 506 "cgesvj.f"
	*info = -4;
#line 507 "cgesvj.f"
	i__1 = -(*info);
#line 507 "cgesvj.f"
	xerbla_("CGESVJ", &i__1, (ftnlen)6);
#line 508 "cgesvj.f"
	return 0;
#line 509 "cgesvj.f"
    }

/*     Initialize the right singular vector matrix. */

#line 513 "cgesvj.f"
    if (rsvec) {
#line 514 "cgesvj.f"
	mvl = *n;
#line 515 "cgesvj.f"
	claset_("A", &mvl, n, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 516 "cgesvj.f"
    } else if (applv) {
#line 517 "cgesvj.f"
	mvl = *mv;
#line 518 "cgesvj.f"
    }
#line 519 "cgesvj.f"
    rsvec = rsvec || applv;

/*     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
/* (!)  If necessary, scale A to protect the largest singular value */
/*     from overflow. It is possible that saving the largest singular */
/*     value destroys the information about the small ones. */
/*     This initial scaling is almost minimal in the sense that the */
/*     goal is to make sure that no column norm overflows, and that */
/*     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
/*     in A are detected, the procedure returns with INFO=-6. */

#line 530 "cgesvj.f"
    skl = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 531 "cgesvj.f"
    noscale = TRUE_;
#line 532 "cgesvj.f"
    goscale = TRUE_;

#line 534 "cgesvj.f"
    if (lower) {
/*        the input matrix is M-by-N lower triangular (trapezoidal) */
#line 536 "cgesvj.f"
	i__1 = *n;
#line 536 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 537 "cgesvj.f"
	    aapp = 0.;
#line 538 "cgesvj.f"
	    aaqq = 1.;
#line 539 "cgesvj.f"
	    i__2 = *m - p + 1;
#line 539 "cgesvj.f"
	    classq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
#line 540 "cgesvj.f"
	    if (aapp > big) {
#line 541 "cgesvj.f"
		*info = -6;
#line 542 "cgesvj.f"
		i__2 = -(*info);
#line 542 "cgesvj.f"
		xerbla_("CGESVJ", &i__2, (ftnlen)6);
#line 543 "cgesvj.f"
		return 0;
#line 544 "cgesvj.f"
	    }
#line 545 "cgesvj.f"
	    aaqq = sqrt(aaqq);
#line 546 "cgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 547 "cgesvj.f"
		sva[p] = aapp * aaqq;
#line 548 "cgesvj.f"
	    } else {
#line 549 "cgesvj.f"
		noscale = FALSE_;
#line 550 "cgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 551 "cgesvj.f"
		if (goscale) {
#line 552 "cgesvj.f"
		    goscale = FALSE_;
#line 553 "cgesvj.f"
		    i__2 = p - 1;
#line 553 "cgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 554 "cgesvj.f"
			sva[q] *= skl;
#line 555 "cgesvj.f"
/* L1873: */
#line 555 "cgesvj.f"
		    }
#line 556 "cgesvj.f"
		}
#line 557 "cgesvj.f"
	    }
#line 558 "cgesvj.f"
/* L1874: */
#line 558 "cgesvj.f"
	}
#line 559 "cgesvj.f"
    } else if (upper) {
/*        the input matrix is M-by-N upper triangular (trapezoidal) */
#line 561 "cgesvj.f"
	i__1 = *n;
#line 561 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 562 "cgesvj.f"
	    aapp = 0.;
#line 563 "cgesvj.f"
	    aaqq = 1.;
#line 564 "cgesvj.f"
	    classq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 565 "cgesvj.f"
	    if (aapp > big) {
#line 566 "cgesvj.f"
		*info = -6;
#line 567 "cgesvj.f"
		i__2 = -(*info);
#line 567 "cgesvj.f"
		xerbla_("CGESVJ", &i__2, (ftnlen)6);
#line 568 "cgesvj.f"
		return 0;
#line 569 "cgesvj.f"
	    }
#line 570 "cgesvj.f"
	    aaqq = sqrt(aaqq);
#line 571 "cgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 572 "cgesvj.f"
		sva[p] = aapp * aaqq;
#line 573 "cgesvj.f"
	    } else {
#line 574 "cgesvj.f"
		noscale = FALSE_;
#line 575 "cgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 576 "cgesvj.f"
		if (goscale) {
#line 577 "cgesvj.f"
		    goscale = FALSE_;
#line 578 "cgesvj.f"
		    i__2 = p - 1;
#line 578 "cgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 579 "cgesvj.f"
			sva[q] *= skl;
#line 580 "cgesvj.f"
/* L2873: */
#line 580 "cgesvj.f"
		    }
#line 581 "cgesvj.f"
		}
#line 582 "cgesvj.f"
	    }
#line 583 "cgesvj.f"
/* L2874: */
#line 583 "cgesvj.f"
	}
#line 584 "cgesvj.f"
    } else {
/*        the input matrix is M-by-N general dense */
#line 586 "cgesvj.f"
	i__1 = *n;
#line 586 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 587 "cgesvj.f"
	    aapp = 0.;
#line 588 "cgesvj.f"
	    aaqq = 1.;
#line 589 "cgesvj.f"
	    classq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 590 "cgesvj.f"
	    if (aapp > big) {
#line 591 "cgesvj.f"
		*info = -6;
#line 592 "cgesvj.f"
		i__2 = -(*info);
#line 592 "cgesvj.f"
		xerbla_("CGESVJ", &i__2, (ftnlen)6);
#line 593 "cgesvj.f"
		return 0;
#line 594 "cgesvj.f"
	    }
#line 595 "cgesvj.f"
	    aaqq = sqrt(aaqq);
#line 596 "cgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 597 "cgesvj.f"
		sva[p] = aapp * aaqq;
#line 598 "cgesvj.f"
	    } else {
#line 599 "cgesvj.f"
		noscale = FALSE_;
#line 600 "cgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 601 "cgesvj.f"
		if (goscale) {
#line 602 "cgesvj.f"
		    goscale = FALSE_;
#line 603 "cgesvj.f"
		    i__2 = p - 1;
#line 603 "cgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 604 "cgesvj.f"
			sva[q] *= skl;
#line 605 "cgesvj.f"
/* L3873: */
#line 605 "cgesvj.f"
		    }
#line 606 "cgesvj.f"
		}
#line 607 "cgesvj.f"
	    }
#line 608 "cgesvj.f"
/* L3874: */
#line 608 "cgesvj.f"
	}
#line 609 "cgesvj.f"
    }

#line 611 "cgesvj.f"
    if (noscale) {
#line 611 "cgesvj.f"
	skl = 1.;
#line 611 "cgesvj.f"
    }

/*     Move the smaller part of the spectrum from the underflow threshold */
/* (!)  Start by determining the position of the nonzero entries of the */
/*     array SVA() relative to ( SFMIN, BIG ). */

#line 617 "cgesvj.f"
    aapp = 0.;
#line 618 "cgesvj.f"
    aaqq = big;
#line 619 "cgesvj.f"
    i__1 = *n;
#line 619 "cgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 620 "cgesvj.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 620 "cgesvj.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 620 "cgesvj.f"
	    aaqq = min(d__1,d__2);
#line 620 "cgesvj.f"
	}
/* Computing MAX */
#line 621 "cgesvj.f"
	d__1 = aapp, d__2 = sva[p];
#line 621 "cgesvj.f"
	aapp = max(d__1,d__2);
#line 622 "cgesvj.f"
/* L4781: */
#line 622 "cgesvj.f"
    }

/* #:) Quick return for zero matrix */

#line 626 "cgesvj.f"
    if (aapp == 0.) {
#line 627 "cgesvj.f"
	if (lsvec) {
#line 627 "cgesvj.f"
	    claset_("G", m, n, &c_b1, &c_b2, &a[a_offset], lda, (ftnlen)1);
#line 627 "cgesvj.f"
	}
#line 628 "cgesvj.f"
	rwork[1] = 1.;
#line 629 "cgesvj.f"
	rwork[2] = 0.;
#line 630 "cgesvj.f"
	rwork[3] = 0.;
#line 631 "cgesvj.f"
	rwork[4] = 0.;
#line 632 "cgesvj.f"
	rwork[5] = 0.;
#line 633 "cgesvj.f"
	rwork[6] = 0.;
#line 634 "cgesvj.f"
	return 0;
#line 635 "cgesvj.f"
    }

/* #:) Quick return for one-column matrix */

#line 639 "cgesvj.f"
    if (*n == 1) {
#line 640 "cgesvj.f"
	if (lsvec) {
#line 640 "cgesvj.f"
	    clascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1]
		    , lda, &ierr, (ftnlen)1);
#line 640 "cgesvj.f"
	}
#line 642 "cgesvj.f"
	rwork[1] = 1. / skl;
#line 643 "cgesvj.f"
	if (sva[1] >= sfmin) {
#line 644 "cgesvj.f"
	    rwork[2] = 1.;
#line 645 "cgesvj.f"
	} else {
#line 646 "cgesvj.f"
	    rwork[2] = 0.;
#line 647 "cgesvj.f"
	}
#line 648 "cgesvj.f"
	rwork[3] = 0.;
#line 649 "cgesvj.f"
	rwork[4] = 0.;
#line 650 "cgesvj.f"
	rwork[5] = 0.;
#line 651 "cgesvj.f"
	rwork[6] = 0.;
#line 652 "cgesvj.f"
	return 0;
#line 653 "cgesvj.f"
    }

/*     Protect small singular values from underflow, and try to */
/*     avoid underflows/overflows in computing Jacobi rotations. */

#line 658 "cgesvj.f"
    sn = sqrt(sfmin / epsln);
#line 659 "cgesvj.f"
    temp1 = sqrt(big / (doublereal) (*n));
#line 660 "cgesvj.f"
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1) {
/* Computing MIN */
#line 662 "cgesvj.f"
	d__1 = big, d__2 = temp1 / aapp;
#line 662 "cgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 665 "cgesvj.f"
    } else if (aaqq <= sn && aapp <= temp1) {
/* Computing MIN */
#line 666 "cgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (aapp * sqrt((doublereal) (*n)));
#line 666 "cgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 669 "cgesvj.f"
    } else if (aaqq >= sn && aapp >= temp1) {
/* Computing MAX */
#line 670 "cgesvj.f"
	d__1 = sn / aaqq, d__2 = temp1 / aapp;
#line 670 "cgesvj.f"
	temp1 = max(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 673 "cgesvj.f"
    } else if (aaqq <= sn && aapp >= temp1) {
/* Computing MIN */
#line 674 "cgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (sqrt((doublereal) (*n)) * aapp);
#line 674 "cgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 677 "cgesvj.f"
    } else {
#line 678 "cgesvj.f"
	temp1 = 1.;
#line 679 "cgesvj.f"
    }

/*     Scale, if necessary */

#line 683 "cgesvj.f"
    if (temp1 != 1.) {
#line 684 "cgesvj.f"
	slascl_("G", &c__0, &c__0, &c_b41, &temp1, n, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 685 "cgesvj.f"
    }
#line 686 "cgesvj.f"
    skl = temp1 * skl;
#line 687 "cgesvj.f"
    if (skl != 1.) {
#line 688 "cgesvj.f"
	clascl_(joba, &c__0, &c__0, &c_b41, &skl, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 689 "cgesvj.f"
	skl = 1. / skl;
#line 690 "cgesvj.f"
    }

/*     Row-cyclic Jacobi SVD algorithm with column pivoting */

#line 694 "cgesvj.f"
    emptsw = *n * (*n - 1) / 2;
#line 695 "cgesvj.f"
    notrot = 0;
#line 697 "cgesvj.f"
    i__1 = *n;
#line 697 "cgesvj.f"
    for (q = 1; q <= i__1; ++q) {
#line 698 "cgesvj.f"
	i__2 = q;
#line 698 "cgesvj.f"
	cwork[i__2].r = 1., cwork[i__2].i = 0.;
#line 699 "cgesvj.f"
/* L1868: */
#line 699 "cgesvj.f"
    }



#line 703 "cgesvj.f"
    swband = 3;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if CGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm CGEJSV. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 711 "cgesvj.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 717 "cgesvj.f"
    nbl = *n / kbl;
#line 718 "cgesvj.f"
    if (nbl * kbl != *n) {
#line 718 "cgesvj.f"
	++nbl;
#line 718 "cgesvj.f"
    }

/* Computing 2nd power */
#line 720 "cgesvj.f"
    i__1 = kbl;
#line 720 "cgesvj.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 723 "cgesvj.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 726 "cgesvj.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */

/* Computing MAX */
#line 734 "cgesvj.f"
    i__1 = 64, i__2 = kbl << 2;
#line 734 "cgesvj.f"
    if ((lower || upper) && *n > max(i__1,i__2)) {
/* [TP] The number of partition levels and the actual partition are */
/*     tuning parameters. */
#line 737 "cgesvj.f"
	n4 = *n / 4;
#line 738 "cgesvj.f"
	n2 = *n / 2;
#line 739 "cgesvj.f"
	n34 = n4 * 3;
#line 740 "cgesvj.f"
	if (applv) {
#line 741 "cgesvj.f"
	    q = 0;
#line 742 "cgesvj.f"
	} else {
#line 743 "cgesvj.f"
	    q = 1;
#line 744 "cgesvj.f"
	}

#line 746 "cgesvj.f"
	if (lower) {

/*     This works very well on lower triangular matrices, in particular */
/*     in the framework of the preconditioned Jacobi SVD (xGEJSV). */
/*     The idea is simple: */
/*     [+ 0 0 0]   Note that Jacobi transformations of [0 0] */
/*     [+ + 0 0]                                       [0 0] */
/*     [+ + x 0]   actually work on [x 0]              [x 0] */
/*     [+ + x x]                    [x x].             [x x] */

#line 756 "cgesvj.f"
	    i__1 = *m - n34;
#line 756 "cgesvj.f"
	    i__2 = *n - n34;
#line 756 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 756 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda,
		     &cwork[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + (
		    n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &
		    cwork[*n + 1], &i__3, &ierr, (ftnlen)1);
#line 761 "cgesvj.f"
	    i__1 = *m - n2;
#line 761 "cgesvj.f"
	    i__2 = n34 - n2;
#line 761 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 761 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &
		    cwork[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 
		    1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &cwork[*n 
		    + 1], &i__3, &ierr, (ftnlen)1);
#line 766 "cgesvj.f"
	    i__1 = *m - n2;
#line 766 "cgesvj.f"
	    i__2 = *n - n2;
#line 766 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 766 "cgesvj.f"
	    cgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], 
		    lda, &cwork[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (
		    n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__3, &ierr, (ftnlen)1);

#line 771 "cgesvj.f"
	    i__1 = *m - n4;
#line 771 "cgesvj.f"
	    i__2 = n2 - n4;
#line 771 "cgesvj.f"
	    i__3 = *lwork - *n;
#line 771 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, &
		    cwork[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 
		    1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n 
		    + 1], &i__3, &ierr, (ftnlen)1);

#line 776 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 776 "cgesvj.f"
	    cgsvj0_(jobv, m, &n4, &a[a_offset], lda, &cwork[1], &sva[1], &mvl,
		     &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*
		    n + 1], &i__1, &ierr, (ftnlen)1);

#line 780 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 780 "cgesvj.f"
	    cgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &cwork[1], &sva[1], 
		    &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);


#line 785 "cgesvj.f"
	} else if (upper) {


#line 788 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 788 "cgesvj.f"
	    cgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &cwork[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 792 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 792 "cgesvj.f"
	    cgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &cwork[n4 
		    + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n + 1],
		     &i__1, &ierr, (ftnlen)1);

#line 797 "cgesvj.f"
	    i__1 = *lwork - *n;
#line 797 "cgesvj.f"
	    cgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &cwork[1], &sva[1]
		    , &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 801 "cgesvj.f"
	    i__1 = n2 + n4;
#line 801 "cgesvj.f"
	    i__2 = *lwork - *n;
#line 801 "cgesvj.f"
	    cgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &cwork[
		    n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n + 1],
		     &i__2, &ierr, (ftnlen)1);
#line 806 "cgesvj.f"
	}

#line 808 "cgesvj.f"
    }

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 812 "cgesvj.f"
    for (i__ = 1; i__ <= 30; ++i__) {

/*     .. go go go ... */

#line 816 "cgesvj.f"
	mxaapq = 0.;
#line 817 "cgesvj.f"
	mxsinj = 0.;
#line 818 "cgesvj.f"
	iswrot = 0;

#line 820 "cgesvj.f"
	notrot = 0;
#line 821 "cgesvj.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 828 "cgesvj.f"
	i__1 = nbl;
#line 828 "cgesvj.f"
	for (ibr = 1; ibr <= i__1; ++ibr) {

#line 830 "cgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 832 "cgesvj.f"
	    i__3 = lkahead, i__4 = nbl - ibr;
#line 832 "cgesvj.f"
	    i__2 = min(i__3,i__4);
#line 832 "cgesvj.f"
	    for (ir1 = 0; ir1 <= i__2; ++ir1) {

#line 834 "cgesvj.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 836 "cgesvj.f"
		i__4 = igl + kbl - 1, i__5 = *n - 1;
#line 836 "cgesvj.f"
		i__3 = min(i__4,i__5);
#line 836 "cgesvj.f"
		for (p = igl; p <= i__3; ++p) {

/*     .. de Rijk's pivoting */

#line 840 "cgesvj.f"
		    i__4 = *n - p + 1;
#line 840 "cgesvj.f"
		    q = isamax_(&i__4, &sva[p], &c__1) + p - 1;
#line 841 "cgesvj.f"
		    if (p != q) {
#line 842 "cgesvj.f"
			cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 843 "cgesvj.f"
			if (rsvec) {
#line 843 "cgesvj.f"
			    cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 843 "cgesvj.f"
			}
#line 845 "cgesvj.f"
			temp1 = sva[p];
#line 846 "cgesvj.f"
			sva[p] = sva[q];
#line 847 "cgesvj.f"
			sva[q] = temp1;
#line 848 "cgesvj.f"
			i__4 = p;
#line 848 "cgesvj.f"
			aapq.r = cwork[i__4].r, aapq.i = cwork[i__4].i;
#line 849 "cgesvj.f"
			i__4 = p;
#line 849 "cgesvj.f"
			i__5 = q;
#line 849 "cgesvj.f"
			cwork[i__4].r = cwork[i__5].r, cwork[i__4].i = cwork[
				i__5].i;
#line 850 "cgesvj.f"
			i__4 = q;
#line 850 "cgesvj.f"
			cwork[i__4].r = aapq.r, cwork[i__4].i = aapq.i;
#line 851 "cgesvj.f"
		    }

#line 853 "cgesvj.f"
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

#line 867 "cgesvj.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 869 "cgesvj.f"
			    sva[p] = scnrm2_(m, &a[p * a_dim1 + 1], &c__1);
#line 870 "cgesvj.f"
			} else {
#line 871 "cgesvj.f"
			    temp1 = 0.;
#line 872 "cgesvj.f"
			    aapp = 1.;
#line 873 "cgesvj.f"
			    classq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 874 "cgesvj.f"
			    sva[p] = temp1 * sqrt(aapp);
#line 875 "cgesvj.f"
			}
#line 876 "cgesvj.f"
			aapp = sva[p];
#line 877 "cgesvj.f"
		    } else {
#line 878 "cgesvj.f"
			aapp = sva[p];
#line 879 "cgesvj.f"
		    }

#line 881 "cgesvj.f"
		    if (aapp > 0.) {

#line 883 "cgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 885 "cgesvj.f"
			i__5 = igl + kbl - 1;
#line 885 "cgesvj.f"
			i__4 = min(i__5,*n);
#line 885 "cgesvj.f"
			for (q = p + 1; q <= i__4; ++q) {

#line 887 "cgesvj.f"
			    aaqq = sva[q];

#line 889 "cgesvj.f"
			    if (aaqq > 0.) {

#line 891 "cgesvj.f"
				aapp0 = aapp;
#line 892 "cgesvj.f"
				if (aaqq >= 1.) {
#line 893 "cgesvj.f"
				    rotok = small * aapp <= aaqq;
#line 894 "cgesvj.f"
				    if (aapp < big / aaqq) {
#line 895 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 895 "cgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 895 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 895 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 897 "cgesvj.f"
				    } else {
#line 898 "cgesvj.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 900 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b41, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 902 "cgesvj.f"
					cdotc_(&z__2, m, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 902 "cgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 902 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 904 "cgesvj.f"
				    }
#line 905 "cgesvj.f"
				} else {
#line 906 "cgesvj.f"
				    rotok = aapp <= aaqq / small;
#line 907 "cgesvj.f"
				    if (aapp > small / aaqq) {
#line 908 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 908 "cgesvj.f"
					z__2.r = z__3.r / aapp, z__2.i = 
						z__3.i / aapp;
#line 908 "cgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 908 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 910 "cgesvj.f"
				    } else {
#line 911 "cgesvj.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 913 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b41, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 916 "cgesvj.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &cwork[*n + 1], &c__1);
#line 916 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 916 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 918 "cgesvj.f"
				    }
#line 919 "cgesvj.f"
				}

/*                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q) */
#line 922 "cgesvj.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 923 "cgesvj.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 923 "cgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 927 "cgesvj.f"
				if (abs(aapq1) > tol) {
#line 928 "cgesvj.f"
				    d__1 = z_abs(&aapq);
#line 928 "cgesvj.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 928 "cgesvj.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 933 "cgesvj.f"
				    if (ir1 == 0) {
#line 934 "cgesvj.f"
					notrot = 0;
#line 935 "cgesvj.f"
					pskipped = 0;
#line 936 "cgesvj.f"
					++iswrot;
#line 937 "cgesvj.f"
				    }

#line 939 "cgesvj.f"
				    if (rotok) {

#line 941 "cgesvj.f"
					aqoap = aaqq / aapp;
#line 942 "cgesvj.f"
					apoaq = aapp / aaqq;
#line 943 "cgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;

#line 945 "cgesvj.f"
					if (abs(theta) > bigtheta) {

#line 947 "cgesvj.f"
					    t = .5 / theta;
#line 948 "cgesvj.f"
					    cs = 1.;
#line 950 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 950 "cgesvj.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 950 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 952 "cgesvj.f"
					    if (rsvec) {
#line 953 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 953 "cgesvj.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 953 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 955 "cgesvj.f"
					    }
/* Computing MAX */
#line 957 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 957 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 959 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 959 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 961 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 961 "cgesvj.f"
					    mxsinj = max(d__1,d__2);

#line 963 "cgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 967 "cgesvj.f"
					    thsign = -d_sign(&c_b41, &aapq1);
#line 968 "cgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 970 "cgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 971 "cgesvj.f"
					    sn = t * cs;

/* Computing MAX */
#line 973 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 973 "cgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 974 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 974 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 976 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 976 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 979 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 979 "cgesvj.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 979 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 981 "cgesvj.f"
					    if (rsvec) {
#line 982 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 982 "cgesvj.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 982 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 984 "cgesvj.f"
					    }
#line 985 "cgesvj.f"
					}
#line 986 "cgesvj.f"
					i__5 = p;
#line 986 "cgesvj.f"
					i__6 = q;
#line 986 "cgesvj.f"
					z__2.r = -cwork[i__6].r, z__2.i = 
						-cwork[i__6].i;
#line 986 "cgesvj.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 986 "cgesvj.f"
					cwork[i__5].r = z__1.r, cwork[i__5].i 
						= z__1.i;

#line 988 "cgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 990 "cgesvj.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 992 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b41, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 995 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b41, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 997 "cgesvj.f"
					z__1.r = -aapq.r, z__1.i = -aapq.i;
#line 997 "cgesvj.f"
					caxpy_(m, &z__1, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 999 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &c_b41, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 1001 "cgesvj.f"
					d__1 = 0., d__2 = 1. - aapq1 * aapq1;
#line 1001 "cgesvj.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 1003 "cgesvj.f"
					mxsinj = max(mxsinj,sfmin);
#line 1004 "cgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 1010 "cgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1010 "cgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1012 "cgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1014 "cgesvj.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 1015 "cgesvj.f"
					} else {
#line 1016 "cgesvj.f"
					    t = 0.;
#line 1017 "cgesvj.f"
					    aaqq = 1.;
#line 1018 "cgesvj.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1020 "cgesvj.f"
					    sva[q] = t * sqrt(aaqq);
#line 1021 "cgesvj.f"
					}
#line 1022 "cgesvj.f"
				    }
#line 1023 "cgesvj.f"
				    if (aapp / aapp0 <= rooteps) {
#line 1024 "cgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1026 "cgesvj.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 1027 "cgesvj.f"
					} else {
#line 1028 "cgesvj.f"
					    t = 0.;
#line 1029 "cgesvj.f"
					    aapp = 1.;
#line 1030 "cgesvj.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1032 "cgesvj.f"
					    aapp = t * sqrt(aapp);
#line 1033 "cgesvj.f"
					}
#line 1034 "cgesvj.f"
					sva[p] = aapp;
#line 1035 "cgesvj.f"
				    }

#line 1037 "cgesvj.f"
				} else {
/*                             A(:,p) and A(:,q) already numerically orthogonal */
#line 1039 "cgesvj.f"
				    if (ir1 == 0) {
#line 1039 "cgesvj.f"
					++notrot;
#line 1039 "cgesvj.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED + 1 */
#line 1041 "cgesvj.f"
				    ++pskipped;
#line 1042 "cgesvj.f"
				}
#line 1043 "cgesvj.f"
			    } else {
/*                          A(:,q) is zero column */
#line 1045 "cgesvj.f"
				if (ir1 == 0) {
#line 1045 "cgesvj.f"
				    ++notrot;
#line 1045 "cgesvj.f"
				}
#line 1046 "cgesvj.f"
				++pskipped;
#line 1047 "cgesvj.f"
			    }

#line 1049 "cgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1051 "cgesvj.f"
				if (ir1 == 0) {
#line 1051 "cgesvj.f"
				    aapp = -aapp;
#line 1051 "cgesvj.f"
				}
#line 1052 "cgesvj.f"
				notrot = 0;
#line 1053 "cgesvj.f"
				goto L2103;
#line 1054 "cgesvj.f"
			    }

#line 1056 "cgesvj.f"
/* L2002: */
#line 1056 "cgesvj.f"
			}
/*     END q-LOOP */

#line 1059 "cgesvj.f"
L2103:
/*     bailed out of q-loop */

#line 1062 "cgesvj.f"
			sva[p] = aapp;

#line 1064 "cgesvj.f"
		    } else {
#line 1065 "cgesvj.f"
			sva[p] = aapp;
#line 1066 "cgesvj.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 1066 "cgesvj.f"
			    i__4 = igl + kbl - 1;
#line 1066 "cgesvj.f"
			    notrot = notrot + min(i__4,*n) - p;
#line 1066 "cgesvj.f"
			}
#line 1068 "cgesvj.f"
		    }

#line 1070 "cgesvj.f"
/* L2001: */
#line 1070 "cgesvj.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 1073 "cgesvj.f"
/* L1002: */
#line 1073 "cgesvj.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 1078 "cgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

#line 1080 "cgesvj.f"
	    i__2 = nbl;
#line 1080 "cgesvj.f"
	    for (jbc = ibr + 1; jbc <= i__2; ++jbc) {

#line 1082 "cgesvj.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 1086 "cgesvj.f"
		ijblsk = 0;
/* Computing MIN */
#line 1087 "cgesvj.f"
		i__4 = igl + kbl - 1;
#line 1087 "cgesvj.f"
		i__3 = min(i__4,*n);
#line 1087 "cgesvj.f"
		for (p = igl; p <= i__3; ++p) {

#line 1089 "cgesvj.f"
		    aapp = sva[p];
#line 1090 "cgesvj.f"
		    if (aapp > 0.) {

#line 1092 "cgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 1094 "cgesvj.f"
			i__5 = jgl + kbl - 1;
#line 1094 "cgesvj.f"
			i__4 = min(i__5,*n);
#line 1094 "cgesvj.f"
			for (q = jgl; q <= i__4; ++q) {

#line 1096 "cgesvj.f"
			    aaqq = sva[q];
#line 1097 "cgesvj.f"
			    if (aaqq > 0.) {
#line 1098 "cgesvj.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 1104 "cgesvj.f"
				if (aaqq >= 1.) {
#line 1105 "cgesvj.f"
				    if (aapp >= aaqq) {
#line 1106 "cgesvj.f"
					rotok = small * aapp <= aaqq;
#line 1107 "cgesvj.f"
				    } else {
#line 1108 "cgesvj.f"
					rotok = small * aaqq <= aapp;
#line 1109 "cgesvj.f"
				    }
#line 1110 "cgesvj.f"
				    if (aapp < big / aaqq) {
#line 1111 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1111 "cgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 1111 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1111 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1113 "cgesvj.f"
				    } else {
#line 1114 "cgesvj.f"
					ccopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 1116 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aapp, &
						c_b41, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1119 "cgesvj.f"
					cdotc_(&z__2, m, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1119 "cgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 1119 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1121 "cgesvj.f"
				    }
#line 1122 "cgesvj.f"
				} else {
#line 1123 "cgesvj.f"
				    if (aapp >= aaqq) {
#line 1124 "cgesvj.f"
					rotok = aapp <= aaqq / small;
#line 1125 "cgesvj.f"
				    } else {
#line 1126 "cgesvj.f"
					rotok = aaqq <= aapp / small;
#line 1127 "cgesvj.f"
				    }
#line 1128 "cgesvj.f"
				    if (aapp > small / aaqq) {
#line 1129 "cgesvj.f"
					cdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1129 "cgesvj.f"
					d__1 = max(aaqq,aapp);
#line 1129 "cgesvj.f"
					z__2.r = z__3.r / d__1, z__2.i = 
						z__3.i / d__1;
#line 1129 "cgesvj.f"
					d__2 = min(aaqq,aapp);
#line 1129 "cgesvj.f"
					z__1.r = z__2.r / d__2, z__1.i = 
						z__2.i / d__2;
#line 1129 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1132 "cgesvj.f"
				    } else {
#line 1133 "cgesvj.f"
					ccopy_(m, &a[q * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 1135 "cgesvj.f"
					clascl_("G", &c__0, &c__0, &aaqq, &
						c_b41, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1138 "cgesvj.f"
					cdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &cwork[*n + 1], &c__1);
#line 1138 "cgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1138 "cgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1140 "cgesvj.f"
				    }
#line 1141 "cgesvj.f"
				}

/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 1144 "cgesvj.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 1145 "cgesvj.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 1145 "cgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 1149 "cgesvj.f"
				if (abs(aapq1) > tol) {
#line 1150 "cgesvj.f"
				    d__1 = z_abs(&aapq);
#line 1150 "cgesvj.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 1150 "cgesvj.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;
#line 1151 "cgesvj.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 1153 "cgesvj.f"
				    pskipped = 0;
#line 1154 "cgesvj.f"
				    ++iswrot;

#line 1156 "cgesvj.f"
				    if (rotok) {

#line 1158 "cgesvj.f"
					aqoap = aaqq / aapp;
#line 1159 "cgesvj.f"
					apoaq = aapp / aaqq;
#line 1160 "cgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 1161 "cgesvj.f"
					if (aaqq > aapp0) {
#line 1161 "cgesvj.f"
					    theta = -theta;
#line 1161 "cgesvj.f"
					}

#line 1163 "cgesvj.f"
					if (abs(theta) > bigtheta) {
#line 1164 "cgesvj.f"
					    t = .5 / theta;
#line 1165 "cgesvj.f"
					    cs = 1.;
#line 1166 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 1166 "cgesvj.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 1166 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 1168 "cgesvj.f"
					    if (rsvec) {
#line 1169 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 1169 "cgesvj.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1169 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 1171 "cgesvj.f"
					    }
/* Computing MAX */
#line 1172 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 1172 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1174 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 1174 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 1176 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 1176 "cgesvj.f"
					    mxsinj = max(d__1,d__2);
#line 1177 "cgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 1181 "cgesvj.f"
					    thsign = -d_sign(&c_b41, &aapq1);
#line 1182 "cgesvj.f"
					    if (aaqq > aapp0) {
#line 1182 "cgesvj.f"
			  thsign = -thsign;
#line 1182 "cgesvj.f"
					    }
#line 1183 "cgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 1185 "cgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 1186 "cgesvj.f"
					    sn = t * cs;
/* Computing MAX */
#line 1187 "cgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 1187 "cgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 1188 "cgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 1188 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1190 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 1190 "cgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 1193 "cgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 1193 "cgesvj.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 1193 "cgesvj.f"
					    crot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 1195 "cgesvj.f"
					    if (rsvec) {
#line 1196 "cgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 1196 "cgesvj.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 1196 "cgesvj.f"
			  crot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 1198 "cgesvj.f"
					    }
#line 1199 "cgesvj.f"
					}
#line 1200 "cgesvj.f"
					i__5 = p;
#line 1200 "cgesvj.f"
					i__6 = q;
#line 1200 "cgesvj.f"
					z__2.r = -cwork[i__6].r, z__2.i = 
						-cwork[i__6].i;
#line 1200 "cgesvj.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 1200 "cgesvj.f"
					cwork[i__5].r = z__1.r, cwork[i__5].i 
						= z__1.i;

#line 1202 "cgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 1204 "cgesvj.f"
					if (aapp > aaqq) {
#line 1205 "cgesvj.f"
					    ccopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &cwork[*n + 1], &
						    c__1);
#line 1207 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b41, m, &c__1, &cwork[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1210 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b41, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1213 "cgesvj.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 1213 "cgesvj.f"
					    caxpy_(m, &z__1, &cwork[*n + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1);
#line 1215 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &c_b41,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1218 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 1218 "cgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 1220 "cgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1221 "cgesvj.f"
					} else {
#line 1222 "cgesvj.f"
					    ccopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &cwork[*n + 1], &
						    c__1);
#line 1224 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b41, m, &c__1, &cwork[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1227 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &aapp, 
						    &c_b41, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1230 "cgesvj.f"
					    d_cnjg(&z__2, &aapq);
#line 1230 "cgesvj.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 1230 "cgesvj.f"
					    caxpy_(m, &z__1, &cwork[*n + 1], &
						    c__1, &a[p * a_dim1 + 1], 
						    &c__1);
#line 1232 "cgesvj.f"
					    clascl_("G", &c__0, &c__0, &c_b41,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1235 "cgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 1235 "cgesvj.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 1237 "cgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1238 "cgesvj.f"
					}
#line 1239 "cgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 1244 "cgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1244 "cgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1246 "cgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1248 "cgesvj.f"
					    sva[q] = scnrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 1249 "cgesvj.f"
					} else {
#line 1250 "cgesvj.f"
					    t = 0.;
#line 1251 "cgesvj.f"
					    aaqq = 1.;
#line 1252 "cgesvj.f"
					    classq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1254 "cgesvj.f"
					    sva[q] = t * sqrt(aaqq);
#line 1255 "cgesvj.f"
					}
#line 1256 "cgesvj.f"
				    }
/* Computing 2nd power */
#line 1257 "cgesvj.f"
				    d__1 = aapp / aapp0;
#line 1257 "cgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1258 "cgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1260 "cgesvj.f"
					    aapp = scnrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 1261 "cgesvj.f"
					} else {
#line 1262 "cgesvj.f"
					    t = 0.;
#line 1263 "cgesvj.f"
					    aapp = 1.;
#line 1264 "cgesvj.f"
					    classq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1266 "cgesvj.f"
					    aapp = t * sqrt(aapp);
#line 1267 "cgesvj.f"
					}
#line 1268 "cgesvj.f"
					sva[p] = aapp;
#line 1269 "cgesvj.f"
				    }
/*              end of OK rotation */
#line 1271 "cgesvj.f"
				} else {
#line 1272 "cgesvj.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1274 "cgesvj.f"
				    ++pskipped;
#line 1275 "cgesvj.f"
				    ++ijblsk;
#line 1276 "cgesvj.f"
				}
#line 1277 "cgesvj.f"
			    } else {
#line 1278 "cgesvj.f"
				++notrot;
#line 1279 "cgesvj.f"
				++pskipped;
#line 1280 "cgesvj.f"
				++ijblsk;
#line 1281 "cgesvj.f"
			    }

#line 1283 "cgesvj.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 1285 "cgesvj.f"
				sva[p] = aapp;
#line 1286 "cgesvj.f"
				notrot = 0;
#line 1287 "cgesvj.f"
				goto L2011;
#line 1288 "cgesvj.f"
			    }
#line 1289 "cgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1291 "cgesvj.f"
				aapp = -aapp;
#line 1292 "cgesvj.f"
				notrot = 0;
#line 1293 "cgesvj.f"
				goto L2203;
#line 1294 "cgesvj.f"
			    }

#line 1296 "cgesvj.f"
/* L2200: */
#line 1296 "cgesvj.f"
			}
/*        end of the q-loop */
#line 1298 "cgesvj.f"
L2203:

#line 1300 "cgesvj.f"
			sva[p] = aapp;

#line 1302 "cgesvj.f"
		    } else {

#line 1304 "cgesvj.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1304 "cgesvj.f"
			    i__4 = jgl + kbl - 1;
#line 1304 "cgesvj.f"
			    notrot = notrot + min(i__4,*n) - jgl + 1;
#line 1304 "cgesvj.f"
			}
#line 1306 "cgesvj.f"
			if (aapp < 0.) {
#line 1306 "cgesvj.f"
			    notrot = 0;
#line 1306 "cgesvj.f"
			}

#line 1308 "cgesvj.f"
		    }

#line 1310 "cgesvj.f"
/* L2100: */
#line 1310 "cgesvj.f"
		}
/*     end of the p-loop */
#line 1312 "cgesvj.f"
/* L2010: */
#line 1312 "cgesvj.f"
	    }
/*     end of the jbc-loop */
#line 1314 "cgesvj.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1316 "cgesvj.f"
	    i__3 = igl + kbl - 1;
#line 1316 "cgesvj.f"
	    i__2 = min(i__3,*n);
#line 1316 "cgesvj.f"
	    for (p = igl; p <= i__2; ++p) {
#line 1317 "cgesvj.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1318 "cgesvj.f"
/* L2012: */
#line 1318 "cgesvj.f"
	    }
/* ** */
#line 1320 "cgesvj.f"
/* L2000: */
#line 1320 "cgesvj.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1324 "cgesvj.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1326 "cgesvj.f"
	    sva[*n] = scnrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 1327 "cgesvj.f"
	} else {
#line 1328 "cgesvj.f"
	    t = 0.;
#line 1329 "cgesvj.f"
	    aapp = 1.;
#line 1330 "cgesvj.f"
	    classq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1331 "cgesvj.f"
	    sva[*n] = t * sqrt(aapp);
#line 1332 "cgesvj.f"
	}

/*     Additional steering devices */

#line 1336 "cgesvj.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1336 "cgesvj.f"
	    swband = i__;
#line 1336 "cgesvj.f"
	}

#line 1339 "cgesvj.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * tol && (
		doublereal) (*n) * mxaapq * mxsinj < tol) {
#line 1341 "cgesvj.f"
	    goto L1994;
#line 1342 "cgesvj.f"
	}

#line 1344 "cgesvj.f"
	if (notrot >= emptsw) {
#line 1344 "cgesvj.f"
	    goto L1994;
#line 1344 "cgesvj.f"
	}

#line 1346 "cgesvj.f"
/* L1993: */
#line 1346 "cgesvj.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 1350 "cgesvj.f"
    *info = 29;
#line 1351 "cgesvj.f"
    goto L1995;

#line 1353 "cgesvj.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 1357 "cgesvj.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1359 "cgesvj.f"
L1995:

/*     Sort the singular values and find how many are above */
/*     the underflow threshold. */

#line 1364 "cgesvj.f"
    n2 = 0;
#line 1365 "cgesvj.f"
    n4 = 0;
#line 1366 "cgesvj.f"
    i__1 = *n - 1;
#line 1366 "cgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 1367 "cgesvj.f"
	i__2 = *n - p + 1;
#line 1367 "cgesvj.f"
	q = isamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1368 "cgesvj.f"
	if (p != q) {
#line 1369 "cgesvj.f"
	    temp1 = sva[p];
#line 1370 "cgesvj.f"
	    sva[p] = sva[q];
#line 1371 "cgesvj.f"
	    sva[q] = temp1;
#line 1372 "cgesvj.f"
	    cswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1373 "cgesvj.f"
	    if (rsvec) {
#line 1373 "cgesvj.f"
		cswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1373 "cgesvj.f"
	    }
#line 1374 "cgesvj.f"
	}
#line 1375 "cgesvj.f"
	if (sva[p] != 0.) {
#line 1376 "cgesvj.f"
	    ++n4;
#line 1377 "cgesvj.f"
	    if (sva[p] * skl > sfmin) {
#line 1377 "cgesvj.f"
		++n2;
#line 1377 "cgesvj.f"
	    }
#line 1378 "cgesvj.f"
	}
#line 1379 "cgesvj.f"
/* L5991: */
#line 1379 "cgesvj.f"
    }
#line 1380 "cgesvj.f"
    if (sva[*n] != 0.) {
#line 1381 "cgesvj.f"
	++n4;
#line 1382 "cgesvj.f"
	if (sva[*n] * skl > sfmin) {
#line 1382 "cgesvj.f"
	    ++n2;
#line 1382 "cgesvj.f"
	}
#line 1383 "cgesvj.f"
    }

/*     Normalize the left singular vectors. */

#line 1387 "cgesvj.f"
    if (lsvec || uctol) {
#line 1388 "cgesvj.f"
	i__1 = n4;
#line 1388 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
/*           CALL CSSCAL( M, ONE / SVA( p ), A( 1, p ), 1 ) */
#line 1390 "cgesvj.f"
	    clascl_("G", &c__0, &c__0, &sva[p], &c_b41, m, &c__1, &a[p * 
		    a_dim1 + 1], m, &ierr, (ftnlen)1);
#line 1391 "cgesvj.f"
/* L1998: */
#line 1391 "cgesvj.f"
	}
#line 1392 "cgesvj.f"
    }

/*     Scale the product of Jacobi rotations. */

#line 1396 "cgesvj.f"
    if (rsvec) {
#line 1397 "cgesvj.f"
	i__1 = *n;
#line 1397 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1398 "cgesvj.f"
	    temp1 = 1. / scnrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
#line 1399 "cgesvj.f"
	    csscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
#line 1400 "cgesvj.f"
/* L2399: */
#line 1400 "cgesvj.f"
	}
#line 1401 "cgesvj.f"
    }

/*     Undo scaling, if necessary (and possible). */
#line 1404 "cgesvj.f"
    if (skl > 1. && sva[1] < big / skl || skl < 1. && sva[max(n2,1)] > sfmin /
	     skl) {
#line 1407 "cgesvj.f"
	i__1 = *n;
#line 1407 "cgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1408 "cgesvj.f"
	    sva[p] = skl * sva[p];
#line 1409 "cgesvj.f"
/* L2400: */
#line 1409 "cgesvj.f"
	}
#line 1410 "cgesvj.f"
	skl = 1.;
#line 1411 "cgesvj.f"
    }

#line 1413 "cgesvj.f"
    rwork[1] = skl;
/*     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
/*     then some of the singular values may overflow or underflow and */
/*     the spectrum is given in this factored representation. */

#line 1418 "cgesvj.f"
    rwork[2] = (doublereal) n4;
/*     N4 is the number of computed nonzero singular values of A. */

#line 1421 "cgesvj.f"
    rwork[3] = (doublereal) n2;
/*     N2 is the number of singular values of A greater than SFMIN. */
/*     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
/*     that may carry some information. */

#line 1426 "cgesvj.f"
    rwork[4] = (doublereal) i__;
/*     i is the index of the last sweep before declaring convergence. */

#line 1429 "cgesvj.f"
    rwork[5] = mxaapq;
/*     MXAAPQ is the largest absolute value of scaled pivots in the */
/*     last sweep */

#line 1433 "cgesvj.f"
    rwork[6] = mxsinj;
/*     MXSINJ is the largest absolute value of the sines of Jacobi angles */
/*     in the last sweep */

#line 1437 "cgesvj.f"
    return 0;
/*     .. */
/*     .. END OF CGESVJ */
/*     .. */
} /* cgesvj_ */

