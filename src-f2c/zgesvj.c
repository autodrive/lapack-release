#line 1 "zgesvj.f"
/* zgesvj.f -- translated by f2c (version 20100827).
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

#line 1 "zgesvj.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b42 = 1.;
static integer c__2 = 2;

/* > \brief <b> ZGESVJ </b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGESVJ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvj.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvj.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvj.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, */
/*                          LDV, CWORK, LWORK, RWORK, LRWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDV, LWORK, LRWORK, M, MV, N */
/*       CHARACTER*1        JOBA, JOBU, JOBV */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ),  V( LDV, * ), CWORK( LWORK ) */
/*       DOUBLE PRECISION   RWORK( LRWORK ),  SVA( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGESVJ computes the singular value decomposition (SVD) of a complex */
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
/* >          = 'U' or 'F': The left singular vectors corresponding to the nonzero */
/* >                 singular values are computed and returned in the leading */
/* >                 columns of A. See more details in the description of A. */
/* >                 The default numerical orthogonality threshold is set to */
/* >                 approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=DLAMCH('E'). */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          If JOBU .EQ. 'U' .OR. JOBU .EQ. 'C': */
/* >                 If INFO .EQ. 0 : */
/* >                 RANKA orthonormal columns of U are returned in the */
/* >                 leading RANKA columns of the array A. Here RANKA <= N */
/* >                 is the number of computed singular values of A that are */
/* >                 above the underflow threshold DLAMCH('S'). The singular */
/* >                 vectors corresponding to underflowed or zero singular */
/* >                 values are not computed. The value of RANKA is returned */
/* >                 in the array RWORK as RANKA=NINT(RWORK(2)). Also see the */
/* >                 descriptions of SVA and RWORK. The computed columns of U */
/* >                 are mutually numerically orthogonal up to approximately */
/* >                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU.EQ.'C'), */
/* >                 see the description of JOBU. */
/* >                 If INFO .GT. 0, */
/* >                 the procedure ZGESVJ did not converge in the given number */
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
/* >                 the procedure ZGESVJ did not converge in the given number */
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
/* >          the procedure ZGESVJ did not converge in the given number of */
/* >          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate. */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* >          MV is INTEGER */
/* >          If JOBV .EQ. 'A', then the product of Jacobi rotations in ZGESVJ */
/* >          is applied to the first MV rows of V. See the description of JOBV. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDV,N) */
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
/* >          CWORK is COMPLEX*16 array, dimension max(1,LWORK). */
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
/* >          RWORK is DOUBLE PRECISION array, dimension max(6,LRWORK). */
/* >          On entry, */
/* >          If JOBU .EQ. 'C' : */
/* >          RWORK(1) = CTOL, where CTOL defines the threshold for convergence. */
/* >                    The process stops if all columns of A are mutually */
/* >                    orthogonal up to CTOL*EPS, EPS=DLAMCH('E'). */
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
/* >                    This is useful information in cases when ZGESVJ did */
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
/* >          > 0 : ZGESVJ did not converge in the maximal allowed number */
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

/* > \ingroup complex16GEcomputational */

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
/* Subroutine */ int zgesvj_(char *joba, char *jobu, char *jobv, integer *m, 
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
    static doublecomplex ompq;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal aapp0, aapq1, temp1, apoaq, aqoap;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal theta, small, sfmin;
    static logical lsvec;
    static doublereal epsln;
    static logical applv, rsvec, uctol;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical lower, upper, rotok;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zgsvj0_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), zgsvj1_(char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ijblsk, swband;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static integer blskip;
    static doublereal mxaapq;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static doublereal thsign;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static doublereal mxsinj;
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static integer emptsw;
    static logical lquery;
    static integer notrot, iswrot, lkahead;
    static logical goscale, noscale;
    static doublereal rootbig, rooteps;
    static integer rowskip;
    static doublereal roottol;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

#line 419 "zgesvj.f"
    /* Parameter adjustments */
#line 419 "zgesvj.f"
    --sva;
#line 419 "zgesvj.f"
    a_dim1 = *lda;
#line 419 "zgesvj.f"
    a_offset = 1 + a_dim1;
#line 419 "zgesvj.f"
    a -= a_offset;
#line 419 "zgesvj.f"
    v_dim1 = *ldv;
#line 419 "zgesvj.f"
    v_offset = 1 + v_dim1;
#line 419 "zgesvj.f"
    v -= v_offset;
#line 419 "zgesvj.f"
    --cwork;
#line 419 "zgesvj.f"
    --rwork;
#line 419 "zgesvj.f"

#line 419 "zgesvj.f"
    /* Function Body */
#line 419 "zgesvj.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 420 "zgesvj.f"
    uctol = lsame_(jobu, "C", (ftnlen)1, (ftnlen)1);
#line 421 "zgesvj.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || lsame_(jobv, "J", (
	    ftnlen)1, (ftnlen)1);
#line 422 "zgesvj.f"
    applv = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 423 "zgesvj.f"
    upper = lsame_(joba, "U", (ftnlen)1, (ftnlen)1);
#line 424 "zgesvj.f"
    lower = lsame_(joba, "L", (ftnlen)1, (ftnlen)1);

#line 426 "zgesvj.f"
    lquery = *lwork == -1 || *lrwork == -1;
#line 427 "zgesvj.f"
    if (! (upper || lower || lsame_(joba, "G", (ftnlen)1, (ftnlen)1))) {
#line 428 "zgesvj.f"
	*info = -1;
#line 429 "zgesvj.f"
    } else if (! (lsvec || uctol || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 430 "zgesvj.f"
	*info = -2;
#line 431 "zgesvj.f"
    } else if (! (rsvec || applv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 432 "zgesvj.f"
	*info = -3;
#line 433 "zgesvj.f"
    } else if (*m < 0) {
#line 434 "zgesvj.f"
	*info = -4;
#line 435 "zgesvj.f"
    } else if (*n < 0 || *n > *m) {
#line 436 "zgesvj.f"
	*info = -5;
#line 437 "zgesvj.f"
    } else if (*lda < *m) {
#line 438 "zgesvj.f"
	*info = -7;
#line 439 "zgesvj.f"
    } else if (*mv < 0) {
#line 440 "zgesvj.f"
	*info = -9;
#line 441 "zgesvj.f"
    } else if (rsvec && *ldv < *n || applv && *ldv < *mv) {
#line 443 "zgesvj.f"
	*info = -11;
#line 444 "zgesvj.f"
    } else if (uctol && rwork[1] <= 1.) {
#line 445 "zgesvj.f"
	*info = -12;
#line 446 "zgesvj.f"
    } else if (*lwork < *m + *n && ! lquery) {
#line 447 "zgesvj.f"
	*info = -13;
#line 448 "zgesvj.f"
    } else if (*lrwork < max(*n,6) && ! lquery) {
#line 449 "zgesvj.f"
	*info = -15;
#line 450 "zgesvj.f"
    } else {
#line 451 "zgesvj.f"
	*info = 0;
#line 452 "zgesvj.f"
    }

/*     #:( */
#line 455 "zgesvj.f"
    if (*info != 0) {
#line 456 "zgesvj.f"
	i__1 = -(*info);
#line 456 "zgesvj.f"
	xerbla_("ZGESVJ", &i__1, (ftnlen)6);
#line 457 "zgesvj.f"
	return 0;
#line 458 "zgesvj.f"
    } else if (lquery) {
#line 459 "zgesvj.f"
	i__1 = *m + *n;
#line 459 "zgesvj.f"
	cwork[1].r = (doublereal) i__1, cwork[1].i = 0.;
#line 460 "zgesvj.f"
	rwork[1] = (doublereal) max(*n,6);
#line 461 "zgesvj.f"
	return 0;
#line 462 "zgesvj.f"
    }

/* #:) Quick return for void matrix */

#line 466 "zgesvj.f"
    if (*m == 0 || *n == 0) {
#line 466 "zgesvj.f"
	return 0;
#line 466 "zgesvj.f"
    }

/*     Set numerical parameters */
/*     The stopping criterion for Jacobi rotations is */

/*     max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS */

/*     where EPS is the round-off and CTOL is defined as follows: */

#line 475 "zgesvj.f"
    if (uctol) {
/*        ... user controlled */
#line 477 "zgesvj.f"
	ctol = rwork[1];
#line 478 "zgesvj.f"
    } else {
/*        ... default */
#line 480 "zgesvj.f"
	if (lsvec || rsvec || applv) {
#line 481 "zgesvj.f"
	    ctol = sqrt((doublereal) (*m));
#line 482 "zgesvj.f"
	} else {
#line 483 "zgesvj.f"
	    ctol = (doublereal) (*m);
#line 484 "zgesvj.f"
	}
#line 485 "zgesvj.f"
    }
/*     ... and the machine dependent parameters are */
/* [!]  (Make sure that SLAMCH() works properly on the target machine.) */

#line 489 "zgesvj.f"
    epsln = dlamch_("Epsilon", (ftnlen)7);
#line 490 "zgesvj.f"
    rooteps = sqrt(epsln);
#line 491 "zgesvj.f"
    sfmin = dlamch_("SafeMinimum", (ftnlen)11);
#line 492 "zgesvj.f"
    rootsfmin = sqrt(sfmin);
#line 493 "zgesvj.f"
    small = sfmin / epsln;
#line 494 "zgesvj.f"
    big = dlamch_("Overflow", (ftnlen)8);
/*     BIG         = ONE    / SFMIN */
#line 496 "zgesvj.f"
    rootbig = 1. / rootsfmin;
/*      LARGE = BIG / SQRT( DBLE( M*N ) ) */
#line 498 "zgesvj.f"
    bigtheta = 1. / rooteps;

#line 500 "zgesvj.f"
    tol = ctol * epsln;
#line 501 "zgesvj.f"
    roottol = sqrt(tol);

#line 503 "zgesvj.f"
    if ((doublereal) (*m) * epsln >= 1.) {
#line 504 "zgesvj.f"
	*info = -4;
#line 505 "zgesvj.f"
	i__1 = -(*info);
#line 505 "zgesvj.f"
	xerbla_("ZGESVJ", &i__1, (ftnlen)6);
#line 506 "zgesvj.f"
	return 0;
#line 507 "zgesvj.f"
    }

/*     Initialize the right singular vector matrix. */

#line 511 "zgesvj.f"
    if (rsvec) {
#line 512 "zgesvj.f"
	mvl = *n;
#line 513 "zgesvj.f"
	zlaset_("A", &mvl, n, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 514 "zgesvj.f"
    } else if (applv) {
#line 515 "zgesvj.f"
	mvl = *mv;
#line 516 "zgesvj.f"
    }
#line 517 "zgesvj.f"
    rsvec = rsvec || applv;

/*     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N ) */
/* (!)  If necessary, scale A to protect the largest singular value */
/*     from overflow. It is possible that saving the largest singular */
/*     value destroys the information about the small ones. */
/*     This initial scaling is almost minimal in the sense that the */
/*     goal is to make sure that no column norm overflows, and that */
/*     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries */
/*     in A are detected, the procedure returns with INFO=-6. */

#line 528 "zgesvj.f"
    skl = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 529 "zgesvj.f"
    noscale = TRUE_;
#line 530 "zgesvj.f"
    goscale = TRUE_;

#line 532 "zgesvj.f"
    if (lower) {
/*        the input matrix is M-by-N lower triangular (trapezoidal) */
#line 534 "zgesvj.f"
	i__1 = *n;
#line 534 "zgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 535 "zgesvj.f"
	    aapp = 0.;
#line 536 "zgesvj.f"
	    aaqq = 1.;
#line 537 "zgesvj.f"
	    i__2 = *m - p + 1;
#line 537 "zgesvj.f"
	    zlassq_(&i__2, &a[p + p * a_dim1], &c__1, &aapp, &aaqq);
#line 538 "zgesvj.f"
	    if (aapp > big) {
#line 539 "zgesvj.f"
		*info = -6;
#line 540 "zgesvj.f"
		i__2 = -(*info);
#line 540 "zgesvj.f"
		xerbla_("ZGESVJ", &i__2, (ftnlen)6);
#line 541 "zgesvj.f"
		return 0;
#line 542 "zgesvj.f"
	    }
#line 543 "zgesvj.f"
	    aaqq = sqrt(aaqq);
#line 544 "zgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 545 "zgesvj.f"
		sva[p] = aapp * aaqq;
#line 546 "zgesvj.f"
	    } else {
#line 547 "zgesvj.f"
		noscale = FALSE_;
#line 548 "zgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 549 "zgesvj.f"
		if (goscale) {
#line 550 "zgesvj.f"
		    goscale = FALSE_;
#line 551 "zgesvj.f"
		    i__2 = p - 1;
#line 551 "zgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 552 "zgesvj.f"
			sva[q] *= skl;
#line 553 "zgesvj.f"
/* L1873: */
#line 553 "zgesvj.f"
		    }
#line 554 "zgesvj.f"
		}
#line 555 "zgesvj.f"
	    }
#line 556 "zgesvj.f"
/* L1874: */
#line 556 "zgesvj.f"
	}
#line 557 "zgesvj.f"
    } else if (upper) {
/*        the input matrix is M-by-N upper triangular (trapezoidal) */
#line 559 "zgesvj.f"
	i__1 = *n;
#line 559 "zgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 560 "zgesvj.f"
	    aapp = 0.;
#line 561 "zgesvj.f"
	    aaqq = 1.;
#line 562 "zgesvj.f"
	    zlassq_(&p, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 563 "zgesvj.f"
	    if (aapp > big) {
#line 564 "zgesvj.f"
		*info = -6;
#line 565 "zgesvj.f"
		i__2 = -(*info);
#line 565 "zgesvj.f"
		xerbla_("ZGESVJ", &i__2, (ftnlen)6);
#line 566 "zgesvj.f"
		return 0;
#line 567 "zgesvj.f"
	    }
#line 568 "zgesvj.f"
	    aaqq = sqrt(aaqq);
#line 569 "zgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 570 "zgesvj.f"
		sva[p] = aapp * aaqq;
#line 571 "zgesvj.f"
	    } else {
#line 572 "zgesvj.f"
		noscale = FALSE_;
#line 573 "zgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 574 "zgesvj.f"
		if (goscale) {
#line 575 "zgesvj.f"
		    goscale = FALSE_;
#line 576 "zgesvj.f"
		    i__2 = p - 1;
#line 576 "zgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 577 "zgesvj.f"
			sva[q] *= skl;
#line 578 "zgesvj.f"
/* L2873: */
#line 578 "zgesvj.f"
		    }
#line 579 "zgesvj.f"
		}
#line 580 "zgesvj.f"
	    }
#line 581 "zgesvj.f"
/* L2874: */
#line 581 "zgesvj.f"
	}
#line 582 "zgesvj.f"
    } else {
/*        the input matrix is M-by-N general dense */
#line 584 "zgesvj.f"
	i__1 = *n;
#line 584 "zgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 585 "zgesvj.f"
	    aapp = 0.;
#line 586 "zgesvj.f"
	    aaqq = 1.;
#line 587 "zgesvj.f"
	    zlassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 588 "zgesvj.f"
	    if (aapp > big) {
#line 589 "zgesvj.f"
		*info = -6;
#line 590 "zgesvj.f"
		i__2 = -(*info);
#line 590 "zgesvj.f"
		xerbla_("ZGESVJ", &i__2, (ftnlen)6);
#line 591 "zgesvj.f"
		return 0;
#line 592 "zgesvj.f"
	    }
#line 593 "zgesvj.f"
	    aaqq = sqrt(aaqq);
#line 594 "zgesvj.f"
	    if (aapp < big / aaqq && noscale) {
#line 595 "zgesvj.f"
		sva[p] = aapp * aaqq;
#line 596 "zgesvj.f"
	    } else {
#line 597 "zgesvj.f"
		noscale = FALSE_;
#line 598 "zgesvj.f"
		sva[p] = aapp * (aaqq * skl);
#line 599 "zgesvj.f"
		if (goscale) {
#line 600 "zgesvj.f"
		    goscale = FALSE_;
#line 601 "zgesvj.f"
		    i__2 = p - 1;
#line 601 "zgesvj.f"
		    for (q = 1; q <= i__2; ++q) {
#line 602 "zgesvj.f"
			sva[q] *= skl;
#line 603 "zgesvj.f"
/* L3873: */
#line 603 "zgesvj.f"
		    }
#line 604 "zgesvj.f"
		}
#line 605 "zgesvj.f"
	    }
#line 606 "zgesvj.f"
/* L3874: */
#line 606 "zgesvj.f"
	}
#line 607 "zgesvj.f"
    }

#line 609 "zgesvj.f"
    if (noscale) {
#line 609 "zgesvj.f"
	skl = 1.;
#line 609 "zgesvj.f"
    }

/*     Move the smaller part of the spectrum from the underflow threshold */
/* (!)  Start by determining the position of the nonzero entries of the */
/*     array SVA() relative to ( SFMIN, BIG ). */

#line 615 "zgesvj.f"
    aapp = 0.;
#line 616 "zgesvj.f"
    aaqq = big;
#line 617 "zgesvj.f"
    i__1 = *n;
#line 617 "zgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 618 "zgesvj.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 618 "zgesvj.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 618 "zgesvj.f"
	    aaqq = min(d__1,d__2);
#line 618 "zgesvj.f"
	}
/* Computing MAX */
#line 619 "zgesvj.f"
	d__1 = aapp, d__2 = sva[p];
#line 619 "zgesvj.f"
	aapp = max(d__1,d__2);
#line 620 "zgesvj.f"
/* L4781: */
#line 620 "zgesvj.f"
    }

/* #:) Quick return for zero matrix */

#line 624 "zgesvj.f"
    if (aapp == 0.) {
#line 625 "zgesvj.f"
	if (lsvec) {
#line 625 "zgesvj.f"
	    zlaset_("G", m, n, &c_b1, &c_b2, &a[a_offset], lda, (ftnlen)1);
#line 625 "zgesvj.f"
	}
#line 626 "zgesvj.f"
	rwork[1] = 1.;
#line 627 "zgesvj.f"
	rwork[2] = 0.;
#line 628 "zgesvj.f"
	rwork[3] = 0.;
#line 629 "zgesvj.f"
	rwork[4] = 0.;
#line 630 "zgesvj.f"
	rwork[5] = 0.;
#line 631 "zgesvj.f"
	rwork[6] = 0.;
#line 632 "zgesvj.f"
	return 0;
#line 633 "zgesvj.f"
    }

/* #:) Quick return for one-column matrix */

#line 637 "zgesvj.f"
    if (*n == 1) {
#line 638 "zgesvj.f"
	if (lsvec) {
#line 638 "zgesvj.f"
	    zlascl_("G", &c__0, &c__0, &sva[1], &skl, m, &c__1, &a[a_dim1 + 1]
		    , lda, &ierr, (ftnlen)1);
#line 638 "zgesvj.f"
	}
#line 640 "zgesvj.f"
	rwork[1] = 1. / skl;
#line 641 "zgesvj.f"
	if (sva[1] >= sfmin) {
#line 642 "zgesvj.f"
	    rwork[2] = 1.;
#line 643 "zgesvj.f"
	} else {
#line 644 "zgesvj.f"
	    rwork[2] = 0.;
#line 645 "zgesvj.f"
	}
#line 646 "zgesvj.f"
	rwork[3] = 0.;
#line 647 "zgesvj.f"
	rwork[4] = 0.;
#line 648 "zgesvj.f"
	rwork[5] = 0.;
#line 649 "zgesvj.f"
	rwork[6] = 0.;
#line 650 "zgesvj.f"
	return 0;
#line 651 "zgesvj.f"
    }

/*     Protect small singular values from underflow, and try to */
/*     avoid underflows/overflows in computing Jacobi rotations. */

#line 656 "zgesvj.f"
    sn = sqrt(sfmin / epsln);
#line 657 "zgesvj.f"
    temp1 = sqrt(big / (doublereal) (*n));
#line 658 "zgesvj.f"
    if (aapp <= sn || aaqq >= temp1 || sn <= aaqq && aapp <= temp1) {
/* Computing MIN */
#line 660 "zgesvj.f"
	d__1 = big, d__2 = temp1 / aapp;
#line 660 "zgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 663 "zgesvj.f"
    } else if (aaqq <= sn && aapp <= temp1) {
/* Computing MIN */
#line 664 "zgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (aapp * sqrt((doublereal) (*n)));
#line 664 "zgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 667 "zgesvj.f"
    } else if (aaqq >= sn && aapp >= temp1) {
/* Computing MAX */
#line 668 "zgesvj.f"
	d__1 = sn / aaqq, d__2 = temp1 / aapp;
#line 668 "zgesvj.f"
	temp1 = max(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 671 "zgesvj.f"
    } else if (aaqq <= sn && aapp >= temp1) {
/* Computing MIN */
#line 672 "zgesvj.f"
	d__1 = sn / aaqq, d__2 = big / (sqrt((doublereal) (*n)) * aapp);
#line 672 "zgesvj.f"
	temp1 = min(d__1,d__2);
/*         AAQQ  = AAQQ*TEMP1 */
/*         AAPP  = AAPP*TEMP1 */
#line 675 "zgesvj.f"
    } else {
#line 676 "zgesvj.f"
	temp1 = 1.;
#line 677 "zgesvj.f"
    }

/*     Scale, if necessary */

#line 681 "zgesvj.f"
    if (temp1 != 1.) {
#line 682 "zgesvj.f"
	dlascl_("G", &c__0, &c__0, &c_b42, &temp1, n, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 683 "zgesvj.f"
    }
#line 684 "zgesvj.f"
    skl = temp1 * skl;
#line 685 "zgesvj.f"
    if (skl != 1.) {
#line 686 "zgesvj.f"
	zlascl_(joba, &c__0, &c__0, &c_b42, &skl, m, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 687 "zgesvj.f"
	skl = 1. / skl;
#line 688 "zgesvj.f"
    }

/*     Row-cyclic Jacobi SVD algorithm with column pivoting */

#line 692 "zgesvj.f"
    emptsw = *n * (*n - 1) / 2;
#line 693 "zgesvj.f"
    notrot = 0;
#line 695 "zgesvj.f"
    i__1 = *n;
#line 695 "zgesvj.f"
    for (q = 1; q <= i__1; ++q) {
#line 696 "zgesvj.f"
	i__2 = q;
#line 696 "zgesvj.f"
	cwork[i__2].r = 1., cwork[i__2].i = 0.;
#line 697 "zgesvj.f"
/* L1868: */
#line 697 "zgesvj.f"
    }



#line 701 "zgesvj.f"
    swband = 3;
/* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
/*     if ZGESVJ is used as a computational routine in the preconditioned */
/*     Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure */
/*     works on pivots inside a band-like region around the diagonal. */
/*     The boundaries are determined dynamically, based on the number of */
/*     pivots above a threshold. */

#line 709 "zgesvj.f"
    kbl = min(8,*n);
/* [TP] KBL is a tuning parameter that defines the tile size in the */
/*     tiling of the p-q loops of pivot pairs. In general, an optimal */
/*     value of KBL depends on the matrix dimensions and on the */
/*     parameters of the computer's memory. */

#line 715 "zgesvj.f"
    nbl = *n / kbl;
#line 716 "zgesvj.f"
    if (nbl * kbl != *n) {
#line 716 "zgesvj.f"
	++nbl;
#line 716 "zgesvj.f"
    }

/* Computing 2nd power */
#line 718 "zgesvj.f"
    i__1 = kbl;
#line 718 "zgesvj.f"
    blskip = i__1 * i__1;
/* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */

#line 721 "zgesvj.f"
    rowskip = min(5,kbl);
/* [TP] ROWSKIP is a tuning parameter. */

#line 724 "zgesvj.f"
    lkahead = 1;
/* [TP] LKAHEAD is a tuning parameter. */

/*     Quasi block transformations, using the lower (upper) triangular */
/*     structure of the input matrix. The quasi-block-cycling usually */
/*     invokes cubic convergence. Big part of this cycle is done inside */
/*     canonical subspaces of dimensions less than M. */

/* Computing MAX */
#line 732 "zgesvj.f"
    i__1 = 64, i__2 = kbl << 2;
#line 732 "zgesvj.f"
    if ((lower || upper) && *n > max(i__1,i__2)) {
/* [TP] The number of partition levels and the actual partition are */
/*     tuning parameters. */
#line 735 "zgesvj.f"
	n4 = *n / 4;
#line 736 "zgesvj.f"
	n2 = *n / 2;
#line 737 "zgesvj.f"
	n34 = n4 * 3;
#line 738 "zgesvj.f"
	if (applv) {
#line 739 "zgesvj.f"
	    q = 0;
#line 740 "zgesvj.f"
	} else {
#line 741 "zgesvj.f"
	    q = 1;
#line 742 "zgesvj.f"
	}

#line 744 "zgesvj.f"
	if (lower) {

/*     This works very well on lower triangular matrices, in particular */
/*     in the framework of the preconditioned Jacobi SVD (xGEJSV). */
/*     The idea is simple: */
/*     [+ 0 0 0]   Note that Jacobi transformations of [0 0] */
/*     [+ + 0 0]                                       [0 0] */
/*     [+ + x 0]   actually work on [x 0]              [x 0] */
/*     [+ + x x]                    [x x].             [x x] */

#line 754 "zgesvj.f"
	    i__1 = *m - n34;
#line 754 "zgesvj.f"
	    i__2 = *n - n34;
#line 754 "zgesvj.f"
	    i__3 = *lwork - *n;
#line 754 "zgesvj.f"
	    zgsvj0_(jobv, &i__1, &i__2, &a[n34 + 1 + (n34 + 1) * a_dim1], lda,
		     &cwork[n34 + 1], &sva[n34 + 1], &mvl, &v[n34 * q + 1 + (
		    n34 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &
		    cwork[*n + 1], &i__3, &ierr, (ftnlen)1);
#line 759 "zgesvj.f"
	    i__1 = *m - n2;
#line 759 "zgesvj.f"
	    i__2 = n34 - n2;
#line 759 "zgesvj.f"
	    i__3 = *lwork - *n;
#line 759 "zgesvj.f"
	    zgsvj0_(jobv, &i__1, &i__2, &a[n2 + 1 + (n2 + 1) * a_dim1], lda, &
		    cwork[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 
		    1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__2, &cwork[*n 
		    + 1], &i__3, &ierr, (ftnlen)1);
#line 764 "zgesvj.f"
	    i__1 = *m - n2;
#line 764 "zgesvj.f"
	    i__2 = *n - n2;
#line 764 "zgesvj.f"
	    i__3 = *lwork - *n;
#line 764 "zgesvj.f"
	    zgsvj1_(jobv, &i__1, &i__2, &n4, &a[n2 + 1 + (n2 + 1) * a_dim1], 
		    lda, &cwork[n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (
		    n2 + 1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__3, &ierr, (ftnlen)1);
#line 769 "zgesvj.f"
	    i__1 = *m - n4;
#line 769 "zgesvj.f"
	    i__2 = n2 - n4;
#line 769 "zgesvj.f"
	    i__3 = *lwork - *n;
#line 769 "zgesvj.f"
	    zgsvj0_(jobv, &i__1, &i__2, &a[n4 + 1 + (n4 + 1) * a_dim1], lda, &
		    cwork[n4 + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 
		    1) * v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n 
		    + 1], &i__3, &ierr, (ftnlen)1);

#line 774 "zgesvj.f"
	    i__1 = *lwork - *n;
#line 774 "zgesvj.f"
	    zgsvj0_(jobv, m, &n4, &a[a_offset], lda, &cwork[1], &sva[1], &mvl,
		     &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*
		    n + 1], &i__1, &ierr, (ftnlen)1);

#line 778 "zgesvj.f"
	    i__1 = *lwork - *n;
#line 778 "zgesvj.f"
	    zgsvj1_(jobv, m, &n2, &n4, &a[a_offset], lda, &cwork[1], &sva[1], 
		    &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);


#line 783 "zgesvj.f"
	} else if (upper) {


#line 786 "zgesvj.f"
	    i__1 = *lwork - *n;
#line 786 "zgesvj.f"
	    zgsvj0_(jobv, &n4, &n4, &a[a_offset], lda, &cwork[1], &sva[1], &
		    mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__2, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 790 "zgesvj.f"
	    i__1 = *lwork - *n;
#line 790 "zgesvj.f"
	    zgsvj0_(jobv, &n2, &n4, &a[(n4 + 1) * a_dim1 + 1], lda, &cwork[n4 
		    + 1], &sva[n4 + 1], &mvl, &v[n4 * q + 1 + (n4 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n + 1],
		     &i__1, &ierr, (ftnlen)1);

#line 795 "zgesvj.f"
	    i__1 = *lwork - *n;
#line 795 "zgesvj.f"
	    zgsvj1_(jobv, &n2, &n2, &n4, &a[a_offset], lda, &cwork[1], &sva[1]
		    , &mvl, &v[v_offset], ldv, &epsln, &sfmin, &tol, &c__1, &
		    cwork[*n + 1], &i__1, &ierr, (ftnlen)1);

#line 799 "zgesvj.f"
	    i__1 = n2 + n4;
#line 799 "zgesvj.f"
	    i__2 = *lwork - *n;
#line 799 "zgesvj.f"
	    zgsvj0_(jobv, &i__1, &n4, &a[(n2 + 1) * a_dim1 + 1], lda, &cwork[
		    n2 + 1], &sva[n2 + 1], &mvl, &v[n2 * q + 1 + (n2 + 1) * 
		    v_dim1], ldv, &epsln, &sfmin, &tol, &c__1, &cwork[*n + 1],
		     &i__2, &ierr, (ftnlen)1);
#line 804 "zgesvj.f"
	}

#line 806 "zgesvj.f"
    }

/*     .. Row-cyclic pivot strategy with de Rijk's pivoting .. */

#line 810 "zgesvj.f"
    for (i__ = 1; i__ <= 30; ++i__) {

/*     .. go go go ... */

#line 814 "zgesvj.f"
	mxaapq = 0.;
#line 815 "zgesvj.f"
	mxsinj = 0.;
#line 816 "zgesvj.f"
	iswrot = 0;

#line 818 "zgesvj.f"
	notrot = 0;
#line 819 "zgesvj.f"
	pskipped = 0;

/*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
/*     1 <= p < q <= N. This is the first step toward a blocked implementation */
/*     of the rotations. New implementation, based on block transformations, */
/*     is under development. */

#line 826 "zgesvj.f"
	i__1 = nbl;
#line 826 "zgesvj.f"
	for (ibr = 1; ibr <= i__1; ++ibr) {

#line 828 "zgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

/* Computing MIN */
#line 830 "zgesvj.f"
	    i__3 = lkahead, i__4 = nbl - ibr;
#line 830 "zgesvj.f"
	    i__2 = min(i__3,i__4);
#line 830 "zgesvj.f"
	    for (ir1 = 0; ir1 <= i__2; ++ir1) {

#line 832 "zgesvj.f"
		igl += ir1 * kbl;

/* Computing MIN */
#line 834 "zgesvj.f"
		i__4 = igl + kbl - 1, i__5 = *n - 1;
#line 834 "zgesvj.f"
		i__3 = min(i__4,i__5);
#line 834 "zgesvj.f"
		for (p = igl; p <= i__3; ++p) {

/*     .. de Rijk's pivoting */

#line 838 "zgesvj.f"
		    i__4 = *n - p + 1;
#line 838 "zgesvj.f"
		    q = idamax_(&i__4, &sva[p], &c__1) + p - 1;
#line 839 "zgesvj.f"
		    if (p != q) {
#line 840 "zgesvj.f"
			zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 
				1], &c__1);
#line 841 "zgesvj.f"
			if (rsvec) {
#line 841 "zgesvj.f"
			    zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				    v_dim1 + 1], &c__1);
#line 841 "zgesvj.f"
			}
#line 843 "zgesvj.f"
			temp1 = sva[p];
#line 844 "zgesvj.f"
			sva[p] = sva[q];
#line 845 "zgesvj.f"
			sva[q] = temp1;
#line 846 "zgesvj.f"
			i__4 = p;
#line 846 "zgesvj.f"
			aapq.r = cwork[i__4].r, aapq.i = cwork[i__4].i;
#line 847 "zgesvj.f"
			i__4 = p;
#line 847 "zgesvj.f"
			i__5 = q;
#line 847 "zgesvj.f"
			cwork[i__4].r = cwork[i__5].r, cwork[i__4].i = cwork[
				i__5].i;
#line 848 "zgesvj.f"
			i__4 = q;
#line 848 "zgesvj.f"
			cwork[i__4].r = aapq.r, cwork[i__4].i = aapq.i;
#line 849 "zgesvj.f"
		    }

#line 851 "zgesvj.f"
		    if (ir1 == 0) {

/*        Column norms are periodically updated by explicit */
/*        norm computation. */
/* [!]     Caveat: */
/*        Unfortunately, some BLAS implementations compute DZNRM2(M,A(1,p),1) */
/*        as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to */
/*        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
/*        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
/*        Hence, DZNRM2 cannot be trusted, not even in the case when */
/*        the true norm is far from the under(over)flow boundaries. */
/*        If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF */
/*        below should be replaced with "AAPP = DZNRM2( M, A(1,p), 1 )". */

#line 865 "zgesvj.f"
			if (sva[p] < rootbig && sva[p] > rootsfmin) {
#line 867 "zgesvj.f"
			    sva[p] = dznrm2_(m, &a[p * a_dim1 + 1], &c__1);
#line 868 "zgesvj.f"
			} else {
#line 869 "zgesvj.f"
			    temp1 = 0.;
#line 870 "zgesvj.f"
			    aapp = 1.;
#line 871 "zgesvj.f"
			    zlassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, &
				    aapp);
#line 872 "zgesvj.f"
			    sva[p] = temp1 * sqrt(aapp);
#line 873 "zgesvj.f"
			}
#line 874 "zgesvj.f"
			aapp = sva[p];
#line 875 "zgesvj.f"
		    } else {
#line 876 "zgesvj.f"
			aapp = sva[p];
#line 877 "zgesvj.f"
		    }

#line 879 "zgesvj.f"
		    if (aapp > 0.) {

#line 881 "zgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 883 "zgesvj.f"
			i__5 = igl + kbl - 1;
#line 883 "zgesvj.f"
			i__4 = min(i__5,*n);
#line 883 "zgesvj.f"
			for (q = p + 1; q <= i__4; ++q) {

#line 885 "zgesvj.f"
			    aaqq = sva[q];

#line 887 "zgesvj.f"
			    if (aaqq > 0.) {

#line 889 "zgesvj.f"
				aapp0 = aapp;
#line 890 "zgesvj.f"
				if (aaqq >= 1.) {
#line 891 "zgesvj.f"
				    rotok = small * aapp <= aaqq;
#line 892 "zgesvj.f"
				    if (aapp < big / aaqq) {
#line 893 "zgesvj.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 893 "zgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 893 "zgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 893 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 895 "zgesvj.f"
				    } else {
#line 896 "zgesvj.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 898 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b42, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 900 "zgesvj.f"
					zdotc_(&z__2, m, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 900 "zgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 900 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 902 "zgesvj.f"
				    }
#line 903 "zgesvj.f"
				} else {
#line 904 "zgesvj.f"
				    rotok = aapp <= aaqq / small;
#line 905 "zgesvj.f"
				    if (aapp > small / aaqq) {
#line 906 "zgesvj.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 906 "zgesvj.f"
					z__2.r = z__3.r / aapp, z__2.i = 
						z__3.i / aapp;
#line 906 "zgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 906 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 908 "zgesvj.f"
				    } else {
#line 909 "zgesvj.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 911 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b42, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 914 "zgesvj.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &cwork[*n + 1], &c__1);
#line 914 "zgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 914 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 916 "zgesvj.f"
				    }
#line 917 "zgesvj.f"
				}

/*                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q) */
#line 921 "zgesvj.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 922 "zgesvj.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 922 "zgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 926 "zgesvj.f"
				if (abs(aapq1) > tol) {
#line 927 "zgesvj.f"
				    d__1 = z_abs(&aapq);
#line 927 "zgesvj.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 927 "zgesvj.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;

/*           .. rotate */
/* [RTD]      ROTATED = ROTATED + ONE */

#line 932 "zgesvj.f"
				    if (ir1 == 0) {
#line 933 "zgesvj.f"
					notrot = 0;
#line 934 "zgesvj.f"
					pskipped = 0;
#line 935 "zgesvj.f"
					++iswrot;
#line 936 "zgesvj.f"
				    }

#line 938 "zgesvj.f"
				    if (rotok) {

#line 940 "zgesvj.f"
					aqoap = aaqq / aapp;
#line 941 "zgesvj.f"
					apoaq = aapp / aaqq;
#line 942 "zgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;

#line 944 "zgesvj.f"
					if (abs(theta) > bigtheta) {

#line 946 "zgesvj.f"
					    t = .5 / theta;
#line 947 "zgesvj.f"
					    cs = 1.;
#line 949 "zgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 949 "zgesvj.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 949 "zgesvj.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 951 "zgesvj.f"
					    if (rsvec) {
#line 952 "zgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 952 "zgesvj.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 952 "zgesvj.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 954 "zgesvj.f"
					    }
/* Computing MAX */
#line 956 "zgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 956 "zgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 958 "zgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 958 "zgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 960 "zgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 960 "zgesvj.f"
					    mxsinj = max(d__1,d__2);

#line 962 "zgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 966 "zgesvj.f"
					    thsign = -d_sign(&c_b42, &aapq1);
#line 967 "zgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 969 "zgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 970 "zgesvj.f"
					    sn = t * cs;

/* Computing MAX */
#line 972 "zgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 972 "zgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 973 "zgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 973 "zgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 975 "zgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 975 "zgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 978 "zgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 978 "zgesvj.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 978 "zgesvj.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 980 "zgesvj.f"
					    if (rsvec) {
#line 981 "zgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 981 "zgesvj.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 981 "zgesvj.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 983 "zgesvj.f"
					    }
#line 984 "zgesvj.f"
					}
#line 985 "zgesvj.f"
					i__5 = p;
#line 985 "zgesvj.f"
					i__6 = q;
#line 985 "zgesvj.f"
					z__2.r = -cwork[i__6].r, z__2.i = 
						-cwork[i__6].i;
#line 985 "zgesvj.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 985 "zgesvj.f"
					cwork[i__5].r = z__1.r, cwork[i__5].i 
						= z__1.i;

#line 987 "zgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 989 "zgesvj.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 991 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b42, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 994 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b42, m, &c__1, &a[q * 
						a_dim1 + 1], lda, &ierr, (
						ftnlen)1);
#line 996 "zgesvj.f"
					z__1.r = -aapq.r, z__1.i = -aapq.i;
#line 996 "zgesvj.f"
					zaxpy_(m, &z__1, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 998 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &c_b42, &
						aaqq, m, &c__1, &a[q * a_dim1 
						+ 1], lda, &ierr, (ftnlen)1);
/* Computing MAX */
#line 1000 "zgesvj.f"
					d__1 = 0., d__2 = 1. - aapq1 * aapq1;
#line 1000 "zgesvj.f"
					sva[q] = aaqq * sqrt((max(d__1,d__2)))
						;
#line 1002 "zgesvj.f"
					mxsinj = max(mxsinj,sfmin);
#line 1003 "zgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           recompute SVA(q), SVA(p). */

/* Computing 2nd power */
#line 1009 "zgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1009 "zgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1011 "zgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1013 "zgesvj.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 1014 "zgesvj.f"
					} else {
#line 1015 "zgesvj.f"
					    t = 0.;
#line 1016 "zgesvj.f"
					    aaqq = 1.;
#line 1017 "zgesvj.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1019 "zgesvj.f"
					    sva[q] = t * sqrt(aaqq);
#line 1020 "zgesvj.f"
					}
#line 1021 "zgesvj.f"
				    }
#line 1022 "zgesvj.f"
				    if (aapp / aapp0 <= rooteps) {
#line 1023 "zgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1025 "zgesvj.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 1026 "zgesvj.f"
					} else {
#line 1027 "zgesvj.f"
					    t = 0.;
#line 1028 "zgesvj.f"
					    aapp = 1.;
#line 1029 "zgesvj.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1031 "zgesvj.f"
					    aapp = t * sqrt(aapp);
#line 1032 "zgesvj.f"
					}
#line 1033 "zgesvj.f"
					sva[p] = aapp;
#line 1034 "zgesvj.f"
				    }

#line 1036 "zgesvj.f"
				} else {
/*                             A(:,p) and A(:,q) already numerically orthogonal */
#line 1038 "zgesvj.f"
				    if (ir1 == 0) {
#line 1038 "zgesvj.f"
					++notrot;
#line 1038 "zgesvj.f"
				    }
/* [RTD]      SKIPPED  = SKIPPED + 1 */
#line 1040 "zgesvj.f"
				    ++pskipped;
#line 1041 "zgesvj.f"
				}
#line 1042 "zgesvj.f"
			    } else {
/*                          A(:,q) is zero column */
#line 1044 "zgesvj.f"
				if (ir1 == 0) {
#line 1044 "zgesvj.f"
				    ++notrot;
#line 1044 "zgesvj.f"
				}
#line 1045 "zgesvj.f"
				++pskipped;
#line 1046 "zgesvj.f"
			    }

#line 1048 "zgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1050 "zgesvj.f"
				if (ir1 == 0) {
#line 1050 "zgesvj.f"
				    aapp = -aapp;
#line 1050 "zgesvj.f"
				}
#line 1051 "zgesvj.f"
				notrot = 0;
#line 1052 "zgesvj.f"
				goto L2103;
#line 1053 "zgesvj.f"
			    }

#line 1055 "zgesvj.f"
/* L2002: */
#line 1055 "zgesvj.f"
			}
/*     END q-LOOP */

#line 1058 "zgesvj.f"
L2103:
/*     bailed out of q-loop */

#line 1061 "zgesvj.f"
			sva[p] = aapp;

#line 1063 "zgesvj.f"
		    } else {
#line 1064 "zgesvj.f"
			sva[p] = aapp;
#line 1065 "zgesvj.f"
			if (ir1 == 0 && aapp == 0.) {
/* Computing MIN */
#line 1065 "zgesvj.f"
			    i__4 = igl + kbl - 1;
#line 1065 "zgesvj.f"
			    notrot = notrot + min(i__4,*n) - p;
#line 1065 "zgesvj.f"
			}
#line 1067 "zgesvj.f"
		    }

#line 1069 "zgesvj.f"
/* L2001: */
#line 1069 "zgesvj.f"
		}
/*     end of the p-loop */
/*     end of doing the block ( ibr, ibr ) */
#line 1072 "zgesvj.f"
/* L1002: */
#line 1072 "zgesvj.f"
	    }
/*     end of ir1-loop */

/* ... go to the off diagonal blocks */

#line 1077 "zgesvj.f"
	    igl = (ibr - 1) * kbl + 1;

#line 1079 "zgesvj.f"
	    i__2 = nbl;
#line 1079 "zgesvj.f"
	    for (jbc = ibr + 1; jbc <= i__2; ++jbc) {

#line 1081 "zgesvj.f"
		jgl = (jbc - 1) * kbl + 1;

/*        doing the block at ( ibr, jbc ) */

#line 1085 "zgesvj.f"
		ijblsk = 0;
/* Computing MIN */
#line 1086 "zgesvj.f"
		i__4 = igl + kbl - 1;
#line 1086 "zgesvj.f"
		i__3 = min(i__4,*n);
#line 1086 "zgesvj.f"
		for (p = igl; p <= i__3; ++p) {

#line 1088 "zgesvj.f"
		    aapp = sva[p];
#line 1089 "zgesvj.f"
		    if (aapp > 0.) {

#line 1091 "zgesvj.f"
			pskipped = 0;

/* Computing MIN */
#line 1093 "zgesvj.f"
			i__5 = jgl + kbl - 1;
#line 1093 "zgesvj.f"
			i__4 = min(i__5,*n);
#line 1093 "zgesvj.f"
			for (q = jgl; q <= i__4; ++q) {

#line 1095 "zgesvj.f"
			    aaqq = sva[q];
#line 1096 "zgesvj.f"
			    if (aaqq > 0.) {
#line 1097 "zgesvj.f"
				aapp0 = aapp;

/*     .. M x 2 Jacobi SVD .. */

/*        Safe Gram matrix computation */

#line 1103 "zgesvj.f"
				if (aaqq >= 1.) {
#line 1104 "zgesvj.f"
				    if (aapp >= aaqq) {
#line 1105 "zgesvj.f"
					rotok = small * aapp <= aaqq;
#line 1106 "zgesvj.f"
				    } else {
#line 1107 "zgesvj.f"
					rotok = small * aaqq <= aapp;
#line 1108 "zgesvj.f"
				    }
#line 1109 "zgesvj.f"
				    if (aapp < big / aaqq) {
#line 1110 "zgesvj.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1110 "zgesvj.f"
					z__2.r = z__3.r / aaqq, z__2.i = 
						z__3.i / aaqq;
#line 1110 "zgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1110 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1112 "zgesvj.f"
				    } else {
#line 1113 "zgesvj.f"
					zcopy_(m, &a[p * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 1115 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &aapp, &
						c_b42, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1118 "zgesvj.f"
					zdotc_(&z__2, m, &cwork[*n + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1118 "zgesvj.f"
					z__1.r = z__2.r / aaqq, z__1.i = 
						z__2.i / aaqq;
#line 1118 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1120 "zgesvj.f"
				    }
#line 1121 "zgesvj.f"
				} else {
#line 1122 "zgesvj.f"
				    if (aapp >= aaqq) {
#line 1123 "zgesvj.f"
					rotok = aapp <= aaqq / small;
#line 1124 "zgesvj.f"
				    } else {
#line 1125 "zgesvj.f"
					rotok = aaqq <= aapp / small;
#line 1126 "zgesvj.f"
				    }
#line 1127 "zgesvj.f"
				    if (aapp > small / aaqq) {
#line 1128 "zgesvj.f"
					zdotc_(&z__3, m, &a[p * a_dim1 + 1], &
						c__1, &a[q * a_dim1 + 1], &
						c__1);
#line 1128 "zgesvj.f"
					d__1 = max(aaqq,aapp);
#line 1128 "zgesvj.f"
					z__2.r = z__3.r / d__1, z__2.i = 
						z__3.i / d__1;
#line 1128 "zgesvj.f"
					d__2 = min(aaqq,aapp);
#line 1128 "zgesvj.f"
					z__1.r = z__2.r / d__2, z__1.i = 
						z__2.i / d__2;
#line 1128 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1131 "zgesvj.f"
				    } else {
#line 1132 "zgesvj.f"
					zcopy_(m, &a[q * a_dim1 + 1], &c__1, &
						cwork[*n + 1], &c__1);
#line 1134 "zgesvj.f"
					zlascl_("G", &c__0, &c__0, &aaqq, &
						c_b42, m, &c__1, &cwork[*n + 
						1], lda, &ierr, (ftnlen)1);
#line 1137 "zgesvj.f"
					zdotc_(&z__2, m, &a[p * a_dim1 + 1], &
						c__1, &cwork[*n + 1], &c__1);
#line 1137 "zgesvj.f"
					z__1.r = z__2.r / aapp, z__1.i = 
						z__2.i / aapp;
#line 1137 "zgesvj.f"
					aapq.r = z__1.r, aapq.i = z__1.i;
#line 1139 "zgesvj.f"
				    }
#line 1140 "zgesvj.f"
				}

/*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
#line 1144 "zgesvj.f"
				aapq1 = -z_abs(&aapq);
/* Computing MAX */
#line 1145 "zgesvj.f"
				d__1 = mxaapq, d__2 = -aapq1;
#line 1145 "zgesvj.f"
				mxaapq = max(d__1,d__2);

/*        TO rotate or NOT to rotate, THAT is the question ... */

#line 1149 "zgesvj.f"
				if (abs(aapq1) > tol) {
#line 1150 "zgesvj.f"
				    d__1 = z_abs(&aapq);
#line 1150 "zgesvj.f"
				    z__1.r = aapq.r / d__1, z__1.i = aapq.i / 
					    d__1;
#line 1150 "zgesvj.f"
				    ompq.r = z__1.r, ompq.i = z__1.i;
#line 1151 "zgesvj.f"
				    notrot = 0;
/* [RTD]      ROTATED  = ROTATED + 1 */
#line 1153 "zgesvj.f"
				    pskipped = 0;
#line 1154 "zgesvj.f"
				    ++iswrot;

#line 1156 "zgesvj.f"
				    if (rotok) {

#line 1158 "zgesvj.f"
					aqoap = aaqq / aapp;
#line 1159 "zgesvj.f"
					apoaq = aapp / aaqq;
#line 1160 "zgesvj.f"
					theta = (d__1 = aqoap - apoaq, abs(
						d__1)) * -.5 / aapq1;
#line 1161 "zgesvj.f"
					if (aaqq > aapp0) {
#line 1161 "zgesvj.f"
					    theta = -theta;
#line 1161 "zgesvj.f"
					}

#line 1163 "zgesvj.f"
					if (abs(theta) > bigtheta) {
#line 1164 "zgesvj.f"
					    t = .5 / theta;
#line 1165 "zgesvj.f"
					    cs = 1.;
#line 1166 "zgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 1166 "zgesvj.f"
					    z__1.r = t * z__2.r, z__1.i = t * 
						    z__2.i;
#line 1166 "zgesvj.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 1168 "zgesvj.f"
					    if (rsvec) {
#line 1169 "zgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 1169 "zgesvj.f"
			  z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1169 "zgesvj.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 1171 "zgesvj.f"
					    }
/* Computing MAX */
#line 1172 "zgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 1172 "zgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1174 "zgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 1174 "zgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));
/* Computing MAX */
#line 1176 "zgesvj.f"
					    d__1 = mxsinj, d__2 = abs(t);
#line 1176 "zgesvj.f"
					    mxsinj = max(d__1,d__2);
#line 1177 "zgesvj.f"
					} else {

/*                 .. choose correct signum for THETA and rotate */

#line 1181 "zgesvj.f"
					    thsign = -d_sign(&c_b42, &aapq1);
#line 1182 "zgesvj.f"
					    if (aaqq > aapp0) {
#line 1182 "zgesvj.f"
			  thsign = -thsign;
#line 1182 "zgesvj.f"
					    }
#line 1183 "zgesvj.f"
					    t = 1. / (theta + thsign * sqrt(
						    theta * theta + 1.));
#line 1185 "zgesvj.f"
					    cs = sqrt(1. / (t * t + 1.));
#line 1186 "zgesvj.f"
					    sn = t * cs;
/* Computing MAX */
#line 1187 "zgesvj.f"
					    d__1 = mxsinj, d__2 = abs(sn);
#line 1187 "zgesvj.f"
					    mxsinj = max(d__1,d__2);
/* Computing MAX */
#line 1188 "zgesvj.f"
					    d__1 = 0., d__2 = t * apoaq * 
						    aapq1 + 1.;
#line 1188 "zgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
/* Computing MAX */
#line 1190 "zgesvj.f"
					    d__1 = 0., d__2 = 1. - t * aqoap *
						     aapq1;
#line 1190 "zgesvj.f"
					    aapp *= sqrt((max(d__1,d__2)));

#line 1193 "zgesvj.f"
					    d_cnjg(&z__2, &ompq);
#line 1193 "zgesvj.f"
					    z__1.r = sn * z__2.r, z__1.i = sn 
						    * z__2.i;
#line 1193 "zgesvj.f"
					    zrot_(m, &a[p * a_dim1 + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1, &cs, &z__1);
#line 1195 "zgesvj.f"
					    if (rsvec) {
#line 1196 "zgesvj.f"
			  d_cnjg(&z__2, &ompq);
#line 1196 "zgesvj.f"
			  z__1.r = sn * z__2.r, z__1.i = sn * z__2.i;
#line 1196 "zgesvj.f"
			  zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * 
				  v_dim1 + 1], &c__1, &cs, &z__1);
#line 1198 "zgesvj.f"
					    }
#line 1199 "zgesvj.f"
					}
#line 1200 "zgesvj.f"
					i__5 = p;
#line 1200 "zgesvj.f"
					i__6 = q;
#line 1200 "zgesvj.f"
					z__2.r = -cwork[i__6].r, z__2.i = 
						-cwork[i__6].i;
#line 1200 "zgesvj.f"
					z__1.r = z__2.r * ompq.r - z__2.i * 
						ompq.i, z__1.i = z__2.r * 
						ompq.i + z__2.i * ompq.r;
#line 1200 "zgesvj.f"
					cwork[i__5].r = z__1.r, cwork[i__5].i 
						= z__1.i;

#line 1202 "zgesvj.f"
				    } else {
/*              .. have to use modified Gram-Schmidt like transformation */
#line 1204 "zgesvj.f"
					if (aapp > aaqq) {
#line 1205 "zgesvj.f"
					    zcopy_(m, &a[p * a_dim1 + 1], &
						    c__1, &cwork[*n + 1], &
						    c__1);
#line 1207 "zgesvj.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &cwork[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1210 "zgesvj.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1213 "zgesvj.f"
					    z__1.r = -aapq.r, z__1.i = 
						    -aapq.i;
#line 1213 "zgesvj.f"
					    zaxpy_(m, &z__1, &cwork[*n + 1], &
						    c__1, &a[q * a_dim1 + 1], 
						    &c__1);
#line 1215 "zgesvj.f"
					    zlascl_("G", &c__0, &c__0, &c_b42,
						     &aaqq, m, &c__1, &a[q * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1218 "zgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 1218 "zgesvj.f"
					    sva[q] = aaqq * sqrt((max(d__1,
						    d__2)));
#line 1220 "zgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1221 "zgesvj.f"
					} else {
#line 1222 "zgesvj.f"
					    zcopy_(m, &a[q * a_dim1 + 1], &
						    c__1, &cwork[*n + 1], &
						    c__1);
#line 1224 "zgesvj.f"
					    zlascl_("G", &c__0, &c__0, &aaqq, 
						    &c_b42, m, &c__1, &cwork[*
						    n + 1], lda, &ierr, (
						    ftnlen)1);
#line 1227 "zgesvj.f"
					    zlascl_("G", &c__0, &c__0, &aapp, 
						    &c_b42, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
#line 1230 "zgesvj.f"
					    d_cnjg(&z__2, &aapq);
#line 1230 "zgesvj.f"
					    z__1.r = -z__2.r, z__1.i = 
						    -z__2.i;
#line 1230 "zgesvj.f"
					    zaxpy_(m, &z__1, &cwork[*n + 1], &
						    c__1, &a[p * a_dim1 + 1], 
						    &c__1);
#line 1232 "zgesvj.f"
					    zlascl_("G", &c__0, &c__0, &c_b42,
						     &aapp, m, &c__1, &a[p * 
						    a_dim1 + 1], lda, &ierr, (
						    ftnlen)1);
/* Computing MAX */
#line 1235 "zgesvj.f"
					    d__1 = 0., d__2 = 1. - aapq1 * 
						    aapq1;
#line 1235 "zgesvj.f"
					    sva[p] = aapp * sqrt((max(d__1,
						    d__2)));
#line 1237 "zgesvj.f"
					    mxsinj = max(mxsinj,sfmin);
#line 1238 "zgesvj.f"
					}
#line 1239 "zgesvj.f"
				    }
/*           END IF ROTOK THEN ... ELSE */

/*           In the case of cancellation in updating SVA(q), SVA(p) */
/*           .. recompute SVA(q), SVA(p) */
/* Computing 2nd power */
#line 1244 "zgesvj.f"
				    d__1 = sva[q] / aaqq;
#line 1244 "zgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1246 "zgesvj.f"
					if (aaqq < rootbig && aaqq > 
						rootsfmin) {
#line 1248 "zgesvj.f"
					    sva[q] = dznrm2_(m, &a[q * a_dim1 
						    + 1], &c__1);
#line 1249 "zgesvj.f"
					} else {
#line 1250 "zgesvj.f"
					    t = 0.;
#line 1251 "zgesvj.f"
					    aaqq = 1.;
#line 1252 "zgesvj.f"
					    zlassq_(m, &a[q * a_dim1 + 1], &
						    c__1, &t, &aaqq);
#line 1254 "zgesvj.f"
					    sva[q] = t * sqrt(aaqq);
#line 1255 "zgesvj.f"
					}
#line 1256 "zgesvj.f"
				    }
/* Computing 2nd power */
#line 1257 "zgesvj.f"
				    d__1 = aapp / aapp0;
#line 1257 "zgesvj.f"
				    if (d__1 * d__1 <= rooteps) {
#line 1258 "zgesvj.f"
					if (aapp < rootbig && aapp > 
						rootsfmin) {
#line 1260 "zgesvj.f"
					    aapp = dznrm2_(m, &a[p * a_dim1 + 
						    1], &c__1);
#line 1261 "zgesvj.f"
					} else {
#line 1262 "zgesvj.f"
					    t = 0.;
#line 1263 "zgesvj.f"
					    aapp = 1.;
#line 1264 "zgesvj.f"
					    zlassq_(m, &a[p * a_dim1 + 1], &
						    c__1, &t, &aapp);
#line 1266 "zgesvj.f"
					    aapp = t * sqrt(aapp);
#line 1267 "zgesvj.f"
					}
#line 1268 "zgesvj.f"
					sva[p] = aapp;
#line 1269 "zgesvj.f"
				    }
/*              end of OK rotation */
#line 1271 "zgesvj.f"
				} else {
#line 1272 "zgesvj.f"
				    ++notrot;
/* [RTD]      SKIPPED  = SKIPPED  + 1 */
#line 1274 "zgesvj.f"
				    ++pskipped;
#line 1275 "zgesvj.f"
				    ++ijblsk;
#line 1276 "zgesvj.f"
				}
#line 1277 "zgesvj.f"
			    } else {
#line 1278 "zgesvj.f"
				++notrot;
#line 1279 "zgesvj.f"
				++pskipped;
#line 1280 "zgesvj.f"
				++ijblsk;
#line 1281 "zgesvj.f"
			    }

#line 1283 "zgesvj.f"
			    if (i__ <= swband && ijblsk >= blskip) {
#line 1285 "zgesvj.f"
				sva[p] = aapp;
#line 1286 "zgesvj.f"
				notrot = 0;
#line 1287 "zgesvj.f"
				goto L2011;
#line 1288 "zgesvj.f"
			    }
#line 1289 "zgesvj.f"
			    if (i__ <= swband && pskipped > rowskip) {
#line 1291 "zgesvj.f"
				aapp = -aapp;
#line 1292 "zgesvj.f"
				notrot = 0;
#line 1293 "zgesvj.f"
				goto L2203;
#line 1294 "zgesvj.f"
			    }

#line 1296 "zgesvj.f"
/* L2200: */
#line 1296 "zgesvj.f"
			}
/*        end of the q-loop */
#line 1298 "zgesvj.f"
L2203:

#line 1300 "zgesvj.f"
			sva[p] = aapp;

#line 1302 "zgesvj.f"
		    } else {

#line 1304 "zgesvj.f"
			if (aapp == 0.) {
/* Computing MIN */
#line 1304 "zgesvj.f"
			    i__4 = jgl + kbl - 1;
#line 1304 "zgesvj.f"
			    notrot = notrot + min(i__4,*n) - jgl + 1;
#line 1304 "zgesvj.f"
			}
#line 1306 "zgesvj.f"
			if (aapp < 0.) {
#line 1306 "zgesvj.f"
			    notrot = 0;
#line 1306 "zgesvj.f"
			}

#line 1308 "zgesvj.f"
		    }

#line 1310 "zgesvj.f"
/* L2100: */
#line 1310 "zgesvj.f"
		}
/*     end of the p-loop */
#line 1312 "zgesvj.f"
/* L2010: */
#line 1312 "zgesvj.f"
	    }
/*     end of the jbc-loop */
#line 1314 "zgesvj.f"
L2011:
/* 2011 bailed out of the jbc-loop */
/* Computing MIN */
#line 1316 "zgesvj.f"
	    i__3 = igl + kbl - 1;
#line 1316 "zgesvj.f"
	    i__2 = min(i__3,*n);
#line 1316 "zgesvj.f"
	    for (p = igl; p <= i__2; ++p) {
#line 1317 "zgesvj.f"
		sva[p] = (d__1 = sva[p], abs(d__1));
#line 1318 "zgesvj.f"
/* L2012: */
#line 1318 "zgesvj.f"
	    }
/* ** */
#line 1320 "zgesvj.f"
/* L2000: */
#line 1320 "zgesvj.f"
	}
/* 2000 :: end of the ibr-loop */

/*     .. update SVA(N) */
#line 1324 "zgesvj.f"
	if (sva[*n] < rootbig && sva[*n] > rootsfmin) {
#line 1326 "zgesvj.f"
	    sva[*n] = dznrm2_(m, &a[*n * a_dim1 + 1], &c__1);
#line 1327 "zgesvj.f"
	} else {
#line 1328 "zgesvj.f"
	    t = 0.;
#line 1329 "zgesvj.f"
	    aapp = 1.;
#line 1330 "zgesvj.f"
	    zlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
#line 1331 "zgesvj.f"
	    sva[*n] = t * sqrt(aapp);
#line 1332 "zgesvj.f"
	}

/*     Additional steering devices */

#line 1336 "zgesvj.f"
	if (i__ < swband && (mxaapq <= roottol || iswrot <= *n)) {
#line 1336 "zgesvj.f"
	    swband = i__;
#line 1336 "zgesvj.f"
	}

#line 1339 "zgesvj.f"
	if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * tol && (
		doublereal) (*n) * mxaapq * mxsinj < tol) {
#line 1341 "zgesvj.f"
	    goto L1994;
#line 1342 "zgesvj.f"
	}

#line 1344 "zgesvj.f"
	if (notrot >= emptsw) {
#line 1344 "zgesvj.f"
	    goto L1994;
#line 1344 "zgesvj.f"
	}

#line 1346 "zgesvj.f"
/* L1993: */
#line 1346 "zgesvj.f"
    }
/*     end i=1:NSWEEP loop */

/* #:( Reaching this point means that the procedure has not converged. */
#line 1350 "zgesvj.f"
    *info = 29;
#line 1351 "zgesvj.f"
    goto L1995;

#line 1353 "zgesvj.f"
L1994:
/* #:) Reaching this point means numerical convergence after the i-th */
/*     sweep. */

#line 1357 "zgesvj.f"
    *info = 0;
/* #:) INFO = 0 confirms successful iterations. */
#line 1359 "zgesvj.f"
L1995:

/*     Sort the singular values and find how many are above */
/*     the underflow threshold. */

#line 1364 "zgesvj.f"
    n2 = 0;
#line 1365 "zgesvj.f"
    n4 = 0;
#line 1366 "zgesvj.f"
    i__1 = *n - 1;
#line 1366 "zgesvj.f"
    for (p = 1; p <= i__1; ++p) {
#line 1367 "zgesvj.f"
	i__2 = *n - p + 1;
#line 1367 "zgesvj.f"
	q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
#line 1368 "zgesvj.f"
	if (p != q) {
#line 1369 "zgesvj.f"
	    temp1 = sva[p];
#line 1370 "zgesvj.f"
	    sva[p] = sva[q];
#line 1371 "zgesvj.f"
	    sva[q] = temp1;
#line 1372 "zgesvj.f"
	    zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
#line 1373 "zgesvj.f"
	    if (rsvec) {
#line 1373 "zgesvj.f"
		zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &
			c__1);
#line 1373 "zgesvj.f"
	    }
#line 1374 "zgesvj.f"
	}
#line 1375 "zgesvj.f"
	if (sva[p] != 0.) {
#line 1376 "zgesvj.f"
	    ++n4;
#line 1377 "zgesvj.f"
	    if (sva[p] * skl > sfmin) {
#line 1377 "zgesvj.f"
		++n2;
#line 1377 "zgesvj.f"
	    }
#line 1378 "zgesvj.f"
	}
#line 1379 "zgesvj.f"
/* L5991: */
#line 1379 "zgesvj.f"
    }
#line 1380 "zgesvj.f"
    if (sva[*n] != 0.) {
#line 1381 "zgesvj.f"
	++n4;
#line 1382 "zgesvj.f"
	if (sva[*n] * skl > sfmin) {
#line 1382 "zgesvj.f"
	    ++n2;
#line 1382 "zgesvj.f"
	}
#line 1383 "zgesvj.f"
    }

/*     Normalize the left singular vectors. */

#line 1387 "zgesvj.f"
    if (lsvec || uctol) {
#line 1388 "zgesvj.f"
	i__1 = n4;
#line 1388 "zgesvj.f"
	for (p = 1; p <= i__1; ++p) {
/*            CALL ZDSCAL( M, ONE / SVA( p ), A( 1, p ), 1 ) */
#line 1390 "zgesvj.f"
	    zlascl_("G", &c__0, &c__0, &sva[p], &c_b42, m, &c__1, &a[p * 
		    a_dim1 + 1], m, &ierr, (ftnlen)1);
#line 1391 "zgesvj.f"
/* L1998: */
#line 1391 "zgesvj.f"
	}
#line 1392 "zgesvj.f"
    }

/*     Scale the product of Jacobi rotations. */

#line 1396 "zgesvj.f"
    if (rsvec) {
#line 1397 "zgesvj.f"
	i__1 = *n;
#line 1397 "zgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1398 "zgesvj.f"
	    temp1 = 1. / dznrm2_(&mvl, &v[p * v_dim1 + 1], &c__1);
#line 1399 "zgesvj.f"
	    zdscal_(&mvl, &temp1, &v[p * v_dim1 + 1], &c__1);
#line 1400 "zgesvj.f"
/* L2399: */
#line 1400 "zgesvj.f"
	}
#line 1401 "zgesvj.f"
    }

/*     Undo scaling, if necessary (and possible). */
#line 1404 "zgesvj.f"
    if (skl > 1. && sva[1] < big / skl || skl < 1. && sva[max(n2,1)] > sfmin /
	     skl) {
#line 1407 "zgesvj.f"
	i__1 = *n;
#line 1407 "zgesvj.f"
	for (p = 1; p <= i__1; ++p) {
#line 1408 "zgesvj.f"
	    sva[p] = skl * sva[p];
#line 1409 "zgesvj.f"
/* L2400: */
#line 1409 "zgesvj.f"
	}
#line 1410 "zgesvj.f"
	skl = 1.;
#line 1411 "zgesvj.f"
    }

#line 1413 "zgesvj.f"
    rwork[1] = skl;
/*     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE */
/*     then some of the singular values may overflow or underflow and */
/*     the spectrum is given in this factored representation. */

#line 1418 "zgesvj.f"
    rwork[2] = (doublereal) n4;
/*     N4 is the number of computed nonzero singular values of A. */

#line 1421 "zgesvj.f"
    rwork[3] = (doublereal) n2;
/*     N2 is the number of singular values of A greater than SFMIN. */
/*     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers */
/*     that may carry some information. */

#line 1426 "zgesvj.f"
    rwork[4] = (doublereal) i__;
/*     i is the index of the last sweep before declaring convergence. */

#line 1429 "zgesvj.f"
    rwork[5] = mxaapq;
/*     MXAAPQ is the largest absolute value of scaled pivots in the */
/*     last sweep */

#line 1433 "zgesvj.f"
    rwork[6] = mxsinj;
/*     MXSINJ is the largest absolute value of the sines of Jacobi angles */
/*     in the last sweep */

#line 1437 "zgesvj.f"
    return 0;
/*     .. */
/*     .. END OF ZGESVJ */
/*     .. */
} /* zgesvj_ */

