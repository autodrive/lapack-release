#line 1 "dsytrd_2stage.f"
/* dsytrd_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsytrd_2stage.f"
/* Table of constant values */

static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;

/* > \brief \b DSYTRD_2STAGE */

/*  @generated from zhetrd_2stage.f, fortran z -> d, Sun Nov  6 19:34:06 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrd_
2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrd_
2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrd_
2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU, */
/*                                 HOUS2, LHOUS2, WORK, LWORK, INFO ) */

/*       IMPLICIT NONE */

/*      .. Scalar Arguments .. */
/*       CHARACTER          VECT, UPLO */
/*       INTEGER            N, LDA, LWORK, LHOUS2, INFO */
/*      .. */
/*      .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), */
/*                          HOUS2( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRD_2STAGE reduces a real symmetric matrix A to real symmetric */
/* > tridiagonal form T by a orthogonal similarity transformation: */
/* > Q1**T Q2**T* A * Q2 * Q1 = T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          = 'N':  No need for the Housholder representation, */
/* >                  in particular for the second stage (Band to */
/* >                  tridiagonal) and thus LHOUS2 is of size max(1, 4*N); */
/* >          = 'V':  the Householder representation is needed to */
/* >                  either generate Q1 Q2 or to apply Q1 Q2, */
/* >                  then LHOUS2 is to be queried and computed. */
/* >                  (NOT AVAILABLE IN THIS RELEASE). */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if UPLO = 'U', the band superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          internal band-diagonal matrix AB, and the elements above */
/* >          the KD superdiagonal, with the array TAU, represent the orthogonal */
/* >          matrix Q1 as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and band subdiagonal of A are over- */
/* >          written by the corresponding elements of the internal band-diagonal */
/* >          matrix AB, and the elements below the KD subdiagonal, with */
/* >          the array TAU, represent the orthogonal matrix Q1 as a product */
/* >          of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (N-KD) */
/* >          The scalar factors of the elementary reflectors of */
/* >          the first stage (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[out] HOUS2 */
/* > \verbatim */
/* >          HOUS2 is DOUBLE PRECISION array, dimension LHOUS2, that */
/* >          store the Householder representation of the stage2 */
/* >          band to tridiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] LHOUS2 */
/* > \verbatim */
/* >          LHOUS2 is INTEGER */
/* >          The dimension of the array HOUS2. LHOUS2 = MAX(1, dimension) */
/* >          If LWORK = -1, or LHOUS2=-1, */
/* >          then a query is assumed; the routine */
/* >          only calculates the optimal size of the HOUS2 array, returns */
/* >          this value as the first entry of the HOUS2 array, and no error */
/* >          message related to LHOUS2 is issued by XERBLA. */
/* >          LHOUS2 = MAX(1, dimension) where */
/* >          dimension = 4*N if VECT='N' */
/* >          not available now if VECT='H' */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK = MAX(1, dimension) */
/* >          If LWORK = -1, or LHOUS2=-1, */
/* >          then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* >          LWORK = MAX(1, dimension) where */
/* >          dimension   = max(stage1,stage2) + (KD+1)*N */
/* >                      = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                        + max(2*KD*KD, KD*NTHREADS) */
/* >                        + (KD+1)*N */
/* >          where KD is the blocking size of the reduction, */
/* >          FACTOPTNB is the blocking used by the QR or LQ */
/* >          algorithm, usually FACTOPTNB=128 is a good choice */
/* >          NTHREADS is the number of threads used when */
/* >          openMP compilation is enabled, otherwise =1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleSYcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Implemented by Azzam Haidar. */
/* > */
/* >  All details are available on technical report, SC11, SC13 papers. */
/* > */
/* >  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* >  Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* >  using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* >  of 2011 International Conference for High Performance Computing, */
/* >  Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* >  Article 8 , 11 pages. */
/* >  http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* >  A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* >  An improved parallel singular value algorithm and its implementation */
/* >  for multicore hardware, In Proceedings of 2013 International Conference */
/* >  for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* >  Denver, Colorado, USA, 2013. */
/* >  Article 90, 12 pages. */
/* >  http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* >  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* >  A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* >  calculations based on fine-grained memory aware tasks. */
/* >  International Journal of High Performance Computing Applications. */
/* >  Volume 28 Issue 2, Pages 196-209, May 2014. */
/* >  http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsytrd_2stage__(char *vect, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *d__, doublereal *e, 
	doublereal *tau, doublereal *hous2, integer *lhous2, doublereal *work,
	 integer *lwork, integer *info, ftnlen vect_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer ib, kd, ldab, lwrk;
    extern /* Subroutine */ int dsytrd_sb2st__(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer wpos;
    extern /* Subroutine */ int dsytrd_sy2sb__(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer abpos, lhmin, lwmin;
    static logical wantq, upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 263 "dsytrd_2stage.f"
    /* Parameter adjustments */
#line 263 "dsytrd_2stage.f"
    a_dim1 = *lda;
#line 263 "dsytrd_2stage.f"
    a_offset = 1 + a_dim1;
#line 263 "dsytrd_2stage.f"
    a -= a_offset;
#line 263 "dsytrd_2stage.f"
    --d__;
#line 263 "dsytrd_2stage.f"
    --e;
#line 263 "dsytrd_2stage.f"
    --tau;
#line 263 "dsytrd_2stage.f"
    --hous2;
#line 263 "dsytrd_2stage.f"
    --work;
#line 263 "dsytrd_2stage.f"

#line 263 "dsytrd_2stage.f"
    /* Function Body */
#line 263 "dsytrd_2stage.f"
    *info = 0;
#line 264 "dsytrd_2stage.f"
    wantq = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 265 "dsytrd_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 266 "dsytrd_2stage.f"
    lquery = *lwork == -1 || *lhous2 == -1;

/*     Determine the block size, the workspace size and the hous size. */

#line 270 "dsytrd_2stage.f"
    kd = ilaenv_(&c__17, "DSYTRD_2STAGE", vect, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 271 "dsytrd_2stage.f"
    ib = ilaenv_(&c__18, "DSYTRD_2STAGE", vect, n, &kd, &c_n1, &c_n1, (ftnlen)
	    13, (ftnlen)1);
#line 272 "dsytrd_2stage.f"
    lhmin = ilaenv_(&c__19, "DSYTRD_2STAGE", vect, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
#line 273 "dsytrd_2stage.f"
    lwmin = ilaenv_(&c__20, "DSYTRD_2STAGE", vect, n, &kd, &ib, &c_n1, (
	    ftnlen)13, (ftnlen)1);
/*      WRITE(*,*),'DSYTRD_2STAGE N KD UPLO LHMIN LWMIN ',N, KD, UPLO, */
/*     $            LHMIN, LWMIN */

#line 277 "dsytrd_2stage.f"
    if (! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 278 "dsytrd_2stage.f"
	*info = -1;
#line 279 "dsytrd_2stage.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 280 "dsytrd_2stage.f"
	*info = -2;
#line 281 "dsytrd_2stage.f"
    } else if (*n < 0) {
#line 282 "dsytrd_2stage.f"
	*info = -3;
#line 283 "dsytrd_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 284 "dsytrd_2stage.f"
	*info = -5;
#line 285 "dsytrd_2stage.f"
    } else if (*lhous2 < lhmin && ! lquery) {
#line 286 "dsytrd_2stage.f"
	*info = -10;
#line 287 "dsytrd_2stage.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 288 "dsytrd_2stage.f"
	*info = -12;
#line 289 "dsytrd_2stage.f"
    }

#line 291 "dsytrd_2stage.f"
    if (*info == 0) {
#line 292 "dsytrd_2stage.f"
	hous2[1] = (doublereal) lhmin;
#line 293 "dsytrd_2stage.f"
	work[1] = (doublereal) lwmin;
#line 294 "dsytrd_2stage.f"
    }

#line 296 "dsytrd_2stage.f"
    if (*info != 0) {
#line 297 "dsytrd_2stage.f"
	i__1 = -(*info);
#line 297 "dsytrd_2stage.f"
	xerbla_("DSYTRD_2STAGE", &i__1, (ftnlen)13);
#line 298 "dsytrd_2stage.f"
	return 0;
#line 299 "dsytrd_2stage.f"
    } else if (lquery) {
#line 300 "dsytrd_2stage.f"
	return 0;
#line 301 "dsytrd_2stage.f"
    }

/*     Quick return if possible */

#line 305 "dsytrd_2stage.f"
    if (*n == 0) {
#line 306 "dsytrd_2stage.f"
	work[1] = 1.;
#line 307 "dsytrd_2stage.f"
	return 0;
#line 308 "dsytrd_2stage.f"
    }

/*     Determine pointer position */

#line 312 "dsytrd_2stage.f"
    ldab = kd + 1;
#line 313 "dsytrd_2stage.f"
    lwrk = *lwork - ldab * *n;
#line 314 "dsytrd_2stage.f"
    abpos = 1;
#line 315 "dsytrd_2stage.f"
    wpos = abpos + ldab * *n;
#line 316 "dsytrd_2stage.f"
    dsytrd_sy2sb__(uplo, n, &kd, &a[a_offset], lda, &work[abpos], &ldab, &tau[
	    1], &work[wpos], &lwrk, info, (ftnlen)1);
#line 318 "dsytrd_2stage.f"
    if (*info != 0) {
#line 319 "dsytrd_2stage.f"
	i__1 = -(*info);
#line 319 "dsytrd_2stage.f"
	xerbla_("DSYTRD_SY2SB", &i__1, (ftnlen)12);
#line 320 "dsytrd_2stage.f"
	return 0;
#line 321 "dsytrd_2stage.f"
    }
#line 322 "dsytrd_2stage.f"
    dsytrd_sb2st__("Y", vect, uplo, n, &kd, &work[abpos], &ldab, &d__[1], &e[
	    1], &hous2[1], lhous2, &work[wpos], &lwrk, info, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1);
#line 325 "dsytrd_2stage.f"
    if (*info != 0) {
#line 326 "dsytrd_2stage.f"
	i__1 = -(*info);
#line 326 "dsytrd_2stage.f"
	xerbla_("DSYTRD_SB2ST", &i__1, (ftnlen)12);
#line 327 "dsytrd_2stage.f"
	return 0;
#line 328 "dsytrd_2stage.f"
    }


#line 331 "dsytrd_2stage.f"
    hous2[1] = (doublereal) lhmin;
#line 332 "dsytrd_2stage.f"
    work[1] = (doublereal) lwmin;
#line 333 "dsytrd_2stage.f"
    return 0;

/*     End of DSYTRD_2STAGE */

} /* dsytrd_2stage__ */

