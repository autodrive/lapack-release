#line 1 "checon_3.f"
/* checon_3.f -- translated by f2c (version 20100827).
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

#line 1 "checon_3.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHECON_3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHECON_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/checon_
3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/checon_
3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/checon_
3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHECON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, */
/*                            WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       COMPLEX            A( LDA, * ), E ( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > CHECON_3 estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian matrix A using the factorization */
/* > computed by CHETRF_RK or CHETRF_BK: */
/* > */
/* >    A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > This routine uses BLAS3 solver CHETRS_3. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix: */
/* >          = 'U':  Upper triangular, form is A = P*U*D*(U**H)*(P**T); */
/* >          = 'L':  Lower triangular, form is A = P*L*D*(L**H)*(P**T). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          Diagonal of the block diagonal matrix D and factors U or L */
/* >          as computed by CHETRF_RK and CHETRF_BK: */
/* >            a) ONLY diagonal elements of the Hermitian block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                should be provided on entry in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N) */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the Hermitian block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          NOTE: For 1-by-1 diagonal block D(k), where */
/* >          1 <= k <= N, the element E(k) is not referenced in both */
/* >          UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by CHETRF_RK or CHETRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* >          estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
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

/* > \date June 2017 */

/* > \ingroup complexHEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > \verbatim */
/* > */
/* >  June 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int checon_3__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *e, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublecomplex *work, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int chetrs_3__(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 217 "checon_3.f"
    /* Parameter adjustments */
#line 217 "checon_3.f"
    a_dim1 = *lda;
#line 217 "checon_3.f"
    a_offset = 1 + a_dim1;
#line 217 "checon_3.f"
    a -= a_offset;
#line 217 "checon_3.f"
    --e;
#line 217 "checon_3.f"
    --ipiv;
#line 217 "checon_3.f"
    --work;
#line 217 "checon_3.f"

#line 217 "checon_3.f"
    /* Function Body */
#line 217 "checon_3.f"
    *info = 0;
#line 218 "checon_3.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 219 "checon_3.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 220 "checon_3.f"
	*info = -1;
#line 221 "checon_3.f"
    } else if (*n < 0) {
#line 222 "checon_3.f"
	*info = -2;
#line 223 "checon_3.f"
    } else if (*lda < max(1,*n)) {
#line 224 "checon_3.f"
	*info = -4;
#line 225 "checon_3.f"
    } else if (*anorm < 0.) {
#line 226 "checon_3.f"
	*info = -7;
#line 227 "checon_3.f"
    }
#line 228 "checon_3.f"
    if (*info != 0) {
#line 229 "checon_3.f"
	i__1 = -(*info);
#line 229 "checon_3.f"
	xerbla_("CHECON_3", &i__1, (ftnlen)8);
#line 230 "checon_3.f"
	return 0;
#line 231 "checon_3.f"
    }

/*     Quick return if possible */

#line 235 "checon_3.f"
    *rcond = 0.;
#line 236 "checon_3.f"
    if (*n == 0) {
#line 237 "checon_3.f"
	*rcond = 1.;
#line 238 "checon_3.f"
	return 0;
#line 239 "checon_3.f"
    } else if (*anorm <= 0.) {
#line 240 "checon_3.f"
	return 0;
#line 241 "checon_3.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 245 "checon_3.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 249 "checon_3.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 250 "checon_3.f"
	    i__1 = i__ + i__ * a_dim1;
#line 250 "checon_3.f"
	    if (ipiv[i__] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 250 "checon_3.f"
		return 0;
#line 250 "checon_3.f"
	    }
#line 252 "checon_3.f"
	}
#line 253 "checon_3.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 257 "checon_3.f"
	i__1 = *n;
#line 257 "checon_3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "checon_3.f"
	    i__2 = i__ + i__ * a_dim1;
#line 258 "checon_3.f"
	    if (ipiv[i__] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 258 "checon_3.f"
		return 0;
#line 258 "checon_3.f"
	    }
#line 260 "checon_3.f"
	}
#line 261 "checon_3.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 265 "checon_3.f"
    kase = 0;
#line 266 "checon_3.f"
L30:
#line 267 "checon_3.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 268 "checon_3.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**H) or inv(U*D*U**H). */

#line 272 "checon_3.f"
	chetrs_3__(uplo, n, &c__1, &a[a_offset], lda, &e[1], &ipiv[1], &work[
		1], n, info, (ftnlen)1);
#line 273 "checon_3.f"
	goto L30;
#line 274 "checon_3.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 278 "checon_3.f"
    if (ainvnm != 0.) {
#line 278 "checon_3.f"
	*rcond = 1. / ainvnm / *anorm;
#line 278 "checon_3.f"
    }

#line 281 "checon_3.f"
    return 0;

/*     End of CHECON_3 */

} /* checon_3__ */

