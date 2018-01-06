#line 1 "zsytrs_aa.f"
/* zsytrs_aa.f -- translated by f2c (version 20100827).
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

#line 1 "zsytrs_aa.f"
/* Table of constant values */

static doublecomplex c_b9 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTRS_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRS_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/*                             WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRS_AA solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A using the factorization A = U*T*U**T or */
/* > A = L*T*L**T computed by ZSYTRF_AA. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*T*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*T*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          Details of factors computed by ZSYTRF_AA. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges as computed by ZSYTRF_AA. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER, LWORK >= MAX(1,3*N-2). */
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

/* > \date November 2017 */

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsytrs_aa__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer k, kp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zgtsv_(integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, integer *), ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), zlacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 169 "zsytrs_aa.f"
    /* Parameter adjustments */
#line 169 "zsytrs_aa.f"
    a_dim1 = *lda;
#line 169 "zsytrs_aa.f"
    a_offset = 1 + a_dim1;
#line 169 "zsytrs_aa.f"
    a -= a_offset;
#line 169 "zsytrs_aa.f"
    --ipiv;
#line 169 "zsytrs_aa.f"
    b_dim1 = *ldb;
#line 169 "zsytrs_aa.f"
    b_offset = 1 + b_dim1;
#line 169 "zsytrs_aa.f"
    b -= b_offset;
#line 169 "zsytrs_aa.f"
    --work;
#line 169 "zsytrs_aa.f"

#line 169 "zsytrs_aa.f"
    /* Function Body */
#line 169 "zsytrs_aa.f"
    *info = 0;
#line 170 "zsytrs_aa.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 171 "zsytrs_aa.f"
    lquery = *lwork == -1;
#line 172 "zsytrs_aa.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 173 "zsytrs_aa.f"
	*info = -1;
#line 174 "zsytrs_aa.f"
    } else if (*n < 0) {
#line 175 "zsytrs_aa.f"
	*info = -2;
#line 176 "zsytrs_aa.f"
    } else if (*nrhs < 0) {
#line 177 "zsytrs_aa.f"
	*info = -3;
#line 178 "zsytrs_aa.f"
    } else if (*lda < max(1,*n)) {
#line 179 "zsytrs_aa.f"
	*info = -5;
#line 180 "zsytrs_aa.f"
    } else if (*ldb < max(1,*n)) {
#line 181 "zsytrs_aa.f"
	*info = -8;
#line 182 "zsytrs_aa.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 182 "zsytrs_aa.f"
	i__1 = 1, i__2 = *n * 3 - 2;
#line 182 "zsytrs_aa.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 183 "zsytrs_aa.f"
	    *info = -10;
#line 184 "zsytrs_aa.f"
	}
#line 184 "zsytrs_aa.f"
    }
#line 185 "zsytrs_aa.f"
    if (*info != 0) {
#line 186 "zsytrs_aa.f"
	i__1 = -(*info);
#line 186 "zsytrs_aa.f"
	xerbla_("ZSYTRS_AA", &i__1, (ftnlen)9);
#line 187 "zsytrs_aa.f"
	return 0;
#line 188 "zsytrs_aa.f"
    } else if (lquery) {
#line 189 "zsytrs_aa.f"
	lwkopt = *n * 3 - 2;
#line 190 "zsytrs_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 191 "zsytrs_aa.f"
	return 0;
#line 192 "zsytrs_aa.f"
    }

/*     Quick return if possible */

#line 196 "zsytrs_aa.f"
    if (*n == 0 || *nrhs == 0) {
#line 196 "zsytrs_aa.f"
	return 0;
#line 196 "zsytrs_aa.f"
    }

#line 199 "zsytrs_aa.f"
    if (upper) {

/*        Solve A*X = B, where A = U*T*U**T. */

/*        Pivot, P**T * B */

#line 205 "zsytrs_aa.f"
	i__1 = *n;
#line 205 "zsytrs_aa.f"
	for (k = 1; k <= i__1; ++k) {
#line 206 "zsytrs_aa.f"
	    kp = ipiv[k];
#line 207 "zsytrs_aa.f"
	    if (kp != k) {
#line 207 "zsytrs_aa.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 207 "zsytrs_aa.f"
	    }
#line 209 "zsytrs_aa.f"
	}

/*        Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 213 "zsytrs_aa.f"
	i__1 = *n - 1;
#line 213 "zsytrs_aa.f"
	ztrsm_("L", "U", "T", "U", &i__1, nrhs, &c_b9, &a[(a_dim1 << 1) + 1], 
		lda, &b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Compute T \ B -> B   [ T \ (U \P**T * B) ] */

#line 218 "zsytrs_aa.f"
	i__1 = *lda + 1;
#line 218 "zsytrs_aa.f"
	zlacpy_("F", &c__1, n, &a[a_dim1 + 1], &i__1, &work[*n], &c__1, (
		ftnlen)1);
#line 219 "zsytrs_aa.f"
	if (*n > 1) {
#line 220 "zsytrs_aa.f"
	    i__1 = *n - 1;
#line 220 "zsytrs_aa.f"
	    i__2 = *lda + 1;
#line 220 "zsytrs_aa.f"
	    zlacpy_("F", &c__1, &i__1, &a[(a_dim1 << 1) + 1], &i__2, &work[1],
		     &c__1, (ftnlen)1);
#line 221 "zsytrs_aa.f"
	    i__1 = *n - 1;
#line 221 "zsytrs_aa.f"
	    i__2 = *lda + 1;
#line 221 "zsytrs_aa.f"
	    zlacpy_("F", &c__1, &i__1, &a[(a_dim1 << 1) + 1], &i__2, &work[*n 
		    * 2], &c__1, (ftnlen)1);
#line 222 "zsytrs_aa.f"
	}
#line 223 "zsytrs_aa.f"
	zgtsv_(n, nrhs, &work[1], &work[*n], &work[*n * 2], &b[b_offset], ldb,
		 info);

/*        Compute (U**T \ B) -> B   [ U**T \ (T \ (U \P**T * B) ) ] */

#line 228 "zsytrs_aa.f"
	i__1 = *n - 1;
#line 228 "zsytrs_aa.f"
	ztrsm_("L", "U", "N", "U", &i__1, nrhs, &c_b9, &a[(a_dim1 << 1) + 1], 
		lda, &b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Pivot, P * B  [ P * (U**T \ (T \ (U \P**T * B) )) ] */

#line 233 "zsytrs_aa.f"
	for (k = *n; k >= 1; --k) {
#line 234 "zsytrs_aa.f"
	    kp = ipiv[k];
#line 235 "zsytrs_aa.f"
	    if (kp != k) {
#line 235 "zsytrs_aa.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 235 "zsytrs_aa.f"
	    }
#line 237 "zsytrs_aa.f"
	}

#line 239 "zsytrs_aa.f"
    } else {

/*        Solve A*X = B, where A = L*T*L**T. */

/*        Pivot, P**T * B */

#line 245 "zsytrs_aa.f"
	i__1 = *n;
#line 245 "zsytrs_aa.f"
	for (k = 1; k <= i__1; ++k) {
#line 246 "zsytrs_aa.f"
	    kp = ipiv[k];
#line 247 "zsytrs_aa.f"
	    if (kp != k) {
#line 247 "zsytrs_aa.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 247 "zsytrs_aa.f"
	    }
#line 249 "zsytrs_aa.f"
	}

/*        Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 253 "zsytrs_aa.f"
	i__1 = *n - 1;
#line 253 "zsytrs_aa.f"
	ztrsm_("L", "L", "N", "U", &i__1, nrhs, &c_b9, &a[a_dim1 + 2], lda, &
		b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		1);

/*        Compute T \ B -> B   [ T \ (L \P**T * B) ] */

#line 258 "zsytrs_aa.f"
	i__1 = *lda + 1;
#line 258 "zsytrs_aa.f"
	zlacpy_("F", &c__1, n, &a[a_dim1 + 1], &i__1, &work[*n], &c__1, (
		ftnlen)1);
#line 259 "zsytrs_aa.f"
	if (*n > 1) {
#line 260 "zsytrs_aa.f"
	    i__1 = *n - 1;
#line 260 "zsytrs_aa.f"
	    i__2 = *lda + 1;
#line 260 "zsytrs_aa.f"
	    zlacpy_("F", &c__1, &i__1, &a[a_dim1 + 2], &i__2, &work[1], &c__1,
		     (ftnlen)1);
#line 261 "zsytrs_aa.f"
	    i__1 = *n - 1;
#line 261 "zsytrs_aa.f"
	    i__2 = *lda + 1;
#line 261 "zsytrs_aa.f"
	    zlacpy_("F", &c__1, &i__1, &a[a_dim1 + 2], &i__2, &work[*n * 2], &
		    c__1, (ftnlen)1);
#line 262 "zsytrs_aa.f"
	}
#line 263 "zsytrs_aa.f"
	zgtsv_(n, nrhs, &work[1], &work[*n], &work[*n * 2], &b[b_offset], ldb,
		 info);

/*        Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ] */

#line 268 "zsytrs_aa.f"
	i__1 = *n - 1;
#line 268 "zsytrs_aa.f"
	ztrsm_("L", "L", "T", "U", &i__1, nrhs, &c_b9, &a[a_dim1 + 2], lda, &
		b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		1);

/*        Pivot, P * B  [ P * (L**T \ (T \ (L \P**T * B) )) ] */

#line 273 "zsytrs_aa.f"
	for (k = *n; k >= 1; --k) {
#line 274 "zsytrs_aa.f"
	    kp = ipiv[k];
#line 275 "zsytrs_aa.f"
	    if (kp != k) {
#line 275 "zsytrs_aa.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 275 "zsytrs_aa.f"
	    }
#line 277 "zsytrs_aa.f"
	}

#line 279 "zsytrs_aa.f"
    }

#line 281 "zsytrs_aa.f"
    return 0;

/*     End of ZSYTRS_AA */

} /* zsytrs_aa__ */

