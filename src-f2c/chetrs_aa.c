#line 1 "chetrs_aa.f"
/* chetrs_aa.f -- translated by f2c (version 20100827).
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

#line 1 "chetrs_aa.f"
/* Table of constant values */

static doublecomplex c_b9 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHETRS_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRS_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/*                             WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRS_AA solves a system of linear equations A*X = B with a complex */
/* > hermitian matrix A using the factorization A = U*T*U**H or */
/* > A = L*T*L**H computed by CHETRF_AA. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*T*U**H; */
/* >          = 'L':  Lower triangular, form is A = L*T*L**H. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          Details of factors computed by CHETRF_AA. */
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
/* >          Details of the interchanges as computed by CHETRF_AA. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
/* Subroutine */ int chetrs_aa__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer k, kp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cgtsv_(integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, integer *), ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , clacpy_(char *, integer *, integer *, doublecomplex *, integer *
	    , doublecomplex *, integer *, ftnlen), xerbla_(char *, integer *, 
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

#line 169 "chetrs_aa.f"
    /* Parameter adjustments */
#line 169 "chetrs_aa.f"
    a_dim1 = *lda;
#line 169 "chetrs_aa.f"
    a_offset = 1 + a_dim1;
#line 169 "chetrs_aa.f"
    a -= a_offset;
#line 169 "chetrs_aa.f"
    --ipiv;
#line 169 "chetrs_aa.f"
    b_dim1 = *ldb;
#line 169 "chetrs_aa.f"
    b_offset = 1 + b_dim1;
#line 169 "chetrs_aa.f"
    b -= b_offset;
#line 169 "chetrs_aa.f"
    --work;
#line 169 "chetrs_aa.f"

#line 169 "chetrs_aa.f"
    /* Function Body */
#line 169 "chetrs_aa.f"
    *info = 0;
#line 170 "chetrs_aa.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 171 "chetrs_aa.f"
    lquery = *lwork == -1;
#line 172 "chetrs_aa.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 173 "chetrs_aa.f"
	*info = -1;
#line 174 "chetrs_aa.f"
    } else if (*n < 0) {
#line 175 "chetrs_aa.f"
	*info = -2;
#line 176 "chetrs_aa.f"
    } else if (*nrhs < 0) {
#line 177 "chetrs_aa.f"
	*info = -3;
#line 178 "chetrs_aa.f"
    } else if (*lda < max(1,*n)) {
#line 179 "chetrs_aa.f"
	*info = -5;
#line 180 "chetrs_aa.f"
    } else if (*ldb < max(1,*n)) {
#line 181 "chetrs_aa.f"
	*info = -8;
#line 182 "chetrs_aa.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 182 "chetrs_aa.f"
	i__1 = 1, i__2 = *n * 3 - 2;
#line 182 "chetrs_aa.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 183 "chetrs_aa.f"
	    *info = -10;
#line 184 "chetrs_aa.f"
	}
#line 184 "chetrs_aa.f"
    }
#line 185 "chetrs_aa.f"
    if (*info != 0) {
#line 186 "chetrs_aa.f"
	i__1 = -(*info);
#line 186 "chetrs_aa.f"
	xerbla_("CHETRS_AA", &i__1, (ftnlen)9);
#line 187 "chetrs_aa.f"
	return 0;
#line 188 "chetrs_aa.f"
    } else if (lquery) {
#line 189 "chetrs_aa.f"
	lwkopt = *n * 3 - 2;
#line 190 "chetrs_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 191 "chetrs_aa.f"
	return 0;
#line 192 "chetrs_aa.f"
    }

/*     Quick return if possible */

#line 196 "chetrs_aa.f"
    if (*n == 0 || *nrhs == 0) {
#line 196 "chetrs_aa.f"
	return 0;
#line 196 "chetrs_aa.f"
    }

#line 199 "chetrs_aa.f"
    if (upper) {

/*        Solve A*X = B, where A = U*T*U**T. */

/*        P**T * B */

#line 205 "chetrs_aa.f"
	k = 1;
#line 206 "chetrs_aa.f"
	while(k <= *n) {
#line 207 "chetrs_aa.f"
	    kp = ipiv[k];
#line 208 "chetrs_aa.f"
	    if (kp != k) {
#line 208 "chetrs_aa.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 208 "chetrs_aa.f"
	    }
#line 210 "chetrs_aa.f"
	    ++k;
#line 211 "chetrs_aa.f"
	}

/*        Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 215 "chetrs_aa.f"
	i__1 = *n - 1;
#line 215 "chetrs_aa.f"
	ctrsm_("L", "U", "C", "U", &i__1, nrhs, &c_b9, &a[(a_dim1 << 1) + 1], 
		lda, &b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Compute T \ B -> B   [ T \ (U \P**T * B) ] */

#line 220 "chetrs_aa.f"
	i__1 = *lda + 1;
#line 220 "chetrs_aa.f"
	clacpy_("F", &c__1, n, &a[a_dim1 + 1], &i__1, &work[*n], &c__1, (
		ftnlen)1);
#line 221 "chetrs_aa.f"
	if (*n > 1) {
#line 222 "chetrs_aa.f"
	    i__1 = *n - 1;
#line 222 "chetrs_aa.f"
	    i__2 = *lda + 1;
#line 222 "chetrs_aa.f"
	    clacpy_("F", &c__1, &i__1, &a[(a_dim1 << 1) + 1], &i__2, &work[*n 
		    * 2], &c__1, (ftnlen)1);
#line 223 "chetrs_aa.f"
	    i__1 = *n - 1;
#line 223 "chetrs_aa.f"
	    i__2 = *lda + 1;
#line 223 "chetrs_aa.f"
	    clacpy_("F", &c__1, &i__1, &a[(a_dim1 << 1) + 1], &i__2, &work[1],
		     &c__1, (ftnlen)1);
#line 224 "chetrs_aa.f"
	    i__1 = *n - 1;
#line 224 "chetrs_aa.f"
	    clacgv_(&i__1, &work[1], &c__1);
#line 225 "chetrs_aa.f"
	}
#line 226 "chetrs_aa.f"
	cgtsv_(n, nrhs, &work[1], &work[*n], &work[*n * 2], &b[b_offset], ldb,
		 info);

/*        Compute (U**T \ B) -> B   [ U**T \ (T \ (U \P**T * B) ) ] */

#line 231 "chetrs_aa.f"
	i__1 = *n - 1;
#line 231 "chetrs_aa.f"
	ctrsm_("L", "U", "N", "U", &i__1, nrhs, &c_b9, &a[(a_dim1 << 1) + 1], 
		lda, &b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Pivot, P * B  [ P * (U**T \ (T \ (U \P**T * B) )) ] */

#line 236 "chetrs_aa.f"
	k = *n;
#line 237 "chetrs_aa.f"
	while(k >= 1) {
#line 238 "chetrs_aa.f"
	    kp = ipiv[k];
#line 239 "chetrs_aa.f"
	    if (kp != k) {
#line 239 "chetrs_aa.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 239 "chetrs_aa.f"
	    }
#line 241 "chetrs_aa.f"
	    --k;
#line 242 "chetrs_aa.f"
	}

#line 244 "chetrs_aa.f"
    } else {

/*        Solve A*X = B, where A = L*T*L**T. */

/*        Pivot, P**T * B */

#line 250 "chetrs_aa.f"
	k = 1;
#line 251 "chetrs_aa.f"
	while(k <= *n) {
#line 252 "chetrs_aa.f"
	    kp = ipiv[k];
#line 253 "chetrs_aa.f"
	    if (kp != k) {
#line 253 "chetrs_aa.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 253 "chetrs_aa.f"
	    }
#line 255 "chetrs_aa.f"
	    ++k;
#line 256 "chetrs_aa.f"
	}

/*        Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 260 "chetrs_aa.f"
	i__1 = *n - 1;
#line 260 "chetrs_aa.f"
	ctrsm_("L", "L", "N", "U", &i__1, nrhs, &c_b9, &a[a_dim1 + 2], lda, &
		b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		1);

/*        Compute T \ B -> B   [ T \ (L \P**T * B) ] */

#line 265 "chetrs_aa.f"
	i__1 = *lda + 1;
#line 265 "chetrs_aa.f"
	clacpy_("F", &c__1, n, &a[a_dim1 + 1], &i__1, &work[*n], &c__1, (
		ftnlen)1);
#line 266 "chetrs_aa.f"
	if (*n > 1) {
#line 267 "chetrs_aa.f"
	    i__1 = *n - 1;
#line 267 "chetrs_aa.f"
	    i__2 = *lda + 1;
#line 267 "chetrs_aa.f"
	    clacpy_("F", &c__1, &i__1, &a[a_dim1 + 2], &i__2, &work[1], &c__1,
		     (ftnlen)1);
#line 268 "chetrs_aa.f"
	    i__1 = *n - 1;
#line 268 "chetrs_aa.f"
	    i__2 = *lda + 1;
#line 268 "chetrs_aa.f"
	    clacpy_("F", &c__1, &i__1, &a[a_dim1 + 2], &i__2, &work[*n * 2], &
		    c__1, (ftnlen)1);
#line 269 "chetrs_aa.f"
	    i__1 = *n - 1;
#line 269 "chetrs_aa.f"
	    clacgv_(&i__1, &work[*n * 2], &c__1);
#line 270 "chetrs_aa.f"
	}
#line 271 "chetrs_aa.f"
	cgtsv_(n, nrhs, &work[1], &work[*n], &work[*n * 2], &b[b_offset], ldb,
		 info);

/*        Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ] */

#line 276 "chetrs_aa.f"
	i__1 = *n - 1;
#line 276 "chetrs_aa.f"
	ctrsm_("L", "L", "C", "U", &i__1, nrhs, &c_b9, &a[a_dim1 + 2], lda, &
		b[b_dim1 + 2], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		1);

/*        Pivot, P * B  [ P * (L**T \ (T \ (L \P**T * B) )) ] */

#line 281 "chetrs_aa.f"
	k = *n;
#line 282 "chetrs_aa.f"
	while(k >= 1) {
#line 283 "chetrs_aa.f"
	    kp = ipiv[k];
#line 284 "chetrs_aa.f"
	    if (kp != k) {
#line 284 "chetrs_aa.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 284 "chetrs_aa.f"
	    }
#line 286 "chetrs_aa.f"
	    --k;
#line 287 "chetrs_aa.f"
	}

#line 289 "chetrs_aa.f"
    }

#line 291 "chetrs_aa.f"
    return 0;

/*     End of CHETRS_AA */

} /* chetrs_aa__ */

