#line 1 "dsytri_3.f"
/* dsytri_3.f -- translated by f2c (version 20100827).
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

#line 1 "dsytri_3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b DSYTRI_3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRI_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri_
3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri_
3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri_
3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, */
/*                            INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > DSYTRI_3 computes the inverse of a real symmetric indefinite */
/* > matrix A using the factorization computed by DSYTRF_RK or DSYTRF_BK: */
/* > */
/* >     A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > DSYTRI_3 sets the leading dimension of the workspace  before calling */
/* > DSYTRI_3X that actually computes the inverse.  This is the blocked */
/* > version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix. */
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
/* >          On entry, diagonal of the block diagonal matrix D and */
/* >          factors U or L as computed by DSYTRF_RK and DSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                should be provided on entry in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, if INFO = 0, the symmetric inverse of the original */
/* >          matrix. */
/* >             If UPLO = 'U': the upper triangular part of the inverse */
/* >             is formed and the part of A below the diagonal is not */
/* >             referenced; */
/* >             If UPLO = 'L': the lower triangular part of the inverse */
/* >             is formed and the part of A above the diagonal is not */
/* >             referenced. */
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
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not refernced; */
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
/* >          as determined by DSYTRF_RK or DSYTRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N+NB+1)*(NB+3). */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK. LWORK >= (N+NB+1)*(NB+3). */
/* > */
/* >          If LDWORK = -1, then a workspace query is assumed; */
/* >          the routine only calculates the optimal size of the optimal */
/* >          size of the WORK array, returns this value as the first */
/* >          entry of the WORK array, and no error message related to */
/* >          LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its */
/* >               inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int dsytri_3__(char *uplo, integer *n, doublereal *a, 
	integer *lda, doublereal *e, integer *ipiv, doublereal *work, integer 
	*lwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int dsytri_3x__(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
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

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 208 "dsytri_3.f"
    /* Parameter adjustments */
#line 208 "dsytri_3.f"
    a_dim1 = *lda;
#line 208 "dsytri_3.f"
    a_offset = 1 + a_dim1;
#line 208 "dsytri_3.f"
    a -= a_offset;
#line 208 "dsytri_3.f"
    --e;
#line 208 "dsytri_3.f"
    --ipiv;
#line 208 "dsytri_3.f"
    --work;
#line 208 "dsytri_3.f"

#line 208 "dsytri_3.f"
    /* Function Body */
#line 208 "dsytri_3.f"
    *info = 0;
#line 209 "dsytri_3.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 210 "dsytri_3.f"
    lquery = *lwork == -1;

/*     Determine the block size */

/* Computing MAX */
#line 214 "dsytri_3.f"
    i__1 = 1, i__2 = ilaenv_(&c__1, "DSYTRI_3", uplo, n, &c_n1, &c_n1, &c_n1, 
	    (ftnlen)8, (ftnlen)1);
#line 214 "dsytri_3.f"
    nb = max(i__1,i__2);
#line 215 "dsytri_3.f"
    lwkopt = (*n + nb + 1) * (nb + 3);

#line 217 "dsytri_3.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 218 "dsytri_3.f"
	*info = -1;
#line 219 "dsytri_3.f"
    } else if (*n < 0) {
#line 220 "dsytri_3.f"
	*info = -2;
#line 221 "dsytri_3.f"
    } else if (*lda < max(1,*n)) {
#line 222 "dsytri_3.f"
	*info = -4;
#line 223 "dsytri_3.f"
    } else if (*lwork < lwkopt && ! lquery) {
#line 224 "dsytri_3.f"
	*info = -8;
#line 225 "dsytri_3.f"
    }

#line 227 "dsytri_3.f"
    if (*info != 0) {
#line 228 "dsytri_3.f"
	i__1 = -(*info);
#line 228 "dsytri_3.f"
	xerbla_("DSYTRI_3", &i__1, (ftnlen)8);
#line 229 "dsytri_3.f"
	return 0;
#line 230 "dsytri_3.f"
    } else if (lquery) {
#line 231 "dsytri_3.f"
	work[1] = (doublereal) lwkopt;
#line 232 "dsytri_3.f"
	return 0;
#line 233 "dsytri_3.f"
    }

/*     Quick return if possible */

#line 237 "dsytri_3.f"
    if (*n == 0) {
#line 237 "dsytri_3.f"
	return 0;
#line 237 "dsytri_3.f"
    }

#line 240 "dsytri_3.f"
    dsytri_3x__(uplo, n, &a[a_offset], lda, &e[1], &ipiv[1], &work[1], &nb, 
	    info, (ftnlen)1);

#line 242 "dsytri_3.f"
    work[1] = (doublereal) lwkopt;

#line 244 "dsytri_3.f"
    return 0;

/*     End of DSYTRI_3 */

} /* dsytri_3__ */

