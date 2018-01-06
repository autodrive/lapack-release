#line 1 "clatsqr.f"
/* clatsqr.f -- translated by f2c (version 20100827).
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

#line 1 "clatsqr.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, */
/*                           LWORK, INFO) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, M, N, MB, NB, LDT, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX           A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLATSQR computes a blocked Tall-Skinny QR factorization of */
/* > an M-by-N matrix A, where M >= N: */
/* > A = Q * R . */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The row block size to be used in the blocked QR. */
/* >          MB > N. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The column block size to be used in the blocked QR. */
/* >          N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal */
/* >          of the array contain the N-by-N upper triangular matrix R; */
/* >          the elements below the diagonal represent Q by the columns */
/* >          of blocked V (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, */
/* >          dimension (LDT, N * Number_of_row_blocks) */
/* >          where Number_of_row_blocks = CEIL((M-N)/(MB-N)) */
/* >          The blocked upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks. */
/* >          See Further Details below. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >         (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          The dimension of the array WORK.  LWORK >= NB*N. */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* >   Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out subdiagonal entries of a block of MB rows of A: */
/* >   Q(1) zeros out the subdiagonal entries of rows 1:MB of A */
/* >   Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A */
/* >   Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N] of A */
/* >   . . . */
/* > */
/* > Q(1) is computed by GEQRT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GEQRT. */
/* > */
/* > Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Householder vectors */
/* > stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N). */
/* > The last Q(k) may use fewer rows. */
/* > For more information see Further Details in TPQRT. */
/* > */
/* > For more details of the overall algorithm, see the description of */
/* > Sequential TSQR in Section 2.2 of [1]. */
/* > */
/* > [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations, */
/* >     J. Demmel, L. Grigori, M. Hoemmen, J. Langou, */
/* >     SIAM J. Sci. Comput, vol. 34, no. 1, 2012 */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clatsqr_(integer *m, integer *n, integer *mb, integer *
	nb, doublecomplex *a, integer *lda, doublecomplex *t, integer *ldt, 
	doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ii, kk, ctr;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cgeqrt_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ctpqrt_(
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *);
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. -- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. EXTERNAL FUNCTIONS .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*     .. */
/*     .. EXECUTABLE STATEMENTS .. */

/*     TEST THE INPUT ARGUMENTS */

#line 183 "clatsqr.f"
    /* Parameter adjustments */
#line 183 "clatsqr.f"
    a_dim1 = *lda;
#line 183 "clatsqr.f"
    a_offset = 1 + a_dim1;
#line 183 "clatsqr.f"
    a -= a_offset;
#line 183 "clatsqr.f"
    t_dim1 = *ldt;
#line 183 "clatsqr.f"
    t_offset = 1 + t_dim1;
#line 183 "clatsqr.f"
    t -= t_offset;
#line 183 "clatsqr.f"
    --work;
#line 183 "clatsqr.f"

#line 183 "clatsqr.f"
    /* Function Body */
#line 183 "clatsqr.f"
    *info = 0;

#line 185 "clatsqr.f"
    lquery = *lwork == -1;

#line 187 "clatsqr.f"
    if (*m < 0) {
#line 188 "clatsqr.f"
	*info = -1;
#line 189 "clatsqr.f"
    } else if (*n < 0 || *m < *n) {
#line 190 "clatsqr.f"
	*info = -2;
#line 191 "clatsqr.f"
    } else if (*mb <= *n) {
#line 192 "clatsqr.f"
	*info = -3;
#line 193 "clatsqr.f"
    } else if (*nb < 1 || *nb > *n && *n > 0) {
#line 194 "clatsqr.f"
	*info = -4;
#line 195 "clatsqr.f"
    } else if (*lda < max(1,*m)) {
#line 196 "clatsqr.f"
	*info = -5;
#line 197 "clatsqr.f"
    } else if (*ldt < *nb) {
#line 198 "clatsqr.f"
	*info = -8;
#line 199 "clatsqr.f"
    } else if (*lwork < *n * *nb && ! lquery) {
#line 200 "clatsqr.f"
	*info = -10;
#line 201 "clatsqr.f"
    }
#line 202 "clatsqr.f"
    if (*info == 0) {
#line 203 "clatsqr.f"
	i__1 = *nb * *n;
#line 203 "clatsqr.f"
	work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 204 "clatsqr.f"
    }
#line 205 "clatsqr.f"
    if (*info != 0) {
#line 206 "clatsqr.f"
	i__1 = -(*info);
#line 206 "clatsqr.f"
	xerbla_("CLATSQR", &i__1, (ftnlen)7);
#line 207 "clatsqr.f"
	return 0;
#line 208 "clatsqr.f"
    } else if (lquery) {
#line 209 "clatsqr.f"
	return 0;
#line 210 "clatsqr.f"
    }

/*     Quick return if possible */

#line 214 "clatsqr.f"
    if (min(*m,*n) == 0) {
#line 215 "clatsqr.f"
	return 0;
#line 216 "clatsqr.f"
    }

/*     The QR Decomposition */

#line 220 "clatsqr.f"
    if (*mb <= *n || *mb >= *m) {
#line 221 "clatsqr.f"
	cgeqrt_(m, n, nb, &a[a_offset], lda, &t[t_offset], ldt, &work[1], 
		info);
#line 222 "clatsqr.f"
	return 0;
#line 223 "clatsqr.f"
    }
#line 224 "clatsqr.f"
    kk = (*m - *n) % (*mb - *n);
#line 225 "clatsqr.f"
    ii = *m - kk + 1;

/*      Compute the QR factorization of the first block A(1:MB,1:N) */

#line 229 "clatsqr.f"
    cgeqrt_(mb, n, nb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &work[1], info)
	    ;
#line 230 "clatsqr.f"
    ctr = 1;

#line 232 "clatsqr.f"
    i__1 = ii - *mb + *n;
#line 232 "clatsqr.f"
    i__2 = *mb - *n;
#line 232 "clatsqr.f"
    for (i__ = *mb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*      Compute the QR factorization of the current block A(I:I+MB-N,1:N) */

#line 236 "clatsqr.f"
	i__3 = *mb - *n;
#line 236 "clatsqr.f"
	ctpqrt_(&i__3, n, &c__0, nb, &a[a_dim1 + 1], lda, &a[i__ + a_dim1], 
		lda, &t[(ctr * *n + 1) * t_dim1 + 1], ldt, &work[1], info);
#line 239 "clatsqr.f"
	++ctr;
#line 240 "clatsqr.f"
    }

/*      Compute the QR factorization of the last block A(II:M,1:N) */

#line 244 "clatsqr.f"
    if (ii <= *m) {
#line 245 "clatsqr.f"
	ctpqrt_(&kk, n, &c__0, nb, &a[a_dim1 + 1], lda, &a[ii + a_dim1], lda, 
		&t[(ctr * *n + 1) * t_dim1 + 1], ldt, &work[1], info);
#line 248 "clatsqr.f"
    }

#line 250 "clatsqr.f"
    i__2 = *n * *nb;
#line 250 "clatsqr.f"
    work[1].r = (doublereal) i__2, work[1].i = 0.;
#line 251 "clatsqr.f"
    return 0;

/*     End of CLATSQR */

} /* clatsqr_ */

