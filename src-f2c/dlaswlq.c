#line 1 "dlaswlq.f"
/* dlaswlq.f -- translated by f2c (version 20100827).
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

#line 1 "dlaswlq.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASWLQ( M, N, MB, NB, A, LDA, T, LDT, WORK, */
/*                            LWORK, INFO) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, M, N, MB, NB, LDT, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION  A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >          DLASWLQ computes a blocked Short-Wide LQ factorization of a */
/* >          M-by-N matrix A, where N >= M: */
/* >          A = L * Q */
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
/* >          The number of columns of the matrix A.  N >= M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The row block size to be used in the blocked QR. */
/* >          M >= MB >= 1 */
/* > \endverbatim */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The column block size to be used in the blocked QR. */
/* >          NB > M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and below the diagonal */
/* >          of the array contain the N-by-N lower triangular matrix L; */
/* >          the elements above the diagonal represent Q by the rows */
/* >          of blocked V (see Further Details). */
/* > */
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
/* >          T is DOUBLE PRECISION array, */
/* >          dimension (LDT, N * Number_of_row_blocks) */
/* >          where Number_of_row_blocks = CEIL((N-M)/(NB-M)) */
/* >          The blocked upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks. */
/* >          See Further Details below. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= MB. */
/* > \endverbatim */
/* > */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >         (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > */
/* > \endverbatim */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          The dimension of the array WORK.  LWORK >= MB*M. */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > */
/* > \endverbatim */
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
/* > Short-Wide LQ (SWLQ) performs LQ by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* >   Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out upper diagonal entries of a block of NB rows of A: */
/* >   Q(1) zeros out the upper diagonal entries of rows 1:NB of A */
/* >   Q(2) zeros out the bottom MB-N rows of rows [1:M,NB+1:2*NB-M] of A */
/* >   Q(3) zeros out the bottom MB-N rows of rows [1:M,2*NB-M+1:3*NB-2*M] of A */
/* >   . . . */
/* > */
/* > Q(1) is computed by GELQT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GELQT. */
/* > */
/* > Q(i) for i>1 is computed by TPLQT, which represents Q(i) by Householder vectors */
/* > stored in columns [(i-1)*(NB-M)+M+1:i*(NB-M)+M] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*M+1:i*M). */
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
/* Subroutine */ int dlaswlq_(integer *m, integer *n, integer *mb, integer *
	nb, doublereal *a, integer *lda, doublereal *t, integer *ldt, 
	doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ii, kk, ctr;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgelqt_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), dtplqt_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *);
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. -- */
/*     June 2017 */

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

#line 184 "dlaswlq.f"
    /* Parameter adjustments */
#line 184 "dlaswlq.f"
    a_dim1 = *lda;
#line 184 "dlaswlq.f"
    a_offset = 1 + a_dim1;
#line 184 "dlaswlq.f"
    a -= a_offset;
#line 184 "dlaswlq.f"
    t_dim1 = *ldt;
#line 184 "dlaswlq.f"
    t_offset = 1 + t_dim1;
#line 184 "dlaswlq.f"
    t -= t_offset;
#line 184 "dlaswlq.f"
    --work;
#line 184 "dlaswlq.f"

#line 184 "dlaswlq.f"
    /* Function Body */
#line 184 "dlaswlq.f"
    *info = 0;

#line 186 "dlaswlq.f"
    lquery = *lwork == -1;

#line 188 "dlaswlq.f"
    if (*m < 0) {
#line 189 "dlaswlq.f"
	*info = -1;
#line 190 "dlaswlq.f"
    } else if (*n < 0 || *n < *m) {
#line 191 "dlaswlq.f"
	*info = -2;
#line 192 "dlaswlq.f"
    } else if (*mb < 1 || *mb > *m && *m > 0) {
#line 193 "dlaswlq.f"
	*info = -3;
#line 194 "dlaswlq.f"
    } else if (*nb <= *m) {
#line 195 "dlaswlq.f"
	*info = -4;
#line 196 "dlaswlq.f"
    } else if (*lda < max(1,*m)) {
#line 197 "dlaswlq.f"
	*info = -5;
#line 198 "dlaswlq.f"
    } else if (*ldt < *mb) {
#line 199 "dlaswlq.f"
	*info = -8;
#line 200 "dlaswlq.f"
    } else if (*lwork < *m * *mb && ! lquery) {
#line 201 "dlaswlq.f"
	*info = -10;
#line 202 "dlaswlq.f"
    }
#line 203 "dlaswlq.f"
    if (*info == 0) {
#line 204 "dlaswlq.f"
	work[1] = (doublereal) (*mb * *m);
#line 205 "dlaswlq.f"
    }

#line 207 "dlaswlq.f"
    if (*info != 0) {
#line 208 "dlaswlq.f"
	i__1 = -(*info);
#line 208 "dlaswlq.f"
	xerbla_("DLASWLQ", &i__1, (ftnlen)7);
#line 209 "dlaswlq.f"
	return 0;
#line 210 "dlaswlq.f"
    } else if (lquery) {
#line 211 "dlaswlq.f"
	return 0;
#line 212 "dlaswlq.f"
    }

/*     Quick return if possible */

#line 216 "dlaswlq.f"
    if (min(*m,*n) == 0) {
#line 217 "dlaswlq.f"
	return 0;
#line 218 "dlaswlq.f"
    }

/*     The LQ Decomposition */

#line 222 "dlaswlq.f"
    if (*m >= *n || *nb <= *m || *nb >= *n) {
#line 223 "dlaswlq.f"
	dgelqt_(m, n, mb, &a[a_offset], lda, &t[t_offset], ldt, &work[1], 
		info);
#line 224 "dlaswlq.f"
	return 0;
#line 225 "dlaswlq.f"
    }

#line 227 "dlaswlq.f"
    kk = (*n - *m) % (*nb - *m);
#line 228 "dlaswlq.f"
    ii = *n - kk + 1;

/*      Compute the LQ factorization of the first block A(1:M,1:NB) */

#line 232 "dlaswlq.f"
    dgelqt_(m, nb, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &work[1], info)
	    ;
#line 233 "dlaswlq.f"
    ctr = 1;

#line 235 "dlaswlq.f"
    i__1 = ii - *nb + *m;
#line 235 "dlaswlq.f"
    i__2 = *nb - *m;
#line 235 "dlaswlq.f"
    for (i__ = *nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*      Compute the QR factorization of the current block A(1:M,I:I+NB-M) */

#line 239 "dlaswlq.f"
	i__3 = *nb - *m;
#line 239 "dlaswlq.f"
	dtplqt_(m, &i__3, &c__0, mb, &a[a_dim1 + 1], lda, &a[i__ * a_dim1 + 1]
		, lda, &t[(ctr * *m + 1) * t_dim1 + 1], ldt, &work[1], info);
#line 242 "dlaswlq.f"
	++ctr;
#line 243 "dlaswlq.f"
    }

/*     Compute the QR factorization of the last block A(1:M,II:N) */

#line 247 "dlaswlq.f"
    if (ii <= *n) {
#line 248 "dlaswlq.f"
	dtplqt_(m, &kk, &c__0, mb, &a[a_dim1 + 1], lda, &a[ii * a_dim1 + 1], 
		lda, &t[(ctr * *m + 1) * t_dim1 + 1], ldt, &work[1], info);
#line 251 "dlaswlq.f"
    }

#line 253 "dlaswlq.f"
    work[1] = (doublereal) (*m * *mb);
#line 254 "dlaswlq.f"
    return 0;

/*     End of DLASWLQ */

} /* dlaswlq_ */

