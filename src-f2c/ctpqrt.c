#line 1 "ctpqrt.f"
/* ctpqrt.f -- translated by f2c (version 20100827).
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

#line 1 "ctpqrt.f"
/* > \brief \b CTPQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPQRT computes a blocked QR factorization of a complex */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B, and the order of the */
/* >          triangular matrix A. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the upper trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the upper triangular N-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first M-L rows */
/* >          are rectangular, and the last L rows are upper trapezoidal. */
/* >          On exit, B contains the pentagonal matrix V.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,N) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See Further Details. */
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
/* >          WORK is COMPLEX array, dimension (NB*N) */
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

/* > \date November 2013 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a (N+M)-by-N matrix */
/* > */
/* >               C = [ A ] */
/* >                   [ B ] */
/* > */
/* >  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ]  <- (M-L)-by-N rectangular */
/* >                   [ B2 ]  <-     L-by-N upper trapezoidal. */
/* > */
/* >  The upper trapezoidal matrix B2 consists of the first L rows of a */
/* >  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is upper triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* >               C = [ A ]  <- upper triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               V = [ V1 ] <- (M-L)-by-N rectangular */
/* >                   [ V2 ] <-     L-by-N upper trapezoidal. */
/* > */
/* >  The columns of V represent the vectors which define the H(i)'s. */
/* > */
/* >  The number of blocks is B = ceiling(N/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctpqrt_(integer *m, integer *n, integer *l, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *t, integer *ldt, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, ib, lb, mb, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ctprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    ctpqrt2_(integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *);


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 217 "ctpqrt.f"
    /* Parameter adjustments */
#line 217 "ctpqrt.f"
    a_dim1 = *lda;
#line 217 "ctpqrt.f"
    a_offset = 1 + a_dim1;
#line 217 "ctpqrt.f"
    a -= a_offset;
#line 217 "ctpqrt.f"
    b_dim1 = *ldb;
#line 217 "ctpqrt.f"
    b_offset = 1 + b_dim1;
#line 217 "ctpqrt.f"
    b -= b_offset;
#line 217 "ctpqrt.f"
    t_dim1 = *ldt;
#line 217 "ctpqrt.f"
    t_offset = 1 + t_dim1;
#line 217 "ctpqrt.f"
    t -= t_offset;
#line 217 "ctpqrt.f"
    --work;
#line 217 "ctpqrt.f"

#line 217 "ctpqrt.f"
    /* Function Body */
#line 217 "ctpqrt.f"
    *info = 0;
#line 218 "ctpqrt.f"
    if (*m < 0) {
#line 219 "ctpqrt.f"
	*info = -1;
#line 220 "ctpqrt.f"
    } else if (*n < 0) {
#line 221 "ctpqrt.f"
	*info = -2;
#line 222 "ctpqrt.f"
    } else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0) {
#line 223 "ctpqrt.f"
	*info = -3;
#line 224 "ctpqrt.f"
    } else if (*nb < 1 || *nb > *n && *n > 0) {
#line 225 "ctpqrt.f"
	*info = -4;
#line 226 "ctpqrt.f"
    } else if (*lda < max(1,*n)) {
#line 227 "ctpqrt.f"
	*info = -6;
#line 228 "ctpqrt.f"
    } else if (*ldb < max(1,*m)) {
#line 229 "ctpqrt.f"
	*info = -8;
#line 230 "ctpqrt.f"
    } else if (*ldt < *nb) {
#line 231 "ctpqrt.f"
	*info = -10;
#line 232 "ctpqrt.f"
    }
#line 233 "ctpqrt.f"
    if (*info != 0) {
#line 234 "ctpqrt.f"
	i__1 = -(*info);
#line 234 "ctpqrt.f"
	xerbla_("CTPQRT", &i__1, (ftnlen)6);
#line 235 "ctpqrt.f"
	return 0;
#line 236 "ctpqrt.f"
    }

/*     Quick return if possible */

#line 240 "ctpqrt.f"
    if (*m == 0 || *n == 0) {
#line 240 "ctpqrt.f"
	return 0;
#line 240 "ctpqrt.f"
    }

#line 242 "ctpqrt.f"
    i__1 = *n;
#line 242 "ctpqrt.f"
    i__2 = *nb;
#line 242 "ctpqrt.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*     Compute the QR factorization of the current block */

/* Computing MIN */
#line 246 "ctpqrt.f"
	i__3 = *n - i__ + 1;
#line 246 "ctpqrt.f"
	ib = min(i__3,*nb);
/* Computing MIN */
#line 247 "ctpqrt.f"
	i__3 = *m - *l + i__ + ib - 1;
#line 247 "ctpqrt.f"
	mb = min(i__3,*m);
#line 248 "ctpqrt.f"
	if (i__ >= *l) {
#line 249 "ctpqrt.f"
	    lb = 0;
#line 250 "ctpqrt.f"
	} else {
#line 251 "ctpqrt.f"
	    lb = mb - *m + *l - i__ + 1;
#line 252 "ctpqrt.f"
	}

#line 254 "ctpqrt.f"
	ctpqrt2_(&mb, &ib, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 
		+ 1], ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);

/*     Update by applying H**H to B(:,I+IB:N) from the left */

#line 259 "ctpqrt.f"
	if (i__ + ib <= *n) {
#line 260 "ctpqrt.f"
	    i__3 = *n - i__ - ib + 1;
#line 260 "ctpqrt.f"
	    ctprfb_("L", "C", "F", "C", &mb, &i__3, &ib, &lb, &b[i__ * b_dim1 
		    + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) 
		    * a_dim1], lda, &b[(i__ + ib) * b_dim1 + 1], ldb, &work[1]
		    , &ib, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 264 "ctpqrt.f"
	}
#line 265 "ctpqrt.f"
    }
#line 266 "ctpqrt.f"
    return 0;

/*     End of CTPQRT */

} /* ctpqrt_ */

