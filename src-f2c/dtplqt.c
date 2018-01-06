#line 1 "dtplqt.f"
/* dtplqt.f -- translated by f2c (version 20100827).
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

#line 1 "dtplqt.f"
/* > \brief \b DTPLQT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtplqt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtplqt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtplqt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPLQT( M, N, L, MB, A, LDA, B, LDB, T, LDT, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, LDB, LDT, N, M, L, MB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION  A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPLQT computes a blocked LQ factorization of a real */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B, and the order of the */
/* >          triangular matrix A. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the lower trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The block size to be used in the blocked QR.  M >= MB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the lower triangular N-by-N matrix A. */
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the lower triangular matrix L. */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first N-L columns */
/* >          are rectangular, and the last L columns are lower trapezoidal. */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
/* >          The lower triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= MB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MB*M) */
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

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a M-by-(M+N) matrix */
/* > */
/* >               C = [ A ] [ B ] */
/* > */
/* > */
/* >  where A is an lower triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a M-by-(N-L) rectangular matrix B1 on left of a M-by-L */
/* >  upper trapezoidal matrix B2: */
/* >          [ B ] = [ B1 ] [ B2 ] */
/* >                   [ B1 ]  <- M-by-(N-L) rectangular */
/* >                   [ B2 ]  <-     M-by-L upper trapezoidal. */
/* > */
/* >  The lower trapezoidal matrix B2 consists of the first L columns of a */
/* >  N-by-N lower triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is lower triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th row */
/* >  above the diagonal (of A) in the M-by-(M+N) input matrix C */
/* >            [ C ] = [ A ] [ B ] */
/* >                   [ A ]  <- lower triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* >            [ W ] = [ I ] [ V ] */
/* >                   [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* >            [ V ] = [ V1 ] [ V2 ] */
/* >                   [ V1 ] <- M-by-(N-L) rectangular */
/* >                   [ V2 ] <-     M-by-L lower trapezoidal. */
/* > */
/* >  The rows of V represent the vectors which define the H(i)'s. */
/* > */
/* >  The number of blocks is B = ceiling(M/MB), where each */
/* >  block is of order MB except for the last block, which is of order */
/* >  IB = M - (M-1)*MB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The MB-by-MB (and IB-by-IB */
/* >  for the last block) T's are stored in the MB-by-N matrix T as */
/* > */
/* >               T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtplqt_(integer *m, integer *n, integer *l, integer *mb, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	t, integer *ldt, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, nb, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dtprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtplqt2_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 217 "dtplqt.f"
    /* Parameter adjustments */
#line 217 "dtplqt.f"
    a_dim1 = *lda;
#line 217 "dtplqt.f"
    a_offset = 1 + a_dim1;
#line 217 "dtplqt.f"
    a -= a_offset;
#line 217 "dtplqt.f"
    b_dim1 = *ldb;
#line 217 "dtplqt.f"
    b_offset = 1 + b_dim1;
#line 217 "dtplqt.f"
    b -= b_offset;
#line 217 "dtplqt.f"
    t_dim1 = *ldt;
#line 217 "dtplqt.f"
    t_offset = 1 + t_dim1;
#line 217 "dtplqt.f"
    t -= t_offset;
#line 217 "dtplqt.f"
    --work;
#line 217 "dtplqt.f"

#line 217 "dtplqt.f"
    /* Function Body */
#line 217 "dtplqt.f"
    *info = 0;
#line 218 "dtplqt.f"
    if (*m < 0) {
#line 219 "dtplqt.f"
	*info = -1;
#line 220 "dtplqt.f"
    } else if (*n < 0) {
#line 221 "dtplqt.f"
	*info = -2;
#line 222 "dtplqt.f"
    } else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0) {
#line 223 "dtplqt.f"
	*info = -3;
#line 224 "dtplqt.f"
    } else if (*mb < 1 || *mb > *m && *m > 0) {
#line 225 "dtplqt.f"
	*info = -4;
#line 226 "dtplqt.f"
    } else if (*lda < max(1,*m)) {
#line 227 "dtplqt.f"
	*info = -6;
#line 228 "dtplqt.f"
    } else if (*ldb < max(1,*m)) {
#line 229 "dtplqt.f"
	*info = -8;
#line 230 "dtplqt.f"
    } else if (*ldt < *mb) {
#line 231 "dtplqt.f"
	*info = -10;
#line 232 "dtplqt.f"
    }
#line 233 "dtplqt.f"
    if (*info != 0) {
#line 234 "dtplqt.f"
	i__1 = -(*info);
#line 234 "dtplqt.f"
	xerbla_("DTPLQT", &i__1, (ftnlen)6);
#line 235 "dtplqt.f"
	return 0;
#line 236 "dtplqt.f"
    }

/*     Quick return if possible */

#line 240 "dtplqt.f"
    if (*m == 0 || *n == 0) {
#line 240 "dtplqt.f"
	return 0;
#line 240 "dtplqt.f"
    }

#line 242 "dtplqt.f"
    i__1 = *m;
#line 242 "dtplqt.f"
    i__2 = *mb;
#line 242 "dtplqt.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*     Compute the QR factorization of the current block */

/* Computing MIN */
#line 246 "dtplqt.f"
	i__3 = *m - i__ + 1;
#line 246 "dtplqt.f"
	ib = min(i__3,*mb);
/* Computing MIN */
#line 247 "dtplqt.f"
	i__3 = *n - *l + i__ + ib - 1;
#line 247 "dtplqt.f"
	nb = min(i__3,*n);
#line 248 "dtplqt.f"
	if (i__ >= *l) {
#line 249 "dtplqt.f"
	    lb = 0;
#line 250 "dtplqt.f"
	} else {
#line 251 "dtplqt.f"
	    lb = nb - *n + *l - i__ + 1;
#line 252 "dtplqt.f"
	}

#line 254 "dtplqt.f"
	dtplqt2_(&ib, &nb, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ + b_dim1],
		 ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);

/*     Update by applying H**T to B(I+IB:M,:) from the right */

#line 259 "dtplqt.f"
	if (i__ + ib <= *m) {
#line 260 "dtplqt.f"
	    i__3 = *m - i__ - ib + 1;
#line 260 "dtplqt.f"
	    i__4 = *m - i__ - ib + 1;
#line 260 "dtplqt.f"
	    dtprfb_("R", "N", "F", "R", &i__3, &nb, &ib, &lb, &b[i__ + b_dim1]
		    , ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + ib + i__ * 
		    a_dim1], lda, &b[i__ + ib + b_dim1], ldb, &work[1], &i__4,
		     (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 264 "dtplqt.f"
	}
#line 265 "dtplqt.f"
    }
#line 266 "dtplqt.f"
    return 0;

/*     End of DTPLQT */

} /* dtplqt_ */

