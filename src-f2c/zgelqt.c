#line 1 "zgelqt.f"
/* zgelqt.f -- translated by f2c (version 20100827).
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

#line 1 "zgelqt.f"
/* > \brief \b ZGELQT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelqt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelqt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelqt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, MB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGELQT computes a blocked LQ factorization of a complex M-by-N matrix A */
/* > using the compact WY representation of Q. */
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
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The block size to be used in the blocked QR.  MIN(M,N) >= MB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the M-by-MIN(M,N) lower trapezoidal matrix L (L is */
/* >          lower triangular if M <= N); the elements above the diagonal */
/* >          are the rows of V. */
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
/* >          T is COMPLEX*16 array, dimension (LDT,MIN(M,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
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
/* >          WORK is COMPLEX*16 array, dimension (MB*N) */
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

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix V stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* >               V = (  1  v1 v1 v1 v1 ) */
/* >                   (     1  v2 v2 v2 ) */
/* >                   (         1 v3 v3 ) */
/* > */
/* > */
/* >  where the vi's represent the vectors which define H(i), which are returned */
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A. */
/* >  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgelqt_(integer *m, integer *n, integer *mb, 
	doublecomplex *a, integer *lda, doublecomplex *t, integer *ldt, 
	doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), zgelqt3_(integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *)
	    ;


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

#line 167 "zgelqt.f"
    /* Parameter adjustments */
#line 167 "zgelqt.f"
    a_dim1 = *lda;
#line 167 "zgelqt.f"
    a_offset = 1 + a_dim1;
#line 167 "zgelqt.f"
    a -= a_offset;
#line 167 "zgelqt.f"
    t_dim1 = *ldt;
#line 167 "zgelqt.f"
    t_offset = 1 + t_dim1;
#line 167 "zgelqt.f"
    t -= t_offset;
#line 167 "zgelqt.f"
    --work;
#line 167 "zgelqt.f"

#line 167 "zgelqt.f"
    /* Function Body */
#line 167 "zgelqt.f"
    *info = 0;
#line 168 "zgelqt.f"
    if (*m < 0) {
#line 169 "zgelqt.f"
	*info = -1;
#line 170 "zgelqt.f"
    } else if (*n < 0) {
#line 171 "zgelqt.f"
	*info = -2;
#line 172 "zgelqt.f"
    } else if (*mb < 1 || *mb > min(*m,*n) && min(*m,*n) > 0) {
#line 173 "zgelqt.f"
	*info = -3;
#line 174 "zgelqt.f"
    } else if (*lda < max(1,*m)) {
#line 175 "zgelqt.f"
	*info = -5;
#line 176 "zgelqt.f"
    } else if (*ldt < *mb) {
#line 177 "zgelqt.f"
	*info = -7;
#line 178 "zgelqt.f"
    }
#line 179 "zgelqt.f"
    if (*info != 0) {
#line 180 "zgelqt.f"
	i__1 = -(*info);
#line 180 "zgelqt.f"
	xerbla_("ZGELQT", &i__1, (ftnlen)6);
#line 181 "zgelqt.f"
	return 0;
#line 182 "zgelqt.f"
    }

/*     Quick return if possible */

#line 186 "zgelqt.f"
    k = min(*m,*n);
#line 187 "zgelqt.f"
    if (k == 0) {
#line 187 "zgelqt.f"
	return 0;
#line 187 "zgelqt.f"
    }

/*     Blocked loop of length K */

#line 191 "zgelqt.f"
    i__1 = k;
#line 191 "zgelqt.f"
    i__2 = *mb;
#line 191 "zgelqt.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 192 "zgelqt.f"
	i__3 = k - i__ + 1;
#line 192 "zgelqt.f"
	ib = min(i__3,*mb);

/*     Compute the LQ factorization of the current block A(I:M,I:I+IB-1) */

#line 196 "zgelqt.f"
	i__3 = *n - i__ + 1;
#line 196 "zgelqt.f"
	zgelqt3_(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 + 1]
		, ldt, &iinfo);
#line 197 "zgelqt.f"
	if (i__ + ib <= *m) {

/*     Update by applying H**T to A(I:M,I+IB:N) from the right */

#line 201 "zgelqt.f"
	    i__3 = *m - i__ - ib + 1;
#line 201 "zgelqt.f"
	    i__4 = *n - i__ + 1;
#line 201 "zgelqt.f"
	    i__5 = *m - i__ - ib + 1;
#line 201 "zgelqt.f"
	    zlarfb_("R", "N", "F", "R", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + ib + 
		    i__ * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
#line 204 "zgelqt.f"
	}
#line 205 "zgelqt.f"
    }
#line 206 "zgelqt.f"
    return 0;

/*     End of ZGELQT */

} /* zgelqt_ */

