#line 1 "dgeqrt.f"
/* dgeqrt.f -- translated by f2c (version 20100827).
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

#line 1 "dgeqrt.f"
/* > \brief \b DGEQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQRT computes a blocked QR factorization of a real M-by-N matrix A */
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
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if M >= N); the elements below the diagonal */
/* >          are the columns of V. */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,MIN(M,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (NB*N) */
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

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix V stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* >               V = (  1       ) */
/* >                   ( v1  1    ) */
/* >                   ( v1 v2  1 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  where the vi's represent the vectors which define H(i), which are returned */
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A. */
/* > */
/* >  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-K matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgeqrt_(integer *m, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *t, integer *ldt, doublereal *work, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgeqrt2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dgeqrt3_(integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 171 "dgeqrt.f"
    /* Parameter adjustments */
#line 171 "dgeqrt.f"
    a_dim1 = *lda;
#line 171 "dgeqrt.f"
    a_offset = 1 + a_dim1;
#line 171 "dgeqrt.f"
    a -= a_offset;
#line 171 "dgeqrt.f"
    t_dim1 = *ldt;
#line 171 "dgeqrt.f"
    t_offset = 1 + t_dim1;
#line 171 "dgeqrt.f"
    t -= t_offset;
#line 171 "dgeqrt.f"
    --work;
#line 171 "dgeqrt.f"

#line 171 "dgeqrt.f"
    /* Function Body */
#line 171 "dgeqrt.f"
    *info = 0;
#line 172 "dgeqrt.f"
    if (*m < 0) {
#line 173 "dgeqrt.f"
	*info = -1;
#line 174 "dgeqrt.f"
    } else if (*n < 0) {
#line 175 "dgeqrt.f"
	*info = -2;
#line 176 "dgeqrt.f"
    } else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
#line 177 "dgeqrt.f"
	*info = -3;
#line 178 "dgeqrt.f"
    } else if (*lda < max(1,*m)) {
#line 179 "dgeqrt.f"
	*info = -5;
#line 180 "dgeqrt.f"
    } else if (*ldt < *nb) {
#line 181 "dgeqrt.f"
	*info = -7;
#line 182 "dgeqrt.f"
    }
#line 183 "dgeqrt.f"
    if (*info != 0) {
#line 184 "dgeqrt.f"
	i__1 = -(*info);
#line 184 "dgeqrt.f"
	xerbla_("DGEQRT", &i__1, (ftnlen)6);
#line 185 "dgeqrt.f"
	return 0;
#line 186 "dgeqrt.f"
    }

/*     Quick return if possible */

#line 190 "dgeqrt.f"
    k = min(*m,*n);
#line 191 "dgeqrt.f"
    if (k == 0) {
#line 191 "dgeqrt.f"
	return 0;
#line 191 "dgeqrt.f"
    }

/*     Blocked loop of length K */

#line 195 "dgeqrt.f"
    i__1 = k;
#line 195 "dgeqrt.f"
    i__2 = *nb;
#line 195 "dgeqrt.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 196 "dgeqrt.f"
	i__3 = k - i__ + 1;
#line 196 "dgeqrt.f"
	ib = min(i__3,*nb);

/*     Compute the QR factorization of the current block A(I:M,I:I+IB-1) */

#line 200 "dgeqrt.f"
	if (TRUE_) {
#line 201 "dgeqrt.f"
	    i__3 = *m - i__ + 1;
#line 201 "dgeqrt.f"
	    dgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
#line 202 "dgeqrt.f"
	} else {
#line 203 "dgeqrt.f"
	    i__3 = *m - i__ + 1;
#line 203 "dgeqrt.f"
	    dgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
#line 204 "dgeqrt.f"
	}
#line 205 "dgeqrt.f"
	if (i__ + ib <= *n) {

/*     Update by applying H**T to A(I:M,I+IB:N) from the left */

#line 209 "dgeqrt.f"
	    i__3 = *m - i__ + 1;
#line 209 "dgeqrt.f"
	    i__4 = *n - i__ - ib + 1;
#line 209 "dgeqrt.f"
	    i__5 = *n - i__ - ib + 1;
#line 209 "dgeqrt.f"
	    dlarfb_("L", "T", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + 
		    ib) * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
#line 212 "dgeqrt.f"
	}
#line 213 "dgeqrt.f"
    }
#line 214 "dgeqrt.f"
    return 0;

/*     End of DGEQRT */

} /* dgeqrt_ */

