#line 1 "cgeqrt.f"
/* cgeqrt.f -- translated by f2c (version 20100827).
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

#line 1 "cgeqrt.f"
/* > \brief \b CGEQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQRT computes a blocked QR factorization of a complex M-by-N matrix A */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          T is COMPLEX array, dimension (LDT,MIN(M,N)) */
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

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

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
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgeqrt_(integer *m, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *t, integer *ldt, 
	doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), cgeqrt2_(integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    ), cgeqrt3_(integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *);


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

#line 171 "cgeqrt.f"
    /* Parameter adjustments */
#line 171 "cgeqrt.f"
    a_dim1 = *lda;
#line 171 "cgeqrt.f"
    a_offset = 1 + a_dim1;
#line 171 "cgeqrt.f"
    a -= a_offset;
#line 171 "cgeqrt.f"
    t_dim1 = *ldt;
#line 171 "cgeqrt.f"
    t_offset = 1 + t_dim1;
#line 171 "cgeqrt.f"
    t -= t_offset;
#line 171 "cgeqrt.f"
    --work;
#line 171 "cgeqrt.f"

#line 171 "cgeqrt.f"
    /* Function Body */
#line 171 "cgeqrt.f"
    *info = 0;
#line 172 "cgeqrt.f"
    if (*m < 0) {
#line 173 "cgeqrt.f"
	*info = -1;
#line 174 "cgeqrt.f"
    } else if (*n < 0) {
#line 175 "cgeqrt.f"
	*info = -2;
#line 176 "cgeqrt.f"
    } else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
#line 177 "cgeqrt.f"
	*info = -3;
#line 178 "cgeqrt.f"
    } else if (*lda < max(1,*m)) {
#line 179 "cgeqrt.f"
	*info = -5;
#line 180 "cgeqrt.f"
    } else if (*ldt < *nb) {
#line 181 "cgeqrt.f"
	*info = -7;
#line 182 "cgeqrt.f"
    }
#line 183 "cgeqrt.f"
    if (*info != 0) {
#line 184 "cgeqrt.f"
	i__1 = -(*info);
#line 184 "cgeqrt.f"
	xerbla_("CGEQRT", &i__1, (ftnlen)6);
#line 185 "cgeqrt.f"
	return 0;
#line 186 "cgeqrt.f"
    }

/*     Quick return if possible */

#line 190 "cgeqrt.f"
    k = min(*m,*n);
#line 191 "cgeqrt.f"
    if (k == 0) {
#line 191 "cgeqrt.f"
	return 0;
#line 191 "cgeqrt.f"
    }

/*     Blocked loop of length K */

#line 195 "cgeqrt.f"
    i__1 = k;
#line 195 "cgeqrt.f"
    i__2 = *nb;
#line 195 "cgeqrt.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 196 "cgeqrt.f"
	i__3 = k - i__ + 1;
#line 196 "cgeqrt.f"
	ib = min(i__3,*nb);

/*     Compute the QR factorization of the current block A(I:M,I:I+IB-1) */

#line 200 "cgeqrt.f"
	if (TRUE_) {
#line 201 "cgeqrt.f"
	    i__3 = *m - i__ + 1;
#line 201 "cgeqrt.f"
	    cgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
#line 202 "cgeqrt.f"
	} else {
#line 203 "cgeqrt.f"
	    i__3 = *m - i__ + 1;
#line 203 "cgeqrt.f"
	    cgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
#line 204 "cgeqrt.f"
	}
#line 205 "cgeqrt.f"
	if (i__ + ib <= *n) {

/*     Update by applying H**H to A(I:M,I+IB:N) from the left */

#line 209 "cgeqrt.f"
	    i__3 = *m - i__ + 1;
#line 209 "cgeqrt.f"
	    i__4 = *n - i__ - ib + 1;
#line 209 "cgeqrt.f"
	    i__5 = *n - i__ - ib + 1;
#line 209 "cgeqrt.f"
	    clarfb_("L", "C", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + 
		    ib) * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
#line 212 "cgeqrt.f"
	}
#line 213 "cgeqrt.f"
    }
#line 214 "cgeqrt.f"
    return 0;

/*     End of CGEQRT */

} /* cgeqrt_ */

