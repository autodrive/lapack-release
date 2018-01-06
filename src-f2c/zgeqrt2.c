#line 1 "zgeqrt2.f"
/* zgeqrt2.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqrt2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZGEQRT2 computes a QR factorization of a general real or complex matrix using the compact WY re
presentation of Q. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQRT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrt2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrt2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrt2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQRT2( M, N, A, LDA, T, LDT, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER   INFO, LDA, LDT, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16   A( LDA, * ), T( LDT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQRT2 computes a QR factorization of a complex M-by-N matrix A, */
/* > using the compact WY representation of Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the complex M-by-N matrix A.  On exit, the elements on and */
/* >          above the diagonal contain the N-by-N upper triangular matrix R; the */
/* >          elements below the diagonal are the columns of V.  See below for */
/* >          further details. */
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
/* >          T is COMPLEX*16 array, dimension (LDT,N) */
/* >          The N-by-N upper triangular factor of the block reflector. */
/* >          The elements on and above the diagonal contain the block */
/* >          reflector T; the elements below the diagonal are not used. */
/* >          See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16GEcomputational */

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
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The */
/* >  block reflector H is then given by */
/* > */
/* >               H = I - V * T * V**H */
/* > */
/* >  where V**H is the conjugate transpose of V. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqrt2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *t, integer *ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k;
    static doublecomplex aii, alpha;
    extern /* Subroutine */ int zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ztrmv_(char *, char *, char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), zlarfg_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 159 "zgeqrt2.f"
    /* Parameter adjustments */
#line 159 "zgeqrt2.f"
    a_dim1 = *lda;
#line 159 "zgeqrt2.f"
    a_offset = 1 + a_dim1;
#line 159 "zgeqrt2.f"
    a -= a_offset;
#line 159 "zgeqrt2.f"
    t_dim1 = *ldt;
#line 159 "zgeqrt2.f"
    t_offset = 1 + t_dim1;
#line 159 "zgeqrt2.f"
    t -= t_offset;
#line 159 "zgeqrt2.f"

#line 159 "zgeqrt2.f"
    /* Function Body */
#line 159 "zgeqrt2.f"
    *info = 0;
#line 160 "zgeqrt2.f"
    if (*m < 0) {
#line 161 "zgeqrt2.f"
	*info = -1;
#line 162 "zgeqrt2.f"
    } else if (*n < 0) {
#line 163 "zgeqrt2.f"
	*info = -2;
#line 164 "zgeqrt2.f"
    } else if (*lda < max(1,*m)) {
#line 165 "zgeqrt2.f"
	*info = -4;
#line 166 "zgeqrt2.f"
    } else if (*ldt < max(1,*n)) {
#line 167 "zgeqrt2.f"
	*info = -6;
#line 168 "zgeqrt2.f"
    }
#line 169 "zgeqrt2.f"
    if (*info != 0) {
#line 170 "zgeqrt2.f"
	i__1 = -(*info);
#line 170 "zgeqrt2.f"
	xerbla_("ZGEQRT2", &i__1, (ftnlen)7);
#line 171 "zgeqrt2.f"
	return 0;
#line 172 "zgeqrt2.f"
    }

#line 174 "zgeqrt2.f"
    k = min(*m,*n);

#line 176 "zgeqrt2.f"
    i__1 = k;
#line 176 "zgeqrt2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1) */

#line 180 "zgeqrt2.f"
	i__2 = *m - i__ + 1;
/* Computing MIN */
#line 180 "zgeqrt2.f"
	i__3 = i__ + 1;
#line 180 "zgeqrt2.f"
	zlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * a_dim1]
		, &c__1, &t[i__ + t_dim1]);
#line 182 "zgeqrt2.f"
	if (i__ < *n) {

/*           Apply H(i) to A(I:M,I+1:N) from the left */

#line 186 "zgeqrt2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 186 "zgeqrt2.f"
	    aii.r = a[i__2].r, aii.i = a[i__2].i;
#line 187 "zgeqrt2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 187 "zgeqrt2.f"
	    a[i__2].r = 1., a[i__2].i = 0.;

/*           W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)] */

#line 191 "zgeqrt2.f"
	    i__2 = *m - i__ + 1;
#line 191 "zgeqrt2.f"
	    i__3 = *n - i__;
#line 191 "zgeqrt2.f"
	    zgemv_("C", &i__2, &i__3, &c_b1, &a[i__ + (i__ + 1) * a_dim1], 
		    lda, &a[i__ + i__ * a_dim1], &c__1, &c_b2, &t[*n * t_dim1 
		    + 1], &c__1, (ftnlen)1);

/*           A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H */

#line 196 "zgeqrt2.f"
	    d_cnjg(&z__2, &t[i__ + t_dim1]);
#line 196 "zgeqrt2.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 196 "zgeqrt2.f"
	    alpha.r = z__1.r, alpha.i = z__1.i;
#line 197 "zgeqrt2.f"
	    i__2 = *m - i__ + 1;
#line 197 "zgeqrt2.f"
	    i__3 = *n - i__;
#line 197 "zgeqrt2.f"
	    zgerc_(&i__2, &i__3, &alpha, &a[i__ + i__ * a_dim1], &c__1, &t[*n 
		    * t_dim1 + 1], &c__1, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 199 "zgeqrt2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 199 "zgeqrt2.f"
	    a[i__2].r = aii.r, a[i__2].i = aii.i;
#line 200 "zgeqrt2.f"
	}
#line 201 "zgeqrt2.f"
    }

#line 203 "zgeqrt2.f"
    i__1 = *n;
#line 203 "zgeqrt2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 204 "zgeqrt2.f"
	i__2 = i__ + i__ * a_dim1;
#line 204 "zgeqrt2.f"
	aii.r = a[i__2].r, aii.i = a[i__2].i;
#line 205 "zgeqrt2.f"
	i__2 = i__ + i__ * a_dim1;
#line 205 "zgeqrt2.f"
	a[i__2].r = 1., a[i__2].i = 0.;

/*        T(1:I-1,I) := alpha * A(I:M,1:I-1)**H * A(I:M,I) */

#line 209 "zgeqrt2.f"
	i__2 = i__ + t_dim1;
#line 209 "zgeqrt2.f"
	z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
#line 209 "zgeqrt2.f"
	alpha.r = z__1.r, alpha.i = z__1.i;
#line 210 "zgeqrt2.f"
	i__2 = *m - i__ + 1;
#line 210 "zgeqrt2.f"
	i__3 = i__ - 1;
#line 210 "zgeqrt2.f"
	zgemv_("C", &i__2, &i__3, &alpha, &a[i__ + a_dim1], lda, &a[i__ + i__ 
		* a_dim1], &c__1, &c_b2, &t[i__ * t_dim1 + 1], &c__1, (ftnlen)
		1);
#line 212 "zgeqrt2.f"
	i__2 = i__ + i__ * a_dim1;
#line 212 "zgeqrt2.f"
	a[i__2].r = aii.r, a[i__2].i = aii.i;

/*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I) */

#line 216 "zgeqrt2.f"
	i__2 = i__ - 1;
#line 216 "zgeqrt2.f"
	ztrmv_("U", "N", "N", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], 
		&c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*           T(I,I) = tau(I) */

#line 220 "zgeqrt2.f"
	i__2 = i__ + i__ * t_dim1;
#line 220 "zgeqrt2.f"
	i__3 = i__ + t_dim1;
#line 220 "zgeqrt2.f"
	t[i__2].r = t[i__3].r, t[i__2].i = t[i__3].i;
#line 221 "zgeqrt2.f"
	i__2 = i__ + t_dim1;
#line 221 "zgeqrt2.f"
	t[i__2].r = 0., t[i__2].i = 0.;
#line 222 "zgeqrt2.f"
    }

/*     End of ZGEQRT2 */

#line 227 "zgeqrt2.f"
    return 0;
} /* zgeqrt2_ */

