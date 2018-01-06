#line 1 "ctzrqf.f"
/* ctzrqf.f -- translated by f2c (version 20100827).
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

#line 1 "ctzrqf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CTZRQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTZRQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctzrqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctzrqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctzrqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTZRQF( M, N, A, LDA, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine CTZRZF. */
/* > */
/* > CTZRQF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A */
/* > to upper triangular form by means of unitary transformations. */
/* > */
/* > The upper trapezoidal matrix A is factored as */
/* > */
/* >    A = ( R  0 ) * Z, */
/* > */
/* > where Z is an N-by-N unitary matrix and R is an M-by-M upper */
/* > triangular matrix. */
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
/* >          The number of columns of the matrix A.  N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the leading M-by-N upper trapezoidal part of the */
/* >          array A must contain the matrix to be factorized. */
/* >          On exit, the leading M-by-M upper triangular part of A */
/* >          contains the upper triangular matrix R, and elements M+1 to */
/* >          N of the first M rows of A, with the array TAU, represent the */
/* >          unitary matrix Z as a product of M elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX array, dimension (M) */
/* >          The scalar factors of the elementary reflectors. */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The  factorization is obtained by Householder's method.  The kth */
/* >  transformation matrix, Z( k ), whose conjugate transpose is used to */
/* >  introduce zeros into the (m - k + 1)th row of A, is given in the form */
/* > */
/* >     Z( k ) = ( I     0   ), */
/* >              ( 0  T( k ) ) */
/* > */
/* >  where */
/* > */
/* >     T( k ) = I - tau*u( k )*u( k )**H,   u( k ) = (   1    ), */
/* >                                                   (   0    ) */
/* >                                                   ( z( k ) ) */
/* > */
/* >  tau is a scalar and z( k ) is an ( n - m ) element vector. */
/* >  tau and z( k ) are chosen to annihilate the elements of the kth row */
/* >  of X. */
/* > */
/* >  The scalar tau is returned in the kth element of TAU and the vector */
/* >  u( k ) in the kth row of A, such that the elements of z( k ) are */
/* >  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in */
/* >  the upper triangular part of A. */
/* > */
/* >  Z is given by */
/* > */
/* >     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctzrqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k, m1;
    extern /* Subroutine */ int cgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublecomplex alpha;
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), caxpy_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), clarfg_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *), 
	    clacgv_(integer *, doublecomplex *, integer *), xerbla_(char *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 175 "ctzrqf.f"
    /* Parameter adjustments */
#line 175 "ctzrqf.f"
    a_dim1 = *lda;
#line 175 "ctzrqf.f"
    a_offset = 1 + a_dim1;
#line 175 "ctzrqf.f"
    a -= a_offset;
#line 175 "ctzrqf.f"
    --tau;
#line 175 "ctzrqf.f"

#line 175 "ctzrqf.f"
    /* Function Body */
#line 175 "ctzrqf.f"
    *info = 0;
#line 176 "ctzrqf.f"
    if (*m < 0) {
#line 177 "ctzrqf.f"
	*info = -1;
#line 178 "ctzrqf.f"
    } else if (*n < *m) {
#line 179 "ctzrqf.f"
	*info = -2;
#line 180 "ctzrqf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "ctzrqf.f"
	*info = -4;
#line 182 "ctzrqf.f"
    }
#line 183 "ctzrqf.f"
    if (*info != 0) {
#line 184 "ctzrqf.f"
	i__1 = -(*info);
#line 184 "ctzrqf.f"
	xerbla_("CTZRQF", &i__1, (ftnlen)6);
#line 185 "ctzrqf.f"
	return 0;
#line 186 "ctzrqf.f"
    }

/*     Perform the factorization. */

#line 190 "ctzrqf.f"
    if (*m == 0) {
#line 190 "ctzrqf.f"
	return 0;
#line 190 "ctzrqf.f"
    }
#line 192 "ctzrqf.f"
    if (*m == *n) {
#line 193 "ctzrqf.f"
	i__1 = *n;
#line 193 "ctzrqf.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "ctzrqf.f"
	    i__2 = i__;
#line 194 "ctzrqf.f"
	    tau[i__2].r = 0., tau[i__2].i = 0.;
#line 195 "ctzrqf.f"
/* L10: */
#line 195 "ctzrqf.f"
	}
#line 196 "ctzrqf.f"
    } else {
/* Computing MIN */
#line 197 "ctzrqf.f"
	i__1 = *m + 1;
#line 197 "ctzrqf.f"
	m1 = min(i__1,*n);
#line 198 "ctzrqf.f"
	for (k = *m; k >= 1; --k) {

/*           Use a Householder reflection to zero the kth row of A. */
/*           First set up the reflection. */

#line 203 "ctzrqf.f"
	    i__1 = k + k * a_dim1;
#line 203 "ctzrqf.f"
	    d_cnjg(&z__1, &a[k + k * a_dim1]);
#line 203 "ctzrqf.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 204 "ctzrqf.f"
	    i__1 = *n - *m;
#line 204 "ctzrqf.f"
	    clacgv_(&i__1, &a[k + m1 * a_dim1], lda);
#line 205 "ctzrqf.f"
	    i__1 = k + k * a_dim1;
#line 205 "ctzrqf.f"
	    alpha.r = a[i__1].r, alpha.i = a[i__1].i;
#line 206 "ctzrqf.f"
	    i__1 = *n - *m + 1;
#line 206 "ctzrqf.f"
	    clarfg_(&i__1, &alpha, &a[k + m1 * a_dim1], lda, &tau[k]);
#line 207 "ctzrqf.f"
	    i__1 = k + k * a_dim1;
#line 207 "ctzrqf.f"
	    a[i__1].r = alpha.r, a[i__1].i = alpha.i;
#line 208 "ctzrqf.f"
	    i__1 = k;
#line 208 "ctzrqf.f"
	    d_cnjg(&z__1, &tau[k]);
#line 208 "ctzrqf.f"
	    tau[i__1].r = z__1.r, tau[i__1].i = z__1.i;

#line 210 "ctzrqf.f"
	    i__1 = k;
#line 210 "ctzrqf.f"
	    if ((tau[i__1].r != 0. || tau[i__1].i != 0.) && k > 1) {

/*              We now perform the operation  A := A*P( k )**H. */

/*              Use the first ( k - 1 ) elements of TAU to store  a( k ), */
/*              where  a( k ) consists of the first ( k - 1 ) elements of */
/*              the  kth column  of  A.  Also  let  B  denote  the  first */
/*              ( k - 1 ) rows of the last ( n - m ) columns of A. */

#line 219 "ctzrqf.f"
		i__1 = k - 1;
#line 219 "ctzrqf.f"
		ccopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &tau[1], &c__1);

/*              Form   w = a( k ) + B*z( k )  in TAU. */

#line 223 "ctzrqf.f"
		i__1 = k - 1;
#line 223 "ctzrqf.f"
		i__2 = *n - *m;
#line 223 "ctzrqf.f"
		cgemv_("No transpose", &i__1, &i__2, &c_b1, &a[m1 * a_dim1 + 
			1], lda, &a[k + m1 * a_dim1], lda, &c_b1, &tau[1], &
			c__1, (ftnlen)12);

/*              Now form  a( k ) := a( k ) - conjg(tau)*w */
/*              and       B      := B      - conjg(tau)*w*z( k )**H. */

#line 229 "ctzrqf.f"
		i__1 = k - 1;
#line 229 "ctzrqf.f"
		d_cnjg(&z__2, &tau[k]);
#line 229 "ctzrqf.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 229 "ctzrqf.f"
		caxpy_(&i__1, &z__1, &tau[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 231 "ctzrqf.f"
		i__1 = k - 1;
#line 231 "ctzrqf.f"
		i__2 = *n - *m;
#line 231 "ctzrqf.f"
		d_cnjg(&z__2, &tau[k]);
#line 231 "ctzrqf.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 231 "ctzrqf.f"
		cgerc_(&i__1, &i__2, &z__1, &tau[1], &c__1, &a[k + m1 * 
			a_dim1], lda, &a[m1 * a_dim1 + 1], lda);
#line 233 "ctzrqf.f"
	    }
#line 234 "ctzrqf.f"
/* L20: */
#line 234 "ctzrqf.f"
	}
#line 235 "ctzrqf.f"
    }

#line 237 "ctzrqf.f"
    return 0;

/*     End of CTZRQF */

} /* ctzrqf_ */

