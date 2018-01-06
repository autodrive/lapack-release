#line 1 "dtzrqf.f"
/* dtzrqf.f -- translated by f2c (version 20100827).
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

#line 1 "dtzrqf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = 1.;

/* > \brief \b DTZRQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTZRQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtzrqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtzrqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtzrqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTZRQF( M, N, A, LDA, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DTZRZF. */
/* > */
/* > DTZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A */
/* > to upper triangular form by means of orthogonal transformations. */
/* > */
/* > The upper trapezoidal matrix A is factored as */
/* > */
/* >    A = ( R  0 ) * Z, */
/* > */
/* > where Z is an N-by-N orthogonal matrix and R is an M-by-M upper */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the leading M-by-N upper trapezoidal part of the */
/* >          array A must contain the matrix to be factorized. */
/* >          On exit, the leading M-by-M upper triangular part of A */
/* >          contains the upper triangular matrix R, and elements M+1 to */
/* >          N of the first M rows of A, with the array TAU, represent the */
/* >          orthogonal matrix Z as a product of M elementary reflectors. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (M) */
/* >          The scalar factors of the elementary reflectors. */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The factorization is obtained by Householder's method.  The kth */
/* >  transformation matrix, Z( k ), which is used to introduce zeros into */
/* >  the ( m - k + 1 )th row of A, is given in the form */
/* > */
/* >     Z( k ) = ( I     0   ), */
/* >              ( 0  T( k ) ) */
/* > */
/* >  where */
/* > */
/* >     T( k ) = I - tau*u( k )*u( k )**T,   u( k ) = (   1    ), */
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
/* Subroutine */ int dtzrqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, m1;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *), dlarfg_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *), 
	    xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 172 "dtzrqf.f"
    /* Parameter adjustments */
#line 172 "dtzrqf.f"
    a_dim1 = *lda;
#line 172 "dtzrqf.f"
    a_offset = 1 + a_dim1;
#line 172 "dtzrqf.f"
    a -= a_offset;
#line 172 "dtzrqf.f"
    --tau;
#line 172 "dtzrqf.f"

#line 172 "dtzrqf.f"
    /* Function Body */
#line 172 "dtzrqf.f"
    *info = 0;
#line 173 "dtzrqf.f"
    if (*m < 0) {
#line 174 "dtzrqf.f"
	*info = -1;
#line 175 "dtzrqf.f"
    } else if (*n < *m) {
#line 176 "dtzrqf.f"
	*info = -2;
#line 177 "dtzrqf.f"
    } else if (*lda < max(1,*m)) {
#line 178 "dtzrqf.f"
	*info = -4;
#line 179 "dtzrqf.f"
    }
#line 180 "dtzrqf.f"
    if (*info != 0) {
#line 181 "dtzrqf.f"
	i__1 = -(*info);
#line 181 "dtzrqf.f"
	xerbla_("DTZRQF", &i__1, (ftnlen)6);
#line 182 "dtzrqf.f"
	return 0;
#line 183 "dtzrqf.f"
    }

/*     Perform the factorization. */

#line 187 "dtzrqf.f"
    if (*m == 0) {
#line 187 "dtzrqf.f"
	return 0;
#line 187 "dtzrqf.f"
    }
#line 189 "dtzrqf.f"
    if (*m == *n) {
#line 190 "dtzrqf.f"
	i__1 = *n;
#line 190 "dtzrqf.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "dtzrqf.f"
	    tau[i__] = 0.;
#line 192 "dtzrqf.f"
/* L10: */
#line 192 "dtzrqf.f"
	}
#line 193 "dtzrqf.f"
    } else {
/* Computing MIN */
#line 194 "dtzrqf.f"
	i__1 = *m + 1;
#line 194 "dtzrqf.f"
	m1 = min(i__1,*n);
#line 195 "dtzrqf.f"
	for (k = *m; k >= 1; --k) {

/*           Use a Householder reflection to zero the kth row of A. */
/*           First set up the reflection. */

#line 200 "dtzrqf.f"
	    i__1 = *n - *m + 1;
#line 200 "dtzrqf.f"
	    dlarfg_(&i__1, &a[k + k * a_dim1], &a[k + m1 * a_dim1], lda, &tau[
		    k]);

#line 202 "dtzrqf.f"
	    if (tau[k] != 0. && k > 1) {

/*              We now perform the operation  A := A*P( k ). */

/*              Use the first ( k - 1 ) elements of TAU to store  a( k ), */
/*              where  a( k ) consists of the first ( k - 1 ) elements of */
/*              the  kth column  of  A.  Also  let  B  denote  the  first */
/*              ( k - 1 ) rows of the last ( n - m ) columns of A. */

#line 211 "dtzrqf.f"
		i__1 = k - 1;
#line 211 "dtzrqf.f"
		dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &tau[1], &c__1);

/*              Form   w = a( k ) + B*z( k )  in TAU. */

#line 215 "dtzrqf.f"
		i__1 = k - 1;
#line 215 "dtzrqf.f"
		i__2 = *n - *m;
#line 215 "dtzrqf.f"
		dgemv_("No transpose", &i__1, &i__2, &c_b8, &a[m1 * a_dim1 + 
			1], lda, &a[k + m1 * a_dim1], lda, &c_b8, &tau[1], &
			c__1, (ftnlen)12);

/*              Now form  a( k ) := a( k ) - tau*w */
/*              and       B      := B      - tau*w*z( k )**T. */

#line 221 "dtzrqf.f"
		i__1 = k - 1;
#line 221 "dtzrqf.f"
		d__1 = -tau[k];
#line 221 "dtzrqf.f"
		daxpy_(&i__1, &d__1, &tau[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 222 "dtzrqf.f"
		i__1 = k - 1;
#line 222 "dtzrqf.f"
		i__2 = *n - *m;
#line 222 "dtzrqf.f"
		d__1 = -tau[k];
#line 222 "dtzrqf.f"
		dger_(&i__1, &i__2, &d__1, &tau[1], &c__1, &a[k + m1 * a_dim1]
			, lda, &a[m1 * a_dim1 + 1], lda);
#line 224 "dtzrqf.f"
	    }
#line 225 "dtzrqf.f"
/* L20: */
#line 225 "dtzrqf.f"
	}
#line 226 "dtzrqf.f"
    }

#line 228 "dtzrqf.f"
    return 0;

/*     End of DTZRQF */

} /* dtzrqf_ */

