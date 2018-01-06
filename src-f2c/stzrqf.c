#line 1 "stzrqf.f"
/* stzrqf.f -- translated by f2c (version 20100827).
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

#line 1 "stzrqf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = 1.;

/* > \brief \b STZRQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STZRQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stzrqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stzrqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stzrqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STZRQF( M, N, A, LDA, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine STZRZF. */
/* > */
/* > STZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          TAU is REAL array, dimension (M) */
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

/* > \ingroup realOTHERcomputational */

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
/* Subroutine */ int stzrqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, m1;
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), scopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), saxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *), xerbla_(
	    char *, integer *, ftnlen), slarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *);


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

#line 172 "stzrqf.f"
    /* Parameter adjustments */
#line 172 "stzrqf.f"
    a_dim1 = *lda;
#line 172 "stzrqf.f"
    a_offset = 1 + a_dim1;
#line 172 "stzrqf.f"
    a -= a_offset;
#line 172 "stzrqf.f"
    --tau;
#line 172 "stzrqf.f"

#line 172 "stzrqf.f"
    /* Function Body */
#line 172 "stzrqf.f"
    *info = 0;
#line 173 "stzrqf.f"
    if (*m < 0) {
#line 174 "stzrqf.f"
	*info = -1;
#line 175 "stzrqf.f"
    } else if (*n < *m) {
#line 176 "stzrqf.f"
	*info = -2;
#line 177 "stzrqf.f"
    } else if (*lda < max(1,*m)) {
#line 178 "stzrqf.f"
	*info = -4;
#line 179 "stzrqf.f"
    }
#line 180 "stzrqf.f"
    if (*info != 0) {
#line 181 "stzrqf.f"
	i__1 = -(*info);
#line 181 "stzrqf.f"
	xerbla_("STZRQF", &i__1, (ftnlen)6);
#line 182 "stzrqf.f"
	return 0;
#line 183 "stzrqf.f"
    }

/*     Perform the factorization. */

#line 187 "stzrqf.f"
    if (*m == 0) {
#line 187 "stzrqf.f"
	return 0;
#line 187 "stzrqf.f"
    }
#line 189 "stzrqf.f"
    if (*m == *n) {
#line 190 "stzrqf.f"
	i__1 = *n;
#line 190 "stzrqf.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "stzrqf.f"
	    tau[i__] = 0.;
#line 192 "stzrqf.f"
/* L10: */
#line 192 "stzrqf.f"
	}
#line 193 "stzrqf.f"
    } else {
/* Computing MIN */
#line 194 "stzrqf.f"
	i__1 = *m + 1;
#line 194 "stzrqf.f"
	m1 = min(i__1,*n);
#line 195 "stzrqf.f"
	for (k = *m; k >= 1; --k) {

/*           Use a Householder reflection to zero the kth row of A. */
/*           First set up the reflection. */

#line 200 "stzrqf.f"
	    i__1 = *n - *m + 1;
#line 200 "stzrqf.f"
	    slarfg_(&i__1, &a[k + k * a_dim1], &a[k + m1 * a_dim1], lda, &tau[
		    k]);

#line 202 "stzrqf.f"
	    if (tau[k] != 0. && k > 1) {

/*              We now perform the operation  A := A*P( k ). */

/*              Use the first ( k - 1 ) elements of TAU to store  a( k ), */
/*              where  a( k ) consists of the first ( k - 1 ) elements of */
/*              the  kth column  of  A.  Also  let  B  denote  the  first */
/*              ( k - 1 ) rows of the last ( n - m ) columns of A. */

#line 211 "stzrqf.f"
		i__1 = k - 1;
#line 211 "stzrqf.f"
		scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &tau[1], &c__1);

/*              Form   w = a( k ) + B*z( k )  in TAU. */

#line 215 "stzrqf.f"
		i__1 = k - 1;
#line 215 "stzrqf.f"
		i__2 = *n - *m;
#line 215 "stzrqf.f"
		sgemv_("No transpose", &i__1, &i__2, &c_b8, &a[m1 * a_dim1 + 
			1], lda, &a[k + m1 * a_dim1], lda, &c_b8, &tau[1], &
			c__1, (ftnlen)12);

/*              Now form  a( k ) := a( k ) - tau*w */
/*              and       B      := B      - tau*w*z( k )**T. */

#line 221 "stzrqf.f"
		i__1 = k - 1;
#line 221 "stzrqf.f"
		d__1 = -tau[k];
#line 221 "stzrqf.f"
		saxpy_(&i__1, &d__1, &tau[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 222 "stzrqf.f"
		i__1 = k - 1;
#line 222 "stzrqf.f"
		i__2 = *n - *m;
#line 222 "stzrqf.f"
		d__1 = -tau[k];
#line 222 "stzrqf.f"
		sger_(&i__1, &i__2, &d__1, &tau[1], &c__1, &a[k + m1 * a_dim1]
			, lda, &a[m1 * a_dim1 + 1], lda);
#line 224 "stzrqf.f"
	    }
#line 225 "stzrqf.f"
/* L20: */
#line 225 "stzrqf.f"
	}
#line 226 "stzrqf.f"
    }

#line 228 "stzrqf.f"
    return 0;

/*     End of STZRQF */

} /* stzrqf_ */

