#line 1 "dlatrz.f"
/* dlatrz.f -- translated by f2c (version 20100827).
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

#line 1 "dlatrz.f"
/* > \brief \b DLATRZ factors an upper trapezoidal matrix by means of orthogonal transformations. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLATRZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatrz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatrz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatrz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLATRZ( M, N, L, A, LDA, TAU, WORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            L, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATRZ factors the M-by-(M+L) real upper trapezoidal matrix */
/* > [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means */
/* > of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal */
/* > matrix and, R and A1 are M-by-M upper triangular matrices. */
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
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of columns of the matrix A containing the */
/* >          meaningful part of the Householder vectors. N-M >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the leading M-by-N upper trapezoidal part of the */
/* >          array A must contain the matrix to be factorized. */
/* >          On exit, the leading M-by-M upper triangular part of A */
/* >          contains the upper triangular matrix R, and elements N-L+1 to */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (M) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

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
/* >                                                 (   0    ) */
/* >                                                 ( z( k ) ) */
/* > */
/* >  tau is a scalar and z( k ) is an l element vector. tau and z( k ) */
/* >  are chosen to annihilate the elements of the kth row of A2. */
/* > */
/* >  The scalar tau is returned in the kth element of TAU and the vector */
/* >  u( k ) in the kth row of A2, such that the elements of z( k ) are */
/* >  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in */
/* >  the upper triangular part of A1. */
/* > */
/* >  Z is given by */
/* > */
/* >     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlatrz_(integer *m, integer *n, integer *l, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int dlarz_(char *, integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), dlarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *);


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

/*     Quick return if possible */

#line 173 "dlatrz.f"
    /* Parameter adjustments */
#line 173 "dlatrz.f"
    a_dim1 = *lda;
#line 173 "dlatrz.f"
    a_offset = 1 + a_dim1;
#line 173 "dlatrz.f"
    a -= a_offset;
#line 173 "dlatrz.f"
    --tau;
#line 173 "dlatrz.f"
    --work;
#line 173 "dlatrz.f"

#line 173 "dlatrz.f"
    /* Function Body */
#line 173 "dlatrz.f"
    if (*m == 0) {
#line 174 "dlatrz.f"
	return 0;
#line 175 "dlatrz.f"
    } else if (*m == *n) {
#line 176 "dlatrz.f"
	i__1 = *n;
#line 176 "dlatrz.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 177 "dlatrz.f"
	    tau[i__] = 0.;
#line 178 "dlatrz.f"
/* L10: */
#line 178 "dlatrz.f"
	}
#line 179 "dlatrz.f"
	return 0;
#line 180 "dlatrz.f"
    }

#line 182 "dlatrz.f"
    for (i__ = *m; i__ >= 1; --i__) {

/*        Generate elementary reflector H(i) to annihilate */
/*        [ A(i,i) A(i,n-l+1:n) ] */

#line 187 "dlatrz.f"
	i__1 = *l + 1;
#line 187 "dlatrz.f"
	dlarfg_(&i__1, &a[i__ + i__ * a_dim1], &a[i__ + (*n - *l + 1) * 
		a_dim1], lda, &tau[i__]);

/*        Apply H(i) to A(1:i-1,i:n) from the right */

#line 191 "dlatrz.f"
	i__1 = i__ - 1;
#line 191 "dlatrz.f"
	i__2 = *n - i__ + 1;
#line 191 "dlatrz.f"
	dlarz_("Right", &i__1, &i__2, l, &a[i__ + (*n - *l + 1) * a_dim1], 
		lda, &tau[i__], &a[i__ * a_dim1 + 1], lda, &work[1], (ftnlen)
		5);

#line 194 "dlatrz.f"
/* L20: */
#line 194 "dlatrz.f"
    }

#line 196 "dlatrz.f"
    return 0;

/*     End of DLATRZ */

} /* dlatrz_ */

