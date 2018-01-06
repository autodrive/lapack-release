#line 1 "slapll.f"
/* slapll.f -- translated by f2c (version 20100827).
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

#line 1 "slapll.f"
/* > \brief \b SLAPLL measures the linear dependence of two vectors. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAPLL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapll.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapll.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapll.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAPLL( N, X, INCX, Y, INCY, SSMIN ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, INCY, N */
/*       REAL               SSMIN */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given two column vectors X and Y, let */
/* > */
/* >                      A = ( X Y ). */
/* > */
/* > The subroutine first computes the QR factorization of A = Q*R, */
/* > and then computes the SVD of the 2-by-2 upper triangular matrix R. */
/* > The smaller singular value of R is returned in SSMIN, which is used */
/* > as the measurement of the linear dependency of the vectors X and Y. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The length of the vectors X and Y. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, */
/* >                         dimension (1+(N-1)*INCX) */
/* >          On entry, X contains the N-vector X. */
/* >          On exit, X is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array, */
/* >                         dimension (1+(N-1)*INCY) */
/* >          On entry, Y contains the N-vector Y. */
/* >          On exit, Y is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >          The increment between successive elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* >          SSMIN is REAL */
/* >          The smallest singular value of the N-by-2 matrix A = ( X Y ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slapll_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *ssmin)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal c__, a11, a12, a22, tau;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int slas2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal ssmax;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), slarfg_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 138 "slapll.f"
    /* Parameter adjustments */
#line 138 "slapll.f"
    --y;
#line 138 "slapll.f"
    --x;
#line 138 "slapll.f"

#line 138 "slapll.f"
    /* Function Body */
#line 138 "slapll.f"
    if (*n <= 1) {
#line 139 "slapll.f"
	*ssmin = 0.;
#line 140 "slapll.f"
	return 0;
#line 141 "slapll.f"
    }

/*     Compute the QR factorization of the N-by-2 matrix ( X Y ) */

#line 145 "slapll.f"
    slarfg_(n, &x[1], &x[*incx + 1], incx, &tau);
#line 146 "slapll.f"
    a11 = x[1];
#line 147 "slapll.f"
    x[1] = 1.;

#line 149 "slapll.f"
    c__ = -tau * sdot_(n, &x[1], incx, &y[1], incy);
#line 150 "slapll.f"
    saxpy_(n, &c__, &x[1], incx, &y[1], incy);

#line 152 "slapll.f"
    i__1 = *n - 1;
#line 152 "slapll.f"
    slarfg_(&i__1, &y[*incy + 1], &y[(*incy << 1) + 1], incy, &tau);

#line 154 "slapll.f"
    a12 = y[1];
#line 155 "slapll.f"
    a22 = y[*incy + 1];

/*     Compute the SVD of 2-by-2 Upper triangular matrix. */

#line 159 "slapll.f"
    slas2_(&a11, &a12, &a22, ssmin, &ssmax);

#line 161 "slapll.f"
    return 0;

/*     End of SLAPLL */

} /* slapll_ */

