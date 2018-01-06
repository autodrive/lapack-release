#line 1 "slarfy.f"
/* slarfy.f -- translated by f2c (version 20100827).
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

#line 1 "slarfy.f"
/* Table of constant values */

static doublereal c_b2 = 1.;
static doublereal c_b3 = 0.;
static integer c__1 = 1;

/* > \brief \b SLARFY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCV, LDC, N */
/*       REAL               TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFY applies an elementary reflector, or Householder matrix, H, */
/* > to an n x n symmetric matrix C, from both the left and the right. */
/* > */
/* > H is represented in the form */
/* > */
/* >    H = I - tau * v * v' */
/* > */
/* > where  tau  is a scalar and  v  is a vector. */
/* > */
/* > If  tau  is  zero, then  H  is taken to be the unit matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix C is stored. */
/* >          = 'U':  Upper triangle */
/* >          = 'L':  Lower triangle */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of rows and columns of the matrix C.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is REAL array, dimension */
/* >                  (1 + (N-1)*abs(INCV)) */
/* >          The vector v as described above. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* >          INCV is INTEGER */
/* >          The increment between successive elements of v.  INCV must */
/* >          not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL */
/* >          The value tau as described above. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC, N) */
/* >          On entry, the matrix C. */
/* >          On exit, C is overwritten by H * C * H'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C.  LDC >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup single_eig */

/*  ===================================================================== */
/* Subroutine */ int slarfy_(char *uplo, integer *n, doublereal *v, integer *
	incv, doublereal *tau, doublereal *c__, integer *ldc, doublereal *
	work, ftnlen uplo_len)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;

    /* Local variables */
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int ssyr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), ssymv_(char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);


/*  -- LAPACK test routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 143 "slarfy.f"
    /* Parameter adjustments */
#line 143 "slarfy.f"
    --v;
#line 143 "slarfy.f"
    c_dim1 = *ldc;
#line 143 "slarfy.f"
    c_offset = 1 + c_dim1;
#line 143 "slarfy.f"
    c__ -= c_offset;
#line 143 "slarfy.f"
    --work;
#line 143 "slarfy.f"

#line 143 "slarfy.f"
    /* Function Body */
#line 143 "slarfy.f"
    if (*tau == 0.) {
#line 143 "slarfy.f"
	return 0;
#line 143 "slarfy.f"
    }

/*     Form  w:= C * v */

#line 148 "slarfy.f"
    ssymv_(uplo, n, &c_b2, &c__[c_offset], ldc, &v[1], incv, &c_b3, &work[1], 
	    &c__1, (ftnlen)1);

#line 150 "slarfy.f"
    alpha = *tau * -.5 * sdot_(n, &work[1], &c__1, &v[1], incv);
#line 151 "slarfy.f"
    saxpy_(n, &alpha, &v[1], incv, &work[1], &c__1);

/*     C := C - v * w' - w * v' */

#line 155 "slarfy.f"
    d__1 = -(*tau);
#line 155 "slarfy.f"
    ssyr2_(uplo, n, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc, 
	    (ftnlen)1);

#line 157 "slarfy.f"
    return 0;

/*     End of SLARFY */

} /* slarfy_ */

