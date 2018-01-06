#line 1 "clarfy.f"
/* clarfy.f -- translated by f2c (version 20100827).
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

#line 1 "clarfy.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CLARFY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCV, LDC, N */
/*       COMPLEX            TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFY applies an elementary reflector, or Householder matrix, H, */
/* > to an n x n Hermitian matrix C, from both the left and the right. */
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
/* >          Hermitian matrix C is stored. */
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
/* >          V is COMPLEX array, dimension */
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
/* >          TAU is COMPLEX */
/* >          The value tau as described above. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC, N) */
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
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex_eig */

/*  ===================================================================== */
/* Subroutine */ int clarfy_(char *uplo, integer *n, doublecomplex *v, 
	integer *incv, doublecomplex *tau, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    extern /* Subroutine */ int cher2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen);
    static doublecomplex alpha;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int chemv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), caxpy_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);


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

#line 145 "clarfy.f"
    /* Parameter adjustments */
#line 145 "clarfy.f"
    --v;
#line 145 "clarfy.f"
    c_dim1 = *ldc;
#line 145 "clarfy.f"
    c_offset = 1 + c_dim1;
#line 145 "clarfy.f"
    c__ -= c_offset;
#line 145 "clarfy.f"
    --work;
#line 145 "clarfy.f"

#line 145 "clarfy.f"
    /* Function Body */
#line 145 "clarfy.f"
    if (tau->r == 0. && tau->i == 0.) {
#line 145 "clarfy.f"
	return 0;
#line 145 "clarfy.f"
    }

/*     Form  w:= C * v */

#line 150 "clarfy.f"
    chemv_(uplo, n, &c_b1, &c__[c_offset], ldc, &v[1], incv, &c_b2, &work[1], 
	    &c__1, (ftnlen)1);

#line 152 "clarfy.f"
    z__3.r = -.5, z__3.i = -0.;
#line 152 "clarfy.f"
    z__2.r = z__3.r * tau->r - z__3.i * tau->i, z__2.i = z__3.r * tau->i + 
	    z__3.i * tau->r;
#line 152 "clarfy.f"
    cdotc_(&z__4, n, &work[1], &c__1, &v[1], incv);
#line 152 "clarfy.f"
    z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * z__4.i + 
	    z__2.i * z__4.r;
#line 152 "clarfy.f"
    alpha.r = z__1.r, alpha.i = z__1.i;
#line 153 "clarfy.f"
    caxpy_(n, &alpha, &v[1], incv, &work[1], &c__1);

/*     C := C - v * w' - w * v' */

#line 157 "clarfy.f"
    z__1.r = -tau->r, z__1.i = -tau->i;
#line 157 "clarfy.f"
    cher2_(uplo, n, &z__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc, 
	    (ftnlen)1);

#line 159 "clarfy.f"
    return 0;

/*     End of CLARFY */

} /* clarfy_ */

