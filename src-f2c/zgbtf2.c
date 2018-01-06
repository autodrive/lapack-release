#line 1 "zgbtf2.f"
/* zgbtf2.f -- translated by f2c (version 20100827).
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

#line 1 "zgbtf2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZGBTF2 computes the LU factorization of a general band matrix using the unblocked version of th
e algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbtf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbtf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbtf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBTF2 computes an LU factorization of a complex m-by-n band matrix */
/* > A using partial pivoting with row interchanges. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On entry, the matrix A in band storage, in rows KL+1 to */
/* >          2*KL+KU+1; rows 1 to KL of the array need not be set. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl) */
/* > */
/* >          On exit, details of the factorization: U is stored as an */
/* >          upper triangular band matrix with KL+KU superdiagonals in */
/* >          rows 1 to KL+KU+1, and the multipliers used during the */
/* >          factorization are stored in rows KL+KU+2 to 2*KL+KU+1. */
/* >          See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (min(M,N)) */
/* >          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization */
/* >               has been completed, but the factor U is exactly */
/* >               singular, and division by zero will occur if it is used */
/* >               to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16GBcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  M = N = 6, KL = 2, KU = 1: */
/* > */
/* >  On entry:                       On exit: */
/* > */
/* >      *    *    *    +    +    +       *    *    *   u14  u25  u36 */
/* >      *    *    +    +    +    +       *    *   u13  u24  u35  u46 */
/* >      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 */
/* >     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 */
/* >     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   * */
/* >     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    * */
/* > */
/* >  Array elements marked * are not used by the routine; elements marked */
/* >  + need not be set on entry, but are required by the routine to store */
/* >  elements of U, because of fill-in resulting from the row */
/* >  interchanges. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, km, jp, ju, kv;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgeru_(integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zswap_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    extern integer izamax_(integer *, doublecomplex *, integer *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     KV is the number of superdiagonals in the factor U, allowing for */
/*     fill-in. */

#line 186 "zgbtf2.f"
    /* Parameter adjustments */
#line 186 "zgbtf2.f"
    ab_dim1 = *ldab;
#line 186 "zgbtf2.f"
    ab_offset = 1 + ab_dim1;
#line 186 "zgbtf2.f"
    ab -= ab_offset;
#line 186 "zgbtf2.f"
    --ipiv;
#line 186 "zgbtf2.f"

#line 186 "zgbtf2.f"
    /* Function Body */
#line 186 "zgbtf2.f"
    kv = *ku + *kl;

/*     Test the input parameters. */

#line 190 "zgbtf2.f"
    *info = 0;
#line 191 "zgbtf2.f"
    if (*m < 0) {
#line 192 "zgbtf2.f"
	*info = -1;
#line 193 "zgbtf2.f"
    } else if (*n < 0) {
#line 194 "zgbtf2.f"
	*info = -2;
#line 195 "zgbtf2.f"
    } else if (*kl < 0) {
#line 196 "zgbtf2.f"
	*info = -3;
#line 197 "zgbtf2.f"
    } else if (*ku < 0) {
#line 198 "zgbtf2.f"
	*info = -4;
#line 199 "zgbtf2.f"
    } else if (*ldab < *kl + kv + 1) {
#line 200 "zgbtf2.f"
	*info = -6;
#line 201 "zgbtf2.f"
    }
#line 202 "zgbtf2.f"
    if (*info != 0) {
#line 203 "zgbtf2.f"
	i__1 = -(*info);
#line 203 "zgbtf2.f"
	xerbla_("ZGBTF2", &i__1, (ftnlen)6);
#line 204 "zgbtf2.f"
	return 0;
#line 205 "zgbtf2.f"
    }

/*     Quick return if possible */

#line 209 "zgbtf2.f"
    if (*m == 0 || *n == 0) {
#line 209 "zgbtf2.f"
	return 0;
#line 209 "zgbtf2.f"
    }

/*     Gaussian elimination with partial pivoting */

/*     Set fill-in elements in columns KU+2 to KV to zero. */

#line 216 "zgbtf2.f"
    i__1 = min(kv,*n);
#line 216 "zgbtf2.f"
    for (j = *ku + 2; j <= i__1; ++j) {
#line 217 "zgbtf2.f"
	i__2 = *kl;
#line 217 "zgbtf2.f"
	for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
#line 218 "zgbtf2.f"
	    i__3 = i__ + j * ab_dim1;
#line 218 "zgbtf2.f"
	    ab[i__3].r = 0., ab[i__3].i = 0.;
#line 219 "zgbtf2.f"
/* L10: */
#line 219 "zgbtf2.f"
	}
#line 220 "zgbtf2.f"
/* L20: */
#line 220 "zgbtf2.f"
    }

/*     JU is the index of the last column affected by the current stage */
/*     of the factorization. */

#line 225 "zgbtf2.f"
    ju = 1;

#line 227 "zgbtf2.f"
    i__1 = min(*m,*n);
#line 227 "zgbtf2.f"
    for (j = 1; j <= i__1; ++j) {

/*        Set fill-in elements in column J+KV to zero. */

#line 231 "zgbtf2.f"
	if (j + kv <= *n) {
#line 232 "zgbtf2.f"
	    i__2 = *kl;
#line 232 "zgbtf2.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 233 "zgbtf2.f"
		i__3 = i__ + (j + kv) * ab_dim1;
#line 233 "zgbtf2.f"
		ab[i__3].r = 0., ab[i__3].i = 0.;
#line 234 "zgbtf2.f"
/* L30: */
#line 234 "zgbtf2.f"
	    }
#line 235 "zgbtf2.f"
	}

/*        Find pivot and test for singularity. KM is the number of */
/*        subdiagonal elements in the current column. */

/* Computing MIN */
#line 240 "zgbtf2.f"
	i__2 = *kl, i__3 = *m - j;
#line 240 "zgbtf2.f"
	km = min(i__2,i__3);
#line 241 "zgbtf2.f"
	i__2 = km + 1;
#line 241 "zgbtf2.f"
	jp = izamax_(&i__2, &ab[kv + 1 + j * ab_dim1], &c__1);
#line 242 "zgbtf2.f"
	ipiv[j] = jp + j - 1;
#line 243 "zgbtf2.f"
	i__2 = kv + jp + j * ab_dim1;
#line 243 "zgbtf2.f"
	if (ab[i__2].r != 0. || ab[i__2].i != 0.) {
/* Computing MAX */
/* Computing MIN */
#line 244 "zgbtf2.f"
	    i__4 = j + *ku + jp - 1;
#line 244 "zgbtf2.f"
	    i__2 = ju, i__3 = min(i__4,*n);
#line 244 "zgbtf2.f"
	    ju = max(i__2,i__3);

/*           Apply interchange to columns J to JU. */

#line 248 "zgbtf2.f"
	    if (jp != 1) {
#line 248 "zgbtf2.f"
		i__2 = ju - j + 1;
#line 248 "zgbtf2.f"
		i__3 = *ldab - 1;
#line 248 "zgbtf2.f"
		i__4 = *ldab - 1;
#line 248 "zgbtf2.f"
		zswap_(&i__2, &ab[kv + jp + j * ab_dim1], &i__3, &ab[kv + 1 + 
			j * ab_dim1], &i__4);
#line 248 "zgbtf2.f"
	    }
#line 251 "zgbtf2.f"
	    if (km > 0) {

/*              Compute multipliers. */

#line 255 "zgbtf2.f"
		z_div(&z__1, &c_b1, &ab[kv + 1 + j * ab_dim1]);
#line 255 "zgbtf2.f"
		zscal_(&km, &z__1, &ab[kv + 2 + j * ab_dim1], &c__1);

/*              Update trailing submatrix within the band. */

#line 259 "zgbtf2.f"
		if (ju > j) {
#line 259 "zgbtf2.f"
		    i__2 = ju - j;
#line 259 "zgbtf2.f"
		    z__1.r = -1., z__1.i = -0.;
#line 259 "zgbtf2.f"
		    i__3 = *ldab - 1;
#line 259 "zgbtf2.f"
		    i__4 = *ldab - 1;
#line 259 "zgbtf2.f"
		    zgeru_(&km, &i__2, &z__1, &ab[kv + 2 + j * ab_dim1], &
			    c__1, &ab[kv + (j + 1) * ab_dim1], &i__3, &ab[kv 
			    + 1 + (j + 1) * ab_dim1], &i__4);
#line 259 "zgbtf2.f"
		}
#line 263 "zgbtf2.f"
	    }
#line 264 "zgbtf2.f"
	} else {

/*           If pivot is zero, set INFO to the index of the pivot */
/*           unless a zero pivot has already been found. */

#line 269 "zgbtf2.f"
	    if (*info == 0) {
#line 269 "zgbtf2.f"
		*info = j;
#line 269 "zgbtf2.f"
	    }
#line 271 "zgbtf2.f"
	}
#line 272 "zgbtf2.f"
/* L40: */
#line 272 "zgbtf2.f"
    }
#line 273 "zgbtf2.f"
    return 0;

/*     End of ZGBTF2 */

} /* zgbtf2_ */

