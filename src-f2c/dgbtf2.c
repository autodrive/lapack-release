#line 1 "dgbtf2.f"
/* dgbtf2.f -- translated by f2c (version 20100827).
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

#line 1 "dgbtf2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;

/* > \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of th
e algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBTF2 computes an LU factorization of a real m-by-n band matrix A */
/* > using partial pivoting with row interchanges. */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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

/* > \date December 2016 */

/* > \ingroup doubleGBcomputational */

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
/* Subroutine */ int dgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, km, jp, ju, kv;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), dswap_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     KV is the number of superdiagonals in the factor U, allowing for */
/*     fill-in. */

#line 185 "dgbtf2.f"
    /* Parameter adjustments */
#line 185 "dgbtf2.f"
    ab_dim1 = *ldab;
#line 185 "dgbtf2.f"
    ab_offset = 1 + ab_dim1;
#line 185 "dgbtf2.f"
    ab -= ab_offset;
#line 185 "dgbtf2.f"
    --ipiv;
#line 185 "dgbtf2.f"

#line 185 "dgbtf2.f"
    /* Function Body */
#line 185 "dgbtf2.f"
    kv = *ku + *kl;

/*     Test the input parameters. */

#line 189 "dgbtf2.f"
    *info = 0;
#line 190 "dgbtf2.f"
    if (*m < 0) {
#line 191 "dgbtf2.f"
	*info = -1;
#line 192 "dgbtf2.f"
    } else if (*n < 0) {
#line 193 "dgbtf2.f"
	*info = -2;
#line 194 "dgbtf2.f"
    } else if (*kl < 0) {
#line 195 "dgbtf2.f"
	*info = -3;
#line 196 "dgbtf2.f"
    } else if (*ku < 0) {
#line 197 "dgbtf2.f"
	*info = -4;
#line 198 "dgbtf2.f"
    } else if (*ldab < *kl + kv + 1) {
#line 199 "dgbtf2.f"
	*info = -6;
#line 200 "dgbtf2.f"
    }
#line 201 "dgbtf2.f"
    if (*info != 0) {
#line 202 "dgbtf2.f"
	i__1 = -(*info);
#line 202 "dgbtf2.f"
	xerbla_("DGBTF2", &i__1, (ftnlen)6);
#line 203 "dgbtf2.f"
	return 0;
#line 204 "dgbtf2.f"
    }

/*     Quick return if possible */

#line 208 "dgbtf2.f"
    if (*m == 0 || *n == 0) {
#line 208 "dgbtf2.f"
	return 0;
#line 208 "dgbtf2.f"
    }

/*     Gaussian elimination with partial pivoting */

/*     Set fill-in elements in columns KU+2 to KV to zero. */

#line 215 "dgbtf2.f"
    i__1 = min(kv,*n);
#line 215 "dgbtf2.f"
    for (j = *ku + 2; j <= i__1; ++j) {
#line 216 "dgbtf2.f"
	i__2 = *kl;
#line 216 "dgbtf2.f"
	for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
#line 217 "dgbtf2.f"
	    ab[i__ + j * ab_dim1] = 0.;
#line 218 "dgbtf2.f"
/* L10: */
#line 218 "dgbtf2.f"
	}
#line 219 "dgbtf2.f"
/* L20: */
#line 219 "dgbtf2.f"
    }

/*     JU is the index of the last column affected by the current stage */
/*     of the factorization. */

#line 224 "dgbtf2.f"
    ju = 1;

#line 226 "dgbtf2.f"
    i__1 = min(*m,*n);
#line 226 "dgbtf2.f"
    for (j = 1; j <= i__1; ++j) {

/*        Set fill-in elements in column J+KV to zero. */

#line 230 "dgbtf2.f"
	if (j + kv <= *n) {
#line 231 "dgbtf2.f"
	    i__2 = *kl;
#line 231 "dgbtf2.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 232 "dgbtf2.f"
		ab[i__ + (j + kv) * ab_dim1] = 0.;
#line 233 "dgbtf2.f"
/* L30: */
#line 233 "dgbtf2.f"
	    }
#line 234 "dgbtf2.f"
	}

/*        Find pivot and test for singularity. KM is the number of */
/*        subdiagonal elements in the current column. */

/* Computing MIN */
#line 239 "dgbtf2.f"
	i__2 = *kl, i__3 = *m - j;
#line 239 "dgbtf2.f"
	km = min(i__2,i__3);
#line 240 "dgbtf2.f"
	i__2 = km + 1;
#line 240 "dgbtf2.f"
	jp = idamax_(&i__2, &ab[kv + 1 + j * ab_dim1], &c__1);
#line 241 "dgbtf2.f"
	ipiv[j] = jp + j - 1;
#line 242 "dgbtf2.f"
	if (ab[kv + jp + j * ab_dim1] != 0.) {
/* Computing MAX */
/* Computing MIN */
#line 243 "dgbtf2.f"
	    i__4 = j + *ku + jp - 1;
#line 243 "dgbtf2.f"
	    i__2 = ju, i__3 = min(i__4,*n);
#line 243 "dgbtf2.f"
	    ju = max(i__2,i__3);

/*           Apply interchange to columns J to JU. */

#line 247 "dgbtf2.f"
	    if (jp != 1) {
#line 247 "dgbtf2.f"
		i__2 = ju - j + 1;
#line 247 "dgbtf2.f"
		i__3 = *ldab - 1;
#line 247 "dgbtf2.f"
		i__4 = *ldab - 1;
#line 247 "dgbtf2.f"
		dswap_(&i__2, &ab[kv + jp + j * ab_dim1], &i__3, &ab[kv + 1 + 
			j * ab_dim1], &i__4);
#line 247 "dgbtf2.f"
	    }

#line 251 "dgbtf2.f"
	    if (km > 0) {

/*              Compute multipliers. */

#line 255 "dgbtf2.f"
		d__1 = 1. / ab[kv + 1 + j * ab_dim1];
#line 255 "dgbtf2.f"
		dscal_(&km, &d__1, &ab[kv + 2 + j * ab_dim1], &c__1);

/*              Update trailing submatrix within the band. */

#line 259 "dgbtf2.f"
		if (ju > j) {
#line 259 "dgbtf2.f"
		    i__2 = ju - j;
#line 259 "dgbtf2.f"
		    i__3 = *ldab - 1;
#line 259 "dgbtf2.f"
		    i__4 = *ldab - 1;
#line 259 "dgbtf2.f"
		    dger_(&km, &i__2, &c_b9, &ab[kv + 2 + j * ab_dim1], &c__1,
			     &ab[kv + (j + 1) * ab_dim1], &i__3, &ab[kv + 1 + 
			    (j + 1) * ab_dim1], &i__4);
#line 259 "dgbtf2.f"
		}
#line 263 "dgbtf2.f"
	    }
#line 264 "dgbtf2.f"
	} else {

/*           If pivot is zero, set INFO to the index of the pivot */
/*           unless a zero pivot has already been found. */

#line 269 "dgbtf2.f"
	    if (*info == 0) {
#line 269 "dgbtf2.f"
		*info = j;
#line 269 "dgbtf2.f"
	    }
#line 271 "dgbtf2.f"
	}
#line 272 "dgbtf2.f"
/* L40: */
#line 272 "dgbtf2.f"
    }
#line 273 "dgbtf2.f"
    return 0;

/*     End of DGBTF2 */

} /* dgbtf2_ */

