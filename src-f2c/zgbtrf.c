#line 1 "zgbtrf.f"
/* zgbtrf.f -- translated by f2c (version 20100827).
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

#line 1 "zgbtrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c__65 = 65;

/* > \brief \b ZGBTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbtrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbtrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbtrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */

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
/* > ZGBTRF computes an LU factorization of a complex m-by-n band matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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

/* > \date November 2011 */

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
/* >  elements of U because of fill-in resulting from the row interchanges. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, i2, i3, j2, j3, k2, jb, nb, ii, jj, jm, ip, jp, km,
	     ju, kv, nw;
    static doublecomplex temp;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static doublecomplex work13[4160]	/* was [65][64] */, work31[4160]	
	    /* was [65][64] */;
    extern /* Subroutine */ int zgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zcopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zswap_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztrsm_(
	    char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), zgbtf2_(integer *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen), izamax_(integer *, 
	    doublecomplex *, integer *);
    extern /* Subroutine */ int zlaswp_(integer *, doublecomplex *, integer *,
	     integer *, integer *, integer *, integer *);


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     KV is the number of superdiagonals in the factor U, allowing for */
/*     fill-in */

#line 194 "zgbtrf.f"
    /* Parameter adjustments */
#line 194 "zgbtrf.f"
    ab_dim1 = *ldab;
#line 194 "zgbtrf.f"
    ab_offset = 1 + ab_dim1;
#line 194 "zgbtrf.f"
    ab -= ab_offset;
#line 194 "zgbtrf.f"
    --ipiv;
#line 194 "zgbtrf.f"

#line 194 "zgbtrf.f"
    /* Function Body */
#line 194 "zgbtrf.f"
    kv = *ku + *kl;

/*     Test the input parameters. */

#line 198 "zgbtrf.f"
    *info = 0;
#line 199 "zgbtrf.f"
    if (*m < 0) {
#line 200 "zgbtrf.f"
	*info = -1;
#line 201 "zgbtrf.f"
    } else if (*n < 0) {
#line 202 "zgbtrf.f"
	*info = -2;
#line 203 "zgbtrf.f"
    } else if (*kl < 0) {
#line 204 "zgbtrf.f"
	*info = -3;
#line 205 "zgbtrf.f"
    } else if (*ku < 0) {
#line 206 "zgbtrf.f"
	*info = -4;
#line 207 "zgbtrf.f"
    } else if (*ldab < *kl + kv + 1) {
#line 208 "zgbtrf.f"
	*info = -6;
#line 209 "zgbtrf.f"
    }
#line 210 "zgbtrf.f"
    if (*info != 0) {
#line 211 "zgbtrf.f"
	i__1 = -(*info);
#line 211 "zgbtrf.f"
	xerbla_("ZGBTRF", &i__1, (ftnlen)6);
#line 212 "zgbtrf.f"
	return 0;
#line 213 "zgbtrf.f"
    }

/*     Quick return if possible */

#line 217 "zgbtrf.f"
    if (*m == 0 || *n == 0) {
#line 217 "zgbtrf.f"
	return 0;
#line 217 "zgbtrf.f"
    }

/*     Determine the block size for this environment */

#line 222 "zgbtrf.f"
    nb = ilaenv_(&c__1, "ZGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);

/*     The block size must not exceed the limit set by the size of the */
/*     local arrays WORK13 and WORK31. */

#line 227 "zgbtrf.f"
    nb = min(nb,64);

#line 229 "zgbtrf.f"
    if (nb <= 1 || nb > *kl) {

/*        Use unblocked code */

#line 233 "zgbtrf.f"
	zgbtf2_(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
#line 234 "zgbtrf.f"
    } else {

/*        Use blocked code */

/*        Zero the superdiagonal elements of the work array WORK13 */

#line 240 "zgbtrf.f"
	i__1 = nb;
#line 240 "zgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 241 "zgbtrf.f"
	    i__2 = j - 1;
#line 241 "zgbtrf.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 242 "zgbtrf.f"
		i__3 = i__ + j * 65 - 66;
#line 242 "zgbtrf.f"
		work13[i__3].r = 0., work13[i__3].i = 0.;
#line 243 "zgbtrf.f"
/* L10: */
#line 243 "zgbtrf.f"
	    }
#line 244 "zgbtrf.f"
/* L20: */
#line 244 "zgbtrf.f"
	}

/*        Zero the subdiagonal elements of the work array WORK31 */

#line 248 "zgbtrf.f"
	i__1 = nb;
#line 248 "zgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 249 "zgbtrf.f"
	    i__2 = nb;
#line 249 "zgbtrf.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 250 "zgbtrf.f"
		i__3 = i__ + j * 65 - 66;
#line 250 "zgbtrf.f"
		work31[i__3].r = 0., work31[i__3].i = 0.;
#line 251 "zgbtrf.f"
/* L30: */
#line 251 "zgbtrf.f"
	    }
#line 252 "zgbtrf.f"
/* L40: */
#line 252 "zgbtrf.f"
	}

/*        Gaussian elimination with partial pivoting */

/*        Set fill-in elements in columns KU+2 to KV to zero */

#line 258 "zgbtrf.f"
	i__1 = min(kv,*n);
#line 258 "zgbtrf.f"
	for (j = *ku + 2; j <= i__1; ++j) {
#line 259 "zgbtrf.f"
	    i__2 = *kl;
#line 259 "zgbtrf.f"
	    for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
#line 260 "zgbtrf.f"
		i__3 = i__ + j * ab_dim1;
#line 260 "zgbtrf.f"
		ab[i__3].r = 0., ab[i__3].i = 0.;
#line 261 "zgbtrf.f"
/* L50: */
#line 261 "zgbtrf.f"
	    }
#line 262 "zgbtrf.f"
/* L60: */
#line 262 "zgbtrf.f"
	}

/*        JU is the index of the last column affected by the current */
/*        stage of the factorization */

#line 267 "zgbtrf.f"
	ju = 1;

#line 269 "zgbtrf.f"
	i__1 = min(*m,*n);
#line 269 "zgbtrf.f"
	i__2 = nb;
#line 269 "zgbtrf.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 270 "zgbtrf.f"
	    i__3 = nb, i__4 = min(*m,*n) - j + 1;
#line 270 "zgbtrf.f"
	    jb = min(i__3,i__4);

/*           The active part of the matrix is partitioned */

/*              A11   A12   A13 */
/*              A21   A22   A23 */
/*              A31   A32   A33 */

/*           Here A11, A21 and A31 denote the current block of JB columns */
/*           which is about to be factorized. The number of rows in the */
/*           partitioning are JB, I2, I3 respectively, and the numbers */
/*           of columns are JB, J2, J3. The superdiagonal elements of A13 */
/*           and the subdiagonal elements of A31 lie outside the band. */

/* Computing MIN */
#line 284 "zgbtrf.f"
	    i__3 = *kl - jb, i__4 = *m - j - jb + 1;
#line 284 "zgbtrf.f"
	    i2 = min(i__3,i__4);
/* Computing MIN */
#line 285 "zgbtrf.f"
	    i__3 = jb, i__4 = *m - j - *kl + 1;
#line 285 "zgbtrf.f"
	    i3 = min(i__3,i__4);

/*           J2 and J3 are computed after JU has been updated. */

/*           Factorize the current block of JB columns */

#line 291 "zgbtrf.f"
	    i__3 = j + jb - 1;
#line 291 "zgbtrf.f"
	    for (jj = j; jj <= i__3; ++jj) {

/*              Set fill-in elements in column JJ+KV to zero */

#line 295 "zgbtrf.f"
		if (jj + kv <= *n) {
#line 296 "zgbtrf.f"
		    i__4 = *kl;
#line 296 "zgbtrf.f"
		    for (i__ = 1; i__ <= i__4; ++i__) {
#line 297 "zgbtrf.f"
			i__5 = i__ + (jj + kv) * ab_dim1;
#line 297 "zgbtrf.f"
			ab[i__5].r = 0., ab[i__5].i = 0.;
#line 298 "zgbtrf.f"
/* L70: */
#line 298 "zgbtrf.f"
		    }
#line 299 "zgbtrf.f"
		}

/*              Find pivot and test for singularity. KM is the number of */
/*              subdiagonal elements in the current column. */

/* Computing MIN */
#line 304 "zgbtrf.f"
		i__4 = *kl, i__5 = *m - jj;
#line 304 "zgbtrf.f"
		km = min(i__4,i__5);
#line 305 "zgbtrf.f"
		i__4 = km + 1;
#line 305 "zgbtrf.f"
		jp = izamax_(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
#line 306 "zgbtrf.f"
		ipiv[jj] = jp + jj - j;
#line 307 "zgbtrf.f"
		i__4 = kv + jp + jj * ab_dim1;
#line 307 "zgbtrf.f"
		if (ab[i__4].r != 0. || ab[i__4].i != 0.) {
/* Computing MAX */
/* Computing MIN */
#line 308 "zgbtrf.f"
		    i__6 = jj + *ku + jp - 1;
#line 308 "zgbtrf.f"
		    i__4 = ju, i__5 = min(i__6,*n);
#line 308 "zgbtrf.f"
		    ju = max(i__4,i__5);
#line 309 "zgbtrf.f"
		    if (jp != 1) {

/*                    Apply interchange to columns J to J+JB-1 */

#line 313 "zgbtrf.f"
			if (jp + jj - 1 < j + *kl) {

#line 315 "zgbtrf.f"
			    i__4 = *ldab - 1;
#line 315 "zgbtrf.f"
			    i__5 = *ldab - 1;
#line 315 "zgbtrf.f"
			    zswap_(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
				    i__4, &ab[kv + jp + jj - j + j * ab_dim1],
				     &i__5);
#line 317 "zgbtrf.f"
			} else {

/*                       The interchange affects columns J to JJ-1 of A31 */
/*                       which are stored in the work array WORK31 */

#line 322 "zgbtrf.f"
			    i__4 = jj - j;
#line 322 "zgbtrf.f"
			    i__5 = *ldab - 1;
#line 322 "zgbtrf.f"
			    zswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], 
				    &i__5, &work31[jp + jj - j - *kl - 1], &
				    c__65);
#line 324 "zgbtrf.f"
			    i__4 = j + jb - jj;
#line 324 "zgbtrf.f"
			    i__5 = *ldab - 1;
#line 324 "zgbtrf.f"
			    i__6 = *ldab - 1;
#line 324 "zgbtrf.f"
			    zswap_(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
				    ab[kv + jp + jj * ab_dim1], &i__6);
#line 326 "zgbtrf.f"
			}
#line 327 "zgbtrf.f"
		    }

/*                 Compute multipliers */

#line 331 "zgbtrf.f"
		    z_div(&z__1, &c_b1, &ab[kv + 1 + jj * ab_dim1]);
#line 331 "zgbtrf.f"
		    zscal_(&km, &z__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

/*                 Update trailing submatrix within the band and within */
/*                 the current block. JM is the index of the last column */
/*                 which needs to be updated. */

/* Computing MIN */
#line 338 "zgbtrf.f"
		    i__4 = ju, i__5 = j + jb - 1;
#line 338 "zgbtrf.f"
		    jm = min(i__4,i__5);
#line 339 "zgbtrf.f"
		    if (jm > jj) {
#line 339 "zgbtrf.f"
			i__4 = jm - jj;
#line 339 "zgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 339 "zgbtrf.f"
			i__5 = *ldab - 1;
#line 339 "zgbtrf.f"
			i__6 = *ldab - 1;
#line 339 "zgbtrf.f"
			zgeru_(&km, &i__4, &z__1, &ab[kv + 2 + jj * ab_dim1], 
				&c__1, &ab[kv + (jj + 1) * ab_dim1], &i__5, &
				ab[kv + 1 + (jj + 1) * ab_dim1], &i__6);
#line 339 "zgbtrf.f"
		    }
#line 343 "zgbtrf.f"
		} else {

/*                 If pivot is zero, set INFO to the index of the pivot */
/*                 unless a zero pivot has already been found. */

#line 348 "zgbtrf.f"
		    if (*info == 0) {
#line 348 "zgbtrf.f"
			*info = jj;
#line 348 "zgbtrf.f"
		    }
#line 350 "zgbtrf.f"
		}

/*              Copy current column of A31 into the work array WORK31 */

/* Computing MIN */
#line 354 "zgbtrf.f"
		i__4 = jj - j + 1;
#line 354 "zgbtrf.f"
		nw = min(i__4,i3);
#line 355 "zgbtrf.f"
		if (nw > 0) {
#line 355 "zgbtrf.f"
		    zcopy_(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
			    c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
#line 355 "zgbtrf.f"
		}
#line 358 "zgbtrf.f"
/* L80: */
#line 358 "zgbtrf.f"
	    }
#line 359 "zgbtrf.f"
	    if (j + jb <= *n) {

/*              Apply the row interchanges to the other blocks. */

/* Computing MIN */
#line 363 "zgbtrf.f"
		i__3 = ju - j + 1;
#line 363 "zgbtrf.f"
		j2 = min(i__3,kv) - jb;
/* Computing MAX */
#line 364 "zgbtrf.f"
		i__3 = 0, i__4 = ju - j - kv + 1;
#line 364 "zgbtrf.f"
		j3 = max(i__3,i__4);

/*              Use ZLASWP to apply the row interchanges to A12, A22, and */
/*              A32. */

#line 369 "zgbtrf.f"
		i__3 = *ldab - 1;
#line 369 "zgbtrf.f"
		zlaswp_(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
			c__1, &jb, &ipiv[j], &c__1);

/*              Adjust the pivot indices. */

#line 374 "zgbtrf.f"
		i__3 = j + jb - 1;
#line 374 "zgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 375 "zgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 376 "zgbtrf.f"
/* L90: */
#line 376 "zgbtrf.f"
		}

/*              Apply the row interchanges to A13, A23, and A33 */
/*              columnwise. */

#line 381 "zgbtrf.f"
		k2 = j - 1 + jb + j2;
#line 382 "zgbtrf.f"
		i__3 = j3;
#line 382 "zgbtrf.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 383 "zgbtrf.f"
		    jj = k2 + i__;
#line 384 "zgbtrf.f"
		    i__4 = j + jb - 1;
#line 384 "zgbtrf.f"
		    for (ii = j + i__ - 1; ii <= i__4; ++ii) {
#line 385 "zgbtrf.f"
			ip = ipiv[ii];
#line 386 "zgbtrf.f"
			if (ip != ii) {
#line 387 "zgbtrf.f"
			    i__5 = kv + 1 + ii - jj + jj * ab_dim1;
#line 387 "zgbtrf.f"
			    temp.r = ab[i__5].r, temp.i = ab[i__5].i;
#line 388 "zgbtrf.f"
			    i__5 = kv + 1 + ii - jj + jj * ab_dim1;
#line 388 "zgbtrf.f"
			    i__6 = kv + 1 + ip - jj + jj * ab_dim1;
#line 388 "zgbtrf.f"
			    ab[i__5].r = ab[i__6].r, ab[i__5].i = ab[i__6].i;
#line 389 "zgbtrf.f"
			    i__5 = kv + 1 + ip - jj + jj * ab_dim1;
#line 389 "zgbtrf.f"
			    ab[i__5].r = temp.r, ab[i__5].i = temp.i;
#line 390 "zgbtrf.f"
			}
#line 391 "zgbtrf.f"
/* L100: */
#line 391 "zgbtrf.f"
		    }
#line 392 "zgbtrf.f"
/* L110: */
#line 392 "zgbtrf.f"
		}

/*              Update the relevant part of the trailing submatrix */

#line 396 "zgbtrf.f"
		if (j2 > 0) {

/*                 Update A12 */

#line 400 "zgbtrf.f"
		    i__3 = *ldab - 1;
#line 400 "zgbtrf.f"
		    i__4 = *ldab - 1;
#line 400 "zgbtrf.f"
		    ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
			    &c_b1, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv + 
			    1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4, (
			    ftnlen)5, (ftnlen)12, (ftnlen)4);

#line 404 "zgbtrf.f"
		    if (i2 > 0) {

/*                    Update A22 */

#line 408 "zgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 408 "zgbtrf.f"
			i__3 = *ldab - 1;
#line 408 "zgbtrf.f"
			i__4 = *ldab - 1;
#line 408 "zgbtrf.f"
			i__5 = *ldab - 1;
#line 408 "zgbtrf.f"
			zgemm_("No transpose", "No transpose", &i2, &j2, &jb, 
				&z__1, &ab[kv + 1 + jb + j * ab_dim1], &i__3, 
				&ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4, 
				&c_b1, &ab[kv + 1 + (j + jb) * ab_dim1], &
				i__5, (ftnlen)12, (ftnlen)12);
#line 412 "zgbtrf.f"
		    }

#line 414 "zgbtrf.f"
		    if (i3 > 0) {

/*                    Update A32 */

#line 418 "zgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 418 "zgbtrf.f"
			i__3 = *ldab - 1;
#line 418 "zgbtrf.f"
			i__4 = *ldab - 1;
#line 418 "zgbtrf.f"
			zgemm_("No transpose", "No transpose", &i3, &j2, &jb, 
				&z__1, work31, &c__65, &ab[kv + 1 - jb + (j + 
				jb) * ab_dim1], &i__3, &c_b1, &ab[kv + *kl + 
				1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)
				12, (ftnlen)12);
#line 422 "zgbtrf.f"
		    }
#line 423 "zgbtrf.f"
		}

#line 425 "zgbtrf.f"
		if (j3 > 0) {

/*                 Copy the lower triangle of A13 into the work array */
/*                 WORK13 */

#line 430 "zgbtrf.f"
		    i__3 = j3;
#line 430 "zgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 431 "zgbtrf.f"
			i__4 = jb;
#line 431 "zgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 432 "zgbtrf.f"
			    i__5 = ii + jj * 65 - 66;
#line 432 "zgbtrf.f"
			    i__6 = ii - jj + 1 + (jj + j + kv - 1) * ab_dim1;
#line 432 "zgbtrf.f"
			    work13[i__5].r = ab[i__6].r, work13[i__5].i = ab[
				    i__6].i;
#line 433 "zgbtrf.f"
/* L120: */
#line 433 "zgbtrf.f"
			}
#line 434 "zgbtrf.f"
/* L130: */
#line 434 "zgbtrf.f"
		    }

/*                 Update A13 in the work array */

#line 438 "zgbtrf.f"
		    i__3 = *ldab - 1;
#line 438 "zgbtrf.f"
		    ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
			    &c_b1, &ab[kv + 1 + j * ab_dim1], &i__3, work13, &
			    c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
			    4);

#line 442 "zgbtrf.f"
		    if (i2 > 0) {

/*                    Update A23 */

#line 446 "zgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 446 "zgbtrf.f"
			i__3 = *ldab - 1;
#line 446 "zgbtrf.f"
			i__4 = *ldab - 1;
#line 446 "zgbtrf.f"
			zgemm_("No transpose", "No transpose", &i2, &j3, &jb, 
				&z__1, &ab[kv + 1 + jb + j * ab_dim1], &i__3, 
				work13, &c__65, &c_b1, &ab[jb + 1 + (j + kv) *
				 ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
#line 450 "zgbtrf.f"
		    }

#line 452 "zgbtrf.f"
		    if (i3 > 0) {

/*                    Update A33 */

#line 456 "zgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 456 "zgbtrf.f"
			i__3 = *ldab - 1;
#line 456 "zgbtrf.f"
			zgemm_("No transpose", "No transpose", &i3, &j3, &jb, 
				&z__1, work31, &c__65, work13, &c__65, &c_b1, 
				&ab[*kl + 1 + (j + kv) * ab_dim1], &i__3, (
				ftnlen)12, (ftnlen)12);
#line 459 "zgbtrf.f"
		    }

/*                 Copy the lower triangle of A13 back into place */

#line 463 "zgbtrf.f"
		    i__3 = j3;
#line 463 "zgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 464 "zgbtrf.f"
			i__4 = jb;
#line 464 "zgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 465 "zgbtrf.f"
			    i__5 = ii - jj + 1 + (jj + j + kv - 1) * ab_dim1;
#line 465 "zgbtrf.f"
			    i__6 = ii + jj * 65 - 66;
#line 465 "zgbtrf.f"
			    ab[i__5].r = work13[i__6].r, ab[i__5].i = work13[
				    i__6].i;
#line 466 "zgbtrf.f"
/* L140: */
#line 466 "zgbtrf.f"
			}
#line 467 "zgbtrf.f"
/* L150: */
#line 467 "zgbtrf.f"
		    }
#line 468 "zgbtrf.f"
		}
#line 469 "zgbtrf.f"
	    } else {

/*              Adjust the pivot indices. */

#line 473 "zgbtrf.f"
		i__3 = j + jb - 1;
#line 473 "zgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 474 "zgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 475 "zgbtrf.f"
/* L160: */
#line 475 "zgbtrf.f"
		}
#line 476 "zgbtrf.f"
	    }

/*           Partially undo the interchanges in the current block to */
/*           restore the upper triangular form of A31 and copy the upper */
/*           triangle of A31 back into place */

#line 482 "zgbtrf.f"
	    i__3 = j;
#line 482 "zgbtrf.f"
	    for (jj = j + jb - 1; jj >= i__3; --jj) {
#line 483 "zgbtrf.f"
		jp = ipiv[jj] - jj + 1;
#line 484 "zgbtrf.f"
		if (jp != 1) {

/*                 Apply interchange to columns J to JJ-1 */

#line 488 "zgbtrf.f"
		    if (jp + jj - 1 < j + *kl) {

/*                    The interchange does not affect A31 */

#line 492 "zgbtrf.f"
			i__4 = jj - j;
#line 492 "zgbtrf.f"
			i__5 = *ldab - 1;
#line 492 "zgbtrf.f"
			i__6 = *ldab - 1;
#line 492 "zgbtrf.f"
			zswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
				i__6);
#line 494 "zgbtrf.f"
		    } else {

/*                    The interchange does affect A31 */

#line 498 "zgbtrf.f"
			i__4 = jj - j;
#line 498 "zgbtrf.f"
			i__5 = *ldab - 1;
#line 498 "zgbtrf.f"
			zswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &work31[jp + jj - j - *kl - 1], &c__65);
#line 500 "zgbtrf.f"
		    }
#line 501 "zgbtrf.f"
		}

/*              Copy the current column of A31 back into place */

/* Computing MIN */
#line 505 "zgbtrf.f"
		i__4 = i3, i__5 = jj - j + 1;
#line 505 "zgbtrf.f"
		nw = min(i__4,i__5);
#line 506 "zgbtrf.f"
		if (nw > 0) {
#line 506 "zgbtrf.f"
		    zcopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
			    kv + *kl + 1 - jj + j + jj * ab_dim1], &c__1);
#line 506 "zgbtrf.f"
		}
#line 509 "zgbtrf.f"
/* L170: */
#line 509 "zgbtrf.f"
	    }
#line 510 "zgbtrf.f"
/* L180: */
#line 510 "zgbtrf.f"
	}
#line 511 "zgbtrf.f"
    }

#line 513 "zgbtrf.f"
    return 0;

/*     End of ZGBTRF */

} /* zgbtrf_ */

