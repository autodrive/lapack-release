#line 1 "cgbtrf.f"
/* cgbtrf.f -- translated by f2c (version 20100827).
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

#line 1 "cgbtrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c__65 = 65;

/* > \brief \b CGBTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbtrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbtrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbtrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBTRF computes an LU factorization of a complex m-by-n band matrix A */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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

/* > \ingroup complexGBcomputational */

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
/* Subroutine */ int cgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cgeru_(integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ccopy_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), cswap_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    ;
    static doublecomplex work13[4160]	/* was [65][64] */, work31[4160]	
	    /* was [65][64] */;
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    cgbtf2_(integer *, integer *, integer *, integer *, doublecomplex 
	    *, integer *, integer *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int claswp_(integer *, doublecomplex *, integer *,
	     integer *, integer *, integer *, integer *);


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

#line 194 "cgbtrf.f"
    /* Parameter adjustments */
#line 194 "cgbtrf.f"
    ab_dim1 = *ldab;
#line 194 "cgbtrf.f"
    ab_offset = 1 + ab_dim1;
#line 194 "cgbtrf.f"
    ab -= ab_offset;
#line 194 "cgbtrf.f"
    --ipiv;
#line 194 "cgbtrf.f"

#line 194 "cgbtrf.f"
    /* Function Body */
#line 194 "cgbtrf.f"
    kv = *ku + *kl;

/*     Test the input parameters. */

#line 198 "cgbtrf.f"
    *info = 0;
#line 199 "cgbtrf.f"
    if (*m < 0) {
#line 200 "cgbtrf.f"
	*info = -1;
#line 201 "cgbtrf.f"
    } else if (*n < 0) {
#line 202 "cgbtrf.f"
	*info = -2;
#line 203 "cgbtrf.f"
    } else if (*kl < 0) {
#line 204 "cgbtrf.f"
	*info = -3;
#line 205 "cgbtrf.f"
    } else if (*ku < 0) {
#line 206 "cgbtrf.f"
	*info = -4;
#line 207 "cgbtrf.f"
    } else if (*ldab < *kl + kv + 1) {
#line 208 "cgbtrf.f"
	*info = -6;
#line 209 "cgbtrf.f"
    }
#line 210 "cgbtrf.f"
    if (*info != 0) {
#line 211 "cgbtrf.f"
	i__1 = -(*info);
#line 211 "cgbtrf.f"
	xerbla_("CGBTRF", &i__1, (ftnlen)6);
#line 212 "cgbtrf.f"
	return 0;
#line 213 "cgbtrf.f"
    }

/*     Quick return if possible */

#line 217 "cgbtrf.f"
    if (*m == 0 || *n == 0) {
#line 217 "cgbtrf.f"
	return 0;
#line 217 "cgbtrf.f"
    }

/*     Determine the block size for this environment */

#line 222 "cgbtrf.f"
    nb = ilaenv_(&c__1, "CGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);

/*     The block size must not exceed the limit set by the size of the */
/*     local arrays WORK13 and WORK31. */

#line 227 "cgbtrf.f"
    nb = min(nb,64);

#line 229 "cgbtrf.f"
    if (nb <= 1 || nb > *kl) {

/*        Use unblocked code */

#line 233 "cgbtrf.f"
	cgbtf2_(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
#line 234 "cgbtrf.f"
    } else {

/*        Use blocked code */

/*        Zero the superdiagonal elements of the work array WORK13 */

#line 240 "cgbtrf.f"
	i__1 = nb;
#line 240 "cgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 241 "cgbtrf.f"
	    i__2 = j - 1;
#line 241 "cgbtrf.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 242 "cgbtrf.f"
		i__3 = i__ + j * 65 - 66;
#line 242 "cgbtrf.f"
		work13[i__3].r = 0., work13[i__3].i = 0.;
#line 243 "cgbtrf.f"
/* L10: */
#line 243 "cgbtrf.f"
	    }
#line 244 "cgbtrf.f"
/* L20: */
#line 244 "cgbtrf.f"
	}

/*        Zero the subdiagonal elements of the work array WORK31 */

#line 248 "cgbtrf.f"
	i__1 = nb;
#line 248 "cgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 249 "cgbtrf.f"
	    i__2 = nb;
#line 249 "cgbtrf.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 250 "cgbtrf.f"
		i__3 = i__ + j * 65 - 66;
#line 250 "cgbtrf.f"
		work31[i__3].r = 0., work31[i__3].i = 0.;
#line 251 "cgbtrf.f"
/* L30: */
#line 251 "cgbtrf.f"
	    }
#line 252 "cgbtrf.f"
/* L40: */
#line 252 "cgbtrf.f"
	}

/*        Gaussian elimination with partial pivoting */

/*        Set fill-in elements in columns KU+2 to KV to zero */

#line 258 "cgbtrf.f"
	i__1 = min(kv,*n);
#line 258 "cgbtrf.f"
	for (j = *ku + 2; j <= i__1; ++j) {
#line 259 "cgbtrf.f"
	    i__2 = *kl;
#line 259 "cgbtrf.f"
	    for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
#line 260 "cgbtrf.f"
		i__3 = i__ + j * ab_dim1;
#line 260 "cgbtrf.f"
		ab[i__3].r = 0., ab[i__3].i = 0.;
#line 261 "cgbtrf.f"
/* L50: */
#line 261 "cgbtrf.f"
	    }
#line 262 "cgbtrf.f"
/* L60: */
#line 262 "cgbtrf.f"
	}

/*        JU is the index of the last column affected by the current */
/*        stage of the factorization */

#line 267 "cgbtrf.f"
	ju = 1;

#line 269 "cgbtrf.f"
	i__1 = min(*m,*n);
#line 269 "cgbtrf.f"
	i__2 = nb;
#line 269 "cgbtrf.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 270 "cgbtrf.f"
	    i__3 = nb, i__4 = min(*m,*n) - j + 1;
#line 270 "cgbtrf.f"
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
#line 284 "cgbtrf.f"
	    i__3 = *kl - jb, i__4 = *m - j - jb + 1;
#line 284 "cgbtrf.f"
	    i2 = min(i__3,i__4);
/* Computing MIN */
#line 285 "cgbtrf.f"
	    i__3 = jb, i__4 = *m - j - *kl + 1;
#line 285 "cgbtrf.f"
	    i3 = min(i__3,i__4);

/*           J2 and J3 are computed after JU has been updated. */

/*           Factorize the current block of JB columns */

#line 291 "cgbtrf.f"
	    i__3 = j + jb - 1;
#line 291 "cgbtrf.f"
	    for (jj = j; jj <= i__3; ++jj) {

/*              Set fill-in elements in column JJ+KV to zero */

#line 295 "cgbtrf.f"
		if (jj + kv <= *n) {
#line 296 "cgbtrf.f"
		    i__4 = *kl;
#line 296 "cgbtrf.f"
		    for (i__ = 1; i__ <= i__4; ++i__) {
#line 297 "cgbtrf.f"
			i__5 = i__ + (jj + kv) * ab_dim1;
#line 297 "cgbtrf.f"
			ab[i__5].r = 0., ab[i__5].i = 0.;
#line 298 "cgbtrf.f"
/* L70: */
#line 298 "cgbtrf.f"
		    }
#line 299 "cgbtrf.f"
		}

/*              Find pivot and test for singularity. KM is the number of */
/*              subdiagonal elements in the current column. */

/* Computing MIN */
#line 304 "cgbtrf.f"
		i__4 = *kl, i__5 = *m - jj;
#line 304 "cgbtrf.f"
		km = min(i__4,i__5);
#line 305 "cgbtrf.f"
		i__4 = km + 1;
#line 305 "cgbtrf.f"
		jp = icamax_(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
#line 306 "cgbtrf.f"
		ipiv[jj] = jp + jj - j;
#line 307 "cgbtrf.f"
		i__4 = kv + jp + jj * ab_dim1;
#line 307 "cgbtrf.f"
		if (ab[i__4].r != 0. || ab[i__4].i != 0.) {
/* Computing MAX */
/* Computing MIN */
#line 308 "cgbtrf.f"
		    i__6 = jj + *ku + jp - 1;
#line 308 "cgbtrf.f"
		    i__4 = ju, i__5 = min(i__6,*n);
#line 308 "cgbtrf.f"
		    ju = max(i__4,i__5);
#line 309 "cgbtrf.f"
		    if (jp != 1) {

/*                    Apply interchange to columns J to J+JB-1 */

#line 313 "cgbtrf.f"
			if (jp + jj - 1 < j + *kl) {

#line 315 "cgbtrf.f"
			    i__4 = *ldab - 1;
#line 315 "cgbtrf.f"
			    i__5 = *ldab - 1;
#line 315 "cgbtrf.f"
			    cswap_(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
				    i__4, &ab[kv + jp + jj - j + j * ab_dim1],
				     &i__5);
#line 317 "cgbtrf.f"
			} else {

/*                       The interchange affects columns J to JJ-1 of A31 */
/*                       which are stored in the work array WORK31 */

#line 322 "cgbtrf.f"
			    i__4 = jj - j;
#line 322 "cgbtrf.f"
			    i__5 = *ldab - 1;
#line 322 "cgbtrf.f"
			    cswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], 
				    &i__5, &work31[jp + jj - j - *kl - 1], &
				    c__65);
#line 324 "cgbtrf.f"
			    i__4 = j + jb - jj;
#line 324 "cgbtrf.f"
			    i__5 = *ldab - 1;
#line 324 "cgbtrf.f"
			    i__6 = *ldab - 1;
#line 324 "cgbtrf.f"
			    cswap_(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
				    ab[kv + jp + jj * ab_dim1], &i__6);
#line 326 "cgbtrf.f"
			}
#line 327 "cgbtrf.f"
		    }

/*                 Compute multipliers */

#line 331 "cgbtrf.f"
		    z_div(&z__1, &c_b1, &ab[kv + 1 + jj * ab_dim1]);
#line 331 "cgbtrf.f"
		    cscal_(&km, &z__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

/*                 Update trailing submatrix within the band and within */
/*                 the current block. JM is the index of the last column */
/*                 which needs to be updated. */

/* Computing MIN */
#line 338 "cgbtrf.f"
		    i__4 = ju, i__5 = j + jb - 1;
#line 338 "cgbtrf.f"
		    jm = min(i__4,i__5);
#line 339 "cgbtrf.f"
		    if (jm > jj) {
#line 339 "cgbtrf.f"
			i__4 = jm - jj;
#line 339 "cgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 339 "cgbtrf.f"
			i__5 = *ldab - 1;
#line 339 "cgbtrf.f"
			i__6 = *ldab - 1;
#line 339 "cgbtrf.f"
			cgeru_(&km, &i__4, &z__1, &ab[kv + 2 + jj * ab_dim1], 
				&c__1, &ab[kv + (jj + 1) * ab_dim1], &i__5, &
				ab[kv + 1 + (jj + 1) * ab_dim1], &i__6);
#line 339 "cgbtrf.f"
		    }
#line 343 "cgbtrf.f"
		} else {

/*                 If pivot is zero, set INFO to the index of the pivot */
/*                 unless a zero pivot has already been found. */

#line 348 "cgbtrf.f"
		    if (*info == 0) {
#line 348 "cgbtrf.f"
			*info = jj;
#line 348 "cgbtrf.f"
		    }
#line 350 "cgbtrf.f"
		}

/*              Copy current column of A31 into the work array WORK31 */

/* Computing MIN */
#line 354 "cgbtrf.f"
		i__4 = jj - j + 1;
#line 354 "cgbtrf.f"
		nw = min(i__4,i3);
#line 355 "cgbtrf.f"
		if (nw > 0) {
#line 355 "cgbtrf.f"
		    ccopy_(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
			    c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
#line 355 "cgbtrf.f"
		}
#line 358 "cgbtrf.f"
/* L80: */
#line 358 "cgbtrf.f"
	    }
#line 359 "cgbtrf.f"
	    if (j + jb <= *n) {

/*              Apply the row interchanges to the other blocks. */

/* Computing MIN */
#line 363 "cgbtrf.f"
		i__3 = ju - j + 1;
#line 363 "cgbtrf.f"
		j2 = min(i__3,kv) - jb;
/* Computing MAX */
#line 364 "cgbtrf.f"
		i__3 = 0, i__4 = ju - j - kv + 1;
#line 364 "cgbtrf.f"
		j3 = max(i__3,i__4);

/*              Use CLASWP to apply the row interchanges to A12, A22, and */
/*              A32. */

#line 369 "cgbtrf.f"
		i__3 = *ldab - 1;
#line 369 "cgbtrf.f"
		claswp_(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
			c__1, &jb, &ipiv[j], &c__1);

/*              Adjust the pivot indices. */

#line 374 "cgbtrf.f"
		i__3 = j + jb - 1;
#line 374 "cgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 375 "cgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 376 "cgbtrf.f"
/* L90: */
#line 376 "cgbtrf.f"
		}

/*              Apply the row interchanges to A13, A23, and A33 */
/*              columnwise. */

#line 381 "cgbtrf.f"
		k2 = j - 1 + jb + j2;
#line 382 "cgbtrf.f"
		i__3 = j3;
#line 382 "cgbtrf.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 383 "cgbtrf.f"
		    jj = k2 + i__;
#line 384 "cgbtrf.f"
		    i__4 = j + jb - 1;
#line 384 "cgbtrf.f"
		    for (ii = j + i__ - 1; ii <= i__4; ++ii) {
#line 385 "cgbtrf.f"
			ip = ipiv[ii];
#line 386 "cgbtrf.f"
			if (ip != ii) {
#line 387 "cgbtrf.f"
			    i__5 = kv + 1 + ii - jj + jj * ab_dim1;
#line 387 "cgbtrf.f"
			    temp.r = ab[i__5].r, temp.i = ab[i__5].i;
#line 388 "cgbtrf.f"
			    i__5 = kv + 1 + ii - jj + jj * ab_dim1;
#line 388 "cgbtrf.f"
			    i__6 = kv + 1 + ip - jj + jj * ab_dim1;
#line 388 "cgbtrf.f"
			    ab[i__5].r = ab[i__6].r, ab[i__5].i = ab[i__6].i;
#line 389 "cgbtrf.f"
			    i__5 = kv + 1 + ip - jj + jj * ab_dim1;
#line 389 "cgbtrf.f"
			    ab[i__5].r = temp.r, ab[i__5].i = temp.i;
#line 390 "cgbtrf.f"
			}
#line 391 "cgbtrf.f"
/* L100: */
#line 391 "cgbtrf.f"
		    }
#line 392 "cgbtrf.f"
/* L110: */
#line 392 "cgbtrf.f"
		}

/*              Update the relevant part of the trailing submatrix */

#line 396 "cgbtrf.f"
		if (j2 > 0) {

/*                 Update A12 */

#line 400 "cgbtrf.f"
		    i__3 = *ldab - 1;
#line 400 "cgbtrf.f"
		    i__4 = *ldab - 1;
#line 400 "cgbtrf.f"
		    ctrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
			    &c_b1, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv + 
			    1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4, (
			    ftnlen)5, (ftnlen)12, (ftnlen)4);

#line 404 "cgbtrf.f"
		    if (i2 > 0) {

/*                    Update A22 */

#line 408 "cgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 408 "cgbtrf.f"
			i__3 = *ldab - 1;
#line 408 "cgbtrf.f"
			i__4 = *ldab - 1;
#line 408 "cgbtrf.f"
			i__5 = *ldab - 1;
#line 408 "cgbtrf.f"
			cgemm_("No transpose", "No transpose", &i2, &j2, &jb, 
				&z__1, &ab[kv + 1 + jb + j * ab_dim1], &i__3, 
				&ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4, 
				&c_b1, &ab[kv + 1 + (j + jb) * ab_dim1], &
				i__5, (ftnlen)12, (ftnlen)12);
#line 412 "cgbtrf.f"
		    }

#line 414 "cgbtrf.f"
		    if (i3 > 0) {

/*                    Update A32 */

#line 418 "cgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 418 "cgbtrf.f"
			i__3 = *ldab - 1;
#line 418 "cgbtrf.f"
			i__4 = *ldab - 1;
#line 418 "cgbtrf.f"
			cgemm_("No transpose", "No transpose", &i3, &j2, &jb, 
				&z__1, work31, &c__65, &ab[kv + 1 - jb + (j + 
				jb) * ab_dim1], &i__3, &c_b1, &ab[kv + *kl + 
				1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)
				12, (ftnlen)12);
#line 422 "cgbtrf.f"
		    }
#line 423 "cgbtrf.f"
		}

#line 425 "cgbtrf.f"
		if (j3 > 0) {

/*                 Copy the lower triangle of A13 into the work array */
/*                 WORK13 */

#line 430 "cgbtrf.f"
		    i__3 = j3;
#line 430 "cgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 431 "cgbtrf.f"
			i__4 = jb;
#line 431 "cgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 432 "cgbtrf.f"
			    i__5 = ii + jj * 65 - 66;
#line 432 "cgbtrf.f"
			    i__6 = ii - jj + 1 + (jj + j + kv - 1) * ab_dim1;
#line 432 "cgbtrf.f"
			    work13[i__5].r = ab[i__6].r, work13[i__5].i = ab[
				    i__6].i;
#line 433 "cgbtrf.f"
/* L120: */
#line 433 "cgbtrf.f"
			}
#line 434 "cgbtrf.f"
/* L130: */
#line 434 "cgbtrf.f"
		    }

/*                 Update A13 in the work array */

#line 438 "cgbtrf.f"
		    i__3 = *ldab - 1;
#line 438 "cgbtrf.f"
		    ctrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
			    &c_b1, &ab[kv + 1 + j * ab_dim1], &i__3, work13, &
			    c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
			    4);

#line 442 "cgbtrf.f"
		    if (i2 > 0) {

/*                    Update A23 */

#line 446 "cgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 446 "cgbtrf.f"
			i__3 = *ldab - 1;
#line 446 "cgbtrf.f"
			i__4 = *ldab - 1;
#line 446 "cgbtrf.f"
			cgemm_("No transpose", "No transpose", &i2, &j3, &jb, 
				&z__1, &ab[kv + 1 + jb + j * ab_dim1], &i__3, 
				work13, &c__65, &c_b1, &ab[jb + 1 + (j + kv) *
				 ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
#line 450 "cgbtrf.f"
		    }

#line 452 "cgbtrf.f"
		    if (i3 > 0) {

/*                    Update A33 */

#line 456 "cgbtrf.f"
			z__1.r = -1., z__1.i = -0.;
#line 456 "cgbtrf.f"
			i__3 = *ldab - 1;
#line 456 "cgbtrf.f"
			cgemm_("No transpose", "No transpose", &i3, &j3, &jb, 
				&z__1, work31, &c__65, work13, &c__65, &c_b1, 
				&ab[*kl + 1 + (j + kv) * ab_dim1], &i__3, (
				ftnlen)12, (ftnlen)12);
#line 459 "cgbtrf.f"
		    }

/*                 Copy the lower triangle of A13 back into place */

#line 463 "cgbtrf.f"
		    i__3 = j3;
#line 463 "cgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 464 "cgbtrf.f"
			i__4 = jb;
#line 464 "cgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 465 "cgbtrf.f"
			    i__5 = ii - jj + 1 + (jj + j + kv - 1) * ab_dim1;
#line 465 "cgbtrf.f"
			    i__6 = ii + jj * 65 - 66;
#line 465 "cgbtrf.f"
			    ab[i__5].r = work13[i__6].r, ab[i__5].i = work13[
				    i__6].i;
#line 466 "cgbtrf.f"
/* L140: */
#line 466 "cgbtrf.f"
			}
#line 467 "cgbtrf.f"
/* L150: */
#line 467 "cgbtrf.f"
		    }
#line 468 "cgbtrf.f"
		}
#line 469 "cgbtrf.f"
	    } else {

/*              Adjust the pivot indices. */

#line 473 "cgbtrf.f"
		i__3 = j + jb - 1;
#line 473 "cgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 474 "cgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 475 "cgbtrf.f"
/* L160: */
#line 475 "cgbtrf.f"
		}
#line 476 "cgbtrf.f"
	    }

/*           Partially undo the interchanges in the current block to */
/*           restore the upper triangular form of A31 and copy the upper */
/*           triangle of A31 back into place */

#line 482 "cgbtrf.f"
	    i__3 = j;
#line 482 "cgbtrf.f"
	    for (jj = j + jb - 1; jj >= i__3; --jj) {
#line 483 "cgbtrf.f"
		jp = ipiv[jj] - jj + 1;
#line 484 "cgbtrf.f"
		if (jp != 1) {

/*                 Apply interchange to columns J to JJ-1 */

#line 488 "cgbtrf.f"
		    if (jp + jj - 1 < j + *kl) {

/*                    The interchange does not affect A31 */

#line 492 "cgbtrf.f"
			i__4 = jj - j;
#line 492 "cgbtrf.f"
			i__5 = *ldab - 1;
#line 492 "cgbtrf.f"
			i__6 = *ldab - 1;
#line 492 "cgbtrf.f"
			cswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
				i__6);
#line 494 "cgbtrf.f"
		    } else {

/*                    The interchange does affect A31 */

#line 498 "cgbtrf.f"
			i__4 = jj - j;
#line 498 "cgbtrf.f"
			i__5 = *ldab - 1;
#line 498 "cgbtrf.f"
			cswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &work31[jp + jj - j - *kl - 1], &c__65);
#line 500 "cgbtrf.f"
		    }
#line 501 "cgbtrf.f"
		}

/*              Copy the current column of A31 back into place */

/* Computing MIN */
#line 505 "cgbtrf.f"
		i__4 = i3, i__5 = jj - j + 1;
#line 505 "cgbtrf.f"
		nw = min(i__4,i__5);
#line 506 "cgbtrf.f"
		if (nw > 0) {
#line 506 "cgbtrf.f"
		    ccopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
			    kv + *kl + 1 - jj + j + jj * ab_dim1], &c__1);
#line 506 "cgbtrf.f"
		}
#line 509 "cgbtrf.f"
/* L170: */
#line 509 "cgbtrf.f"
	    }
#line 510 "cgbtrf.f"
/* L180: */
#line 510 "cgbtrf.f"
	}
#line 511 "cgbtrf.f"
    }

#line 513 "cgbtrf.f"
    return 0;

/*     End of CGBTRF */

} /* cgbtrf_ */

