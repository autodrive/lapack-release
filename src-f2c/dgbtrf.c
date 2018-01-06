#line 1 "dgbtrf.f"
/* dgbtrf.f -- translated by f2c (version 20100827).
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

#line 1 "dgbtrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__65 = 65;
static doublereal c_b18 = -1.;
static doublereal c_b31 = 1.;

/* > \brief \b DGBTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */

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
/* > DGBTRF computes an LU factorization of a real m-by-n band matrix A */
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
/* >  elements of U because of fill-in resulting from the row interchanges. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, i2, i3, j2, j3, k2, jb, nb, ii, jj, jm, ip, jp, km,
	     ju, kv, nw;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dcopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    dswap_(integer *, doublereal *, integer *, doublereal *, integer *
	    );
    static doublereal work13[4160]	/* was [65][64] */, work31[4160]	
	    /* was [65][64] */;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgbtf2_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaswp_(integer *, doublereal *, integer *, 
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

#line 193 "dgbtrf.f"
    /* Parameter adjustments */
#line 193 "dgbtrf.f"
    ab_dim1 = *ldab;
#line 193 "dgbtrf.f"
    ab_offset = 1 + ab_dim1;
#line 193 "dgbtrf.f"
    ab -= ab_offset;
#line 193 "dgbtrf.f"
    --ipiv;
#line 193 "dgbtrf.f"

#line 193 "dgbtrf.f"
    /* Function Body */
#line 193 "dgbtrf.f"
    kv = *ku + *kl;

/*     Test the input parameters. */

#line 197 "dgbtrf.f"
    *info = 0;
#line 198 "dgbtrf.f"
    if (*m < 0) {
#line 199 "dgbtrf.f"
	*info = -1;
#line 200 "dgbtrf.f"
    } else if (*n < 0) {
#line 201 "dgbtrf.f"
	*info = -2;
#line 202 "dgbtrf.f"
    } else if (*kl < 0) {
#line 203 "dgbtrf.f"
	*info = -3;
#line 204 "dgbtrf.f"
    } else if (*ku < 0) {
#line 205 "dgbtrf.f"
	*info = -4;
#line 206 "dgbtrf.f"
    } else if (*ldab < *kl + kv + 1) {
#line 207 "dgbtrf.f"
	*info = -6;
#line 208 "dgbtrf.f"
    }
#line 209 "dgbtrf.f"
    if (*info != 0) {
#line 210 "dgbtrf.f"
	i__1 = -(*info);
#line 210 "dgbtrf.f"
	xerbla_("DGBTRF", &i__1, (ftnlen)6);
#line 211 "dgbtrf.f"
	return 0;
#line 212 "dgbtrf.f"
    }

/*     Quick return if possible */

#line 216 "dgbtrf.f"
    if (*m == 0 || *n == 0) {
#line 216 "dgbtrf.f"
	return 0;
#line 216 "dgbtrf.f"
    }

/*     Determine the block size for this environment */

#line 221 "dgbtrf.f"
    nb = ilaenv_(&c__1, "DGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);

/*     The block size must not exceed the limit set by the size of the */
/*     local arrays WORK13 and WORK31. */

#line 226 "dgbtrf.f"
    nb = min(nb,64);

#line 228 "dgbtrf.f"
    if (nb <= 1 || nb > *kl) {

/*        Use unblocked code */

#line 232 "dgbtrf.f"
	dgbtf2_(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
#line 233 "dgbtrf.f"
    } else {

/*        Use blocked code */

/*        Zero the superdiagonal elements of the work array WORK13 */

#line 239 "dgbtrf.f"
	i__1 = nb;
#line 239 "dgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 240 "dgbtrf.f"
	    i__2 = j - 1;
#line 240 "dgbtrf.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 241 "dgbtrf.f"
		work13[i__ + j * 65 - 66] = 0.;
#line 242 "dgbtrf.f"
/* L10: */
#line 242 "dgbtrf.f"
	    }
#line 243 "dgbtrf.f"
/* L20: */
#line 243 "dgbtrf.f"
	}

/*        Zero the subdiagonal elements of the work array WORK31 */

#line 247 "dgbtrf.f"
	i__1 = nb;
#line 247 "dgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 248 "dgbtrf.f"
	    i__2 = nb;
#line 248 "dgbtrf.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 249 "dgbtrf.f"
		work31[i__ + j * 65 - 66] = 0.;
#line 250 "dgbtrf.f"
/* L30: */
#line 250 "dgbtrf.f"
	    }
#line 251 "dgbtrf.f"
/* L40: */
#line 251 "dgbtrf.f"
	}

/*        Gaussian elimination with partial pivoting */

/*        Set fill-in elements in columns KU+2 to KV to zero */

#line 257 "dgbtrf.f"
	i__1 = min(kv,*n);
#line 257 "dgbtrf.f"
	for (j = *ku + 2; j <= i__1; ++j) {
#line 258 "dgbtrf.f"
	    i__2 = *kl;
#line 258 "dgbtrf.f"
	    for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
#line 259 "dgbtrf.f"
		ab[i__ + j * ab_dim1] = 0.;
#line 260 "dgbtrf.f"
/* L50: */
#line 260 "dgbtrf.f"
	    }
#line 261 "dgbtrf.f"
/* L60: */
#line 261 "dgbtrf.f"
	}

/*        JU is the index of the last column affected by the current */
/*        stage of the factorization */

#line 266 "dgbtrf.f"
	ju = 1;

#line 268 "dgbtrf.f"
	i__1 = min(*m,*n);
#line 268 "dgbtrf.f"
	i__2 = nb;
#line 268 "dgbtrf.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 269 "dgbtrf.f"
	    i__3 = nb, i__4 = min(*m,*n) - j + 1;
#line 269 "dgbtrf.f"
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
#line 283 "dgbtrf.f"
	    i__3 = *kl - jb, i__4 = *m - j - jb + 1;
#line 283 "dgbtrf.f"
	    i2 = min(i__3,i__4);
/* Computing MIN */
#line 284 "dgbtrf.f"
	    i__3 = jb, i__4 = *m - j - *kl + 1;
#line 284 "dgbtrf.f"
	    i3 = min(i__3,i__4);

/*           J2 and J3 are computed after JU has been updated. */

/*           Factorize the current block of JB columns */

#line 290 "dgbtrf.f"
	    i__3 = j + jb - 1;
#line 290 "dgbtrf.f"
	    for (jj = j; jj <= i__3; ++jj) {

/*              Set fill-in elements in column JJ+KV to zero */

#line 294 "dgbtrf.f"
		if (jj + kv <= *n) {
#line 295 "dgbtrf.f"
		    i__4 = *kl;
#line 295 "dgbtrf.f"
		    for (i__ = 1; i__ <= i__4; ++i__) {
#line 296 "dgbtrf.f"
			ab[i__ + (jj + kv) * ab_dim1] = 0.;
#line 297 "dgbtrf.f"
/* L70: */
#line 297 "dgbtrf.f"
		    }
#line 298 "dgbtrf.f"
		}

/*              Find pivot and test for singularity. KM is the number of */
/*              subdiagonal elements in the current column. */

/* Computing MIN */
#line 303 "dgbtrf.f"
		i__4 = *kl, i__5 = *m - jj;
#line 303 "dgbtrf.f"
		km = min(i__4,i__5);
#line 304 "dgbtrf.f"
		i__4 = km + 1;
#line 304 "dgbtrf.f"
		jp = idamax_(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
#line 305 "dgbtrf.f"
		ipiv[jj] = jp + jj - j;
#line 306 "dgbtrf.f"
		if (ab[kv + jp + jj * ab_dim1] != 0.) {
/* Computing MAX */
/* Computing MIN */
#line 307 "dgbtrf.f"
		    i__6 = jj + *ku + jp - 1;
#line 307 "dgbtrf.f"
		    i__4 = ju, i__5 = min(i__6,*n);
#line 307 "dgbtrf.f"
		    ju = max(i__4,i__5);
#line 308 "dgbtrf.f"
		    if (jp != 1) {

/*                    Apply interchange to columns J to J+JB-1 */

#line 312 "dgbtrf.f"
			if (jp + jj - 1 < j + *kl) {

#line 314 "dgbtrf.f"
			    i__4 = *ldab - 1;
#line 314 "dgbtrf.f"
			    i__5 = *ldab - 1;
#line 314 "dgbtrf.f"
			    dswap_(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
				    i__4, &ab[kv + jp + jj - j + j * ab_dim1],
				     &i__5);
#line 316 "dgbtrf.f"
			} else {

/*                       The interchange affects columns J to JJ-1 of A31 */
/*                       which are stored in the work array WORK31 */

#line 321 "dgbtrf.f"
			    i__4 = jj - j;
#line 321 "dgbtrf.f"
			    i__5 = *ldab - 1;
#line 321 "dgbtrf.f"
			    dswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], 
				    &i__5, &work31[jp + jj - j - *kl - 1], &
				    c__65);
#line 323 "dgbtrf.f"
			    i__4 = j + jb - jj;
#line 323 "dgbtrf.f"
			    i__5 = *ldab - 1;
#line 323 "dgbtrf.f"
			    i__6 = *ldab - 1;
#line 323 "dgbtrf.f"
			    dswap_(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
				    ab[kv + jp + jj * ab_dim1], &i__6);
#line 325 "dgbtrf.f"
			}
#line 326 "dgbtrf.f"
		    }

/*                 Compute multipliers */

#line 330 "dgbtrf.f"
		    d__1 = 1. / ab[kv + 1 + jj * ab_dim1];
#line 330 "dgbtrf.f"
		    dscal_(&km, &d__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

/*                 Update trailing submatrix within the band and within */
/*                 the current block. JM is the index of the last column */
/*                 which needs to be updated. */

/* Computing MIN */
#line 337 "dgbtrf.f"
		    i__4 = ju, i__5 = j + jb - 1;
#line 337 "dgbtrf.f"
		    jm = min(i__4,i__5);
#line 338 "dgbtrf.f"
		    if (jm > jj) {
#line 338 "dgbtrf.f"
			i__4 = jm - jj;
#line 338 "dgbtrf.f"
			i__5 = *ldab - 1;
#line 338 "dgbtrf.f"
			i__6 = *ldab - 1;
#line 338 "dgbtrf.f"
			dger_(&km, &i__4, &c_b18, &ab[kv + 2 + jj * ab_dim1], 
				&c__1, &ab[kv + (jj + 1) * ab_dim1], &i__5, &
				ab[kv + 1 + (jj + 1) * ab_dim1], &i__6);
#line 338 "dgbtrf.f"
		    }
#line 342 "dgbtrf.f"
		} else {

/*                 If pivot is zero, set INFO to the index of the pivot */
/*                 unless a zero pivot has already been found. */

#line 347 "dgbtrf.f"
		    if (*info == 0) {
#line 347 "dgbtrf.f"
			*info = jj;
#line 347 "dgbtrf.f"
		    }
#line 349 "dgbtrf.f"
		}

/*              Copy current column of A31 into the work array WORK31 */

/* Computing MIN */
#line 353 "dgbtrf.f"
		i__4 = jj - j + 1;
#line 353 "dgbtrf.f"
		nw = min(i__4,i3);
#line 354 "dgbtrf.f"
		if (nw > 0) {
#line 354 "dgbtrf.f"
		    dcopy_(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
			    c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
#line 354 "dgbtrf.f"
		}
#line 357 "dgbtrf.f"
/* L80: */
#line 357 "dgbtrf.f"
	    }
#line 358 "dgbtrf.f"
	    if (j + jb <= *n) {

/*              Apply the row interchanges to the other blocks. */

/* Computing MIN */
#line 362 "dgbtrf.f"
		i__3 = ju - j + 1;
#line 362 "dgbtrf.f"
		j2 = min(i__3,kv) - jb;
/* Computing MAX */
#line 363 "dgbtrf.f"
		i__3 = 0, i__4 = ju - j - kv + 1;
#line 363 "dgbtrf.f"
		j3 = max(i__3,i__4);

/*              Use DLASWP to apply the row interchanges to A12, A22, and */
/*              A32. */

#line 368 "dgbtrf.f"
		i__3 = *ldab - 1;
#line 368 "dgbtrf.f"
		dlaswp_(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
			c__1, &jb, &ipiv[j], &c__1);

/*              Adjust the pivot indices. */

#line 373 "dgbtrf.f"
		i__3 = j + jb - 1;
#line 373 "dgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 374 "dgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 375 "dgbtrf.f"
/* L90: */
#line 375 "dgbtrf.f"
		}

/*              Apply the row interchanges to A13, A23, and A33 */
/*              columnwise. */

#line 380 "dgbtrf.f"
		k2 = j - 1 + jb + j2;
#line 381 "dgbtrf.f"
		i__3 = j3;
#line 381 "dgbtrf.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 382 "dgbtrf.f"
		    jj = k2 + i__;
#line 383 "dgbtrf.f"
		    i__4 = j + jb - 1;
#line 383 "dgbtrf.f"
		    for (ii = j + i__ - 1; ii <= i__4; ++ii) {
#line 384 "dgbtrf.f"
			ip = ipiv[ii];
#line 385 "dgbtrf.f"
			if (ip != ii) {
#line 386 "dgbtrf.f"
			    temp = ab[kv + 1 + ii - jj + jj * ab_dim1];
#line 387 "dgbtrf.f"
			    ab[kv + 1 + ii - jj + jj * ab_dim1] = ab[kv + 1 + 
				    ip - jj + jj * ab_dim1];
#line 388 "dgbtrf.f"
			    ab[kv + 1 + ip - jj + jj * ab_dim1] = temp;
#line 389 "dgbtrf.f"
			}
#line 390 "dgbtrf.f"
/* L100: */
#line 390 "dgbtrf.f"
		    }
#line 391 "dgbtrf.f"
/* L110: */
#line 391 "dgbtrf.f"
		}

/*              Update the relevant part of the trailing submatrix */

#line 395 "dgbtrf.f"
		if (j2 > 0) {

/*                 Update A12 */

#line 399 "dgbtrf.f"
		    i__3 = *ldab - 1;
#line 399 "dgbtrf.f"
		    i__4 = *ldab - 1;
#line 399 "dgbtrf.f"
		    dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
			    &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv 
			    + 1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4, 
			    (ftnlen)5, (ftnlen)12, (ftnlen)4);

#line 403 "dgbtrf.f"
		    if (i2 > 0) {

/*                    Update A22 */

#line 407 "dgbtrf.f"
			i__3 = *ldab - 1;
#line 407 "dgbtrf.f"
			i__4 = *ldab - 1;
#line 407 "dgbtrf.f"
			i__5 = *ldab - 1;
#line 407 "dgbtrf.f"
			dgemm_("No transpose", "No transpose", &i2, &j2, &jb, 
				&c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
				 &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4,
				 &c_b31, &ab[kv + 1 + (j + jb) * ab_dim1], &
				i__5, (ftnlen)12, (ftnlen)12);
#line 411 "dgbtrf.f"
		    }

#line 413 "dgbtrf.f"
		    if (i3 > 0) {

/*                    Update A32 */

#line 417 "dgbtrf.f"
			i__3 = *ldab - 1;
#line 417 "dgbtrf.f"
			i__4 = *ldab - 1;
#line 417 "dgbtrf.f"
			dgemm_("No transpose", "No transpose", &i3, &j2, &jb, 
				&c_b18, work31, &c__65, &ab[kv + 1 - jb + (j 
				+ jb) * ab_dim1], &i__3, &c_b31, &ab[kv + *kl 
				+ 1 - jb + (j + jb) * ab_dim1], &i__4, (
				ftnlen)12, (ftnlen)12);
#line 421 "dgbtrf.f"
		    }
#line 422 "dgbtrf.f"
		}

#line 424 "dgbtrf.f"
		if (j3 > 0) {

/*                 Copy the lower triangle of A13 into the work array */
/*                 WORK13 */

#line 429 "dgbtrf.f"
		    i__3 = j3;
#line 429 "dgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 430 "dgbtrf.f"
			i__4 = jb;
#line 430 "dgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 431 "dgbtrf.f"
			    work13[ii + jj * 65 - 66] = ab[ii - jj + 1 + (jj 
				    + j + kv - 1) * ab_dim1];
#line 432 "dgbtrf.f"
/* L120: */
#line 432 "dgbtrf.f"
			}
#line 433 "dgbtrf.f"
/* L130: */
#line 433 "dgbtrf.f"
		    }

/*                 Update A13 in the work array */

#line 437 "dgbtrf.f"
		    i__3 = *ldab - 1;
#line 437 "dgbtrf.f"
		    dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
			    &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, work13, 
			    &c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
			    4);

#line 441 "dgbtrf.f"
		    if (i2 > 0) {

/*                    Update A23 */

#line 445 "dgbtrf.f"
			i__3 = *ldab - 1;
#line 445 "dgbtrf.f"
			i__4 = *ldab - 1;
#line 445 "dgbtrf.f"
			dgemm_("No transpose", "No transpose", &i2, &j3, &jb, 
				&c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
				 work13, &c__65, &c_b31, &ab[jb + 1 + (j + kv)
				 * ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
#line 449 "dgbtrf.f"
		    }

#line 451 "dgbtrf.f"
		    if (i3 > 0) {

/*                    Update A33 */

#line 455 "dgbtrf.f"
			i__3 = *ldab - 1;
#line 455 "dgbtrf.f"
			dgemm_("No transpose", "No transpose", &i3, &j3, &jb, 
				&c_b18, work31, &c__65, work13, &c__65, &
				c_b31, &ab[*kl + 1 + (j + kv) * ab_dim1], &
				i__3, (ftnlen)12, (ftnlen)12);
#line 458 "dgbtrf.f"
		    }

/*                 Copy the lower triangle of A13 back into place */

#line 462 "dgbtrf.f"
		    i__3 = j3;
#line 462 "dgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 463 "dgbtrf.f"
			i__4 = jb;
#line 463 "dgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 464 "dgbtrf.f"
			    ab[ii - jj + 1 + (jj + j + kv - 1) * ab_dim1] = 
				    work13[ii + jj * 65 - 66];
#line 465 "dgbtrf.f"
/* L140: */
#line 465 "dgbtrf.f"
			}
#line 466 "dgbtrf.f"
/* L150: */
#line 466 "dgbtrf.f"
		    }
#line 467 "dgbtrf.f"
		}
#line 468 "dgbtrf.f"
	    } else {

/*              Adjust the pivot indices. */

#line 472 "dgbtrf.f"
		i__3 = j + jb - 1;
#line 472 "dgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 473 "dgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 474 "dgbtrf.f"
/* L160: */
#line 474 "dgbtrf.f"
		}
#line 475 "dgbtrf.f"
	    }

/*           Partially undo the interchanges in the current block to */
/*           restore the upper triangular form of A31 and copy the upper */
/*           triangle of A31 back into place */

#line 481 "dgbtrf.f"
	    i__3 = j;
#line 481 "dgbtrf.f"
	    for (jj = j + jb - 1; jj >= i__3; --jj) {
#line 482 "dgbtrf.f"
		jp = ipiv[jj] - jj + 1;
#line 483 "dgbtrf.f"
		if (jp != 1) {

/*                 Apply interchange to columns J to JJ-1 */

#line 487 "dgbtrf.f"
		    if (jp + jj - 1 < j + *kl) {

/*                    The interchange does not affect A31 */

#line 491 "dgbtrf.f"
			i__4 = jj - j;
#line 491 "dgbtrf.f"
			i__5 = *ldab - 1;
#line 491 "dgbtrf.f"
			i__6 = *ldab - 1;
#line 491 "dgbtrf.f"
			dswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
				i__6);
#line 493 "dgbtrf.f"
		    } else {

/*                    The interchange does affect A31 */

#line 497 "dgbtrf.f"
			i__4 = jj - j;
#line 497 "dgbtrf.f"
			i__5 = *ldab - 1;
#line 497 "dgbtrf.f"
			dswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &work31[jp + jj - j - *kl - 1], &c__65);
#line 499 "dgbtrf.f"
		    }
#line 500 "dgbtrf.f"
		}

/*              Copy the current column of A31 back into place */

/* Computing MIN */
#line 504 "dgbtrf.f"
		i__4 = i3, i__5 = jj - j + 1;
#line 504 "dgbtrf.f"
		nw = min(i__4,i__5);
#line 505 "dgbtrf.f"
		if (nw > 0) {
#line 505 "dgbtrf.f"
		    dcopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
			    kv + *kl + 1 - jj + j + jj * ab_dim1], &c__1);
#line 505 "dgbtrf.f"
		}
#line 508 "dgbtrf.f"
/* L170: */
#line 508 "dgbtrf.f"
	    }
#line 509 "dgbtrf.f"
/* L180: */
#line 509 "dgbtrf.f"
	}
#line 510 "dgbtrf.f"
    }

#line 512 "dgbtrf.f"
    return 0;

/*     End of DGBTRF */

} /* dgbtrf_ */

