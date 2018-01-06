#line 1 "sgbtrf.f"
/* sgbtrf.f -- translated by f2c (version 20100827).
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

#line 1 "sgbtrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__65 = 65;
static doublereal c_b18 = -1.;
static doublereal c_b31 = 1.;

/* > \brief \b SGBTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbtrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbtrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbtrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBTRF computes an LU factorization of a real m-by-n band matrix A */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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

/* > \ingroup realGBcomputational */

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
/* Subroutine */ int sgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, i2, i3, j2, j3, k2, jb, nb, ii, jj, jm, ip, jp, km,
	     ju, kv, nw;
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal work13[4160]	/* was [65][64] */, work31[4160]	
	    /* was [65][64] */;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), strsm_(char *, char *, char *, char *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), sgbtf2_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen), isamax_(integer *, 
	    doublereal *, integer *);
    extern /* Subroutine */ int slaswp_(integer *, doublereal *, integer *, 
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

#line 193 "sgbtrf.f"
    /* Parameter adjustments */
#line 193 "sgbtrf.f"
    ab_dim1 = *ldab;
#line 193 "sgbtrf.f"
    ab_offset = 1 + ab_dim1;
#line 193 "sgbtrf.f"
    ab -= ab_offset;
#line 193 "sgbtrf.f"
    --ipiv;
#line 193 "sgbtrf.f"

#line 193 "sgbtrf.f"
    /* Function Body */
#line 193 "sgbtrf.f"
    kv = *ku + *kl;

/*     Test the input parameters. */

#line 197 "sgbtrf.f"
    *info = 0;
#line 198 "sgbtrf.f"
    if (*m < 0) {
#line 199 "sgbtrf.f"
	*info = -1;
#line 200 "sgbtrf.f"
    } else if (*n < 0) {
#line 201 "sgbtrf.f"
	*info = -2;
#line 202 "sgbtrf.f"
    } else if (*kl < 0) {
#line 203 "sgbtrf.f"
	*info = -3;
#line 204 "sgbtrf.f"
    } else if (*ku < 0) {
#line 205 "sgbtrf.f"
	*info = -4;
#line 206 "sgbtrf.f"
    } else if (*ldab < *kl + kv + 1) {
#line 207 "sgbtrf.f"
	*info = -6;
#line 208 "sgbtrf.f"
    }
#line 209 "sgbtrf.f"
    if (*info != 0) {
#line 210 "sgbtrf.f"
	i__1 = -(*info);
#line 210 "sgbtrf.f"
	xerbla_("SGBTRF", &i__1, (ftnlen)6);
#line 211 "sgbtrf.f"
	return 0;
#line 212 "sgbtrf.f"
    }

/*     Quick return if possible */

#line 216 "sgbtrf.f"
    if (*m == 0 || *n == 0) {
#line 216 "sgbtrf.f"
	return 0;
#line 216 "sgbtrf.f"
    }

/*     Determine the block size for this environment */

#line 221 "sgbtrf.f"
    nb = ilaenv_(&c__1, "SGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);

/*     The block size must not exceed the limit set by the size of the */
/*     local arrays WORK13 and WORK31. */

#line 226 "sgbtrf.f"
    nb = min(nb,64);

#line 228 "sgbtrf.f"
    if (nb <= 1 || nb > *kl) {

/*        Use unblocked code */

#line 232 "sgbtrf.f"
	sgbtf2_(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
#line 233 "sgbtrf.f"
    } else {

/*        Use blocked code */

/*        Zero the superdiagonal elements of the work array WORK13 */

#line 239 "sgbtrf.f"
	i__1 = nb;
#line 239 "sgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 240 "sgbtrf.f"
	    i__2 = j - 1;
#line 240 "sgbtrf.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 241 "sgbtrf.f"
		work13[i__ + j * 65 - 66] = 0.;
#line 242 "sgbtrf.f"
/* L10: */
#line 242 "sgbtrf.f"
	    }
#line 243 "sgbtrf.f"
/* L20: */
#line 243 "sgbtrf.f"
	}

/*        Zero the subdiagonal elements of the work array WORK31 */

#line 247 "sgbtrf.f"
	i__1 = nb;
#line 247 "sgbtrf.f"
	for (j = 1; j <= i__1; ++j) {
#line 248 "sgbtrf.f"
	    i__2 = nb;
#line 248 "sgbtrf.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 249 "sgbtrf.f"
		work31[i__ + j * 65 - 66] = 0.;
#line 250 "sgbtrf.f"
/* L30: */
#line 250 "sgbtrf.f"
	    }
#line 251 "sgbtrf.f"
/* L40: */
#line 251 "sgbtrf.f"
	}

/*        Gaussian elimination with partial pivoting */

/*        Set fill-in elements in columns KU+2 to KV to zero */

#line 257 "sgbtrf.f"
	i__1 = min(kv,*n);
#line 257 "sgbtrf.f"
	for (j = *ku + 2; j <= i__1; ++j) {
#line 258 "sgbtrf.f"
	    i__2 = *kl;
#line 258 "sgbtrf.f"
	    for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
#line 259 "sgbtrf.f"
		ab[i__ + j * ab_dim1] = 0.;
#line 260 "sgbtrf.f"
/* L50: */
#line 260 "sgbtrf.f"
	    }
#line 261 "sgbtrf.f"
/* L60: */
#line 261 "sgbtrf.f"
	}

/*        JU is the index of the last column affected by the current */
/*        stage of the factorization */

#line 266 "sgbtrf.f"
	ju = 1;

#line 268 "sgbtrf.f"
	i__1 = min(*m,*n);
#line 268 "sgbtrf.f"
	i__2 = nb;
#line 268 "sgbtrf.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 269 "sgbtrf.f"
	    i__3 = nb, i__4 = min(*m,*n) - j + 1;
#line 269 "sgbtrf.f"
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
#line 283 "sgbtrf.f"
	    i__3 = *kl - jb, i__4 = *m - j - jb + 1;
#line 283 "sgbtrf.f"
	    i2 = min(i__3,i__4);
/* Computing MIN */
#line 284 "sgbtrf.f"
	    i__3 = jb, i__4 = *m - j - *kl + 1;
#line 284 "sgbtrf.f"
	    i3 = min(i__3,i__4);

/*           J2 and J3 are computed after JU has been updated. */

/*           Factorize the current block of JB columns */

#line 290 "sgbtrf.f"
	    i__3 = j + jb - 1;
#line 290 "sgbtrf.f"
	    for (jj = j; jj <= i__3; ++jj) {

/*              Set fill-in elements in column JJ+KV to zero */

#line 294 "sgbtrf.f"
		if (jj + kv <= *n) {
#line 295 "sgbtrf.f"
		    i__4 = *kl;
#line 295 "sgbtrf.f"
		    for (i__ = 1; i__ <= i__4; ++i__) {
#line 296 "sgbtrf.f"
			ab[i__ + (jj + kv) * ab_dim1] = 0.;
#line 297 "sgbtrf.f"
/* L70: */
#line 297 "sgbtrf.f"
		    }
#line 298 "sgbtrf.f"
		}

/*              Find pivot and test for singularity. KM is the number of */
/*              subdiagonal elements in the current column. */

/* Computing MIN */
#line 303 "sgbtrf.f"
		i__4 = *kl, i__5 = *m - jj;
#line 303 "sgbtrf.f"
		km = min(i__4,i__5);
#line 304 "sgbtrf.f"
		i__4 = km + 1;
#line 304 "sgbtrf.f"
		jp = isamax_(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
#line 305 "sgbtrf.f"
		ipiv[jj] = jp + jj - j;
#line 306 "sgbtrf.f"
		if (ab[kv + jp + jj * ab_dim1] != 0.) {
/* Computing MAX */
/* Computing MIN */
#line 307 "sgbtrf.f"
		    i__6 = jj + *ku + jp - 1;
#line 307 "sgbtrf.f"
		    i__4 = ju, i__5 = min(i__6,*n);
#line 307 "sgbtrf.f"
		    ju = max(i__4,i__5);
#line 308 "sgbtrf.f"
		    if (jp != 1) {

/*                    Apply interchange to columns J to J+JB-1 */

#line 312 "sgbtrf.f"
			if (jp + jj - 1 < j + *kl) {

#line 314 "sgbtrf.f"
			    i__4 = *ldab - 1;
#line 314 "sgbtrf.f"
			    i__5 = *ldab - 1;
#line 314 "sgbtrf.f"
			    sswap_(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
				    i__4, &ab[kv + jp + jj - j + j * ab_dim1],
				     &i__5);
#line 316 "sgbtrf.f"
			} else {

/*                       The interchange affects columns J to JJ-1 of A31 */
/*                       which are stored in the work array WORK31 */

#line 321 "sgbtrf.f"
			    i__4 = jj - j;
#line 321 "sgbtrf.f"
			    i__5 = *ldab - 1;
#line 321 "sgbtrf.f"
			    sswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], 
				    &i__5, &work31[jp + jj - j - *kl - 1], &
				    c__65);
#line 323 "sgbtrf.f"
			    i__4 = j + jb - jj;
#line 323 "sgbtrf.f"
			    i__5 = *ldab - 1;
#line 323 "sgbtrf.f"
			    i__6 = *ldab - 1;
#line 323 "sgbtrf.f"
			    sswap_(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
				    ab[kv + jp + jj * ab_dim1], &i__6);
#line 325 "sgbtrf.f"
			}
#line 326 "sgbtrf.f"
		    }

/*                 Compute multipliers */

#line 330 "sgbtrf.f"
		    d__1 = 1. / ab[kv + 1 + jj * ab_dim1];
#line 330 "sgbtrf.f"
		    sscal_(&km, &d__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

/*                 Update trailing submatrix within the band and within */
/*                 the current block. JM is the index of the last column */
/*                 which needs to be updated. */

/* Computing MIN */
#line 337 "sgbtrf.f"
		    i__4 = ju, i__5 = j + jb - 1;
#line 337 "sgbtrf.f"
		    jm = min(i__4,i__5);
#line 338 "sgbtrf.f"
		    if (jm > jj) {
#line 338 "sgbtrf.f"
			i__4 = jm - jj;
#line 338 "sgbtrf.f"
			i__5 = *ldab - 1;
#line 338 "sgbtrf.f"
			i__6 = *ldab - 1;
#line 338 "sgbtrf.f"
			sger_(&km, &i__4, &c_b18, &ab[kv + 2 + jj * ab_dim1], 
				&c__1, &ab[kv + (jj + 1) * ab_dim1], &i__5, &
				ab[kv + 1 + (jj + 1) * ab_dim1], &i__6);
#line 338 "sgbtrf.f"
		    }
#line 342 "sgbtrf.f"
		} else {

/*                 If pivot is zero, set INFO to the index of the pivot */
/*                 unless a zero pivot has already been found. */

#line 347 "sgbtrf.f"
		    if (*info == 0) {
#line 347 "sgbtrf.f"
			*info = jj;
#line 347 "sgbtrf.f"
		    }
#line 349 "sgbtrf.f"
		}

/*              Copy current column of A31 into the work array WORK31 */

/* Computing MIN */
#line 353 "sgbtrf.f"
		i__4 = jj - j + 1;
#line 353 "sgbtrf.f"
		nw = min(i__4,i3);
#line 354 "sgbtrf.f"
		if (nw > 0) {
#line 354 "sgbtrf.f"
		    scopy_(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
			    c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
#line 354 "sgbtrf.f"
		}
#line 357 "sgbtrf.f"
/* L80: */
#line 357 "sgbtrf.f"
	    }
#line 358 "sgbtrf.f"
	    if (j + jb <= *n) {

/*              Apply the row interchanges to the other blocks. */

/* Computing MIN */
#line 362 "sgbtrf.f"
		i__3 = ju - j + 1;
#line 362 "sgbtrf.f"
		j2 = min(i__3,kv) - jb;
/* Computing MAX */
#line 363 "sgbtrf.f"
		i__3 = 0, i__4 = ju - j - kv + 1;
#line 363 "sgbtrf.f"
		j3 = max(i__3,i__4);

/*              Use SLASWP to apply the row interchanges to A12, A22, and */
/*              A32. */

#line 368 "sgbtrf.f"
		i__3 = *ldab - 1;
#line 368 "sgbtrf.f"
		slaswp_(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
			c__1, &jb, &ipiv[j], &c__1);

/*              Adjust the pivot indices. */

#line 373 "sgbtrf.f"
		i__3 = j + jb - 1;
#line 373 "sgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 374 "sgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 375 "sgbtrf.f"
/* L90: */
#line 375 "sgbtrf.f"
		}

/*              Apply the row interchanges to A13, A23, and A33 */
/*              columnwise. */

#line 380 "sgbtrf.f"
		k2 = j - 1 + jb + j2;
#line 381 "sgbtrf.f"
		i__3 = j3;
#line 381 "sgbtrf.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 382 "sgbtrf.f"
		    jj = k2 + i__;
#line 383 "sgbtrf.f"
		    i__4 = j + jb - 1;
#line 383 "sgbtrf.f"
		    for (ii = j + i__ - 1; ii <= i__4; ++ii) {
#line 384 "sgbtrf.f"
			ip = ipiv[ii];
#line 385 "sgbtrf.f"
			if (ip != ii) {
#line 386 "sgbtrf.f"
			    temp = ab[kv + 1 + ii - jj + jj * ab_dim1];
#line 387 "sgbtrf.f"
			    ab[kv + 1 + ii - jj + jj * ab_dim1] = ab[kv + 1 + 
				    ip - jj + jj * ab_dim1];
#line 388 "sgbtrf.f"
			    ab[kv + 1 + ip - jj + jj * ab_dim1] = temp;
#line 389 "sgbtrf.f"
			}
#line 390 "sgbtrf.f"
/* L100: */
#line 390 "sgbtrf.f"
		    }
#line 391 "sgbtrf.f"
/* L110: */
#line 391 "sgbtrf.f"
		}

/*              Update the relevant part of the trailing submatrix */

#line 395 "sgbtrf.f"
		if (j2 > 0) {

/*                 Update A12 */

#line 399 "sgbtrf.f"
		    i__3 = *ldab - 1;
#line 399 "sgbtrf.f"
		    i__4 = *ldab - 1;
#line 399 "sgbtrf.f"
		    strsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
			    &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv 
			    + 1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4, 
			    (ftnlen)5, (ftnlen)12, (ftnlen)4);

#line 403 "sgbtrf.f"
		    if (i2 > 0) {

/*                    Update A22 */

#line 407 "sgbtrf.f"
			i__3 = *ldab - 1;
#line 407 "sgbtrf.f"
			i__4 = *ldab - 1;
#line 407 "sgbtrf.f"
			i__5 = *ldab - 1;
#line 407 "sgbtrf.f"
			sgemm_("No transpose", "No transpose", &i2, &j2, &jb, 
				&c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
				 &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4,
				 &c_b31, &ab[kv + 1 + (j + jb) * ab_dim1], &
				i__5, (ftnlen)12, (ftnlen)12);
#line 411 "sgbtrf.f"
		    }

#line 413 "sgbtrf.f"
		    if (i3 > 0) {

/*                    Update A32 */

#line 417 "sgbtrf.f"
			i__3 = *ldab - 1;
#line 417 "sgbtrf.f"
			i__4 = *ldab - 1;
#line 417 "sgbtrf.f"
			sgemm_("No transpose", "No transpose", &i3, &j2, &jb, 
				&c_b18, work31, &c__65, &ab[kv + 1 - jb + (j 
				+ jb) * ab_dim1], &i__3, &c_b31, &ab[kv + *kl 
				+ 1 - jb + (j + jb) * ab_dim1], &i__4, (
				ftnlen)12, (ftnlen)12);
#line 421 "sgbtrf.f"
		    }
#line 422 "sgbtrf.f"
		}

#line 424 "sgbtrf.f"
		if (j3 > 0) {

/*                 Copy the lower triangle of A13 into the work array */
/*                 WORK13 */

#line 429 "sgbtrf.f"
		    i__3 = j3;
#line 429 "sgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 430 "sgbtrf.f"
			i__4 = jb;
#line 430 "sgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 431 "sgbtrf.f"
			    work13[ii + jj * 65 - 66] = ab[ii - jj + 1 + (jj 
				    + j + kv - 1) * ab_dim1];
#line 432 "sgbtrf.f"
/* L120: */
#line 432 "sgbtrf.f"
			}
#line 433 "sgbtrf.f"
/* L130: */
#line 433 "sgbtrf.f"
		    }

/*                 Update A13 in the work array */

#line 437 "sgbtrf.f"
		    i__3 = *ldab - 1;
#line 437 "sgbtrf.f"
		    strsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
			    &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, work13, 
			    &c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
			    4);

#line 441 "sgbtrf.f"
		    if (i2 > 0) {

/*                    Update A23 */

#line 445 "sgbtrf.f"
			i__3 = *ldab - 1;
#line 445 "sgbtrf.f"
			i__4 = *ldab - 1;
#line 445 "sgbtrf.f"
			sgemm_("No transpose", "No transpose", &i2, &j3, &jb, 
				&c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
				 work13, &c__65, &c_b31, &ab[jb + 1 + (j + kv)
				 * ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
#line 449 "sgbtrf.f"
		    }

#line 451 "sgbtrf.f"
		    if (i3 > 0) {

/*                    Update A33 */

#line 455 "sgbtrf.f"
			i__3 = *ldab - 1;
#line 455 "sgbtrf.f"
			sgemm_("No transpose", "No transpose", &i3, &j3, &jb, 
				&c_b18, work31, &c__65, work13, &c__65, &
				c_b31, &ab[*kl + 1 + (j + kv) * ab_dim1], &
				i__3, (ftnlen)12, (ftnlen)12);
#line 458 "sgbtrf.f"
		    }

/*                 Copy the lower triangle of A13 back into place */

#line 462 "sgbtrf.f"
		    i__3 = j3;
#line 462 "sgbtrf.f"
		    for (jj = 1; jj <= i__3; ++jj) {
#line 463 "sgbtrf.f"
			i__4 = jb;
#line 463 "sgbtrf.f"
			for (ii = jj; ii <= i__4; ++ii) {
#line 464 "sgbtrf.f"
			    ab[ii - jj + 1 + (jj + j + kv - 1) * ab_dim1] = 
				    work13[ii + jj * 65 - 66];
#line 465 "sgbtrf.f"
/* L140: */
#line 465 "sgbtrf.f"
			}
#line 466 "sgbtrf.f"
/* L150: */
#line 466 "sgbtrf.f"
		    }
#line 467 "sgbtrf.f"
		}
#line 468 "sgbtrf.f"
	    } else {

/*              Adjust the pivot indices. */

#line 472 "sgbtrf.f"
		i__3 = j + jb - 1;
#line 472 "sgbtrf.f"
		for (i__ = j; i__ <= i__3; ++i__) {
#line 473 "sgbtrf.f"
		    ipiv[i__] = ipiv[i__] + j - 1;
#line 474 "sgbtrf.f"
/* L160: */
#line 474 "sgbtrf.f"
		}
#line 475 "sgbtrf.f"
	    }

/*           Partially undo the interchanges in the current block to */
/*           restore the upper triangular form of A31 and copy the upper */
/*           triangle of A31 back into place */

#line 481 "sgbtrf.f"
	    i__3 = j;
#line 481 "sgbtrf.f"
	    for (jj = j + jb - 1; jj >= i__3; --jj) {
#line 482 "sgbtrf.f"
		jp = ipiv[jj] - jj + 1;
#line 483 "sgbtrf.f"
		if (jp != 1) {

/*                 Apply interchange to columns J to JJ-1 */

#line 487 "sgbtrf.f"
		    if (jp + jj - 1 < j + *kl) {

/*                    The interchange does not affect A31 */

#line 491 "sgbtrf.f"
			i__4 = jj - j;
#line 491 "sgbtrf.f"
			i__5 = *ldab - 1;
#line 491 "sgbtrf.f"
			i__6 = *ldab - 1;
#line 491 "sgbtrf.f"
			sswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
				i__6);
#line 493 "sgbtrf.f"
		    } else {

/*                    The interchange does affect A31 */

#line 497 "sgbtrf.f"
			i__4 = jj - j;
#line 497 "sgbtrf.f"
			i__5 = *ldab - 1;
#line 497 "sgbtrf.f"
			sswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
				i__5, &work31[jp + jj - j - *kl - 1], &c__65);
#line 499 "sgbtrf.f"
		    }
#line 500 "sgbtrf.f"
		}

/*              Copy the current column of A31 back into place */

/* Computing MIN */
#line 504 "sgbtrf.f"
		i__4 = i3, i__5 = jj - j + 1;
#line 504 "sgbtrf.f"
		nw = min(i__4,i__5);
#line 505 "sgbtrf.f"
		if (nw > 0) {
#line 505 "sgbtrf.f"
		    scopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
			    kv + *kl + 1 - jj + j + jj * ab_dim1], &c__1);
#line 505 "sgbtrf.f"
		}
#line 508 "sgbtrf.f"
/* L170: */
#line 508 "sgbtrf.f"
	    }
#line 509 "sgbtrf.f"
/* L180: */
#line 509 "sgbtrf.f"
	}
#line 510 "sgbtrf.f"
    }

#line 512 "sgbtrf.f"
    return 0;

/*     End of SGBTRF */

} /* sgbtrf_ */

