#line 1 "ssytrf_aa_2stage.f"
/* ssytrf_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssytrf_aa_2stage.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b12 = 1.;
static doublereal c_b13 = 0.;
static doublereal c_b21 = -1.;

/* > \brief \b SSYTRF_AA_2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRF_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrf_
aa_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrf_
aa_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrf_
aa_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE SSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
/*                                   IPIV2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LTB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IPIV2( * ) */
/*       REAL               A( LDA, * ), TB( * ), WORK( * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRF_AA_2STAGE computes the factorization of a real symmetric matrix A */
/* > using the Aasen's algorithm.  The form of the factorization is */
/* > */
/* >    A = U*T*U**T  or  A = L*T*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a symmetric band matrix with the */
/* > bandwidth of NB (NB is internally selected and stored in TB( 1 ), and T is */
/* > LU factorized with partial pivoting). */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, L is stored below (or above) the subdiaonal blocks, */
/* >          when UPLO  is 'L' (or 'U'). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TB */
/* > \verbatim */
/* >          TB is REAL array, dimension (LTB) */
/* >          On exit, details of the LU factorization of the band matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LTB */
/* > \verbatim */
/* >          The size of the array TB. LTB >= 4*N, internally */
/* >          used to select NB such that LTB >= (3*NB+1)*N. */
/* > */
/* >          If LTB = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of LTB, */
/* >          returns this value as the first entry of TB, and */
/* >          no error message related to LTB is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On exit, it contains the details of the interchanges, i.e., */
/* >          the row and column k of A were interchanged with the */
/* >          row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV2 */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On exit, it contains the details of the interchanges, i.e., */
/* >          the row and column k of T were interchanged with the */
/* >          row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL workspace of size LWORK */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          The size of WORK. LWORK >= N, internally used to select NB */
/* >          such that LWORK >= N*NB. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of the WORK array, */
/* >          returns this value as the first entry of the WORK array, and */
/* >          no error message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, band LU factorization failed on i-th column */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssytrf_aa_2stage__(char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *tb, integer *ltb, integer *ipiv, integer *
	ipiv2, doublereal *work, integer *lwork, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, i1, i2, jb, kb, nb, td, nt;
    static doublereal piv;
    static integer ldtb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), strsm_(char *, char *, char *, char *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgbtrf_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    sgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), slacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), slaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    static logical tquery, wquery;
    extern /* Subroutine */ int ssygst_(integer *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*     .. Parameters .. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 205 "ssytrf_aa_2stage.f"
    /* Parameter adjustments */
#line 205 "ssytrf_aa_2stage.f"
    a_dim1 = *lda;
#line 205 "ssytrf_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 205 "ssytrf_aa_2stage.f"
    a -= a_offset;
#line 205 "ssytrf_aa_2stage.f"
    --tb;
#line 205 "ssytrf_aa_2stage.f"
    --ipiv;
#line 205 "ssytrf_aa_2stage.f"
    --ipiv2;
#line 205 "ssytrf_aa_2stage.f"
    --work;
#line 205 "ssytrf_aa_2stage.f"

#line 205 "ssytrf_aa_2stage.f"
    /* Function Body */
#line 205 "ssytrf_aa_2stage.f"
    *info = 0;
#line 206 "ssytrf_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "ssytrf_aa_2stage.f"
    wquery = *lwork == -1;
#line 208 "ssytrf_aa_2stage.f"
    tquery = *ltb == -1;
#line 209 "ssytrf_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 210 "ssytrf_aa_2stage.f"
	*info = -1;
#line 211 "ssytrf_aa_2stage.f"
    } else if (*n < 0) {
#line 212 "ssytrf_aa_2stage.f"
	*info = -2;
#line 213 "ssytrf_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 214 "ssytrf_aa_2stage.f"
	*info = -4;
#line 215 "ssytrf_aa_2stage.f"
    } else if (*ltb < *n << 2 && ! tquery) {
#line 216 "ssytrf_aa_2stage.f"
	*info = -6;
#line 217 "ssytrf_aa_2stage.f"
    } else if (*lwork < *n && ! wquery) {
#line 218 "ssytrf_aa_2stage.f"
	*info = -10;
#line 219 "ssytrf_aa_2stage.f"
    }

#line 221 "ssytrf_aa_2stage.f"
    if (*info != 0) {
#line 222 "ssytrf_aa_2stage.f"
	i__1 = -(*info);
#line 222 "ssytrf_aa_2stage.f"
	xerbla_("SSYTRF_AA_2STAGE", &i__1, (ftnlen)16);
#line 223 "ssytrf_aa_2stage.f"
	return 0;
#line 224 "ssytrf_aa_2stage.f"
    }

/*     Answer the query */

#line 228 "ssytrf_aa_2stage.f"
    nb = ilaenv_(&c__1, "SSYTRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)16, (ftnlen)1);
#line 229 "ssytrf_aa_2stage.f"
    if (*info == 0) {
#line 230 "ssytrf_aa_2stage.f"
	if (tquery) {
#line 231 "ssytrf_aa_2stage.f"
	    tb[1] = (doublereal) ((nb * 3 + 1) * *n);
#line 232 "ssytrf_aa_2stage.f"
	}
#line 233 "ssytrf_aa_2stage.f"
	if (wquery) {
#line 234 "ssytrf_aa_2stage.f"
	    work[1] = (doublereal) (*n * nb);
#line 235 "ssytrf_aa_2stage.f"
	}
#line 236 "ssytrf_aa_2stage.f"
    }
#line 237 "ssytrf_aa_2stage.f"
    if (tquery || wquery) {
#line 238 "ssytrf_aa_2stage.f"
	return 0;
#line 239 "ssytrf_aa_2stage.f"
    }

/*     Quick return */

#line 243 "ssytrf_aa_2stage.f"
    if (*n == 0) {
#line 244 "ssytrf_aa_2stage.f"
	return 0;
#line 245 "ssytrf_aa_2stage.f"
    }

/*     Determine the number of the block size */

#line 249 "ssytrf_aa_2stage.f"
    ldtb = *ltb / *n;
#line 250 "ssytrf_aa_2stage.f"
    if (ldtb < nb * 3 + 1) {
#line 251 "ssytrf_aa_2stage.f"
	nb = (ldtb - 1) / 3;
#line 252 "ssytrf_aa_2stage.f"
    }
#line 253 "ssytrf_aa_2stage.f"
    if (*lwork < nb * *n) {
#line 254 "ssytrf_aa_2stage.f"
	nb = *lwork / *n;
#line 255 "ssytrf_aa_2stage.f"
    }

/*     Determine the number of the block columns */

#line 259 "ssytrf_aa_2stage.f"
    nt = (*n + nb - 1) / nb;
#line 260 "ssytrf_aa_2stage.f"
    td = nb << 1;
#line 261 "ssytrf_aa_2stage.f"
    kb = min(nb,*n);

/*     Initialize vectors/matrices */

#line 265 "ssytrf_aa_2stage.f"
    i__1 = kb;
#line 265 "ssytrf_aa_2stage.f"
    for (j = 1; j <= i__1; ++j) {
#line 266 "ssytrf_aa_2stage.f"
	ipiv[j] = j;
#line 267 "ssytrf_aa_2stage.f"
    }

/*     Save NB */

#line 271 "ssytrf_aa_2stage.f"
    tb[1] = (doublereal) nb;

#line 273 "ssytrf_aa_2stage.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

#line 279 "ssytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 279 "ssytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 283 "ssytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 283 "ssytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 284 "ssytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 284 "ssytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 285 "ssytrf_aa_2stage.f"
		if (i__ == 1) {
/*                 H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J) */
#line 287 "ssytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 288 "ssytrf_aa_2stage.f"
			jb = nb + kb;
#line 289 "ssytrf_aa_2stage.f"
		    } else {
#line 290 "ssytrf_aa_2stage.f"
			jb = nb << 1;
#line 291 "ssytrf_aa_2stage.f"
		    }
#line 292 "ssytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 292 "ssytrf_aa_2stage.f"
		    sgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &
			    c_b12, &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(
			    i__ - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &
			    c_b13, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)11);
#line 297 "ssytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 299 "ssytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 300 "ssytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 301 "ssytrf_aa_2stage.f"
		    } else {
#line 302 "ssytrf_aa_2stage.f"
			jb = nb * 3;
#line 303 "ssytrf_aa_2stage.f"
		    }
#line 304 "ssytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 304 "ssytrf_aa_2stage.f"
		    sgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &
			    c_b12, &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &
			    i__3, &a[(i__ - 2) * nb + 1 + (j * nb + 1) * 
			    a_dim1], lda, &c_b13, &work[i__ * nb + 1], n, (
			    ftnlen)11, (ftnlen)11);
#line 310 "ssytrf_aa_2stage.f"
		}
#line 311 "ssytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 315 "ssytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 315 "ssytrf_aa_2stage.f"
	    slacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 317 "ssytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = U(1:J,J)'*H(1:J) */
#line 319 "ssytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 319 "ssytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 319 "ssytrf_aa_2stage.f"
		sgemm_("Transpose", "NoTranspose", &kb, &kb, &i__2, &c_b21, &
			a[(j * nb + 1) * a_dim1 + 1], lda, &work[nb + 1], n, &
			c_b12, &tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)9, 
			(ftnlen)11);
/*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
#line 325 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 325 "ssytrf_aa_2stage.f"
		sgemm_("Transpose", "NoTranspose", &kb, &nb, &kb, &c_b12, &a[(
			j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td 
			+ nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b13, &work[
			1], n, (ftnlen)9, (ftnlen)11);
#line 330 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 330 "ssytrf_aa_2stage.f"
		sgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &c_b21, &
			work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &c_b12, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)11);
#line 335 "ssytrf_aa_2stage.f"
	    }
#line 336 "ssytrf_aa_2stage.f"
	    if (j > 0) {
#line 337 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 337 "ssytrf_aa_2stage.f"
		ssygst_(&c__1, "Upper", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 340 "ssytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 344 "ssytrf_aa_2stage.f"
	    i__2 = kb;
#line 344 "ssytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 345 "ssytrf_aa_2stage.f"
		i__3 = kb;
#line 345 "ssytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 346 "ssytrf_aa_2stage.f"
		    tb[td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb] = tb[
			    td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb];
#line 348 "ssytrf_aa_2stage.f"
		}
#line 349 "ssytrf_aa_2stage.f"
	    }

#line 351 "ssytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 352 "ssytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 356 "ssytrf_aa_2stage.f"
		    if (j == 1) {
#line 357 "ssytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 357 "ssytrf_aa_2stage.f"
			sgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, &
				c_b12, &tb[td + 1 + j * nb * ldtb], &i__2, &a[
				(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
				lda, &c_b13, &work[j * nb + 1], n, (ftnlen)11,
				 (ftnlen)11);
#line 362 "ssytrf_aa_2stage.f"
		    } else {
#line 363 "ssytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 363 "ssytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 363 "ssytrf_aa_2stage.f"
			sgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, 
				&c_b12, &tb[td + nb + 1 + (j - 1) * nb * ldtb]
				, &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
				a_dim1], lda, &c_b13, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)11);
#line 369 "ssytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 373 "ssytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 373 "ssytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 373 "ssytrf_aa_2stage.f"
		    sgemm_("Transpose", "NoTranspose", &nb, &i__2, &i__3, &
			    c_b21, &work[nb + 1], n, &a[((j + 1) * nb + 1) * 
			    a_dim1 + 1], lda, &c_b12, &a[j * nb + 1 + ((j + 1)
			     * nb + 1) * a_dim1], lda, (ftnlen)9, (ftnlen)11);
#line 378 "ssytrf_aa_2stage.f"
		}

/*              Copy panel to workspace to call SGETRF */

#line 382 "ssytrf_aa_2stage.f"
		i__2 = nb;
#line 382 "ssytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 383 "ssytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 383 "ssytrf_aa_2stage.f"
		    scopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1]
			    , lda, &work[(k - 1) * *n + 1], &c__1);
#line 386 "ssytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 390 "ssytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 390 "ssytrf_aa_2stage.f"
		sgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], &
			iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Copy panel back */

#line 399 "ssytrf_aa_2stage.f"
		i__2 = nb;
#line 399 "ssytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 400 "ssytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 400 "ssytrf_aa_2stage.f"
		    scopy_(&i__3, &work[(k - 1) * *n + 1], &c__1, &a[j * nb + 
			    k + ((j + 1) * nb + 1) * a_dim1], lda);
#line 403 "ssytrf_aa_2stage.f"
		}

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 407 "ssytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 407 "ssytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 408 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 408 "ssytrf_aa_2stage.f"
		slaset_("Full", &kb, &nb, &c_b13, &c_b13, &tb[td + nb + 1 + j 
			* nb * ldtb], &i__2, (ftnlen)4);
#line 410 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 410 "ssytrf_aa_2stage.f"
		slacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)5);
#line 413 "ssytrf_aa_2stage.f"
		if (j > 0) {
#line 414 "ssytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 414 "ssytrf_aa_2stage.f"
		    strsm_("R", "U", "N", "U", &kb, &nb, &c_b12, &a[(j - 1) * 
			    nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 417 "ssytrf_aa_2stage.f"
		}

/*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
/*              updates */

#line 422 "ssytrf_aa_2stage.f"
		i__2 = nb;
#line 422 "ssytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 423 "ssytrf_aa_2stage.f"
		    i__3 = kb;
#line 423 "ssytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 424 "ssytrf_aa_2stage.f"
			tb[td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1) * 
				ldtb] = tb[td + nb + i__ - k + 1 + (j * nb + 
				k - 1) * ldtb];
#line 426 "ssytrf_aa_2stage.f"
		    }
#line 427 "ssytrf_aa_2stage.f"
		}
#line 428 "ssytrf_aa_2stage.f"
		slaset_("Lower", &kb, &nb, &c_b13, &c_b12, &a[j * nb + 1 + ((
			j + 1) * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 433 "ssytrf_aa_2stage.f"
		i__2 = kb;
#line 433 "ssytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 435 "ssytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 437 "ssytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 438 "ssytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 439 "ssytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 441 "ssytrf_aa_2stage.f"
			i__3 = k - 1;
#line 441 "ssytrf_aa_2stage.f"
			sswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], &
				c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], &
				c__1);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 444 "ssytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 444 "ssytrf_aa_2stage.f"
			sswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[i1 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 447 "ssytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 447 "ssytrf_aa_2stage.f"
			sswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[i2 
				+ (i2 + 1) * a_dim1], lda);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 450 "ssytrf_aa_2stage.f"
			piv = a[i1 + i1 * a_dim1];
#line 451 "ssytrf_aa_2stage.f"
			a[i1 + i1 * a_dim1] = a[i2 + i2 * a_dim1];
#line 452 "ssytrf_aa_2stage.f"
			a[i2 + i2 * a_dim1] = piv;
/*                    > Apply pivots to previous columns of L */
#line 454 "ssytrf_aa_2stage.f"
			if (j > 0) {
#line 455 "ssytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 455 "ssytrf_aa_2stage.f"
			    sswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * 
				    a_dim1 + 1], &c__1);
#line 457 "ssytrf_aa_2stage.f"
			}
#line 458 "ssytrf_aa_2stage.f"
		    }
#line 459 "ssytrf_aa_2stage.f"
		}
#line 460 "ssytrf_aa_2stage.f"
	    }
#line 461 "ssytrf_aa_2stage.f"
	}
#line 462 "ssytrf_aa_2stage.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 468 "ssytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 468 "ssytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 472 "ssytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 472 "ssytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 473 "ssytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 473 "ssytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 474 "ssytrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
#line 476 "ssytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 477 "ssytrf_aa_2stage.f"
			jb = nb + kb;
#line 478 "ssytrf_aa_2stage.f"
		    } else {
#line 479 "ssytrf_aa_2stage.f"
			jb = nb << 1;
#line 480 "ssytrf_aa_2stage.f"
		    }
#line 481 "ssytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 481 "ssytrf_aa_2stage.f"
		    sgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b12, 
			    &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[j * nb + 
			    1 + ((i__ - 1) * nb + 1) * a_dim1], lda, &c_b13, &
			    work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)9);
#line 486 "ssytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
#line 488 "ssytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 489 "ssytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 490 "ssytrf_aa_2stage.f"
		    } else {
#line 491 "ssytrf_aa_2stage.f"
			jb = nb * 3;
#line 492 "ssytrf_aa_2stage.f"
		    }
#line 493 "ssytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 493 "ssytrf_aa_2stage.f"
		    sgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b12, 
			    &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, &
			    a[j * nb + 1 + ((i__ - 2) * nb + 1) * a_dim1], 
			    lda, &c_b13, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)9);
#line 499 "ssytrf_aa_2stage.f"
		}
#line 500 "ssytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 504 "ssytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 504 "ssytrf_aa_2stage.f"
	    slacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 506 "ssytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = L(J,1:J)*H(1:J) */
#line 508 "ssytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 508 "ssytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 508 "ssytrf_aa_2stage.f"
		sgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &c_b21, 
			&a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &
			c_b12, &tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)11,
			 (ftnlen)11);
/*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
#line 514 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 514 "ssytrf_aa_2stage.f"
		sgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b12, &
			a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[
			td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b13, &
			work[1], n, (ftnlen)11, (ftnlen)11);
#line 519 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 519 "ssytrf_aa_2stage.f"
		sgemm_("NoTranspose", "Transpose", &kb, &kb, &nb, &c_b21, &
			work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
			a_dim1], lda, &c_b12, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)9);
#line 524 "ssytrf_aa_2stage.f"
	    }
#line 525 "ssytrf_aa_2stage.f"
	    if (j > 0) {
#line 526 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 526 "ssytrf_aa_2stage.f"
		ssygst_(&c__1, "Lower", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 529 "ssytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 533 "ssytrf_aa_2stage.f"
	    i__2 = kb;
#line 533 "ssytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 534 "ssytrf_aa_2stage.f"
		i__3 = kb;
#line 534 "ssytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 535 "ssytrf_aa_2stage.f"
		    tb[td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb] = tb[
			    td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb];
#line 537 "ssytrf_aa_2stage.f"
		}
#line 538 "ssytrf_aa_2stage.f"
	    }

#line 540 "ssytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 541 "ssytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 545 "ssytrf_aa_2stage.f"
		    if (j == 1) {
#line 546 "ssytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 546 "ssytrf_aa_2stage.f"
			sgemm_("NoTranspose", "Transpose", &kb, &kb, &kb, &
				c_b12, &tb[td + 1 + j * nb * ldtb], &i__2, &a[
				j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
				lda, &c_b13, &work[j * nb + 1], n, (ftnlen)11,
				 (ftnlen)9);
#line 551 "ssytrf_aa_2stage.f"
		    } else {
#line 552 "ssytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 552 "ssytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 552 "ssytrf_aa_2stage.f"
			sgemm_("NoTranspose", "Transpose", &kb, &kb, &i__2, &
				c_b12, &tb[td + nb + 1 + (j - 1) * nb * ldtb],
				 &i__3, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
				a_dim1], lda, &c_b13, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)9);
#line 558 "ssytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 562 "ssytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 562 "ssytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 562 "ssytrf_aa_2stage.f"
		    sgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, &
			    c_b21, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[
			    nb + 1], n, &c_b12, &a[(j + 1) * nb + 1 + (j * nb 
			    + 1) * a_dim1], lda, (ftnlen)11, (ftnlen)11);
#line 567 "ssytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 571 "ssytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 571 "ssytrf_aa_2stage.f"
		sgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 580 "ssytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 580 "ssytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 581 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 581 "ssytrf_aa_2stage.f"
		slaset_("Full", &kb, &nb, &c_b13, &c_b13, &tb[td + nb + 1 + j 
			* nb * ldtb], &i__2, (ftnlen)4);
#line 583 "ssytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 583 "ssytrf_aa_2stage.f"
		slacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) 
			* a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &
			i__2, (ftnlen)5);
#line 586 "ssytrf_aa_2stage.f"
		if (j > 0) {
#line 587 "ssytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 587 "ssytrf_aa_2stage.f"
		    strsm_("R", "L", "T", "U", &kb, &nb, &c_b12, &a[j * nb + 
			    1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + 
			    nb + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 590 "ssytrf_aa_2stage.f"
		}

/*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
/*              updates */

#line 595 "ssytrf_aa_2stage.f"
		i__2 = nb;
#line 595 "ssytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 596 "ssytrf_aa_2stage.f"
		    i__3 = kb;
#line 596 "ssytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 597 "ssytrf_aa_2stage.f"
			tb[td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1) * 
				ldtb] = tb[td + nb + i__ - k + 1 + (j * nb + 
				k - 1) * ldtb];
#line 599 "ssytrf_aa_2stage.f"
		    }
#line 600 "ssytrf_aa_2stage.f"
		}
#line 601 "ssytrf_aa_2stage.f"
		slaset_("Upper", &kb, &nb, &c_b13, &c_b12, &a[(j + 1) * nb + 
			1 + (j * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 606 "ssytrf_aa_2stage.f"
		i__2 = kb;
#line 606 "ssytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 608 "ssytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 610 "ssytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 611 "ssytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 612 "ssytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 614 "ssytrf_aa_2stage.f"
			i__3 = k - 1;
#line 614 "ssytrf_aa_2stage.f"
			sswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], 
				lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], 
				lda);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 617 "ssytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 617 "ssytrf_aa_2stage.f"
			sswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ (i1 + 1) * a_dim1], lda);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 620 "ssytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 620 "ssytrf_aa_2stage.f"
			sswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 623 "ssytrf_aa_2stage.f"
			piv = a[i1 + i1 * a_dim1];
#line 624 "ssytrf_aa_2stage.f"
			a[i1 + i1 * a_dim1] = a[i2 + i2 * a_dim1];
#line 625 "ssytrf_aa_2stage.f"
			a[i2 + i2 * a_dim1] = piv;
/*                    > Apply pivots to previous columns of L */
#line 627 "ssytrf_aa_2stage.f"
			if (j > 0) {
#line 628 "ssytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 628 "ssytrf_aa_2stage.f"
			    sswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + 
				    a_dim1], lda);
#line 630 "ssytrf_aa_2stage.f"
			}
#line 631 "ssytrf_aa_2stage.f"
		    }
#line 632 "ssytrf_aa_2stage.f"
		}

/*              Apply pivots to previous columns of L */

/*               CALL SLASWP( J*NB, A( 1, 1 ), LDA, */
/*     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
#line 638 "ssytrf_aa_2stage.f"
	    }
#line 639 "ssytrf_aa_2stage.f"
	}
#line 640 "ssytrf_aa_2stage.f"
    }

/*     Factor the band matrix */
#line 643 "ssytrf_aa_2stage.f"
    sgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);

/*     End of SSYTRF_AA_2STAGE */

#line 647 "ssytrf_aa_2stage.f"
    return 0;
} /* ssytrf_aa_2stage__ */

