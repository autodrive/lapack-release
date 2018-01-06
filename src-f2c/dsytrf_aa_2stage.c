#line 1 "dsytrf_aa_2stage.f"
/* dsytrf_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "dsytrf_aa_2stage.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b12 = 1.;
static doublereal c_b13 = 0.;
static doublereal c_b21 = -1.;

/* > \brief \b DSYTRF_AA_2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRF_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrf_
aa_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrf_
aa_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrf_
aa_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE DSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
/*                                   IPIV2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LTB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IPIV2( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), TB( * ), WORK( * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRF_AA_2STAGE computes the factorization of a real symmetric matrix A */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          TB is DOUBLE PRECISION array, dimension (LTB) */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION workspace of size LWORK */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsytrf_aa_2stage__(char *uplo, integer *n, doublereal *a,
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
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), dtrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dgbtrf_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dsygst_(integer *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static logical tquery, wquery;


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

#line 205 "dsytrf_aa_2stage.f"
    /* Parameter adjustments */
#line 205 "dsytrf_aa_2stage.f"
    a_dim1 = *lda;
#line 205 "dsytrf_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 205 "dsytrf_aa_2stage.f"
    a -= a_offset;
#line 205 "dsytrf_aa_2stage.f"
    --tb;
#line 205 "dsytrf_aa_2stage.f"
    --ipiv;
#line 205 "dsytrf_aa_2stage.f"
    --ipiv2;
#line 205 "dsytrf_aa_2stage.f"
    --work;
#line 205 "dsytrf_aa_2stage.f"

#line 205 "dsytrf_aa_2stage.f"
    /* Function Body */
#line 205 "dsytrf_aa_2stage.f"
    *info = 0;
#line 206 "dsytrf_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "dsytrf_aa_2stage.f"
    wquery = *lwork == -1;
#line 208 "dsytrf_aa_2stage.f"
    tquery = *ltb == -1;
#line 209 "dsytrf_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 210 "dsytrf_aa_2stage.f"
	*info = -1;
#line 211 "dsytrf_aa_2stage.f"
    } else if (*n < 0) {
#line 212 "dsytrf_aa_2stage.f"
	*info = -2;
#line 213 "dsytrf_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 214 "dsytrf_aa_2stage.f"
	*info = -4;
#line 215 "dsytrf_aa_2stage.f"
    } else if (*ltb < *n << 2 && ! tquery) {
#line 216 "dsytrf_aa_2stage.f"
	*info = -6;
#line 217 "dsytrf_aa_2stage.f"
    } else if (*lwork < *n && ! wquery) {
#line 218 "dsytrf_aa_2stage.f"
	*info = -10;
#line 219 "dsytrf_aa_2stage.f"
    }

#line 221 "dsytrf_aa_2stage.f"
    if (*info != 0) {
#line 222 "dsytrf_aa_2stage.f"
	i__1 = -(*info);
#line 222 "dsytrf_aa_2stage.f"
	xerbla_("DSYTRF_AA_2STAGE", &i__1, (ftnlen)16);
#line 223 "dsytrf_aa_2stage.f"
	return 0;
#line 224 "dsytrf_aa_2stage.f"
    }

/*     Answer the query */

#line 228 "dsytrf_aa_2stage.f"
    nb = ilaenv_(&c__1, "DSYTRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)16, (ftnlen)1);
#line 229 "dsytrf_aa_2stage.f"
    if (*info == 0) {
#line 230 "dsytrf_aa_2stage.f"
	if (tquery) {
#line 231 "dsytrf_aa_2stage.f"
	    tb[1] = (doublereal) ((nb * 3 + 1) * *n);
#line 232 "dsytrf_aa_2stage.f"
	}
#line 233 "dsytrf_aa_2stage.f"
	if (wquery) {
#line 234 "dsytrf_aa_2stage.f"
	    work[1] = (doublereal) (*n * nb);
#line 235 "dsytrf_aa_2stage.f"
	}
#line 236 "dsytrf_aa_2stage.f"
    }
#line 237 "dsytrf_aa_2stage.f"
    if (tquery || wquery) {
#line 238 "dsytrf_aa_2stage.f"
	return 0;
#line 239 "dsytrf_aa_2stage.f"
    }

/*     Quick return */

#line 243 "dsytrf_aa_2stage.f"
    if (*n == 0) {
#line 244 "dsytrf_aa_2stage.f"
	return 0;
#line 245 "dsytrf_aa_2stage.f"
    }

/*     Determine the number of the block size */

#line 249 "dsytrf_aa_2stage.f"
    ldtb = *ltb / *n;
#line 250 "dsytrf_aa_2stage.f"
    if (ldtb < nb * 3 + 1) {
#line 251 "dsytrf_aa_2stage.f"
	nb = (ldtb - 1) / 3;
#line 252 "dsytrf_aa_2stage.f"
    }
#line 253 "dsytrf_aa_2stage.f"
    if (*lwork < nb * *n) {
#line 254 "dsytrf_aa_2stage.f"
	nb = *lwork / *n;
#line 255 "dsytrf_aa_2stage.f"
    }

/*     Determine the number of the block columns */

#line 259 "dsytrf_aa_2stage.f"
    nt = (*n + nb - 1) / nb;
#line 260 "dsytrf_aa_2stage.f"
    td = nb << 1;
#line 261 "dsytrf_aa_2stage.f"
    kb = min(nb,*n);

/*     Initialize vectors/matrices */

#line 265 "dsytrf_aa_2stage.f"
    i__1 = kb;
#line 265 "dsytrf_aa_2stage.f"
    for (j = 1; j <= i__1; ++j) {
#line 266 "dsytrf_aa_2stage.f"
	ipiv[j] = j;
#line 267 "dsytrf_aa_2stage.f"
    }

/*     Save NB */

#line 271 "dsytrf_aa_2stage.f"
    tb[1] = (doublereal) nb;

#line 273 "dsytrf_aa_2stage.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

#line 279 "dsytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 279 "dsytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 283 "dsytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 283 "dsytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 284 "dsytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 284 "dsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 285 "dsytrf_aa_2stage.f"
		if (i__ == 1) {
/*                 H(I,J) = T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 287 "dsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 288 "dsytrf_aa_2stage.f"
			jb = nb + kb;
#line 289 "dsytrf_aa_2stage.f"
		    } else {
#line 290 "dsytrf_aa_2stage.f"
			jb = nb << 1;
#line 291 "dsytrf_aa_2stage.f"
		    }
#line 292 "dsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 292 "dsytrf_aa_2stage.f"
		    dgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &
			    c_b12, &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(
			    i__ - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &
			    c_b13, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)11);
#line 297 "dsytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 299 "dsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 300 "dsytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 301 "dsytrf_aa_2stage.f"
		    } else {
#line 302 "dsytrf_aa_2stage.f"
			jb = nb * 3;
#line 303 "dsytrf_aa_2stage.f"
		    }
#line 304 "dsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 304 "dsytrf_aa_2stage.f"
		    dgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &
			    c_b12, &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &
			    i__3, &a[(i__ - 2) * nb + 1 + (j * nb + 1) * 
			    a_dim1], lda, &c_b13, &work[i__ * nb + 1], n, (
			    ftnlen)11, (ftnlen)11);
#line 310 "dsytrf_aa_2stage.f"
		}
#line 311 "dsytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 315 "dsytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 315 "dsytrf_aa_2stage.f"
	    dlacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 317 "dsytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = U(1:J,J)'*H(1:J) */
#line 319 "dsytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 319 "dsytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 319 "dsytrf_aa_2stage.f"
		dgemm_("Transpose", "NoTranspose", &kb, &kb, &i__2, &c_b21, &
			a[(j * nb + 1) * a_dim1 + 1], lda, &work[nb + 1], n, &
			c_b12, &tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)9, 
			(ftnlen)11);
/*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
#line 325 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 325 "dsytrf_aa_2stage.f"
		dgemm_("Transpose", "NoTranspose", &kb, &nb, &kb, &c_b12, &a[(
			j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td 
			+ nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b13, &work[
			1], n, (ftnlen)9, (ftnlen)11);
#line 330 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 330 "dsytrf_aa_2stage.f"
		dgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &c_b21, &
			work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &c_b12, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)11);
#line 335 "dsytrf_aa_2stage.f"
	    }
#line 336 "dsytrf_aa_2stage.f"
	    if (j > 0) {
#line 337 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 337 "dsytrf_aa_2stage.f"
		dsygst_(&c__1, "Upper", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 340 "dsytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 344 "dsytrf_aa_2stage.f"
	    i__2 = kb;
#line 344 "dsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 345 "dsytrf_aa_2stage.f"
		i__3 = kb;
#line 345 "dsytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 346 "dsytrf_aa_2stage.f"
		    tb[td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb] = tb[
			    td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb];
#line 348 "dsytrf_aa_2stage.f"
		}
#line 349 "dsytrf_aa_2stage.f"
	    }

#line 351 "dsytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 352 "dsytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 356 "dsytrf_aa_2stage.f"
		    if (j == 1) {
#line 357 "dsytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 357 "dsytrf_aa_2stage.f"
			dgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, &
				c_b12, &tb[td + 1 + j * nb * ldtb], &i__2, &a[
				(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
				lda, &c_b13, &work[j * nb + 1], n, (ftnlen)11,
				 (ftnlen)11);
#line 362 "dsytrf_aa_2stage.f"
		    } else {
#line 363 "dsytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 363 "dsytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 363 "dsytrf_aa_2stage.f"
			dgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, 
				&c_b12, &tb[td + nb + 1 + (j - 1) * nb * ldtb]
				, &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
				a_dim1], lda, &c_b13, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)11);
#line 369 "dsytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 373 "dsytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 373 "dsytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 373 "dsytrf_aa_2stage.f"
		    dgemm_("Transpose", "NoTranspose", &nb, &i__2, &i__3, &
			    c_b21, &work[nb + 1], n, &a[((j + 1) * nb + 1) * 
			    a_dim1 + 1], lda, &c_b12, &a[j * nb + 1 + ((j + 1)
			     * nb + 1) * a_dim1], lda, (ftnlen)9, (ftnlen)11);
#line 378 "dsytrf_aa_2stage.f"
		}

/*              Copy panel to workspace to call DGETRF */

#line 382 "dsytrf_aa_2stage.f"
		i__2 = nb;
#line 382 "dsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 383 "dsytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 383 "dsytrf_aa_2stage.f"
		    dcopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1]
			    , lda, &work[(k - 1) * *n + 1], &c__1);
#line 386 "dsytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 390 "dsytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 390 "dsytrf_aa_2stage.f"
		dgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], &
			iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Copy panel back */

#line 399 "dsytrf_aa_2stage.f"
		i__2 = nb;
#line 399 "dsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 400 "dsytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 400 "dsytrf_aa_2stage.f"
		    dcopy_(&i__3, &work[(k - 1) * *n + 1], &c__1, &a[j * nb + 
			    k + ((j + 1) * nb + 1) * a_dim1], lda);
#line 403 "dsytrf_aa_2stage.f"
		}

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 407 "dsytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 407 "dsytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 408 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 408 "dsytrf_aa_2stage.f"
		dlaset_("Full", &kb, &nb, &c_b13, &c_b13, &tb[td + nb + 1 + j 
			* nb * ldtb], &i__2, (ftnlen)4);
#line 410 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 410 "dsytrf_aa_2stage.f"
		dlacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)5);
#line 413 "dsytrf_aa_2stage.f"
		if (j > 0) {
#line 414 "dsytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 414 "dsytrf_aa_2stage.f"
		    dtrsm_("R", "U", "N", "U", &kb, &nb, &c_b12, &a[(j - 1) * 
			    nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 417 "dsytrf_aa_2stage.f"
		}

/*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
/*              updates */

#line 422 "dsytrf_aa_2stage.f"
		i__2 = nb;
#line 422 "dsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 423 "dsytrf_aa_2stage.f"
		    i__3 = kb;
#line 423 "dsytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 424 "dsytrf_aa_2stage.f"
			tb[td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1) * 
				ldtb] = tb[td + nb + i__ - k + 1 + (j * nb + 
				k - 1) * ldtb];
#line 426 "dsytrf_aa_2stage.f"
		    }
#line 427 "dsytrf_aa_2stage.f"
		}
#line 428 "dsytrf_aa_2stage.f"
		dlaset_("Lower", &kb, &nb, &c_b13, &c_b12, &a[j * nb + 1 + ((
			j + 1) * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 433 "dsytrf_aa_2stage.f"
		i__2 = kb;
#line 433 "dsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 435 "dsytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 437 "dsytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 438 "dsytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 439 "dsytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 441 "dsytrf_aa_2stage.f"
			i__3 = k - 1;
#line 441 "dsytrf_aa_2stage.f"
			dswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], &
				c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], &
				c__1);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 444 "dsytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 444 "dsytrf_aa_2stage.f"
			dswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[i1 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 447 "dsytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 447 "dsytrf_aa_2stage.f"
			dswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[i2 
				+ (i2 + 1) * a_dim1], lda);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 450 "dsytrf_aa_2stage.f"
			piv = a[i1 + i1 * a_dim1];
#line 451 "dsytrf_aa_2stage.f"
			a[i1 + i1 * a_dim1] = a[i2 + i2 * a_dim1];
#line 452 "dsytrf_aa_2stage.f"
			a[i2 + i2 * a_dim1] = piv;
/*                    > Apply pivots to previous columns of L */
#line 454 "dsytrf_aa_2stage.f"
			if (j > 0) {
#line 455 "dsytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 455 "dsytrf_aa_2stage.f"
			    dswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * 
				    a_dim1 + 1], &c__1);
#line 457 "dsytrf_aa_2stage.f"
			}
#line 458 "dsytrf_aa_2stage.f"
		    }
#line 459 "dsytrf_aa_2stage.f"
		}
#line 460 "dsytrf_aa_2stage.f"
	    }
#line 461 "dsytrf_aa_2stage.f"
	}
#line 462 "dsytrf_aa_2stage.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 468 "dsytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 468 "dsytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 472 "dsytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 472 "dsytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 473 "dsytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 473 "dsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 474 "dsytrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
#line 476 "dsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 477 "dsytrf_aa_2stage.f"
			jb = nb + kb;
#line 478 "dsytrf_aa_2stage.f"
		    } else {
#line 479 "dsytrf_aa_2stage.f"
			jb = nb << 1;
#line 480 "dsytrf_aa_2stage.f"
		    }
#line 481 "dsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 481 "dsytrf_aa_2stage.f"
		    dgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b12, 
			    &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[j * nb + 
			    1 + ((i__ - 1) * nb + 1) * a_dim1], lda, &c_b13, &
			    work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)9);
#line 486 "dsytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
#line 488 "dsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 489 "dsytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 490 "dsytrf_aa_2stage.f"
		    } else {
#line 491 "dsytrf_aa_2stage.f"
			jb = nb * 3;
#line 492 "dsytrf_aa_2stage.f"
		    }
#line 493 "dsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 493 "dsytrf_aa_2stage.f"
		    dgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b12, 
			    &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, &
			    a[j * nb + 1 + ((i__ - 2) * nb + 1) * a_dim1], 
			    lda, &c_b13, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)9);
#line 499 "dsytrf_aa_2stage.f"
		}
#line 500 "dsytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 504 "dsytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 504 "dsytrf_aa_2stage.f"
	    dlacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 506 "dsytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = L(J,1:J)*H(1:J) */
#line 508 "dsytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 508 "dsytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 508 "dsytrf_aa_2stage.f"
		dgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &c_b21, 
			&a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &
			c_b12, &tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)11,
			 (ftnlen)11);
/*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
#line 514 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 514 "dsytrf_aa_2stage.f"
		dgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b12, &
			a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[
			td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b13, &
			work[1], n, (ftnlen)11, (ftnlen)11);
#line 519 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 519 "dsytrf_aa_2stage.f"
		dgemm_("NoTranspose", "Transpose", &kb, &kb, &nb, &c_b21, &
			work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
			a_dim1], lda, &c_b12, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)9);
#line 524 "dsytrf_aa_2stage.f"
	    }
#line 525 "dsytrf_aa_2stage.f"
	    if (j > 0) {
#line 526 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 526 "dsytrf_aa_2stage.f"
		dsygst_(&c__1, "Lower", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 529 "dsytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 533 "dsytrf_aa_2stage.f"
	    i__2 = kb;
#line 533 "dsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 534 "dsytrf_aa_2stage.f"
		i__3 = kb;
#line 534 "dsytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 535 "dsytrf_aa_2stage.f"
		    tb[td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb] = tb[
			    td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb];
#line 537 "dsytrf_aa_2stage.f"
		}
#line 538 "dsytrf_aa_2stage.f"
	    }

#line 540 "dsytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 541 "dsytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 545 "dsytrf_aa_2stage.f"
		    if (j == 1) {
#line 546 "dsytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 546 "dsytrf_aa_2stage.f"
			dgemm_("NoTranspose", "Transpose", &kb, &kb, &kb, &
				c_b12, &tb[td + 1 + j * nb * ldtb], &i__2, &a[
				j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
				lda, &c_b13, &work[j * nb + 1], n, (ftnlen)11,
				 (ftnlen)9);
#line 551 "dsytrf_aa_2stage.f"
		    } else {
#line 552 "dsytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 552 "dsytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 552 "dsytrf_aa_2stage.f"
			dgemm_("NoTranspose", "Transpose", &kb, &kb, &i__2, &
				c_b12, &tb[td + nb + 1 + (j - 1) * nb * ldtb],
				 &i__3, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
				a_dim1], lda, &c_b13, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)9);
#line 558 "dsytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 562 "dsytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 562 "dsytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 562 "dsytrf_aa_2stage.f"
		    dgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, &
			    c_b21, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[
			    nb + 1], n, &c_b12, &a[(j + 1) * nb + 1 + (j * nb 
			    + 1) * a_dim1], lda, (ftnlen)11, (ftnlen)11);
#line 567 "dsytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 571 "dsytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 571 "dsytrf_aa_2stage.f"
		dgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 580 "dsytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 580 "dsytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 581 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 581 "dsytrf_aa_2stage.f"
		dlaset_("Full", &kb, &nb, &c_b13, &c_b13, &tb[td + nb + 1 + j 
			* nb * ldtb], &i__2, (ftnlen)4);
#line 583 "dsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 583 "dsytrf_aa_2stage.f"
		dlacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) 
			* a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &
			i__2, (ftnlen)5);
#line 586 "dsytrf_aa_2stage.f"
		if (j > 0) {
#line 587 "dsytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 587 "dsytrf_aa_2stage.f"
		    dtrsm_("R", "L", "T", "U", &kb, &nb, &c_b12, &a[j * nb + 
			    1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + 
			    nb + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 590 "dsytrf_aa_2stage.f"
		}

/*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
/*              updates */

#line 595 "dsytrf_aa_2stage.f"
		i__2 = nb;
#line 595 "dsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 596 "dsytrf_aa_2stage.f"
		    i__3 = kb;
#line 596 "dsytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 597 "dsytrf_aa_2stage.f"
			tb[td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1) * 
				ldtb] = tb[td + nb + i__ - k + 1 + (j * nb + 
				k - 1) * ldtb];
#line 599 "dsytrf_aa_2stage.f"
		    }
#line 600 "dsytrf_aa_2stage.f"
		}
#line 601 "dsytrf_aa_2stage.f"
		dlaset_("Upper", &kb, &nb, &c_b13, &c_b12, &a[(j + 1) * nb + 
			1 + (j * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 606 "dsytrf_aa_2stage.f"
		i__2 = kb;
#line 606 "dsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 608 "dsytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 610 "dsytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 611 "dsytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 612 "dsytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 614 "dsytrf_aa_2stage.f"
			i__3 = k - 1;
#line 614 "dsytrf_aa_2stage.f"
			dswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], 
				lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], 
				lda);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 617 "dsytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 617 "dsytrf_aa_2stage.f"
			dswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ (i1 + 1) * a_dim1], lda);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 620 "dsytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 620 "dsytrf_aa_2stage.f"
			dswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 623 "dsytrf_aa_2stage.f"
			piv = a[i1 + i1 * a_dim1];
#line 624 "dsytrf_aa_2stage.f"
			a[i1 + i1 * a_dim1] = a[i2 + i2 * a_dim1];
#line 625 "dsytrf_aa_2stage.f"
			a[i2 + i2 * a_dim1] = piv;
/*                    > Apply pivots to previous columns of L */
#line 627 "dsytrf_aa_2stage.f"
			if (j > 0) {
#line 628 "dsytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 628 "dsytrf_aa_2stage.f"
			    dswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + 
				    a_dim1], lda);
#line 630 "dsytrf_aa_2stage.f"
			}
#line 631 "dsytrf_aa_2stage.f"
		    }
#line 632 "dsytrf_aa_2stage.f"
		}

/*              Apply pivots to previous columns of L */

/*               CALL DLASWP( J*NB, A( 1, 1 ), LDA, */
/*     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
#line 638 "dsytrf_aa_2stage.f"
	    }
#line 639 "dsytrf_aa_2stage.f"
	}
#line 640 "dsytrf_aa_2stage.f"
    }

/*     Factor the band matrix */
#line 643 "dsytrf_aa_2stage.f"
    dgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);

/*     End of DSYTRF_AA_2STAGE */

#line 647 "dsytrf_aa_2stage.f"
    return 0;
} /* dsytrf_aa_2stage__ */

