#line 1 "zsytrf_aa_2stage.f"
/* zsytrf_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zsytrf_aa_2stage.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZSYTRF_AA_2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRF_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrf_
aa_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrf_
aa_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrf_
aa_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE ZSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
/*                                   IPIV2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LTB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IPIV2( * ) */
/*       COMPLEX*16         A( LDA, * ), TB( * ), WORK( * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRF_AA_2STAGE computes the factorization of a complex symmetric matrix A */
/* > using the Aasen's algorithm.  The form of the factorization is */
/* > */
/* >    A = U*T*U**T  or  A = L*T*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a complex symmetric band matrix with the */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the hermitian matrix A.  If UPLO = 'U', the leading */
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
/* >          TB is COMPLEX*16 array, dimension (LTB) */
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
/* >          WORK is COMPLEX*16 workspace of size LWORK */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsytrf_aa_2stage__(char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublecomplex *tb, integer *ltb, integer *ipiv, 
	integer *ipiv2, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, k, i1, i2, jb, kb, nb, td, nt;
    static doublecomplex piv;
    static integer ldtb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ztrsm_(char *, char *, 
	    char *, char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgbtrf_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, integer *), 
	    zgetrf_(integer *, integer *, doublecomplex *, integer *, integer 
	    *, integer *), zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
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

#line 205 "zsytrf_aa_2stage.f"
    /* Parameter adjustments */
#line 205 "zsytrf_aa_2stage.f"
    a_dim1 = *lda;
#line 205 "zsytrf_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 205 "zsytrf_aa_2stage.f"
    a -= a_offset;
#line 205 "zsytrf_aa_2stage.f"
    --tb;
#line 205 "zsytrf_aa_2stage.f"
    --ipiv;
#line 205 "zsytrf_aa_2stage.f"
    --ipiv2;
#line 205 "zsytrf_aa_2stage.f"
    --work;
#line 205 "zsytrf_aa_2stage.f"

#line 205 "zsytrf_aa_2stage.f"
    /* Function Body */
#line 205 "zsytrf_aa_2stage.f"
    *info = 0;
#line 206 "zsytrf_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "zsytrf_aa_2stage.f"
    wquery = *lwork == -1;
#line 208 "zsytrf_aa_2stage.f"
    tquery = *ltb == -1;
#line 209 "zsytrf_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 210 "zsytrf_aa_2stage.f"
	*info = -1;
#line 211 "zsytrf_aa_2stage.f"
    } else if (*n < 0) {
#line 212 "zsytrf_aa_2stage.f"
	*info = -2;
#line 213 "zsytrf_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 214 "zsytrf_aa_2stage.f"
	*info = -4;
#line 215 "zsytrf_aa_2stage.f"
    } else if (*ltb < *n << 2 && ! tquery) {
#line 216 "zsytrf_aa_2stage.f"
	*info = -6;
#line 217 "zsytrf_aa_2stage.f"
    } else if (*lwork < *n && ! wquery) {
#line 218 "zsytrf_aa_2stage.f"
	*info = -10;
#line 219 "zsytrf_aa_2stage.f"
    }

#line 221 "zsytrf_aa_2stage.f"
    if (*info != 0) {
#line 222 "zsytrf_aa_2stage.f"
	i__1 = -(*info);
#line 222 "zsytrf_aa_2stage.f"
	xerbla_("ZSYTRF_AA_2STAGE", &i__1, (ftnlen)16);
#line 223 "zsytrf_aa_2stage.f"
	return 0;
#line 224 "zsytrf_aa_2stage.f"
    }

/*     Answer the query */

#line 228 "zsytrf_aa_2stage.f"
    nb = ilaenv_(&c__1, "ZSYTRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)16, (ftnlen)1);
#line 229 "zsytrf_aa_2stage.f"
    if (*info == 0) {
#line 230 "zsytrf_aa_2stage.f"
	if (tquery) {
#line 231 "zsytrf_aa_2stage.f"
	    i__1 = (nb * 3 + 1) * *n;
#line 231 "zsytrf_aa_2stage.f"
	    tb[1].r = (doublereal) i__1, tb[1].i = 0.;
#line 232 "zsytrf_aa_2stage.f"
	}
#line 233 "zsytrf_aa_2stage.f"
	if (wquery) {
#line 234 "zsytrf_aa_2stage.f"
	    i__1 = *n * nb;
#line 234 "zsytrf_aa_2stage.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 235 "zsytrf_aa_2stage.f"
	}
#line 236 "zsytrf_aa_2stage.f"
    }
#line 237 "zsytrf_aa_2stage.f"
    if (tquery || wquery) {
#line 238 "zsytrf_aa_2stage.f"
	return 0;
#line 239 "zsytrf_aa_2stage.f"
    }

/*     Quick return */

#line 243 "zsytrf_aa_2stage.f"
    if (*n == 0) {
#line 244 "zsytrf_aa_2stage.f"
	return 0;
#line 245 "zsytrf_aa_2stage.f"
    }

/*     Determine the number of the block size */

#line 249 "zsytrf_aa_2stage.f"
    ldtb = *ltb / *n;
#line 250 "zsytrf_aa_2stage.f"
    if (ldtb < nb * 3 + 1) {
#line 251 "zsytrf_aa_2stage.f"
	nb = (ldtb - 1) / 3;
#line 252 "zsytrf_aa_2stage.f"
    }
#line 253 "zsytrf_aa_2stage.f"
    if (*lwork < nb * *n) {
#line 254 "zsytrf_aa_2stage.f"
	nb = *lwork / *n;
#line 255 "zsytrf_aa_2stage.f"
    }

/*     Determine the number of the block columns */

#line 259 "zsytrf_aa_2stage.f"
    nt = (*n + nb - 1) / nb;
#line 260 "zsytrf_aa_2stage.f"
    td = nb << 1;
#line 261 "zsytrf_aa_2stage.f"
    kb = min(nb,*n);

/*     Initialize vectors/matrices */

#line 265 "zsytrf_aa_2stage.f"
    i__1 = kb;
#line 265 "zsytrf_aa_2stage.f"
    for (j = 1; j <= i__1; ++j) {
#line 266 "zsytrf_aa_2stage.f"
	ipiv[j] = j;
#line 267 "zsytrf_aa_2stage.f"
    }

/*     Save NB */

#line 271 "zsytrf_aa_2stage.f"
    tb[1].r = (doublereal) nb, tb[1].i = 0.;

#line 273 "zsytrf_aa_2stage.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

#line 279 "zsytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 279 "zsytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 283 "zsytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 283 "zsytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 284 "zsytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 284 "zsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 285 "zsytrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J) */
#line 287 "zsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 288 "zsytrf_aa_2stage.f"
			jb = nb + kb;
#line 289 "zsytrf_aa_2stage.f"
		    } else {
#line 290 "zsytrf_aa_2stage.f"
			jb = nb << 1;
#line 291 "zsytrf_aa_2stage.f"
		    }
#line 292 "zsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 292 "zsytrf_aa_2stage.f"
		    zgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(i__ - 
			    1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, 
			    &work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)11);
#line 297 "zsytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 299 "zsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 300 "zsytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 301 "zsytrf_aa_2stage.f"
		    } else {
#line 302 "zsytrf_aa_2stage.f"
			jb = nb * 3;
#line 303 "zsytrf_aa_2stage.f"
		    }
#line 304 "zsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 304 "zsytrf_aa_2stage.f"
		    zgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, 
			    &a[(i__ - 2) * nb + 1 + (j * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)11);
#line 310 "zsytrf_aa_2stage.f"
		}
#line 311 "zsytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 315 "zsytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 315 "zsytrf_aa_2stage.f"
	    zlacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 317 "zsytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = U(1:J,J)'*H(1:J) */
#line 319 "zsytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 319 "zsytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 319 "zsytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 319 "zsytrf_aa_2stage.f"
		zgemm_("Transpose", "NoTranspose", &kb, &kb, &i__2, &z__1, &a[
			(j * nb + 1) * a_dim1 + 1], lda, &work[nb + 1], n, &
			c_b2, &tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)9, (
			ftnlen)11);
/*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
#line 325 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 325 "zsytrf_aa_2stage.f"
		zgemm_("Transpose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[(
			j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td 
			+ nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &work[
			1], n, (ftnlen)9, (ftnlen)11);
#line 330 "zsytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 330 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 330 "zsytrf_aa_2stage.f"
		zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &z__1, &
			work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)11);
#line 335 "zsytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 339 "zsytrf_aa_2stage.f"
	    i__2 = kb;
#line 339 "zsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 340 "zsytrf_aa_2stage.f"
		i__3 = kb;
#line 340 "zsytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 341 "zsytrf_aa_2stage.f"
		    i__4 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 341 "zsytrf_aa_2stage.f"
		    i__5 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 341 "zsytrf_aa_2stage.f"
		    tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 343 "zsytrf_aa_2stage.f"
		}
#line 344 "zsytrf_aa_2stage.f"
	    }
#line 345 "zsytrf_aa_2stage.f"
	    if (j > 0) {
/*               CALL CHEGST( 1, 'Upper', KB, */
/*     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1, */
/*     $                      A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO ) */
#line 349 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 349 "zsytrf_aa_2stage.f"
		ztrsm_("L", "U", "T", "N", &kb, &kb, &c_b2, &a[(j - 1) * nb + 
			1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 352 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 352 "zsytrf_aa_2stage.f"
		ztrsm_("R", "U", "N", "N", &kb, &kb, &c_b2, &a[(j - 1) * nb + 
			1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 355 "zsytrf_aa_2stage.f"
	    }

#line 357 "zsytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 358 "zsytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 362 "zsytrf_aa_2stage.f"
		    if (j == 1) {
#line 363 "zsytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 363 "zsytrf_aa_2stage.f"
			zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, &
				c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[(
				j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda,
				 &c_b1, &work[j * nb + 1], n, (ftnlen)11, (
				ftnlen)11);
#line 368 "zsytrf_aa_2stage.f"
		    } else {
#line 369 "zsytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 369 "zsytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 369 "zsytrf_aa_2stage.f"
			zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, 
				&c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb],
				 &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)11);
#line 375 "zsytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 379 "zsytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 379 "zsytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 379 "zsytrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 379 "zsytrf_aa_2stage.f"
		    zgemm_("Transpose", "NoTranspose", &nb, &i__2, &i__3, &
			    z__1, &work[nb + 1], n, &a[((j + 1) * nb + 1) * 
			    a_dim1 + 1], lda, &c_b2, &a[j * nb + 1 + ((j + 1) 
			    * nb + 1) * a_dim1], lda, (ftnlen)9, (ftnlen)11);
#line 384 "zsytrf_aa_2stage.f"
		}

/*              Copy panel to workspace to call ZGETRF */

#line 388 "zsytrf_aa_2stage.f"
		i__2 = nb;
#line 388 "zsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 389 "zsytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 389 "zsytrf_aa_2stage.f"
		    zcopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1]
			    , lda, &work[(k - 1) * *n + 1], &c__1);
#line 392 "zsytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 396 "zsytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 396 "zsytrf_aa_2stage.f"
		zgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], &
			iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Copy panel back */

#line 405 "zsytrf_aa_2stage.f"
		i__2 = nb;
#line 405 "zsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 406 "zsytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 406 "zsytrf_aa_2stage.f"
		    zcopy_(&i__3, &work[(k - 1) * *n + 1], &c__1, &a[j * nb + 
			    k + ((j + 1) * nb + 1) * a_dim1], lda);
#line 409 "zsytrf_aa_2stage.f"
		}

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 413 "zsytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 413 "zsytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 414 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 414 "zsytrf_aa_2stage.f"
		zlaset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 416 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 416 "zsytrf_aa_2stage.f"
		zlacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)5);
#line 419 "zsytrf_aa_2stage.f"
		if (j > 0) {
#line 420 "zsytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 420 "zsytrf_aa_2stage.f"
		    ztrsm_("R", "U", "N", "U", &kb, &nb, &c_b2, &a[(j - 1) * 
			    nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 423 "zsytrf_aa_2stage.f"
		}

/*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
/*              updates */

#line 428 "zsytrf_aa_2stage.f"
		i__2 = nb;
#line 428 "zsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 429 "zsytrf_aa_2stage.f"
		    i__3 = kb;
#line 429 "zsytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 430 "zsytrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 430 "zsytrf_aa_2stage.f"
			i__5 = td + nb + i__ - k + 1 + (j * nb + k - 1) * 
				ldtb;
#line 430 "zsytrf_aa_2stage.f"
			tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 432 "zsytrf_aa_2stage.f"
		    }
#line 433 "zsytrf_aa_2stage.f"
		}
#line 434 "zsytrf_aa_2stage.f"
		zlaset_("Lower", &kb, &nb, &c_b1, &c_b2, &a[j * nb + 1 + ((j 
			+ 1) * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 439 "zsytrf_aa_2stage.f"
		i__2 = kb;
#line 439 "zsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 441 "zsytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 443 "zsytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 444 "zsytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 445 "zsytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 447 "zsytrf_aa_2stage.f"
			i__3 = k - 1;
#line 447 "zsytrf_aa_2stage.f"
			zswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], &
				c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], &
				c__1);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 450 "zsytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 450 "zsytrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[i1 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 453 "zsytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 453 "zsytrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[i2 
				+ (i2 + 1) * a_dim1], lda);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 456 "zsytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 456 "zsytrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 457 "zsytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 457 "zsytrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 457 "zsytrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 458 "zsytrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 458 "zsytrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 460 "zsytrf_aa_2stage.f"
			if (j > 0) {
#line 461 "zsytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 461 "zsytrf_aa_2stage.f"
			    zswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * 
				    a_dim1 + 1], &c__1);
#line 463 "zsytrf_aa_2stage.f"
			}
#line 464 "zsytrf_aa_2stage.f"
		    }
#line 465 "zsytrf_aa_2stage.f"
		}
#line 466 "zsytrf_aa_2stage.f"
	    }
#line 467 "zsytrf_aa_2stage.f"
	}
#line 468 "zsytrf_aa_2stage.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 474 "zsytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 474 "zsytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 478 "zsytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 478 "zsytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 479 "zsytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 479 "zsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 480 "zsytrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
#line 482 "zsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 483 "zsytrf_aa_2stage.f"
			jb = nb + kb;
#line 484 "zsytrf_aa_2stage.f"
		    } else {
#line 485 "zsytrf_aa_2stage.f"
			jb = nb << 1;
#line 486 "zsytrf_aa_2stage.f"
		    }
#line 487 "zsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 487 "zsytrf_aa_2stage.f"
		    zgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b2, &
			    tb[td + 1 + i__ * nb * ldtb], &i__3, &a[j * nb + 
			    1 + ((i__ - 1) * nb + 1) * a_dim1], lda, &c_b1, &
			    work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)9);
#line 492 "zsytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
#line 494 "zsytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 495 "zsytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 496 "zsytrf_aa_2stage.f"
		    } else {
#line 497 "zsytrf_aa_2stage.f"
			jb = nb * 3;
#line 498 "zsytrf_aa_2stage.f"
		    }
#line 499 "zsytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 499 "zsytrf_aa_2stage.f"
		    zgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b2, &
			    tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, &
			    a[j * nb + 1 + ((i__ - 2) * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)9);
#line 505 "zsytrf_aa_2stage.f"
		}
#line 506 "zsytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 510 "zsytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 510 "zsytrf_aa_2stage.f"
	    zlacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 512 "zsytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = L(J,1:J)*H(1:J) */
#line 514 "zsytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 514 "zsytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 514 "zsytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 514 "zsytrf_aa_2stage.f"
		zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &z__1, &
			a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &c_b2, 
			&tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)11, (
			ftnlen)11);
/*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
#line 520 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 520 "zsytrf_aa_2stage.f"
		zgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[
			j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[
			td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &
			work[1], n, (ftnlen)11, (ftnlen)11);
#line 525 "zsytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 525 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 525 "zsytrf_aa_2stage.f"
		zgemm_("NoTranspose", "Transpose", &kb, &kb, &nb, &z__1, &
			work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
			a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)9);
#line 530 "zsytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 534 "zsytrf_aa_2stage.f"
	    i__2 = kb;
#line 534 "zsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 535 "zsytrf_aa_2stage.f"
		i__3 = kb;
#line 535 "zsytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 536 "zsytrf_aa_2stage.f"
		    i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 536 "zsytrf_aa_2stage.f"
		    i__5 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 536 "zsytrf_aa_2stage.f"
		    tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 538 "zsytrf_aa_2stage.f"
		}
#line 539 "zsytrf_aa_2stage.f"
	    }
#line 540 "zsytrf_aa_2stage.f"
	    if (j > 0) {
/*               CALL CHEGST( 1, 'Lower', KB, */
/*     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1, */
/*     $                      A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO ) */
#line 544 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 544 "zsytrf_aa_2stage.f"
		ztrsm_("L", "L", "N", "N", &kb, &kb, &c_b2, &a[j * nb + 1 + ((
			j - 1) * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 547 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 547 "zsytrf_aa_2stage.f"
		ztrsm_("R", "L", "T", "N", &kb, &kb, &c_b2, &a[j * nb + 1 + ((
			j - 1) * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 550 "zsytrf_aa_2stage.f"
	    }

/*           Symmetrize T(J,J) */

#line 554 "zsytrf_aa_2stage.f"
	    i__2 = kb;
#line 554 "zsytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 555 "zsytrf_aa_2stage.f"
		i__3 = kb;
#line 555 "zsytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 556 "zsytrf_aa_2stage.f"
		    i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 556 "zsytrf_aa_2stage.f"
		    i__5 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 556 "zsytrf_aa_2stage.f"
		    tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 558 "zsytrf_aa_2stage.f"
		}
#line 559 "zsytrf_aa_2stage.f"
	    }

#line 561 "zsytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 562 "zsytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 566 "zsytrf_aa_2stage.f"
		    if (j == 1) {
#line 567 "zsytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 567 "zsytrf_aa_2stage.f"
			zgemm_("NoTranspose", "Transpose", &kb, &kb, &kb, &
				c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[
				j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
				lda, &c_b1, &work[j * nb + 1], n, (ftnlen)11, 
				(ftnlen)9);
#line 572 "zsytrf_aa_2stage.f"
		    } else {
#line 573 "zsytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 573 "zsytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 573 "zsytrf_aa_2stage.f"
			zgemm_("NoTranspose", "Transpose", &kb, &kb, &i__2, &
				c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb], 
				&i__3, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)9);
#line 579 "zsytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 583 "zsytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 583 "zsytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 583 "zsytrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 583 "zsytrf_aa_2stage.f"
		    zgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, &
			    z__1, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[
			    nb + 1], n, &c_b2, &a[(j + 1) * nb + 1 + (j * nb 
			    + 1) * a_dim1], lda, (ftnlen)11, (ftnlen)11);
#line 588 "zsytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 592 "zsytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 592 "zsytrf_aa_2stage.f"
		zgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 601 "zsytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 601 "zsytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 602 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 602 "zsytrf_aa_2stage.f"
		zlaset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 604 "zsytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 604 "zsytrf_aa_2stage.f"
		zlacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) 
			* a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &
			i__2, (ftnlen)5);
#line 607 "zsytrf_aa_2stage.f"
		if (j > 0) {
#line 608 "zsytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 608 "zsytrf_aa_2stage.f"
		    ztrsm_("R", "L", "T", "U", &kb, &nb, &c_b2, &a[j * nb + 1 
			    + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 611 "zsytrf_aa_2stage.f"
		}

/*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
/*              updates */

#line 616 "zsytrf_aa_2stage.f"
		i__2 = nb;
#line 616 "zsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 617 "zsytrf_aa_2stage.f"
		    i__3 = kb;
#line 617 "zsytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 618 "zsytrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 618 "zsytrf_aa_2stage.f"
			i__5 = td + nb + i__ - k + 1 + (j * nb + k - 1) * 
				ldtb;
#line 618 "zsytrf_aa_2stage.f"
			tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 620 "zsytrf_aa_2stage.f"
		    }
#line 621 "zsytrf_aa_2stage.f"
		}
#line 622 "zsytrf_aa_2stage.f"
		zlaset_("Upper", &kb, &nb, &c_b1, &c_b2, &a[(j + 1) * nb + 1 
			+ (j * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 627 "zsytrf_aa_2stage.f"
		i__2 = kb;
#line 627 "zsytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 629 "zsytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 631 "zsytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 632 "zsytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 633 "zsytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 635 "zsytrf_aa_2stage.f"
			i__3 = k - 1;
#line 635 "zsytrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], 
				lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], 
				lda);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 638 "zsytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 638 "zsytrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ (i1 + 1) * a_dim1], lda);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 641 "zsytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 641 "zsytrf_aa_2stage.f"
			zswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 644 "zsytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 644 "zsytrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 645 "zsytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 645 "zsytrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 645 "zsytrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 646 "zsytrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 646 "zsytrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 648 "zsytrf_aa_2stage.f"
			if (j > 0) {
#line 649 "zsytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 649 "zsytrf_aa_2stage.f"
			    zswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + 
				    a_dim1], lda);
#line 651 "zsytrf_aa_2stage.f"
			}
#line 652 "zsytrf_aa_2stage.f"
		    }
#line 653 "zsytrf_aa_2stage.f"
		}

/*              Apply pivots to previous columns of L */

/*               CALL ZLASWP( J*NB, A( 1, 1 ), LDA, */
/*     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
#line 659 "zsytrf_aa_2stage.f"
	    }
#line 660 "zsytrf_aa_2stage.f"
	}
#line 661 "zsytrf_aa_2stage.f"
    }

/*     Factor the band matrix */
#line 664 "zsytrf_aa_2stage.f"
    zgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);

/*     End of ZSYTRF_AA_2STAGE */

#line 668 "zsytrf_aa_2stage.f"
    return 0;
} /* zsytrf_aa_2stage__ */

