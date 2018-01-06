#line 1 "csytrf_aa_2stage.f"
/* csytrf_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "csytrf_aa_2stage.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CSYTRF_AA_2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRF_AA_2STAGE + dependencies */
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

/*      SUBROUTINE CSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
/*                                   IPIV2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LTB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IPIV2( * ) */
/*       COMPLEX            A( LDA, * ), TB( * ), WORK( * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRF_AA_2STAGE computes the factorization of a complex symmetric matrix A */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          TB is COMPLEX array, dimension (LTB) */
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
/* >          WORK is COMPLEX workspace of size LWORK */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csytrf_aa_2stage__(char *uplo, integer *n, doublecomplex 
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
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ctrsm_(char *, char *, 
	    char *, char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int cgbtrf_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, integer *), 
	    cgetrf_(integer *, integer *, doublecomplex *, integer *, integer 
	    *, integer *), clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

#line 205 "csytrf_aa_2stage.f"
    /* Parameter adjustments */
#line 205 "csytrf_aa_2stage.f"
    a_dim1 = *lda;
#line 205 "csytrf_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 205 "csytrf_aa_2stage.f"
    a -= a_offset;
#line 205 "csytrf_aa_2stage.f"
    --tb;
#line 205 "csytrf_aa_2stage.f"
    --ipiv;
#line 205 "csytrf_aa_2stage.f"
    --ipiv2;
#line 205 "csytrf_aa_2stage.f"
    --work;
#line 205 "csytrf_aa_2stage.f"

#line 205 "csytrf_aa_2stage.f"
    /* Function Body */
#line 205 "csytrf_aa_2stage.f"
    *info = 0;
#line 206 "csytrf_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "csytrf_aa_2stage.f"
    wquery = *lwork == -1;
#line 208 "csytrf_aa_2stage.f"
    tquery = *ltb == -1;
#line 209 "csytrf_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 210 "csytrf_aa_2stage.f"
	*info = -1;
#line 211 "csytrf_aa_2stage.f"
    } else if (*n < 0) {
#line 212 "csytrf_aa_2stage.f"
	*info = -2;
#line 213 "csytrf_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 214 "csytrf_aa_2stage.f"
	*info = -4;
#line 215 "csytrf_aa_2stage.f"
    } else if (*ltb < *n << 2 && ! tquery) {
#line 216 "csytrf_aa_2stage.f"
	*info = -6;
#line 217 "csytrf_aa_2stage.f"
    } else if (*lwork < *n && ! wquery) {
#line 218 "csytrf_aa_2stage.f"
	*info = -10;
#line 219 "csytrf_aa_2stage.f"
    }

#line 221 "csytrf_aa_2stage.f"
    if (*info != 0) {
#line 222 "csytrf_aa_2stage.f"
	i__1 = -(*info);
#line 222 "csytrf_aa_2stage.f"
	xerbla_("CSYTRF_AA_2STAGE", &i__1, (ftnlen)16);
#line 223 "csytrf_aa_2stage.f"
	return 0;
#line 224 "csytrf_aa_2stage.f"
    }

/*     Answer the query */

#line 228 "csytrf_aa_2stage.f"
    nb = ilaenv_(&c__1, "CSYTRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)16, (ftnlen)1);
#line 229 "csytrf_aa_2stage.f"
    if (*info == 0) {
#line 230 "csytrf_aa_2stage.f"
	if (tquery) {
#line 231 "csytrf_aa_2stage.f"
	    i__1 = (nb * 3 + 1) * *n;
#line 231 "csytrf_aa_2stage.f"
	    tb[1].r = (doublereal) i__1, tb[1].i = 0.;
#line 232 "csytrf_aa_2stage.f"
	}
#line 233 "csytrf_aa_2stage.f"
	if (wquery) {
#line 234 "csytrf_aa_2stage.f"
	    i__1 = *n * nb;
#line 234 "csytrf_aa_2stage.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 235 "csytrf_aa_2stage.f"
	}
#line 236 "csytrf_aa_2stage.f"
    }
#line 237 "csytrf_aa_2stage.f"
    if (tquery || wquery) {
#line 238 "csytrf_aa_2stage.f"
	return 0;
#line 239 "csytrf_aa_2stage.f"
    }

/*     Quick return */

#line 243 "csytrf_aa_2stage.f"
    if (*n == 0) {
#line 244 "csytrf_aa_2stage.f"
	return 0;
#line 245 "csytrf_aa_2stage.f"
    }

/*     Determine the number of the block size */

#line 249 "csytrf_aa_2stage.f"
    ldtb = *ltb / *n;
#line 250 "csytrf_aa_2stage.f"
    if (ldtb < nb * 3 + 1) {
#line 251 "csytrf_aa_2stage.f"
	nb = (ldtb - 1) / 3;
#line 252 "csytrf_aa_2stage.f"
    }
#line 253 "csytrf_aa_2stage.f"
    if (*lwork < nb * *n) {
#line 254 "csytrf_aa_2stage.f"
	nb = *lwork / *n;
#line 255 "csytrf_aa_2stage.f"
    }

/*     Determine the number of the block columns */

#line 259 "csytrf_aa_2stage.f"
    nt = (*n + nb - 1) / nb;
#line 260 "csytrf_aa_2stage.f"
    td = nb << 1;
#line 261 "csytrf_aa_2stage.f"
    kb = min(nb,*n);

/*     Initialize vectors/matrices */

#line 265 "csytrf_aa_2stage.f"
    i__1 = kb;
#line 265 "csytrf_aa_2stage.f"
    for (j = 1; j <= i__1; ++j) {
#line 266 "csytrf_aa_2stage.f"
	ipiv[j] = j;
#line 267 "csytrf_aa_2stage.f"
    }

/*     Save NB */

#line 271 "csytrf_aa_2stage.f"
    tb[1].r = (doublereal) nb, tb[1].i = 0.;

#line 273 "csytrf_aa_2stage.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

#line 279 "csytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 279 "csytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 283 "csytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 283 "csytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 284 "csytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 284 "csytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 285 "csytrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J) */
#line 287 "csytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 288 "csytrf_aa_2stage.f"
			jb = nb + kb;
#line 289 "csytrf_aa_2stage.f"
		    } else {
#line 290 "csytrf_aa_2stage.f"
			jb = nb << 1;
#line 291 "csytrf_aa_2stage.f"
		    }
#line 292 "csytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 292 "csytrf_aa_2stage.f"
		    cgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(i__ - 
			    1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, 
			    &work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)11);
#line 297 "csytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 299 "csytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 300 "csytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 301 "csytrf_aa_2stage.f"
		    } else {
#line 302 "csytrf_aa_2stage.f"
			jb = nb * 3;
#line 303 "csytrf_aa_2stage.f"
		    }
#line 304 "csytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 304 "csytrf_aa_2stage.f"
		    cgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, 
			    &a[(i__ - 2) * nb + 1 + (j * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)11);
#line 310 "csytrf_aa_2stage.f"
		}
#line 311 "csytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 315 "csytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 315 "csytrf_aa_2stage.f"
	    clacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 317 "csytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = U(1:J,J)'*H(1:J) */
#line 319 "csytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 319 "csytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 319 "csytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 319 "csytrf_aa_2stage.f"
		cgemm_("Transpose", "NoTranspose", &kb, &kb, &i__2, &z__1, &a[
			(j * nb + 1) * a_dim1 + 1], lda, &work[nb + 1], n, &
			c_b2, &tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)9, (
			ftnlen)11);
/*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
#line 325 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 325 "csytrf_aa_2stage.f"
		cgemm_("Transpose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[(
			j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td 
			+ nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &work[
			1], n, (ftnlen)9, (ftnlen)11);
#line 330 "csytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 330 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 330 "csytrf_aa_2stage.f"
		cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &z__1, &
			work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)11);
#line 335 "csytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 339 "csytrf_aa_2stage.f"
	    i__2 = kb;
#line 339 "csytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 340 "csytrf_aa_2stage.f"
		i__3 = kb;
#line 340 "csytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 341 "csytrf_aa_2stage.f"
		    i__4 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 341 "csytrf_aa_2stage.f"
		    i__5 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 341 "csytrf_aa_2stage.f"
		    tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 343 "csytrf_aa_2stage.f"
		}
#line 344 "csytrf_aa_2stage.f"
	    }
#line 345 "csytrf_aa_2stage.f"
	    if (j > 0) {
/*               CALL CHEGST( 1, 'Upper', KB, */
/*     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1, */
/*     $                      A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO ) */
#line 349 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 349 "csytrf_aa_2stage.f"
		ctrsm_("L", "U", "T", "N", &kb, &kb, &c_b2, &a[(j - 1) * nb + 
			1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 352 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 352 "csytrf_aa_2stage.f"
		ctrsm_("R", "U", "N", "N", &kb, &kb, &c_b2, &a[(j - 1) * nb + 
			1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 355 "csytrf_aa_2stage.f"
	    }

#line 357 "csytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 358 "csytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 362 "csytrf_aa_2stage.f"
		    if (j == 1) {
#line 363 "csytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 363 "csytrf_aa_2stage.f"
			cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, &
				c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[(
				j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda,
				 &c_b1, &work[j * nb + 1], n, (ftnlen)11, (
				ftnlen)11);
#line 368 "csytrf_aa_2stage.f"
		    } else {
#line 369 "csytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 369 "csytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 369 "csytrf_aa_2stage.f"
			cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, 
				&c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb],
				 &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)11);
#line 375 "csytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 379 "csytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 379 "csytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 379 "csytrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 379 "csytrf_aa_2stage.f"
		    cgemm_("Transpose", "NoTranspose", &nb, &i__2, &i__3, &
			    z__1, &work[nb + 1], n, &a[((j + 1) * nb + 1) * 
			    a_dim1 + 1], lda, &c_b2, &a[j * nb + 1 + ((j + 1) 
			    * nb + 1) * a_dim1], lda, (ftnlen)9, (ftnlen)11);
#line 384 "csytrf_aa_2stage.f"
		}

/*              Copy panel to workspace to call CGETRF */

#line 388 "csytrf_aa_2stage.f"
		i__2 = nb;
#line 388 "csytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 389 "csytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 389 "csytrf_aa_2stage.f"
		    ccopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1]
			    , lda, &work[(k - 1) * *n + 1], &c__1);
#line 392 "csytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 396 "csytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 396 "csytrf_aa_2stage.f"
		cgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], &
			iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Copy panel back */

#line 405 "csytrf_aa_2stage.f"
		i__2 = nb;
#line 405 "csytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 406 "csytrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 406 "csytrf_aa_2stage.f"
		    ccopy_(&i__3, &work[(k - 1) * *n + 1], &c__1, &a[j * nb + 
			    k + ((j + 1) * nb + 1) * a_dim1], lda);
#line 409 "csytrf_aa_2stage.f"
		}

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 413 "csytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 413 "csytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 414 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 414 "csytrf_aa_2stage.f"
		claset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 416 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 416 "csytrf_aa_2stage.f"
		clacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)5);
#line 419 "csytrf_aa_2stage.f"
		if (j > 0) {
#line 420 "csytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 420 "csytrf_aa_2stage.f"
		    ctrsm_("R", "U", "N", "U", &kb, &nb, &c_b2, &a[(j - 1) * 
			    nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 423 "csytrf_aa_2stage.f"
		}

/*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
/*              updates */

#line 428 "csytrf_aa_2stage.f"
		i__2 = nb;
#line 428 "csytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 429 "csytrf_aa_2stage.f"
		    i__3 = kb;
#line 429 "csytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 430 "csytrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 430 "csytrf_aa_2stage.f"
			i__5 = td + nb + i__ - k + 1 + (j * nb + k - 1) * 
				ldtb;
#line 430 "csytrf_aa_2stage.f"
			tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 432 "csytrf_aa_2stage.f"
		    }
#line 433 "csytrf_aa_2stage.f"
		}
#line 434 "csytrf_aa_2stage.f"
		claset_("Lower", &kb, &nb, &c_b1, &c_b2, &a[j * nb + 1 + ((j 
			+ 1) * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 439 "csytrf_aa_2stage.f"
		i__2 = kb;
#line 439 "csytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 441 "csytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 443 "csytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 444 "csytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 445 "csytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 447 "csytrf_aa_2stage.f"
			i__3 = k - 1;
#line 447 "csytrf_aa_2stage.f"
			cswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], &
				c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], &
				c__1);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 450 "csytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 450 "csytrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[i1 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 453 "csytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 453 "csytrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[i2 
				+ (i2 + 1) * a_dim1], lda);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 456 "csytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 456 "csytrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 457 "csytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 457 "csytrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 457 "csytrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 458 "csytrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 458 "csytrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 460 "csytrf_aa_2stage.f"
			if (j > 0) {
#line 461 "csytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 461 "csytrf_aa_2stage.f"
			    cswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * 
				    a_dim1 + 1], &c__1);
#line 463 "csytrf_aa_2stage.f"
			}
#line 464 "csytrf_aa_2stage.f"
		    }
#line 465 "csytrf_aa_2stage.f"
		}
#line 466 "csytrf_aa_2stage.f"
	    }
#line 467 "csytrf_aa_2stage.f"
	}
#line 468 "csytrf_aa_2stage.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 474 "csytrf_aa_2stage.f"
	i__1 = nt - 1;
#line 474 "csytrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 478 "csytrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 478 "csytrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 479 "csytrf_aa_2stage.f"
	    i__2 = j - 1;
#line 479 "csytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 480 "csytrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
#line 482 "csytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 483 "csytrf_aa_2stage.f"
			jb = nb + kb;
#line 484 "csytrf_aa_2stage.f"
		    } else {
#line 485 "csytrf_aa_2stage.f"
			jb = nb << 1;
#line 486 "csytrf_aa_2stage.f"
		    }
#line 487 "csytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 487 "csytrf_aa_2stage.f"
		    cgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b2, &
			    tb[td + 1 + i__ * nb * ldtb], &i__3, &a[j * nb + 
			    1 + ((i__ - 1) * nb + 1) * a_dim1], lda, &c_b1, &
			    work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)9);
#line 492 "csytrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
#line 494 "csytrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 495 "csytrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 496 "csytrf_aa_2stage.f"
		    } else {
#line 497 "csytrf_aa_2stage.f"
			jb = nb * 3;
#line 498 "csytrf_aa_2stage.f"
		    }
#line 499 "csytrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 499 "csytrf_aa_2stage.f"
		    cgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b2, &
			    tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, &
			    a[j * nb + 1 + ((i__ - 2) * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)9);
#line 505 "csytrf_aa_2stage.f"
		}
#line 506 "csytrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 510 "csytrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 510 "csytrf_aa_2stage.f"
	    clacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 512 "csytrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = L(J,1:J)*H(1:J) */
#line 514 "csytrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 514 "csytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 514 "csytrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 514 "csytrf_aa_2stage.f"
		cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &z__1, &
			a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &c_b2, 
			&tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)11, (
			ftnlen)11);
/*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
#line 520 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 520 "csytrf_aa_2stage.f"
		cgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[
			j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[
			td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &
			work[1], n, (ftnlen)11, (ftnlen)11);
#line 525 "csytrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 525 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 525 "csytrf_aa_2stage.f"
		cgemm_("NoTranspose", "Transpose", &kb, &kb, &nb, &z__1, &
			work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
			a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)9);
#line 530 "csytrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 534 "csytrf_aa_2stage.f"
	    i__2 = kb;
#line 534 "csytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 535 "csytrf_aa_2stage.f"
		i__3 = kb;
#line 535 "csytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 536 "csytrf_aa_2stage.f"
		    i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 536 "csytrf_aa_2stage.f"
		    i__5 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 536 "csytrf_aa_2stage.f"
		    tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 538 "csytrf_aa_2stage.f"
		}
#line 539 "csytrf_aa_2stage.f"
	    }
#line 540 "csytrf_aa_2stage.f"
	    if (j > 0) {
/*               CALL CHEGST( 1, 'Lower', KB, */
/*     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1, */
/*     $                      A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO ) */
#line 544 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 544 "csytrf_aa_2stage.f"
		ctrsm_("L", "L", "N", "N", &kb, &kb, &c_b2, &a[j * nb + 1 + ((
			j - 1) * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 547 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 547 "csytrf_aa_2stage.f"
		ctrsm_("R", "L", "T", "N", &kb, &kb, &c_b2, &a[j * nb + 1 + ((
			j - 1) * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb *
			 ldtb], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 550 "csytrf_aa_2stage.f"
	    }

/*           Symmetrize T(J,J) */

#line 554 "csytrf_aa_2stage.f"
	    i__2 = kb;
#line 554 "csytrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 555 "csytrf_aa_2stage.f"
		i__3 = kb;
#line 555 "csytrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 556 "csytrf_aa_2stage.f"
		    i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 556 "csytrf_aa_2stage.f"
		    i__5 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 556 "csytrf_aa_2stage.f"
		    tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 558 "csytrf_aa_2stage.f"
		}
#line 559 "csytrf_aa_2stage.f"
	    }

#line 561 "csytrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 562 "csytrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 566 "csytrf_aa_2stage.f"
		    if (j == 1) {
#line 567 "csytrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 567 "csytrf_aa_2stage.f"
			cgemm_("NoTranspose", "Transpose", &kb, &kb, &kb, &
				c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[
				j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
				lda, &c_b1, &work[j * nb + 1], n, (ftnlen)11, 
				(ftnlen)9);
#line 572 "csytrf_aa_2stage.f"
		    } else {
#line 573 "csytrf_aa_2stage.f"
			i__2 = nb + kb;
#line 573 "csytrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 573 "csytrf_aa_2stage.f"
			cgemm_("NoTranspose", "Transpose", &kb, &kb, &i__2, &
				c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb], 
				&i__3, &a[j * nb + 1 + ((j - 2) * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)9);
#line 579 "csytrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 583 "csytrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 583 "csytrf_aa_2stage.f"
		    i__3 = j * nb;
#line 583 "csytrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 583 "csytrf_aa_2stage.f"
		    cgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, &
			    z__1, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[
			    nb + 1], n, &c_b2, &a[(j + 1) * nb + 1 + (j * nb 
			    + 1) * a_dim1], lda, (ftnlen)11, (ftnlen)11);
#line 588 "csytrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 592 "csytrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 592 "csytrf_aa_2stage.f"
		cgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 601 "csytrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 601 "csytrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 602 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 602 "csytrf_aa_2stage.f"
		claset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 604 "csytrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 604 "csytrf_aa_2stage.f"
		clacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) 
			* a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &
			i__2, (ftnlen)5);
#line 607 "csytrf_aa_2stage.f"
		if (j > 0) {
#line 608 "csytrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 608 "csytrf_aa_2stage.f"
		    ctrsm_("R", "L", "T", "U", &kb, &nb, &c_b2, &a[j * nb + 1 
			    + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 611 "csytrf_aa_2stage.f"
		}

/*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
/*              updates */

#line 616 "csytrf_aa_2stage.f"
		i__2 = nb;
#line 616 "csytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 617 "csytrf_aa_2stage.f"
		    i__3 = kb;
#line 617 "csytrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 618 "csytrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 618 "csytrf_aa_2stage.f"
			i__5 = td + nb + i__ - k + 1 + (j * nb + k - 1) * 
				ldtb;
#line 618 "csytrf_aa_2stage.f"
			tb[i__4].r = tb[i__5].r, tb[i__4].i = tb[i__5].i;
#line 620 "csytrf_aa_2stage.f"
		    }
#line 621 "csytrf_aa_2stage.f"
		}
#line 622 "csytrf_aa_2stage.f"
		claset_("Upper", &kb, &nb, &c_b1, &c_b2, &a[(j + 1) * nb + 1 
			+ (j * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 627 "csytrf_aa_2stage.f"
		i__2 = kb;
#line 627 "csytrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 629 "csytrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 631 "csytrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 632 "csytrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 633 "csytrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 635 "csytrf_aa_2stage.f"
			i__3 = k - 1;
#line 635 "csytrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], 
				lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], 
				lda);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 638 "csytrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 638 "csytrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ (i1 + 1) * a_dim1], lda);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 641 "csytrf_aa_2stage.f"
			i__3 = *n - i2;
#line 641 "csytrf_aa_2stage.f"
			cswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 644 "csytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 644 "csytrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 645 "csytrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 645 "csytrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 645 "csytrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 646 "csytrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 646 "csytrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 648 "csytrf_aa_2stage.f"
			if (j > 0) {
#line 649 "csytrf_aa_2stage.f"
			    i__3 = j * nb;
#line 649 "csytrf_aa_2stage.f"
			    cswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + 
				    a_dim1], lda);
#line 651 "csytrf_aa_2stage.f"
			}
#line 652 "csytrf_aa_2stage.f"
		    }
#line 653 "csytrf_aa_2stage.f"
		}

/*              Apply pivots to previous columns of L */

/*               CALL CLASWP( J*NB, A( 1, 1 ), LDA, */
/*     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
#line 659 "csytrf_aa_2stage.f"
	    }
#line 660 "csytrf_aa_2stage.f"
	}
#line 661 "csytrf_aa_2stage.f"
    }

/*     Factor the band matrix */
#line 664 "csytrf_aa_2stage.f"
    cgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);

/*     End of CSYTRF_AA_2STAGE */

#line 668 "csytrf_aa_2stage.f"
    return 0;
} /* csytrf_aa_2stage__ */

