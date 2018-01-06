#line 1 "zhetrf_aa_2stage.f"
/* zhetrf_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrf_aa_2stage.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZHETRF_AA_2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRF_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrf_
aa_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrf_
aa_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrf_
aa_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE ZHETRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
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
/* > ZHETRF_AA_2STAGE computes the factorization of a double hermitian matrix A */
/* > using the Aasen's algorithm.  The form of the factorization is */
/* > */
/* >    A = U*T*U**T  or  A = L*T*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a hermitian band matrix with the */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhetrf_aa_2stage__(char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublecomplex *tb, integer *ltb, integer *ipiv, 
	integer *ipiv2, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    , zgbtrf_(integer *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *), zgetrf_(
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    integer *), zlacpy_(char *, integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, ftnlen), zlaset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), zhegst_(integer *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
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

#line 206 "zhetrf_aa_2stage.f"
    /* Parameter adjustments */
#line 206 "zhetrf_aa_2stage.f"
    a_dim1 = *lda;
#line 206 "zhetrf_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 206 "zhetrf_aa_2stage.f"
    a -= a_offset;
#line 206 "zhetrf_aa_2stage.f"
    --tb;
#line 206 "zhetrf_aa_2stage.f"
    --ipiv;
#line 206 "zhetrf_aa_2stage.f"
    --ipiv2;
#line 206 "zhetrf_aa_2stage.f"
    --work;
#line 206 "zhetrf_aa_2stage.f"

#line 206 "zhetrf_aa_2stage.f"
    /* Function Body */
#line 206 "zhetrf_aa_2stage.f"
    *info = 0;
#line 207 "zhetrf_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 208 "zhetrf_aa_2stage.f"
    wquery = *lwork == -1;
#line 209 "zhetrf_aa_2stage.f"
    tquery = *ltb == -1;
#line 210 "zhetrf_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 211 "zhetrf_aa_2stage.f"
	*info = -1;
#line 212 "zhetrf_aa_2stage.f"
    } else if (*n < 0) {
#line 213 "zhetrf_aa_2stage.f"
	*info = -2;
#line 214 "zhetrf_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 215 "zhetrf_aa_2stage.f"
	*info = -4;
#line 216 "zhetrf_aa_2stage.f"
    } else if (*ltb < *n << 2 && ! tquery) {
#line 217 "zhetrf_aa_2stage.f"
	*info = -6;
#line 218 "zhetrf_aa_2stage.f"
    } else if (*lwork < *n && ! wquery) {
#line 219 "zhetrf_aa_2stage.f"
	*info = -10;
#line 220 "zhetrf_aa_2stage.f"
    }

#line 222 "zhetrf_aa_2stage.f"
    if (*info != 0) {
#line 223 "zhetrf_aa_2stage.f"
	i__1 = -(*info);
#line 223 "zhetrf_aa_2stage.f"
	xerbla_("ZHETRF_AA_2STAGE", &i__1, (ftnlen)16);
#line 224 "zhetrf_aa_2stage.f"
	return 0;
#line 225 "zhetrf_aa_2stage.f"
    }

/*     Answer the query */

#line 229 "zhetrf_aa_2stage.f"
    nb = ilaenv_(&c__1, "ZHETRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)16, (ftnlen)1);
#line 230 "zhetrf_aa_2stage.f"
    if (*info == 0) {
#line 231 "zhetrf_aa_2stage.f"
	if (tquery) {
#line 232 "zhetrf_aa_2stage.f"
	    i__1 = (nb * 3 + 1) * *n;
#line 232 "zhetrf_aa_2stage.f"
	    tb[1].r = (doublereal) i__1, tb[1].i = 0.;
#line 233 "zhetrf_aa_2stage.f"
	}
#line 234 "zhetrf_aa_2stage.f"
	if (wquery) {
#line 235 "zhetrf_aa_2stage.f"
	    i__1 = *n * nb;
#line 235 "zhetrf_aa_2stage.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 236 "zhetrf_aa_2stage.f"
	}
#line 237 "zhetrf_aa_2stage.f"
    }
#line 238 "zhetrf_aa_2stage.f"
    if (tquery || wquery) {
#line 239 "zhetrf_aa_2stage.f"
	return 0;
#line 240 "zhetrf_aa_2stage.f"
    }

/*     Quick return */

#line 244 "zhetrf_aa_2stage.f"
    if (*n == 0) {
#line 245 "zhetrf_aa_2stage.f"
	return 0;
#line 246 "zhetrf_aa_2stage.f"
    }

/*     Determine the number of the block size */

#line 250 "zhetrf_aa_2stage.f"
    ldtb = *ltb / *n;
#line 251 "zhetrf_aa_2stage.f"
    if (ldtb < nb * 3 + 1) {
#line 252 "zhetrf_aa_2stage.f"
	nb = (ldtb - 1) / 3;
#line 253 "zhetrf_aa_2stage.f"
    }
#line 254 "zhetrf_aa_2stage.f"
    if (*lwork < nb * *n) {
#line 255 "zhetrf_aa_2stage.f"
	nb = *lwork / *n;
#line 256 "zhetrf_aa_2stage.f"
    }

/*     Determine the number of the block columns */

#line 260 "zhetrf_aa_2stage.f"
    nt = (*n + nb - 1) / nb;
#line 261 "zhetrf_aa_2stage.f"
    td = nb << 1;
#line 262 "zhetrf_aa_2stage.f"
    kb = min(nb,*n);

/*     Initialize vectors/matrices */

#line 266 "zhetrf_aa_2stage.f"
    i__1 = kb;
#line 266 "zhetrf_aa_2stage.f"
    for (j = 1; j <= i__1; ++j) {
#line 267 "zhetrf_aa_2stage.f"
	ipiv[j] = j;
#line 268 "zhetrf_aa_2stage.f"
    }

/*     Save NB */

#line 272 "zhetrf_aa_2stage.f"
    tb[1].r = (doublereal) nb, tb[1].i = 0.;

#line 274 "zhetrf_aa_2stage.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

#line 280 "zhetrf_aa_2stage.f"
	i__1 = nt - 1;
#line 280 "zhetrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 284 "zhetrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 284 "zhetrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 285 "zhetrf_aa_2stage.f"
	    i__2 = j - 1;
#line 285 "zhetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 286 "zhetrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J) */
#line 288 "zhetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 289 "zhetrf_aa_2stage.f"
			jb = nb + kb;
#line 290 "zhetrf_aa_2stage.f"
		    } else {
#line 291 "zhetrf_aa_2stage.f"
			jb = nb << 1;
#line 292 "zhetrf_aa_2stage.f"
		    }
#line 293 "zhetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 293 "zhetrf_aa_2stage.f"
		    zgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(i__ - 
			    1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, 
			    &work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)11);
#line 298 "zhetrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 300 "zhetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 301 "zhetrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 302 "zhetrf_aa_2stage.f"
		    } else {
#line 303 "zhetrf_aa_2stage.f"
			jb = nb * 3;
#line 304 "zhetrf_aa_2stage.f"
		    }
#line 305 "zhetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 305 "zhetrf_aa_2stage.f"
		    zgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, 
			    &a[(i__ - 2) * nb + 1 + (j * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)11);
#line 311 "zhetrf_aa_2stage.f"
		}
#line 312 "zhetrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 316 "zhetrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 316 "zhetrf_aa_2stage.f"
	    zlacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 318 "zhetrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = U(1:J,J)'*H(1:J) */
#line 320 "zhetrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 320 "zhetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 320 "zhetrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 320 "zhetrf_aa_2stage.f"
		zgemm_("Conjugate transpose", "NoTranspose", &kb, &kb, &i__2, 
			&z__1, &a[(j * nb + 1) * a_dim1 + 1], lda, &work[nb + 
			1], n, &c_b2, &tb[td + 1 + j * nb * ldtb], &i__3, (
			ftnlen)19, (ftnlen)11);
/*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
#line 326 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 326 "zhetrf_aa_2stage.f"
		zgemm_("Conjugate transpose", "NoTranspose", &kb, &nb, &kb, &
			c_b2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
			lda, &tb[td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &
			c_b1, &work[1], n, (ftnlen)19, (ftnlen)11);
#line 331 "zhetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 331 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 331 "zhetrf_aa_2stage.f"
		zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &z__1, &
			work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)11);
#line 336 "zhetrf_aa_2stage.f"
	    }
#line 337 "zhetrf_aa_2stage.f"
	    if (j > 0) {
#line 338 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 338 "zhetrf_aa_2stage.f"
		zhegst_(&c__1, "Upper", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 341 "zhetrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 345 "zhetrf_aa_2stage.f"
	    i__2 = kb;
#line 345 "zhetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 346 "zhetrf_aa_2stage.f"
		i__3 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 346 "zhetrf_aa_2stage.f"
		i__4 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 346 "zhetrf_aa_2stage.f"
		d__1 = tb[i__4].r;
#line 346 "zhetrf_aa_2stage.f"
		tb[i__3].r = d__1, tb[i__3].i = 0.;
#line 348 "zhetrf_aa_2stage.f"
		i__3 = kb;
#line 348 "zhetrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 349 "zhetrf_aa_2stage.f"
		    i__4 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 349 "zhetrf_aa_2stage.f"
		    d_cnjg(&z__1, &tb[td - (k - (i__ + 1)) + (j * nb + k - 1) 
			    * ldtb]);
#line 349 "zhetrf_aa_2stage.f"
		    tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 351 "zhetrf_aa_2stage.f"
		}
#line 352 "zhetrf_aa_2stage.f"
	    }

#line 354 "zhetrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 355 "zhetrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 359 "zhetrf_aa_2stage.f"
		    if (j == 1) {
#line 360 "zhetrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 360 "zhetrf_aa_2stage.f"
			zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, &
				c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[(
				j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda,
				 &c_b1, &work[j * nb + 1], n, (ftnlen)11, (
				ftnlen)11);
#line 365 "zhetrf_aa_2stage.f"
		    } else {
#line 366 "zhetrf_aa_2stage.f"
			i__2 = nb + kb;
#line 366 "zhetrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 366 "zhetrf_aa_2stage.f"
			zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, 
				&c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb],
				 &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)11);
#line 372 "zhetrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 376 "zhetrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 376 "zhetrf_aa_2stage.f"
		    i__3 = j * nb;
#line 376 "zhetrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 376 "zhetrf_aa_2stage.f"
		    zgemm_("Conjugate transpose", "NoTranspose", &nb, &i__2, &
			    i__3, &z__1, &work[nb + 1], n, &a[((j + 1) * nb + 
			    1) * a_dim1 + 1], lda, &c_b2, &a[j * nb + 1 + ((j 
			    + 1) * nb + 1) * a_dim1], lda, (ftnlen)19, (
			    ftnlen)11);
#line 381 "zhetrf_aa_2stage.f"
		}

/*              Copy panel to workspace to call ZGETRF */

#line 385 "zhetrf_aa_2stage.f"
		i__2 = nb;
#line 385 "zhetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 386 "zhetrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 386 "zhetrf_aa_2stage.f"
		    zcopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1]
			    , lda, &work[(k - 1) * *n + 1], &c__1);
#line 389 "zhetrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 393 "zhetrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 393 "zhetrf_aa_2stage.f"
		zgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], &
			iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Copy panel back */

#line 402 "zhetrf_aa_2stage.f"
		i__2 = nb;
#line 402 "zhetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {

/*                  Copy only L-factor */

#line 406 "zhetrf_aa_2stage.f"
		    i__3 = *n - k - (j + 1) * nb;
#line 406 "zhetrf_aa_2stage.f"
		    zcopy_(&i__3, &work[k + 1 + (k - 1) * *n], &c__1, &a[j * 
			    nb + k + ((j + 1) * nb + k + 1) * a_dim1], lda);

/*                  Transpose U-factor to be copied back into T(J+1, J) */

#line 412 "zhetrf_aa_2stage.f"
		    zlacgv_(&k, &work[(k - 1) * *n + 1], &c__1);
#line 413 "zhetrf_aa_2stage.f"
		}

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 417 "zhetrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 417 "zhetrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 418 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 418 "zhetrf_aa_2stage.f"
		zlaset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 420 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 420 "zhetrf_aa_2stage.f"
		zlacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)5);
#line 423 "zhetrf_aa_2stage.f"
		if (j > 0) {
#line 424 "zhetrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 424 "zhetrf_aa_2stage.f"
		    ztrsm_("R", "U", "N", "U", &kb, &nb, &c_b2, &a[(j - 1) * 
			    nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 427 "zhetrf_aa_2stage.f"
		}

/*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
/*              updates */

#line 432 "zhetrf_aa_2stage.f"
		i__2 = nb;
#line 432 "zhetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 433 "zhetrf_aa_2stage.f"
		    i__3 = kb;
#line 433 "zhetrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 434 "zhetrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 434 "zhetrf_aa_2stage.f"
			d_cnjg(&z__1, &tb[td + nb + i__ - k + 1 + (j * nb + k 
				- 1) * ldtb]);
#line 434 "zhetrf_aa_2stage.f"
			tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 436 "zhetrf_aa_2stage.f"
		    }
#line 437 "zhetrf_aa_2stage.f"
		}
#line 438 "zhetrf_aa_2stage.f"
		zlaset_("Lower", &kb, &nb, &c_b1, &c_b2, &a[j * nb + 1 + ((j 
			+ 1) * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 443 "zhetrf_aa_2stage.f"
		i__2 = kb;
#line 443 "zhetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 445 "zhetrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 447 "zhetrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 448 "zhetrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 449 "zhetrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 451 "zhetrf_aa_2stage.f"
			i__3 = k - 1;
#line 451 "zhetrf_aa_2stage.f"
			zswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], &
				c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], &
				c__1);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 454 "zhetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 454 "zhetrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[i1 
				+ 1 + i2 * a_dim1], &c__1);
#line 456 "zhetrf_aa_2stage.f"
			i__3 = i2 - i1;
#line 456 "zhetrf_aa_2stage.f"
			zlacgv_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda);
#line 457 "zhetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 457 "zhetrf_aa_2stage.f"
			zlacgv_(&i__3, &a[i1 + 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 459 "zhetrf_aa_2stage.f"
			i__3 = *n - i2;
#line 459 "zhetrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[i2 
				+ (i2 + 1) * a_dim1], lda);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 462 "zhetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 462 "zhetrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 463 "zhetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 463 "zhetrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 463 "zhetrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 464 "zhetrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 464 "zhetrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 466 "zhetrf_aa_2stage.f"
			if (j > 0) {
#line 467 "zhetrf_aa_2stage.f"
			    i__3 = j * nb;
#line 467 "zhetrf_aa_2stage.f"
			    zswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * 
				    a_dim1 + 1], &c__1);
#line 469 "zhetrf_aa_2stage.f"
			}
#line 470 "zhetrf_aa_2stage.f"
		    }
#line 471 "zhetrf_aa_2stage.f"
		}
#line 472 "zhetrf_aa_2stage.f"
	    }
#line 473 "zhetrf_aa_2stage.f"
	}
#line 474 "zhetrf_aa_2stage.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 480 "zhetrf_aa_2stage.f"
	i__1 = nt - 1;
#line 480 "zhetrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 484 "zhetrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 484 "zhetrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 485 "zhetrf_aa_2stage.f"
	    i__2 = j - 1;
#line 485 "zhetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 486 "zhetrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
#line 488 "zhetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 489 "zhetrf_aa_2stage.f"
			jb = nb + kb;
#line 490 "zhetrf_aa_2stage.f"
		    } else {
#line 491 "zhetrf_aa_2stage.f"
			jb = nb << 1;
#line 492 "zhetrf_aa_2stage.f"
		    }
#line 493 "zhetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 493 "zhetrf_aa_2stage.f"
		    zgemm_("NoTranspose", "Conjugate transpose", &nb, &kb, &
			    jb, &c_b2, &tb[td + 1 + i__ * nb * ldtb], &i__3, &
			    a[j * nb + 1 + ((i__ - 1) * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)19);
#line 498 "zhetrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
#line 500 "zhetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 501 "zhetrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 502 "zhetrf_aa_2stage.f"
		    } else {
#line 503 "zhetrf_aa_2stage.f"
			jb = nb * 3;
#line 504 "zhetrf_aa_2stage.f"
		    }
#line 505 "zhetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 505 "zhetrf_aa_2stage.f"
		    zgemm_("NoTranspose", "Conjugate transpose", &nb, &kb, &
			    jb, &c_b2, &tb[td + nb + 1 + (i__ - 1) * nb * 
			    ldtb], &i__3, &a[j * nb + 1 + ((i__ - 2) * nb + 1)
			     * a_dim1], lda, &c_b1, &work[i__ * nb + 1], n, (
			    ftnlen)11, (ftnlen)19);
#line 511 "zhetrf_aa_2stage.f"
		}
#line 512 "zhetrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 516 "zhetrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 516 "zhetrf_aa_2stage.f"
	    zlacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 518 "zhetrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = L(J,1:J)*H(1:J) */
#line 520 "zhetrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 520 "zhetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 520 "zhetrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 520 "zhetrf_aa_2stage.f"
		zgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &z__1, &
			a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &c_b2, 
			&tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)11, (
			ftnlen)11);
/*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
#line 526 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 526 "zhetrf_aa_2stage.f"
		zgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[
			j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[
			td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &
			work[1], n, (ftnlen)11, (ftnlen)11);
#line 531 "zhetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 531 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 531 "zhetrf_aa_2stage.f"
		zgemm_("NoTranspose", "Conjugate transpose", &kb, &kb, &nb, &
			z__1, &work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) 
			* a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)19);
#line 536 "zhetrf_aa_2stage.f"
	    }
#line 537 "zhetrf_aa_2stage.f"
	    if (j > 0) {
#line 538 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 538 "zhetrf_aa_2stage.f"
		zhegst_(&c__1, "Lower", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 541 "zhetrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 545 "zhetrf_aa_2stage.f"
	    i__2 = kb;
#line 545 "zhetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 546 "zhetrf_aa_2stage.f"
		i__3 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 546 "zhetrf_aa_2stage.f"
		i__4 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 546 "zhetrf_aa_2stage.f"
		d__1 = tb[i__4].r;
#line 546 "zhetrf_aa_2stage.f"
		tb[i__3].r = d__1, tb[i__3].i = 0.;
#line 548 "zhetrf_aa_2stage.f"
		i__3 = kb;
#line 548 "zhetrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 549 "zhetrf_aa_2stage.f"
		    i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 549 "zhetrf_aa_2stage.f"
		    d_cnjg(&z__1, &tb[td + (k - i__) + 1 + (j * nb + i__ - 1) 
			    * ldtb]);
#line 549 "zhetrf_aa_2stage.f"
		    tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 551 "zhetrf_aa_2stage.f"
		}
#line 552 "zhetrf_aa_2stage.f"
	    }

#line 554 "zhetrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 555 "zhetrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 559 "zhetrf_aa_2stage.f"
		    if (j == 1) {
#line 560 "zhetrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 560 "zhetrf_aa_2stage.f"
			zgemm_("NoTranspose", "Conjugate transpose", &kb, &kb,
				 &kb, &c_b2, &tb[td + 1 + j * nb * ldtb], &
				i__2, &a[j * nb + 1 + ((j - 1) * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)19);
#line 565 "zhetrf_aa_2stage.f"
		    } else {
#line 566 "zhetrf_aa_2stage.f"
			i__2 = nb + kb;
#line 566 "zhetrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 566 "zhetrf_aa_2stage.f"
			zgemm_("NoTranspose", "Conjugate transpose", &kb, &kb,
				 &i__2, &c_b2, &tb[td + nb + 1 + (j - 1) * nb 
				* ldtb], &i__3, &a[j * nb + 1 + ((j - 2) * nb 
				+ 1) * a_dim1], lda, &c_b1, &work[j * nb + 1],
				 n, (ftnlen)11, (ftnlen)19);
#line 572 "zhetrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 576 "zhetrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 576 "zhetrf_aa_2stage.f"
		    i__3 = j * nb;
#line 576 "zhetrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 576 "zhetrf_aa_2stage.f"
		    zgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, &
			    z__1, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[
			    nb + 1], n, &c_b2, &a[(j + 1) * nb + 1 + (j * nb 
			    + 1) * a_dim1], lda, (ftnlen)11, (ftnlen)11);
#line 581 "zhetrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 585 "zhetrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 585 "zhetrf_aa_2stage.f"
		zgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 594 "zhetrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 594 "zhetrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 595 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 595 "zhetrf_aa_2stage.f"
		zlaset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 597 "zhetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 597 "zhetrf_aa_2stage.f"
		zlacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) 
			* a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &
			i__2, (ftnlen)5);
#line 600 "zhetrf_aa_2stage.f"
		if (j > 0) {
#line 601 "zhetrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 601 "zhetrf_aa_2stage.f"
		    ztrsm_("R", "L", "C", "U", &kb, &nb, &c_b2, &a[j * nb + 1 
			    + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 604 "zhetrf_aa_2stage.f"
		}

/*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
/*              updates */

#line 609 "zhetrf_aa_2stage.f"
		i__2 = nb;
#line 609 "zhetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 610 "zhetrf_aa_2stage.f"
		    i__3 = kb;
#line 610 "zhetrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 611 "zhetrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 611 "zhetrf_aa_2stage.f"
			d_cnjg(&z__1, &tb[td + nb + i__ - k + 1 + (j * nb + k 
				- 1) * ldtb]);
#line 611 "zhetrf_aa_2stage.f"
			tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 613 "zhetrf_aa_2stage.f"
		    }
#line 614 "zhetrf_aa_2stage.f"
		}
#line 615 "zhetrf_aa_2stage.f"
		zlaset_("Upper", &kb, &nb, &c_b1, &c_b2, &a[(j + 1) * nb + 1 
			+ (j * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 620 "zhetrf_aa_2stage.f"
		i__2 = kb;
#line 620 "zhetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 622 "zhetrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 624 "zhetrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 625 "zhetrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 626 "zhetrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 628 "zhetrf_aa_2stage.f"
			i__3 = k - 1;
#line 628 "zhetrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], 
				lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], 
				lda);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 631 "zhetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 631 "zhetrf_aa_2stage.f"
			zswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ (i1 + 1) * a_dim1], lda);
#line 633 "zhetrf_aa_2stage.f"
			i__3 = i2 - i1;
#line 633 "zhetrf_aa_2stage.f"
			zlacgv_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1);
#line 634 "zhetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 634 "zhetrf_aa_2stage.f"
			zlacgv_(&i__3, &a[i2 + (i1 + 1) * a_dim1], lda);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 636 "zhetrf_aa_2stage.f"
			i__3 = *n - i2;
#line 636 "zhetrf_aa_2stage.f"
			zswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 639 "zhetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 639 "zhetrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 640 "zhetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 640 "zhetrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 640 "zhetrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 641 "zhetrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 641 "zhetrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 643 "zhetrf_aa_2stage.f"
			if (j > 0) {
#line 644 "zhetrf_aa_2stage.f"
			    i__3 = j * nb;
#line 644 "zhetrf_aa_2stage.f"
			    zswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + 
				    a_dim1], lda);
#line 646 "zhetrf_aa_2stage.f"
			}
#line 647 "zhetrf_aa_2stage.f"
		    }
#line 648 "zhetrf_aa_2stage.f"
		}

/*              Apply pivots to previous columns of L */

/*               CALL ZLASWP( J*NB, A( 1, 1 ), LDA, */
/*     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
#line 654 "zhetrf_aa_2stage.f"
	    }
#line 655 "zhetrf_aa_2stage.f"
	}
#line 656 "zhetrf_aa_2stage.f"
    }

/*     Factor the band matrix */
#line 659 "zhetrf_aa_2stage.f"
    zgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);

/*     End of ZHETRF_AA_2STAGE */

#line 663 "zhetrf_aa_2stage.f"
    return 0;
} /* zhetrf_aa_2stage__ */

