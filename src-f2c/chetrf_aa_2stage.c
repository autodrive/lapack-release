#line 1 "chetrf_aa_2stage.f"
/* chetrf_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "chetrf_aa_2stage.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CHETRF_AA_2STAGE */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRF_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf_
aa_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf_
aa_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf_
aa_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE CHETRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
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
/* > CHETRF_AA_2STAGE computes the factorization of a real hermitian matrix A */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int chetrf_aa_2stage__(char *uplo, integer *n, doublecomplex 
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
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , cgbtrf_(integer *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *), cgetrf_(
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    integer *), clacpy_(char *, integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, ftnlen), claset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int chegst_(integer *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);
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

#line 207 "chetrf_aa_2stage.f"
    /* Parameter adjustments */
#line 207 "chetrf_aa_2stage.f"
    a_dim1 = *lda;
#line 207 "chetrf_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 207 "chetrf_aa_2stage.f"
    a -= a_offset;
#line 207 "chetrf_aa_2stage.f"
    --tb;
#line 207 "chetrf_aa_2stage.f"
    --ipiv;
#line 207 "chetrf_aa_2stage.f"
    --ipiv2;
#line 207 "chetrf_aa_2stage.f"
    --work;
#line 207 "chetrf_aa_2stage.f"

#line 207 "chetrf_aa_2stage.f"
    /* Function Body */
#line 207 "chetrf_aa_2stage.f"
    *info = 0;
#line 208 "chetrf_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 209 "chetrf_aa_2stage.f"
    wquery = *lwork == -1;
#line 210 "chetrf_aa_2stage.f"
    tquery = *ltb == -1;
#line 211 "chetrf_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 212 "chetrf_aa_2stage.f"
	*info = -1;
#line 213 "chetrf_aa_2stage.f"
    } else if (*n < 0) {
#line 214 "chetrf_aa_2stage.f"
	*info = -2;
#line 215 "chetrf_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 216 "chetrf_aa_2stage.f"
	*info = -4;
#line 217 "chetrf_aa_2stage.f"
    } else if (*ltb < *n << 2 && ! tquery) {
#line 218 "chetrf_aa_2stage.f"
	*info = -6;
#line 219 "chetrf_aa_2stage.f"
    } else if (*lwork < *n && ! wquery) {
#line 220 "chetrf_aa_2stage.f"
	*info = -10;
#line 221 "chetrf_aa_2stage.f"
    }

#line 223 "chetrf_aa_2stage.f"
    if (*info != 0) {
#line 224 "chetrf_aa_2stage.f"
	i__1 = -(*info);
#line 224 "chetrf_aa_2stage.f"
	xerbla_("CHETRF_AA_2STAGE", &i__1, (ftnlen)16);
#line 225 "chetrf_aa_2stage.f"
	return 0;
#line 226 "chetrf_aa_2stage.f"
    }

/*     Answer the query */

#line 230 "chetrf_aa_2stage.f"
    nb = ilaenv_(&c__1, "CHETRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1, (
	    ftnlen)16, (ftnlen)1);
#line 231 "chetrf_aa_2stage.f"
    if (*info == 0) {
#line 232 "chetrf_aa_2stage.f"
	if (tquery) {
#line 233 "chetrf_aa_2stage.f"
	    i__1 = (nb * 3 + 1) * *n;
#line 233 "chetrf_aa_2stage.f"
	    tb[1].r = (doublereal) i__1, tb[1].i = 0.;
#line 234 "chetrf_aa_2stage.f"
	}
#line 235 "chetrf_aa_2stage.f"
	if (wquery) {
#line 236 "chetrf_aa_2stage.f"
	    i__1 = *n * nb;
#line 236 "chetrf_aa_2stage.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 237 "chetrf_aa_2stage.f"
	}
#line 238 "chetrf_aa_2stage.f"
    }
#line 239 "chetrf_aa_2stage.f"
    if (tquery || wquery) {
#line 240 "chetrf_aa_2stage.f"
	return 0;
#line 241 "chetrf_aa_2stage.f"
    }

/*     Quick return */

#line 245 "chetrf_aa_2stage.f"
    if (*n == 0) {
#line 246 "chetrf_aa_2stage.f"
	return 0;
#line 247 "chetrf_aa_2stage.f"
    }

/*     Determine the number of the block size */

#line 251 "chetrf_aa_2stage.f"
    ldtb = *ltb / *n;
#line 252 "chetrf_aa_2stage.f"
    if (ldtb < nb * 3 + 1) {
#line 253 "chetrf_aa_2stage.f"
	nb = (ldtb - 1) / 3;
#line 254 "chetrf_aa_2stage.f"
    }
#line 255 "chetrf_aa_2stage.f"
    if (*lwork < nb * *n) {
#line 256 "chetrf_aa_2stage.f"
	nb = *lwork / *n;
#line 257 "chetrf_aa_2stage.f"
    }

/*     Determine the number of the block columns */

#line 261 "chetrf_aa_2stage.f"
    nt = (*n + nb - 1) / nb;
#line 262 "chetrf_aa_2stage.f"
    td = nb << 1;
#line 263 "chetrf_aa_2stage.f"
    kb = min(nb,*n);

/*     Initialize vectors/matrices */

#line 267 "chetrf_aa_2stage.f"
    i__1 = kb;
#line 267 "chetrf_aa_2stage.f"
    for (j = 1; j <= i__1; ++j) {
#line 268 "chetrf_aa_2stage.f"
	ipiv[j] = j;
#line 269 "chetrf_aa_2stage.f"
    }

/*     Save NB */

#line 273 "chetrf_aa_2stage.f"
    tb[1].r = (doublereal) nb, tb[1].i = 0.;

#line 275 "chetrf_aa_2stage.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

#line 281 "chetrf_aa_2stage.f"
	i__1 = nt - 1;
#line 281 "chetrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 285 "chetrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 285 "chetrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 286 "chetrf_aa_2stage.f"
	    i__2 = j - 1;
#line 286 "chetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 287 "chetrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J) */
#line 289 "chetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 290 "chetrf_aa_2stage.f"
			jb = nb + kb;
#line 291 "chetrf_aa_2stage.f"
		    } else {
#line 292 "chetrf_aa_2stage.f"
			jb = nb << 1;
#line 293 "chetrf_aa_2stage.f"
		    }
#line 294 "chetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 294 "chetrf_aa_2stage.f"
		    cgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(i__ - 
			    1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, 
			    &work[i__ * nb + 1], n, (ftnlen)11, (ftnlen)11);
#line 299 "chetrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
#line 301 "chetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 302 "chetrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 303 "chetrf_aa_2stage.f"
		    } else {
#line 304 "chetrf_aa_2stage.f"
			jb = nb * 3;
#line 305 "chetrf_aa_2stage.f"
		    }
#line 306 "chetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 306 "chetrf_aa_2stage.f"
		    cgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2,
			     &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, 
			    &a[(i__ - 2) * nb + 1 + (j * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)11);
#line 312 "chetrf_aa_2stage.f"
		}
#line 313 "chetrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 317 "chetrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 317 "chetrf_aa_2stage.f"
	    clacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 319 "chetrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = U(1:J,J)'*H(1:J) */
#line 321 "chetrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 321 "chetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 321 "chetrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 321 "chetrf_aa_2stage.f"
		cgemm_("Conjugate transpose", "NoTranspose", &kb, &kb, &i__2, 
			&z__1, &a[(j * nb + 1) * a_dim1 + 1], lda, &work[nb + 
			1], n, &c_b2, &tb[td + 1 + j * nb * ldtb], &i__3, (
			ftnlen)19, (ftnlen)11);
/*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
#line 327 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 327 "chetrf_aa_2stage.f"
		cgemm_("Conjugate transpose", "NoTranspose", &kb, &nb, &kb, &
			c_b2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
			lda, &tb[td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &
			c_b1, &work[1], n, (ftnlen)19, (ftnlen)11);
#line 332 "chetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 332 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 332 "chetrf_aa_2stage.f"
		cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &z__1, &
			work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)11);
#line 337 "chetrf_aa_2stage.f"
	    }
#line 338 "chetrf_aa_2stage.f"
	    if (j > 0) {
#line 339 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 339 "chetrf_aa_2stage.f"
		chegst_(&c__1, "Upper", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 342 "chetrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 346 "chetrf_aa_2stage.f"
	    i__2 = kb;
#line 346 "chetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 347 "chetrf_aa_2stage.f"
		i__3 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 347 "chetrf_aa_2stage.f"
		i__4 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 347 "chetrf_aa_2stage.f"
		d__1 = tb[i__4].r;
#line 347 "chetrf_aa_2stage.f"
		tb[i__3].r = d__1, tb[i__3].i = 0.;
#line 349 "chetrf_aa_2stage.f"
		i__3 = kb;
#line 349 "chetrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 350 "chetrf_aa_2stage.f"
		    i__4 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
#line 350 "chetrf_aa_2stage.f"
		    d_cnjg(&z__1, &tb[td - (k - (i__ + 1)) + (j * nb + k - 1) 
			    * ldtb]);
#line 350 "chetrf_aa_2stage.f"
		    tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 352 "chetrf_aa_2stage.f"
		}
#line 353 "chetrf_aa_2stage.f"
	    }

#line 355 "chetrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 356 "chetrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 360 "chetrf_aa_2stage.f"
		    if (j == 1) {
#line 361 "chetrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 361 "chetrf_aa_2stage.f"
			cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, &
				c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[(
				j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda,
				 &c_b1, &work[j * nb + 1], n, (ftnlen)11, (
				ftnlen)11);
#line 366 "chetrf_aa_2stage.f"
		    } else {
#line 367 "chetrf_aa_2stage.f"
			i__2 = nb + kb;
#line 367 "chetrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 367 "chetrf_aa_2stage.f"
			cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, 
				&c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb],
				 &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)11);
#line 373 "chetrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 377 "chetrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 377 "chetrf_aa_2stage.f"
		    i__3 = j * nb;
#line 377 "chetrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 377 "chetrf_aa_2stage.f"
		    cgemm_("Conjugate transpose", "NoTranspose", &nb, &i__2, &
			    i__3, &z__1, &work[nb + 1], n, &a[((j + 1) * nb + 
			    1) * a_dim1 + 1], lda, &c_b2, &a[j * nb + 1 + ((j 
			    + 1) * nb + 1) * a_dim1], lda, (ftnlen)19, (
			    ftnlen)11);
#line 382 "chetrf_aa_2stage.f"
		}

/*              Copy panel to workspace to call CGETRF */

#line 386 "chetrf_aa_2stage.f"
		i__2 = nb;
#line 386 "chetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 387 "chetrf_aa_2stage.f"
		    i__3 = *n - (j + 1) * nb;
#line 387 "chetrf_aa_2stage.f"
		    ccopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1]
			    , lda, &work[(k - 1) * *n + 1], &c__1);
#line 390 "chetrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 394 "chetrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 394 "chetrf_aa_2stage.f"
		cgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], &
			iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Copy panel back */

#line 403 "chetrf_aa_2stage.f"
		i__2 = nb;
#line 403 "chetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {

/*                  Copy only L-factor */

#line 407 "chetrf_aa_2stage.f"
		    i__3 = *n - k - (j + 1) * nb;
#line 407 "chetrf_aa_2stage.f"
		    ccopy_(&i__3, &work[k + 1 + (k - 1) * *n], &c__1, &a[j * 
			    nb + k + ((j + 1) * nb + k + 1) * a_dim1], lda);

/*                  Transpose U-factor to be copied back into T(J+1, J) */

#line 413 "chetrf_aa_2stage.f"
		    clacgv_(&k, &work[(k - 1) * *n + 1], &c__1);
#line 414 "chetrf_aa_2stage.f"
		}

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 418 "chetrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 418 "chetrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 419 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 419 "chetrf_aa_2stage.f"
		claset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 421 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 421 "chetrf_aa_2stage.f"
		clacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)5);
#line 424 "chetrf_aa_2stage.f"
		if (j > 0) {
#line 425 "chetrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 425 "chetrf_aa_2stage.f"
		    ctrsm_("R", "U", "N", "U", &kb, &nb, &c_b2, &a[(j - 1) * 
			    nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 428 "chetrf_aa_2stage.f"
		}

/*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
/*              updates */

#line 433 "chetrf_aa_2stage.f"
		i__2 = nb;
#line 433 "chetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 434 "chetrf_aa_2stage.f"
		    i__3 = kb;
#line 434 "chetrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 435 "chetrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 435 "chetrf_aa_2stage.f"
			d_cnjg(&z__1, &tb[td + nb + i__ - k + 1 + (j * nb + k 
				- 1) * ldtb]);
#line 435 "chetrf_aa_2stage.f"
			tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 437 "chetrf_aa_2stage.f"
		    }
#line 438 "chetrf_aa_2stage.f"
		}
#line 439 "chetrf_aa_2stage.f"
		claset_("Lower", &kb, &nb, &c_b1, &c_b2, &a[j * nb + 1 + ((j 
			+ 1) * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 444 "chetrf_aa_2stage.f"
		i__2 = kb;
#line 444 "chetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 446 "chetrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 448 "chetrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 449 "chetrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 450 "chetrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 452 "chetrf_aa_2stage.f"
			i__3 = k - 1;
#line 452 "chetrf_aa_2stage.f"
			cswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], &
				c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], &
				c__1);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 455 "chetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 455 "chetrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[i1 
				+ 1 + i2 * a_dim1], &c__1);
#line 457 "chetrf_aa_2stage.f"
			i__3 = i2 - i1;
#line 457 "chetrf_aa_2stage.f"
			clacgv_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda);
#line 458 "chetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 458 "chetrf_aa_2stage.f"
			clacgv_(&i__3, &a[i1 + 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 460 "chetrf_aa_2stage.f"
			i__3 = *n - i2;
#line 460 "chetrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[i2 
				+ (i2 + 1) * a_dim1], lda);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 463 "chetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 463 "chetrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 464 "chetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 464 "chetrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 464 "chetrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 465 "chetrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 465 "chetrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 467 "chetrf_aa_2stage.f"
			if (j > 0) {
#line 468 "chetrf_aa_2stage.f"
			    i__3 = j * nb;
#line 468 "chetrf_aa_2stage.f"
			    cswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * 
				    a_dim1 + 1], &c__1);
#line 470 "chetrf_aa_2stage.f"
			}
#line 471 "chetrf_aa_2stage.f"
		    }
#line 472 "chetrf_aa_2stage.f"
		}
#line 473 "chetrf_aa_2stage.f"
	    }
#line 474 "chetrf_aa_2stage.f"
	}
#line 475 "chetrf_aa_2stage.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 481 "chetrf_aa_2stage.f"
	i__1 = nt - 1;
#line 481 "chetrf_aa_2stage.f"
	for (j = 0; j <= i__1; ++j) {

/*           Generate Jth column of W and H */

/* Computing MIN */
#line 485 "chetrf_aa_2stage.f"
	    i__2 = nb, i__3 = *n - j * nb;
#line 485 "chetrf_aa_2stage.f"
	    kb = min(i__2,i__3);
#line 486 "chetrf_aa_2stage.f"
	    i__2 = j - 1;
#line 486 "chetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 487 "chetrf_aa_2stage.f"
		if (i__ == 1) {
/*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
#line 489 "chetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 490 "chetrf_aa_2stage.f"
			jb = nb + kb;
#line 491 "chetrf_aa_2stage.f"
		    } else {
#line 492 "chetrf_aa_2stage.f"
			jb = nb << 1;
#line 493 "chetrf_aa_2stage.f"
		    }
#line 494 "chetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 494 "chetrf_aa_2stage.f"
		    cgemm_("NoTranspose", "Conjugate transpose", &nb, &kb, &
			    jb, &c_b2, &tb[td + 1 + i__ * nb * ldtb], &i__3, &
			    a[j * nb + 1 + ((i__ - 1) * nb + 1) * a_dim1], 
			    lda, &c_b1, &work[i__ * nb + 1], n, (ftnlen)11, (
			    ftnlen)19);
#line 499 "chetrf_aa_2stage.f"
		} else {
/*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
#line 501 "chetrf_aa_2stage.f"
		    if (i__ == j - 1) {
#line 502 "chetrf_aa_2stage.f"
			jb = (nb << 1) + kb;
#line 503 "chetrf_aa_2stage.f"
		    } else {
#line 504 "chetrf_aa_2stage.f"
			jb = nb * 3;
#line 505 "chetrf_aa_2stage.f"
		    }
#line 506 "chetrf_aa_2stage.f"
		    i__3 = ldtb - 1;
#line 506 "chetrf_aa_2stage.f"
		    cgemm_("NoTranspose", "Conjugate transpose", &nb, &kb, &
			    jb, &c_b2, &tb[td + nb + 1 + (i__ - 1) * nb * 
			    ldtb], &i__3, &a[j * nb + 1 + ((i__ - 2) * nb + 1)
			     * a_dim1], lda, &c_b1, &work[i__ * nb + 1], n, (
			    ftnlen)11, (ftnlen)19);
#line 512 "chetrf_aa_2stage.f"
		}
#line 513 "chetrf_aa_2stage.f"
	    }

/*           Compute T(J,J) */

#line 517 "chetrf_aa_2stage.f"
	    i__2 = ldtb - 1;
#line 517 "chetrf_aa_2stage.f"
	    clacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1],
		     lda, &tb[td + 1 + j * nb * ldtb], &i__2, (ftnlen)5);
#line 519 "chetrf_aa_2stage.f"
	    if (j > 1) {
/*              T(J,J) = L(J,1:J)*H(1:J) */
#line 521 "chetrf_aa_2stage.f"
		i__2 = (j - 1) * nb;
#line 521 "chetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 521 "chetrf_aa_2stage.f"
		i__3 = ldtb - 1;
#line 521 "chetrf_aa_2stage.f"
		cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &z__1, &
			a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &c_b2, 
			&tb[td + 1 + j * nb * ldtb], &i__3, (ftnlen)11, (
			ftnlen)11);
/*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
#line 527 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 527 "chetrf_aa_2stage.f"
		cgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[
			j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[
			td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &
			work[1], n, (ftnlen)11, (ftnlen)11);
#line 532 "chetrf_aa_2stage.f"
		z__1.r = -1., z__1.i = -0.;
#line 532 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 532 "chetrf_aa_2stage.f"
		cgemm_("NoTranspose", "Conjugate transpose", &kb, &kb, &nb, &
			z__1, &work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) 
			* a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], &
			i__2, (ftnlen)11, (ftnlen)19);
#line 537 "chetrf_aa_2stage.f"
	    }
#line 538 "chetrf_aa_2stage.f"
	    if (j > 0) {
#line 539 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 539 "chetrf_aa_2stage.f"
		chegst_(&c__1, "Lower", &kb, &tb[td + 1 + j * nb * ldtb], &
			i__2, &a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], 
			lda, &iinfo, (ftnlen)5);
#line 542 "chetrf_aa_2stage.f"
	    }

/*           Expand T(J,J) into full format */

#line 546 "chetrf_aa_2stage.f"
	    i__2 = kb;
#line 546 "chetrf_aa_2stage.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 547 "chetrf_aa_2stage.f"
		i__3 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 547 "chetrf_aa_2stage.f"
		i__4 = td + 1 + (j * nb + i__ - 1) * ldtb;
#line 547 "chetrf_aa_2stage.f"
		d__1 = tb[i__4].r;
#line 547 "chetrf_aa_2stage.f"
		tb[i__3].r = d__1, tb[i__3].i = 0.;
#line 549 "chetrf_aa_2stage.f"
		i__3 = kb;
#line 549 "chetrf_aa_2stage.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 550 "chetrf_aa_2stage.f"
		    i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
#line 550 "chetrf_aa_2stage.f"
		    d_cnjg(&z__1, &tb[td + (k - i__) + 1 + (j * nb + i__ - 1) 
			    * ldtb]);
#line 550 "chetrf_aa_2stage.f"
		    tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 552 "chetrf_aa_2stage.f"
		}
#line 553 "chetrf_aa_2stage.f"
	    }

#line 555 "chetrf_aa_2stage.f"
	    if (j < nt - 1) {
#line 556 "chetrf_aa_2stage.f"
		if (j > 0) {

/*                 Compute H(J,J) */

#line 560 "chetrf_aa_2stage.f"
		    if (j == 1) {
#line 561 "chetrf_aa_2stage.f"
			i__2 = ldtb - 1;
#line 561 "chetrf_aa_2stage.f"
			cgemm_("NoTranspose", "Conjugate transpose", &kb, &kb,
				 &kb, &c_b2, &tb[td + 1 + j * nb * ldtb], &
				i__2, &a[j * nb + 1 + ((j - 1) * nb + 1) * 
				a_dim1], lda, &c_b1, &work[j * nb + 1], n, (
				ftnlen)11, (ftnlen)19);
#line 566 "chetrf_aa_2stage.f"
		    } else {
#line 567 "chetrf_aa_2stage.f"
			i__2 = nb + kb;
#line 567 "chetrf_aa_2stage.f"
			i__3 = ldtb - 1;
#line 567 "chetrf_aa_2stage.f"
			cgemm_("NoTranspose", "Conjugate transpose", &kb, &kb,
				 &i__2, &c_b2, &tb[td + nb + 1 + (j - 1) * nb 
				* ldtb], &i__3, &a[j * nb + 1 + ((j - 2) * nb 
				+ 1) * a_dim1], lda, &c_b1, &work[j * nb + 1],
				 n, (ftnlen)11, (ftnlen)19);
#line 573 "chetrf_aa_2stage.f"
		    }

/*                 Update with the previous column */

#line 577 "chetrf_aa_2stage.f"
		    i__2 = *n - (j + 1) * nb;
#line 577 "chetrf_aa_2stage.f"
		    i__3 = j * nb;
#line 577 "chetrf_aa_2stage.f"
		    z__1.r = -1., z__1.i = -0.;
#line 577 "chetrf_aa_2stage.f"
		    cgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, &
			    z__1, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[
			    nb + 1], n, &c_b2, &a[(j + 1) * nb + 1 + (j * nb 
			    + 1) * a_dim1], lda, (ftnlen)11, (ftnlen)11);
#line 582 "chetrf_aa_2stage.f"
		}

/*              Factorize panel */

#line 586 "chetrf_aa_2stage.f"
		i__2 = *n - (j + 1) * nb;
#line 586 "chetrf_aa_2stage.f"
		cgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * 
			a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
/*               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
/*                  INFO = IINFO+(J+1)*NB */
/*               END IF */

/*              Compute T(J+1, J), zero out for GEMM update */

/* Computing MIN */
#line 595 "chetrf_aa_2stage.f"
		i__2 = nb, i__3 = *n - (j + 1) * nb;
#line 595 "chetrf_aa_2stage.f"
		kb = min(i__2,i__3);
#line 596 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 596 "chetrf_aa_2stage.f"
		claset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * 
			nb * ldtb], &i__2, (ftnlen)4);
#line 598 "chetrf_aa_2stage.f"
		i__2 = ldtb - 1;
#line 598 "chetrf_aa_2stage.f"
		clacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) 
			* a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &
			i__2, (ftnlen)5);
#line 601 "chetrf_aa_2stage.f"
		if (j > 0) {
#line 602 "chetrf_aa_2stage.f"
		    i__2 = ldtb - 1;
#line 602 "chetrf_aa_2stage.f"
		    ctrsm_("R", "L", "C", "U", &kb, &nb, &c_b2, &a[j * nb + 1 
			    + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + nb 
			    + 1 + j * nb * ldtb], &i__2, (ftnlen)1, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 605 "chetrf_aa_2stage.f"
		}

/*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
/*              updates */

#line 610 "chetrf_aa_2stage.f"
		i__2 = nb;
#line 610 "chetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
#line 611 "chetrf_aa_2stage.f"
		    i__3 = kb;
#line 611 "chetrf_aa_2stage.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 612 "chetrf_aa_2stage.f"
			i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1)
				 * ldtb;
#line 612 "chetrf_aa_2stage.f"
			d_cnjg(&z__1, &tb[td + nb + i__ - k + 1 + (j * nb + k 
				- 1) * ldtb]);
#line 612 "chetrf_aa_2stage.f"
			tb[i__4].r = z__1.r, tb[i__4].i = z__1.i;
#line 614 "chetrf_aa_2stage.f"
		    }
#line 615 "chetrf_aa_2stage.f"
		}
#line 616 "chetrf_aa_2stage.f"
		claset_("Upper", &kb, &nb, &c_b1, &c_b2, &a[(j + 1) * nb + 1 
			+ (j * nb + 1) * a_dim1], lda, (ftnlen)5);

/*              Apply pivots to trailing submatrix of A */

#line 621 "chetrf_aa_2stage.f"
		i__2 = kb;
#line 621 "chetrf_aa_2stage.f"
		for (k = 1; k <= i__2; ++k) {
/*                 > Adjust ipiv */
#line 623 "chetrf_aa_2stage.f"
		    ipiv[(j + 1) * nb + k] += (j + 1) * nb;

#line 625 "chetrf_aa_2stage.f"
		    i1 = (j + 1) * nb + k;
#line 626 "chetrf_aa_2stage.f"
		    i2 = ipiv[(j + 1) * nb + k];
#line 627 "chetrf_aa_2stage.f"
		    if (i1 != i2) {
/*                    > Apply pivots to previous columns of L */
#line 629 "chetrf_aa_2stage.f"
			i__3 = k - 1;
#line 629 "chetrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], 
				lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], 
				lda);
/*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
#line 632 "chetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 632 "chetrf_aa_2stage.f"
			cswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ (i1 + 1) * a_dim1], lda);
#line 634 "chetrf_aa_2stage.f"
			i__3 = i2 - i1;
#line 634 "chetrf_aa_2stage.f"
			clacgv_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1);
#line 635 "chetrf_aa_2stage.f"
			i__3 = i2 - i1 - 1;
#line 635 "chetrf_aa_2stage.f"
			clacgv_(&i__3, &a[i2 + (i1 + 1) * a_dim1], lda);
/*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
#line 637 "chetrf_aa_2stage.f"
			i__3 = *n - i2;
#line 637 "chetrf_aa_2stage.f"
			cswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[i2 
				+ 1 + i2 * a_dim1], &c__1);
/*                    > Swap A(I1, I1) with A(I2, I2) */
#line 640 "chetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 640 "chetrf_aa_2stage.f"
			piv.r = a[i__3].r, piv.i = a[i__3].i;
#line 641 "chetrf_aa_2stage.f"
			i__3 = i1 + i1 * a_dim1;
#line 641 "chetrf_aa_2stage.f"
			i__4 = i2 + i2 * a_dim1;
#line 641 "chetrf_aa_2stage.f"
			a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 642 "chetrf_aa_2stage.f"
			i__3 = i2 + i2 * a_dim1;
#line 642 "chetrf_aa_2stage.f"
			a[i__3].r = piv.r, a[i__3].i = piv.i;
/*                    > Apply pivots to previous columns of L */
#line 644 "chetrf_aa_2stage.f"
			if (j > 0) {
#line 645 "chetrf_aa_2stage.f"
			    i__3 = j * nb;
#line 645 "chetrf_aa_2stage.f"
			    cswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + 
				    a_dim1], lda);
#line 647 "chetrf_aa_2stage.f"
			}
#line 648 "chetrf_aa_2stage.f"
		    }
#line 649 "chetrf_aa_2stage.f"
		}

/*              Apply pivots to previous columns of L */

/*               CALL CLASWP( J*NB, A( 1, 1 ), LDA, */
/*     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
#line 655 "chetrf_aa_2stage.f"
	    }
#line 656 "chetrf_aa_2stage.f"
	}
#line 657 "chetrf_aa_2stage.f"
    }

/*     Factor the band matrix */
#line 660 "chetrf_aa_2stage.f"
    cgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);

/*     End of CHETRF_AA_2STAGE */

#line 664 "chetrf_aa_2stage.f"
    return 0;
} /* chetrf_aa_2stage__ */

