#line 1 "zlahef_aa.f"
/* zlahef_aa.f -- translated by f2c (version 20100827).
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

#line 1 "zlahef_aa.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLAHEF_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAHEF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahef_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahef_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahef_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAHEF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
/*                             H, LDH, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER    UPLO */
/*       INTEGER      J1, M, NB, LDA, LDH, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER      IPIV( * ) */
/*       COMPLEX*16   A( LDA, * ), H( LDH, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAHEF_AA factorizes a panel of a complex hermitian matrix A using */
/* > the Aasen's algorithm. The panel consists of a set of NB rows of A */
/* > when UPLO is U, or a set of NB columns when UPLO is L. */
/* > */
/* > In order to factorize the panel, the Aasen's algorithm requires the */
/* > last row, or column, of the previous panel. The first row, or column, */
/* > of A is set to be the first row, or column, of an identity matrix, */
/* > which is used to factorize the first panel. */
/* > */
/* > The resulting J-th row of U, or J-th column of L, is stored in the */
/* > (J-1)-th row, or column, of A (without the unit diagonals), while */
/* > the diagonal and subdiagonal of A are overwritten by those of T. */
/* > */
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
/* > \param[in] J1 */
/* > \verbatim */
/* >          J1 is INTEGER */
/* >          The location of the first row, or column, of the panel */
/* >          within the submatrix of A, passed to this routine, e.g., */
/* >          when called by ZHETRF_AA, for the first panel, J1 is 1, */
/* >          while for the remaining panels, J1 is 2. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The dimension of the submatrix. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The dimension of the panel to be facotorized. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,M) for */
/* >          the first panel, while dimension (LDA,M+1) for the */
/* >          remaining panels. */
/* > */
/* >          On entry, A contains the last row, or column, of */
/* >          the previous panel, and the trailing submatrix of A */
/* >          to be factorized, except for the first panel, only */
/* >          the panel is passed. */
/* > */
/* >          On exit, the leading panel is factorized. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the row and column interchanges, */
/* >          the row and column k were interchanged with the row and */
/* >          column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX*16 workspace, dimension (LDH,NB). */
/* > */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the workspace H. LDH >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 workspace, dimension (M). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization */
/* >                has been completed, but the block diagonal matrix D is */
/* >                exactly singular, and division by zero will occur if it */
/* >                is used to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zlahef_aa__(char *uplo, integer *j1, integer *m, integer 
	*nb, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *
	h__, integer *ldh, doublecomplex *work, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, h_dim1, h_offset, i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, i1, k1, i2;
    static doublecomplex piv, alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zlacgv_(
	    integer *, doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);


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

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 195 "zlahef_aa.f"
    /* Parameter adjustments */
#line 195 "zlahef_aa.f"
    a_dim1 = *lda;
#line 195 "zlahef_aa.f"
    a_offset = 1 + a_dim1;
#line 195 "zlahef_aa.f"
    a -= a_offset;
#line 195 "zlahef_aa.f"
    --ipiv;
#line 195 "zlahef_aa.f"
    h_dim1 = *ldh;
#line 195 "zlahef_aa.f"
    h_offset = 1 + h_dim1;
#line 195 "zlahef_aa.f"
    h__ -= h_offset;
#line 195 "zlahef_aa.f"
    --work;
#line 195 "zlahef_aa.f"

#line 195 "zlahef_aa.f"
    /* Function Body */
#line 195 "zlahef_aa.f"
    *info = 0;
#line 196 "zlahef_aa.f"
    j = 1;

/*     K1 is the first column of the panel to be factorized */
/*     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks */

#line 201 "zlahef_aa.f"
    k1 = 2 - *j1 + 1;

#line 203 "zlahef_aa.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        ..................................................... */
/*        Factorize A as U**T*D*U using the upper triangle of A */
/*        ..................................................... */

#line 209 "zlahef_aa.f"
L10:
#line 210 "zlahef_aa.f"
	if (j > min(*m,*nb)) {
#line 210 "zlahef_aa.f"
	    goto L20;
#line 210 "zlahef_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from ZHETRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 218 "zlahef_aa.f"
	k = *j1 + j - 1;

/*        H(J:N, J) := A(J, J:N) - H(J:N, 1:(J-1)) * L(J1:(J-1), J), */
/*         where H(J:N, J) has been initialized to be A(J, J:N) */

#line 223 "zlahef_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 231 "zlahef_aa.f"
	    i__1 = j - k1;
#line 231 "zlahef_aa.f"
	    zlacgv_(&i__1, &a[j * a_dim1 + 1], &c__1);
#line 232 "zlahef_aa.f"
	    i__1 = *m - j + 1;
#line 232 "zlahef_aa.f"
	    i__2 = j - k1;
#line 232 "zlahef_aa.f"
	    z__1.r = -1., z__1.i = -0.;
#line 232 "zlahef_aa.f"
	    zgemv_("No transpose", &i__1, &i__2, &z__1, &h__[j + k1 * h_dim1],
		     ldh, &a[j * a_dim1 + 1], &c__1, &c_b2, &h__[j + j * 
		    h_dim1], &c__1, (ftnlen)12);
#line 236 "zlahef_aa.f"
	    i__1 = j - k1;
#line 236 "zlahef_aa.f"
	    zlacgv_(&i__1, &a[j * a_dim1 + 1], &c__1);
#line 237 "zlahef_aa.f"
	}

/*        Copy H(i:n, i) into WORK */

#line 241 "zlahef_aa.f"
	i__1 = *m - j + 1;
#line 241 "zlahef_aa.f"
	zcopy_(&i__1, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 243 "zlahef_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J-1, J:N) * T(J-1,J), */
/*            where A(J-1, J) stores T(J-1, J) and A(J-2, J:N) stores U(J-1, J:N) */

#line 248 "zlahef_aa.f"
	    d_cnjg(&z__2, &a[k - 1 + j * a_dim1]);
#line 248 "zlahef_aa.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 248 "zlahef_aa.f"
	    alpha.r = z__1.r, alpha.i = z__1.i;
#line 249 "zlahef_aa.f"
	    i__1 = *m - j + 1;
#line 249 "zlahef_aa.f"
	    zaxpy_(&i__1, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &
		    c__1);
#line 250 "zlahef_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 254 "zlahef_aa.f"
	i__1 = k + j * a_dim1;
#line 254 "zlahef_aa.f"
	d__1 = work[1].r;
#line 254 "zlahef_aa.f"
	a[i__1].r = d__1, a[i__1].i = 0.;

#line 256 "zlahef_aa.f"
	if (j < *m) {

/*           Compute WORK(2:N) = T(J, J) L(J, (J+1):N) */
/*            where A(J, J) stores T(J, J) and A(J-1, (J+1):N) stores U(J, (J+1):N) */

#line 261 "zlahef_aa.f"
	    if (k > 1) {
#line 262 "zlahef_aa.f"
		i__1 = k + j * a_dim1;
#line 262 "zlahef_aa.f"
		z__1.r = -a[i__1].r, z__1.i = -a[i__1].i;
#line 262 "zlahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 263 "zlahef_aa.f"
		i__1 = *m - j;
#line 263 "zlahef_aa.f"
		zaxpy_(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, &
			work[2], &c__1);
#line 265 "zlahef_aa.f"
	    }

/*           Find max(|WORK(2:n)|) */

#line 269 "zlahef_aa.f"
	    i__1 = *m - j;
#line 269 "zlahef_aa.f"
	    i2 = izamax_(&i__1, &work[2], &c__1) + 1;
#line 270 "zlahef_aa.f"
	    i__1 = i2;
#line 270 "zlahef_aa.f"
	    piv.r = work[i__1].r, piv.i = work[i__1].i;

/*           Apply hermitian pivot */

#line 274 "zlahef_aa.f"
	    if (i2 != 2 && (piv.r != 0. || piv.i != 0.)) {

/*              Swap WORK(I1) and WORK(I2) */

#line 278 "zlahef_aa.f"
		i1 = 2;
#line 279 "zlahef_aa.f"
		i__1 = i2;
#line 279 "zlahef_aa.f"
		i__2 = i1;
#line 279 "zlahef_aa.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 280 "zlahef_aa.f"
		i__1 = i1;
#line 280 "zlahef_aa.f"
		work[i__1].r = piv.r, work[i__1].i = piv.i;

/*              Swap A(I1, I1+1:N) with A(I1+1:N, I2) */

#line 284 "zlahef_aa.f"
		i1 = i1 + j - 1;
#line 285 "zlahef_aa.f"
		i2 = i2 + j - 1;
#line 286 "zlahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 286 "zlahef_aa.f"
		zswap_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda, &a[*
			j1 + i1 + i2 * a_dim1], &c__1);
#line 288 "zlahef_aa.f"
		i__1 = i2 - i1;
#line 288 "zlahef_aa.f"
		zlacgv_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda);
#line 289 "zlahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 289 "zlahef_aa.f"
		zlacgv_(&i__1, &a[*j1 + i1 + i2 * a_dim1], &c__1);

/*              Swap A(I1, I2+1:N) with A(I2, I2+1:N) */

#line 293 "zlahef_aa.f"
		i__1 = *m - i2;
#line 293 "zlahef_aa.f"
		zswap_(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda, &a[*
			j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);

/*              Swap A(I1, I1) with A(I2,I2) */

#line 298 "zlahef_aa.f"
		i__1 = i1 + *j1 - 1 + i1 * a_dim1;
#line 298 "zlahef_aa.f"
		piv.r = a[i__1].r, piv.i = a[i__1].i;
#line 299 "zlahef_aa.f"
		i__1 = *j1 + i1 - 1 + i1 * a_dim1;
#line 299 "zlahef_aa.f"
		i__2 = *j1 + i2 - 1 + i2 * a_dim1;
#line 299 "zlahef_aa.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 300 "zlahef_aa.f"
		i__1 = *j1 + i2 - 1 + i2 * a_dim1;
#line 300 "zlahef_aa.f"
		a[i__1].r = piv.r, a[i__1].i = piv.i;

/*              Swap H(I1, 1:J1) with H(I2, 1:J1) */

#line 304 "zlahef_aa.f"
		i__1 = i1 - 1;
#line 304 "zlahef_aa.f"
		zswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 305 "zlahef_aa.f"
		ipiv[i1] = i2;

#line 307 "zlahef_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 312 "zlahef_aa.f"
		    i__1 = i1 - k1 + 1;
#line 312 "zlahef_aa.f"
		    zswap_(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 
			    + 1], &c__1);
#line 314 "zlahef_aa.f"
		}
#line 315 "zlahef_aa.f"
	    } else {
#line 316 "zlahef_aa.f"
		ipiv[j + 1] = j + 1;
#line 317 "zlahef_aa.f"
	    }

/*           Set A(J, J+1) = T(J, J+1) */

#line 321 "zlahef_aa.f"
	    i__1 = k + (j + 1) * a_dim1;
#line 321 "zlahef_aa.f"
	    a[i__1].r = work[2].r, a[i__1].i = work[2].i;
#line 322 "zlahef_aa.f"
	    i__1 = k + j * a_dim1;
#line 322 "zlahef_aa.f"
	    i__2 = k + (j + 1) * a_dim1;
#line 322 "zlahef_aa.f"
	    if (a[i__1].r == 0. && a[i__1].i == 0. && (j == *m || a[i__2].r ==
		     0. && a[i__2].i == 0.)) {
#line 324 "zlahef_aa.f"
		if (*info == 0) {
#line 325 "zlahef_aa.f"
		    *info = j;
#line 326 "zlahef_aa.f"
		}
#line 327 "zlahef_aa.f"
	    }

#line 329 "zlahef_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:N, J+1) into H(J:N, J), */

#line 333 "zlahef_aa.f"
		i__1 = *m - j;
#line 333 "zlahef_aa.f"
		zcopy_(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda, &h__[j + 1 + 
			(j + 1) * h_dim1], &c__1);
#line 335 "zlahef_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1) */

#line 340 "zlahef_aa.f"
	    i__1 = k + (j + 1) * a_dim1;
#line 340 "zlahef_aa.f"
	    if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 341 "zlahef_aa.f"
		z_div(&z__1, &c_b2, &a[k + (j + 1) * a_dim1]);
#line 341 "zlahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 342 "zlahef_aa.f"
		i__1 = *m - j - 1;
#line 342 "zlahef_aa.f"
		zcopy_(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
#line 343 "zlahef_aa.f"
		i__1 = *m - j - 1;
#line 343 "zlahef_aa.f"
		zscal_(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
#line 344 "zlahef_aa.f"
	    } else {
#line 345 "zlahef_aa.f"
		i__1 = *m - j - 1;
#line 345 "zlahef_aa.f"
		zlaset_("Full", &c__1, &i__1, &c_b1, &c_b1, &a[k + (j + 2) * 
			a_dim1], lda, (ftnlen)4);
#line 347 "zlahef_aa.f"
	    }
#line 348 "zlahef_aa.f"
	} else {
#line 349 "zlahef_aa.f"
	    i__1 = k + j * a_dim1;
#line 349 "zlahef_aa.f"
	    if (a[i__1].r == 0. && a[i__1].i == 0. && *info == 0) {
#line 350 "zlahef_aa.f"
		*info = j;
#line 351 "zlahef_aa.f"
	    }
#line 352 "zlahef_aa.f"
	}
#line 353 "zlahef_aa.f"
	++j;
#line 354 "zlahef_aa.f"
	goto L10;
#line 355 "zlahef_aa.f"
L20:

#line 357 "zlahef_aa.f"
	;
#line 357 "zlahef_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 363 "zlahef_aa.f"
L30:
#line 364 "zlahef_aa.f"
	if (j > min(*m,*nb)) {
#line 364 "zlahef_aa.f"
	    goto L40;
#line 364 "zlahef_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from ZHETRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 372 "zlahef_aa.f"
	k = *j1 + j - 1;

/*        H(J:N, J) := A(J:N, J) - H(J:N, 1:(J-1)) * L(J, J1:(J-1))^T, */
/*         where H(J:N, J) has been initialized to be A(J:N, J) */

#line 377 "zlahef_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 385 "zlahef_aa.f"
	    i__1 = j - k1;
#line 385 "zlahef_aa.f"
	    zlacgv_(&i__1, &a[j + a_dim1], lda);
#line 386 "zlahef_aa.f"
	    i__1 = *m - j + 1;
#line 386 "zlahef_aa.f"
	    i__2 = j - k1;
#line 386 "zlahef_aa.f"
	    z__1.r = -1., z__1.i = -0.;
#line 386 "zlahef_aa.f"
	    zgemv_("No transpose", &i__1, &i__2, &z__1, &h__[j + k1 * h_dim1],
		     ldh, &a[j + a_dim1], lda, &c_b2, &h__[j + j * h_dim1], &
		    c__1, (ftnlen)12);
#line 390 "zlahef_aa.f"
	    i__1 = j - k1;
#line 390 "zlahef_aa.f"
	    zlacgv_(&i__1, &a[j + a_dim1], lda);
#line 391 "zlahef_aa.f"
	}

/*        Copy H(J:N, J) into WORK */

#line 395 "zlahef_aa.f"
	i__1 = *m - j + 1;
#line 395 "zlahef_aa.f"
	zcopy_(&i__1, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 397 "zlahef_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J:N, J-1) * T(J-1,J), */
/*            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */

#line 402 "zlahef_aa.f"
	    d_cnjg(&z__2, &a[j + (k - 1) * a_dim1]);
#line 402 "zlahef_aa.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 402 "zlahef_aa.f"
	    alpha.r = z__1.r, alpha.i = z__1.i;
#line 403 "zlahef_aa.f"
	    i__1 = *m - j + 1;
#line 403 "zlahef_aa.f"
	    zaxpy_(&i__1, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 404 "zlahef_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 408 "zlahef_aa.f"
	i__1 = j + k * a_dim1;
#line 408 "zlahef_aa.f"
	d__1 = work[1].r;
#line 408 "zlahef_aa.f"
	a[i__1].r = d__1, a[i__1].i = 0.;

#line 410 "zlahef_aa.f"
	if (j < *m) {

/*           Compute WORK(2:N) = T(J, J) L((J+1):N, J) */
/*            where A(J, J) = T(J, J) and A((J+1):N, J-1) = L((J+1):N, J) */

#line 415 "zlahef_aa.f"
	    if (k > 1) {
#line 416 "zlahef_aa.f"
		i__1 = j + k * a_dim1;
#line 416 "zlahef_aa.f"
		z__1.r = -a[i__1].r, z__1.i = -a[i__1].i;
#line 416 "zlahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 417 "zlahef_aa.f"
		i__1 = *m - j;
#line 417 "zlahef_aa.f"
		zaxpy_(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, &
			work[2], &c__1);
#line 419 "zlahef_aa.f"
	    }

/*           Find max(|WORK(2:n)|) */

#line 423 "zlahef_aa.f"
	    i__1 = *m - j;
#line 423 "zlahef_aa.f"
	    i2 = izamax_(&i__1, &work[2], &c__1) + 1;
#line 424 "zlahef_aa.f"
	    i__1 = i2;
#line 424 "zlahef_aa.f"
	    piv.r = work[i__1].r, piv.i = work[i__1].i;

/*           Apply hermitian pivot */

#line 428 "zlahef_aa.f"
	    if (i2 != 2 && (piv.r != 0. || piv.i != 0.)) {

/*              Swap WORK(I1) and WORK(I2) */

#line 432 "zlahef_aa.f"
		i1 = 2;
#line 433 "zlahef_aa.f"
		i__1 = i2;
#line 433 "zlahef_aa.f"
		i__2 = i1;
#line 433 "zlahef_aa.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 434 "zlahef_aa.f"
		i__1 = i1;
#line 434 "zlahef_aa.f"
		work[i__1].r = piv.r, work[i__1].i = piv.i;

/*              Swap A(I1+1:N, I1) with A(I2, I1+1:N) */

#line 438 "zlahef_aa.f"
		i1 = i1 + j - 1;
#line 439 "zlahef_aa.f"
		i2 = i2 + j - 1;
#line 440 "zlahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 440 "zlahef_aa.f"
		zswap_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + (*j1 + i1) * a_dim1], lda);
#line 442 "zlahef_aa.f"
		i__1 = i2 - i1;
#line 442 "zlahef_aa.f"
		zlacgv_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1);
#line 443 "zlahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 443 "zlahef_aa.f"
		zlacgv_(&i__1, &a[i2 + (*j1 + i1) * a_dim1], lda);

/*              Swap A(I2+1:N, I1) with A(I2+1:N, I2) */

#line 447 "zlahef_aa.f"
		i__1 = *m - i2;
#line 447 "zlahef_aa.f"
		zswap_(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);

/*              Swap A(I1, I1) with A(I2, I2) */

#line 452 "zlahef_aa.f"
		i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
#line 452 "zlahef_aa.f"
		piv.r = a[i__1].r, piv.i = a[i__1].i;
#line 453 "zlahef_aa.f"
		i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
#line 453 "zlahef_aa.f"
		i__2 = i2 + (*j1 + i2 - 1) * a_dim1;
#line 453 "zlahef_aa.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 454 "zlahef_aa.f"
		i__1 = i2 + (*j1 + i2 - 1) * a_dim1;
#line 454 "zlahef_aa.f"
		a[i__1].r = piv.r, a[i__1].i = piv.i;

/*              Swap H(I1, I1:J1) with H(I2, I2:J1) */

#line 458 "zlahef_aa.f"
		i__1 = i1 - 1;
#line 458 "zlahef_aa.f"
		zswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 459 "zlahef_aa.f"
		ipiv[i1] = i2;

#line 461 "zlahef_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 466 "zlahef_aa.f"
		    i__1 = i1 - k1 + 1;
#line 466 "zlahef_aa.f"
		    zswap_(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
#line 468 "zlahef_aa.f"
		}
#line 469 "zlahef_aa.f"
	    } else {
#line 470 "zlahef_aa.f"
		ipiv[j + 1] = j + 1;
#line 471 "zlahef_aa.f"
	    }

/*           Set A(J+1, J) = T(J+1, J) */

#line 475 "zlahef_aa.f"
	    i__1 = j + 1 + k * a_dim1;
#line 475 "zlahef_aa.f"
	    a[i__1].r = work[2].r, a[i__1].i = work[2].i;
#line 476 "zlahef_aa.f"
	    i__1 = j + k * a_dim1;
#line 476 "zlahef_aa.f"
	    i__2 = j + 1 + k * a_dim1;
#line 476 "zlahef_aa.f"
	    if (a[i__1].r == 0. && a[i__1].i == 0. && (j == *m || a[i__2].r ==
		     0. && a[i__2].i == 0.)) {
#line 478 "zlahef_aa.f"
		if (*info == 0) {
#line 478 "zlahef_aa.f"
		    *info = j;
#line 478 "zlahef_aa.f"
		}
#line 480 "zlahef_aa.f"
	    }

#line 482 "zlahef_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:N, J+1) into H(J+1:N, J), */

#line 486 "zlahef_aa.f"
		i__1 = *m - j;
#line 486 "zlahef_aa.f"
		zcopy_(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1, &h__[j + 1 
			+ (j + 1) * h_dim1], &c__1);
#line 488 "zlahef_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1) */

#line 493 "zlahef_aa.f"
	    i__1 = j + 1 + k * a_dim1;
#line 493 "zlahef_aa.f"
	    if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 494 "zlahef_aa.f"
		z_div(&z__1, &c_b2, &a[j + 1 + k * a_dim1]);
#line 494 "zlahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 495 "zlahef_aa.f"
		i__1 = *m - j - 1;
#line 495 "zlahef_aa.f"
		zcopy_(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], &c__1);
#line 496 "zlahef_aa.f"
		i__1 = *m - j - 1;
#line 496 "zlahef_aa.f"
		zscal_(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
#line 497 "zlahef_aa.f"
	    } else {
#line 498 "zlahef_aa.f"
		i__1 = *m - j - 1;
#line 498 "zlahef_aa.f"
		zlaset_("Full", &i__1, &c__1, &c_b1, &c_b1, &a[j + 2 + k * 
			a_dim1], lda, (ftnlen)4);
#line 500 "zlahef_aa.f"
	    }
#line 501 "zlahef_aa.f"
	} else {
#line 502 "zlahef_aa.f"
	    i__1 = j + k * a_dim1;
#line 502 "zlahef_aa.f"
	    if (a[i__1].r == 0. && a[i__1].i == 0. && j == *m && *info == 0) {
#line 502 "zlahef_aa.f"
		*info = j;
#line 502 "zlahef_aa.f"
	    }
#line 504 "zlahef_aa.f"
	}
#line 505 "zlahef_aa.f"
	++j;
#line 506 "zlahef_aa.f"
	goto L30;
#line 507 "zlahef_aa.f"
L40:
#line 508 "zlahef_aa.f"
	;
#line 508 "zlahef_aa.f"
    }
#line 509 "zlahef_aa.f"
    return 0;

/*     End of ZLAHEF_AA */

} /* zlahef_aa__ */

