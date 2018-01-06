#line 1 "dlasyf_aa.f"
/* dlasyf_aa.f -- translated by f2c (version 20100827).
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

#line 1 "dlasyf_aa.f"
/* Table of constant values */

static doublereal c_b6 = -1.;
static integer c__1 = 1;
static doublereal c_b8 = 1.;
static doublereal c_b22 = 0.;

/* > \brief \b DLASYF_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASYF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasyf_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasyf_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasyf_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASYF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
/*                             H, LDH, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            J1, M, NB, LDA, LDH */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), H( LDH, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATRF_AA factorizes a panel of a real symmetric matrix A using */
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
/* >          when called by DSYTRF_AA, for the first panel, J1 is 1, */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,M) for */
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
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (M) */
/* >          Details of the row and column interchanges, */
/* >          the row and column k were interchanged with the row and */
/* >          column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION workspace, dimension (LDH,NB). */
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
/* >          WORK is DOUBLE PRECISION workspace, dimension (M). */
/* > \endverbatim */
/* > */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasyf_aa__(char *uplo, integer *j1, integer *m, integer 
	*nb, doublereal *a, integer *lda, integer *ipiv, doublereal *h__, 
	integer *ldh, doublereal *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, h_dim1, h_offset, i__1;

    /* Local variables */
    static integer j, k, i1, k1, i2, mj;
    static doublereal piv, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);


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

#line 186 "dlasyf_aa.f"
    /* Parameter adjustments */
#line 186 "dlasyf_aa.f"
    a_dim1 = *lda;
#line 186 "dlasyf_aa.f"
    a_offset = 1 + a_dim1;
#line 186 "dlasyf_aa.f"
    a -= a_offset;
#line 186 "dlasyf_aa.f"
    --ipiv;
#line 186 "dlasyf_aa.f"
    h_dim1 = *ldh;
#line 186 "dlasyf_aa.f"
    h_offset = 1 + h_dim1;
#line 186 "dlasyf_aa.f"
    h__ -= h_offset;
#line 186 "dlasyf_aa.f"
    --work;
#line 186 "dlasyf_aa.f"

#line 186 "dlasyf_aa.f"
    /* Function Body */
#line 186 "dlasyf_aa.f"
    j = 1;

/*     K1 is the first column of the panel to be factorized */
/*     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks */

#line 191 "dlasyf_aa.f"
    k1 = 2 - *j1 + 1;

#line 193 "dlasyf_aa.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        ..................................................... */
/*        Factorize A as U**T*D*U using the upper triangle of A */
/*        ..................................................... */

#line 199 "dlasyf_aa.f"
L10:
#line 200 "dlasyf_aa.f"
	if (j > min(*m,*nb)) {
#line 200 "dlasyf_aa.f"
	    goto L20;
#line 200 "dlasyf_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from DSYTRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 208 "dlasyf_aa.f"
	k = *j1 + j - 1;
#line 209 "dlasyf_aa.f"
	if (j == *m) {

/*            Only need to compute T(J, J) */

#line 213 "dlasyf_aa.f"
	    mj = 1;
#line 214 "dlasyf_aa.f"
	} else {
#line 215 "dlasyf_aa.f"
	    mj = *m - j + 1;
#line 216 "dlasyf_aa.f"
	}

/*        H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J), */
/*         where H(J:M, J) has been initialized to be A(J, J:M) */

#line 221 "dlasyf_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 229 "dlasyf_aa.f"
	    i__1 = j - k1;
#line 229 "dlasyf_aa.f"
	    dgemv_("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], 
		    ldh, &a[j * a_dim1 + 1], &c__1, &c_b8, &h__[j + j * 
		    h_dim1], &c__1, (ftnlen)12);
#line 233 "dlasyf_aa.f"
	}

/*        Copy H(i:M, i) into WORK */

#line 237 "dlasyf_aa.f"
	dcopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 239 "dlasyf_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J-1, J:M) * T(J-1,J), */
/*            where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M) */

#line 244 "dlasyf_aa.f"
	    alpha = -a[k - 1 + j * a_dim1];
#line 245 "dlasyf_aa.f"
	    daxpy_(&mj, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &c__1);
#line 246 "dlasyf_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 250 "dlasyf_aa.f"
	a[k + j * a_dim1] = work[1];

#line 252 "dlasyf_aa.f"
	if (j < *m) {

/*           Compute WORK(2:M) = T(J, J) L(J, (J+1):M) */
/*            where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M) */

#line 257 "dlasyf_aa.f"
	    if (k > 1) {
#line 258 "dlasyf_aa.f"
		alpha = -a[k + j * a_dim1];
#line 259 "dlasyf_aa.f"
		i__1 = *m - j;
#line 259 "dlasyf_aa.f"
		daxpy_(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, &
			work[2], &c__1);
#line 261 "dlasyf_aa.f"
	    }

/*           Find max(|WORK(2:M)|) */

#line 265 "dlasyf_aa.f"
	    i__1 = *m - j;
#line 265 "dlasyf_aa.f"
	    i2 = idamax_(&i__1, &work[2], &c__1) + 1;
#line 266 "dlasyf_aa.f"
	    piv = work[i2];

/*           Apply symmetric pivot */

#line 270 "dlasyf_aa.f"
	    if (i2 != 2 && piv != 0.) {

/*              Swap WORK(I1) and WORK(I2) */

#line 274 "dlasyf_aa.f"
		i1 = 2;
#line 275 "dlasyf_aa.f"
		work[i2] = work[i1];
#line 276 "dlasyf_aa.f"
		work[i1] = piv;

/*              Swap A(I1, I1+1:M) with A(I1+1:M, I2) */

#line 280 "dlasyf_aa.f"
		i1 = i1 + j - 1;
#line 281 "dlasyf_aa.f"
		i2 = i2 + j - 1;
#line 282 "dlasyf_aa.f"
		i__1 = i2 - i1 - 1;
#line 282 "dlasyf_aa.f"
		dswap_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda, &a[*
			j1 + i1 + i2 * a_dim1], &c__1);

/*              Swap A(I1, I2+1:M) with A(I2, I2+1:M) */

#line 287 "dlasyf_aa.f"
		i__1 = *m - i2;
#line 287 "dlasyf_aa.f"
		dswap_(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda, &a[*
			j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);

/*              Swap A(I1, I1) with A(I2,I2) */

#line 292 "dlasyf_aa.f"
		piv = a[i1 + *j1 - 1 + i1 * a_dim1];
#line 293 "dlasyf_aa.f"
		a[*j1 + i1 - 1 + i1 * a_dim1] = a[*j1 + i2 - 1 + i2 * a_dim1];
#line 294 "dlasyf_aa.f"
		a[*j1 + i2 - 1 + i2 * a_dim1] = piv;

/*              Swap H(I1, 1:J1) with H(I2, 1:J1) */

#line 298 "dlasyf_aa.f"
		i__1 = i1 - 1;
#line 298 "dlasyf_aa.f"
		dswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 299 "dlasyf_aa.f"
		ipiv[i1] = i2;

#line 301 "dlasyf_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 306 "dlasyf_aa.f"
		    i__1 = i1 - k1 + 1;
#line 306 "dlasyf_aa.f"
		    dswap_(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 
			    + 1], &c__1);
#line 308 "dlasyf_aa.f"
		}
#line 309 "dlasyf_aa.f"
	    } else {
#line 310 "dlasyf_aa.f"
		ipiv[j + 1] = j + 1;
#line 311 "dlasyf_aa.f"
	    }

/*           Set A(J, J+1) = T(J, J+1) */

#line 315 "dlasyf_aa.f"
	    a[k + (j + 1) * a_dim1] = work[2];

#line 317 "dlasyf_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:M, J+1) into H(J:M, J), */

#line 321 "dlasyf_aa.f"
		i__1 = *m - j;
#line 321 "dlasyf_aa.f"
		dcopy_(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda, &h__[j + 1 + 
			(j + 1) * h_dim1], &c__1);
#line 323 "dlasyf_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */

#line 328 "dlasyf_aa.f"
	    if (a[k + (j + 1) * a_dim1] != 0.) {
#line 329 "dlasyf_aa.f"
		alpha = 1. / a[k + (j + 1) * a_dim1];
#line 330 "dlasyf_aa.f"
		i__1 = *m - j - 1;
#line 330 "dlasyf_aa.f"
		dcopy_(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
#line 331 "dlasyf_aa.f"
		i__1 = *m - j - 1;
#line 331 "dlasyf_aa.f"
		dscal_(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
#line 332 "dlasyf_aa.f"
	    } else {
#line 333 "dlasyf_aa.f"
		i__1 = *m - j - 1;
#line 333 "dlasyf_aa.f"
		dlaset_("Full", &c__1, &i__1, &c_b22, &c_b22, &a[k + (j + 2) *
			 a_dim1], lda, (ftnlen)4);
#line 335 "dlasyf_aa.f"
	    }
#line 336 "dlasyf_aa.f"
	}
#line 337 "dlasyf_aa.f"
	++j;
#line 338 "dlasyf_aa.f"
	goto L10;
#line 339 "dlasyf_aa.f"
L20:

#line 341 "dlasyf_aa.f"
	;
#line 341 "dlasyf_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 347 "dlasyf_aa.f"
L30:
#line 348 "dlasyf_aa.f"
	if (j > min(*m,*nb)) {
#line 348 "dlasyf_aa.f"
	    goto L40;
#line 348 "dlasyf_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from DSYTRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 356 "dlasyf_aa.f"
	k = *j1 + j - 1;
#line 357 "dlasyf_aa.f"
	if (j == *m) {

/*            Only need to compute T(J, J) */

#line 361 "dlasyf_aa.f"
	    mj = 1;
#line 362 "dlasyf_aa.f"
	} else {
#line 363 "dlasyf_aa.f"
	    mj = *m - j + 1;
#line 364 "dlasyf_aa.f"
	}

/*        H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T, */
/*         where H(J:M, J) has been initialized to be A(J:M, J) */

#line 369 "dlasyf_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 377 "dlasyf_aa.f"
	    i__1 = j - k1;
#line 377 "dlasyf_aa.f"
	    dgemv_("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], 
		    ldh, &a[j + a_dim1], lda, &c_b8, &h__[j + j * h_dim1], &
		    c__1, (ftnlen)12);
#line 381 "dlasyf_aa.f"
	}

/*        Copy H(J:M, J) into WORK */

#line 385 "dlasyf_aa.f"
	dcopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 387 "dlasyf_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J:M, J-1) * T(J-1,J), */
/*            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */

#line 392 "dlasyf_aa.f"
	    alpha = -a[j + (k - 1) * a_dim1];
#line 393 "dlasyf_aa.f"
	    daxpy_(&mj, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 394 "dlasyf_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 398 "dlasyf_aa.f"
	a[j + k * a_dim1] = work[1];

#line 400 "dlasyf_aa.f"
	if (j < *m) {

/*           Compute WORK(2:M) = T(J, J) L((J+1):M, J) */
/*            where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J) */

#line 405 "dlasyf_aa.f"
	    if (k > 1) {
#line 406 "dlasyf_aa.f"
		alpha = -a[j + k * a_dim1];
#line 407 "dlasyf_aa.f"
		i__1 = *m - j;
#line 407 "dlasyf_aa.f"
		daxpy_(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, &
			work[2], &c__1);
#line 409 "dlasyf_aa.f"
	    }

/*           Find max(|WORK(2:M)|) */

#line 413 "dlasyf_aa.f"
	    i__1 = *m - j;
#line 413 "dlasyf_aa.f"
	    i2 = idamax_(&i__1, &work[2], &c__1) + 1;
#line 414 "dlasyf_aa.f"
	    piv = work[i2];

/*           Apply symmetric pivot */

#line 418 "dlasyf_aa.f"
	    if (i2 != 2 && piv != 0.) {

/*              Swap WORK(I1) and WORK(I2) */

#line 422 "dlasyf_aa.f"
		i1 = 2;
#line 423 "dlasyf_aa.f"
		work[i2] = work[i1];
#line 424 "dlasyf_aa.f"
		work[i1] = piv;

/*              Swap A(I1+1:M, I1) with A(I2, I1+1:M) */

#line 428 "dlasyf_aa.f"
		i1 = i1 + j - 1;
#line 429 "dlasyf_aa.f"
		i2 = i2 + j - 1;
#line 430 "dlasyf_aa.f"
		i__1 = i2 - i1 - 1;
#line 430 "dlasyf_aa.f"
		dswap_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + (*j1 + i1) * a_dim1], lda);

/*              Swap A(I2+1:M, I1) with A(I2+1:M, I2) */

#line 435 "dlasyf_aa.f"
		i__1 = *m - i2;
#line 435 "dlasyf_aa.f"
		dswap_(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);

/*              Swap A(I1, I1) with A(I2, I2) */

#line 440 "dlasyf_aa.f"
		piv = a[i1 + (*j1 + i1 - 1) * a_dim1];
#line 441 "dlasyf_aa.f"
		a[i1 + (*j1 + i1 - 1) * a_dim1] = a[i2 + (*j1 + i2 - 1) * 
			a_dim1];
#line 442 "dlasyf_aa.f"
		a[i2 + (*j1 + i2 - 1) * a_dim1] = piv;

/*              Swap H(I1, I1:J1) with H(I2, I2:J1) */

#line 446 "dlasyf_aa.f"
		i__1 = i1 - 1;
#line 446 "dlasyf_aa.f"
		dswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 447 "dlasyf_aa.f"
		ipiv[i1] = i2;

#line 449 "dlasyf_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 454 "dlasyf_aa.f"
		    i__1 = i1 - k1 + 1;
#line 454 "dlasyf_aa.f"
		    dswap_(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
#line 456 "dlasyf_aa.f"
		}
#line 457 "dlasyf_aa.f"
	    } else {
#line 458 "dlasyf_aa.f"
		ipiv[j + 1] = j + 1;
#line 459 "dlasyf_aa.f"
	    }

/*           Set A(J+1, J) = T(J+1, J) */

#line 463 "dlasyf_aa.f"
	    a[j + 1 + k * a_dim1] = work[2];

#line 465 "dlasyf_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:M, J+1) into H(J+1:M, J), */

#line 469 "dlasyf_aa.f"
		i__1 = *m - j;
#line 469 "dlasyf_aa.f"
		dcopy_(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1, &h__[j + 1 
			+ (j + 1) * h_dim1], &c__1);
#line 471 "dlasyf_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */

#line 476 "dlasyf_aa.f"
	    if (a[j + 1 + k * a_dim1] != 0.) {
#line 477 "dlasyf_aa.f"
		alpha = 1. / a[j + 1 + k * a_dim1];
#line 478 "dlasyf_aa.f"
		i__1 = *m - j - 1;
#line 478 "dlasyf_aa.f"
		dcopy_(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], &c__1);
#line 479 "dlasyf_aa.f"
		i__1 = *m - j - 1;
#line 479 "dlasyf_aa.f"
		dscal_(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
#line 480 "dlasyf_aa.f"
	    } else {
#line 481 "dlasyf_aa.f"
		i__1 = *m - j - 1;
#line 481 "dlasyf_aa.f"
		dlaset_("Full", &i__1, &c__1, &c_b22, &c_b22, &a[j + 2 + k * 
			a_dim1], lda, (ftnlen)4);
#line 483 "dlasyf_aa.f"
	    }
#line 484 "dlasyf_aa.f"
	}
#line 485 "dlasyf_aa.f"
	++j;
#line 486 "dlasyf_aa.f"
	goto L30;
#line 487 "dlasyf_aa.f"
L40:
#line 488 "dlasyf_aa.f"
	;
#line 488 "dlasyf_aa.f"
    }
#line 489 "dlasyf_aa.f"
    return 0;

/*     End of DLASYF_AA */

} /* dlasyf_aa__ */

