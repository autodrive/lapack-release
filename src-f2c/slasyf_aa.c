#line 1 "slasyf_aa.f"
/* slasyf_aa.f -- translated by f2c (version 20100827).
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

#line 1 "slasyf_aa.f"
/* Table of constant values */

static doublereal c_b6 = -1.;
static integer c__1 = 1;
static doublereal c_b8 = 1.;
static doublereal c_b22 = 0.;

/* > \brief \b SLASYF_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASYF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASYF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
/*                             H, LDH, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            J1, M, NB, LDA, LDH */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), H( LDH, * ), WORK( * ) */
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
/* >          when called by SSYTRF_AA, for the first panel, J1 is 1, */
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
/* >          A is REAL array, dimension (LDA,M) for */
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
/* >          H is REAL workspace, dimension (LDH,NB). */
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
/* >          WORK is REAL workspace, dimension (M). */
/* > \endverbatim */
/* > */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int slasyf_aa__(char *uplo, integer *j1, integer *m, integer 
	*nb, doublereal *a, integer *lda, integer *ipiv, doublereal *h__, 
	integer *ldh, doublereal *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, h_dim1, h_offset, i__1;

    /* Local variables */
    static integer j, k, i1, k1, i2, mj;
    static doublereal piv, alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), scopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), sswap_(integer *, doublereal 
	    *, integer *, doublereal *, integer *), saxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
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

#line 186 "slasyf_aa.f"
    /* Parameter adjustments */
#line 186 "slasyf_aa.f"
    a_dim1 = *lda;
#line 186 "slasyf_aa.f"
    a_offset = 1 + a_dim1;
#line 186 "slasyf_aa.f"
    a -= a_offset;
#line 186 "slasyf_aa.f"
    --ipiv;
#line 186 "slasyf_aa.f"
    h_dim1 = *ldh;
#line 186 "slasyf_aa.f"
    h_offset = 1 + h_dim1;
#line 186 "slasyf_aa.f"
    h__ -= h_offset;
#line 186 "slasyf_aa.f"
    --work;
#line 186 "slasyf_aa.f"

#line 186 "slasyf_aa.f"
    /* Function Body */
#line 186 "slasyf_aa.f"
    j = 1;

/*     K1 is the first column of the panel to be factorized */
/*     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks */

#line 191 "slasyf_aa.f"
    k1 = 2 - *j1 + 1;

#line 193 "slasyf_aa.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        ..................................................... */
/*        Factorize A as U**T*D*U using the upper triangle of A */
/*        ..................................................... */

#line 199 "slasyf_aa.f"
L10:
#line 200 "slasyf_aa.f"
	if (j > min(*m,*nb)) {
#line 200 "slasyf_aa.f"
	    goto L20;
#line 200 "slasyf_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from SSYTRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 208 "slasyf_aa.f"
	k = *j1 + j - 1;
#line 209 "slasyf_aa.f"
	if (j == *m) {

/*            Only need to compute T(J, J) */

#line 213 "slasyf_aa.f"
	    mj = 1;
#line 214 "slasyf_aa.f"
	} else {
#line 215 "slasyf_aa.f"
	    mj = *m - j + 1;
#line 216 "slasyf_aa.f"
	}

/*        H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J), */
/*         where H(J:M, J) has been initialized to be A(J, J:M) */

#line 221 "slasyf_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 229 "slasyf_aa.f"
	    i__1 = j - k1;
#line 229 "slasyf_aa.f"
	    sgemv_("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], 
		    ldh, &a[j * a_dim1 + 1], &c__1, &c_b8, &h__[j + j * 
		    h_dim1], &c__1, (ftnlen)12);
#line 233 "slasyf_aa.f"
	}

/*        Copy H(i:M, i) into WORK */

#line 237 "slasyf_aa.f"
	scopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 239 "slasyf_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J-1, J:M) * T(J-1,J), */
/*            where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M) */

#line 244 "slasyf_aa.f"
	    alpha = -a[k - 1 + j * a_dim1];
#line 245 "slasyf_aa.f"
	    saxpy_(&mj, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &c__1);
#line 246 "slasyf_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 250 "slasyf_aa.f"
	a[k + j * a_dim1] = work[1];

#line 252 "slasyf_aa.f"
	if (j < *m) {

/*           Compute WORK(2:M) = T(J, J) L(J, (J+1):M) */
/*            where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M) */

#line 257 "slasyf_aa.f"
	    if (k > 1) {
#line 258 "slasyf_aa.f"
		alpha = -a[k + j * a_dim1];
#line 259 "slasyf_aa.f"
		i__1 = *m - j;
#line 259 "slasyf_aa.f"
		saxpy_(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, &
			work[2], &c__1);
#line 261 "slasyf_aa.f"
	    }

/*           Find max(|WORK(2:M)|) */

#line 265 "slasyf_aa.f"
	    i__1 = *m - j;
#line 265 "slasyf_aa.f"
	    i2 = isamax_(&i__1, &work[2], &c__1) + 1;
#line 266 "slasyf_aa.f"
	    piv = work[i2];

/*           Apply symmetric pivot */

#line 270 "slasyf_aa.f"
	    if (i2 != 2 && piv != 0.) {

/*              Swap WORK(I1) and WORK(I2) */

#line 274 "slasyf_aa.f"
		i1 = 2;
#line 275 "slasyf_aa.f"
		work[i2] = work[i1];
#line 276 "slasyf_aa.f"
		work[i1] = piv;

/*              Swap A(I1, I1+1:M) with A(I1+1:M, I2) */

#line 280 "slasyf_aa.f"
		i1 = i1 + j - 1;
#line 281 "slasyf_aa.f"
		i2 = i2 + j - 1;
#line 282 "slasyf_aa.f"
		i__1 = i2 - i1 - 1;
#line 282 "slasyf_aa.f"
		sswap_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda, &a[*
			j1 + i1 + i2 * a_dim1], &c__1);

/*              Swap A(I1, I2+1:M) with A(I2, I2+1:M) */

#line 287 "slasyf_aa.f"
		i__1 = *m - i2;
#line 287 "slasyf_aa.f"
		sswap_(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda, &a[*
			j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);

/*              Swap A(I1, I1) with A(I2,I2) */

#line 292 "slasyf_aa.f"
		piv = a[i1 + *j1 - 1 + i1 * a_dim1];
#line 293 "slasyf_aa.f"
		a[*j1 + i1 - 1 + i1 * a_dim1] = a[*j1 + i2 - 1 + i2 * a_dim1];
#line 294 "slasyf_aa.f"
		a[*j1 + i2 - 1 + i2 * a_dim1] = piv;

/*              Swap H(I1, 1:J1) with H(I2, 1:J1) */

#line 298 "slasyf_aa.f"
		i__1 = i1 - 1;
#line 298 "slasyf_aa.f"
		sswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 299 "slasyf_aa.f"
		ipiv[i1] = i2;

#line 301 "slasyf_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 306 "slasyf_aa.f"
		    i__1 = i1 - k1 + 1;
#line 306 "slasyf_aa.f"
		    sswap_(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 
			    + 1], &c__1);
#line 308 "slasyf_aa.f"
		}
#line 309 "slasyf_aa.f"
	    } else {
#line 310 "slasyf_aa.f"
		ipiv[j + 1] = j + 1;
#line 311 "slasyf_aa.f"
	    }

/*           Set A(J, J+1) = T(J, J+1) */

#line 315 "slasyf_aa.f"
	    a[k + (j + 1) * a_dim1] = work[2];

#line 317 "slasyf_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:M, J+1) into H(J:M, J), */

#line 321 "slasyf_aa.f"
		i__1 = *m - j;
#line 321 "slasyf_aa.f"
		scopy_(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda, &h__[j + 1 + 
			(j + 1) * h_dim1], &c__1);
#line 323 "slasyf_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */

#line 328 "slasyf_aa.f"
	    if (a[k + (j + 1) * a_dim1] != 0.) {
#line 329 "slasyf_aa.f"
		alpha = 1. / a[k + (j + 1) * a_dim1];
#line 330 "slasyf_aa.f"
		i__1 = *m - j - 1;
#line 330 "slasyf_aa.f"
		scopy_(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
#line 331 "slasyf_aa.f"
		i__1 = *m - j - 1;
#line 331 "slasyf_aa.f"
		sscal_(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
#line 332 "slasyf_aa.f"
	    } else {
#line 333 "slasyf_aa.f"
		i__1 = *m - j - 1;
#line 333 "slasyf_aa.f"
		slaset_("Full", &c__1, &i__1, &c_b22, &c_b22, &a[k + (j + 2) *
			 a_dim1], lda, (ftnlen)4);
#line 335 "slasyf_aa.f"
	    }
#line 336 "slasyf_aa.f"
	}
#line 337 "slasyf_aa.f"
	++j;
#line 338 "slasyf_aa.f"
	goto L10;
#line 339 "slasyf_aa.f"
L20:

#line 341 "slasyf_aa.f"
	;
#line 341 "slasyf_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 347 "slasyf_aa.f"
L30:
#line 348 "slasyf_aa.f"
	if (j > min(*m,*nb)) {
#line 348 "slasyf_aa.f"
	    goto L40;
#line 348 "slasyf_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from SSYTRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 356 "slasyf_aa.f"
	k = *j1 + j - 1;
#line 357 "slasyf_aa.f"
	if (j == *m) {

/*            Only need to compute T(J, J) */

#line 361 "slasyf_aa.f"
	    mj = 1;
#line 362 "slasyf_aa.f"
	} else {
#line 363 "slasyf_aa.f"
	    mj = *m - j + 1;
#line 364 "slasyf_aa.f"
	}

/*        H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T, */
/*         where H(J:M, J) has been initialized to be A(J:M, J) */

#line 369 "slasyf_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 377 "slasyf_aa.f"
	    i__1 = j - k1;
#line 377 "slasyf_aa.f"
	    sgemv_("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], 
		    ldh, &a[j + a_dim1], lda, &c_b8, &h__[j + j * h_dim1], &
		    c__1, (ftnlen)12);
#line 381 "slasyf_aa.f"
	}

/*        Copy H(J:M, J) into WORK */

#line 385 "slasyf_aa.f"
	scopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 387 "slasyf_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J:M, J-1) * T(J-1,J), */
/*            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */

#line 392 "slasyf_aa.f"
	    alpha = -a[j + (k - 1) * a_dim1];
#line 393 "slasyf_aa.f"
	    saxpy_(&mj, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 394 "slasyf_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 398 "slasyf_aa.f"
	a[j + k * a_dim1] = work[1];

#line 400 "slasyf_aa.f"
	if (j < *m) {

/*           Compute WORK(2:M) = T(J, J) L((J+1):M, J) */
/*            where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J) */

#line 405 "slasyf_aa.f"
	    if (k > 1) {
#line 406 "slasyf_aa.f"
		alpha = -a[j + k * a_dim1];
#line 407 "slasyf_aa.f"
		i__1 = *m - j;
#line 407 "slasyf_aa.f"
		saxpy_(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, &
			work[2], &c__1);
#line 409 "slasyf_aa.f"
	    }

/*           Find max(|WORK(2:M)|) */

#line 413 "slasyf_aa.f"
	    i__1 = *m - j;
#line 413 "slasyf_aa.f"
	    i2 = isamax_(&i__1, &work[2], &c__1) + 1;
#line 414 "slasyf_aa.f"
	    piv = work[i2];

/*           Apply symmetric pivot */

#line 418 "slasyf_aa.f"
	    if (i2 != 2 && piv != 0.) {

/*              Swap WORK(I1) and WORK(I2) */

#line 422 "slasyf_aa.f"
		i1 = 2;
#line 423 "slasyf_aa.f"
		work[i2] = work[i1];
#line 424 "slasyf_aa.f"
		work[i1] = piv;

/*              Swap A(I1+1:M, I1) with A(I2, I1+1:M) */

#line 428 "slasyf_aa.f"
		i1 = i1 + j - 1;
#line 429 "slasyf_aa.f"
		i2 = i2 + j - 1;
#line 430 "slasyf_aa.f"
		i__1 = i2 - i1 - 1;
#line 430 "slasyf_aa.f"
		sswap_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + (*j1 + i1) * a_dim1], lda);

/*              Swap A(I2+1:M, I1) with A(I2+1:M, I2) */

#line 435 "slasyf_aa.f"
		i__1 = *m - i2;
#line 435 "slasyf_aa.f"
		sswap_(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);

/*              Swap A(I1, I1) with A(I2, I2) */

#line 440 "slasyf_aa.f"
		piv = a[i1 + (*j1 + i1 - 1) * a_dim1];
#line 441 "slasyf_aa.f"
		a[i1 + (*j1 + i1 - 1) * a_dim1] = a[i2 + (*j1 + i2 - 1) * 
			a_dim1];
#line 442 "slasyf_aa.f"
		a[i2 + (*j1 + i2 - 1) * a_dim1] = piv;

/*              Swap H(I1, I1:J1) with H(I2, I2:J1) */

#line 446 "slasyf_aa.f"
		i__1 = i1 - 1;
#line 446 "slasyf_aa.f"
		sswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 447 "slasyf_aa.f"
		ipiv[i1] = i2;

#line 449 "slasyf_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 454 "slasyf_aa.f"
		    i__1 = i1 - k1 + 1;
#line 454 "slasyf_aa.f"
		    sswap_(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
#line 456 "slasyf_aa.f"
		}
#line 457 "slasyf_aa.f"
	    } else {
#line 458 "slasyf_aa.f"
		ipiv[j + 1] = j + 1;
#line 459 "slasyf_aa.f"
	    }

/*           Set A(J+1, J) = T(J+1, J) */

#line 463 "slasyf_aa.f"
	    a[j + 1 + k * a_dim1] = work[2];

#line 465 "slasyf_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:M, J+1) into H(J+1:M, J), */

#line 469 "slasyf_aa.f"
		i__1 = *m - j;
#line 469 "slasyf_aa.f"
		scopy_(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1, &h__[j + 1 
			+ (j + 1) * h_dim1], &c__1);
#line 471 "slasyf_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */

#line 476 "slasyf_aa.f"
	    if (a[j + 1 + k * a_dim1] != 0.) {
#line 477 "slasyf_aa.f"
		alpha = 1. / a[j + 1 + k * a_dim1];
#line 478 "slasyf_aa.f"
		i__1 = *m - j - 1;
#line 478 "slasyf_aa.f"
		scopy_(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], &c__1);
#line 479 "slasyf_aa.f"
		i__1 = *m - j - 1;
#line 479 "slasyf_aa.f"
		sscal_(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
#line 480 "slasyf_aa.f"
	    } else {
#line 481 "slasyf_aa.f"
		i__1 = *m - j - 1;
#line 481 "slasyf_aa.f"
		slaset_("Full", &i__1, &c__1, &c_b22, &c_b22, &a[j + 2 + k * 
			a_dim1], lda, (ftnlen)4);
#line 483 "slasyf_aa.f"
	    }
#line 484 "slasyf_aa.f"
	}
#line 485 "slasyf_aa.f"
	++j;
#line 486 "slasyf_aa.f"
	goto L30;
#line 487 "slasyf_aa.f"
L40:
#line 488 "slasyf_aa.f"
	;
#line 488 "slasyf_aa.f"
    }
#line 489 "slasyf_aa.f"
    return 0;

/*     End of SLASYF_AA */

} /* slasyf_aa__ */

