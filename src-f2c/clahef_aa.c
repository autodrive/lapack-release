#line 1 "clahef_aa.f"
/* clahef_aa.f -- translated by f2c (version 20100827).
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

#line 1 "clahef_aa.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLAHEF_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAHEF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAHEF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
/*                             H, LDH, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER    UPLO */
/*       INTEGER      J1, M, NB, LDA, LDH */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER      IPIV( * ) */
/*       COMPLEX      A( LDA, * ), H( LDH, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAHEF_AA factorizes a panel of a complex hermitian matrix A using */
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
/* >          when called by CHETRF_AA, for the first panel, J1 is 1, */
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
/* >          A is COMPLEX array, dimension (LDA,M) for */
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
/* >          H is COMPLEX workspace, dimension (LDH,NB). */
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
/* >          WORK is COMPLEX workspace, dimension (M). */
/* > \endverbatim */
/* > */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int clahef_aa__(char *uplo, integer *j1, integer *m, integer 
	*nb, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *
	h__, integer *ldh, doublecomplex *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, h_dim1, h_offset, i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, i1, k1, i2, mj;
    static doublecomplex piv, alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), clacgv_(
	    integer *, doublecomplex *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
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

#line 186 "clahef_aa.f"
    /* Parameter adjustments */
#line 186 "clahef_aa.f"
    a_dim1 = *lda;
#line 186 "clahef_aa.f"
    a_offset = 1 + a_dim1;
#line 186 "clahef_aa.f"
    a -= a_offset;
#line 186 "clahef_aa.f"
    --ipiv;
#line 186 "clahef_aa.f"
    h_dim1 = *ldh;
#line 186 "clahef_aa.f"
    h_offset = 1 + h_dim1;
#line 186 "clahef_aa.f"
    h__ -= h_offset;
#line 186 "clahef_aa.f"
    --work;
#line 186 "clahef_aa.f"

#line 186 "clahef_aa.f"
    /* Function Body */
#line 186 "clahef_aa.f"
    j = 1;

/*     K1 is the first column of the panel to be factorized */
/*     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks */

#line 191 "clahef_aa.f"
    k1 = 2 - *j1 + 1;

#line 193 "clahef_aa.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        ..................................................... */
/*        Factorize A as U**T*D*U using the upper triangle of A */
/*        ..................................................... */

#line 199 "clahef_aa.f"
L10:
#line 200 "clahef_aa.f"
	if (j > min(*m,*nb)) {
#line 200 "clahef_aa.f"
	    goto L20;
#line 200 "clahef_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from CHETRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 208 "clahef_aa.f"
	k = *j1 + j - 1;
#line 209 "clahef_aa.f"
	if (j == *m) {

/*            Only need to compute T(J, J) */

#line 213 "clahef_aa.f"
	    mj = 1;
#line 214 "clahef_aa.f"
	} else {
#line 215 "clahef_aa.f"
	    mj = *m - j + 1;
#line 216 "clahef_aa.f"
	}

/*        H(J:N, J) := A(J, J:N) - H(J:N, 1:(J-1)) * L(J1:(J-1), J), */
/*         where H(J:N, J) has been initialized to be A(J, J:N) */

#line 221 "clahef_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 229 "clahef_aa.f"
	    i__1 = j - k1;
#line 229 "clahef_aa.f"
	    clacgv_(&i__1, &a[j * a_dim1 + 1], &c__1);
#line 230 "clahef_aa.f"
	    i__1 = j - k1;
#line 230 "clahef_aa.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "clahef_aa.f"
	    cgemv_("No transpose", &mj, &i__1, &z__1, &h__[j + k1 * h_dim1], 
		    ldh, &a[j * a_dim1 + 1], &c__1, &c_b2, &h__[j + j * 
		    h_dim1], &c__1, (ftnlen)12);
#line 234 "clahef_aa.f"
	    i__1 = j - k1;
#line 234 "clahef_aa.f"
	    clacgv_(&i__1, &a[j * a_dim1 + 1], &c__1);
#line 235 "clahef_aa.f"
	}

/*        Copy H(i:n, i) into WORK */

#line 239 "clahef_aa.f"
	ccopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 241 "clahef_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J-1, J:N) * T(J-1,J), */
/*            where A(J-1, J) stores T(J-1, J) and A(J-2, J:N) stores U(J-1, J:N) */

#line 246 "clahef_aa.f"
	    d_cnjg(&z__2, &a[k - 1 + j * a_dim1]);
#line 246 "clahef_aa.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 246 "clahef_aa.f"
	    alpha.r = z__1.r, alpha.i = z__1.i;
#line 247 "clahef_aa.f"
	    caxpy_(&mj, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &c__1);
#line 248 "clahef_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 252 "clahef_aa.f"
	i__1 = k + j * a_dim1;
#line 252 "clahef_aa.f"
	d__1 = work[1].r;
#line 252 "clahef_aa.f"
	a[i__1].r = d__1, a[i__1].i = 0.;

#line 254 "clahef_aa.f"
	if (j < *m) {

/*           Compute WORK(2:N) = T(J, J) L(J, (J+1):N) */
/*            where A(J, J) stores T(J, J) and A(J-1, (J+1):N) stores U(J, (J+1):N) */

#line 259 "clahef_aa.f"
	    if (k > 1) {
#line 260 "clahef_aa.f"
		i__1 = k + j * a_dim1;
#line 260 "clahef_aa.f"
		z__1.r = -a[i__1].r, z__1.i = -a[i__1].i;
#line 260 "clahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 261 "clahef_aa.f"
		i__1 = *m - j;
#line 261 "clahef_aa.f"
		caxpy_(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, &
			work[2], &c__1);
#line 263 "clahef_aa.f"
	    }

/*           Find max(|WORK(2:n)|) */

#line 267 "clahef_aa.f"
	    i__1 = *m - j;
#line 267 "clahef_aa.f"
	    i2 = icamax_(&i__1, &work[2], &c__1) + 1;
#line 268 "clahef_aa.f"
	    i__1 = i2;
#line 268 "clahef_aa.f"
	    piv.r = work[i__1].r, piv.i = work[i__1].i;

/*           Apply hermitian pivot */

#line 272 "clahef_aa.f"
	    if (i2 != 2 && (piv.r != 0. || piv.i != 0.)) {

/*              Swap WORK(I1) and WORK(I2) */

#line 276 "clahef_aa.f"
		i1 = 2;
#line 277 "clahef_aa.f"
		i__1 = i2;
#line 277 "clahef_aa.f"
		i__2 = i1;
#line 277 "clahef_aa.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 278 "clahef_aa.f"
		i__1 = i1;
#line 278 "clahef_aa.f"
		work[i__1].r = piv.r, work[i__1].i = piv.i;

/*              Swap A(I1, I1+1:N) with A(I1+1:N, I2) */

#line 282 "clahef_aa.f"
		i1 = i1 + j - 1;
#line 283 "clahef_aa.f"
		i2 = i2 + j - 1;
#line 284 "clahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 284 "clahef_aa.f"
		cswap_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda, &a[*
			j1 + i1 + i2 * a_dim1], &c__1);
#line 286 "clahef_aa.f"
		i__1 = i2 - i1;
#line 286 "clahef_aa.f"
		clacgv_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda);
#line 287 "clahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 287 "clahef_aa.f"
		clacgv_(&i__1, &a[*j1 + i1 + i2 * a_dim1], &c__1);

/*              Swap A(I1, I2+1:N) with A(I2, I2+1:N) */

#line 291 "clahef_aa.f"
		i__1 = *m - i2;
#line 291 "clahef_aa.f"
		cswap_(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda, &a[*
			j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);

/*              Swap A(I1, I1) with A(I2,I2) */

#line 296 "clahef_aa.f"
		i__1 = i1 + *j1 - 1 + i1 * a_dim1;
#line 296 "clahef_aa.f"
		piv.r = a[i__1].r, piv.i = a[i__1].i;
#line 297 "clahef_aa.f"
		i__1 = *j1 + i1 - 1 + i1 * a_dim1;
#line 297 "clahef_aa.f"
		i__2 = *j1 + i2 - 1 + i2 * a_dim1;
#line 297 "clahef_aa.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 298 "clahef_aa.f"
		i__1 = *j1 + i2 - 1 + i2 * a_dim1;
#line 298 "clahef_aa.f"
		a[i__1].r = piv.r, a[i__1].i = piv.i;

/*              Swap H(I1, 1:J1) with H(I2, 1:J1) */

#line 302 "clahef_aa.f"
		i__1 = i1 - 1;
#line 302 "clahef_aa.f"
		cswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 303 "clahef_aa.f"
		ipiv[i1] = i2;

#line 305 "clahef_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 310 "clahef_aa.f"
		    i__1 = i1 - k1 + 1;
#line 310 "clahef_aa.f"
		    cswap_(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 
			    + 1], &c__1);
#line 312 "clahef_aa.f"
		}
#line 313 "clahef_aa.f"
	    } else {
#line 314 "clahef_aa.f"
		ipiv[j + 1] = j + 1;
#line 315 "clahef_aa.f"
	    }

/*           Set A(J, J+1) = T(J, J+1) */

#line 319 "clahef_aa.f"
	    i__1 = k + (j + 1) * a_dim1;
#line 319 "clahef_aa.f"
	    a[i__1].r = work[2].r, a[i__1].i = work[2].i;

#line 321 "clahef_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:N, J+1) into H(J:N, J), */

#line 325 "clahef_aa.f"
		i__1 = *m - j;
#line 325 "clahef_aa.f"
		ccopy_(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda, &h__[j + 1 + 
			(j + 1) * h_dim1], &c__1);
#line 327 "clahef_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1) */

#line 332 "clahef_aa.f"
	    i__1 = k + (j + 1) * a_dim1;
#line 332 "clahef_aa.f"
	    if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 333 "clahef_aa.f"
		z_div(&z__1, &c_b2, &a[k + (j + 1) * a_dim1]);
#line 333 "clahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 334 "clahef_aa.f"
		i__1 = *m - j - 1;
#line 334 "clahef_aa.f"
		ccopy_(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
#line 335 "clahef_aa.f"
		i__1 = *m - j - 1;
#line 335 "clahef_aa.f"
		cscal_(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
#line 336 "clahef_aa.f"
	    } else {
#line 337 "clahef_aa.f"
		i__1 = *m - j - 1;
#line 337 "clahef_aa.f"
		claset_("Full", &c__1, &i__1, &c_b1, &c_b1, &a[k + (j + 2) * 
			a_dim1], lda, (ftnlen)4);
#line 339 "clahef_aa.f"
	    }
#line 340 "clahef_aa.f"
	}
#line 341 "clahef_aa.f"
	++j;
#line 342 "clahef_aa.f"
	goto L10;
#line 343 "clahef_aa.f"
L20:

#line 345 "clahef_aa.f"
	;
#line 345 "clahef_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

#line 351 "clahef_aa.f"
L30:
#line 352 "clahef_aa.f"
	if (j > min(*m,*nb)) {
#line 352 "clahef_aa.f"
	    goto L40;
#line 352 "clahef_aa.f"
	}

/*        K is the column to be factorized */
/*         when being called from CHETRF_AA, */
/*         > for the first block column, J1 is 1, hence J1+J-1 is J, */
/*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */

#line 360 "clahef_aa.f"
	k = *j1 + j - 1;
#line 361 "clahef_aa.f"
	if (j == *m) {

/*            Only need to compute T(J, J) */

#line 365 "clahef_aa.f"
	    mj = 1;
#line 366 "clahef_aa.f"
	} else {
#line 367 "clahef_aa.f"
	    mj = *m - j + 1;
#line 368 "clahef_aa.f"
	}

/*        H(J:N, J) := A(J:N, J) - H(J:N, 1:(J-1)) * L(J, J1:(J-1))^T, */
/*         where H(J:N, J) has been initialized to be A(J:N, J) */

#line 373 "clahef_aa.f"
	if (k > 2) {

/*        K is the column to be factorized */
/*         > for the first block column, K is J, skipping the first two */
/*           columns */
/*         > for the rest of the columns, K is J+1, skipping only the */
/*           first column */

#line 381 "clahef_aa.f"
	    i__1 = j - k1;
#line 381 "clahef_aa.f"
	    clacgv_(&i__1, &a[j + a_dim1], lda);
#line 382 "clahef_aa.f"
	    i__1 = j - k1;
#line 382 "clahef_aa.f"
	    z__1.r = -1., z__1.i = -0.;
#line 382 "clahef_aa.f"
	    cgemv_("No transpose", &mj, &i__1, &z__1, &h__[j + k1 * h_dim1], 
		    ldh, &a[j + a_dim1], lda, &c_b2, &h__[j + j * h_dim1], &
		    c__1, (ftnlen)12);
#line 386 "clahef_aa.f"
	    i__1 = j - k1;
#line 386 "clahef_aa.f"
	    clacgv_(&i__1, &a[j + a_dim1], lda);
#line 387 "clahef_aa.f"
	}

/*        Copy H(J:N, J) into WORK */

#line 391 "clahef_aa.f"
	ccopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);

#line 393 "clahef_aa.f"
	if (j > k1) {

/*           Compute WORK := WORK - L(J:N, J-1) * T(J-1,J), */
/*            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */

#line 398 "clahef_aa.f"
	    d_cnjg(&z__2, &a[j + (k - 1) * a_dim1]);
#line 398 "clahef_aa.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 398 "clahef_aa.f"
	    alpha.r = z__1.r, alpha.i = z__1.i;
#line 399 "clahef_aa.f"
	    caxpy_(&mj, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 400 "clahef_aa.f"
	}

/*        Set A(J, J) = T(J, J) */

#line 404 "clahef_aa.f"
	i__1 = j + k * a_dim1;
#line 404 "clahef_aa.f"
	d__1 = work[1].r;
#line 404 "clahef_aa.f"
	a[i__1].r = d__1, a[i__1].i = 0.;

#line 406 "clahef_aa.f"
	if (j < *m) {

/*           Compute WORK(2:N) = T(J, J) L((J+1):N, J) */
/*            where A(J, J) = T(J, J) and A((J+1):N, J-1) = L((J+1):N, J) */

#line 411 "clahef_aa.f"
	    if (k > 1) {
#line 412 "clahef_aa.f"
		i__1 = j + k * a_dim1;
#line 412 "clahef_aa.f"
		z__1.r = -a[i__1].r, z__1.i = -a[i__1].i;
#line 412 "clahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 413 "clahef_aa.f"
		i__1 = *m - j;
#line 413 "clahef_aa.f"
		caxpy_(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, &
			work[2], &c__1);
#line 415 "clahef_aa.f"
	    }

/*           Find max(|WORK(2:n)|) */

#line 419 "clahef_aa.f"
	    i__1 = *m - j;
#line 419 "clahef_aa.f"
	    i2 = icamax_(&i__1, &work[2], &c__1) + 1;
#line 420 "clahef_aa.f"
	    i__1 = i2;
#line 420 "clahef_aa.f"
	    piv.r = work[i__1].r, piv.i = work[i__1].i;

/*           Apply hermitian pivot */

#line 424 "clahef_aa.f"
	    if (i2 != 2 && (piv.r != 0. || piv.i != 0.)) {

/*              Swap WORK(I1) and WORK(I2) */

#line 428 "clahef_aa.f"
		i1 = 2;
#line 429 "clahef_aa.f"
		i__1 = i2;
#line 429 "clahef_aa.f"
		i__2 = i1;
#line 429 "clahef_aa.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 430 "clahef_aa.f"
		i__1 = i1;
#line 430 "clahef_aa.f"
		work[i__1].r = piv.r, work[i__1].i = piv.i;

/*              Swap A(I1+1:N, I1) with A(I2, I1+1:N) */

#line 434 "clahef_aa.f"
		i1 = i1 + j - 1;
#line 435 "clahef_aa.f"
		i2 = i2 + j - 1;
#line 436 "clahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 436 "clahef_aa.f"
		cswap_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + (*j1 + i1) * a_dim1], lda);
#line 438 "clahef_aa.f"
		i__1 = i2 - i1;
#line 438 "clahef_aa.f"
		clacgv_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1);
#line 439 "clahef_aa.f"
		i__1 = i2 - i1 - 1;
#line 439 "clahef_aa.f"
		clacgv_(&i__1, &a[i2 + (*j1 + i1) * a_dim1], lda);

/*              Swap A(I2+1:N, I1) with A(I2+1:N, I2) */

#line 443 "clahef_aa.f"
		i__1 = *m - i2;
#line 443 "clahef_aa.f"
		cswap_(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[
			i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);

/*              Swap A(I1, I1) with A(I2, I2) */

#line 448 "clahef_aa.f"
		i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
#line 448 "clahef_aa.f"
		piv.r = a[i__1].r, piv.i = a[i__1].i;
#line 449 "clahef_aa.f"
		i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
#line 449 "clahef_aa.f"
		i__2 = i2 + (*j1 + i2 - 1) * a_dim1;
#line 449 "clahef_aa.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 450 "clahef_aa.f"
		i__1 = i2 + (*j1 + i2 - 1) * a_dim1;
#line 450 "clahef_aa.f"
		a[i__1].r = piv.r, a[i__1].i = piv.i;

/*              Swap H(I1, I1:J1) with H(I2, I2:J1) */

#line 454 "clahef_aa.f"
		i__1 = i1 - 1;
#line 454 "clahef_aa.f"
		cswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
#line 455 "clahef_aa.f"
		ipiv[i1] = i2;

#line 457 "clahef_aa.f"
		if (i1 > k1 - 1) {

/*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
/*                  skipping the first column */

#line 462 "clahef_aa.f"
		    i__1 = i1 - k1 + 1;
#line 462 "clahef_aa.f"
		    cswap_(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
#line 464 "clahef_aa.f"
		}
#line 465 "clahef_aa.f"
	    } else {
#line 466 "clahef_aa.f"
		ipiv[j + 1] = j + 1;
#line 467 "clahef_aa.f"
	    }

/*           Set A(J+1, J) = T(J+1, J) */

#line 471 "clahef_aa.f"
	    i__1 = j + 1 + k * a_dim1;
#line 471 "clahef_aa.f"
	    a[i__1].r = work[2].r, a[i__1].i = work[2].i;

#line 473 "clahef_aa.f"
	    if (j < *nb) {

/*              Copy A(J+1:N, J+1) into H(J+1:N, J), */

#line 477 "clahef_aa.f"
		i__1 = *m - j;
#line 477 "clahef_aa.f"
		ccopy_(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1, &h__[j + 1 
			+ (j + 1) * h_dim1], &c__1);
#line 479 "clahef_aa.f"
	    }

/*           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1), */
/*            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1) */

#line 484 "clahef_aa.f"
	    i__1 = j + 1 + k * a_dim1;
#line 484 "clahef_aa.f"
	    if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 485 "clahef_aa.f"
		z_div(&z__1, &c_b2, &a[j + 1 + k * a_dim1]);
#line 485 "clahef_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 486 "clahef_aa.f"
		i__1 = *m - j - 1;
#line 486 "clahef_aa.f"
		ccopy_(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], &c__1);
#line 487 "clahef_aa.f"
		i__1 = *m - j - 1;
#line 487 "clahef_aa.f"
		cscal_(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
#line 488 "clahef_aa.f"
	    } else {
#line 489 "clahef_aa.f"
		i__1 = *m - j - 1;
#line 489 "clahef_aa.f"
		claset_("Full", &i__1, &c__1, &c_b1, &c_b1, &a[j + 2 + k * 
			a_dim1], lda, (ftnlen)4);
#line 491 "clahef_aa.f"
	    }
#line 492 "clahef_aa.f"
	}
#line 493 "clahef_aa.f"
	++j;
#line 494 "clahef_aa.f"
	goto L30;
#line 495 "clahef_aa.f"
L40:
#line 496 "clahef_aa.f"
	;
#line 496 "clahef_aa.f"
    }
#line 497 "clahef_aa.f"
    return 0;

/*     End of CLAHEF_AA */

} /* clahef_aa__ */

