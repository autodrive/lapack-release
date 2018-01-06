#line 1 "clauu2.f"
/* clauu2.f -- translated by f2c (version 20100827).
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

#line 1 "clauu2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLAUU2 computes the product UUH or LHL, where U and L are upper or lower triangular matrices (u
nblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAUU2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauu2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauu2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauu2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAUU2( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAUU2 computes the product U * U**H or L**H * L, where the triangular */
/* > factor U or L is stored in the upper or lower triangular part of */
/* > the array A. */
/* > */
/* > If UPLO = 'U' or 'u' then the upper triangle of the result is stored, */
/* > overwriting the factor U in A. */
/* > If UPLO = 'L' or 'l' then the lower triangle of the result is stored, */
/* > overwriting the factor L in A. */
/* > */
/* > This is the unblocked form of the algorithm, calling Level 2 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the triangular factor stored in the array A */
/* >          is upper or lower triangular: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the triangular factor U or L.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the triangular factor U or L. */
/* >          On exit, if UPLO = 'U', the upper triangle of A is */
/* >          overwritten with the upper triangle of the product U * U**H; */
/* >          if UPLO = 'L', the lower triangle of A is overwritten with */
/* >          the lower triangle of the product L**H * L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clauu2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Local variables */
    static integer i__;
    static doublereal aii;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , csscal_(integer *, doublereal *, doublecomplex *, integer *), 
	    xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
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

#line 144 "clauu2.f"
    /* Parameter adjustments */
#line 144 "clauu2.f"
    a_dim1 = *lda;
#line 144 "clauu2.f"
    a_offset = 1 + a_dim1;
#line 144 "clauu2.f"
    a -= a_offset;
#line 144 "clauu2.f"

#line 144 "clauu2.f"
    /* Function Body */
#line 144 "clauu2.f"
    *info = 0;
#line 145 "clauu2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 146 "clauu2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 147 "clauu2.f"
	*info = -1;
#line 148 "clauu2.f"
    } else if (*n < 0) {
#line 149 "clauu2.f"
	*info = -2;
#line 150 "clauu2.f"
    } else if (*lda < max(1,*n)) {
#line 151 "clauu2.f"
	*info = -4;
#line 152 "clauu2.f"
    }
#line 153 "clauu2.f"
    if (*info != 0) {
#line 154 "clauu2.f"
	i__1 = -(*info);
#line 154 "clauu2.f"
	xerbla_("CLAUU2", &i__1, (ftnlen)6);
#line 155 "clauu2.f"
	return 0;
#line 156 "clauu2.f"
    }

/*     Quick return if possible */

#line 160 "clauu2.f"
    if (*n == 0) {
#line 160 "clauu2.f"
	return 0;
#line 160 "clauu2.f"
    }

#line 163 "clauu2.f"
    if (upper) {

/*        Compute the product U * U**H. */

#line 167 "clauu2.f"
	i__1 = *n;
#line 167 "clauu2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 168 "clauu2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 168 "clauu2.f"
	    aii = a[i__2].r;
#line 169 "clauu2.f"
	    if (i__ < *n) {
#line 170 "clauu2.f"
		i__2 = i__ + i__ * a_dim1;
#line 170 "clauu2.f"
		i__3 = *n - i__;
#line 170 "clauu2.f"
		cdotc_(&z__1, &i__3, &a[i__ + (i__ + 1) * a_dim1], lda, &a[
			i__ + (i__ + 1) * a_dim1], lda);
#line 170 "clauu2.f"
		d__1 = aii * aii + z__1.r;
#line 170 "clauu2.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 172 "clauu2.f"
		i__2 = *n - i__;
#line 172 "clauu2.f"
		clacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 173 "clauu2.f"
		i__2 = i__ - 1;
#line 173 "clauu2.f"
		i__3 = *n - i__;
#line 173 "clauu2.f"
		z__1.r = aii, z__1.i = 0.;
#line 173 "clauu2.f"
		cgemv_("No transpose", &i__2, &i__3, &c_b1, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &a[i__ + (i__ + 1) * a_dim1], lda, &
			z__1, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 176 "clauu2.f"
		i__2 = *n - i__;
#line 176 "clauu2.f"
		clacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 177 "clauu2.f"
	    } else {
#line 178 "clauu2.f"
		csscal_(&i__, &aii, &a[i__ * a_dim1 + 1], &c__1);
#line 179 "clauu2.f"
	    }
#line 180 "clauu2.f"
/* L10: */
#line 180 "clauu2.f"
	}

#line 182 "clauu2.f"
    } else {

/*        Compute the product L**H * L. */

#line 186 "clauu2.f"
	i__1 = *n;
#line 186 "clauu2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 187 "clauu2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 187 "clauu2.f"
	    aii = a[i__2].r;
#line 188 "clauu2.f"
	    if (i__ < *n) {
#line 189 "clauu2.f"
		i__2 = i__ + i__ * a_dim1;
#line 189 "clauu2.f"
		i__3 = *n - i__;
#line 189 "clauu2.f"
		cdotc_(&z__1, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[
			i__ + 1 + i__ * a_dim1], &c__1);
#line 189 "clauu2.f"
		d__1 = aii * aii + z__1.r;
#line 189 "clauu2.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 191 "clauu2.f"
		i__2 = i__ - 1;
#line 191 "clauu2.f"
		clacgv_(&i__2, &a[i__ + a_dim1], lda);
#line 192 "clauu2.f"
		i__2 = *n - i__;
#line 192 "clauu2.f"
		i__3 = i__ - 1;
#line 192 "clauu2.f"
		z__1.r = aii, z__1.i = 0.;
#line 192 "clauu2.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b1, &a[i__ + 1 
			+ a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &
			z__1, &a[i__ + a_dim1], lda, (ftnlen)19);
#line 195 "clauu2.f"
		i__2 = i__ - 1;
#line 195 "clauu2.f"
		clacgv_(&i__2, &a[i__ + a_dim1], lda);
#line 196 "clauu2.f"
	    } else {
#line 197 "clauu2.f"
		csscal_(&i__, &aii, &a[i__ + a_dim1], lda);
#line 198 "clauu2.f"
	    }
#line 199 "clauu2.f"
/* L20: */
#line 199 "clauu2.f"
	}
#line 200 "clauu2.f"
    }

#line 202 "clauu2.f"
    return 0;

/*     End of CLAUU2 */

} /* clauu2_ */

