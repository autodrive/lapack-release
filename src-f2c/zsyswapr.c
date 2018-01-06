#line 1 "zsyswapr.f"
/* zsyswapr.f -- translated by f2c (version 20100827).
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

#line 1 "zsyswapr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZSYSWAPR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYSWAPR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyswap
r.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyswap
r.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyswap
r.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYSWAPR( UPLO, N, A, LDA, I1, I2) */

/*       .. Scalar Arguments .. */
/*       CHARACTER        UPLO */
/*       INTEGER          I1, I2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16       A( LDA, N ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYSWAPR applies an elementary permutation on the rows and the columns of */
/* > a symmetric matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
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
/* >          On entry, the NB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by ZSYTRF. */
/* > */
/* >          On exit, if INFO = 0, the (symmetric) inverse of the original */
/* >          matrix.  If UPLO = 'U', the upper triangular part of the */
/* >          inverse is formed and the part of A below the diagonal is not */
/* >          referenced; if UPLO = 'L' the lower triangular part of the */
/* >          inverse is formed and the part of A above the diagonal is */
/* >          not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] I1 */
/* > \verbatim */
/* >          I1 is INTEGER */
/* >          Index of the first row to swap */
/* > \endverbatim */
/* > */
/* > \param[in] I2 */
/* > \verbatim */
/* >          I2 is INTEGER */
/* >          Index of the second row to swap */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16SYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zsyswapr_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *i1, integer *i2, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    static doublecomplex tmp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 134 "zsyswapr.f"
    /* Parameter adjustments */
#line 134 "zsyswapr.f"
    a_dim1 = *lda;
#line 134 "zsyswapr.f"
    a_offset = 1 + a_dim1;
#line 134 "zsyswapr.f"
    a -= a_offset;
#line 134 "zsyswapr.f"

#line 134 "zsyswapr.f"
    /* Function Body */
#line 134 "zsyswapr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 135 "zsyswapr.f"
    if (upper) {

/*         UPPER */
/*         first swap */
/*          - swap column I1 and I2 from I1 to I1-1 */
#line 140 "zsyswapr.f"
	i__1 = *i1 - 1;
#line 140 "zsyswapr.f"
	zswap_(&i__1, &a[*i1 * a_dim1 + 1], &c__1, &a[*i2 * a_dim1 + 1], &
		c__1);

/*          second swap : */
/*          - swap A(I1,I1) and A(I2,I2) */
/*          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1 */
#line 145 "zsyswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 145 "zsyswapr.f"
	tmp.r = a[i__1].r, tmp.i = a[i__1].i;
#line 146 "zsyswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 146 "zsyswapr.f"
	i__2 = *i2 + *i2 * a_dim1;
#line 146 "zsyswapr.f"
	a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 147 "zsyswapr.f"
	i__1 = *i2 + *i2 * a_dim1;
#line 147 "zsyswapr.f"
	a[i__1].r = tmp.r, a[i__1].i = tmp.i;

#line 149 "zsyswapr.f"
	i__1 = *i2 - *i1 - 1;
#line 149 "zsyswapr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 150 "zsyswapr.f"
	    i__2 = *i1 + (*i1 + i__) * a_dim1;
#line 150 "zsyswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 151 "zsyswapr.f"
	    i__2 = *i1 + (*i1 + i__) * a_dim1;
#line 151 "zsyswapr.f"
	    i__3 = *i1 + i__ + *i2 * a_dim1;
#line 151 "zsyswapr.f"
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 152 "zsyswapr.f"
	    i__2 = *i1 + i__ + *i2 * a_dim1;
#line 152 "zsyswapr.f"
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
#line 153 "zsyswapr.f"
	}

/*          third swap */
/*          - swap row I1 and I2 from I2+1 to N */
#line 157 "zsyswapr.f"
	i__1 = *n;
#line 157 "zsyswapr.f"
	for (i__ = *i2 + 1; i__ <= i__1; ++i__) {
#line 158 "zsyswapr.f"
	    i__2 = *i1 + i__ * a_dim1;
#line 158 "zsyswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 159 "zsyswapr.f"
	    i__2 = *i1 + i__ * a_dim1;
#line 159 "zsyswapr.f"
	    i__3 = *i2 + i__ * a_dim1;
#line 159 "zsyswapr.f"
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 160 "zsyswapr.f"
	    i__2 = *i2 + i__ * a_dim1;
#line 160 "zsyswapr.f"
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
#line 161 "zsyswapr.f"
	}

#line 163 "zsyswapr.f"
    } else {

/*         LOWER */
/*         first swap */
/*          - swap row I1 and I2 from I1 to I1-1 */
#line 168 "zsyswapr.f"
	i__1 = *i1 - 1;
#line 168 "zsyswapr.f"
	zswap_(&i__1, &a[*i1 + a_dim1], lda, &a[*i2 + a_dim1], lda);

/*         second swap : */
/*          - swap A(I1,I1) and A(I2,I2) */
/*          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1 */
#line 173 "zsyswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 173 "zsyswapr.f"
	tmp.r = a[i__1].r, tmp.i = a[i__1].i;
#line 174 "zsyswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 174 "zsyswapr.f"
	i__2 = *i2 + *i2 * a_dim1;
#line 174 "zsyswapr.f"
	a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 175 "zsyswapr.f"
	i__1 = *i2 + *i2 * a_dim1;
#line 175 "zsyswapr.f"
	a[i__1].r = tmp.r, a[i__1].i = tmp.i;

#line 177 "zsyswapr.f"
	i__1 = *i2 - *i1 - 1;
#line 177 "zsyswapr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "zsyswapr.f"
	    i__2 = *i1 + i__ + *i1 * a_dim1;
#line 178 "zsyswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 179 "zsyswapr.f"
	    i__2 = *i1 + i__ + *i1 * a_dim1;
#line 179 "zsyswapr.f"
	    i__3 = *i2 + (*i1 + i__) * a_dim1;
#line 179 "zsyswapr.f"
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 180 "zsyswapr.f"
	    i__2 = *i2 + (*i1 + i__) * a_dim1;
#line 180 "zsyswapr.f"
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
#line 181 "zsyswapr.f"
	}

/*         third swap */
/*          - swap col I1 and I2 from I2+1 to N */
#line 185 "zsyswapr.f"
	i__1 = *n;
#line 185 "zsyswapr.f"
	for (i__ = *i2 + 1; i__ <= i__1; ++i__) {
#line 186 "zsyswapr.f"
	    i__2 = i__ + *i1 * a_dim1;
#line 186 "zsyswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 187 "zsyswapr.f"
	    i__2 = i__ + *i1 * a_dim1;
#line 187 "zsyswapr.f"
	    i__3 = i__ + *i2 * a_dim1;
#line 187 "zsyswapr.f"
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 188 "zsyswapr.f"
	    i__2 = i__ + *i2 * a_dim1;
#line 188 "zsyswapr.f"
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
#line 189 "zsyswapr.f"
	}

#line 191 "zsyswapr.f"
    }
#line 192 "zsyswapr.f"
    return 0;
} /* zsyswapr_ */

