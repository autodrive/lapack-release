#line 1 "cheswapr.f"
/* cheswapr.f -- translated by f2c (version 20100827).
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

#line 1 "cheswapr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHESWAPR applies an elementary permutation on the rows and columns of a Hermitian matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHESWAPR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheswap
r.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheswap
r.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheswap
r.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHESWAPR( UPLO, N, A, LDA, I1, I2) */

/*       .. Scalar Arguments .. */
/*       CHARACTER        UPLO */
/*       INTEGER          I1, I2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX          A( LDA, N ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHESWAPR applies an elementary permutation on the rows and the columns of */
/* > a hermitian matrix. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the NB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by CSYTRF. */
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

/* > \ingroup complexHEauxiliary */

/*  ===================================================================== */
/* Subroutine */ int cheswapr_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *i1, integer *i2, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublecomplex tmp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper;


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

#line 134 "cheswapr.f"
    /* Parameter adjustments */
#line 134 "cheswapr.f"
    a_dim1 = *lda;
#line 134 "cheswapr.f"
    a_offset = 1 + a_dim1;
#line 134 "cheswapr.f"
    a -= a_offset;
#line 134 "cheswapr.f"

#line 134 "cheswapr.f"
    /* Function Body */
#line 134 "cheswapr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 135 "cheswapr.f"
    if (upper) {

/*         UPPER */
/*         first swap */
/*          - swap column I1 and I2 from I1 to I1-1 */
#line 140 "cheswapr.f"
	i__1 = *i1 - 1;
#line 140 "cheswapr.f"
	cswap_(&i__1, &a[*i1 * a_dim1 + 1], &c__1, &a[*i2 * a_dim1 + 1], &
		c__1);

/*          second swap : */
/*          - swap A(I1,I1) and A(I2,I2) */
/*          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1 */
/*          - swap A(I2,I1) and A(I1,I2) */
#line 147 "cheswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 147 "cheswapr.f"
	tmp.r = a[i__1].r, tmp.i = a[i__1].i;
#line 148 "cheswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 148 "cheswapr.f"
	i__2 = *i2 + *i2 * a_dim1;
#line 148 "cheswapr.f"
	a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 149 "cheswapr.f"
	i__1 = *i2 + *i2 * a_dim1;
#line 149 "cheswapr.f"
	a[i__1].r = tmp.r, a[i__1].i = tmp.i;

#line 151 "cheswapr.f"
	i__1 = *i2 - *i1 - 1;
#line 151 "cheswapr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "cheswapr.f"
	    i__2 = *i1 + (*i1 + i__) * a_dim1;
#line 152 "cheswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 153 "cheswapr.f"
	    i__2 = *i1 + (*i1 + i__) * a_dim1;
#line 153 "cheswapr.f"
	    d_cnjg(&z__1, &a[*i1 + i__ + *i2 * a_dim1]);
#line 153 "cheswapr.f"
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 154 "cheswapr.f"
	    i__2 = *i1 + i__ + *i2 * a_dim1;
#line 154 "cheswapr.f"
	    d_cnjg(&z__1, &tmp);
#line 154 "cheswapr.f"
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 155 "cheswapr.f"
	}

#line 157 "cheswapr.f"
	i__1 = *i1 + *i2 * a_dim1;
#line 157 "cheswapr.f"
	d_cnjg(&z__1, &a[*i1 + *i2 * a_dim1]);
#line 157 "cheswapr.f"
	a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*          third swap */
/*          - swap row I1 and I2 from I2+1 to N */
#line 162 "cheswapr.f"
	i__1 = *n;
#line 162 "cheswapr.f"
	for (i__ = *i2 + 1; i__ <= i__1; ++i__) {
#line 163 "cheswapr.f"
	    i__2 = *i1 + i__ * a_dim1;
#line 163 "cheswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 164 "cheswapr.f"
	    i__2 = *i1 + i__ * a_dim1;
#line 164 "cheswapr.f"
	    i__3 = *i2 + i__ * a_dim1;
#line 164 "cheswapr.f"
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 165 "cheswapr.f"
	    i__2 = *i2 + i__ * a_dim1;
#line 165 "cheswapr.f"
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
#line 166 "cheswapr.f"
	}

#line 168 "cheswapr.f"
    } else {

/*         LOWER */
/*         first swap */
/*          - swap row I1 and I2 from 1 to I1-1 */
#line 173 "cheswapr.f"
	i__1 = *i1 - 1;
#line 173 "cheswapr.f"
	cswap_(&i__1, &a[*i1 + a_dim1], lda, &a[*i2 + a_dim1], lda);

/*         second swap : */
/*          - swap A(I1,I1) and A(I2,I2) */
/*          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1 */
/*          - swap A(I2,I1) and A(I1,I2) */
#line 180 "cheswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 180 "cheswapr.f"
	tmp.r = a[i__1].r, tmp.i = a[i__1].i;
#line 181 "cheswapr.f"
	i__1 = *i1 + *i1 * a_dim1;
#line 181 "cheswapr.f"
	i__2 = *i2 + *i2 * a_dim1;
#line 181 "cheswapr.f"
	a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 182 "cheswapr.f"
	i__1 = *i2 + *i2 * a_dim1;
#line 182 "cheswapr.f"
	a[i__1].r = tmp.r, a[i__1].i = tmp.i;

#line 184 "cheswapr.f"
	i__1 = *i2 - *i1 - 1;
#line 184 "cheswapr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "cheswapr.f"
	    i__2 = *i1 + i__ + *i1 * a_dim1;
#line 185 "cheswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 186 "cheswapr.f"
	    i__2 = *i1 + i__ + *i1 * a_dim1;
#line 186 "cheswapr.f"
	    d_cnjg(&z__1, &a[*i2 + (*i1 + i__) * a_dim1]);
#line 186 "cheswapr.f"
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 187 "cheswapr.f"
	    i__2 = *i2 + (*i1 + i__) * a_dim1;
#line 187 "cheswapr.f"
	    d_cnjg(&z__1, &tmp);
#line 187 "cheswapr.f"
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 188 "cheswapr.f"
	}

#line 190 "cheswapr.f"
	i__1 = *i2 + *i1 * a_dim1;
#line 190 "cheswapr.f"
	d_cnjg(&z__1, &a[*i2 + *i1 * a_dim1]);
#line 190 "cheswapr.f"
	a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*         third swap */
/*          - swap col I1 and I2 from I2+1 to N */
#line 194 "cheswapr.f"
	i__1 = *n;
#line 194 "cheswapr.f"
	for (i__ = *i2 + 1; i__ <= i__1; ++i__) {
#line 195 "cheswapr.f"
	    i__2 = i__ + *i1 * a_dim1;
#line 195 "cheswapr.f"
	    tmp.r = a[i__2].r, tmp.i = a[i__2].i;
#line 196 "cheswapr.f"
	    i__2 = i__ + *i1 * a_dim1;
#line 196 "cheswapr.f"
	    i__3 = i__ + *i2 * a_dim1;
#line 196 "cheswapr.f"
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 197 "cheswapr.f"
	    i__2 = i__ + *i2 * a_dim1;
#line 197 "cheswapr.f"
	    a[i__2].r = tmp.r, a[i__2].i = tmp.i;
#line 198 "cheswapr.f"
	}

#line 200 "cheswapr.f"
    }
#line 202 "cheswapr.f"
    return 0;
} /* cheswapr_ */

