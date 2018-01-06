#line 1 "slauum.f"
/* slauum.f -- translated by f2c (version 20100827).
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

#line 1 "slauum.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b15 = 1.;

/* > \brief \b SLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (b
locked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAUUM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slauum.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slauum.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slauum.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAUUM( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAUUM computes the product U * U**T or L**T * L, where the triangular */
/* > factor U or L is stored in the upper or lower triangular part of */
/* > the array A. */
/* > */
/* > If UPLO = 'U' or 'u' then the upper triangle of the result is stored, */
/* > overwriting the factor U in A. */
/* > If UPLO = 'L' or 'l' then the lower triangle of the result is stored, */
/* > overwriting the factor L in A. */
/* > */
/* > This is the blocked form of the algorithm, calling Level 3 BLAS. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the triangular factor U or L. */
/* >          On exit, if UPLO = 'U', the upper triangle of A is */
/* >          overwritten with the upper triangle of the product U * U**T; */
/* >          if UPLO = 'L', the lower triangle of A is overwritten with */
/* >          the lower triangle of the product L**T * L. */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slauum_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), ssyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     slauu2_(char *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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

#line 143 "slauum.f"
    /* Parameter adjustments */
#line 143 "slauum.f"
    a_dim1 = *lda;
#line 143 "slauum.f"
    a_offset = 1 + a_dim1;
#line 143 "slauum.f"
    a -= a_offset;
#line 143 "slauum.f"

#line 143 "slauum.f"
    /* Function Body */
#line 143 "slauum.f"
    *info = 0;
#line 144 "slauum.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 145 "slauum.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 146 "slauum.f"
	*info = -1;
#line 147 "slauum.f"
    } else if (*n < 0) {
#line 148 "slauum.f"
	*info = -2;
#line 149 "slauum.f"
    } else if (*lda < max(1,*n)) {
#line 150 "slauum.f"
	*info = -4;
#line 151 "slauum.f"
    }
#line 152 "slauum.f"
    if (*info != 0) {
#line 153 "slauum.f"
	i__1 = -(*info);
#line 153 "slauum.f"
	xerbla_("SLAUUM", &i__1, (ftnlen)6);
#line 154 "slauum.f"
	return 0;
#line 155 "slauum.f"
    }

/*     Quick return if possible */

#line 159 "slauum.f"
    if (*n == 0) {
#line 159 "slauum.f"
	return 0;
#line 159 "slauum.f"
    }

/*     Determine the block size for this environment. */

#line 164 "slauum.f"
    nb = ilaenv_(&c__1, "SLAUUM", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 166 "slauum.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 170 "slauum.f"
	slauu2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
#line 171 "slauum.f"
    } else {

/*        Use blocked code */

#line 175 "slauum.f"
	if (upper) {

/*           Compute the product U * U**T. */

#line 179 "slauum.f"
	    i__1 = *n;
#line 179 "slauum.f"
	    i__2 = nb;
#line 179 "slauum.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 180 "slauum.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 180 "slauum.f"
		ib = min(i__3,i__4);
#line 181 "slauum.f"
		i__3 = i__ - 1;
#line 181 "slauum.f"
		strmm_("Right", "Upper", "Transpose", "Non-unit", &i__3, &ib, 
			&c_b15, &a[i__ + i__ * a_dim1], lda, &a[i__ * a_dim1 
			+ 1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8)
			;
#line 184 "slauum.f"
		slauu2_("Upper", &ib, &a[i__ + i__ * a_dim1], lda, info, (
			ftnlen)5);
#line 185 "slauum.f"
		if (i__ + ib <= *n) {
#line 186 "slauum.f"
		    i__3 = i__ - 1;
#line 186 "slauum.f"
		    i__4 = *n - i__ - ib + 1;
#line 186 "slauum.f"
		    sgemm_("No transpose", "Transpose", &i__3, &ib, &i__4, &
			    c_b15, &a[(i__ + ib) * a_dim1 + 1], lda, &a[i__ + 
			    (i__ + ib) * a_dim1], lda, &c_b15, &a[i__ * 
			    a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 189 "slauum.f"
		    i__3 = *n - i__ - ib + 1;
#line 189 "slauum.f"
		    ssyrk_("Upper", "No transpose", &ib, &i__3, &c_b15, &a[
			    i__ + (i__ + ib) * a_dim1], lda, &c_b15, &a[i__ + 
			    i__ * a_dim1], lda, (ftnlen)5, (ftnlen)12);
#line 192 "slauum.f"
		}
#line 193 "slauum.f"
/* L10: */
#line 193 "slauum.f"
	    }
#line 194 "slauum.f"
	} else {

/*           Compute the product L**T * L. */

#line 198 "slauum.f"
	    i__2 = *n;
#line 198 "slauum.f"
	    i__1 = nb;
#line 198 "slauum.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 199 "slauum.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 199 "slauum.f"
		ib = min(i__3,i__4);
#line 200 "slauum.f"
		i__3 = i__ - 1;
#line 200 "slauum.f"
		strmm_("Left", "Lower", "Transpose", "Non-unit", &ib, &i__3, &
			c_b15, &a[i__ + i__ * a_dim1], lda, &a[i__ + a_dim1], 
			lda, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 202 "slauum.f"
		slauu2_("Lower", &ib, &a[i__ + i__ * a_dim1], lda, info, (
			ftnlen)5);
#line 203 "slauum.f"
		if (i__ + ib <= *n) {
#line 204 "slauum.f"
		    i__3 = i__ - 1;
#line 204 "slauum.f"
		    i__4 = *n - i__ - ib + 1;
#line 204 "slauum.f"
		    sgemm_("Transpose", "No transpose", &ib, &i__3, &i__4, &
			    c_b15, &a[i__ + ib + i__ * a_dim1], lda, &a[i__ + 
			    ib + a_dim1], lda, &c_b15, &a[i__ + a_dim1], lda, 
			    (ftnlen)9, (ftnlen)12);
#line 207 "slauum.f"
		    i__3 = *n - i__ - ib + 1;
#line 207 "slauum.f"
		    ssyrk_("Lower", "Transpose", &ib, &i__3, &c_b15, &a[i__ + 
			    ib + i__ * a_dim1], lda, &c_b15, &a[i__ + i__ * 
			    a_dim1], lda, (ftnlen)5, (ftnlen)9);
#line 209 "slauum.f"
		}
#line 210 "slauum.f"
/* L20: */
#line 210 "slauum.f"
	    }
#line 211 "slauum.f"
	}
#line 212 "slauum.f"
    }

#line 214 "slauum.f"
    return 0;

/*     End of SLAUUM */

} /* slauum_ */

