#line 1 "clauum.f"
/* clauum.f -- translated by f2c (version 20100827).
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

#line 1 "clauum.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b21 = 1.;

/* > \brief \b CLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (b
locked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAUUM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauum.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauum.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauum.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAUUM( UPLO, N, A, LDA, INFO ) */

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
/* > CLAUUM computes the product U * U**H or L**H * L, where the triangular */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clauum_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, nb;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cherk_(char *, char *, integer *, 
	    integer *, doublereal *, doublecomplex *, integer *, doublereal *,
	     doublecomplex *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int clauu2_(char *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 145 "clauum.f"
    /* Parameter adjustments */
#line 145 "clauum.f"
    a_dim1 = *lda;
#line 145 "clauum.f"
    a_offset = 1 + a_dim1;
#line 145 "clauum.f"
    a -= a_offset;
#line 145 "clauum.f"

#line 145 "clauum.f"
    /* Function Body */
#line 145 "clauum.f"
    *info = 0;
#line 146 "clauum.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 147 "clauum.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 148 "clauum.f"
	*info = -1;
#line 149 "clauum.f"
    } else if (*n < 0) {
#line 150 "clauum.f"
	*info = -2;
#line 151 "clauum.f"
    } else if (*lda < max(1,*n)) {
#line 152 "clauum.f"
	*info = -4;
#line 153 "clauum.f"
    }
#line 154 "clauum.f"
    if (*info != 0) {
#line 155 "clauum.f"
	i__1 = -(*info);
#line 155 "clauum.f"
	xerbla_("CLAUUM", &i__1, (ftnlen)6);
#line 156 "clauum.f"
	return 0;
#line 157 "clauum.f"
    }

/*     Quick return if possible */

#line 161 "clauum.f"
    if (*n == 0) {
#line 161 "clauum.f"
	return 0;
#line 161 "clauum.f"
    }

/*     Determine the block size for this environment. */

#line 166 "clauum.f"
    nb = ilaenv_(&c__1, "CLAUUM", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 168 "clauum.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 172 "clauum.f"
	clauu2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
#line 173 "clauum.f"
    } else {

/*        Use blocked code */

#line 177 "clauum.f"
	if (upper) {

/*           Compute the product U * U**H. */

#line 181 "clauum.f"
	    i__1 = *n;
#line 181 "clauum.f"
	    i__2 = nb;
#line 181 "clauum.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 182 "clauum.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 182 "clauum.f"
		ib = min(i__3,i__4);
#line 183 "clauum.f"
		i__3 = i__ - 1;
#line 183 "clauum.f"
		ctrmm_("Right", "Upper", "Conjugate transpose", "Non-unit", &
			i__3, &ib, &c_b1, &a[i__ + i__ * a_dim1], lda, &a[i__ 
			* a_dim1 + 1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)19, 
			(ftnlen)8);
#line 186 "clauum.f"
		clauu2_("Upper", &ib, &a[i__ + i__ * a_dim1], lda, info, (
			ftnlen)5);
#line 187 "clauum.f"
		if (i__ + ib <= *n) {
#line 188 "clauum.f"
		    i__3 = i__ - 1;
#line 188 "clauum.f"
		    i__4 = *n - i__ - ib + 1;
#line 188 "clauum.f"
		    cgemm_("No transpose", "Conjugate transpose", &i__3, &ib, 
			    &i__4, &c_b1, &a[(i__ + ib) * a_dim1 + 1], lda, &
			    a[i__ + (i__ + ib) * a_dim1], lda, &c_b1, &a[i__ *
			     a_dim1 + 1], lda, (ftnlen)12, (ftnlen)19);
#line 192 "clauum.f"
		    i__3 = *n - i__ - ib + 1;
#line 192 "clauum.f"
		    cherk_("Upper", "No transpose", &ib, &i__3, &c_b21, &a[
			    i__ + (i__ + ib) * a_dim1], lda, &c_b21, &a[i__ + 
			    i__ * a_dim1], lda, (ftnlen)5, (ftnlen)12);
#line 195 "clauum.f"
		}
#line 196 "clauum.f"
/* L10: */
#line 196 "clauum.f"
	    }
#line 197 "clauum.f"
	} else {

/*           Compute the product L**H * L. */

#line 201 "clauum.f"
	    i__2 = *n;
#line 201 "clauum.f"
	    i__1 = nb;
#line 201 "clauum.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 202 "clauum.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 202 "clauum.f"
		ib = min(i__3,i__4);
#line 203 "clauum.f"
		i__3 = i__ - 1;
#line 203 "clauum.f"
		ctrmm_("Left", "Lower", "Conjugate transpose", "Non-unit", &
			ib, &i__3, &c_b1, &a[i__ + i__ * a_dim1], lda, &a[i__ 
			+ a_dim1], lda, (ftnlen)4, (ftnlen)5, (ftnlen)19, (
			ftnlen)8);
#line 206 "clauum.f"
		clauu2_("Lower", &ib, &a[i__ + i__ * a_dim1], lda, info, (
			ftnlen)5);
#line 207 "clauum.f"
		if (i__ + ib <= *n) {
#line 208 "clauum.f"
		    i__3 = i__ - 1;
#line 208 "clauum.f"
		    i__4 = *n - i__ - ib + 1;
#line 208 "clauum.f"
		    cgemm_("Conjugate transpose", "No transpose", &ib, &i__3, 
			    &i__4, &c_b1, &a[i__ + ib + i__ * a_dim1], lda, &
			    a[i__ + ib + a_dim1], lda, &c_b1, &a[i__ + a_dim1]
			    , lda, (ftnlen)19, (ftnlen)12);
#line 211 "clauum.f"
		    i__3 = *n - i__ - ib + 1;
#line 211 "clauum.f"
		    cherk_("Lower", "Conjugate transpose", &ib, &i__3, &c_b21,
			     &a[i__ + ib + i__ * a_dim1], lda, &c_b21, &a[i__ 
			    + i__ * a_dim1], lda, (ftnlen)5, (ftnlen)19);
#line 214 "clauum.f"
		}
#line 215 "clauum.f"
/* L20: */
#line 215 "clauum.f"
	    }
#line 216 "clauum.f"
	}
#line 217 "clauum.f"
    }

#line 219 "clauum.f"
    return 0;

/*     End of CLAUUM */

} /* clauum_ */

