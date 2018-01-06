#line 1 "cpotri.f"
/* cpotri.f -- translated by f2c (version 20100827).
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

#line 1 "cpotri.f"
/* > \brief \b CPOTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPOTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPOTRI( UPLO, N, A, LDA, INFO ) */

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
/* > CPOTRI computes the inverse of a complex Hermitian positive definite */
/* > matrix A using the Cholesky factorization A = U**H*U or A = L*L**H */
/* > computed by CPOTRF. */
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
/* >          On entry, the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H, as computed by */
/* >          CPOTRF. */
/* >          On exit, the upper or lower triangle of the (Hermitian) */
/* >          inverse of A, overwriting the input factor U or L. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the (i,i) element of the factor U or L is */
/* >                zero, and the inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpotri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), clauum_(
	    char *, integer *, doublecomplex *, integer *, integer *, ftnlen),
	     ctrtri_(char *, char *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 127 "cpotri.f"
    /* Parameter adjustments */
#line 127 "cpotri.f"
    a_dim1 = *lda;
#line 127 "cpotri.f"
    a_offset = 1 + a_dim1;
#line 127 "cpotri.f"
    a -= a_offset;
#line 127 "cpotri.f"

#line 127 "cpotri.f"
    /* Function Body */
#line 127 "cpotri.f"
    *info = 0;
#line 128 "cpotri.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 129 "cpotri.f"
	*info = -1;
#line 130 "cpotri.f"
    } else if (*n < 0) {
#line 131 "cpotri.f"
	*info = -2;
#line 132 "cpotri.f"
    } else if (*lda < max(1,*n)) {
#line 133 "cpotri.f"
	*info = -4;
#line 134 "cpotri.f"
    }
#line 135 "cpotri.f"
    if (*info != 0) {
#line 136 "cpotri.f"
	i__1 = -(*info);
#line 136 "cpotri.f"
	xerbla_("CPOTRI", &i__1, (ftnlen)6);
#line 137 "cpotri.f"
	return 0;
#line 138 "cpotri.f"
    }

/*     Quick return if possible */

#line 142 "cpotri.f"
    if (*n == 0) {
#line 142 "cpotri.f"
	return 0;
#line 142 "cpotri.f"
    }

/*     Invert the triangular Cholesky factor U or L. */

#line 147 "cpotri.f"
    ctrtri_(uplo, "Non-unit", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)
	    8);
#line 148 "cpotri.f"
    if (*info > 0) {
#line 148 "cpotri.f"
	return 0;
#line 148 "cpotri.f"
    }

/*     Form inv(U) * inv(U)**H or inv(L)**H * inv(L). */

#line 153 "cpotri.f"
    clauum_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);

#line 155 "cpotri.f"
    return 0;

/*     End of CPOTRI */

} /* cpotri_ */

