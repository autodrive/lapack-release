#line 1 "zpotri.f"
/* zpotri.f -- translated by f2c (version 20100827).
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

#line 1 "zpotri.f"
/* > \brief \b ZPOTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPOTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpotri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpotri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpotri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPOTRI( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPOTRI computes the inverse of a complex Hermitian positive definite */
/* > matrix A using the Cholesky factorization A = U**H*U or A = L*L**H */
/* > computed by ZPOTRF. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H, as computed by */
/* >          ZPOTRF. */
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

/* > \date November 2011 */

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpotri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlauum_(
	    char *, integer *, doublecomplex *, integer *, integer *, ftnlen),
	     ztrtri_(char *, char *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 127 "zpotri.f"
    /* Parameter adjustments */
#line 127 "zpotri.f"
    a_dim1 = *lda;
#line 127 "zpotri.f"
    a_offset = 1 + a_dim1;
#line 127 "zpotri.f"
    a -= a_offset;
#line 127 "zpotri.f"

#line 127 "zpotri.f"
    /* Function Body */
#line 127 "zpotri.f"
    *info = 0;
#line 128 "zpotri.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 129 "zpotri.f"
	*info = -1;
#line 130 "zpotri.f"
    } else if (*n < 0) {
#line 131 "zpotri.f"
	*info = -2;
#line 132 "zpotri.f"
    } else if (*lda < max(1,*n)) {
#line 133 "zpotri.f"
	*info = -4;
#line 134 "zpotri.f"
    }
#line 135 "zpotri.f"
    if (*info != 0) {
#line 136 "zpotri.f"
	i__1 = -(*info);
#line 136 "zpotri.f"
	xerbla_("ZPOTRI", &i__1, (ftnlen)6);
#line 137 "zpotri.f"
	return 0;
#line 138 "zpotri.f"
    }

/*     Quick return if possible */

#line 142 "zpotri.f"
    if (*n == 0) {
#line 142 "zpotri.f"
	return 0;
#line 142 "zpotri.f"
    }

/*     Invert the triangular Cholesky factor U or L. */

#line 147 "zpotri.f"
    ztrtri_(uplo, "Non-unit", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)
	    8);
#line 148 "zpotri.f"
    if (*info > 0) {
#line 148 "zpotri.f"
	return 0;
#line 148 "zpotri.f"
    }

/*     Form inv(U) * inv(U)**H or inv(L)**H * inv(L). */

#line 153 "zpotri.f"
    zlauum_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);

#line 155 "zpotri.f"
    return 0;

/*     End of ZPOTRI */

} /* zpotri_ */

