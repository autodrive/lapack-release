#line 1 "zsycon.f"
/* zsycon.f -- translated by f2c (version 20100827).
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

#line 1 "zsycon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZSYCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsycon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsycon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsycon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex symmetric matrix A using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by ZSYTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* >          estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsycon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int zsytrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 171 "zsycon.f"
    /* Parameter adjustments */
#line 171 "zsycon.f"
    a_dim1 = *lda;
#line 171 "zsycon.f"
    a_offset = 1 + a_dim1;
#line 171 "zsycon.f"
    a -= a_offset;
#line 171 "zsycon.f"
    --ipiv;
#line 171 "zsycon.f"
    --work;
#line 171 "zsycon.f"

#line 171 "zsycon.f"
    /* Function Body */
#line 171 "zsycon.f"
    *info = 0;
#line 172 "zsycon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 173 "zsycon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "zsycon.f"
	*info = -1;
#line 175 "zsycon.f"
    } else if (*n < 0) {
#line 176 "zsycon.f"
	*info = -2;
#line 177 "zsycon.f"
    } else if (*lda < max(1,*n)) {
#line 178 "zsycon.f"
	*info = -4;
#line 179 "zsycon.f"
    } else if (*anorm < 0.) {
#line 180 "zsycon.f"
	*info = -6;
#line 181 "zsycon.f"
    }
#line 182 "zsycon.f"
    if (*info != 0) {
#line 183 "zsycon.f"
	i__1 = -(*info);
#line 183 "zsycon.f"
	xerbla_("ZSYCON", &i__1, (ftnlen)6);
#line 184 "zsycon.f"
	return 0;
#line 185 "zsycon.f"
    }

/*     Quick return if possible */

#line 189 "zsycon.f"
    *rcond = 0.;
#line 190 "zsycon.f"
    if (*n == 0) {
#line 191 "zsycon.f"
	*rcond = 1.;
#line 192 "zsycon.f"
	return 0;
#line 193 "zsycon.f"
    } else if (*anorm <= 0.) {
#line 194 "zsycon.f"
	return 0;
#line 195 "zsycon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 199 "zsycon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 203 "zsycon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 204 "zsycon.f"
	    i__1 = i__ + i__ * a_dim1;
#line 204 "zsycon.f"
	    if (ipiv[i__] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 204 "zsycon.f"
		return 0;
#line 204 "zsycon.f"
	    }
#line 206 "zsycon.f"
/* L10: */
#line 206 "zsycon.f"
	}
#line 207 "zsycon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 211 "zsycon.f"
	i__1 = *n;
#line 211 "zsycon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "zsycon.f"
	    i__2 = i__ + i__ * a_dim1;
#line 212 "zsycon.f"
	    if (ipiv[i__] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 212 "zsycon.f"
		return 0;
#line 212 "zsycon.f"
	    }
#line 214 "zsycon.f"
/* L20: */
#line 214 "zsycon.f"
	}
#line 215 "zsycon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 219 "zsycon.f"
    kase = 0;
#line 220 "zsycon.f"
L30:
#line 221 "zsycon.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 222 "zsycon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 226 "zsycon.f"
	zsytrs_(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], n, 
		info, (ftnlen)1);
#line 227 "zsycon.f"
	goto L30;
#line 228 "zsycon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 232 "zsycon.f"
    if (ainvnm != 0.) {
#line 232 "zsycon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 232 "zsycon.f"
    }

#line 235 "zsycon.f"
    return 0;

/*     End of ZSYCON */

} /* zsycon_ */

