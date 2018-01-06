#line 1 "ssycon.f"
/* ssycon.f -- translated by f2c (version 20100827).
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

#line 1 "ssycon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SSYCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssycon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssycon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssycon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric matrix A using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by SSYTRF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by SSYTRF. */
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
/* >          as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* >          estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssycon_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int ssytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
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

#line 176 "ssycon.f"
    /* Parameter adjustments */
#line 176 "ssycon.f"
    a_dim1 = *lda;
#line 176 "ssycon.f"
    a_offset = 1 + a_dim1;
#line 176 "ssycon.f"
    a -= a_offset;
#line 176 "ssycon.f"
    --ipiv;
#line 176 "ssycon.f"
    --work;
#line 176 "ssycon.f"
    --iwork;
#line 176 "ssycon.f"

#line 176 "ssycon.f"
    /* Function Body */
#line 176 "ssycon.f"
    *info = 0;
#line 177 "ssycon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 178 "ssycon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 179 "ssycon.f"
	*info = -1;
#line 180 "ssycon.f"
    } else if (*n < 0) {
#line 181 "ssycon.f"
	*info = -2;
#line 182 "ssycon.f"
    } else if (*lda < max(1,*n)) {
#line 183 "ssycon.f"
	*info = -4;
#line 184 "ssycon.f"
    } else if (*anorm < 0.) {
#line 185 "ssycon.f"
	*info = -6;
#line 186 "ssycon.f"
    }
#line 187 "ssycon.f"
    if (*info != 0) {
#line 188 "ssycon.f"
	i__1 = -(*info);
#line 188 "ssycon.f"
	xerbla_("SSYCON", &i__1, (ftnlen)6);
#line 189 "ssycon.f"
	return 0;
#line 190 "ssycon.f"
    }

/*     Quick return if possible */

#line 194 "ssycon.f"
    *rcond = 0.;
#line 195 "ssycon.f"
    if (*n == 0) {
#line 196 "ssycon.f"
	*rcond = 1.;
#line 197 "ssycon.f"
	return 0;
#line 198 "ssycon.f"
    } else if (*anorm <= 0.) {
#line 199 "ssycon.f"
	return 0;
#line 200 "ssycon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 204 "ssycon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 208 "ssycon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 209 "ssycon.f"
	    if (ipiv[i__] > 0 && a[i__ + i__ * a_dim1] == 0.) {
#line 209 "ssycon.f"
		return 0;
#line 209 "ssycon.f"
	    }
#line 211 "ssycon.f"
/* L10: */
#line 211 "ssycon.f"
	}
#line 212 "ssycon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 216 "ssycon.f"
	i__1 = *n;
#line 216 "ssycon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 217 "ssycon.f"
	    if (ipiv[i__] > 0 && a[i__ + i__ * a_dim1] == 0.) {
#line 217 "ssycon.f"
		return 0;
#line 217 "ssycon.f"
	    }
#line 219 "ssycon.f"
/* L20: */
#line 219 "ssycon.f"
	}
#line 220 "ssycon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 224 "ssycon.f"
    kase = 0;
#line 225 "ssycon.f"
L30:
#line 226 "ssycon.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 227 "ssycon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 231 "ssycon.f"
	ssytrs_(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], n, 
		info, (ftnlen)1);
#line 232 "ssycon.f"
	goto L30;
#line 233 "ssycon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 237 "ssycon.f"
    if (ainvnm != 0.) {
#line 237 "ssycon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 237 "ssycon.f"
    }

#line 240 "ssycon.f"
    return 0;

/*     End of SSYCON */

} /* ssycon_ */

