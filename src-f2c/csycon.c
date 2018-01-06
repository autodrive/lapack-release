#line 1 "csycon.f"
/* csycon.f -- translated by f2c (version 20100827).
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

#line 1 "csycon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CSYCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csycon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csycon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csycon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex symmetric matrix A using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by CSYTRF. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSYTRF. */
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
/* >          as determined by CSYTRF. */
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
/* >          WORK is COMPLEX array, dimension (2*N) */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csycon_(char *uplo, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csytrs_(char *, integer *, integer *, 
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

#line 171 "csycon.f"
    /* Parameter adjustments */
#line 171 "csycon.f"
    a_dim1 = *lda;
#line 171 "csycon.f"
    a_offset = 1 + a_dim1;
#line 171 "csycon.f"
    a -= a_offset;
#line 171 "csycon.f"
    --ipiv;
#line 171 "csycon.f"
    --work;
#line 171 "csycon.f"

#line 171 "csycon.f"
    /* Function Body */
#line 171 "csycon.f"
    *info = 0;
#line 172 "csycon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 173 "csycon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "csycon.f"
	*info = -1;
#line 175 "csycon.f"
    } else if (*n < 0) {
#line 176 "csycon.f"
	*info = -2;
#line 177 "csycon.f"
    } else if (*lda < max(1,*n)) {
#line 178 "csycon.f"
	*info = -4;
#line 179 "csycon.f"
    } else if (*anorm < 0.) {
#line 180 "csycon.f"
	*info = -6;
#line 181 "csycon.f"
    }
#line 182 "csycon.f"
    if (*info != 0) {
#line 183 "csycon.f"
	i__1 = -(*info);
#line 183 "csycon.f"
	xerbla_("CSYCON", &i__1, (ftnlen)6);
#line 184 "csycon.f"
	return 0;
#line 185 "csycon.f"
    }

/*     Quick return if possible */

#line 189 "csycon.f"
    *rcond = 0.;
#line 190 "csycon.f"
    if (*n == 0) {
#line 191 "csycon.f"
	*rcond = 1.;
#line 192 "csycon.f"
	return 0;
#line 193 "csycon.f"
    } else if (*anorm <= 0.) {
#line 194 "csycon.f"
	return 0;
#line 195 "csycon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 199 "csycon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 203 "csycon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 204 "csycon.f"
	    i__1 = i__ + i__ * a_dim1;
#line 204 "csycon.f"
	    if (ipiv[i__] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 204 "csycon.f"
		return 0;
#line 204 "csycon.f"
	    }
#line 206 "csycon.f"
/* L10: */
#line 206 "csycon.f"
	}
#line 207 "csycon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 211 "csycon.f"
	i__1 = *n;
#line 211 "csycon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "csycon.f"
	    i__2 = i__ + i__ * a_dim1;
#line 212 "csycon.f"
	    if (ipiv[i__] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 212 "csycon.f"
		return 0;
#line 212 "csycon.f"
	    }
#line 214 "csycon.f"
/* L20: */
#line 214 "csycon.f"
	}
#line 215 "csycon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 219 "csycon.f"
    kase = 0;
#line 220 "csycon.f"
L30:
#line 221 "csycon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 222 "csycon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 226 "csycon.f"
	csytrs_(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], n, 
		info, (ftnlen)1);
#line 227 "csycon.f"
	goto L30;
#line 228 "csycon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 232 "csycon.f"
    if (ainvnm != 0.) {
#line 232 "csycon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 232 "csycon.f"
    }

#line 235 "csycon.f"
    return 0;

/*     End of CSYCON */

} /* csycon_ */

