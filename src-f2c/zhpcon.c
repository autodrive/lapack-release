#line 1 "zhpcon.f"
/* zhpcon.f -- translated by f2c (version 20100827).
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

#line 1 "zhpcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPCON estimates the reciprocal of the condition number of a complex */
/* > Hermitian packed matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by ZHPTRF. */
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
/* >          = 'U':  Upper triangular, form is A = U*D*U**H; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZHPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZHPTRF. */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhpcon_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, ip, kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int zhptrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 161 "zhpcon.f"
    /* Parameter adjustments */
#line 161 "zhpcon.f"
    --work;
#line 161 "zhpcon.f"
    --ipiv;
#line 161 "zhpcon.f"
    --ap;
#line 161 "zhpcon.f"

#line 161 "zhpcon.f"
    /* Function Body */
#line 161 "zhpcon.f"
    *info = 0;
#line 162 "zhpcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "zhpcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "zhpcon.f"
	*info = -1;
#line 165 "zhpcon.f"
    } else if (*n < 0) {
#line 166 "zhpcon.f"
	*info = -2;
#line 167 "zhpcon.f"
    } else if (*anorm < 0.) {
#line 168 "zhpcon.f"
	*info = -5;
#line 169 "zhpcon.f"
    }
#line 170 "zhpcon.f"
    if (*info != 0) {
#line 171 "zhpcon.f"
	i__1 = -(*info);
#line 171 "zhpcon.f"
	xerbla_("ZHPCON", &i__1, (ftnlen)6);
#line 172 "zhpcon.f"
	return 0;
#line 173 "zhpcon.f"
    }

/*     Quick return if possible */

#line 177 "zhpcon.f"
    *rcond = 0.;
#line 178 "zhpcon.f"
    if (*n == 0) {
#line 179 "zhpcon.f"
	*rcond = 1.;
#line 180 "zhpcon.f"
	return 0;
#line 181 "zhpcon.f"
    } else if (*anorm <= 0.) {
#line 182 "zhpcon.f"
	return 0;
#line 183 "zhpcon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 187 "zhpcon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 191 "zhpcon.f"
	ip = *n * (*n + 1) / 2;
#line 192 "zhpcon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 193 "zhpcon.f"
	    i__1 = ip;
#line 193 "zhpcon.f"
	    if (ipiv[i__] > 0 && (ap[i__1].r == 0. && ap[i__1].i == 0.)) {
#line 193 "zhpcon.f"
		return 0;
#line 193 "zhpcon.f"
	    }
#line 195 "zhpcon.f"
	    ip -= i__;
#line 196 "zhpcon.f"
/* L10: */
#line 196 "zhpcon.f"
	}
#line 197 "zhpcon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 201 "zhpcon.f"
	ip = 1;
#line 202 "zhpcon.f"
	i__1 = *n;
#line 202 "zhpcon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "zhpcon.f"
	    i__2 = ip;
#line 203 "zhpcon.f"
	    if (ipiv[i__] > 0 && (ap[i__2].r == 0. && ap[i__2].i == 0.)) {
#line 203 "zhpcon.f"
		return 0;
#line 203 "zhpcon.f"
	    }
#line 205 "zhpcon.f"
	    ip = ip + *n - i__ + 1;
#line 206 "zhpcon.f"
/* L20: */
#line 206 "zhpcon.f"
	}
#line 207 "zhpcon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 211 "zhpcon.f"
    kase = 0;
#line 212 "zhpcon.f"
L30:
#line 213 "zhpcon.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 214 "zhpcon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**H) or inv(U*D*U**H). */

#line 218 "zhpcon.f"
	zhptrs_(uplo, n, &c__1, &ap[1], &ipiv[1], &work[1], n, info, (ftnlen)
		1);
#line 219 "zhpcon.f"
	goto L30;
#line 220 "zhpcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 224 "zhpcon.f"
    if (ainvnm != 0.) {
#line 224 "zhpcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 224 "zhpcon.f"
    }

#line 227 "zhpcon.f"
    return 0;

/*     End of ZHPCON */

} /* zhpcon_ */

