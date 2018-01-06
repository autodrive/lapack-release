#line 1 "cspcon.f"
/* cspcon.f -- translated by f2c (version 20100827).
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

#line 1 "cspcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CSPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cspcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cspcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cspcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex symmetric packed matrix A using the */
/* > factorization A = U*D*U**T or A = L*D*L**T computed by CSPTRF. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by CSPTRF. */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cspcon_(char *uplo, integer *n, doublecomplex *ap, 
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
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csptrs_(char *, integer *, integer *, 
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

#line 161 "cspcon.f"
    /* Parameter adjustments */
#line 161 "cspcon.f"
    --work;
#line 161 "cspcon.f"
    --ipiv;
#line 161 "cspcon.f"
    --ap;
#line 161 "cspcon.f"

#line 161 "cspcon.f"
    /* Function Body */
#line 161 "cspcon.f"
    *info = 0;
#line 162 "cspcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "cspcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "cspcon.f"
	*info = -1;
#line 165 "cspcon.f"
    } else if (*n < 0) {
#line 166 "cspcon.f"
	*info = -2;
#line 167 "cspcon.f"
    } else if (*anorm < 0.) {
#line 168 "cspcon.f"
	*info = -5;
#line 169 "cspcon.f"
    }
#line 170 "cspcon.f"
    if (*info != 0) {
#line 171 "cspcon.f"
	i__1 = -(*info);
#line 171 "cspcon.f"
	xerbla_("CSPCON", &i__1, (ftnlen)6);
#line 172 "cspcon.f"
	return 0;
#line 173 "cspcon.f"
    }

/*     Quick return if possible */

#line 177 "cspcon.f"
    *rcond = 0.;
#line 178 "cspcon.f"
    if (*n == 0) {
#line 179 "cspcon.f"
	*rcond = 1.;
#line 180 "cspcon.f"
	return 0;
#line 181 "cspcon.f"
    } else if (*anorm <= 0.) {
#line 182 "cspcon.f"
	return 0;
#line 183 "cspcon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 187 "cspcon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 191 "cspcon.f"
	ip = *n * (*n + 1) / 2;
#line 192 "cspcon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 193 "cspcon.f"
	    i__1 = ip;
#line 193 "cspcon.f"
	    if (ipiv[i__] > 0 && (ap[i__1].r == 0. && ap[i__1].i == 0.)) {
#line 193 "cspcon.f"
		return 0;
#line 193 "cspcon.f"
	    }
#line 195 "cspcon.f"
	    ip -= i__;
#line 196 "cspcon.f"
/* L10: */
#line 196 "cspcon.f"
	}
#line 197 "cspcon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 201 "cspcon.f"
	ip = 1;
#line 202 "cspcon.f"
	i__1 = *n;
#line 202 "cspcon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "cspcon.f"
	    i__2 = ip;
#line 203 "cspcon.f"
	    if (ipiv[i__] > 0 && (ap[i__2].r == 0. && ap[i__2].i == 0.)) {
#line 203 "cspcon.f"
		return 0;
#line 203 "cspcon.f"
	    }
#line 205 "cspcon.f"
	    ip = ip + *n - i__ + 1;
#line 206 "cspcon.f"
/* L20: */
#line 206 "cspcon.f"
	}
#line 207 "cspcon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 211 "cspcon.f"
    kase = 0;
#line 212 "cspcon.f"
L30:
#line 213 "cspcon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 214 "cspcon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 218 "cspcon.f"
	csptrs_(uplo, n, &c__1, &ap[1], &ipiv[1], &work[1], n, info, (ftnlen)
		1);
#line 219 "cspcon.f"
	goto L30;
#line 220 "cspcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 224 "cspcon.f"
    if (ainvnm != 0.) {
#line 224 "cspcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 224 "cspcon.f"
    }

#line 227 "cspcon.f"
    return 0;

/*     End of CSPCON */

} /* cspcon_ */

