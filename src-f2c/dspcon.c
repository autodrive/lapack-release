#line 1 "dspcon.f"
/* dspcon.f -- translated by f2c (version 20100827).
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

#line 1 "dspcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric packed matrix A using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by DSPTRF. */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSPTRF. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dspcon_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer 
	*iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ip, kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int dsptrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
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

#line 168 "dspcon.f"
    /* Parameter adjustments */
#line 168 "dspcon.f"
    --iwork;
#line 168 "dspcon.f"
    --work;
#line 168 "dspcon.f"
    --ipiv;
#line 168 "dspcon.f"
    --ap;
#line 168 "dspcon.f"

#line 168 "dspcon.f"
    /* Function Body */
#line 168 "dspcon.f"
    *info = 0;
#line 169 "dspcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "dspcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 171 "dspcon.f"
	*info = -1;
#line 172 "dspcon.f"
    } else if (*n < 0) {
#line 173 "dspcon.f"
	*info = -2;
#line 174 "dspcon.f"
    } else if (*anorm < 0.) {
#line 175 "dspcon.f"
	*info = -5;
#line 176 "dspcon.f"
    }
#line 177 "dspcon.f"
    if (*info != 0) {
#line 178 "dspcon.f"
	i__1 = -(*info);
#line 178 "dspcon.f"
	xerbla_("DSPCON", &i__1, (ftnlen)6);
#line 179 "dspcon.f"
	return 0;
#line 180 "dspcon.f"
    }

/*     Quick return if possible */

#line 184 "dspcon.f"
    *rcond = 0.;
#line 185 "dspcon.f"
    if (*n == 0) {
#line 186 "dspcon.f"
	*rcond = 1.;
#line 187 "dspcon.f"
	return 0;
#line 188 "dspcon.f"
    } else if (*anorm <= 0.) {
#line 189 "dspcon.f"
	return 0;
#line 190 "dspcon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 194 "dspcon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 198 "dspcon.f"
	ip = *n * (*n + 1) / 2;
#line 199 "dspcon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 200 "dspcon.f"
	    if (ipiv[i__] > 0 && ap[ip] == 0.) {
#line 200 "dspcon.f"
		return 0;
#line 200 "dspcon.f"
	    }
#line 202 "dspcon.f"
	    ip -= i__;
#line 203 "dspcon.f"
/* L10: */
#line 203 "dspcon.f"
	}
#line 204 "dspcon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 208 "dspcon.f"
	ip = 1;
#line 209 "dspcon.f"
	i__1 = *n;
#line 209 "dspcon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 210 "dspcon.f"
	    if (ipiv[i__] > 0 && ap[ip] == 0.) {
#line 210 "dspcon.f"
		return 0;
#line 210 "dspcon.f"
	    }
#line 212 "dspcon.f"
	    ip = ip + *n - i__ + 1;
#line 213 "dspcon.f"
/* L20: */
#line 213 "dspcon.f"
	}
#line 214 "dspcon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 218 "dspcon.f"
    kase = 0;
#line 219 "dspcon.f"
L30:
#line 220 "dspcon.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 221 "dspcon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 225 "dspcon.f"
	dsptrs_(uplo, n, &c__1, &ap[1], &ipiv[1], &work[1], n, info, (ftnlen)
		1);
#line 226 "dspcon.f"
	goto L30;
#line 227 "dspcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 231 "dspcon.f"
    if (ainvnm != 0.) {
#line 231 "dspcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 231 "dspcon.f"
    }

#line 234 "dspcon.f"
    return 0;

/*     End of DSPCON */

} /* dspcon_ */

