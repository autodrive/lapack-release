#line 1 "zspcon.f"
/* zspcon.f -- translated by f2c (version 20100827).
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

#line 1 "zspcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZSPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO ) */

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
/* > ZSPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex symmetric packed matrix A using the */
/* > factorization A = U*D*U**T or A = L*D*L**T computed by ZSPTRF. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSPTRF. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zspcon_(char *uplo, integer *n, doublecomplex *ap, 
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
    extern /* Subroutine */ int zsptrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 161 "zspcon.f"
    /* Parameter adjustments */
#line 161 "zspcon.f"
    --work;
#line 161 "zspcon.f"
    --ipiv;
#line 161 "zspcon.f"
    --ap;
#line 161 "zspcon.f"

#line 161 "zspcon.f"
    /* Function Body */
#line 161 "zspcon.f"
    *info = 0;
#line 162 "zspcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "zspcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "zspcon.f"
	*info = -1;
#line 165 "zspcon.f"
    } else if (*n < 0) {
#line 166 "zspcon.f"
	*info = -2;
#line 167 "zspcon.f"
    } else if (*anorm < 0.) {
#line 168 "zspcon.f"
	*info = -5;
#line 169 "zspcon.f"
    }
#line 170 "zspcon.f"
    if (*info != 0) {
#line 171 "zspcon.f"
	i__1 = -(*info);
#line 171 "zspcon.f"
	xerbla_("ZSPCON", &i__1, (ftnlen)6);
#line 172 "zspcon.f"
	return 0;
#line 173 "zspcon.f"
    }

/*     Quick return if possible */

#line 177 "zspcon.f"
    *rcond = 0.;
#line 178 "zspcon.f"
    if (*n == 0) {
#line 179 "zspcon.f"
	*rcond = 1.;
#line 180 "zspcon.f"
	return 0;
#line 181 "zspcon.f"
    } else if (*anorm <= 0.) {
#line 182 "zspcon.f"
	return 0;
#line 183 "zspcon.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 187 "zspcon.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 191 "zspcon.f"
	ip = *n * (*n + 1) / 2;
#line 192 "zspcon.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 193 "zspcon.f"
	    i__1 = ip;
#line 193 "zspcon.f"
	    if (ipiv[i__] > 0 && (ap[i__1].r == 0. && ap[i__1].i == 0.)) {
#line 193 "zspcon.f"
		return 0;
#line 193 "zspcon.f"
	    }
#line 195 "zspcon.f"
	    ip -= i__;
#line 196 "zspcon.f"
/* L10: */
#line 196 "zspcon.f"
	}
#line 197 "zspcon.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 201 "zspcon.f"
	ip = 1;
#line 202 "zspcon.f"
	i__1 = *n;
#line 202 "zspcon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "zspcon.f"
	    i__2 = ip;
#line 203 "zspcon.f"
	    if (ipiv[i__] > 0 && (ap[i__2].r == 0. && ap[i__2].i == 0.)) {
#line 203 "zspcon.f"
		return 0;
#line 203 "zspcon.f"
	    }
#line 205 "zspcon.f"
	    ip = ip + *n - i__ + 1;
#line 206 "zspcon.f"
/* L20: */
#line 206 "zspcon.f"
	}
#line 207 "zspcon.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 211 "zspcon.f"
    kase = 0;
#line 212 "zspcon.f"
L30:
#line 213 "zspcon.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 214 "zspcon.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 218 "zspcon.f"
	zsptrs_(uplo, n, &c__1, &ap[1], &ipiv[1], &work[1], n, info, (ftnlen)
		1);
#line 219 "zspcon.f"
	goto L30;
#line 220 "zspcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 224 "zspcon.f"
    if (ainvnm != 0.) {
#line 224 "zspcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 224 "zspcon.f"
    }

#line 227 "zspcon.f"
    return 0;

/*     End of ZSPCON */

} /* zspcon_ */

