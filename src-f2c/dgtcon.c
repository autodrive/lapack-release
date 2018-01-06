#line 1 "dgtcon.f"
/* dgtcon.f -- translated by f2c (version 20100827).
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

#line 1 "dgtcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGTCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGTCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), DL( * ), DU( * ), DU2( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTCON estimates the reciprocal of the condition number of a real */
/* > tridiagonal matrix A using the LU factorization as computed by */
/* > DGTTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies whether the 1-norm condition number or the */
/* >          infinity-norm condition number is required: */
/* >          = '1' or 'O':  1-norm; */
/* >          = 'I':         Infinity-norm. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A as computed by DGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is DOUBLE PRECISION array, dimension (N-2) */
/* >          The (n-2) elements of the second superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
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

/* > \date September 2012 */

/* > \ingroup doubleGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgtcon_(char *norm, integer *n, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, kase, kase1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    static logical onenrm;
    extern /* Subroutine */ int dgttrs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

/*     Test the input arguments. */

#line 189 "dgtcon.f"
    /* Parameter adjustments */
#line 189 "dgtcon.f"
    --iwork;
#line 189 "dgtcon.f"
    --work;
#line 189 "dgtcon.f"
    --ipiv;
#line 189 "dgtcon.f"
    --du2;
#line 189 "dgtcon.f"
    --du;
#line 189 "dgtcon.f"
    --d__;
#line 189 "dgtcon.f"
    --dl;
#line 189 "dgtcon.f"

#line 189 "dgtcon.f"
    /* Function Body */
#line 189 "dgtcon.f"
    *info = 0;
#line 190 "dgtcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 191 "dgtcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 192 "dgtcon.f"
	*info = -1;
#line 193 "dgtcon.f"
    } else if (*n < 0) {
#line 194 "dgtcon.f"
	*info = -2;
#line 195 "dgtcon.f"
    } else if (*anorm < 0.) {
#line 196 "dgtcon.f"
	*info = -8;
#line 197 "dgtcon.f"
    }
#line 198 "dgtcon.f"
    if (*info != 0) {
#line 199 "dgtcon.f"
	i__1 = -(*info);
#line 199 "dgtcon.f"
	xerbla_("DGTCON", &i__1, (ftnlen)6);
#line 200 "dgtcon.f"
	return 0;
#line 201 "dgtcon.f"
    }

/*     Quick return if possible */

#line 205 "dgtcon.f"
    *rcond = 0.;
#line 206 "dgtcon.f"
    if (*n == 0) {
#line 207 "dgtcon.f"
	*rcond = 1.;
#line 208 "dgtcon.f"
	return 0;
#line 209 "dgtcon.f"
    } else if (*anorm == 0.) {
#line 210 "dgtcon.f"
	return 0;
#line 211 "dgtcon.f"
    }

/*     Check that D(1:N) is non-zero. */

#line 215 "dgtcon.f"
    i__1 = *n;
#line 215 "dgtcon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "dgtcon.f"
	if (d__[i__] == 0.) {
#line 216 "dgtcon.f"
	    return 0;
#line 216 "dgtcon.f"
	}
#line 218 "dgtcon.f"
/* L10: */
#line 218 "dgtcon.f"
    }

#line 220 "dgtcon.f"
    ainvnm = 0.;
#line 221 "dgtcon.f"
    if (onenrm) {
#line 222 "dgtcon.f"
	kase1 = 1;
#line 223 "dgtcon.f"
    } else {
#line 224 "dgtcon.f"
	kase1 = 2;
#line 225 "dgtcon.f"
    }
#line 226 "dgtcon.f"
    kase = 0;
#line 227 "dgtcon.f"
L20:
#line 228 "dgtcon.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 229 "dgtcon.f"
    if (kase != 0) {
#line 230 "dgtcon.f"
	if (kase == kase1) {

/*           Multiply by inv(U)*inv(L). */

#line 234 "dgtcon.f"
	    dgttrs_("No transpose", n, &c__1, &dl[1], &d__[1], &du[1], &du2[1]
		    , &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 236 "dgtcon.f"
	} else {

/*           Multiply by inv(L**T)*inv(U**T). */

#line 240 "dgtcon.f"
	    dgttrs_("Transpose", n, &c__1, &dl[1], &d__[1], &du[1], &du2[1], &
		    ipiv[1], &work[1], n, info, (ftnlen)9);
#line 242 "dgtcon.f"
	}
#line 243 "dgtcon.f"
	goto L20;
#line 244 "dgtcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 248 "dgtcon.f"
    if (ainvnm != 0.) {
#line 248 "dgtcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 248 "dgtcon.f"
    }

#line 251 "dgtcon.f"
    return 0;

/*     End of DGTCON */

} /* dgtcon_ */

