#line 1 "sgtcon.f"
/* sgtcon.f -- translated by f2c (version 20100827).
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

#line 1 "sgtcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SGTCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGTCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               D( * ), DL( * ), DU( * ), DU2( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTCON estimates the reciprocal of the condition number of a real */
/* > tridiagonal matrix A using the LU factorization as computed by */
/* > SGTTRF. */
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
/* >          DL is REAL array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A as computed by SGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is REAL array, dimension (N-2) */
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
/* >          ANORM is REAL */
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
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

/* > \date September 2012 */

/* > \ingroup realGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgtcon_(char *norm, integer *n, doublereal *dl, 
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
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    static logical onenrm;
    extern /* Subroutine */ int sgttrs_(char *, integer *, integer *, 
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

#line 189 "sgtcon.f"
    /* Parameter adjustments */
#line 189 "sgtcon.f"
    --iwork;
#line 189 "sgtcon.f"
    --work;
#line 189 "sgtcon.f"
    --ipiv;
#line 189 "sgtcon.f"
    --du2;
#line 189 "sgtcon.f"
    --du;
#line 189 "sgtcon.f"
    --d__;
#line 189 "sgtcon.f"
    --dl;
#line 189 "sgtcon.f"

#line 189 "sgtcon.f"
    /* Function Body */
#line 189 "sgtcon.f"
    *info = 0;
#line 190 "sgtcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 191 "sgtcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 192 "sgtcon.f"
	*info = -1;
#line 193 "sgtcon.f"
    } else if (*n < 0) {
#line 194 "sgtcon.f"
	*info = -2;
#line 195 "sgtcon.f"
    } else if (*anorm < 0.) {
#line 196 "sgtcon.f"
	*info = -8;
#line 197 "sgtcon.f"
    }
#line 198 "sgtcon.f"
    if (*info != 0) {
#line 199 "sgtcon.f"
	i__1 = -(*info);
#line 199 "sgtcon.f"
	xerbla_("SGTCON", &i__1, (ftnlen)6);
#line 200 "sgtcon.f"
	return 0;
#line 201 "sgtcon.f"
    }

/*     Quick return if possible */

#line 205 "sgtcon.f"
    *rcond = 0.;
#line 206 "sgtcon.f"
    if (*n == 0) {
#line 207 "sgtcon.f"
	*rcond = 1.;
#line 208 "sgtcon.f"
	return 0;
#line 209 "sgtcon.f"
    } else if (*anorm == 0.) {
#line 210 "sgtcon.f"
	return 0;
#line 211 "sgtcon.f"
    }

/*     Check that D(1:N) is non-zero. */

#line 215 "sgtcon.f"
    i__1 = *n;
#line 215 "sgtcon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "sgtcon.f"
	if (d__[i__] == 0.) {
#line 216 "sgtcon.f"
	    return 0;
#line 216 "sgtcon.f"
	}
#line 218 "sgtcon.f"
/* L10: */
#line 218 "sgtcon.f"
    }

#line 220 "sgtcon.f"
    ainvnm = 0.;
#line 221 "sgtcon.f"
    if (onenrm) {
#line 222 "sgtcon.f"
	kase1 = 1;
#line 223 "sgtcon.f"
    } else {
#line 224 "sgtcon.f"
	kase1 = 2;
#line 225 "sgtcon.f"
    }
#line 226 "sgtcon.f"
    kase = 0;
#line 227 "sgtcon.f"
L20:
#line 228 "sgtcon.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 229 "sgtcon.f"
    if (kase != 0) {
#line 230 "sgtcon.f"
	if (kase == kase1) {

/*           Multiply by inv(U)*inv(L). */

#line 234 "sgtcon.f"
	    sgttrs_("No transpose", n, &c__1, &dl[1], &d__[1], &du[1], &du2[1]
		    , &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 236 "sgtcon.f"
	} else {

/*           Multiply by inv(L**T)*inv(U**T). */

#line 240 "sgtcon.f"
	    sgttrs_("Transpose", n, &c__1, &dl[1], &d__[1], &du[1], &du2[1], &
		    ipiv[1], &work[1], n, info, (ftnlen)9);
#line 242 "sgtcon.f"
	}
#line 243 "sgtcon.f"
	goto L20;
#line 244 "sgtcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 248 "sgtcon.f"
    if (ainvnm != 0.) {
#line 248 "sgtcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 248 "sgtcon.f"
    }

#line 251 "sgtcon.f"
    return 0;

/*     End of SGTCON */

} /* sgtcon_ */

