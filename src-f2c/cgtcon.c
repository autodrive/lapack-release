#line 1 "cgtcon.f"
/* cgtcon.f -- translated by f2c (version 20100827).
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

#line 1 "cgtcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CGTCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGTCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            D( * ), DL( * ), DU( * ), DU2( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTCON estimates the reciprocal of the condition number of a complex */
/* > tridiagonal matrix A using the LU factorization as computed by */
/* > CGTTRF. */
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
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A as computed by CGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX array, dimension (N-2) */
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

/* > \date September 2012 */

/* > \ingroup complexGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgtcon_(char *norm, integer *n, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, 
	integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, kase, kase1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    static logical onenrm;
    extern /* Subroutine */ int cgttrs_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, doublecomplex *, integer *, integer *, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 187 "cgtcon.f"
    /* Parameter adjustments */
#line 187 "cgtcon.f"
    --work;
#line 187 "cgtcon.f"
    --ipiv;
#line 187 "cgtcon.f"
    --du2;
#line 187 "cgtcon.f"
    --du;
#line 187 "cgtcon.f"
    --d__;
#line 187 "cgtcon.f"
    --dl;
#line 187 "cgtcon.f"

#line 187 "cgtcon.f"
    /* Function Body */
#line 187 "cgtcon.f"
    *info = 0;
#line 188 "cgtcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 189 "cgtcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 190 "cgtcon.f"
	*info = -1;
#line 191 "cgtcon.f"
    } else if (*n < 0) {
#line 192 "cgtcon.f"
	*info = -2;
#line 193 "cgtcon.f"
    } else if (*anorm < 0.) {
#line 194 "cgtcon.f"
	*info = -8;
#line 195 "cgtcon.f"
    }
#line 196 "cgtcon.f"
    if (*info != 0) {
#line 197 "cgtcon.f"
	i__1 = -(*info);
#line 197 "cgtcon.f"
	xerbla_("CGTCON", &i__1, (ftnlen)6);
#line 198 "cgtcon.f"
	return 0;
#line 199 "cgtcon.f"
    }

/*     Quick return if possible */

#line 203 "cgtcon.f"
    *rcond = 0.;
#line 204 "cgtcon.f"
    if (*n == 0) {
#line 205 "cgtcon.f"
	*rcond = 1.;
#line 206 "cgtcon.f"
	return 0;
#line 207 "cgtcon.f"
    } else if (*anorm == 0.) {
#line 208 "cgtcon.f"
	return 0;
#line 209 "cgtcon.f"
    }

/*     Check that D(1:N) is non-zero. */

#line 213 "cgtcon.f"
    i__1 = *n;
#line 213 "cgtcon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "cgtcon.f"
	i__2 = i__;
#line 214 "cgtcon.f"
	if (d__[i__2].r == 0. && d__[i__2].i == 0.) {
#line 214 "cgtcon.f"
	    return 0;
#line 214 "cgtcon.f"
	}
#line 216 "cgtcon.f"
/* L10: */
#line 216 "cgtcon.f"
    }

#line 218 "cgtcon.f"
    ainvnm = 0.;
#line 219 "cgtcon.f"
    if (onenrm) {
#line 220 "cgtcon.f"
	kase1 = 1;
#line 221 "cgtcon.f"
    } else {
#line 222 "cgtcon.f"
	kase1 = 2;
#line 223 "cgtcon.f"
    }
#line 224 "cgtcon.f"
    kase = 0;
#line 225 "cgtcon.f"
L20:
#line 226 "cgtcon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 227 "cgtcon.f"
    if (kase != 0) {
#line 228 "cgtcon.f"
	if (kase == kase1) {

/*           Multiply by inv(U)*inv(L). */

#line 232 "cgtcon.f"
	    cgttrs_("No transpose", n, &c__1, &dl[1], &d__[1], &du[1], &du2[1]
		    , &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 234 "cgtcon.f"
	} else {

/*           Multiply by inv(L**H)*inv(U**H). */

#line 238 "cgtcon.f"
	    cgttrs_("Conjugate transpose", n, &c__1, &dl[1], &d__[1], &du[1], 
		    &du2[1], &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 240 "cgtcon.f"
	}
#line 241 "cgtcon.f"
	goto L20;
#line 242 "cgtcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 246 "cgtcon.f"
    if (ainvnm != 0.) {
#line 246 "cgtcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 246 "cgtcon.f"
    }

#line 249 "cgtcon.f"
    return 0;

/*     End of CGTCON */

} /* cgtcon_ */

