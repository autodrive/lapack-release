#line 1 "zla_gbrcond_x.f"
/* zla_gbrcond_x.f -- translated by f2c (version 20100827).
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

#line 1 "zla_gbrcond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_GBRCOND_X computes the infinity norm condition number of op(A)*diag(x) for general banded m
atrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_GBRCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gbr
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gbr
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gbr
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_GBRCOND_X( TRANS, N, KL, KU, AB, */
/*                                                LDAB, AFB, LDAFB, IPIV, */
/*                                                X, INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), */
/*      $                   X( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */



/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLA_GBRCOND_X Computes the infinity norm condition number of */
/* >    op(A) * diag(X) where X is a COMPLEX*16 vector. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >     Specifies the form of the system of equations: */
/* >       = 'N':  A * X = B     (No transpose) */
/* >       = 'T':  A**T * X = B  (Transpose) */
/* >       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >     The number of linear equations, i.e., the order of the */
/* >     matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >     The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >     The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >     On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >     The j-th column of A is stored in the j-th column of the */
/* >     array AB as follows: */
/* >     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl) */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >     The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* >          AFB is COMPLEX*16 array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by ZGBTRF.  U is stored as an upper triangular */
/* >     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* >     and the multipliers used during the factorization are stored */
/* >     in rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >     The pivot indices from the factorization A = P*L*U */
/* >     as computed by ZGBTRF; row i of the matrix was interchanged */
/* >     with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (N) */
/* >     The vector X in the formula op(A) * diag(X). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >       = 0:  Successful exit. */
/* >     i > 0:  The ith argument is invalid. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (2*N). */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N). */
/* >     Workspace. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16GBcomputational */

/*  ===================================================================== */
doublereal zla_gbrcond_x__(char *trans, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, 
	integer *ipiv, doublecomplex *x, integer *info, doublecomplex *work, 
	doublereal *rwork, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, kd, ke;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int zgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static logical notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */


/*  ===================================================================== */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 203 "zla_gbrcond_x.f"
    /* Parameter adjustments */
#line 203 "zla_gbrcond_x.f"
    ab_dim1 = *ldab;
#line 203 "zla_gbrcond_x.f"
    ab_offset = 1 + ab_dim1;
#line 203 "zla_gbrcond_x.f"
    ab -= ab_offset;
#line 203 "zla_gbrcond_x.f"
    afb_dim1 = *ldafb;
#line 203 "zla_gbrcond_x.f"
    afb_offset = 1 + afb_dim1;
#line 203 "zla_gbrcond_x.f"
    afb -= afb_offset;
#line 203 "zla_gbrcond_x.f"
    --ipiv;
#line 203 "zla_gbrcond_x.f"
    --x;
#line 203 "zla_gbrcond_x.f"
    --work;
#line 203 "zla_gbrcond_x.f"
    --rwork;
#line 203 "zla_gbrcond_x.f"

#line 203 "zla_gbrcond_x.f"
    /* Function Body */
#line 203 "zla_gbrcond_x.f"
    ret_val = 0.;

#line 205 "zla_gbrcond_x.f"
    *info = 0;
#line 206 "zla_gbrcond_x.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 207 "zla_gbrcond_x.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 209 "zla_gbrcond_x.f"
	*info = -1;
#line 210 "zla_gbrcond_x.f"
    } else if (*n < 0) {
#line 211 "zla_gbrcond_x.f"
	*info = -2;
#line 212 "zla_gbrcond_x.f"
    } else if (*kl < 0 || *kl > *n - 1) {
#line 213 "zla_gbrcond_x.f"
	*info = -3;
#line 214 "zla_gbrcond_x.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 215 "zla_gbrcond_x.f"
	*info = -4;
#line 216 "zla_gbrcond_x.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 217 "zla_gbrcond_x.f"
	*info = -6;
#line 218 "zla_gbrcond_x.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 219 "zla_gbrcond_x.f"
	*info = -8;
#line 220 "zla_gbrcond_x.f"
    }
#line 221 "zla_gbrcond_x.f"
    if (*info != 0) {
#line 222 "zla_gbrcond_x.f"
	i__1 = -(*info);
#line 222 "zla_gbrcond_x.f"
	xerbla_("ZLA_GBRCOND_X", &i__1, (ftnlen)13);
#line 223 "zla_gbrcond_x.f"
	return ret_val;
#line 224 "zla_gbrcond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 228 "zla_gbrcond_x.f"
    kd = *ku + 1;
#line 229 "zla_gbrcond_x.f"
    ke = *kl + 1;
#line 230 "zla_gbrcond_x.f"
    anorm = 0.;
#line 231 "zla_gbrcond_x.f"
    if (notrans) {
#line 232 "zla_gbrcond_x.f"
	i__1 = *n;
#line 232 "zla_gbrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "zla_gbrcond_x.f"
	    tmp = 0.;
/* Computing MAX */
#line 234 "zla_gbrcond_x.f"
	    i__2 = i__ - *kl;
/* Computing MIN */
#line 234 "zla_gbrcond_x.f"
	    i__4 = i__ + *ku;
#line 234 "zla_gbrcond_x.f"
	    i__3 = min(i__4,*n);
#line 234 "zla_gbrcond_x.f"
	    for (j = max(i__2,1); j <= i__3; ++j) {
#line 235 "zla_gbrcond_x.f"
		i__2 = kd + i__ - j + j * ab_dim1;
#line 235 "zla_gbrcond_x.f"
		i__4 = j;
#line 235 "zla_gbrcond_x.f"
		z__2.r = ab[i__2].r * x[i__4].r - ab[i__2].i * x[i__4].i, 
			z__2.i = ab[i__2].r * x[i__4].i + ab[i__2].i * x[i__4]
			.r;
#line 235 "zla_gbrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 235 "zla_gbrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 236 "zla_gbrcond_x.f"
	    }
#line 237 "zla_gbrcond_x.f"
	    rwork[i__] = tmp;
#line 238 "zla_gbrcond_x.f"
	    anorm = max(anorm,tmp);
#line 239 "zla_gbrcond_x.f"
	}
#line 240 "zla_gbrcond_x.f"
    } else {
#line 241 "zla_gbrcond_x.f"
	i__1 = *n;
#line 241 "zla_gbrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "zla_gbrcond_x.f"
	    tmp = 0.;
/* Computing MAX */
#line 243 "zla_gbrcond_x.f"
	    i__3 = i__ - *kl;
/* Computing MIN */
#line 243 "zla_gbrcond_x.f"
	    i__4 = i__ + *ku;
#line 243 "zla_gbrcond_x.f"
	    i__2 = min(i__4,*n);
#line 243 "zla_gbrcond_x.f"
	    for (j = max(i__3,1); j <= i__2; ++j) {
#line 244 "zla_gbrcond_x.f"
		i__3 = ke - i__ + j + i__ * ab_dim1;
#line 244 "zla_gbrcond_x.f"
		i__4 = j;
#line 244 "zla_gbrcond_x.f"
		z__2.r = ab[i__3].r * x[i__4].r - ab[i__3].i * x[i__4].i, 
			z__2.i = ab[i__3].r * x[i__4].i + ab[i__3].i * x[i__4]
			.r;
#line 244 "zla_gbrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 244 "zla_gbrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 245 "zla_gbrcond_x.f"
	    }
#line 246 "zla_gbrcond_x.f"
	    rwork[i__] = tmp;
#line 247 "zla_gbrcond_x.f"
	    anorm = max(anorm,tmp);
#line 248 "zla_gbrcond_x.f"
	}
#line 249 "zla_gbrcond_x.f"
    }

/*     Quick return if possible. */

#line 253 "zla_gbrcond_x.f"
    if (*n == 0) {
#line 254 "zla_gbrcond_x.f"
	ret_val = 1.;
#line 255 "zla_gbrcond_x.f"
	return ret_val;
#line 256 "zla_gbrcond_x.f"
    } else if (anorm == 0.) {
#line 257 "zla_gbrcond_x.f"
	return ret_val;
#line 258 "zla_gbrcond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 262 "zla_gbrcond_x.f"
    ainvnm = 0.;

#line 264 "zla_gbrcond_x.f"
    kase = 0;
#line 265 "zla_gbrcond_x.f"
L10:
#line 266 "zla_gbrcond_x.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 267 "zla_gbrcond_x.f"
    if (kase != 0) {
#line 268 "zla_gbrcond_x.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 272 "zla_gbrcond_x.f"
	    i__1 = *n;
#line 272 "zla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "zla_gbrcond_x.f"
		i__2 = i__;
#line 273 "zla_gbrcond_x.f"
		i__3 = i__;
#line 273 "zla_gbrcond_x.f"
		i__4 = i__;
#line 273 "zla_gbrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 273 "zla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 274 "zla_gbrcond_x.f"
	    }

#line 276 "zla_gbrcond_x.f"
	    if (notrans) {
#line 277 "zla_gbrcond_x.f"
		zgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 279 "zla_gbrcond_x.f"
	    } else {
#line 280 "zla_gbrcond_x.f"
		zgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info, (
			ftnlen)19);
#line 282 "zla_gbrcond_x.f"
	    }

/*           Multiply by inv(X). */

#line 286 "zla_gbrcond_x.f"
	    i__1 = *n;
#line 286 "zla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "zla_gbrcond_x.f"
		i__2 = i__;
#line 287 "zla_gbrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 287 "zla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 288 "zla_gbrcond_x.f"
	    }
#line 289 "zla_gbrcond_x.f"
	} else {

/*           Multiply by inv(X**H). */

#line 293 "zla_gbrcond_x.f"
	    i__1 = *n;
#line 293 "zla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 294 "zla_gbrcond_x.f"
		i__2 = i__;
#line 294 "zla_gbrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 294 "zla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 295 "zla_gbrcond_x.f"
	    }

#line 297 "zla_gbrcond_x.f"
	    if (notrans) {
#line 298 "zla_gbrcond_x.f"
		zgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info, (
			ftnlen)19);
#line 300 "zla_gbrcond_x.f"
	    } else {
#line 301 "zla_gbrcond_x.f"
		zgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 303 "zla_gbrcond_x.f"
	    }

/*           Multiply by R. */

#line 307 "zla_gbrcond_x.f"
	    i__1 = *n;
#line 307 "zla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "zla_gbrcond_x.f"
		i__2 = i__;
#line 308 "zla_gbrcond_x.f"
		i__3 = i__;
#line 308 "zla_gbrcond_x.f"
		i__4 = i__;
#line 308 "zla_gbrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 308 "zla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 309 "zla_gbrcond_x.f"
	    }
#line 310 "zla_gbrcond_x.f"
	}
#line 311 "zla_gbrcond_x.f"
	goto L10;
#line 312 "zla_gbrcond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 316 "zla_gbrcond_x.f"
    if (ainvnm != 0.) {
#line 316 "zla_gbrcond_x.f"
	ret_val = 1. / ainvnm;
#line 316 "zla_gbrcond_x.f"
    }

#line 319 "zla_gbrcond_x.f"
    return ret_val;

} /* zla_gbrcond_x__ */

