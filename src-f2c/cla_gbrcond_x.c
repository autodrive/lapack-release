#line 1 "cla_gbrcond_x.f"
/* cla_gbrcond_x.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gbrcond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_GBRCOND_X computes the infinity norm condition number of op(A)*diag(x) for general banded m
atrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GBRCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbr
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbr
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbr
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_GBRCOND_X( TRANS, N, KL, KU, AB, LDAB, AFB, */
/*                                    LDAFB, IPIV, X, INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), */
/*      $                   X( * ) */
/*       REAL               RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLA_GBRCOND_X Computes the infinity norm condition number of */
/* >    op(A) * diag(X) where X is a COMPLEX vector. */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          AFB is COMPLEX array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by CGBTRF.  U is stored as an upper triangular */
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
/* >     as computed by CGBTRF; row i of the matrix was interchanged */
/* >     with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (N) */
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
/* >          WORK is COMPLEX array, dimension (2*N). */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N). */
/* >     Workspace. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
doublereal cla_gbrcond_x__(char *trans, integer *n, integer *kl, integer *ku, 
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
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen), cgbtrs_(char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static doublereal ainvnm;
    static logical notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

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

#line 201 "cla_gbrcond_x.f"
    /* Parameter adjustments */
#line 201 "cla_gbrcond_x.f"
    ab_dim1 = *ldab;
#line 201 "cla_gbrcond_x.f"
    ab_offset = 1 + ab_dim1;
#line 201 "cla_gbrcond_x.f"
    ab -= ab_offset;
#line 201 "cla_gbrcond_x.f"
    afb_dim1 = *ldafb;
#line 201 "cla_gbrcond_x.f"
    afb_offset = 1 + afb_dim1;
#line 201 "cla_gbrcond_x.f"
    afb -= afb_offset;
#line 201 "cla_gbrcond_x.f"
    --ipiv;
#line 201 "cla_gbrcond_x.f"
    --x;
#line 201 "cla_gbrcond_x.f"
    --work;
#line 201 "cla_gbrcond_x.f"
    --rwork;
#line 201 "cla_gbrcond_x.f"

#line 201 "cla_gbrcond_x.f"
    /* Function Body */
#line 201 "cla_gbrcond_x.f"
    ret_val = 0.;

#line 203 "cla_gbrcond_x.f"
    *info = 0;
#line 204 "cla_gbrcond_x.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 205 "cla_gbrcond_x.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 207 "cla_gbrcond_x.f"
	*info = -1;
#line 208 "cla_gbrcond_x.f"
    } else if (*n < 0) {
#line 209 "cla_gbrcond_x.f"
	*info = -2;
#line 210 "cla_gbrcond_x.f"
    } else if (*kl < 0 || *kl > *n - 1) {
#line 211 "cla_gbrcond_x.f"
	*info = -3;
#line 212 "cla_gbrcond_x.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 213 "cla_gbrcond_x.f"
	*info = -4;
#line 214 "cla_gbrcond_x.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 215 "cla_gbrcond_x.f"
	*info = -6;
#line 216 "cla_gbrcond_x.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 217 "cla_gbrcond_x.f"
	*info = -8;
#line 218 "cla_gbrcond_x.f"
    }
#line 219 "cla_gbrcond_x.f"
    if (*info != 0) {
#line 220 "cla_gbrcond_x.f"
	i__1 = -(*info);
#line 220 "cla_gbrcond_x.f"
	xerbla_("CLA_GBRCOND_X", &i__1, (ftnlen)13);
#line 221 "cla_gbrcond_x.f"
	return ret_val;
#line 222 "cla_gbrcond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 226 "cla_gbrcond_x.f"
    kd = *ku + 1;
#line 227 "cla_gbrcond_x.f"
    ke = *kl + 1;
#line 228 "cla_gbrcond_x.f"
    anorm = 0.;
#line 229 "cla_gbrcond_x.f"
    if (notrans) {
#line 230 "cla_gbrcond_x.f"
	i__1 = *n;
#line 230 "cla_gbrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 231 "cla_gbrcond_x.f"
	    tmp = 0.;
/* Computing MAX */
#line 232 "cla_gbrcond_x.f"
	    i__2 = i__ - *kl;
/* Computing MIN */
#line 232 "cla_gbrcond_x.f"
	    i__4 = i__ + *ku;
#line 232 "cla_gbrcond_x.f"
	    i__3 = min(i__4,*n);
#line 232 "cla_gbrcond_x.f"
	    for (j = max(i__2,1); j <= i__3; ++j) {
#line 233 "cla_gbrcond_x.f"
		i__2 = kd + i__ - j + j * ab_dim1;
#line 233 "cla_gbrcond_x.f"
		i__4 = j;
#line 233 "cla_gbrcond_x.f"
		z__2.r = ab[i__2].r * x[i__4].r - ab[i__2].i * x[i__4].i, 
			z__2.i = ab[i__2].r * x[i__4].i + ab[i__2].i * x[i__4]
			.r;
#line 233 "cla_gbrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 233 "cla_gbrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 234 "cla_gbrcond_x.f"
	    }
#line 235 "cla_gbrcond_x.f"
	    rwork[i__] = tmp;
#line 236 "cla_gbrcond_x.f"
	    anorm = max(anorm,tmp);
#line 237 "cla_gbrcond_x.f"
	}
#line 238 "cla_gbrcond_x.f"
    } else {
#line 239 "cla_gbrcond_x.f"
	i__1 = *n;
#line 239 "cla_gbrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "cla_gbrcond_x.f"
	    tmp = 0.;
/* Computing MAX */
#line 241 "cla_gbrcond_x.f"
	    i__3 = i__ - *kl;
/* Computing MIN */
#line 241 "cla_gbrcond_x.f"
	    i__4 = i__ + *ku;
#line 241 "cla_gbrcond_x.f"
	    i__2 = min(i__4,*n);
#line 241 "cla_gbrcond_x.f"
	    for (j = max(i__3,1); j <= i__2; ++j) {
#line 242 "cla_gbrcond_x.f"
		i__3 = ke - i__ + j + i__ * ab_dim1;
#line 242 "cla_gbrcond_x.f"
		i__4 = j;
#line 242 "cla_gbrcond_x.f"
		z__2.r = ab[i__3].r * x[i__4].r - ab[i__3].i * x[i__4].i, 
			z__2.i = ab[i__3].r * x[i__4].i + ab[i__3].i * x[i__4]
			.r;
#line 242 "cla_gbrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 242 "cla_gbrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 243 "cla_gbrcond_x.f"
	    }
#line 244 "cla_gbrcond_x.f"
	    rwork[i__] = tmp;
#line 245 "cla_gbrcond_x.f"
	    anorm = max(anorm,tmp);
#line 246 "cla_gbrcond_x.f"
	}
#line 247 "cla_gbrcond_x.f"
    }

/*     Quick return if possible. */

#line 251 "cla_gbrcond_x.f"
    if (*n == 0) {
#line 252 "cla_gbrcond_x.f"
	ret_val = 1.;
#line 253 "cla_gbrcond_x.f"
	return ret_val;
#line 254 "cla_gbrcond_x.f"
    } else if (anorm == 0.) {
#line 255 "cla_gbrcond_x.f"
	return ret_val;
#line 256 "cla_gbrcond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 260 "cla_gbrcond_x.f"
    ainvnm = 0.;

#line 262 "cla_gbrcond_x.f"
    kase = 0;
#line 263 "cla_gbrcond_x.f"
L10:
#line 264 "cla_gbrcond_x.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 265 "cla_gbrcond_x.f"
    if (kase != 0) {
#line 266 "cla_gbrcond_x.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 270 "cla_gbrcond_x.f"
	    i__1 = *n;
#line 270 "cla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "cla_gbrcond_x.f"
		i__2 = i__;
#line 271 "cla_gbrcond_x.f"
		i__3 = i__;
#line 271 "cla_gbrcond_x.f"
		i__4 = i__;
#line 271 "cla_gbrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 271 "cla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 272 "cla_gbrcond_x.f"
	    }

#line 274 "cla_gbrcond_x.f"
	    if (notrans) {
#line 275 "cla_gbrcond_x.f"
		cgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 277 "cla_gbrcond_x.f"
	    } else {
#line 278 "cla_gbrcond_x.f"
		cgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info, (
			ftnlen)19);
#line 280 "cla_gbrcond_x.f"
	    }

/*           Multiply by inv(X). */

#line 284 "cla_gbrcond_x.f"
	    i__1 = *n;
#line 284 "cla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "cla_gbrcond_x.f"
		i__2 = i__;
#line 285 "cla_gbrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 285 "cla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 286 "cla_gbrcond_x.f"
	    }
#line 287 "cla_gbrcond_x.f"
	} else {

/*           Multiply by inv(X**H). */

#line 291 "cla_gbrcond_x.f"
	    i__1 = *n;
#line 291 "cla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 292 "cla_gbrcond_x.f"
		i__2 = i__;
#line 292 "cla_gbrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 292 "cla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 293 "cla_gbrcond_x.f"
	    }

#line 295 "cla_gbrcond_x.f"
	    if (notrans) {
#line 296 "cla_gbrcond_x.f"
		cgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info, (
			ftnlen)19);
#line 298 "cla_gbrcond_x.f"
	    } else {
#line 299 "cla_gbrcond_x.f"
		cgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 301 "cla_gbrcond_x.f"
	    }

/*           Multiply by R. */

#line 305 "cla_gbrcond_x.f"
	    i__1 = *n;
#line 305 "cla_gbrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "cla_gbrcond_x.f"
		i__2 = i__;
#line 306 "cla_gbrcond_x.f"
		i__3 = i__;
#line 306 "cla_gbrcond_x.f"
		i__4 = i__;
#line 306 "cla_gbrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 306 "cla_gbrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 307 "cla_gbrcond_x.f"
	    }
#line 308 "cla_gbrcond_x.f"
	}
#line 309 "cla_gbrcond_x.f"
	goto L10;
#line 310 "cla_gbrcond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 314 "cla_gbrcond_x.f"
    if (ainvnm != 0.) {
#line 314 "cla_gbrcond_x.f"
	ret_val = 1. / ainvnm;
#line 314 "cla_gbrcond_x.f"
    }

#line 317 "cla_gbrcond_x.f"
    return ret_val;

} /* cla_gbrcond_x__ */

