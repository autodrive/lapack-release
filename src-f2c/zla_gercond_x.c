#line 1 "zla_gercond_x.f"
/* zla_gercond_x.f -- translated by f2c (version 20100827).
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

#line 1 "zla_gercond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_GERCOND_X computes the infinity norm condition number of op(A)*diag(x) for general matrices
. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_GERCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_ger
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_ger
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_ger
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_GERCOND_X( TRANS, N, A, LDA, AF, */
/*                                                LDAF, IPIV, X, INFO, */
/*                                                WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLA_GERCOND_X computes the infinity norm condition number of */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >     The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by ZGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >     The pivot indices from the factorization A = P*L*U */
/* >     as computed by ZGETRF; row i of the matrix was interchanged */
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

/* > \date December 2016 */

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
doublereal zla_gercond_x__(char *trans, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *
	x, integer *info, doublecomplex *work, doublereal *rwork, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int zgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    static logical notrans;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 184 "zla_gercond_x.f"
    /* Parameter adjustments */
#line 184 "zla_gercond_x.f"
    a_dim1 = *lda;
#line 184 "zla_gercond_x.f"
    a_offset = 1 + a_dim1;
#line 184 "zla_gercond_x.f"
    a -= a_offset;
#line 184 "zla_gercond_x.f"
    af_dim1 = *ldaf;
#line 184 "zla_gercond_x.f"
    af_offset = 1 + af_dim1;
#line 184 "zla_gercond_x.f"
    af -= af_offset;
#line 184 "zla_gercond_x.f"
    --ipiv;
#line 184 "zla_gercond_x.f"
    --x;
#line 184 "zla_gercond_x.f"
    --work;
#line 184 "zla_gercond_x.f"
    --rwork;
#line 184 "zla_gercond_x.f"

#line 184 "zla_gercond_x.f"
    /* Function Body */
#line 184 "zla_gercond_x.f"
    ret_val = 0.;

#line 186 "zla_gercond_x.f"
    *info = 0;
#line 187 "zla_gercond_x.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 188 "zla_gercond_x.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 190 "zla_gercond_x.f"
	*info = -1;
#line 191 "zla_gercond_x.f"
    } else if (*n < 0) {
#line 192 "zla_gercond_x.f"
	*info = -2;
#line 193 "zla_gercond_x.f"
    } else if (*lda < max(1,*n)) {
#line 194 "zla_gercond_x.f"
	*info = -4;
#line 195 "zla_gercond_x.f"
    } else if (*ldaf < max(1,*n)) {
#line 196 "zla_gercond_x.f"
	*info = -6;
#line 197 "zla_gercond_x.f"
    }
#line 198 "zla_gercond_x.f"
    if (*info != 0) {
#line 199 "zla_gercond_x.f"
	i__1 = -(*info);
#line 199 "zla_gercond_x.f"
	xerbla_("ZLA_GERCOND_X", &i__1, (ftnlen)13);
#line 200 "zla_gercond_x.f"
	return ret_val;
#line 201 "zla_gercond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 205 "zla_gercond_x.f"
    anorm = 0.;
#line 206 "zla_gercond_x.f"
    if (notrans) {
#line 207 "zla_gercond_x.f"
	i__1 = *n;
#line 207 "zla_gercond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 208 "zla_gercond_x.f"
	    tmp = 0.;
#line 209 "zla_gercond_x.f"
	    i__2 = *n;
#line 209 "zla_gercond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 210 "zla_gercond_x.f"
		i__3 = i__ + j * a_dim1;
#line 210 "zla_gercond_x.f"
		i__4 = j;
#line 210 "zla_gercond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 210 "zla_gercond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 210 "zla_gercond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 211 "zla_gercond_x.f"
	    }
#line 212 "zla_gercond_x.f"
	    rwork[i__] = tmp;
#line 213 "zla_gercond_x.f"
	    anorm = max(anorm,tmp);
#line 214 "zla_gercond_x.f"
	}
#line 215 "zla_gercond_x.f"
    } else {
#line 216 "zla_gercond_x.f"
	i__1 = *n;
#line 216 "zla_gercond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 217 "zla_gercond_x.f"
	    tmp = 0.;
#line 218 "zla_gercond_x.f"
	    i__2 = *n;
#line 218 "zla_gercond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 219 "zla_gercond_x.f"
		i__3 = j + i__ * a_dim1;
#line 219 "zla_gercond_x.f"
		i__4 = j;
#line 219 "zla_gercond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 219 "zla_gercond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 219 "zla_gercond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 220 "zla_gercond_x.f"
	    }
#line 221 "zla_gercond_x.f"
	    rwork[i__] = tmp;
#line 222 "zla_gercond_x.f"
	    anorm = max(anorm,tmp);
#line 223 "zla_gercond_x.f"
	}
#line 224 "zla_gercond_x.f"
    }

/*     Quick return if possible. */

#line 228 "zla_gercond_x.f"
    if (*n == 0) {
#line 229 "zla_gercond_x.f"
	ret_val = 1.;
#line 230 "zla_gercond_x.f"
	return ret_val;
#line 231 "zla_gercond_x.f"
    } else if (anorm == 0.) {
#line 232 "zla_gercond_x.f"
	return ret_val;
#line 233 "zla_gercond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 237 "zla_gercond_x.f"
    ainvnm = 0.;

#line 239 "zla_gercond_x.f"
    kase = 0;
#line 240 "zla_gercond_x.f"
L10:
#line 241 "zla_gercond_x.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 242 "zla_gercond_x.f"
    if (kase != 0) {
#line 243 "zla_gercond_x.f"
	if (kase == 2) {
/*           Multiply by R. */
#line 245 "zla_gercond_x.f"
	    i__1 = *n;
#line 245 "zla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "zla_gercond_x.f"
		i__2 = i__;
#line 246 "zla_gercond_x.f"
		i__3 = i__;
#line 246 "zla_gercond_x.f"
		i__4 = i__;
#line 246 "zla_gercond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 246 "zla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 247 "zla_gercond_x.f"
	    }

#line 249 "zla_gercond_x.f"
	    if (notrans) {
#line 250 "zla_gercond_x.f"
		zgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 252 "zla_gercond_x.f"
	    } else {
#line 253 "zla_gercond_x.f"
		zgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 255 "zla_gercond_x.f"
	    }

/*           Multiply by inv(X). */

#line 259 "zla_gercond_x.f"
	    i__1 = *n;
#line 259 "zla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 260 "zla_gercond_x.f"
		i__2 = i__;
#line 260 "zla_gercond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 260 "zla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 261 "zla_gercond_x.f"
	    }
#line 262 "zla_gercond_x.f"
	} else {

/*           Multiply by inv(X**H). */

#line 266 "zla_gercond_x.f"
	    i__1 = *n;
#line 266 "zla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "zla_gercond_x.f"
		i__2 = i__;
#line 267 "zla_gercond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 267 "zla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 268 "zla_gercond_x.f"
	    }

#line 270 "zla_gercond_x.f"
	    if (notrans) {
#line 271 "zla_gercond_x.f"
		zgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 273 "zla_gercond_x.f"
	    } else {
#line 274 "zla_gercond_x.f"
		zgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 276 "zla_gercond_x.f"
	    }

/*           Multiply by R. */

#line 280 "zla_gercond_x.f"
	    i__1 = *n;
#line 280 "zla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 281 "zla_gercond_x.f"
		i__2 = i__;
#line 281 "zla_gercond_x.f"
		i__3 = i__;
#line 281 "zla_gercond_x.f"
		i__4 = i__;
#line 281 "zla_gercond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 281 "zla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 282 "zla_gercond_x.f"
	    }
#line 283 "zla_gercond_x.f"
	}
#line 284 "zla_gercond_x.f"
	goto L10;
#line 285 "zla_gercond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 289 "zla_gercond_x.f"
    if (ainvnm != 0.) {
#line 289 "zla_gercond_x.f"
	ret_val = 1. / ainvnm;
#line 289 "zla_gercond_x.f"
    }

#line 292 "zla_gercond_x.f"
    return ret_val;

} /* zla_gercond_x__ */

