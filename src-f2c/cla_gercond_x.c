#line 1 "cla_gercond_x.f"
/* cla_gercond_x.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gercond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_GERCOND_X computes the infinity norm condition number of op(A)*diag(x) for general matrices
. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GERCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X, */
/*                                    INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * ) */
/*       REAL               RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* >    CLA_GERCOND_X computes the infinity norm condition number of */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by CGETRF. */
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
/* >     as computed by CGETRF; row i of the matrix was interchanged */
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

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
doublereal cla_gercond_x__(char *trans, integer *n, doublecomplex *a, integer 
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
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen), cgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    static doublereal ainvnm;
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

#line 183 "cla_gercond_x.f"
    /* Parameter adjustments */
#line 183 "cla_gercond_x.f"
    a_dim1 = *lda;
#line 183 "cla_gercond_x.f"
    a_offset = 1 + a_dim1;
#line 183 "cla_gercond_x.f"
    a -= a_offset;
#line 183 "cla_gercond_x.f"
    af_dim1 = *ldaf;
#line 183 "cla_gercond_x.f"
    af_offset = 1 + af_dim1;
#line 183 "cla_gercond_x.f"
    af -= af_offset;
#line 183 "cla_gercond_x.f"
    --ipiv;
#line 183 "cla_gercond_x.f"
    --x;
#line 183 "cla_gercond_x.f"
    --work;
#line 183 "cla_gercond_x.f"
    --rwork;
#line 183 "cla_gercond_x.f"

#line 183 "cla_gercond_x.f"
    /* Function Body */
#line 183 "cla_gercond_x.f"
    ret_val = 0.;

#line 185 "cla_gercond_x.f"
    *info = 0;
#line 186 "cla_gercond_x.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 187 "cla_gercond_x.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 189 "cla_gercond_x.f"
	*info = -1;
#line 190 "cla_gercond_x.f"
    } else if (*n < 0) {
#line 191 "cla_gercond_x.f"
	*info = -2;
#line 192 "cla_gercond_x.f"
    } else if (*lda < max(1,*n)) {
#line 193 "cla_gercond_x.f"
	*info = -4;
#line 194 "cla_gercond_x.f"
    } else if (*ldaf < max(1,*n)) {
#line 195 "cla_gercond_x.f"
	*info = -6;
#line 196 "cla_gercond_x.f"
    }
#line 197 "cla_gercond_x.f"
    if (*info != 0) {
#line 198 "cla_gercond_x.f"
	i__1 = -(*info);
#line 198 "cla_gercond_x.f"
	xerbla_("CLA_GERCOND_X", &i__1, (ftnlen)13);
#line 199 "cla_gercond_x.f"
	return ret_val;
#line 200 "cla_gercond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 204 "cla_gercond_x.f"
    anorm = 0.;
#line 205 "cla_gercond_x.f"
    if (notrans) {
#line 206 "cla_gercond_x.f"
	i__1 = *n;
#line 206 "cla_gercond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 207 "cla_gercond_x.f"
	    tmp = 0.;
#line 208 "cla_gercond_x.f"
	    i__2 = *n;
#line 208 "cla_gercond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 209 "cla_gercond_x.f"
		i__3 = i__ + j * a_dim1;
#line 209 "cla_gercond_x.f"
		i__4 = j;
#line 209 "cla_gercond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 209 "cla_gercond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 209 "cla_gercond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 210 "cla_gercond_x.f"
	    }
#line 211 "cla_gercond_x.f"
	    rwork[i__] = tmp;
#line 212 "cla_gercond_x.f"
	    anorm = max(anorm,tmp);
#line 213 "cla_gercond_x.f"
	}
#line 214 "cla_gercond_x.f"
    } else {
#line 215 "cla_gercond_x.f"
	i__1 = *n;
#line 215 "cla_gercond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "cla_gercond_x.f"
	    tmp = 0.;
#line 217 "cla_gercond_x.f"
	    i__2 = *n;
#line 217 "cla_gercond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 218 "cla_gercond_x.f"
		i__3 = j + i__ * a_dim1;
#line 218 "cla_gercond_x.f"
		i__4 = j;
#line 218 "cla_gercond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 218 "cla_gercond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 218 "cla_gercond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 219 "cla_gercond_x.f"
	    }
#line 220 "cla_gercond_x.f"
	    rwork[i__] = tmp;
#line 221 "cla_gercond_x.f"
	    anorm = max(anorm,tmp);
#line 222 "cla_gercond_x.f"
	}
#line 223 "cla_gercond_x.f"
    }

/*     Quick return if possible. */

#line 227 "cla_gercond_x.f"
    if (*n == 0) {
#line 228 "cla_gercond_x.f"
	ret_val = 1.;
#line 229 "cla_gercond_x.f"
	return ret_val;
#line 230 "cla_gercond_x.f"
    } else if (anorm == 0.) {
#line 231 "cla_gercond_x.f"
	return ret_val;
#line 232 "cla_gercond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 236 "cla_gercond_x.f"
    ainvnm = 0.;

#line 238 "cla_gercond_x.f"
    kase = 0;
#line 239 "cla_gercond_x.f"
L10:
#line 240 "cla_gercond_x.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 241 "cla_gercond_x.f"
    if (kase != 0) {
#line 242 "cla_gercond_x.f"
	if (kase == 2) {
/*           Multiply by R. */
#line 244 "cla_gercond_x.f"
	    i__1 = *n;
#line 244 "cla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 245 "cla_gercond_x.f"
		i__2 = i__;
#line 245 "cla_gercond_x.f"
		i__3 = i__;
#line 245 "cla_gercond_x.f"
		i__4 = i__;
#line 245 "cla_gercond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 245 "cla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 246 "cla_gercond_x.f"
	    }

#line 248 "cla_gercond_x.f"
	    if (notrans) {
#line 249 "cla_gercond_x.f"
		cgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 251 "cla_gercond_x.f"
	    } else {
#line 252 "cla_gercond_x.f"
		cgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 254 "cla_gercond_x.f"
	    }

/*           Multiply by inv(X). */

#line 258 "cla_gercond_x.f"
	    i__1 = *n;
#line 258 "cla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 259 "cla_gercond_x.f"
		i__2 = i__;
#line 259 "cla_gercond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 259 "cla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 260 "cla_gercond_x.f"
	    }
#line 261 "cla_gercond_x.f"
	} else {

/*           Multiply by inv(X**H). */

#line 265 "cla_gercond_x.f"
	    i__1 = *n;
#line 265 "cla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "cla_gercond_x.f"
		i__2 = i__;
#line 266 "cla_gercond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 266 "cla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 267 "cla_gercond_x.f"
	    }

#line 269 "cla_gercond_x.f"
	    if (notrans) {
#line 270 "cla_gercond_x.f"
		cgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 272 "cla_gercond_x.f"
	    } else {
#line 273 "cla_gercond_x.f"
		cgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 275 "cla_gercond_x.f"
	    }

/*           Multiply by R. */

#line 279 "cla_gercond_x.f"
	    i__1 = *n;
#line 279 "cla_gercond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "cla_gercond_x.f"
		i__2 = i__;
#line 280 "cla_gercond_x.f"
		i__3 = i__;
#line 280 "cla_gercond_x.f"
		i__4 = i__;
#line 280 "cla_gercond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 280 "cla_gercond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 281 "cla_gercond_x.f"
	    }
#line 282 "cla_gercond_x.f"
	}
#line 283 "cla_gercond_x.f"
	goto L10;
#line 284 "cla_gercond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 288 "cla_gercond_x.f"
    if (ainvnm != 0.) {
#line 288 "cla_gercond_x.f"
	ret_val = 1. / ainvnm;
#line 288 "cla_gercond_x.f"
    }

#line 291 "cla_gercond_x.f"
    return ret_val;

} /* cla_gercond_x__ */

