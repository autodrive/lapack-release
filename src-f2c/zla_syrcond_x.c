#line 1 "zla_syrcond_x.f"
/* zla_syrcond_x.f -- translated by f2c (version 20100827).
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

#line 1 "zla_syrcond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_SYRCOND_X computes the infinity norm condition number of op(A)*diag(x) for symmetric indefi
nite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_SYRCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_syr
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_syr
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_syr
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_SYRCOND_X( UPLO, N, A, LDA, AF, */
/*                                                LDAF, IPIV, X, INFO, */
/*                                                WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
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
/* >    ZLA_SYRCOND_X Computes the infinity norm condition number of */
/* >    op(A) * diag(X) where X is a COMPLEX*16 vector. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >       = 'U':  Upper triangle of A is stored; */
/* >       = 'L':  Lower triangle of A is stored. */
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
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by ZSYTRF. */
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
/* >     Details of the interchanges and the block structure of D */
/* >     as determined by ZSYTRF. */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
doublereal zla_syrcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *
	x, integer *info, doublecomplex *work, doublereal *rwork, ftnlen 
	uplo_len)
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
    static logical up;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    static logical upper;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int zsytrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


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

#line 181 "zla_syrcond_x.f"
    /* Parameter adjustments */
#line 181 "zla_syrcond_x.f"
    a_dim1 = *lda;
#line 181 "zla_syrcond_x.f"
    a_offset = 1 + a_dim1;
#line 181 "zla_syrcond_x.f"
    a -= a_offset;
#line 181 "zla_syrcond_x.f"
    af_dim1 = *ldaf;
#line 181 "zla_syrcond_x.f"
    af_offset = 1 + af_dim1;
#line 181 "zla_syrcond_x.f"
    af -= af_offset;
#line 181 "zla_syrcond_x.f"
    --ipiv;
#line 181 "zla_syrcond_x.f"
    --x;
#line 181 "zla_syrcond_x.f"
    --work;
#line 181 "zla_syrcond_x.f"
    --rwork;
#line 181 "zla_syrcond_x.f"

#line 181 "zla_syrcond_x.f"
    /* Function Body */
#line 181 "zla_syrcond_x.f"
    ret_val = 0.;

#line 183 "zla_syrcond_x.f"
    *info = 0;
#line 184 "zla_syrcond_x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 185 "zla_syrcond_x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 186 "zla_syrcond_x.f"
	*info = -1;
#line 187 "zla_syrcond_x.f"
    } else if (*n < 0) {
#line 188 "zla_syrcond_x.f"
	*info = -2;
#line 189 "zla_syrcond_x.f"
    } else if (*lda < max(1,*n)) {
#line 190 "zla_syrcond_x.f"
	*info = -4;
#line 191 "zla_syrcond_x.f"
    } else if (*ldaf < max(1,*n)) {
#line 192 "zla_syrcond_x.f"
	*info = -6;
#line 193 "zla_syrcond_x.f"
    }
#line 194 "zla_syrcond_x.f"
    if (*info != 0) {
#line 195 "zla_syrcond_x.f"
	i__1 = -(*info);
#line 195 "zla_syrcond_x.f"
	xerbla_("ZLA_SYRCOND_X", &i__1, (ftnlen)13);
#line 196 "zla_syrcond_x.f"
	return ret_val;
#line 197 "zla_syrcond_x.f"
    }
#line 198 "zla_syrcond_x.f"
    up = FALSE_;
#line 199 "zla_syrcond_x.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 199 "zla_syrcond_x.f"
	up = TRUE_;
#line 199 "zla_syrcond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 203 "zla_syrcond_x.f"
    anorm = 0.;
#line 204 "zla_syrcond_x.f"
    if (up) {
#line 205 "zla_syrcond_x.f"
	i__1 = *n;
#line 205 "zla_syrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "zla_syrcond_x.f"
	    tmp = 0.;
#line 207 "zla_syrcond_x.f"
	    i__2 = i__;
#line 207 "zla_syrcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 208 "zla_syrcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 208 "zla_syrcond_x.f"
		i__4 = j;
#line 208 "zla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 208 "zla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 208 "zla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 209 "zla_syrcond_x.f"
	    }
#line 210 "zla_syrcond_x.f"
	    i__2 = *n;
#line 210 "zla_syrcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 211 "zla_syrcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 211 "zla_syrcond_x.f"
		i__4 = j;
#line 211 "zla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 211 "zla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 211 "zla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 212 "zla_syrcond_x.f"
	    }
#line 213 "zla_syrcond_x.f"
	    rwork[i__] = tmp;
#line 214 "zla_syrcond_x.f"
	    anorm = max(anorm,tmp);
#line 215 "zla_syrcond_x.f"
	}
#line 216 "zla_syrcond_x.f"
    } else {
#line 217 "zla_syrcond_x.f"
	i__1 = *n;
#line 217 "zla_syrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "zla_syrcond_x.f"
	    tmp = 0.;
#line 219 "zla_syrcond_x.f"
	    i__2 = i__;
#line 219 "zla_syrcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 220 "zla_syrcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 220 "zla_syrcond_x.f"
		i__4 = j;
#line 220 "zla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 220 "zla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 220 "zla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 221 "zla_syrcond_x.f"
	    }
#line 222 "zla_syrcond_x.f"
	    i__2 = *n;
#line 222 "zla_syrcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 223 "zla_syrcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 223 "zla_syrcond_x.f"
		i__4 = j;
#line 223 "zla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 223 "zla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 223 "zla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 224 "zla_syrcond_x.f"
	    }
#line 225 "zla_syrcond_x.f"
	    rwork[i__] = tmp;
#line 226 "zla_syrcond_x.f"
	    anorm = max(anorm,tmp);
#line 227 "zla_syrcond_x.f"
	}
#line 228 "zla_syrcond_x.f"
    }

/*     Quick return if possible. */

#line 232 "zla_syrcond_x.f"
    if (*n == 0) {
#line 233 "zla_syrcond_x.f"
	ret_val = 1.;
#line 234 "zla_syrcond_x.f"
	return ret_val;
#line 235 "zla_syrcond_x.f"
    } else if (anorm == 0.) {
#line 236 "zla_syrcond_x.f"
	return ret_val;
#line 237 "zla_syrcond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 241 "zla_syrcond_x.f"
    ainvnm = 0.;

#line 243 "zla_syrcond_x.f"
    kase = 0;
#line 244 "zla_syrcond_x.f"
L10:
#line 245 "zla_syrcond_x.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 246 "zla_syrcond_x.f"
    if (kase != 0) {
#line 247 "zla_syrcond_x.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 251 "zla_syrcond_x.f"
	    i__1 = *n;
#line 251 "zla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 252 "zla_syrcond_x.f"
		i__2 = i__;
#line 252 "zla_syrcond_x.f"
		i__3 = i__;
#line 252 "zla_syrcond_x.f"
		i__4 = i__;
#line 252 "zla_syrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 252 "zla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 253 "zla_syrcond_x.f"
	    }

#line 255 "zla_syrcond_x.f"
	    if (up) {
#line 256 "zla_syrcond_x.f"
		zsytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 258 "zla_syrcond_x.f"
	    } else {
#line 259 "zla_syrcond_x.f"
		zsytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 261 "zla_syrcond_x.f"
	    }

/*           Multiply by inv(X). */

#line 265 "zla_syrcond_x.f"
	    i__1 = *n;
#line 265 "zla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "zla_syrcond_x.f"
		i__2 = i__;
#line 266 "zla_syrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 266 "zla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 267 "zla_syrcond_x.f"
	    }
#line 268 "zla_syrcond_x.f"
	} else {

/*           Multiply by inv(X**T). */

#line 272 "zla_syrcond_x.f"
	    i__1 = *n;
#line 272 "zla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "zla_syrcond_x.f"
		i__2 = i__;
#line 273 "zla_syrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 273 "zla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 274 "zla_syrcond_x.f"
	    }

#line 276 "zla_syrcond_x.f"
	    if (up) {
#line 277 "zla_syrcond_x.f"
		zsytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 279 "zla_syrcond_x.f"
	    } else {
#line 280 "zla_syrcond_x.f"
		zsytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 282 "zla_syrcond_x.f"
	    }

/*           Multiply by R. */

#line 286 "zla_syrcond_x.f"
	    i__1 = *n;
#line 286 "zla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "zla_syrcond_x.f"
		i__2 = i__;
#line 287 "zla_syrcond_x.f"
		i__3 = i__;
#line 287 "zla_syrcond_x.f"
		i__4 = i__;
#line 287 "zla_syrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 287 "zla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 288 "zla_syrcond_x.f"
	    }
#line 289 "zla_syrcond_x.f"
	}
#line 290 "zla_syrcond_x.f"
	goto L10;
#line 291 "zla_syrcond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 295 "zla_syrcond_x.f"
    if (ainvnm != 0.) {
#line 295 "zla_syrcond_x.f"
	ret_val = 1. / ainvnm;
#line 295 "zla_syrcond_x.f"
    }

#line 298 "zla_syrcond_x.f"
    return ret_val;

} /* zla_syrcond_x__ */

