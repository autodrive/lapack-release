#line 1 "zla_porcond_x.f"
/* zla_porcond_x.f -- translated by f2c (version 20100827).
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

#line 1 "zla_porcond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_PORCOND_X computes the infinity norm condition number of op(A)*diag(x) for Hermitian positi
ve-definite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_PORCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_por
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_por
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_por
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_PORCOND_X( UPLO, N, A, LDA, AF, */
/*                                                LDAF, X, INFO, WORK, */
/*                                                RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLA_PORCOND_X Computes the infinity norm condition number of */
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
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**H*U or A = L*L**H, as computed by ZPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
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

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
doublereal zla_porcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, doublecomplex *x, integer *
	info, doublecomplex *work, doublereal *rwork, ftnlen uplo_len)
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
    extern /* Subroutine */ int zpotrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);


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

#line 171 "zla_porcond_x.f"
    /* Parameter adjustments */
#line 171 "zla_porcond_x.f"
    a_dim1 = *lda;
#line 171 "zla_porcond_x.f"
    a_offset = 1 + a_dim1;
#line 171 "zla_porcond_x.f"
    a -= a_offset;
#line 171 "zla_porcond_x.f"
    af_dim1 = *ldaf;
#line 171 "zla_porcond_x.f"
    af_offset = 1 + af_dim1;
#line 171 "zla_porcond_x.f"
    af -= af_offset;
#line 171 "zla_porcond_x.f"
    --x;
#line 171 "zla_porcond_x.f"
    --work;
#line 171 "zla_porcond_x.f"
    --rwork;
#line 171 "zla_porcond_x.f"

#line 171 "zla_porcond_x.f"
    /* Function Body */
#line 171 "zla_porcond_x.f"
    ret_val = 0.;

#line 173 "zla_porcond_x.f"
    *info = 0;
#line 174 "zla_porcond_x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 175 "zla_porcond_x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 176 "zla_porcond_x.f"
	*info = -1;
#line 177 "zla_porcond_x.f"
    } else if (*n < 0) {
#line 178 "zla_porcond_x.f"
	*info = -2;
#line 179 "zla_porcond_x.f"
    } else if (*lda < max(1,*n)) {
#line 180 "zla_porcond_x.f"
	*info = -4;
#line 181 "zla_porcond_x.f"
    } else if (*ldaf < max(1,*n)) {
#line 182 "zla_porcond_x.f"
	*info = -6;
#line 183 "zla_porcond_x.f"
    }
#line 184 "zla_porcond_x.f"
    if (*info != 0) {
#line 185 "zla_porcond_x.f"
	i__1 = -(*info);
#line 185 "zla_porcond_x.f"
	xerbla_("ZLA_PORCOND_X", &i__1, (ftnlen)13);
#line 186 "zla_porcond_x.f"
	return ret_val;
#line 187 "zla_porcond_x.f"
    }
#line 188 "zla_porcond_x.f"
    up = FALSE_;
#line 189 "zla_porcond_x.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "zla_porcond_x.f"
	up = TRUE_;
#line 189 "zla_porcond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 193 "zla_porcond_x.f"
    anorm = 0.;
#line 194 "zla_porcond_x.f"
    if (up) {
#line 195 "zla_porcond_x.f"
	i__1 = *n;
#line 195 "zla_porcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 196 "zla_porcond_x.f"
	    tmp = 0.;
#line 197 "zla_porcond_x.f"
	    i__2 = i__;
#line 197 "zla_porcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 198 "zla_porcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 198 "zla_porcond_x.f"
		i__4 = j;
#line 198 "zla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 198 "zla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 198 "zla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 199 "zla_porcond_x.f"
	    }
#line 200 "zla_porcond_x.f"
	    i__2 = *n;
#line 200 "zla_porcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 201 "zla_porcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 201 "zla_porcond_x.f"
		i__4 = j;
#line 201 "zla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 201 "zla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 201 "zla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 202 "zla_porcond_x.f"
	    }
#line 203 "zla_porcond_x.f"
	    rwork[i__] = tmp;
#line 204 "zla_porcond_x.f"
	    anorm = max(anorm,tmp);
#line 205 "zla_porcond_x.f"
	}
#line 206 "zla_porcond_x.f"
    } else {
#line 207 "zla_porcond_x.f"
	i__1 = *n;
#line 207 "zla_porcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 208 "zla_porcond_x.f"
	    tmp = 0.;
#line 209 "zla_porcond_x.f"
	    i__2 = i__;
#line 209 "zla_porcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 210 "zla_porcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 210 "zla_porcond_x.f"
		i__4 = j;
#line 210 "zla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 210 "zla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 210 "zla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 211 "zla_porcond_x.f"
	    }
#line 212 "zla_porcond_x.f"
	    i__2 = *n;
#line 212 "zla_porcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 213 "zla_porcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 213 "zla_porcond_x.f"
		i__4 = j;
#line 213 "zla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 213 "zla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 213 "zla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 214 "zla_porcond_x.f"
	    }
#line 215 "zla_porcond_x.f"
	    rwork[i__] = tmp;
#line 216 "zla_porcond_x.f"
	    anorm = max(anorm,tmp);
#line 217 "zla_porcond_x.f"
	}
#line 218 "zla_porcond_x.f"
    }

/*     Quick return if possible. */

#line 222 "zla_porcond_x.f"
    if (*n == 0) {
#line 223 "zla_porcond_x.f"
	ret_val = 1.;
#line 224 "zla_porcond_x.f"
	return ret_val;
#line 225 "zla_porcond_x.f"
    } else if (anorm == 0.) {
#line 226 "zla_porcond_x.f"
	return ret_val;
#line 227 "zla_porcond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 231 "zla_porcond_x.f"
    ainvnm = 0.;

#line 233 "zla_porcond_x.f"
    kase = 0;
#line 234 "zla_porcond_x.f"
L10:
#line 235 "zla_porcond_x.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 236 "zla_porcond_x.f"
    if (kase != 0) {
#line 237 "zla_porcond_x.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 241 "zla_porcond_x.f"
	    i__1 = *n;
#line 241 "zla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "zla_porcond_x.f"
		i__2 = i__;
#line 242 "zla_porcond_x.f"
		i__3 = i__;
#line 242 "zla_porcond_x.f"
		i__4 = i__;
#line 242 "zla_porcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 242 "zla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 243 "zla_porcond_x.f"
	    }

#line 245 "zla_porcond_x.f"
	    if (up) {
#line 246 "zla_porcond_x.f"
		zpotrs_("U", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 248 "zla_porcond_x.f"
	    } else {
#line 249 "zla_porcond_x.f"
		zpotrs_("L", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 251 "zla_porcond_x.f"
	    }

/*           Multiply by inv(X). */

#line 255 "zla_porcond_x.f"
	    i__1 = *n;
#line 255 "zla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 256 "zla_porcond_x.f"
		i__2 = i__;
#line 256 "zla_porcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 256 "zla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 257 "zla_porcond_x.f"
	    }
#line 258 "zla_porcond_x.f"
	} else {

/*           Multiply by inv(X**H). */

#line 262 "zla_porcond_x.f"
	    i__1 = *n;
#line 262 "zla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 263 "zla_porcond_x.f"
		i__2 = i__;
#line 263 "zla_porcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 263 "zla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 264 "zla_porcond_x.f"
	    }

#line 266 "zla_porcond_x.f"
	    if (up) {
#line 267 "zla_porcond_x.f"
		zpotrs_("U", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 269 "zla_porcond_x.f"
	    } else {
#line 270 "zla_porcond_x.f"
		zpotrs_("L", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 272 "zla_porcond_x.f"
	    }

/*           Multiply by R. */

#line 276 "zla_porcond_x.f"
	    i__1 = *n;
#line 276 "zla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 277 "zla_porcond_x.f"
		i__2 = i__;
#line 277 "zla_porcond_x.f"
		i__3 = i__;
#line 277 "zla_porcond_x.f"
		i__4 = i__;
#line 277 "zla_porcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 277 "zla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 278 "zla_porcond_x.f"
	    }
#line 279 "zla_porcond_x.f"
	}
#line 280 "zla_porcond_x.f"
	goto L10;
#line 281 "zla_porcond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 285 "zla_porcond_x.f"
    if (ainvnm != 0.) {
#line 285 "zla_porcond_x.f"
	ret_val = 1. / ainvnm;
#line 285 "zla_porcond_x.f"
    }

#line 288 "zla_porcond_x.f"
    return ret_val;

} /* zla_porcond_x__ */

