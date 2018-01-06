#line 1 "cla_porcond_x.f"
/* cla_porcond_x.f -- translated by f2c (version 20100827).
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

#line 1 "cla_porcond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_PORCOND_X computes the infinity norm condition number of op(A)*diag(x) for Hermitian positi
ve-definite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_PORCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_por
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_por
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_por
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_PORCOND_X( UPLO, N, A, LDA, AF, LDAF, X, INFO, */
/*                                    WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * ) */
/*       REAL               RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLA_PORCOND_X Computes the infinity norm condition number of */
/* >    op(A) * diag(X) where X is a COMPLEX vector. */
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
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**H*U or A = L*L**H, as computed by CPOTRF. */
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

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
doublereal cla_porcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
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
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int cpotrs_(char *, integer *, integer *, 
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

#line 169 "cla_porcond_x.f"
    /* Parameter adjustments */
#line 169 "cla_porcond_x.f"
    a_dim1 = *lda;
#line 169 "cla_porcond_x.f"
    a_offset = 1 + a_dim1;
#line 169 "cla_porcond_x.f"
    a -= a_offset;
#line 169 "cla_porcond_x.f"
    af_dim1 = *ldaf;
#line 169 "cla_porcond_x.f"
    af_offset = 1 + af_dim1;
#line 169 "cla_porcond_x.f"
    af -= af_offset;
#line 169 "cla_porcond_x.f"
    --x;
#line 169 "cla_porcond_x.f"
    --work;
#line 169 "cla_porcond_x.f"
    --rwork;
#line 169 "cla_porcond_x.f"

#line 169 "cla_porcond_x.f"
    /* Function Body */
#line 169 "cla_porcond_x.f"
    ret_val = 0.;

#line 171 "cla_porcond_x.f"
    *info = 0;
#line 172 "cla_porcond_x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 173 "cla_porcond_x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "cla_porcond_x.f"
	*info = -1;
#line 175 "cla_porcond_x.f"
    } else if (*n < 0) {
#line 176 "cla_porcond_x.f"
	*info = -2;
#line 177 "cla_porcond_x.f"
    } else if (*lda < max(1,*n)) {
#line 178 "cla_porcond_x.f"
	*info = -4;
#line 179 "cla_porcond_x.f"
    } else if (*ldaf < max(1,*n)) {
#line 180 "cla_porcond_x.f"
	*info = -6;
#line 181 "cla_porcond_x.f"
    }
#line 182 "cla_porcond_x.f"
    if (*info != 0) {
#line 183 "cla_porcond_x.f"
	i__1 = -(*info);
#line 183 "cla_porcond_x.f"
	xerbla_("CLA_PORCOND_X", &i__1, (ftnlen)13);
#line 184 "cla_porcond_x.f"
	return ret_val;
#line 185 "cla_porcond_x.f"
    }
#line 186 "cla_porcond_x.f"
    up = FALSE_;
#line 187 "cla_porcond_x.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 187 "cla_porcond_x.f"
	up = TRUE_;
#line 187 "cla_porcond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 191 "cla_porcond_x.f"
    anorm = 0.;
#line 192 "cla_porcond_x.f"
    if (up) {
#line 193 "cla_porcond_x.f"
	i__1 = *n;
#line 193 "cla_porcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "cla_porcond_x.f"
	    tmp = 0.;
#line 195 "cla_porcond_x.f"
	    i__2 = i__;
#line 195 "cla_porcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 196 "cla_porcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 196 "cla_porcond_x.f"
		i__4 = j;
#line 196 "cla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 196 "cla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 196 "cla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 197 "cla_porcond_x.f"
	    }
#line 198 "cla_porcond_x.f"
	    i__2 = *n;
#line 198 "cla_porcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 199 "cla_porcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 199 "cla_porcond_x.f"
		i__4 = j;
#line 199 "cla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 199 "cla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 199 "cla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 200 "cla_porcond_x.f"
	    }
#line 201 "cla_porcond_x.f"
	    rwork[i__] = tmp;
#line 202 "cla_porcond_x.f"
	    anorm = max(anorm,tmp);
#line 203 "cla_porcond_x.f"
	}
#line 204 "cla_porcond_x.f"
    } else {
#line 205 "cla_porcond_x.f"
	i__1 = *n;
#line 205 "cla_porcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "cla_porcond_x.f"
	    tmp = 0.;
#line 207 "cla_porcond_x.f"
	    i__2 = i__;
#line 207 "cla_porcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 208 "cla_porcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 208 "cla_porcond_x.f"
		i__4 = j;
#line 208 "cla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 208 "cla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 208 "cla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 209 "cla_porcond_x.f"
	    }
#line 210 "cla_porcond_x.f"
	    i__2 = *n;
#line 210 "cla_porcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 211 "cla_porcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 211 "cla_porcond_x.f"
		i__4 = j;
#line 211 "cla_porcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 211 "cla_porcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 211 "cla_porcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 212 "cla_porcond_x.f"
	    }
#line 213 "cla_porcond_x.f"
	    rwork[i__] = tmp;
#line 214 "cla_porcond_x.f"
	    anorm = max(anorm,tmp);
#line 215 "cla_porcond_x.f"
	}
#line 216 "cla_porcond_x.f"
    }

/*     Quick return if possible. */

#line 220 "cla_porcond_x.f"
    if (*n == 0) {
#line 221 "cla_porcond_x.f"
	ret_val = 1.;
#line 222 "cla_porcond_x.f"
	return ret_val;
#line 223 "cla_porcond_x.f"
    } else if (anorm == 0.) {
#line 224 "cla_porcond_x.f"
	return ret_val;
#line 225 "cla_porcond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 229 "cla_porcond_x.f"
    ainvnm = 0.;

#line 231 "cla_porcond_x.f"
    kase = 0;
#line 232 "cla_porcond_x.f"
L10:
#line 233 "cla_porcond_x.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 234 "cla_porcond_x.f"
    if (kase != 0) {
#line 235 "cla_porcond_x.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 239 "cla_porcond_x.f"
	    i__1 = *n;
#line 239 "cla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "cla_porcond_x.f"
		i__2 = i__;
#line 240 "cla_porcond_x.f"
		i__3 = i__;
#line 240 "cla_porcond_x.f"
		i__4 = i__;
#line 240 "cla_porcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 240 "cla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 241 "cla_porcond_x.f"
	    }

#line 243 "cla_porcond_x.f"
	    if (up) {
#line 244 "cla_porcond_x.f"
		cpotrs_("U", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 246 "cla_porcond_x.f"
	    } else {
#line 247 "cla_porcond_x.f"
		cpotrs_("L", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 249 "cla_porcond_x.f"
	    }

/*           Multiply by inv(X). */

#line 253 "cla_porcond_x.f"
	    i__1 = *n;
#line 253 "cla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 254 "cla_porcond_x.f"
		i__2 = i__;
#line 254 "cla_porcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 254 "cla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 255 "cla_porcond_x.f"
	    }
#line 256 "cla_porcond_x.f"
	} else {

/*           Multiply by inv(X**H). */

#line 260 "cla_porcond_x.f"
	    i__1 = *n;
#line 260 "cla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 261 "cla_porcond_x.f"
		i__2 = i__;
#line 261 "cla_porcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 261 "cla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 262 "cla_porcond_x.f"
	    }

#line 264 "cla_porcond_x.f"
	    if (up) {
#line 265 "cla_porcond_x.f"
		cpotrs_("U", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 267 "cla_porcond_x.f"
	    } else {
#line 268 "cla_porcond_x.f"
		cpotrs_("L", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 270 "cla_porcond_x.f"
	    }

/*           Multiply by R. */

#line 274 "cla_porcond_x.f"
	    i__1 = *n;
#line 274 "cla_porcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "cla_porcond_x.f"
		i__2 = i__;
#line 275 "cla_porcond_x.f"
		i__3 = i__;
#line 275 "cla_porcond_x.f"
		i__4 = i__;
#line 275 "cla_porcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 275 "cla_porcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 276 "cla_porcond_x.f"
	    }
#line 277 "cla_porcond_x.f"
	}
#line 278 "cla_porcond_x.f"
	goto L10;
#line 279 "cla_porcond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 283 "cla_porcond_x.f"
    if (ainvnm != 0.) {
#line 283 "cla_porcond_x.f"
	ret_val = 1. / ainvnm;
#line 283 "cla_porcond_x.f"
    }

#line 286 "cla_porcond_x.f"
    return ret_val;

} /* cla_porcond_x__ */

