#line 1 "cla_syrcond_x.f"
/* cla_syrcond_x.f -- translated by f2c (version 20100827).
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

#line 1 "cla_syrcond_x.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_SYRCOND_X computes the infinity norm condition number of op(A)*diag(x) for symmetric indefi
nite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_SYRCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syr
cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syr
cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syr
cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_SYRCOND_X( UPLO, N, A, LDA, AF, LDAF, IPIV, X, */
/*                                    INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
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
/* >    CLA_SYRCOND_X Computes the infinity norm condition number of */
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
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by CSYTRF. */
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
/* >     as determined by CSYTRF. */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
doublereal cla_syrcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
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
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csytrs_(char *, integer *, integer *, 
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

#line 179 "cla_syrcond_x.f"
    /* Parameter adjustments */
#line 179 "cla_syrcond_x.f"
    a_dim1 = *lda;
#line 179 "cla_syrcond_x.f"
    a_offset = 1 + a_dim1;
#line 179 "cla_syrcond_x.f"
    a -= a_offset;
#line 179 "cla_syrcond_x.f"
    af_dim1 = *ldaf;
#line 179 "cla_syrcond_x.f"
    af_offset = 1 + af_dim1;
#line 179 "cla_syrcond_x.f"
    af -= af_offset;
#line 179 "cla_syrcond_x.f"
    --ipiv;
#line 179 "cla_syrcond_x.f"
    --x;
#line 179 "cla_syrcond_x.f"
    --work;
#line 179 "cla_syrcond_x.f"
    --rwork;
#line 179 "cla_syrcond_x.f"

#line 179 "cla_syrcond_x.f"
    /* Function Body */
#line 179 "cla_syrcond_x.f"
    ret_val = 0.;

#line 181 "cla_syrcond_x.f"
    *info = 0;
#line 182 "cla_syrcond_x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 183 "cla_syrcond_x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 184 "cla_syrcond_x.f"
	*info = -1;
#line 185 "cla_syrcond_x.f"
    } else if (*n < 0) {
#line 186 "cla_syrcond_x.f"
	*info = -2;
#line 187 "cla_syrcond_x.f"
    } else if (*lda < max(1,*n)) {
#line 188 "cla_syrcond_x.f"
	*info = -4;
#line 189 "cla_syrcond_x.f"
    } else if (*ldaf < max(1,*n)) {
#line 190 "cla_syrcond_x.f"
	*info = -6;
#line 191 "cla_syrcond_x.f"
    }
#line 192 "cla_syrcond_x.f"
    if (*info != 0) {
#line 193 "cla_syrcond_x.f"
	i__1 = -(*info);
#line 193 "cla_syrcond_x.f"
	xerbla_("CLA_SYRCOND_X", &i__1, (ftnlen)13);
#line 194 "cla_syrcond_x.f"
	return ret_val;
#line 195 "cla_syrcond_x.f"
    }
#line 196 "cla_syrcond_x.f"
    up = FALSE_;
#line 197 "cla_syrcond_x.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 197 "cla_syrcond_x.f"
	up = TRUE_;
#line 197 "cla_syrcond_x.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 201 "cla_syrcond_x.f"
    anorm = 0.;
#line 202 "cla_syrcond_x.f"
    if (up) {
#line 203 "cla_syrcond_x.f"
	i__1 = *n;
#line 203 "cla_syrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "cla_syrcond_x.f"
	    tmp = 0.;
#line 205 "cla_syrcond_x.f"
	    i__2 = i__;
#line 205 "cla_syrcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 206 "cla_syrcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 206 "cla_syrcond_x.f"
		i__4 = j;
#line 206 "cla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 206 "cla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 206 "cla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 207 "cla_syrcond_x.f"
	    }
#line 208 "cla_syrcond_x.f"
	    i__2 = *n;
#line 208 "cla_syrcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 209 "cla_syrcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 209 "cla_syrcond_x.f"
		i__4 = j;
#line 209 "cla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 209 "cla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 209 "cla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 210 "cla_syrcond_x.f"
	    }
#line 211 "cla_syrcond_x.f"
	    rwork[i__] = tmp;
#line 212 "cla_syrcond_x.f"
	    anorm = max(anorm,tmp);
#line 213 "cla_syrcond_x.f"
	}
#line 214 "cla_syrcond_x.f"
    } else {
#line 215 "cla_syrcond_x.f"
	i__1 = *n;
#line 215 "cla_syrcond_x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "cla_syrcond_x.f"
	    tmp = 0.;
#line 217 "cla_syrcond_x.f"
	    i__2 = i__;
#line 217 "cla_syrcond_x.f"
	    for (j = 1; j <= i__2; ++j) {
#line 218 "cla_syrcond_x.f"
		i__3 = i__ + j * a_dim1;
#line 218 "cla_syrcond_x.f"
		i__4 = j;
#line 218 "cla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 218 "cla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 218 "cla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 219 "cla_syrcond_x.f"
	    }
#line 220 "cla_syrcond_x.f"
	    i__2 = *n;
#line 220 "cla_syrcond_x.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 221 "cla_syrcond_x.f"
		i__3 = j + i__ * a_dim1;
#line 221 "cla_syrcond_x.f"
		i__4 = j;
#line 221 "cla_syrcond_x.f"
		z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4]
			.r;
#line 221 "cla_syrcond_x.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 221 "cla_syrcond_x.f"
		tmp += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
			abs(d__2));
#line 222 "cla_syrcond_x.f"
	    }
#line 223 "cla_syrcond_x.f"
	    rwork[i__] = tmp;
#line 224 "cla_syrcond_x.f"
	    anorm = max(anorm,tmp);
#line 225 "cla_syrcond_x.f"
	}
#line 226 "cla_syrcond_x.f"
    }

/*     Quick return if possible. */

#line 230 "cla_syrcond_x.f"
    if (*n == 0) {
#line 231 "cla_syrcond_x.f"
	ret_val = 1.;
#line 232 "cla_syrcond_x.f"
	return ret_val;
#line 233 "cla_syrcond_x.f"
    } else if (anorm == 0.) {
#line 234 "cla_syrcond_x.f"
	return ret_val;
#line 235 "cla_syrcond_x.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 239 "cla_syrcond_x.f"
    ainvnm = 0.;

#line 241 "cla_syrcond_x.f"
    kase = 0;
#line 242 "cla_syrcond_x.f"
L10:
#line 243 "cla_syrcond_x.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 244 "cla_syrcond_x.f"
    if (kase != 0) {
#line 245 "cla_syrcond_x.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 249 "cla_syrcond_x.f"
	    i__1 = *n;
#line 249 "cla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 250 "cla_syrcond_x.f"
		i__2 = i__;
#line 250 "cla_syrcond_x.f"
		i__3 = i__;
#line 250 "cla_syrcond_x.f"
		i__4 = i__;
#line 250 "cla_syrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 250 "cla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 251 "cla_syrcond_x.f"
	    }

#line 253 "cla_syrcond_x.f"
	    if (up) {
#line 254 "cla_syrcond_x.f"
		csytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 256 "cla_syrcond_x.f"
	    } else {
#line 257 "cla_syrcond_x.f"
		csytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 259 "cla_syrcond_x.f"
	    }

/*           Multiply by inv(X). */

#line 263 "cla_syrcond_x.f"
	    i__1 = *n;
#line 263 "cla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 264 "cla_syrcond_x.f"
		i__2 = i__;
#line 264 "cla_syrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 264 "cla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 265 "cla_syrcond_x.f"
	    }
#line 266 "cla_syrcond_x.f"
	} else {

/*           Multiply by inv(X**T). */

#line 270 "cla_syrcond_x.f"
	    i__1 = *n;
#line 270 "cla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "cla_syrcond_x.f"
		i__2 = i__;
#line 271 "cla_syrcond_x.f"
		z_div(&z__1, &work[i__], &x[i__]);
#line 271 "cla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 272 "cla_syrcond_x.f"
	    }

#line 274 "cla_syrcond_x.f"
	    if (up) {
#line 275 "cla_syrcond_x.f"
		csytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 277 "cla_syrcond_x.f"
	    } else {
#line 278 "cla_syrcond_x.f"
		csytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 280 "cla_syrcond_x.f"
	    }

/*           Multiply by R. */

#line 284 "cla_syrcond_x.f"
	    i__1 = *n;
#line 284 "cla_syrcond_x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "cla_syrcond_x.f"
		i__2 = i__;
#line 285 "cla_syrcond_x.f"
		i__3 = i__;
#line 285 "cla_syrcond_x.f"
		i__4 = i__;
#line 285 "cla_syrcond_x.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 285 "cla_syrcond_x.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 286 "cla_syrcond_x.f"
	    }
#line 287 "cla_syrcond_x.f"
	}
#line 288 "cla_syrcond_x.f"
	goto L10;
#line 289 "cla_syrcond_x.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 293 "cla_syrcond_x.f"
    if (ainvnm != 0.) {
#line 293 "cla_syrcond_x.f"
	ret_val = 1. / ainvnm;
#line 293 "cla_syrcond_x.f"
    }

#line 296 "cla_syrcond_x.f"
    return ret_val;

} /* cla_syrcond_x__ */

