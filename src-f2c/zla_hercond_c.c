#line 1 "zla_hercond_c.f"
/* zla_hercond_c.f -- translated by f2c (version 20100827).
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

#line 1 "zla_hercond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_HERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian i
ndefinite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_HERCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_her
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_her
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_her
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_HERCOND_C( UPLO, N, A, LDA, AF, */
/*                                                LDAF, IPIV, C, CAPPLY, */
/*                                                INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       LOGICAL            CAPPLY */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/*       DOUBLE PRECISION   C ( * ), RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLA_HERCOND_C computes the infinity norm condition number of */
/* >    op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector. */
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
/* >     On entry, the N-by-N matrix A */
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
/* >     obtain the factor U or L as computed by ZHETRF. */
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
/* >     as determined by CHETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N) */
/* >     The vector C in the formula op(A) * inv(diag(C)). */
/* > \endverbatim */
/* > */
/* > \param[in] CAPPLY */
/* > \verbatim */
/* >          CAPPLY is LOGICAL */
/* >     If .TRUE. then access the vector C in the formula above. */
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

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
doublereal zla_hercond_c__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *c__,
	 logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

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
    extern /* Subroutine */ int zhetrs_(char *, integer *, integer *, 
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

#line 188 "zla_hercond_c.f"
    /* Parameter adjustments */
#line 188 "zla_hercond_c.f"
    a_dim1 = *lda;
#line 188 "zla_hercond_c.f"
    a_offset = 1 + a_dim1;
#line 188 "zla_hercond_c.f"
    a -= a_offset;
#line 188 "zla_hercond_c.f"
    af_dim1 = *ldaf;
#line 188 "zla_hercond_c.f"
    af_offset = 1 + af_dim1;
#line 188 "zla_hercond_c.f"
    af -= af_offset;
#line 188 "zla_hercond_c.f"
    --ipiv;
#line 188 "zla_hercond_c.f"
    --c__;
#line 188 "zla_hercond_c.f"
    --work;
#line 188 "zla_hercond_c.f"
    --rwork;
#line 188 "zla_hercond_c.f"

#line 188 "zla_hercond_c.f"
    /* Function Body */
#line 188 "zla_hercond_c.f"
    ret_val = 0.;

#line 190 "zla_hercond_c.f"
    *info = 0;
#line 191 "zla_hercond_c.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 192 "zla_hercond_c.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 193 "zla_hercond_c.f"
	*info = -1;
#line 194 "zla_hercond_c.f"
    } else if (*n < 0) {
#line 195 "zla_hercond_c.f"
	*info = -2;
#line 196 "zla_hercond_c.f"
    } else if (*lda < max(1,*n)) {
#line 197 "zla_hercond_c.f"
	*info = -4;
#line 198 "zla_hercond_c.f"
    } else if (*ldaf < max(1,*n)) {
#line 199 "zla_hercond_c.f"
	*info = -6;
#line 200 "zla_hercond_c.f"
    }
#line 201 "zla_hercond_c.f"
    if (*info != 0) {
#line 202 "zla_hercond_c.f"
	i__1 = -(*info);
#line 202 "zla_hercond_c.f"
	xerbla_("ZLA_HERCOND_C", &i__1, (ftnlen)13);
#line 203 "zla_hercond_c.f"
	return ret_val;
#line 204 "zla_hercond_c.f"
    }
#line 205 "zla_hercond_c.f"
    up = FALSE_;
#line 206 "zla_hercond_c.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 206 "zla_hercond_c.f"
	up = TRUE_;
#line 206 "zla_hercond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 210 "zla_hercond_c.f"
    anorm = 0.;
#line 211 "zla_hercond_c.f"
    if (up) {
#line 212 "zla_hercond_c.f"
	i__1 = *n;
#line 212 "zla_hercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 213 "zla_hercond_c.f"
	    tmp = 0.;
#line 214 "zla_hercond_c.f"
	    if (*capply) {
#line 215 "zla_hercond_c.f"
		i__2 = i__;
#line 215 "zla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 216 "zla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 216 "zla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 217 "zla_hercond_c.f"
		}
#line 218 "zla_hercond_c.f"
		i__2 = *n;
#line 218 "zla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 219 "zla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 219 "zla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 220 "zla_hercond_c.f"
		}
#line 221 "zla_hercond_c.f"
	    } else {
#line 222 "zla_hercond_c.f"
		i__2 = i__;
#line 222 "zla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 223 "zla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 223 "zla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 224 "zla_hercond_c.f"
		}
#line 225 "zla_hercond_c.f"
		i__2 = *n;
#line 225 "zla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 226 "zla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 226 "zla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 227 "zla_hercond_c.f"
		}
#line 228 "zla_hercond_c.f"
	    }
#line 229 "zla_hercond_c.f"
	    rwork[i__] = tmp;
#line 230 "zla_hercond_c.f"
	    anorm = max(anorm,tmp);
#line 231 "zla_hercond_c.f"
	}
#line 232 "zla_hercond_c.f"
    } else {
#line 233 "zla_hercond_c.f"
	i__1 = *n;
#line 233 "zla_hercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 234 "zla_hercond_c.f"
	    tmp = 0.;
#line 235 "zla_hercond_c.f"
	    if (*capply) {
#line 236 "zla_hercond_c.f"
		i__2 = i__;
#line 236 "zla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 237 "zla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 237 "zla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 238 "zla_hercond_c.f"
		}
#line 239 "zla_hercond_c.f"
		i__2 = *n;
#line 239 "zla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 240 "zla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 240 "zla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 241 "zla_hercond_c.f"
		}
#line 242 "zla_hercond_c.f"
	    } else {
#line 243 "zla_hercond_c.f"
		i__2 = i__;
#line 243 "zla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 244 "zla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 244 "zla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 245 "zla_hercond_c.f"
		}
#line 246 "zla_hercond_c.f"
		i__2 = *n;
#line 246 "zla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 247 "zla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 247 "zla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 248 "zla_hercond_c.f"
		}
#line 249 "zla_hercond_c.f"
	    }
#line 250 "zla_hercond_c.f"
	    rwork[i__] = tmp;
#line 251 "zla_hercond_c.f"
	    anorm = max(anorm,tmp);
#line 252 "zla_hercond_c.f"
	}
#line 253 "zla_hercond_c.f"
    }

/*     Quick return if possible. */

#line 257 "zla_hercond_c.f"
    if (*n == 0) {
#line 258 "zla_hercond_c.f"
	ret_val = 1.;
#line 259 "zla_hercond_c.f"
	return ret_val;
#line 260 "zla_hercond_c.f"
    } else if (anorm == 0.) {
#line 261 "zla_hercond_c.f"
	return ret_val;
#line 262 "zla_hercond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 266 "zla_hercond_c.f"
    ainvnm = 0.;

#line 268 "zla_hercond_c.f"
    kase = 0;
#line 269 "zla_hercond_c.f"
L10:
#line 270 "zla_hercond_c.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 271 "zla_hercond_c.f"
    if (kase != 0) {
#line 272 "zla_hercond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 276 "zla_hercond_c.f"
	    i__1 = *n;
#line 276 "zla_hercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 277 "zla_hercond_c.f"
		i__2 = i__;
#line 277 "zla_hercond_c.f"
		i__3 = i__;
#line 277 "zla_hercond_c.f"
		i__4 = i__;
#line 277 "zla_hercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 277 "zla_hercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 278 "zla_hercond_c.f"
	    }

#line 280 "zla_hercond_c.f"
	    if (up) {
#line 281 "zla_hercond_c.f"
		zhetrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 283 "zla_hercond_c.f"
	    } else {
#line 284 "zla_hercond_c.f"
		zhetrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 286 "zla_hercond_c.f"
	    }

/*           Multiply by inv(C). */

#line 290 "zla_hercond_c.f"
	    if (*capply) {
#line 291 "zla_hercond_c.f"
		i__1 = *n;
#line 291 "zla_hercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 292 "zla_hercond_c.f"
		    i__2 = i__;
#line 292 "zla_hercond_c.f"
		    i__3 = i__;
#line 292 "zla_hercond_c.f"
		    i__4 = i__;
#line 292 "zla_hercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 292 "zla_hercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 293 "zla_hercond_c.f"
		}
#line 294 "zla_hercond_c.f"
	    }
#line 295 "zla_hercond_c.f"
	} else {

/*           Multiply by inv(C**H). */

#line 299 "zla_hercond_c.f"
	    if (*capply) {
#line 300 "zla_hercond_c.f"
		i__1 = *n;
#line 300 "zla_hercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 301 "zla_hercond_c.f"
		    i__2 = i__;
#line 301 "zla_hercond_c.f"
		    i__3 = i__;
#line 301 "zla_hercond_c.f"
		    i__4 = i__;
#line 301 "zla_hercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 301 "zla_hercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 302 "zla_hercond_c.f"
		}
#line 303 "zla_hercond_c.f"
	    }

#line 305 "zla_hercond_c.f"
	    if (up) {
#line 306 "zla_hercond_c.f"
		zhetrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 308 "zla_hercond_c.f"
	    } else {
#line 309 "zla_hercond_c.f"
		zhetrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 311 "zla_hercond_c.f"
	    }

/*           Multiply by R. */

#line 315 "zla_hercond_c.f"
	    i__1 = *n;
#line 315 "zla_hercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "zla_hercond_c.f"
		i__2 = i__;
#line 316 "zla_hercond_c.f"
		i__3 = i__;
#line 316 "zla_hercond_c.f"
		i__4 = i__;
#line 316 "zla_hercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 316 "zla_hercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 317 "zla_hercond_c.f"
	    }
#line 318 "zla_hercond_c.f"
	}
#line 319 "zla_hercond_c.f"
	goto L10;
#line 320 "zla_hercond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 324 "zla_hercond_c.f"
    if (ainvnm != 0.) {
#line 324 "zla_hercond_c.f"
	ret_val = 1. / ainvnm;
#line 324 "zla_hercond_c.f"
    }

#line 327 "zla_hercond_c.f"
    return ret_val;

} /* zla_hercond_c__ */

