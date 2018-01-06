#line 1 "zla_syrcond_c.f"
/* zla_syrcond_c.f -- translated by f2c (version 20100827).
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

#line 1 "zla_syrcond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_SYRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for symmetric i
ndefinite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_SYRCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_syr
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_syr
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_syr
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_SYRCOND_C( UPLO, N, A, LDA, AF, */
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
/*       DOUBLE PRECISION   C( * ), RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLA_SYRCOND_C Computes the infinity norm condition number of */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
doublereal zla_syrcond_c__(char *uplo, integer *n, doublecomplex *a, integer *
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

#line 189 "zla_syrcond_c.f"
    /* Parameter adjustments */
#line 189 "zla_syrcond_c.f"
    a_dim1 = *lda;
#line 189 "zla_syrcond_c.f"
    a_offset = 1 + a_dim1;
#line 189 "zla_syrcond_c.f"
    a -= a_offset;
#line 189 "zla_syrcond_c.f"
    af_dim1 = *ldaf;
#line 189 "zla_syrcond_c.f"
    af_offset = 1 + af_dim1;
#line 189 "zla_syrcond_c.f"
    af -= af_offset;
#line 189 "zla_syrcond_c.f"
    --ipiv;
#line 189 "zla_syrcond_c.f"
    --c__;
#line 189 "zla_syrcond_c.f"
    --work;
#line 189 "zla_syrcond_c.f"
    --rwork;
#line 189 "zla_syrcond_c.f"

#line 189 "zla_syrcond_c.f"
    /* Function Body */
#line 189 "zla_syrcond_c.f"
    ret_val = 0.;

#line 191 "zla_syrcond_c.f"
    *info = 0;
#line 192 "zla_syrcond_c.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 193 "zla_syrcond_c.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 194 "zla_syrcond_c.f"
	*info = -1;
#line 195 "zla_syrcond_c.f"
    } else if (*n < 0) {
#line 196 "zla_syrcond_c.f"
	*info = -2;
#line 197 "zla_syrcond_c.f"
    } else if (*lda < max(1,*n)) {
#line 198 "zla_syrcond_c.f"
	*info = -4;
#line 199 "zla_syrcond_c.f"
    } else if (*ldaf < max(1,*n)) {
#line 200 "zla_syrcond_c.f"
	*info = -6;
#line 201 "zla_syrcond_c.f"
    }
#line 202 "zla_syrcond_c.f"
    if (*info != 0) {
#line 203 "zla_syrcond_c.f"
	i__1 = -(*info);
#line 203 "zla_syrcond_c.f"
	xerbla_("ZLA_SYRCOND_C", &i__1, (ftnlen)13);
#line 204 "zla_syrcond_c.f"
	return ret_val;
#line 205 "zla_syrcond_c.f"
    }
#line 206 "zla_syrcond_c.f"
    up = FALSE_;
#line 207 "zla_syrcond_c.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 207 "zla_syrcond_c.f"
	up = TRUE_;
#line 207 "zla_syrcond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 211 "zla_syrcond_c.f"
    anorm = 0.;
#line 212 "zla_syrcond_c.f"
    if (up) {
#line 213 "zla_syrcond_c.f"
	i__1 = *n;
#line 213 "zla_syrcond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "zla_syrcond_c.f"
	    tmp = 0.;
#line 215 "zla_syrcond_c.f"
	    if (*capply) {
#line 216 "zla_syrcond_c.f"
		i__2 = i__;
#line 216 "zla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 217 "zla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 217 "zla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 218 "zla_syrcond_c.f"
		}
#line 219 "zla_syrcond_c.f"
		i__2 = *n;
#line 219 "zla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 220 "zla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 220 "zla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 221 "zla_syrcond_c.f"
		}
#line 222 "zla_syrcond_c.f"
	    } else {
#line 223 "zla_syrcond_c.f"
		i__2 = i__;
#line 223 "zla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 224 "zla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 224 "zla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 225 "zla_syrcond_c.f"
		}
#line 226 "zla_syrcond_c.f"
		i__2 = *n;
#line 226 "zla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 227 "zla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 227 "zla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 228 "zla_syrcond_c.f"
		}
#line 229 "zla_syrcond_c.f"
	    }
#line 230 "zla_syrcond_c.f"
	    rwork[i__] = tmp;
#line 231 "zla_syrcond_c.f"
	    anorm = max(anorm,tmp);
#line 232 "zla_syrcond_c.f"
	}
#line 233 "zla_syrcond_c.f"
    } else {
#line 234 "zla_syrcond_c.f"
	i__1 = *n;
#line 234 "zla_syrcond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "zla_syrcond_c.f"
	    tmp = 0.;
#line 236 "zla_syrcond_c.f"
	    if (*capply) {
#line 237 "zla_syrcond_c.f"
		i__2 = i__;
#line 237 "zla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 238 "zla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 238 "zla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 239 "zla_syrcond_c.f"
		}
#line 240 "zla_syrcond_c.f"
		i__2 = *n;
#line 240 "zla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 241 "zla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 241 "zla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 242 "zla_syrcond_c.f"
		}
#line 243 "zla_syrcond_c.f"
	    } else {
#line 244 "zla_syrcond_c.f"
		i__2 = i__;
#line 244 "zla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 245 "zla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 245 "zla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 246 "zla_syrcond_c.f"
		}
#line 247 "zla_syrcond_c.f"
		i__2 = *n;
#line 247 "zla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 248 "zla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 248 "zla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 249 "zla_syrcond_c.f"
		}
#line 250 "zla_syrcond_c.f"
	    }
#line 251 "zla_syrcond_c.f"
	    rwork[i__] = tmp;
#line 252 "zla_syrcond_c.f"
	    anorm = max(anorm,tmp);
#line 253 "zla_syrcond_c.f"
	}
#line 254 "zla_syrcond_c.f"
    }

/*     Quick return if possible. */

#line 258 "zla_syrcond_c.f"
    if (*n == 0) {
#line 259 "zla_syrcond_c.f"
	ret_val = 1.;
#line 260 "zla_syrcond_c.f"
	return ret_val;
#line 261 "zla_syrcond_c.f"
    } else if (anorm == 0.) {
#line 262 "zla_syrcond_c.f"
	return ret_val;
#line 263 "zla_syrcond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 267 "zla_syrcond_c.f"
    ainvnm = 0.;

#line 269 "zla_syrcond_c.f"
    kase = 0;
#line 270 "zla_syrcond_c.f"
L10:
#line 271 "zla_syrcond_c.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 272 "zla_syrcond_c.f"
    if (kase != 0) {
#line 273 "zla_syrcond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 277 "zla_syrcond_c.f"
	    i__1 = *n;
#line 277 "zla_syrcond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "zla_syrcond_c.f"
		i__2 = i__;
#line 278 "zla_syrcond_c.f"
		i__3 = i__;
#line 278 "zla_syrcond_c.f"
		i__4 = i__;
#line 278 "zla_syrcond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 278 "zla_syrcond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 279 "zla_syrcond_c.f"
	    }

#line 281 "zla_syrcond_c.f"
	    if (up) {
#line 282 "zla_syrcond_c.f"
		zsytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 284 "zla_syrcond_c.f"
	    } else {
#line 285 "zla_syrcond_c.f"
		zsytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 287 "zla_syrcond_c.f"
	    }

/*           Multiply by inv(C). */

#line 291 "zla_syrcond_c.f"
	    if (*capply) {
#line 292 "zla_syrcond_c.f"
		i__1 = *n;
#line 292 "zla_syrcond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 293 "zla_syrcond_c.f"
		    i__2 = i__;
#line 293 "zla_syrcond_c.f"
		    i__3 = i__;
#line 293 "zla_syrcond_c.f"
		    i__4 = i__;
#line 293 "zla_syrcond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 293 "zla_syrcond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 294 "zla_syrcond_c.f"
		}
#line 295 "zla_syrcond_c.f"
	    }
#line 296 "zla_syrcond_c.f"
	} else {

/*           Multiply by inv(C**T). */

#line 300 "zla_syrcond_c.f"
	    if (*capply) {
#line 301 "zla_syrcond_c.f"
		i__1 = *n;
#line 301 "zla_syrcond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "zla_syrcond_c.f"
		    i__2 = i__;
#line 302 "zla_syrcond_c.f"
		    i__3 = i__;
#line 302 "zla_syrcond_c.f"
		    i__4 = i__;
#line 302 "zla_syrcond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 302 "zla_syrcond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 303 "zla_syrcond_c.f"
		}
#line 304 "zla_syrcond_c.f"
	    }

#line 306 "zla_syrcond_c.f"
	    if (up) {
#line 307 "zla_syrcond_c.f"
		zsytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 309 "zla_syrcond_c.f"
	    } else {
#line 310 "zla_syrcond_c.f"
		zsytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 312 "zla_syrcond_c.f"
	    }

/*           Multiply by R. */

#line 316 "zla_syrcond_c.f"
	    i__1 = *n;
#line 316 "zla_syrcond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 317 "zla_syrcond_c.f"
		i__2 = i__;
#line 317 "zla_syrcond_c.f"
		i__3 = i__;
#line 317 "zla_syrcond_c.f"
		i__4 = i__;
#line 317 "zla_syrcond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 317 "zla_syrcond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 318 "zla_syrcond_c.f"
	    }
#line 319 "zla_syrcond_c.f"
	}
#line 320 "zla_syrcond_c.f"
	goto L10;
#line 321 "zla_syrcond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 325 "zla_syrcond_c.f"
    if (ainvnm != 0.) {
#line 325 "zla_syrcond_c.f"
	ret_val = 1. / ainvnm;
#line 325 "zla_syrcond_c.f"
    }

#line 328 "zla_syrcond_c.f"
    return ret_val;

} /* zla_syrcond_c__ */

