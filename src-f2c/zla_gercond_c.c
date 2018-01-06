#line 1 "zla_gercond_c.f"
/* zla_gercond_c.f -- translated by f2c (version 20100827).
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

#line 1 "zla_gercond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_GERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general mat
rices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_GERCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_ger
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_ger
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_ger
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_GERCOND_C( TRANS, N, A, LDA, AF, */
/*                                                LDAF, IPIV, C, CAPPLY, */
/*                                                INFO, WORK, RWORK ) */

/*       .. Scalar Aguments .. */
/*       CHARACTER          TRANS */
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
/* >    ZLA_GERCOND_C computes the infinity norm condition number of */
/* >    op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector. */
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

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
doublereal zla_gercond_c__(char *trans, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *
	c__, logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

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


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Aguments .. */
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
#line 190 "zla_gercond_c.f"
    /* Parameter adjustments */
#line 190 "zla_gercond_c.f"
    a_dim1 = *lda;
#line 190 "zla_gercond_c.f"
    a_offset = 1 + a_dim1;
#line 190 "zla_gercond_c.f"
    a -= a_offset;
#line 190 "zla_gercond_c.f"
    af_dim1 = *ldaf;
#line 190 "zla_gercond_c.f"
    af_offset = 1 + af_dim1;
#line 190 "zla_gercond_c.f"
    af -= af_offset;
#line 190 "zla_gercond_c.f"
    --ipiv;
#line 190 "zla_gercond_c.f"
    --c__;
#line 190 "zla_gercond_c.f"
    --work;
#line 190 "zla_gercond_c.f"
    --rwork;
#line 190 "zla_gercond_c.f"

#line 190 "zla_gercond_c.f"
    /* Function Body */
#line 190 "zla_gercond_c.f"
    ret_val = 0.;

#line 192 "zla_gercond_c.f"
    *info = 0;
#line 193 "zla_gercond_c.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 194 "zla_gercond_c.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 196 "zla_gercond_c.f"
	*info = -1;
#line 197 "zla_gercond_c.f"
    } else if (*n < 0) {
#line 198 "zla_gercond_c.f"
	*info = -2;
#line 199 "zla_gercond_c.f"
    } else if (*lda < max(1,*n)) {
#line 200 "zla_gercond_c.f"
	*info = -4;
#line 201 "zla_gercond_c.f"
    } else if (*ldaf < max(1,*n)) {
#line 202 "zla_gercond_c.f"
	*info = -6;
#line 203 "zla_gercond_c.f"
    }
#line 204 "zla_gercond_c.f"
    if (*info != 0) {
#line 205 "zla_gercond_c.f"
	i__1 = -(*info);
#line 205 "zla_gercond_c.f"
	xerbla_("ZLA_GERCOND_C", &i__1, (ftnlen)13);
#line 206 "zla_gercond_c.f"
	return ret_val;
#line 207 "zla_gercond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 211 "zla_gercond_c.f"
    anorm = 0.;
#line 212 "zla_gercond_c.f"
    if (notrans) {
#line 213 "zla_gercond_c.f"
	i__1 = *n;
#line 213 "zla_gercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "zla_gercond_c.f"
	    tmp = 0.;
#line 215 "zla_gercond_c.f"
	    if (*capply) {
#line 216 "zla_gercond_c.f"
		i__2 = *n;
#line 216 "zla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 217 "zla_gercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 217 "zla_gercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 218 "zla_gercond_c.f"
		}
#line 219 "zla_gercond_c.f"
	    } else {
#line 220 "zla_gercond_c.f"
		i__2 = *n;
#line 220 "zla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 221 "zla_gercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 221 "zla_gercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 222 "zla_gercond_c.f"
		}
#line 223 "zla_gercond_c.f"
	    }
#line 224 "zla_gercond_c.f"
	    rwork[i__] = tmp;
#line 225 "zla_gercond_c.f"
	    anorm = max(anorm,tmp);
#line 226 "zla_gercond_c.f"
	}
#line 227 "zla_gercond_c.f"
    } else {
#line 228 "zla_gercond_c.f"
	i__1 = *n;
#line 228 "zla_gercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "zla_gercond_c.f"
	    tmp = 0.;
#line 230 "zla_gercond_c.f"
	    if (*capply) {
#line 231 "zla_gercond_c.f"
		i__2 = *n;
#line 231 "zla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 232 "zla_gercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 232 "zla_gercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 233 "zla_gercond_c.f"
		}
#line 234 "zla_gercond_c.f"
	    } else {
#line 235 "zla_gercond_c.f"
		i__2 = *n;
#line 235 "zla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 236 "zla_gercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 236 "zla_gercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 237 "zla_gercond_c.f"
		}
#line 238 "zla_gercond_c.f"
	    }
#line 239 "zla_gercond_c.f"
	    rwork[i__] = tmp;
#line 240 "zla_gercond_c.f"
	    anorm = max(anorm,tmp);
#line 241 "zla_gercond_c.f"
	}
#line 242 "zla_gercond_c.f"
    }

/*     Quick return if possible. */

#line 246 "zla_gercond_c.f"
    if (*n == 0) {
#line 247 "zla_gercond_c.f"
	ret_val = 1.;
#line 248 "zla_gercond_c.f"
	return ret_val;
#line 249 "zla_gercond_c.f"
    } else if (anorm == 0.) {
#line 250 "zla_gercond_c.f"
	return ret_val;
#line 251 "zla_gercond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 255 "zla_gercond_c.f"
    ainvnm = 0.;

#line 257 "zla_gercond_c.f"
    kase = 0;
#line 258 "zla_gercond_c.f"
L10:
#line 259 "zla_gercond_c.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 260 "zla_gercond_c.f"
    if (kase != 0) {
#line 261 "zla_gercond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 265 "zla_gercond_c.f"
	    i__1 = *n;
#line 265 "zla_gercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "zla_gercond_c.f"
		i__2 = i__;
#line 266 "zla_gercond_c.f"
		i__3 = i__;
#line 266 "zla_gercond_c.f"
		i__4 = i__;
#line 266 "zla_gercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 266 "zla_gercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 267 "zla_gercond_c.f"
	    }

#line 269 "zla_gercond_c.f"
	    if (notrans) {
#line 270 "zla_gercond_c.f"
		zgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 272 "zla_gercond_c.f"
	    } else {
#line 273 "zla_gercond_c.f"
		zgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 275 "zla_gercond_c.f"
	    }

/*           Multiply by inv(C). */

#line 279 "zla_gercond_c.f"
	    if (*capply) {
#line 280 "zla_gercond_c.f"
		i__1 = *n;
#line 280 "zla_gercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 281 "zla_gercond_c.f"
		    i__2 = i__;
#line 281 "zla_gercond_c.f"
		    i__3 = i__;
#line 281 "zla_gercond_c.f"
		    i__4 = i__;
#line 281 "zla_gercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 281 "zla_gercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 282 "zla_gercond_c.f"
		}
#line 283 "zla_gercond_c.f"
	    }
#line 284 "zla_gercond_c.f"
	} else {

/*           Multiply by inv(C**H). */

#line 288 "zla_gercond_c.f"
	    if (*capply) {
#line 289 "zla_gercond_c.f"
		i__1 = *n;
#line 289 "zla_gercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "zla_gercond_c.f"
		    i__2 = i__;
#line 290 "zla_gercond_c.f"
		    i__3 = i__;
#line 290 "zla_gercond_c.f"
		    i__4 = i__;
#line 290 "zla_gercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 290 "zla_gercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 291 "zla_gercond_c.f"
		}
#line 292 "zla_gercond_c.f"
	    }

#line 294 "zla_gercond_c.f"
	    if (notrans) {
#line 295 "zla_gercond_c.f"
		zgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 297 "zla_gercond_c.f"
	    } else {
#line 298 "zla_gercond_c.f"
		zgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 300 "zla_gercond_c.f"
	    }

/*           Multiply by R. */

#line 304 "zla_gercond_c.f"
	    i__1 = *n;
#line 304 "zla_gercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 305 "zla_gercond_c.f"
		i__2 = i__;
#line 305 "zla_gercond_c.f"
		i__3 = i__;
#line 305 "zla_gercond_c.f"
		i__4 = i__;
#line 305 "zla_gercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 305 "zla_gercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 306 "zla_gercond_c.f"
	    }
#line 307 "zla_gercond_c.f"
	}
#line 308 "zla_gercond_c.f"
	goto L10;
#line 309 "zla_gercond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 313 "zla_gercond_c.f"
    if (ainvnm != 0.) {
#line 313 "zla_gercond_c.f"
	ret_val = 1. / ainvnm;
#line 313 "zla_gercond_c.f"
    }

#line 316 "zla_gercond_c.f"
    return ret_val;

} /* zla_gercond_c__ */

