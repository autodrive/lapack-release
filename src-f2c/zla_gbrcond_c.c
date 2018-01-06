#line 1 "zla_gbrcond_c.f"
/* zla_gbrcond_c.f -- translated by f2c (version 20100827).
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

#line 1 "zla_gbrcond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLA_GBRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general ban
ded matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_GBRCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gbr
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gbr
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gbr
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_GBRCOND_C( TRANS, N, KL, KU, AB, */
/*                                                LDAB, AFB, LDAFB, IPIV, */
/*                                                C, CAPPLY, INFO, WORK, */
/*                                                RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       LOGICAL            CAPPLY */
/*       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ) */
/*       DOUBLE PRECISION   C( * ), RWORK( * ) */



/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLA_GBRCOND_C Computes the infinity norm condition number of */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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
/* >          AFB is COMPLEX*16 array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by ZGBTRF.  U is stored as an upper triangular */
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
/* >     as computed by ZGBTRF; row i of the matrix was interchanged */
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

/* > \ingroup complex16GBcomputational */

/*  ===================================================================== */
doublereal zla_gbrcond_c__(char *trans, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, 
	integer *ipiv, doublereal *c__, logical *capply, integer *info, 
	doublecomplex *work, doublereal *rwork, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, kd, ke;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int zgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static logical notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */


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
#line 210 "zla_gbrcond_c.f"
    /* Parameter adjustments */
#line 210 "zla_gbrcond_c.f"
    ab_dim1 = *ldab;
#line 210 "zla_gbrcond_c.f"
    ab_offset = 1 + ab_dim1;
#line 210 "zla_gbrcond_c.f"
    ab -= ab_offset;
#line 210 "zla_gbrcond_c.f"
    afb_dim1 = *ldafb;
#line 210 "zla_gbrcond_c.f"
    afb_offset = 1 + afb_dim1;
#line 210 "zla_gbrcond_c.f"
    afb -= afb_offset;
#line 210 "zla_gbrcond_c.f"
    --ipiv;
#line 210 "zla_gbrcond_c.f"
    --c__;
#line 210 "zla_gbrcond_c.f"
    --work;
#line 210 "zla_gbrcond_c.f"
    --rwork;
#line 210 "zla_gbrcond_c.f"

#line 210 "zla_gbrcond_c.f"
    /* Function Body */
#line 210 "zla_gbrcond_c.f"
    ret_val = 0.;

#line 212 "zla_gbrcond_c.f"
    *info = 0;
#line 213 "zla_gbrcond_c.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 214 "zla_gbrcond_c.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 216 "zla_gbrcond_c.f"
	*info = -1;
#line 217 "zla_gbrcond_c.f"
    } else if (*n < 0) {
#line 218 "zla_gbrcond_c.f"
	*info = -2;
#line 219 "zla_gbrcond_c.f"
    } else if (*kl < 0 || *kl > *n - 1) {
#line 220 "zla_gbrcond_c.f"
	*info = -3;
#line 221 "zla_gbrcond_c.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 222 "zla_gbrcond_c.f"
	*info = -4;
#line 223 "zla_gbrcond_c.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 224 "zla_gbrcond_c.f"
	*info = -6;
#line 225 "zla_gbrcond_c.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 226 "zla_gbrcond_c.f"
	*info = -8;
#line 227 "zla_gbrcond_c.f"
    }
#line 228 "zla_gbrcond_c.f"
    if (*info != 0) {
#line 229 "zla_gbrcond_c.f"
	i__1 = -(*info);
#line 229 "zla_gbrcond_c.f"
	xerbla_("ZLA_GBRCOND_C", &i__1, (ftnlen)13);
#line 230 "zla_gbrcond_c.f"
	return ret_val;
#line 231 "zla_gbrcond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 235 "zla_gbrcond_c.f"
    anorm = 0.;
#line 236 "zla_gbrcond_c.f"
    kd = *ku + 1;
#line 237 "zla_gbrcond_c.f"
    ke = *kl + 1;
#line 238 "zla_gbrcond_c.f"
    if (notrans) {
#line 239 "zla_gbrcond_c.f"
	i__1 = *n;
#line 239 "zla_gbrcond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "zla_gbrcond_c.f"
	    tmp = 0.;
#line 241 "zla_gbrcond_c.f"
	    if (*capply) {
/* Computing MAX */
#line 242 "zla_gbrcond_c.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 242 "zla_gbrcond_c.f"
		i__4 = i__ + *ku;
#line 242 "zla_gbrcond_c.f"
		i__3 = min(i__4,*n);
#line 242 "zla_gbrcond_c.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 243 "zla_gbrcond_c.f"
		    i__2 = kd + i__ - j + j * ab_dim1;
#line 243 "zla_gbrcond_c.f"
		    tmp += ((d__1 = ab[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    ab[kd + i__ - j + j * ab_dim1]), abs(d__2))) / 
			    c__[j];
#line 244 "zla_gbrcond_c.f"
		}
#line 245 "zla_gbrcond_c.f"
	    } else {
/* Computing MAX */
#line 246 "zla_gbrcond_c.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 246 "zla_gbrcond_c.f"
		i__4 = i__ + *ku;
#line 246 "zla_gbrcond_c.f"
		i__2 = min(i__4,*n);
#line 246 "zla_gbrcond_c.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 247 "zla_gbrcond_c.f"
		    i__3 = kd + i__ - j + j * ab_dim1;
#line 247 "zla_gbrcond_c.f"
		    tmp += (d__1 = ab[i__3].r, abs(d__1)) + (d__2 = d_imag(&
			    ab[kd + i__ - j + j * ab_dim1]), abs(d__2));
#line 248 "zla_gbrcond_c.f"
		}
#line 249 "zla_gbrcond_c.f"
	    }
#line 250 "zla_gbrcond_c.f"
	    rwork[i__] = tmp;
#line 251 "zla_gbrcond_c.f"
	    anorm = max(anorm,tmp);
#line 252 "zla_gbrcond_c.f"
	}
#line 253 "zla_gbrcond_c.f"
    } else {
#line 254 "zla_gbrcond_c.f"
	i__1 = *n;
#line 254 "zla_gbrcond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 255 "zla_gbrcond_c.f"
	    tmp = 0.;
#line 256 "zla_gbrcond_c.f"
	    if (*capply) {
/* Computing MAX */
#line 257 "zla_gbrcond_c.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 257 "zla_gbrcond_c.f"
		i__4 = i__ + *ku;
#line 257 "zla_gbrcond_c.f"
		i__3 = min(i__4,*n);
#line 257 "zla_gbrcond_c.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 258 "zla_gbrcond_c.f"
		    i__2 = ke - i__ + j + i__ * ab_dim1;
#line 258 "zla_gbrcond_c.f"
		    tmp += ((d__1 = ab[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    ab[ke - i__ + j + i__ * ab_dim1]), abs(d__2))) / 
			    c__[j];
#line 259 "zla_gbrcond_c.f"
		}
#line 260 "zla_gbrcond_c.f"
	    } else {
/* Computing MAX */
#line 261 "zla_gbrcond_c.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 261 "zla_gbrcond_c.f"
		i__4 = i__ + *ku;
#line 261 "zla_gbrcond_c.f"
		i__2 = min(i__4,*n);
#line 261 "zla_gbrcond_c.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 262 "zla_gbrcond_c.f"
		    i__3 = ke - i__ + j + i__ * ab_dim1;
#line 262 "zla_gbrcond_c.f"
		    tmp += (d__1 = ab[i__3].r, abs(d__1)) + (d__2 = d_imag(&
			    ab[ke - i__ + j + i__ * ab_dim1]), abs(d__2));
#line 263 "zla_gbrcond_c.f"
		}
#line 264 "zla_gbrcond_c.f"
	    }
#line 265 "zla_gbrcond_c.f"
	    rwork[i__] = tmp;
#line 266 "zla_gbrcond_c.f"
	    anorm = max(anorm,tmp);
#line 267 "zla_gbrcond_c.f"
	}
#line 268 "zla_gbrcond_c.f"
    }

/*     Quick return if possible. */

#line 272 "zla_gbrcond_c.f"
    if (*n == 0) {
#line 273 "zla_gbrcond_c.f"
	ret_val = 1.;
#line 274 "zla_gbrcond_c.f"
	return ret_val;
#line 275 "zla_gbrcond_c.f"
    } else if (anorm == 0.) {
#line 276 "zla_gbrcond_c.f"
	return ret_val;
#line 277 "zla_gbrcond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 281 "zla_gbrcond_c.f"
    ainvnm = 0.;

#line 283 "zla_gbrcond_c.f"
    kase = 0;
#line 284 "zla_gbrcond_c.f"
L10:
#line 285 "zla_gbrcond_c.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 286 "zla_gbrcond_c.f"
    if (kase != 0) {
#line 287 "zla_gbrcond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 291 "zla_gbrcond_c.f"
	    i__1 = *n;
#line 291 "zla_gbrcond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 292 "zla_gbrcond_c.f"
		i__2 = i__;
#line 292 "zla_gbrcond_c.f"
		i__3 = i__;
#line 292 "zla_gbrcond_c.f"
		i__4 = i__;
#line 292 "zla_gbrcond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 292 "zla_gbrcond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 293 "zla_gbrcond_c.f"
	    }

#line 295 "zla_gbrcond_c.f"
	    if (notrans) {
#line 296 "zla_gbrcond_c.f"
		zgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 298 "zla_gbrcond_c.f"
	    } else {
#line 299 "zla_gbrcond_c.f"
		zgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info, (
			ftnlen)19);
#line 301 "zla_gbrcond_c.f"
	    }

/*           Multiply by inv(C). */

#line 305 "zla_gbrcond_c.f"
	    if (*capply) {
#line 306 "zla_gbrcond_c.f"
		i__1 = *n;
#line 306 "zla_gbrcond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 307 "zla_gbrcond_c.f"
		    i__2 = i__;
#line 307 "zla_gbrcond_c.f"
		    i__3 = i__;
#line 307 "zla_gbrcond_c.f"
		    i__4 = i__;
#line 307 "zla_gbrcond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 307 "zla_gbrcond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 308 "zla_gbrcond_c.f"
		}
#line 309 "zla_gbrcond_c.f"
	    }
#line 310 "zla_gbrcond_c.f"
	} else {

/*           Multiply by inv(C**H). */

#line 314 "zla_gbrcond_c.f"
	    if (*capply) {
#line 315 "zla_gbrcond_c.f"
		i__1 = *n;
#line 315 "zla_gbrcond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "zla_gbrcond_c.f"
		    i__2 = i__;
#line 316 "zla_gbrcond_c.f"
		    i__3 = i__;
#line 316 "zla_gbrcond_c.f"
		    i__4 = i__;
#line 316 "zla_gbrcond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 316 "zla_gbrcond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 317 "zla_gbrcond_c.f"
		}
#line 318 "zla_gbrcond_c.f"
	    }

#line 320 "zla_gbrcond_c.f"
	    if (notrans) {
#line 321 "zla_gbrcond_c.f"
		zgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info, (
			ftnlen)19);
#line 323 "zla_gbrcond_c.f"
	    } else {
#line 324 "zla_gbrcond_c.f"
		zgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 326 "zla_gbrcond_c.f"
	    }

/*           Multiply by R. */

#line 330 "zla_gbrcond_c.f"
	    i__1 = *n;
#line 330 "zla_gbrcond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 331 "zla_gbrcond_c.f"
		i__2 = i__;
#line 331 "zla_gbrcond_c.f"
		i__3 = i__;
#line 331 "zla_gbrcond_c.f"
		i__4 = i__;
#line 331 "zla_gbrcond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 331 "zla_gbrcond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 332 "zla_gbrcond_c.f"
	    }
#line 333 "zla_gbrcond_c.f"
	}
#line 334 "zla_gbrcond_c.f"
	goto L10;
#line 335 "zla_gbrcond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 339 "zla_gbrcond_c.f"
    if (ainvnm != 0.) {
#line 339 "zla_gbrcond_c.f"
	ret_val = 1. / ainvnm;
#line 339 "zla_gbrcond_c.f"
    }

#line 342 "zla_gbrcond_c.f"
    return ret_val;

} /* zla_gbrcond_c__ */

