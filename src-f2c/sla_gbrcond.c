#line 1 "sla_gbrcond.f"
/* sla_gbrcond.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gbrcond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLA_GBRCOND estimates the Skeel condition number for a general banded matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GBRCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbr
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbr
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbr
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, */
/*                                  IPIV, CMODE, C, INFO, WORK, IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, LDAB, LDAFB, INFO, KL, KU, CMODE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ), IPIV( * ) */
/*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), */
/*      $                   C( * ) */
/*      .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C) */
/* >    where op2 is determined by CMODE as follows */
/* >    CMODE =  1    op2(C) = C */
/* >    CMODE =  0    op2(C) = I */
/* >    CMODE = -1    op2(C) = inv(C) */
/* >    The Skeel condition number  cond(A) = norminf( |inv(A)||A| ) */
/* >    is computed by computing scaling factors R such that */
/* >    diag(R)*A*op2(C) is row equilibrated and computing the standard */
/* >    infinity-norm condition number. */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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
/* >          AFB is REAL array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by SGBTRF.  U is stored as an upper triangular */
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
/* >     as computed by SGBTRF; row i of the matrix was interchanged */
/* >     with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] CMODE */
/* > \verbatim */
/* >          CMODE is INTEGER */
/* >     Determines op2(C) in the formula op(A) * op2(C) as follows: */
/* >     CMODE =  1    op2(C) = C */
/* >     CMODE =  0    op2(C) = I */
/* >     CMODE = -1    op2(C) = inv(C) */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >     The vector C in the formula op(A) * op2(C). */
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
/* >          WORK is REAL array, dimension (5*N). */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N). */
/* >     Workspace. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
doublereal sla_gbrcond__(char *trans, integer *n, integer *kl, integer *ku, 
	doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, 
	integer *ipiv, integer *cmode, doublereal *c__, integer *info, 
	doublereal *work, integer *iwork, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, j, kd, ke;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int sgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical notrans;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*    .. */

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
/*     .. Executable Statements .. */

#line 208 "sla_gbrcond.f"
    /* Parameter adjustments */
#line 208 "sla_gbrcond.f"
    ab_dim1 = *ldab;
#line 208 "sla_gbrcond.f"
    ab_offset = 1 + ab_dim1;
#line 208 "sla_gbrcond.f"
    ab -= ab_offset;
#line 208 "sla_gbrcond.f"
    afb_dim1 = *ldafb;
#line 208 "sla_gbrcond.f"
    afb_offset = 1 + afb_dim1;
#line 208 "sla_gbrcond.f"
    afb -= afb_offset;
#line 208 "sla_gbrcond.f"
    --ipiv;
#line 208 "sla_gbrcond.f"
    --c__;
#line 208 "sla_gbrcond.f"
    --work;
#line 208 "sla_gbrcond.f"
    --iwork;
#line 208 "sla_gbrcond.f"

#line 208 "sla_gbrcond.f"
    /* Function Body */
#line 208 "sla_gbrcond.f"
    ret_val = 0.;

#line 210 "sla_gbrcond.f"
    *info = 0;
#line 211 "sla_gbrcond.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 212 "sla_gbrcond.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 214 "sla_gbrcond.f"
	*info = -1;
#line 215 "sla_gbrcond.f"
    } else if (*n < 0) {
#line 216 "sla_gbrcond.f"
	*info = -2;
#line 217 "sla_gbrcond.f"
    } else if (*kl < 0 || *kl > *n - 1) {
#line 218 "sla_gbrcond.f"
	*info = -3;
#line 219 "sla_gbrcond.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 220 "sla_gbrcond.f"
	*info = -4;
#line 221 "sla_gbrcond.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 222 "sla_gbrcond.f"
	*info = -6;
#line 223 "sla_gbrcond.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 224 "sla_gbrcond.f"
	*info = -8;
#line 225 "sla_gbrcond.f"
    }
#line 226 "sla_gbrcond.f"
    if (*info != 0) {
#line 227 "sla_gbrcond.f"
	i__1 = -(*info);
#line 227 "sla_gbrcond.f"
	xerbla_("SLA_GBRCOND", &i__1, (ftnlen)11);
#line 228 "sla_gbrcond.f"
	return ret_val;
#line 229 "sla_gbrcond.f"
    }
#line 230 "sla_gbrcond.f"
    if (*n == 0) {
#line 231 "sla_gbrcond.f"
	ret_val = 1.;
#line 232 "sla_gbrcond.f"
	return ret_val;
#line 233 "sla_gbrcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 238 "sla_gbrcond.f"
    kd = *ku + 1;
#line 239 "sla_gbrcond.f"
    ke = *kl + 1;
#line 240 "sla_gbrcond.f"
    if (notrans) {
#line 241 "sla_gbrcond.f"
	i__1 = *n;
#line 241 "sla_gbrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "sla_gbrcond.f"
	    tmp = 0.;
#line 243 "sla_gbrcond.f"
	    if (*cmode == 1) {
/* Computing MAX */
#line 244 "sla_gbrcond.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 244 "sla_gbrcond.f"
		i__4 = i__ + *ku;
#line 244 "sla_gbrcond.f"
		i__3 = min(i__4,*n);
#line 244 "sla_gbrcond.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 245 "sla_gbrcond.f"
		    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] * c__[j], 
			    abs(d__1));
#line 246 "sla_gbrcond.f"
		}
#line 247 "sla_gbrcond.f"
	    } else if (*cmode == 0) {
/* Computing MAX */
#line 248 "sla_gbrcond.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 248 "sla_gbrcond.f"
		i__4 = i__ + *ku;
#line 248 "sla_gbrcond.f"
		i__2 = min(i__4,*n);
#line 248 "sla_gbrcond.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 249 "sla_gbrcond.f"
		    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(d__1));
#line 250 "sla_gbrcond.f"
		}
#line 251 "sla_gbrcond.f"
	    } else {
/* Computing MAX */
#line 252 "sla_gbrcond.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 252 "sla_gbrcond.f"
		i__4 = i__ + *ku;
#line 252 "sla_gbrcond.f"
		i__3 = min(i__4,*n);
#line 252 "sla_gbrcond.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 253 "sla_gbrcond.f"
		    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] / c__[j], 
			    abs(d__1));
#line 254 "sla_gbrcond.f"
		}
#line 255 "sla_gbrcond.f"
	    }
#line 256 "sla_gbrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 257 "sla_gbrcond.f"
	}
#line 258 "sla_gbrcond.f"
    } else {
#line 259 "sla_gbrcond.f"
	i__1 = *n;
#line 259 "sla_gbrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 260 "sla_gbrcond.f"
	    tmp = 0.;
#line 261 "sla_gbrcond.f"
	    if (*cmode == 1) {
/* Computing MAX */
#line 262 "sla_gbrcond.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 262 "sla_gbrcond.f"
		i__4 = i__ + *ku;
#line 262 "sla_gbrcond.f"
		i__2 = min(i__4,*n);
#line 262 "sla_gbrcond.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 263 "sla_gbrcond.f"
		    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] * c__[j], 
			    abs(d__1));
#line 264 "sla_gbrcond.f"
		}
#line 265 "sla_gbrcond.f"
	    } else if (*cmode == 0) {
/* Computing MAX */
#line 266 "sla_gbrcond.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 266 "sla_gbrcond.f"
		i__4 = i__ + *ku;
#line 266 "sla_gbrcond.f"
		i__3 = min(i__4,*n);
#line 266 "sla_gbrcond.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 267 "sla_gbrcond.f"
		    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(d__1)
			    );
#line 268 "sla_gbrcond.f"
		}
#line 269 "sla_gbrcond.f"
	    } else {
/* Computing MAX */
#line 270 "sla_gbrcond.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 270 "sla_gbrcond.f"
		i__4 = i__ + *ku;
#line 270 "sla_gbrcond.f"
		i__2 = min(i__4,*n);
#line 270 "sla_gbrcond.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 271 "sla_gbrcond.f"
		    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] / c__[j], 
			    abs(d__1));
#line 272 "sla_gbrcond.f"
		}
#line 273 "sla_gbrcond.f"
	    }
#line 274 "sla_gbrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 275 "sla_gbrcond.f"
	}
#line 276 "sla_gbrcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 280 "sla_gbrcond.f"
    ainvnm = 0.;
#line 282 "sla_gbrcond.f"
    kase = 0;
#line 283 "sla_gbrcond.f"
L10:
#line 284 "sla_gbrcond.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 285 "sla_gbrcond.f"
    if (kase != 0) {
#line 286 "sla_gbrcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 290 "sla_gbrcond.f"
	    i__1 = *n;
#line 290 "sla_gbrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "sla_gbrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 292 "sla_gbrcond.f"
	    }
#line 294 "sla_gbrcond.f"
	    if (notrans) {
#line 295 "sla_gbrcond.f"
		sgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 297 "sla_gbrcond.f"
	    } else {
#line 298 "sla_gbrcond.f"
		sgbtrs_("Transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)9);
#line 300 "sla_gbrcond.f"
	    }

/*           Multiply by inv(C). */

#line 304 "sla_gbrcond.f"
	    if (*cmode == 1) {
#line 305 "sla_gbrcond.f"
		i__1 = *n;
#line 305 "sla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "sla_gbrcond.f"
		    work[i__] /= c__[i__];
#line 307 "sla_gbrcond.f"
		}
#line 308 "sla_gbrcond.f"
	    } else if (*cmode == -1) {
#line 309 "sla_gbrcond.f"
		i__1 = *n;
#line 309 "sla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 310 "sla_gbrcond.f"
		    work[i__] *= c__[i__];
#line 311 "sla_gbrcond.f"
		}
#line 312 "sla_gbrcond.f"
	    }
#line 313 "sla_gbrcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 317 "sla_gbrcond.f"
	    if (*cmode == 1) {
#line 318 "sla_gbrcond.f"
		i__1 = *n;
#line 318 "sla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 319 "sla_gbrcond.f"
		    work[i__] /= c__[i__];
#line 320 "sla_gbrcond.f"
		}
#line 321 "sla_gbrcond.f"
	    } else if (*cmode == -1) {
#line 322 "sla_gbrcond.f"
		i__1 = *n;
#line 322 "sla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 323 "sla_gbrcond.f"
		    work[i__] *= c__[i__];
#line 324 "sla_gbrcond.f"
		}
#line 325 "sla_gbrcond.f"
	    }
#line 327 "sla_gbrcond.f"
	    if (notrans) {
#line 328 "sla_gbrcond.f"
		sgbtrs_("Transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)9);
#line 330 "sla_gbrcond.f"
	    } else {
#line 331 "sla_gbrcond.f"
		sgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 333 "sla_gbrcond.f"
	    }

/*           Multiply by R. */

#line 337 "sla_gbrcond.f"
	    i__1 = *n;
#line 337 "sla_gbrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 338 "sla_gbrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 339 "sla_gbrcond.f"
	    }
#line 340 "sla_gbrcond.f"
	}
#line 341 "sla_gbrcond.f"
	goto L10;
#line 342 "sla_gbrcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 346 "sla_gbrcond.f"
    if (ainvnm != 0.) {
#line 346 "sla_gbrcond.f"
	ret_val = 1. / ainvnm;
#line 346 "sla_gbrcond.f"
    }

#line 349 "sla_gbrcond.f"
    return ret_val;

} /* sla_gbrcond__ */

