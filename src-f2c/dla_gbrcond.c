#line 1 "dla_gbrcond.f"
/* dla_gbrcond.f -- translated by f2c (version 20100827).
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

#line 1 "dla_gbrcond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLA_GBRCOND estimates the Skeel condition number for a general banded matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_GBRCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_gbr
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_gbr
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_gbr
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, */
/*                                              AFB, LDAFB, IPIV, CMODE, C, */
/*                                              INFO, WORK, IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, LDAB, LDAFB, INFO, KL, KU, CMODE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ), IPIV( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), */
/*      $                   C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C) */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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
/* >          AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by DGBTRF.  U is stored as an upper triangular */
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
/* >     as computed by DGBTRF; row i of the matrix was interchanged */
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
/* >          C is DOUBLE PRECISION array, dimension (N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (5*N). */
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

/* > \date September 2012 */

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
doublereal dla_gbrcond__(char *trans, integer *n, integer *kl, integer *ku, 
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
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen), dgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal ainvnm;
    static logical notrans;


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
/*     .. Executable Statements .. */

#line 210 "dla_gbrcond.f"
    /* Parameter adjustments */
#line 210 "dla_gbrcond.f"
    ab_dim1 = *ldab;
#line 210 "dla_gbrcond.f"
    ab_offset = 1 + ab_dim1;
#line 210 "dla_gbrcond.f"
    ab -= ab_offset;
#line 210 "dla_gbrcond.f"
    afb_dim1 = *ldafb;
#line 210 "dla_gbrcond.f"
    afb_offset = 1 + afb_dim1;
#line 210 "dla_gbrcond.f"
    afb -= afb_offset;
#line 210 "dla_gbrcond.f"
    --ipiv;
#line 210 "dla_gbrcond.f"
    --c__;
#line 210 "dla_gbrcond.f"
    --work;
#line 210 "dla_gbrcond.f"
    --iwork;
#line 210 "dla_gbrcond.f"

#line 210 "dla_gbrcond.f"
    /* Function Body */
#line 210 "dla_gbrcond.f"
    ret_val = 0.;

#line 212 "dla_gbrcond.f"
    *info = 0;
#line 213 "dla_gbrcond.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 214 "dla_gbrcond.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 216 "dla_gbrcond.f"
	*info = -1;
#line 217 "dla_gbrcond.f"
    } else if (*n < 0) {
#line 218 "dla_gbrcond.f"
	*info = -2;
#line 219 "dla_gbrcond.f"
    } else if (*kl < 0 || *kl > *n - 1) {
#line 220 "dla_gbrcond.f"
	*info = -3;
#line 221 "dla_gbrcond.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 222 "dla_gbrcond.f"
	*info = -4;
#line 223 "dla_gbrcond.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 224 "dla_gbrcond.f"
	*info = -6;
#line 225 "dla_gbrcond.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 226 "dla_gbrcond.f"
	*info = -8;
#line 227 "dla_gbrcond.f"
    }
#line 228 "dla_gbrcond.f"
    if (*info != 0) {
#line 229 "dla_gbrcond.f"
	i__1 = -(*info);
#line 229 "dla_gbrcond.f"
	xerbla_("DLA_GBRCOND", &i__1, (ftnlen)11);
#line 230 "dla_gbrcond.f"
	return ret_val;
#line 231 "dla_gbrcond.f"
    }
#line 232 "dla_gbrcond.f"
    if (*n == 0) {
#line 233 "dla_gbrcond.f"
	ret_val = 1.;
#line 234 "dla_gbrcond.f"
	return ret_val;
#line 235 "dla_gbrcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 240 "dla_gbrcond.f"
    kd = *ku + 1;
#line 241 "dla_gbrcond.f"
    ke = *kl + 1;
#line 242 "dla_gbrcond.f"
    if (notrans) {
#line 243 "dla_gbrcond.f"
	i__1 = *n;
#line 243 "dla_gbrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "dla_gbrcond.f"
	    tmp = 0.;
#line 245 "dla_gbrcond.f"
	    if (*cmode == 1) {
/* Computing MAX */
#line 246 "dla_gbrcond.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 246 "dla_gbrcond.f"
		i__4 = i__ + *ku;
#line 246 "dla_gbrcond.f"
		i__3 = min(i__4,*n);
#line 246 "dla_gbrcond.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 247 "dla_gbrcond.f"
		    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] * c__[j], 
			    abs(d__1));
#line 248 "dla_gbrcond.f"
		}
#line 249 "dla_gbrcond.f"
	    } else if (*cmode == 0) {
/* Computing MAX */
#line 250 "dla_gbrcond.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 250 "dla_gbrcond.f"
		i__4 = i__ + *ku;
#line 250 "dla_gbrcond.f"
		i__2 = min(i__4,*n);
#line 250 "dla_gbrcond.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 251 "dla_gbrcond.f"
		    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(d__1));
#line 252 "dla_gbrcond.f"
		}
#line 253 "dla_gbrcond.f"
	    } else {
/* Computing MAX */
#line 254 "dla_gbrcond.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 254 "dla_gbrcond.f"
		i__4 = i__ + *ku;
#line 254 "dla_gbrcond.f"
		i__3 = min(i__4,*n);
#line 254 "dla_gbrcond.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 255 "dla_gbrcond.f"
		    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] / c__[j], 
			    abs(d__1));
#line 256 "dla_gbrcond.f"
		}
#line 257 "dla_gbrcond.f"
	    }
#line 258 "dla_gbrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 259 "dla_gbrcond.f"
	}
#line 260 "dla_gbrcond.f"
    } else {
#line 261 "dla_gbrcond.f"
	i__1 = *n;
#line 261 "dla_gbrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "dla_gbrcond.f"
	    tmp = 0.;
#line 263 "dla_gbrcond.f"
	    if (*cmode == 1) {
/* Computing MAX */
#line 264 "dla_gbrcond.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 264 "dla_gbrcond.f"
		i__4 = i__ + *ku;
#line 264 "dla_gbrcond.f"
		i__2 = min(i__4,*n);
#line 264 "dla_gbrcond.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 265 "dla_gbrcond.f"
		    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] * c__[j], 
			    abs(d__1));
#line 266 "dla_gbrcond.f"
		}
#line 267 "dla_gbrcond.f"
	    } else if (*cmode == 0) {
/* Computing MAX */
#line 268 "dla_gbrcond.f"
		i__2 = i__ - *kl;
/* Computing MIN */
#line 268 "dla_gbrcond.f"
		i__4 = i__ + *ku;
#line 268 "dla_gbrcond.f"
		i__3 = min(i__4,*n);
#line 268 "dla_gbrcond.f"
		for (j = max(i__2,1); j <= i__3; ++j) {
#line 269 "dla_gbrcond.f"
		    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(d__1)
			    );
#line 270 "dla_gbrcond.f"
		}
#line 271 "dla_gbrcond.f"
	    } else {
/* Computing MAX */
#line 272 "dla_gbrcond.f"
		i__3 = i__ - *kl;
/* Computing MIN */
#line 272 "dla_gbrcond.f"
		i__4 = i__ + *ku;
#line 272 "dla_gbrcond.f"
		i__2 = min(i__4,*n);
#line 272 "dla_gbrcond.f"
		for (j = max(i__3,1); j <= i__2; ++j) {
#line 273 "dla_gbrcond.f"
		    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] / c__[j], 
			    abs(d__1));
#line 274 "dla_gbrcond.f"
		}
#line 275 "dla_gbrcond.f"
	    }
#line 276 "dla_gbrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 277 "dla_gbrcond.f"
	}
#line 278 "dla_gbrcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 282 "dla_gbrcond.f"
    ainvnm = 0.;
#line 284 "dla_gbrcond.f"
    kase = 0;
#line 285 "dla_gbrcond.f"
L10:
#line 286 "dla_gbrcond.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 287 "dla_gbrcond.f"
    if (kase != 0) {
#line 288 "dla_gbrcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 292 "dla_gbrcond.f"
	    i__1 = *n;
#line 292 "dla_gbrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 293 "dla_gbrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 294 "dla_gbrcond.f"
	    }
#line 296 "dla_gbrcond.f"
	    if (notrans) {
#line 297 "dla_gbrcond.f"
		dgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 299 "dla_gbrcond.f"
	    } else {
#line 300 "dla_gbrcond.f"
		dgbtrs_("Transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)9);
#line 302 "dla_gbrcond.f"
	    }

/*           Multiply by inv(C). */

#line 306 "dla_gbrcond.f"
	    if (*cmode == 1) {
#line 307 "dla_gbrcond.f"
		i__1 = *n;
#line 307 "dla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "dla_gbrcond.f"
		    work[i__] /= c__[i__];
#line 309 "dla_gbrcond.f"
		}
#line 310 "dla_gbrcond.f"
	    } else if (*cmode == -1) {
#line 311 "dla_gbrcond.f"
		i__1 = *n;
#line 311 "dla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 312 "dla_gbrcond.f"
		    work[i__] *= c__[i__];
#line 313 "dla_gbrcond.f"
		}
#line 314 "dla_gbrcond.f"
	    }
#line 315 "dla_gbrcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 319 "dla_gbrcond.f"
	    if (*cmode == 1) {
#line 320 "dla_gbrcond.f"
		i__1 = *n;
#line 320 "dla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 321 "dla_gbrcond.f"
		    work[i__] /= c__[i__];
#line 322 "dla_gbrcond.f"
		}
#line 323 "dla_gbrcond.f"
	    } else if (*cmode == -1) {
#line 324 "dla_gbrcond.f"
		i__1 = *n;
#line 324 "dla_gbrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 325 "dla_gbrcond.f"
		    work[i__] *= c__[i__];
#line 326 "dla_gbrcond.f"
		}
#line 327 "dla_gbrcond.f"
	    }
#line 329 "dla_gbrcond.f"
	    if (notrans) {
#line 330 "dla_gbrcond.f"
		dgbtrs_("Transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)9);
#line 332 "dla_gbrcond.f"
	    } else {
#line 333 "dla_gbrcond.f"
		dgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info, (ftnlen)12);
#line 335 "dla_gbrcond.f"
	    }

/*           Multiply by R. */

#line 339 "dla_gbrcond.f"
	    i__1 = *n;
#line 339 "dla_gbrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 340 "dla_gbrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 341 "dla_gbrcond.f"
	    }
#line 342 "dla_gbrcond.f"
	}
#line 343 "dla_gbrcond.f"
	goto L10;
#line 344 "dla_gbrcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 348 "dla_gbrcond.f"
    if (ainvnm != 0.) {
#line 348 "dla_gbrcond.f"
	ret_val = 1. / ainvnm;
#line 348 "dla_gbrcond.f"
    }

#line 351 "dla_gbrcond.f"
    return ret_val;

} /* dla_gbrcond__ */

