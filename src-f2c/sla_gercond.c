#line 1 "sla_gercond.f"
/* sla_gercond.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gercond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLA_GERCOND estimates the Skeel condition number for a general matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GERCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_ger
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_ger
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_ger
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_GERCOND ( TRANS, N, A, LDA, AF, LDAF, IPIV, */
/*                                   CMODE, C, INFO, WORK, IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, LDA, LDAF, INFO, CMODE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), */
/*      $                   C( * ) */
/*      .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLA_GERCOND estimates the Skeel condition number of op(A) * op2(C) */
/* >    where op2 is determined by CMODE as follows */
/* >    CMODE =  1    op2(C) = C */
/* >    CMODE =  0    op2(C) = I */
/* >    CMODE = -1    op2(C) = inv(C) */
/* >    The Skeel condition number cond(A) = norminf( |inv(A)||A| ) */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          AF is REAL array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by SGETRF. */
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
/* >     as computed by SGETRF; row i of the matrix was interchanged */
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
/* >          WORK is REAL array, dimension (3*N). */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N). */
/* >     Workspace.2 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
doublereal sla_gercond__(char *trans, integer *n, doublereal *a, integer *lda,
	 doublereal *af, integer *ldaf, integer *ipiv, integer *cmode, 
	doublereal *c__, integer *info, doublereal *work, integer *iwork, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int sgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 190 "sla_gercond.f"
    /* Parameter adjustments */
#line 190 "sla_gercond.f"
    a_dim1 = *lda;
#line 190 "sla_gercond.f"
    a_offset = 1 + a_dim1;
#line 190 "sla_gercond.f"
    a -= a_offset;
#line 190 "sla_gercond.f"
    af_dim1 = *ldaf;
#line 190 "sla_gercond.f"
    af_offset = 1 + af_dim1;
#line 190 "sla_gercond.f"
    af -= af_offset;
#line 190 "sla_gercond.f"
    --ipiv;
#line 190 "sla_gercond.f"
    --c__;
#line 190 "sla_gercond.f"
    --work;
#line 190 "sla_gercond.f"
    --iwork;
#line 190 "sla_gercond.f"

#line 190 "sla_gercond.f"
    /* Function Body */
#line 190 "sla_gercond.f"
    ret_val = 0.;

#line 192 "sla_gercond.f"
    *info = 0;
#line 193 "sla_gercond.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 194 "sla_gercond.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 196 "sla_gercond.f"
	*info = -1;
#line 197 "sla_gercond.f"
    } else if (*n < 0) {
#line 198 "sla_gercond.f"
	*info = -2;
#line 199 "sla_gercond.f"
    } else if (*lda < max(1,*n)) {
#line 200 "sla_gercond.f"
	*info = -4;
#line 201 "sla_gercond.f"
    } else if (*ldaf < max(1,*n)) {
#line 202 "sla_gercond.f"
	*info = -6;
#line 203 "sla_gercond.f"
    }
#line 204 "sla_gercond.f"
    if (*info != 0) {
#line 205 "sla_gercond.f"
	i__1 = -(*info);
#line 205 "sla_gercond.f"
	xerbla_("SLA_GERCOND", &i__1, (ftnlen)11);
#line 206 "sla_gercond.f"
	return ret_val;
#line 207 "sla_gercond.f"
    }
#line 208 "sla_gercond.f"
    if (*n == 0) {
#line 209 "sla_gercond.f"
	ret_val = 1.;
#line 210 "sla_gercond.f"
	return ret_val;
#line 211 "sla_gercond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 216 "sla_gercond.f"
    if (notrans) {
#line 217 "sla_gercond.f"
	i__1 = *n;
#line 217 "sla_gercond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "sla_gercond.f"
	    tmp = 0.;
#line 219 "sla_gercond.f"
	    if (*cmode == 1) {
#line 220 "sla_gercond.f"
		i__2 = *n;
#line 220 "sla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 221 "sla_gercond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 222 "sla_gercond.f"
		}
#line 223 "sla_gercond.f"
	    } else if (*cmode == 0) {
#line 224 "sla_gercond.f"
		i__2 = *n;
#line 224 "sla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 225 "sla_gercond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 226 "sla_gercond.f"
		}
#line 227 "sla_gercond.f"
	    } else {
#line 228 "sla_gercond.f"
		i__2 = *n;
#line 228 "sla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 229 "sla_gercond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 230 "sla_gercond.f"
		}
#line 231 "sla_gercond.f"
	    }
#line 232 "sla_gercond.f"
	    work[(*n << 1) + i__] = tmp;
#line 233 "sla_gercond.f"
	}
#line 234 "sla_gercond.f"
    } else {
#line 235 "sla_gercond.f"
	i__1 = *n;
#line 235 "sla_gercond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "sla_gercond.f"
	    tmp = 0.;
#line 237 "sla_gercond.f"
	    if (*cmode == 1) {
#line 238 "sla_gercond.f"
		i__2 = *n;
#line 238 "sla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 239 "sla_gercond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 240 "sla_gercond.f"
		}
#line 241 "sla_gercond.f"
	    } else if (*cmode == 0) {
#line 242 "sla_gercond.f"
		i__2 = *n;
#line 242 "sla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 243 "sla_gercond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 244 "sla_gercond.f"
		}
#line 245 "sla_gercond.f"
	    } else {
#line 246 "sla_gercond.f"
		i__2 = *n;
#line 246 "sla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 247 "sla_gercond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 248 "sla_gercond.f"
		}
#line 249 "sla_gercond.f"
	    }
#line 250 "sla_gercond.f"
	    work[(*n << 1) + i__] = tmp;
#line 251 "sla_gercond.f"
	}
#line 252 "sla_gercond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 256 "sla_gercond.f"
    ainvnm = 0.;
#line 258 "sla_gercond.f"
    kase = 0;
#line 259 "sla_gercond.f"
L10:
#line 260 "sla_gercond.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 261 "sla_gercond.f"
    if (kase != 0) {
#line 262 "sla_gercond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 266 "sla_gercond.f"
	    i__1 = *n;
#line 266 "sla_gercond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "sla_gercond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 268 "sla_gercond.f"
	    }
#line 270 "sla_gercond.f"
	    if (notrans) {
#line 271 "sla_gercond.f"
		sgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 273 "sla_gercond.f"
	    } else {
#line 274 "sla_gercond.f"
		sgetrs_("Transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1],
			 &work[1], n, info, (ftnlen)9);
#line 276 "sla_gercond.f"
	    }

/*           Multiply by inv(C). */

#line 280 "sla_gercond.f"
	    if (*cmode == 1) {
#line 281 "sla_gercond.f"
		i__1 = *n;
#line 281 "sla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 282 "sla_gercond.f"
		    work[i__] /= c__[i__];
#line 283 "sla_gercond.f"
		}
#line 284 "sla_gercond.f"
	    } else if (*cmode == -1) {
#line 285 "sla_gercond.f"
		i__1 = *n;
#line 285 "sla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 286 "sla_gercond.f"
		    work[i__] *= c__[i__];
#line 287 "sla_gercond.f"
		}
#line 288 "sla_gercond.f"
	    }
#line 289 "sla_gercond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 293 "sla_gercond.f"
	    if (*cmode == 1) {
#line 294 "sla_gercond.f"
		i__1 = *n;
#line 294 "sla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 295 "sla_gercond.f"
		    work[i__] /= c__[i__];
#line 296 "sla_gercond.f"
		}
#line 297 "sla_gercond.f"
	    } else if (*cmode == -1) {
#line 298 "sla_gercond.f"
		i__1 = *n;
#line 298 "sla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "sla_gercond.f"
		    work[i__] *= c__[i__];
#line 300 "sla_gercond.f"
		}
#line 301 "sla_gercond.f"
	    }
#line 303 "sla_gercond.f"
	    if (notrans) {
#line 304 "sla_gercond.f"
		sgetrs_("Transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1],
			 &work[1], n, info, (ftnlen)9);
#line 306 "sla_gercond.f"
	    } else {
#line 307 "sla_gercond.f"
		sgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 309 "sla_gercond.f"
	    }

/*           Multiply by R. */

#line 313 "sla_gercond.f"
	    i__1 = *n;
#line 313 "sla_gercond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 314 "sla_gercond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 315 "sla_gercond.f"
	    }
#line 316 "sla_gercond.f"
	}
#line 317 "sla_gercond.f"
	goto L10;
#line 318 "sla_gercond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 322 "sla_gercond.f"
    if (ainvnm != 0.) {
#line 322 "sla_gercond.f"
	ret_val = 1. / ainvnm;
#line 322 "sla_gercond.f"
    }

#line 325 "sla_gercond.f"
    return ret_val;

} /* sla_gercond__ */

