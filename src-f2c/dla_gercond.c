#line 1 "dla_gercond.f"
/* dla_gercond.f -- translated by f2c (version 20100827).
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

#line 1 "dla_gercond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLA_GERCOND estimates the Skeel condition number for a general matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_GERCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_ger
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_ger
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_ger
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLA_GERCOND ( TRANS, N, A, LDA, AF, */
/*                                               LDAF, IPIV, CMODE, C, */
/*                                               INFO, WORK, IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            N, LDA, LDAF, INFO, CMODE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * ), */
/*      $                   C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLA_GERCOND estimates the Skeel condition number of op(A) * op2(C) */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          AF is DOUBLE PRECISION array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by DGETRF. */
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
/* >     as computed by DGETRF; row i of the matrix was interchanged */
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
/* >          WORK is DOUBLE PRECISION array, dimension (3*N). */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
doublereal dla_gercond__(char *trans, integer *n, doublereal *a, integer *lda,
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
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
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

#line 192 "dla_gercond.f"
    /* Parameter adjustments */
#line 192 "dla_gercond.f"
    a_dim1 = *lda;
#line 192 "dla_gercond.f"
    a_offset = 1 + a_dim1;
#line 192 "dla_gercond.f"
    a -= a_offset;
#line 192 "dla_gercond.f"
    af_dim1 = *ldaf;
#line 192 "dla_gercond.f"
    af_offset = 1 + af_dim1;
#line 192 "dla_gercond.f"
    af -= af_offset;
#line 192 "dla_gercond.f"
    --ipiv;
#line 192 "dla_gercond.f"
    --c__;
#line 192 "dla_gercond.f"
    --work;
#line 192 "dla_gercond.f"
    --iwork;
#line 192 "dla_gercond.f"

#line 192 "dla_gercond.f"
    /* Function Body */
#line 192 "dla_gercond.f"
    ret_val = 0.;

#line 194 "dla_gercond.f"
    *info = 0;
#line 195 "dla_gercond.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 196 "dla_gercond.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 198 "dla_gercond.f"
	*info = -1;
#line 199 "dla_gercond.f"
    } else if (*n < 0) {
#line 200 "dla_gercond.f"
	*info = -2;
#line 201 "dla_gercond.f"
    } else if (*lda < max(1,*n)) {
#line 202 "dla_gercond.f"
	*info = -4;
#line 203 "dla_gercond.f"
    } else if (*ldaf < max(1,*n)) {
#line 204 "dla_gercond.f"
	*info = -6;
#line 205 "dla_gercond.f"
    }
#line 206 "dla_gercond.f"
    if (*info != 0) {
#line 207 "dla_gercond.f"
	i__1 = -(*info);
#line 207 "dla_gercond.f"
	xerbla_("DLA_GERCOND", &i__1, (ftnlen)11);
#line 208 "dla_gercond.f"
	return ret_val;
#line 209 "dla_gercond.f"
    }
#line 210 "dla_gercond.f"
    if (*n == 0) {
#line 211 "dla_gercond.f"
	ret_val = 1.;
#line 212 "dla_gercond.f"
	return ret_val;
#line 213 "dla_gercond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 218 "dla_gercond.f"
    if (notrans) {
#line 219 "dla_gercond.f"
	i__1 = *n;
#line 219 "dla_gercond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "dla_gercond.f"
	    tmp = 0.;
#line 221 "dla_gercond.f"
	    if (*cmode == 1) {
#line 222 "dla_gercond.f"
		i__2 = *n;
#line 222 "dla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 223 "dla_gercond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 224 "dla_gercond.f"
		}
#line 225 "dla_gercond.f"
	    } else if (*cmode == 0) {
#line 226 "dla_gercond.f"
		i__2 = *n;
#line 226 "dla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 227 "dla_gercond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 228 "dla_gercond.f"
		}
#line 229 "dla_gercond.f"
	    } else {
#line 230 "dla_gercond.f"
		i__2 = *n;
#line 230 "dla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 231 "dla_gercond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 232 "dla_gercond.f"
		}
#line 233 "dla_gercond.f"
	    }
#line 234 "dla_gercond.f"
	    work[(*n << 1) + i__] = tmp;
#line 235 "dla_gercond.f"
	}
#line 236 "dla_gercond.f"
    } else {
#line 237 "dla_gercond.f"
	i__1 = *n;
#line 237 "dla_gercond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 238 "dla_gercond.f"
	    tmp = 0.;
#line 239 "dla_gercond.f"
	    if (*cmode == 1) {
#line 240 "dla_gercond.f"
		i__2 = *n;
#line 240 "dla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 241 "dla_gercond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 242 "dla_gercond.f"
		}
#line 243 "dla_gercond.f"
	    } else if (*cmode == 0) {
#line 244 "dla_gercond.f"
		i__2 = *n;
#line 244 "dla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 245 "dla_gercond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 246 "dla_gercond.f"
		}
#line 247 "dla_gercond.f"
	    } else {
#line 248 "dla_gercond.f"
		i__2 = *n;
#line 248 "dla_gercond.f"
		for (j = 1; j <= i__2; ++j) {
#line 249 "dla_gercond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 250 "dla_gercond.f"
		}
#line 251 "dla_gercond.f"
	    }
#line 252 "dla_gercond.f"
	    work[(*n << 1) + i__] = tmp;
#line 253 "dla_gercond.f"
	}
#line 254 "dla_gercond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 258 "dla_gercond.f"
    ainvnm = 0.;
#line 260 "dla_gercond.f"
    kase = 0;
#line 261 "dla_gercond.f"
L10:
#line 262 "dla_gercond.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 263 "dla_gercond.f"
    if (kase != 0) {
#line 264 "dla_gercond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 268 "dla_gercond.f"
	    i__1 = *n;
#line 268 "dla_gercond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "dla_gercond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 270 "dla_gercond.f"
	    }
#line 272 "dla_gercond.f"
	    if (notrans) {
#line 273 "dla_gercond.f"
		dgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 275 "dla_gercond.f"
	    } else {
#line 276 "dla_gercond.f"
		dgetrs_("Transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1],
			 &work[1], n, info, (ftnlen)9);
#line 278 "dla_gercond.f"
	    }

/*           Multiply by inv(C). */

#line 282 "dla_gercond.f"
	    if (*cmode == 1) {
#line 283 "dla_gercond.f"
		i__1 = *n;
#line 283 "dla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 284 "dla_gercond.f"
		    work[i__] /= c__[i__];
#line 285 "dla_gercond.f"
		}
#line 286 "dla_gercond.f"
	    } else if (*cmode == -1) {
#line 287 "dla_gercond.f"
		i__1 = *n;
#line 287 "dla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 288 "dla_gercond.f"
		    work[i__] *= c__[i__];
#line 289 "dla_gercond.f"
		}
#line 290 "dla_gercond.f"
	    }
#line 291 "dla_gercond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 295 "dla_gercond.f"
	    if (*cmode == 1) {
#line 296 "dla_gercond.f"
		i__1 = *n;
#line 296 "dla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "dla_gercond.f"
		    work[i__] /= c__[i__];
#line 298 "dla_gercond.f"
		}
#line 299 "dla_gercond.f"
	    } else if (*cmode == -1) {
#line 300 "dla_gercond.f"
		i__1 = *n;
#line 300 "dla_gercond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 301 "dla_gercond.f"
		    work[i__] *= c__[i__];
#line 302 "dla_gercond.f"
		}
#line 303 "dla_gercond.f"
	    }
#line 305 "dla_gercond.f"
	    if (notrans) {
#line 306 "dla_gercond.f"
		dgetrs_("Transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1],
			 &work[1], n, info, (ftnlen)9);
#line 308 "dla_gercond.f"
	    } else {
#line 309 "dla_gercond.f"
		dgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 311 "dla_gercond.f"
	    }

/*           Multiply by R. */

#line 315 "dla_gercond.f"
	    i__1 = *n;
#line 315 "dla_gercond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "dla_gercond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 317 "dla_gercond.f"
	    }
#line 318 "dla_gercond.f"
	}
#line 319 "dla_gercond.f"
	goto L10;
#line 320 "dla_gercond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 324 "dla_gercond.f"
    if (ainvnm != 0.) {
#line 324 "dla_gercond.f"
	ret_val = 1. / ainvnm;
#line 324 "dla_gercond.f"
    }

#line 327 "dla_gercond.f"
    return ret_val;

} /* dla_gercond__ */

