#line 1 "sla_porcond.f"
/* sla_porcond.f -- translated by f2c (version 20100827).
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

#line 1 "sla_porcond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLA_PORCOND estimates the Skeel condition number for a symmetric positive-definite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_PORCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_por
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_por
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_por
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C, */
/*                                  INFO, WORK, IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LDAF, INFO, CMODE */
/*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), */
/*      $                   C( * ) */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C) */
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
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**T*U or A = L*L**T, as computed by SPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
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
/* >     Workspace. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realPOcomputational */

/*  ===================================================================== */
doublereal sla_porcond__(char *uplo, integer *n, doublereal *a, integer *lda, 
	doublereal *af, integer *ldaf, integer *cmode, doublereal *c__, 
	integer *info, doublereal *work, integer *iwork, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, j;
    static logical up;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int spotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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
/*     .. Array Arguments .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 181 "sla_porcond.f"
    /* Parameter adjustments */
#line 181 "sla_porcond.f"
    a_dim1 = *lda;
#line 181 "sla_porcond.f"
    a_offset = 1 + a_dim1;
#line 181 "sla_porcond.f"
    a -= a_offset;
#line 181 "sla_porcond.f"
    af_dim1 = *ldaf;
#line 181 "sla_porcond.f"
    af_offset = 1 + af_dim1;
#line 181 "sla_porcond.f"
    af -= af_offset;
#line 181 "sla_porcond.f"
    --c__;
#line 181 "sla_porcond.f"
    --work;
#line 181 "sla_porcond.f"
    --iwork;
#line 181 "sla_porcond.f"

#line 181 "sla_porcond.f"
    /* Function Body */
#line 181 "sla_porcond.f"
    ret_val = 0.;

#line 183 "sla_porcond.f"
    *info = 0;
#line 184 "sla_porcond.f"
    if (*n < 0) {
#line 185 "sla_porcond.f"
	*info = -2;
#line 186 "sla_porcond.f"
    }
#line 187 "sla_porcond.f"
    if (*info != 0) {
#line 188 "sla_porcond.f"
	i__1 = -(*info);
#line 188 "sla_porcond.f"
	xerbla_("SLA_PORCOND", &i__1, (ftnlen)11);
#line 189 "sla_porcond.f"
	return ret_val;
#line 190 "sla_porcond.f"
    }
#line 192 "sla_porcond.f"
    if (*n == 0) {
#line 193 "sla_porcond.f"
	ret_val = 1.;
#line 194 "sla_porcond.f"
	return ret_val;
#line 195 "sla_porcond.f"
    }
#line 196 "sla_porcond.f"
    up = FALSE_;
#line 197 "sla_porcond.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 197 "sla_porcond.f"
	up = TRUE_;
#line 197 "sla_porcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 202 "sla_porcond.f"
    if (up) {
#line 203 "sla_porcond.f"
	i__1 = *n;
#line 203 "sla_porcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "sla_porcond.f"
	    tmp = 0.;
#line 205 "sla_porcond.f"
	    if (*cmode == 1) {
#line 206 "sla_porcond.f"
		i__2 = i__;
#line 206 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 207 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 208 "sla_porcond.f"
		}
#line 209 "sla_porcond.f"
		i__2 = *n;
#line 209 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 210 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 211 "sla_porcond.f"
		}
#line 212 "sla_porcond.f"
	    } else if (*cmode == 0) {
#line 213 "sla_porcond.f"
		i__2 = i__;
#line 213 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 214 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 215 "sla_porcond.f"
		}
#line 216 "sla_porcond.f"
		i__2 = *n;
#line 216 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 217 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 218 "sla_porcond.f"
		}
#line 219 "sla_porcond.f"
	    } else {
#line 220 "sla_porcond.f"
		i__2 = i__;
#line 220 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 221 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 222 "sla_porcond.f"
		}
#line 223 "sla_porcond.f"
		i__2 = *n;
#line 223 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 224 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 225 "sla_porcond.f"
		}
#line 226 "sla_porcond.f"
	    }
#line 227 "sla_porcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 228 "sla_porcond.f"
	}
#line 229 "sla_porcond.f"
    } else {
#line 230 "sla_porcond.f"
	i__1 = *n;
#line 230 "sla_porcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 231 "sla_porcond.f"
	    tmp = 0.;
#line 232 "sla_porcond.f"
	    if (*cmode == 1) {
#line 233 "sla_porcond.f"
		i__2 = i__;
#line 233 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 234 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 235 "sla_porcond.f"
		}
#line 236 "sla_porcond.f"
		i__2 = *n;
#line 236 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 237 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 238 "sla_porcond.f"
		}
#line 239 "sla_porcond.f"
	    } else if (*cmode == 0) {
#line 240 "sla_porcond.f"
		i__2 = i__;
#line 240 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 241 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 242 "sla_porcond.f"
		}
#line 243 "sla_porcond.f"
		i__2 = *n;
#line 243 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 244 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 245 "sla_porcond.f"
		}
#line 246 "sla_porcond.f"
	    } else {
#line 247 "sla_porcond.f"
		i__2 = i__;
#line 247 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 248 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 249 "sla_porcond.f"
		}
#line 250 "sla_porcond.f"
		i__2 = *n;
#line 250 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 251 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 252 "sla_porcond.f"
		}
#line 253 "sla_porcond.f"
	    }
#line 254 "sla_porcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 255 "sla_porcond.f"
	}
#line 256 "sla_porcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 260 "sla_porcond.f"
    ainvnm = 0.;
#line 262 "sla_porcond.f"
    kase = 0;
#line 263 "sla_porcond.f"
L10:
#line 264 "sla_porcond.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 265 "sla_porcond.f"
    if (kase != 0) {
#line 266 "sla_porcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 270 "sla_porcond.f"
	    i__1 = *n;
#line 270 "sla_porcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "sla_porcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 272 "sla_porcond.f"
	    }
#line 274 "sla_porcond.f"
	    if (up) {
#line 275 "sla_porcond.f"
		spotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 276 "sla_porcond.f"
	    } else {
#line 277 "sla_porcond.f"
		spotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 278 "sla_porcond.f"
	    }

/*           Multiply by inv(C). */

#line 282 "sla_porcond.f"
	    if (*cmode == 1) {
#line 283 "sla_porcond.f"
		i__1 = *n;
#line 283 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 284 "sla_porcond.f"
		    work[i__] /= c__[i__];
#line 285 "sla_porcond.f"
		}
#line 286 "sla_porcond.f"
	    } else if (*cmode == -1) {
#line 287 "sla_porcond.f"
		i__1 = *n;
#line 287 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 288 "sla_porcond.f"
		    work[i__] *= c__[i__];
#line 289 "sla_porcond.f"
		}
#line 290 "sla_porcond.f"
	    }
#line 291 "sla_porcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 295 "sla_porcond.f"
	    if (*cmode == 1) {
#line 296 "sla_porcond.f"
		i__1 = *n;
#line 296 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "sla_porcond.f"
		    work[i__] /= c__[i__];
#line 298 "sla_porcond.f"
		}
#line 299 "sla_porcond.f"
	    } else if (*cmode == -1) {
#line 300 "sla_porcond.f"
		i__1 = *n;
#line 300 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 301 "sla_porcond.f"
		    work[i__] *= c__[i__];
#line 302 "sla_porcond.f"
		}
#line 303 "sla_porcond.f"
	    }
#line 305 "sla_porcond.f"
	    if (up) {
#line 306 "sla_porcond.f"
		spotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 307 "sla_porcond.f"
	    } else {
#line 308 "sla_porcond.f"
		spotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 309 "sla_porcond.f"
	    }

/*           Multiply by R. */

#line 313 "sla_porcond.f"
	    i__1 = *n;
#line 313 "sla_porcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 314 "sla_porcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 315 "sla_porcond.f"
	    }
#line 316 "sla_porcond.f"
	}
#line 317 "sla_porcond.f"
	goto L10;
#line 318 "sla_porcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 322 "sla_porcond.f"
    if (ainvnm != 0.) {
#line 322 "sla_porcond.f"
	ret_val = 1. / ainvnm;
#line 322 "sla_porcond.f"
    }

#line 325 "sla_porcond.f"
    return ret_val;

} /* sla_porcond__ */

