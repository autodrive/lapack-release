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

/* > \date December 2016 */

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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 180 "sla_porcond.f"
    /* Parameter adjustments */
#line 180 "sla_porcond.f"
    a_dim1 = *lda;
#line 180 "sla_porcond.f"
    a_offset = 1 + a_dim1;
#line 180 "sla_porcond.f"
    a -= a_offset;
#line 180 "sla_porcond.f"
    af_dim1 = *ldaf;
#line 180 "sla_porcond.f"
    af_offset = 1 + af_dim1;
#line 180 "sla_porcond.f"
    af -= af_offset;
#line 180 "sla_porcond.f"
    --c__;
#line 180 "sla_porcond.f"
    --work;
#line 180 "sla_porcond.f"
    --iwork;
#line 180 "sla_porcond.f"

#line 180 "sla_porcond.f"
    /* Function Body */
#line 180 "sla_porcond.f"
    ret_val = 0.;

#line 182 "sla_porcond.f"
    *info = 0;
#line 183 "sla_porcond.f"
    if (*n < 0) {
#line 184 "sla_porcond.f"
	*info = -2;
#line 185 "sla_porcond.f"
    }
#line 186 "sla_porcond.f"
    if (*info != 0) {
#line 187 "sla_porcond.f"
	i__1 = -(*info);
#line 187 "sla_porcond.f"
	xerbla_("SLA_PORCOND", &i__1, (ftnlen)11);
#line 188 "sla_porcond.f"
	return ret_val;
#line 189 "sla_porcond.f"
    }
#line 191 "sla_porcond.f"
    if (*n == 0) {
#line 192 "sla_porcond.f"
	ret_val = 1.;
#line 193 "sla_porcond.f"
	return ret_val;
#line 194 "sla_porcond.f"
    }
#line 195 "sla_porcond.f"
    up = FALSE_;
#line 196 "sla_porcond.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "sla_porcond.f"
	up = TRUE_;
#line 196 "sla_porcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 201 "sla_porcond.f"
    if (up) {
#line 202 "sla_porcond.f"
	i__1 = *n;
#line 202 "sla_porcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "sla_porcond.f"
	    tmp = 0.;
#line 204 "sla_porcond.f"
	    if (*cmode == 1) {
#line 205 "sla_porcond.f"
		i__2 = i__;
#line 205 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 206 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 207 "sla_porcond.f"
		}
#line 208 "sla_porcond.f"
		i__2 = *n;
#line 208 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 209 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 210 "sla_porcond.f"
		}
#line 211 "sla_porcond.f"
	    } else if (*cmode == 0) {
#line 212 "sla_porcond.f"
		i__2 = i__;
#line 212 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 213 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 214 "sla_porcond.f"
		}
#line 215 "sla_porcond.f"
		i__2 = *n;
#line 215 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 216 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 217 "sla_porcond.f"
		}
#line 218 "sla_porcond.f"
	    } else {
#line 219 "sla_porcond.f"
		i__2 = i__;
#line 219 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 220 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 221 "sla_porcond.f"
		}
#line 222 "sla_porcond.f"
		i__2 = *n;
#line 222 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 223 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 224 "sla_porcond.f"
		}
#line 225 "sla_porcond.f"
	    }
#line 226 "sla_porcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 227 "sla_porcond.f"
	}
#line 228 "sla_porcond.f"
    } else {
#line 229 "sla_porcond.f"
	i__1 = *n;
#line 229 "sla_porcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 230 "sla_porcond.f"
	    tmp = 0.;
#line 231 "sla_porcond.f"
	    if (*cmode == 1) {
#line 232 "sla_porcond.f"
		i__2 = i__;
#line 232 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 233 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 234 "sla_porcond.f"
		}
#line 235 "sla_porcond.f"
		i__2 = *n;
#line 235 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 236 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 237 "sla_porcond.f"
		}
#line 238 "sla_porcond.f"
	    } else if (*cmode == 0) {
#line 239 "sla_porcond.f"
		i__2 = i__;
#line 239 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 240 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 241 "sla_porcond.f"
		}
#line 242 "sla_porcond.f"
		i__2 = *n;
#line 242 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 243 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 244 "sla_porcond.f"
		}
#line 245 "sla_porcond.f"
	    } else {
#line 246 "sla_porcond.f"
		i__2 = i__;
#line 246 "sla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 247 "sla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 248 "sla_porcond.f"
		}
#line 249 "sla_porcond.f"
		i__2 = *n;
#line 249 "sla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 250 "sla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 251 "sla_porcond.f"
		}
#line 252 "sla_porcond.f"
	    }
#line 253 "sla_porcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 254 "sla_porcond.f"
	}
#line 255 "sla_porcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 259 "sla_porcond.f"
    ainvnm = 0.;
#line 261 "sla_porcond.f"
    kase = 0;
#line 262 "sla_porcond.f"
L10:
#line 263 "sla_porcond.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 264 "sla_porcond.f"
    if (kase != 0) {
#line 265 "sla_porcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 269 "sla_porcond.f"
	    i__1 = *n;
#line 269 "sla_porcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 270 "sla_porcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 271 "sla_porcond.f"
	    }
#line 273 "sla_porcond.f"
	    if (up) {
#line 274 "sla_porcond.f"
		spotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 275 "sla_porcond.f"
	    } else {
#line 276 "sla_porcond.f"
		spotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 277 "sla_porcond.f"
	    }

/*           Multiply by inv(C). */

#line 281 "sla_porcond.f"
	    if (*cmode == 1) {
#line 282 "sla_porcond.f"
		i__1 = *n;
#line 282 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 283 "sla_porcond.f"
		    work[i__] /= c__[i__];
#line 284 "sla_porcond.f"
		}
#line 285 "sla_porcond.f"
	    } else if (*cmode == -1) {
#line 286 "sla_porcond.f"
		i__1 = *n;
#line 286 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "sla_porcond.f"
		    work[i__] *= c__[i__];
#line 288 "sla_porcond.f"
		}
#line 289 "sla_porcond.f"
	    }
#line 290 "sla_porcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 294 "sla_porcond.f"
	    if (*cmode == 1) {
#line 295 "sla_porcond.f"
		i__1 = *n;
#line 295 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 296 "sla_porcond.f"
		    work[i__] /= c__[i__];
#line 297 "sla_porcond.f"
		}
#line 298 "sla_porcond.f"
	    } else if (*cmode == -1) {
#line 299 "sla_porcond.f"
		i__1 = *n;
#line 299 "sla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 300 "sla_porcond.f"
		    work[i__] *= c__[i__];
#line 301 "sla_porcond.f"
		}
#line 302 "sla_porcond.f"
	    }
#line 304 "sla_porcond.f"
	    if (up) {
#line 305 "sla_porcond.f"
		spotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 306 "sla_porcond.f"
	    } else {
#line 307 "sla_porcond.f"
		spotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 308 "sla_porcond.f"
	    }

/*           Multiply by R. */

#line 312 "sla_porcond.f"
	    i__1 = *n;
#line 312 "sla_porcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "sla_porcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 314 "sla_porcond.f"
	    }
#line 315 "sla_porcond.f"
	}
#line 316 "sla_porcond.f"
	goto L10;
#line 317 "sla_porcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 321 "sla_porcond.f"
    if (ainvnm != 0.) {
#line 321 "sla_porcond.f"
	ret_val = 1. / ainvnm;
#line 321 "sla_porcond.f"
    }

#line 324 "sla_porcond.f"
    return ret_val;

} /* sla_porcond__ */

