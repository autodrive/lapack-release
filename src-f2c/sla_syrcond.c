#line 1 "sla_syrcond.f"
/* sla_syrcond.f -- translated by f2c (version 20100827).
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

#line 1 "sla_syrcond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLA_SYRCOND estimates the Skeel condition number for a symmetric indefinite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_SYRCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syr
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syr
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syr
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, CMODE, */
/*                                  C, INFO, WORK, IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LDAF, INFO, CMODE */
/*       .. */
/*       .. Array Arguments */
/*       INTEGER            IWORK( * ), IPIV( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C) */
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
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by SSYTRF. */
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
/* >     as determined by SSYTRF. */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
doublereal sla_syrcond__(char *uplo, integer *n, doublereal *a, integer *lda, 
	doublereal *af, integer *ldaf, integer *ipiv, integer *cmode, 
	doublereal *c__, integer *info, doublereal *work, integer *iwork, 
	ftnlen uplo_len)
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
	     integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    static char normin[1];
    static doublereal smlnum;
    extern /* Subroutine */ int ssytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments */
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

#line 188 "sla_syrcond.f"
    /* Parameter adjustments */
#line 188 "sla_syrcond.f"
    a_dim1 = *lda;
#line 188 "sla_syrcond.f"
    a_offset = 1 + a_dim1;
#line 188 "sla_syrcond.f"
    a -= a_offset;
#line 188 "sla_syrcond.f"
    af_dim1 = *ldaf;
#line 188 "sla_syrcond.f"
    af_offset = 1 + af_dim1;
#line 188 "sla_syrcond.f"
    af -= af_offset;
#line 188 "sla_syrcond.f"
    --ipiv;
#line 188 "sla_syrcond.f"
    --c__;
#line 188 "sla_syrcond.f"
    --work;
#line 188 "sla_syrcond.f"
    --iwork;
#line 188 "sla_syrcond.f"

#line 188 "sla_syrcond.f"
    /* Function Body */
#line 188 "sla_syrcond.f"
    ret_val = 0.;

#line 190 "sla_syrcond.f"
    *info = 0;
#line 191 "sla_syrcond.f"
    if (*n < 0) {
#line 192 "sla_syrcond.f"
	*info = -2;
#line 193 "sla_syrcond.f"
    } else if (*lda < max(1,*n)) {
#line 194 "sla_syrcond.f"
	*info = -4;
#line 195 "sla_syrcond.f"
    } else if (*ldaf < max(1,*n)) {
#line 196 "sla_syrcond.f"
	*info = -6;
#line 197 "sla_syrcond.f"
    }
#line 198 "sla_syrcond.f"
    if (*info != 0) {
#line 199 "sla_syrcond.f"
	i__1 = -(*info);
#line 199 "sla_syrcond.f"
	xerbla_("SLA_SYRCOND", &i__1, (ftnlen)11);
#line 200 "sla_syrcond.f"
	return ret_val;
#line 201 "sla_syrcond.f"
    }
#line 202 "sla_syrcond.f"
    if (*n == 0) {
#line 203 "sla_syrcond.f"
	ret_val = 1.;
#line 204 "sla_syrcond.f"
	return ret_val;
#line 205 "sla_syrcond.f"
    }
#line 206 "sla_syrcond.f"
    up = FALSE_;
#line 207 "sla_syrcond.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 207 "sla_syrcond.f"
	up = TRUE_;
#line 207 "sla_syrcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 212 "sla_syrcond.f"
    if (up) {
#line 213 "sla_syrcond.f"
	i__1 = *n;
#line 213 "sla_syrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "sla_syrcond.f"
	    tmp = 0.;
#line 215 "sla_syrcond.f"
	    if (*cmode == 1) {
#line 216 "sla_syrcond.f"
		i__2 = i__;
#line 216 "sla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 217 "sla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 218 "sla_syrcond.f"
		}
#line 219 "sla_syrcond.f"
		i__2 = *n;
#line 219 "sla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 220 "sla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 221 "sla_syrcond.f"
		}
#line 222 "sla_syrcond.f"
	    } else if (*cmode == 0) {
#line 223 "sla_syrcond.f"
		i__2 = i__;
#line 223 "sla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 224 "sla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 225 "sla_syrcond.f"
		}
#line 226 "sla_syrcond.f"
		i__2 = *n;
#line 226 "sla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 227 "sla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 228 "sla_syrcond.f"
		}
#line 229 "sla_syrcond.f"
	    } else {
#line 230 "sla_syrcond.f"
		i__2 = i__;
#line 230 "sla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 231 "sla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 232 "sla_syrcond.f"
		}
#line 233 "sla_syrcond.f"
		i__2 = *n;
#line 233 "sla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 234 "sla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 235 "sla_syrcond.f"
		}
#line 236 "sla_syrcond.f"
	    }
#line 237 "sla_syrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 238 "sla_syrcond.f"
	}
#line 239 "sla_syrcond.f"
    } else {
#line 240 "sla_syrcond.f"
	i__1 = *n;
#line 240 "sla_syrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "sla_syrcond.f"
	    tmp = 0.;
#line 242 "sla_syrcond.f"
	    if (*cmode == 1) {
#line 243 "sla_syrcond.f"
		i__2 = i__;
#line 243 "sla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 244 "sla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 245 "sla_syrcond.f"
		}
#line 246 "sla_syrcond.f"
		i__2 = *n;
#line 246 "sla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 247 "sla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 248 "sla_syrcond.f"
		}
#line 249 "sla_syrcond.f"
	    } else if (*cmode == 0) {
#line 250 "sla_syrcond.f"
		i__2 = i__;
#line 250 "sla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 251 "sla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 252 "sla_syrcond.f"
		}
#line 253 "sla_syrcond.f"
		i__2 = *n;
#line 253 "sla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 254 "sla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 255 "sla_syrcond.f"
		}
#line 256 "sla_syrcond.f"
	    } else {
#line 257 "sla_syrcond.f"
		i__2 = i__;
#line 257 "sla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 258 "sla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 259 "sla_syrcond.f"
		}
#line 260 "sla_syrcond.f"
		i__2 = *n;
#line 260 "sla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 261 "sla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 262 "sla_syrcond.f"
		}
#line 263 "sla_syrcond.f"
	    }
#line 264 "sla_syrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 265 "sla_syrcond.f"
	}
#line 266 "sla_syrcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 270 "sla_syrcond.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 271 "sla_syrcond.f"
    ainvnm = 0.;
#line 272 "sla_syrcond.f"
    *(unsigned char *)normin = 'N';
#line 274 "sla_syrcond.f"
    kase = 0;
#line 275 "sla_syrcond.f"
L10:
#line 276 "sla_syrcond.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 277 "sla_syrcond.f"
    if (kase != 0) {
#line 278 "sla_syrcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 282 "sla_syrcond.f"
	    i__1 = *n;
#line 282 "sla_syrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 283 "sla_syrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 284 "sla_syrcond.f"
	    }
#line 286 "sla_syrcond.f"
	    if (up) {
#line 287 "sla_syrcond.f"
		ssytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 288 "sla_syrcond.f"
	    } else {
#line 289 "sla_syrcond.f"
		ssytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 290 "sla_syrcond.f"
	    }

/*           Multiply by inv(C). */

#line 294 "sla_syrcond.f"
	    if (*cmode == 1) {
#line 295 "sla_syrcond.f"
		i__1 = *n;
#line 295 "sla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 296 "sla_syrcond.f"
		    work[i__] /= c__[i__];
#line 297 "sla_syrcond.f"
		}
#line 298 "sla_syrcond.f"
	    } else if (*cmode == -1) {
#line 299 "sla_syrcond.f"
		i__1 = *n;
#line 299 "sla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 300 "sla_syrcond.f"
		    work[i__] *= c__[i__];
#line 301 "sla_syrcond.f"
		}
#line 302 "sla_syrcond.f"
	    }
#line 303 "sla_syrcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 307 "sla_syrcond.f"
	    if (*cmode == 1) {
#line 308 "sla_syrcond.f"
		i__1 = *n;
#line 308 "sla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 309 "sla_syrcond.f"
		    work[i__] /= c__[i__];
#line 310 "sla_syrcond.f"
		}
#line 311 "sla_syrcond.f"
	    } else if (*cmode == -1) {
#line 312 "sla_syrcond.f"
		i__1 = *n;
#line 312 "sla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "sla_syrcond.f"
		    work[i__] *= c__[i__];
#line 314 "sla_syrcond.f"
		}
#line 315 "sla_syrcond.f"
	    }
#line 317 "sla_syrcond.f"
	    if (up) {
#line 318 "sla_syrcond.f"
		ssytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 319 "sla_syrcond.f"
	    } else {
#line 320 "sla_syrcond.f"
		ssytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 321 "sla_syrcond.f"
	    }

/*           Multiply by R. */

#line 325 "sla_syrcond.f"
	    i__1 = *n;
#line 325 "sla_syrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 326 "sla_syrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 327 "sla_syrcond.f"
	    }
#line 328 "sla_syrcond.f"
	}

#line 330 "sla_syrcond.f"
	goto L10;
#line 331 "sla_syrcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 335 "sla_syrcond.f"
    if (ainvnm != 0.) {
#line 335 "sla_syrcond.f"
	ret_val = 1. / ainvnm;
#line 335 "sla_syrcond.f"
    }

#line 338 "sla_syrcond.f"
    return ret_val;

} /* sla_syrcond__ */

