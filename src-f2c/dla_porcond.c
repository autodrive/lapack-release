#line 1 "dla_porcond.f"
/* dla_porcond.f -- translated by f2c (version 20100827).
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

#line 1 "dla_porcond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLA_PORCOND estimates the Skeel condition number for a symmetric positive-definite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_PORCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_por
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_por
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_por
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, */
/*                                              CMODE, C, INFO, WORK, */
/*                                              IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LDAF, INFO, CMODE */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * ), */
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
/* >    DLA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C) */
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
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**T*U or A = L*L**T, as computed by DPOTRF. */
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

/* > \date December 2016 */

/* > \ingroup doublePOcomputational */

/*  ===================================================================== */
doublereal dla_porcond__(char *uplo, integer *n, doublereal *a, integer *lda, 
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
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
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

#line 182 "dla_porcond.f"
    /* Parameter adjustments */
#line 182 "dla_porcond.f"
    a_dim1 = *lda;
#line 182 "dla_porcond.f"
    a_offset = 1 + a_dim1;
#line 182 "dla_porcond.f"
    a -= a_offset;
#line 182 "dla_porcond.f"
    af_dim1 = *ldaf;
#line 182 "dla_porcond.f"
    af_offset = 1 + af_dim1;
#line 182 "dla_porcond.f"
    af -= af_offset;
#line 182 "dla_porcond.f"
    --c__;
#line 182 "dla_porcond.f"
    --work;
#line 182 "dla_porcond.f"
    --iwork;
#line 182 "dla_porcond.f"

#line 182 "dla_porcond.f"
    /* Function Body */
#line 182 "dla_porcond.f"
    ret_val = 0.;

#line 184 "dla_porcond.f"
    *info = 0;
#line 185 "dla_porcond.f"
    if (*n < 0) {
#line 186 "dla_porcond.f"
	*info = -2;
#line 187 "dla_porcond.f"
    }
#line 188 "dla_porcond.f"
    if (*info != 0) {
#line 189 "dla_porcond.f"
	i__1 = -(*info);
#line 189 "dla_porcond.f"
	xerbla_("DLA_PORCOND", &i__1, (ftnlen)11);
#line 190 "dla_porcond.f"
	return ret_val;
#line 191 "dla_porcond.f"
    }
#line 193 "dla_porcond.f"
    if (*n == 0) {
#line 194 "dla_porcond.f"
	ret_val = 1.;
#line 195 "dla_porcond.f"
	return ret_val;
#line 196 "dla_porcond.f"
    }
#line 197 "dla_porcond.f"
    up = FALSE_;
#line 198 "dla_porcond.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 198 "dla_porcond.f"
	up = TRUE_;
#line 198 "dla_porcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 203 "dla_porcond.f"
    if (up) {
#line 204 "dla_porcond.f"
	i__1 = *n;
#line 204 "dla_porcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 205 "dla_porcond.f"
	    tmp = 0.;
#line 206 "dla_porcond.f"
	    if (*cmode == 1) {
#line 207 "dla_porcond.f"
		i__2 = i__;
#line 207 "dla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 208 "dla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 209 "dla_porcond.f"
		}
#line 210 "dla_porcond.f"
		i__2 = *n;
#line 210 "dla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 211 "dla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 212 "dla_porcond.f"
		}
#line 213 "dla_porcond.f"
	    } else if (*cmode == 0) {
#line 214 "dla_porcond.f"
		i__2 = i__;
#line 214 "dla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 215 "dla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 216 "dla_porcond.f"
		}
#line 217 "dla_porcond.f"
		i__2 = *n;
#line 217 "dla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 218 "dla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 219 "dla_porcond.f"
		}
#line 220 "dla_porcond.f"
	    } else {
#line 221 "dla_porcond.f"
		i__2 = i__;
#line 221 "dla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 222 "dla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 223 "dla_porcond.f"
		}
#line 224 "dla_porcond.f"
		i__2 = *n;
#line 224 "dla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 225 "dla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 226 "dla_porcond.f"
		}
#line 227 "dla_porcond.f"
	    }
#line 228 "dla_porcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 229 "dla_porcond.f"
	}
#line 230 "dla_porcond.f"
    } else {
#line 231 "dla_porcond.f"
	i__1 = *n;
#line 231 "dla_porcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "dla_porcond.f"
	    tmp = 0.;
#line 233 "dla_porcond.f"
	    if (*cmode == 1) {
#line 234 "dla_porcond.f"
		i__2 = i__;
#line 234 "dla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 235 "dla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 236 "dla_porcond.f"
		}
#line 237 "dla_porcond.f"
		i__2 = *n;
#line 237 "dla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 238 "dla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 239 "dla_porcond.f"
		}
#line 240 "dla_porcond.f"
	    } else if (*cmode == 0) {
#line 241 "dla_porcond.f"
		i__2 = i__;
#line 241 "dla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 242 "dla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 243 "dla_porcond.f"
		}
#line 244 "dla_porcond.f"
		i__2 = *n;
#line 244 "dla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 245 "dla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 246 "dla_porcond.f"
		}
#line 247 "dla_porcond.f"
	    } else {
#line 248 "dla_porcond.f"
		i__2 = i__;
#line 248 "dla_porcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 249 "dla_porcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 250 "dla_porcond.f"
		}
#line 251 "dla_porcond.f"
		i__2 = *n;
#line 251 "dla_porcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 252 "dla_porcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 253 "dla_porcond.f"
		}
#line 254 "dla_porcond.f"
	    }
#line 255 "dla_porcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 256 "dla_porcond.f"
	}
#line 257 "dla_porcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 261 "dla_porcond.f"
    ainvnm = 0.;
#line 263 "dla_porcond.f"
    kase = 0;
#line 264 "dla_porcond.f"
L10:
#line 265 "dla_porcond.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 266 "dla_porcond.f"
    if (kase != 0) {
#line 267 "dla_porcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 271 "dla_porcond.f"
	    i__1 = *n;
#line 271 "dla_porcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 272 "dla_porcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 273 "dla_porcond.f"
	    }
#line 275 "dla_porcond.f"
	    if (up) {
#line 276 "dla_porcond.f"
		dpotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 277 "dla_porcond.f"
	    } else {
#line 278 "dla_porcond.f"
		dpotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 279 "dla_porcond.f"
	    }

/*           Multiply by inv(C). */

#line 283 "dla_porcond.f"
	    if (*cmode == 1) {
#line 284 "dla_porcond.f"
		i__1 = *n;
#line 284 "dla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "dla_porcond.f"
		    work[i__] /= c__[i__];
#line 286 "dla_porcond.f"
		}
#line 287 "dla_porcond.f"
	    } else if (*cmode == -1) {
#line 288 "dla_porcond.f"
		i__1 = *n;
#line 288 "dla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "dla_porcond.f"
		    work[i__] *= c__[i__];
#line 290 "dla_porcond.f"
		}
#line 291 "dla_porcond.f"
	    }
#line 292 "dla_porcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 296 "dla_porcond.f"
	    if (*cmode == 1) {
#line 297 "dla_porcond.f"
		i__1 = *n;
#line 297 "dla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "dla_porcond.f"
		    work[i__] /= c__[i__];
#line 299 "dla_porcond.f"
		}
#line 300 "dla_porcond.f"
	    } else if (*cmode == -1) {
#line 301 "dla_porcond.f"
		i__1 = *n;
#line 301 "dla_porcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "dla_porcond.f"
		    work[i__] *= c__[i__];
#line 303 "dla_porcond.f"
		}
#line 304 "dla_porcond.f"
	    }
#line 306 "dla_porcond.f"
	    if (up) {
#line 307 "dla_porcond.f"
		dpotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 308 "dla_porcond.f"
	    } else {
#line 309 "dla_porcond.f"
		dpotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)5);
#line 310 "dla_porcond.f"
	    }

/*           Multiply by R. */

#line 314 "dla_porcond.f"
	    i__1 = *n;
#line 314 "dla_porcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 315 "dla_porcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 316 "dla_porcond.f"
	    }
#line 317 "dla_porcond.f"
	}
#line 318 "dla_porcond.f"
	goto L10;
#line 319 "dla_porcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 323 "dla_porcond.f"
    if (ainvnm != 0.) {
#line 323 "dla_porcond.f"
	ret_val = 1. / ainvnm;
#line 323 "dla_porcond.f"
    }

#line 326 "dla_porcond.f"
    return ret_val;

} /* dla_porcond__ */

