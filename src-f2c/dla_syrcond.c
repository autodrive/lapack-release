#line 1 "dla_syrcond.f"
/* dla_syrcond.f -- translated by f2c (version 20100827).
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

#line 1 "dla_syrcond.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLA_SYRCOND estimates the Skeel condition number for a symmetric indefinite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_SYRCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_syr
cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_syr
cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_syr
cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, */
/*                                              IPIV, CMODE, C, INFO, WORK, */
/*                                              IWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LDAF, INFO, CMODE */
/*       .. */
/*       .. Array Arguments */
/*       INTEGER            IWORK( * ), IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C) */
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
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by DSYTRF. */
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
/* >     as determined by DSYTRF. */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
doublereal dla_syrcond__(char *uplo, integer *n, doublereal *a, integer *lda, 
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
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    static char normin[1];
    static doublereal smlnum;
    extern /* Subroutine */ int dsytrs_(char *, integer *, integer *, 
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

#line 190 "dla_syrcond.f"
    /* Parameter adjustments */
#line 190 "dla_syrcond.f"
    a_dim1 = *lda;
#line 190 "dla_syrcond.f"
    a_offset = 1 + a_dim1;
#line 190 "dla_syrcond.f"
    a -= a_offset;
#line 190 "dla_syrcond.f"
    af_dim1 = *ldaf;
#line 190 "dla_syrcond.f"
    af_offset = 1 + af_dim1;
#line 190 "dla_syrcond.f"
    af -= af_offset;
#line 190 "dla_syrcond.f"
    --ipiv;
#line 190 "dla_syrcond.f"
    --c__;
#line 190 "dla_syrcond.f"
    --work;
#line 190 "dla_syrcond.f"
    --iwork;
#line 190 "dla_syrcond.f"

#line 190 "dla_syrcond.f"
    /* Function Body */
#line 190 "dla_syrcond.f"
    ret_val = 0.;

#line 192 "dla_syrcond.f"
    *info = 0;
#line 193 "dla_syrcond.f"
    if (*n < 0) {
#line 194 "dla_syrcond.f"
	*info = -2;
#line 195 "dla_syrcond.f"
    } else if (*lda < max(1,*n)) {
#line 196 "dla_syrcond.f"
	*info = -4;
#line 197 "dla_syrcond.f"
    } else if (*ldaf < max(1,*n)) {
#line 198 "dla_syrcond.f"
	*info = -6;
#line 199 "dla_syrcond.f"
    }
#line 200 "dla_syrcond.f"
    if (*info != 0) {
#line 201 "dla_syrcond.f"
	i__1 = -(*info);
#line 201 "dla_syrcond.f"
	xerbla_("DLA_SYRCOND", &i__1, (ftnlen)11);
#line 202 "dla_syrcond.f"
	return ret_val;
#line 203 "dla_syrcond.f"
    }
#line 204 "dla_syrcond.f"
    if (*n == 0) {
#line 205 "dla_syrcond.f"
	ret_val = 1.;
#line 206 "dla_syrcond.f"
	return ret_val;
#line 207 "dla_syrcond.f"
    }
#line 208 "dla_syrcond.f"
    up = FALSE_;
#line 209 "dla_syrcond.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 209 "dla_syrcond.f"
	up = TRUE_;
#line 209 "dla_syrcond.f"
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

#line 214 "dla_syrcond.f"
    if (up) {
#line 215 "dla_syrcond.f"
	i__1 = *n;
#line 215 "dla_syrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "dla_syrcond.f"
	    tmp = 0.;
#line 217 "dla_syrcond.f"
	    if (*cmode == 1) {
#line 218 "dla_syrcond.f"
		i__2 = i__;
#line 218 "dla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 219 "dla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 220 "dla_syrcond.f"
		}
#line 221 "dla_syrcond.f"
		i__2 = *n;
#line 221 "dla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 222 "dla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 223 "dla_syrcond.f"
		}
#line 224 "dla_syrcond.f"
	    } else if (*cmode == 0) {
#line 225 "dla_syrcond.f"
		i__2 = i__;
#line 225 "dla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 226 "dla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 227 "dla_syrcond.f"
		}
#line 228 "dla_syrcond.f"
		i__2 = *n;
#line 228 "dla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 229 "dla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 230 "dla_syrcond.f"
		}
#line 231 "dla_syrcond.f"
	    } else {
#line 232 "dla_syrcond.f"
		i__2 = i__;
#line 232 "dla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 233 "dla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 234 "dla_syrcond.f"
		}
#line 235 "dla_syrcond.f"
		i__2 = *n;
#line 235 "dla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 236 "dla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 237 "dla_syrcond.f"
		}
#line 238 "dla_syrcond.f"
	    }
#line 239 "dla_syrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 240 "dla_syrcond.f"
	}
#line 241 "dla_syrcond.f"
    } else {
#line 242 "dla_syrcond.f"
	i__1 = *n;
#line 242 "dla_syrcond.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 243 "dla_syrcond.f"
	    tmp = 0.;
#line 244 "dla_syrcond.f"
	    if (*cmode == 1) {
#line 245 "dla_syrcond.f"
		i__2 = i__;
#line 245 "dla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 246 "dla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], abs(d__1));
#line 247 "dla_syrcond.f"
		}
#line 248 "dla_syrcond.f"
		i__2 = *n;
#line 248 "dla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 249 "dla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], abs(d__1));
#line 250 "dla_syrcond.f"
		}
#line 251 "dla_syrcond.f"
	    } else if (*cmode == 0) {
#line 252 "dla_syrcond.f"
		i__2 = i__;
#line 252 "dla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 253 "dla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 254 "dla_syrcond.f"
		}
#line 255 "dla_syrcond.f"
		i__2 = *n;
#line 255 "dla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 256 "dla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 257 "dla_syrcond.f"
		}
#line 258 "dla_syrcond.f"
	    } else {
#line 259 "dla_syrcond.f"
		i__2 = i__;
#line 259 "dla_syrcond.f"
		for (j = 1; j <= i__2; ++j) {
#line 260 "dla_syrcond.f"
		    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], abs(d__1));
#line 261 "dla_syrcond.f"
		}
#line 262 "dla_syrcond.f"
		i__2 = *n;
#line 262 "dla_syrcond.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 263 "dla_syrcond.f"
		    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], abs(d__1));
#line 264 "dla_syrcond.f"
		}
#line 265 "dla_syrcond.f"
	    }
#line 266 "dla_syrcond.f"
	    work[(*n << 1) + i__] = tmp;
#line 267 "dla_syrcond.f"
	}
#line 268 "dla_syrcond.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 272 "dla_syrcond.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 273 "dla_syrcond.f"
    ainvnm = 0.;
#line 274 "dla_syrcond.f"
    *(unsigned char *)normin = 'N';
#line 276 "dla_syrcond.f"
    kase = 0;
#line 277 "dla_syrcond.f"
L10:
#line 278 "dla_syrcond.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 279 "dla_syrcond.f"
    if (kase != 0) {
#line 280 "dla_syrcond.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 284 "dla_syrcond.f"
	    i__1 = *n;
#line 284 "dla_syrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "dla_syrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 286 "dla_syrcond.f"
	    }
#line 288 "dla_syrcond.f"
	    if (up) {
#line 289 "dla_syrcond.f"
		dsytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 290 "dla_syrcond.f"
	    } else {
#line 291 "dla_syrcond.f"
		dsytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 292 "dla_syrcond.f"
	    }

/*           Multiply by inv(C). */

#line 296 "dla_syrcond.f"
	    if (*cmode == 1) {
#line 297 "dla_syrcond.f"
		i__1 = *n;
#line 297 "dla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "dla_syrcond.f"
		    work[i__] /= c__[i__];
#line 299 "dla_syrcond.f"
		}
#line 300 "dla_syrcond.f"
	    } else if (*cmode == -1) {
#line 301 "dla_syrcond.f"
		i__1 = *n;
#line 301 "dla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "dla_syrcond.f"
		    work[i__] *= c__[i__];
#line 303 "dla_syrcond.f"
		}
#line 304 "dla_syrcond.f"
	    }
#line 305 "dla_syrcond.f"
	} else {

/*           Multiply by inv(C**T). */

#line 309 "dla_syrcond.f"
	    if (*cmode == 1) {
#line 310 "dla_syrcond.f"
		i__1 = *n;
#line 310 "dla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 311 "dla_syrcond.f"
		    work[i__] /= c__[i__];
#line 312 "dla_syrcond.f"
		}
#line 313 "dla_syrcond.f"
	    } else if (*cmode == -1) {
#line 314 "dla_syrcond.f"
		i__1 = *n;
#line 314 "dla_syrcond.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 315 "dla_syrcond.f"
		    work[i__] *= c__[i__];
#line 316 "dla_syrcond.f"
		}
#line 317 "dla_syrcond.f"
	    }
#line 319 "dla_syrcond.f"
	    if (up) {
#line 320 "dla_syrcond.f"
		dsytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 321 "dla_syrcond.f"
	    } else {
#line 322 "dla_syrcond.f"
		dsytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 323 "dla_syrcond.f"
	    }

/*           Multiply by R. */

#line 327 "dla_syrcond.f"
	    i__1 = *n;
#line 327 "dla_syrcond.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 328 "dla_syrcond.f"
		work[i__] *= work[(*n << 1) + i__];
#line 329 "dla_syrcond.f"
	    }
#line 330 "dla_syrcond.f"
	}

#line 332 "dla_syrcond.f"
	goto L10;
#line 333 "dla_syrcond.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 337 "dla_syrcond.f"
    if (ainvnm != 0.) {
#line 337 "dla_syrcond.f"
	ret_val = 1. / ainvnm;
#line 337 "dla_syrcond.f"
    }

#line 340 "dla_syrcond.f"
    return ret_val;

} /* dla_syrcond__ */

