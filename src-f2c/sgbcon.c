#line 1 "sgbcon.f"
/* sgbcon.f -- translated by f2c (version 20100827).
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

#line 1 "sgbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SGBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, KL, KU, LDAB, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBCON estimates the reciprocal of the condition number of a real */
/* > general band matrix A, in either the 1-norm or the infinity-norm, */
/* > using the LU factorization computed by SGBTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as */
/* >    RCOND = 1 / ( norm(A) * norm(inv(A)) ). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies whether the 1-norm condition number or the */
/* >          infinity-norm condition number is required: */
/* >          = '1' or 'O':  1-norm; */
/* >          = 'I':         Infinity-norm. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by SGBTRF.  U is stored as an upper triangular band */
/* >          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
/* >          the multipliers used during the factorization are stored in */
/* >          rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= N, row i of the matrix was */
/* >          interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
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
/* Subroutine */ int sgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublereal *work, integer *iwork, integer *info, 
	ftnlen norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t;
    static integer kd, lm, jp, ix, kase;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical lnoti;
    extern /* Subroutine */ int srscl_(integer *, doublereal *, doublereal *, 
	    integer *), saxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *), slacn2_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal ainvnm;
    extern /* Subroutine */ int slatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static logical onenrm;
    static char normin[1];
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
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

/*     Test the input parameters. */

#line 195 "sgbcon.f"
    /* Parameter adjustments */
#line 195 "sgbcon.f"
    ab_dim1 = *ldab;
#line 195 "sgbcon.f"
    ab_offset = 1 + ab_dim1;
#line 195 "sgbcon.f"
    ab -= ab_offset;
#line 195 "sgbcon.f"
    --ipiv;
#line 195 "sgbcon.f"
    --work;
#line 195 "sgbcon.f"
    --iwork;
#line 195 "sgbcon.f"

#line 195 "sgbcon.f"
    /* Function Body */
#line 195 "sgbcon.f"
    *info = 0;
#line 196 "sgbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 197 "sgbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 198 "sgbcon.f"
	*info = -1;
#line 199 "sgbcon.f"
    } else if (*n < 0) {
#line 200 "sgbcon.f"
	*info = -2;
#line 201 "sgbcon.f"
    } else if (*kl < 0) {
#line 202 "sgbcon.f"
	*info = -3;
#line 203 "sgbcon.f"
    } else if (*ku < 0) {
#line 204 "sgbcon.f"
	*info = -4;
#line 205 "sgbcon.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 206 "sgbcon.f"
	*info = -6;
#line 207 "sgbcon.f"
    } else if (*anorm < 0.) {
#line 208 "sgbcon.f"
	*info = -8;
#line 209 "sgbcon.f"
    }
#line 210 "sgbcon.f"
    if (*info != 0) {
#line 211 "sgbcon.f"
	i__1 = -(*info);
#line 211 "sgbcon.f"
	xerbla_("SGBCON", &i__1, (ftnlen)6);
#line 212 "sgbcon.f"
	return 0;
#line 213 "sgbcon.f"
    }

/*     Quick return if possible */

#line 217 "sgbcon.f"
    *rcond = 0.;
#line 218 "sgbcon.f"
    if (*n == 0) {
#line 219 "sgbcon.f"
	*rcond = 1.;
#line 220 "sgbcon.f"
	return 0;
#line 221 "sgbcon.f"
    } else if (*anorm == 0.) {
#line 222 "sgbcon.f"
	return 0;
#line 223 "sgbcon.f"
    }

#line 225 "sgbcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 229 "sgbcon.f"
    ainvnm = 0.;
#line 230 "sgbcon.f"
    *(unsigned char *)normin = 'N';
#line 231 "sgbcon.f"
    if (onenrm) {
#line 232 "sgbcon.f"
	kase1 = 1;
#line 233 "sgbcon.f"
    } else {
#line 234 "sgbcon.f"
	kase1 = 2;
#line 235 "sgbcon.f"
    }
#line 236 "sgbcon.f"
    kd = *kl + *ku + 1;
#line 237 "sgbcon.f"
    lnoti = *kl > 0;
#line 238 "sgbcon.f"
    kase = 0;
#line 239 "sgbcon.f"
L10:
#line 240 "sgbcon.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 241 "sgbcon.f"
    if (kase != 0) {
#line 242 "sgbcon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 246 "sgbcon.f"
	    if (lnoti) {
#line 247 "sgbcon.f"
		i__1 = *n - 1;
#line 247 "sgbcon.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 248 "sgbcon.f"
		    i__2 = *kl, i__3 = *n - j;
#line 248 "sgbcon.f"
		    lm = min(i__2,i__3);
#line 249 "sgbcon.f"
		    jp = ipiv[j];
#line 250 "sgbcon.f"
		    t = work[jp];
#line 251 "sgbcon.f"
		    if (jp != j) {
#line 252 "sgbcon.f"
			work[jp] = work[j];
#line 253 "sgbcon.f"
			work[j] = t;
#line 254 "sgbcon.f"
		    }
#line 255 "sgbcon.f"
		    d__1 = -t;
#line 255 "sgbcon.f"
		    saxpy_(&lm, &d__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 256 "sgbcon.f"
/* L20: */
#line 256 "sgbcon.f"
		}
#line 257 "sgbcon.f"
	    }

/*           Multiply by inv(U). */

#line 261 "sgbcon.f"
	    i__1 = *kl + *ku;
#line 261 "sgbcon.f"
	    slatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
		    ab[ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 
		    1], info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 264 "sgbcon.f"
	} else {

/*           Multiply by inv(U**T). */

#line 268 "sgbcon.f"
	    i__1 = *kl + *ku;
#line 268 "sgbcon.f"
	    slatbs_("Upper", "Transpose", "Non-unit", normin, n, &i__1, &ab[
		    ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**T). */

#line 274 "sgbcon.f"
	    if (lnoti) {
#line 275 "sgbcon.f"
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 276 "sgbcon.f"
		    i__1 = *kl, i__2 = *n - j;
#line 276 "sgbcon.f"
		    lm = min(i__1,i__2);
#line 277 "sgbcon.f"
		    work[j] -= sdot_(&lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 279 "sgbcon.f"
		    jp = ipiv[j];
#line 280 "sgbcon.f"
		    if (jp != j) {
#line 281 "sgbcon.f"
			t = work[jp];
#line 282 "sgbcon.f"
			work[jp] = work[j];
#line 283 "sgbcon.f"
			work[j] = t;
#line 284 "sgbcon.f"
		    }
#line 285 "sgbcon.f"
/* L30: */
#line 285 "sgbcon.f"
		}
#line 286 "sgbcon.f"
	    }
#line 287 "sgbcon.f"
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

#line 291 "sgbcon.f"
	*(unsigned char *)normin = 'Y';
#line 292 "sgbcon.f"
	if (scale != 1.) {
#line 293 "sgbcon.f"
	    ix = isamax_(n, &work[1], &c__1);
#line 294 "sgbcon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 294 "sgbcon.f"
		goto L40;
#line 294 "sgbcon.f"
	    }
#line 296 "sgbcon.f"
	    srscl_(n, &scale, &work[1], &c__1);
#line 297 "sgbcon.f"
	}
#line 298 "sgbcon.f"
	goto L10;
#line 299 "sgbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 303 "sgbcon.f"
    if (ainvnm != 0.) {
#line 303 "sgbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 303 "sgbcon.f"
    }

#line 306 "sgbcon.f"
L40:
#line 307 "sgbcon.f"
    return 0;

/*     End of SGBCON */

} /* sgbcon_ */

