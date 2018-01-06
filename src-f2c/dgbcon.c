#line 1 "dgbcon.f"
/* dgbcon.f -- translated by f2c (version 20100827).
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

#line 1 "dgbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, KL, KU, LDAB, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBCON estimates the reciprocal of the condition number of a real */
/* > general band matrix A, in either the 1-norm or the infinity-norm, */
/* > using the LU factorization computed by DGBTRF. */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by DGBTRF.  U is stored as an upper triangular band */
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
/* >          ANORM is DOUBLE PRECISION */
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
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

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgbcon_(char *norm, integer *n, integer *kl, integer *ku,
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
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lnoti;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dlacn2_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
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

#line 195 "dgbcon.f"
    /* Parameter adjustments */
#line 195 "dgbcon.f"
    ab_dim1 = *ldab;
#line 195 "dgbcon.f"
    ab_offset = 1 + ab_dim1;
#line 195 "dgbcon.f"
    ab -= ab_offset;
#line 195 "dgbcon.f"
    --ipiv;
#line 195 "dgbcon.f"
    --work;
#line 195 "dgbcon.f"
    --iwork;
#line 195 "dgbcon.f"

#line 195 "dgbcon.f"
    /* Function Body */
#line 195 "dgbcon.f"
    *info = 0;
#line 196 "dgbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 197 "dgbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 198 "dgbcon.f"
	*info = -1;
#line 199 "dgbcon.f"
    } else if (*n < 0) {
#line 200 "dgbcon.f"
	*info = -2;
#line 201 "dgbcon.f"
    } else if (*kl < 0) {
#line 202 "dgbcon.f"
	*info = -3;
#line 203 "dgbcon.f"
    } else if (*ku < 0) {
#line 204 "dgbcon.f"
	*info = -4;
#line 205 "dgbcon.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 206 "dgbcon.f"
	*info = -6;
#line 207 "dgbcon.f"
    } else if (*anorm < 0.) {
#line 208 "dgbcon.f"
	*info = -8;
#line 209 "dgbcon.f"
    }
#line 210 "dgbcon.f"
    if (*info != 0) {
#line 211 "dgbcon.f"
	i__1 = -(*info);
#line 211 "dgbcon.f"
	xerbla_("DGBCON", &i__1, (ftnlen)6);
#line 212 "dgbcon.f"
	return 0;
#line 213 "dgbcon.f"
    }

/*     Quick return if possible */

#line 217 "dgbcon.f"
    *rcond = 0.;
#line 218 "dgbcon.f"
    if (*n == 0) {
#line 219 "dgbcon.f"
	*rcond = 1.;
#line 220 "dgbcon.f"
	return 0;
#line 221 "dgbcon.f"
    } else if (*anorm == 0.) {
#line 222 "dgbcon.f"
	return 0;
#line 223 "dgbcon.f"
    }

#line 225 "dgbcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 229 "dgbcon.f"
    ainvnm = 0.;
#line 230 "dgbcon.f"
    *(unsigned char *)normin = 'N';
#line 231 "dgbcon.f"
    if (onenrm) {
#line 232 "dgbcon.f"
	kase1 = 1;
#line 233 "dgbcon.f"
    } else {
#line 234 "dgbcon.f"
	kase1 = 2;
#line 235 "dgbcon.f"
    }
#line 236 "dgbcon.f"
    kd = *kl + *ku + 1;
#line 237 "dgbcon.f"
    lnoti = *kl > 0;
#line 238 "dgbcon.f"
    kase = 0;
#line 239 "dgbcon.f"
L10:
#line 240 "dgbcon.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 241 "dgbcon.f"
    if (kase != 0) {
#line 242 "dgbcon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 246 "dgbcon.f"
	    if (lnoti) {
#line 247 "dgbcon.f"
		i__1 = *n - 1;
#line 247 "dgbcon.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 248 "dgbcon.f"
		    i__2 = *kl, i__3 = *n - j;
#line 248 "dgbcon.f"
		    lm = min(i__2,i__3);
#line 249 "dgbcon.f"
		    jp = ipiv[j];
#line 250 "dgbcon.f"
		    t = work[jp];
#line 251 "dgbcon.f"
		    if (jp != j) {
#line 252 "dgbcon.f"
			work[jp] = work[j];
#line 253 "dgbcon.f"
			work[j] = t;
#line 254 "dgbcon.f"
		    }
#line 255 "dgbcon.f"
		    d__1 = -t;
#line 255 "dgbcon.f"
		    daxpy_(&lm, &d__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 256 "dgbcon.f"
/* L20: */
#line 256 "dgbcon.f"
		}
#line 257 "dgbcon.f"
	    }

/*           Multiply by inv(U). */

#line 261 "dgbcon.f"
	    i__1 = *kl + *ku;
#line 261 "dgbcon.f"
	    dlatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
		    ab[ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 
		    1], info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 264 "dgbcon.f"
	} else {

/*           Multiply by inv(U**T). */

#line 268 "dgbcon.f"
	    i__1 = *kl + *ku;
#line 268 "dgbcon.f"
	    dlatbs_("Upper", "Transpose", "Non-unit", normin, n, &i__1, &ab[
		    ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**T). */

#line 274 "dgbcon.f"
	    if (lnoti) {
#line 275 "dgbcon.f"
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 276 "dgbcon.f"
		    i__1 = *kl, i__2 = *n - j;
#line 276 "dgbcon.f"
		    lm = min(i__1,i__2);
#line 277 "dgbcon.f"
		    work[j] -= ddot_(&lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 279 "dgbcon.f"
		    jp = ipiv[j];
#line 280 "dgbcon.f"
		    if (jp != j) {
#line 281 "dgbcon.f"
			t = work[jp];
#line 282 "dgbcon.f"
			work[jp] = work[j];
#line 283 "dgbcon.f"
			work[j] = t;
#line 284 "dgbcon.f"
		    }
#line 285 "dgbcon.f"
/* L30: */
#line 285 "dgbcon.f"
		}
#line 286 "dgbcon.f"
	    }
#line 287 "dgbcon.f"
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

#line 291 "dgbcon.f"
	*(unsigned char *)normin = 'Y';
#line 292 "dgbcon.f"
	if (scale != 1.) {
#line 293 "dgbcon.f"
	    ix = idamax_(n, &work[1], &c__1);
#line 294 "dgbcon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 294 "dgbcon.f"
		goto L40;
#line 294 "dgbcon.f"
	    }
#line 296 "dgbcon.f"
	    drscl_(n, &scale, &work[1], &c__1);
#line 297 "dgbcon.f"
	}
#line 298 "dgbcon.f"
	goto L10;
#line 299 "dgbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 303 "dgbcon.f"
    if (ainvnm != 0.) {
#line 303 "dgbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 303 "dgbcon.f"
    }

#line 306 "dgbcon.f"
L40:
#line 307 "dgbcon.f"
    return 0;

/*     End of DGBCON */

} /* dgbcon_ */

