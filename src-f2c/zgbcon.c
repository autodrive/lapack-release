#line 1 "zgbcon.f"
/* zgbcon.f -- translated by f2c (version 20100827).
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

#line 1 "zgbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZGBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, KL, KU, LDAB, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBCON estimates the reciprocal of the condition number of a complex */
/* > general band matrix A, in either the 1-norm or the infinity-norm, */
/* > using the LU factorization computed by ZGBTRF. */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by ZGBTRF.  U is stored as an upper triangular band */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup complex16GBcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer j;
    static doublecomplex t;
    static integer kd, lm, jp, ix, kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical lnoti;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zlacn2_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical onenrm;
    extern /* Subroutine */ int zlatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), zdrscl_(integer *, doublereal *, doublecomplex *, 
	    integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 205 "zgbcon.f"
    /* Parameter adjustments */
#line 205 "zgbcon.f"
    ab_dim1 = *ldab;
#line 205 "zgbcon.f"
    ab_offset = 1 + ab_dim1;
#line 205 "zgbcon.f"
    ab -= ab_offset;
#line 205 "zgbcon.f"
    --ipiv;
#line 205 "zgbcon.f"
    --work;
#line 205 "zgbcon.f"
    --rwork;
#line 205 "zgbcon.f"

#line 205 "zgbcon.f"
    /* Function Body */
#line 205 "zgbcon.f"
    *info = 0;
#line 206 "zgbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 207 "zgbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 208 "zgbcon.f"
	*info = -1;
#line 209 "zgbcon.f"
    } else if (*n < 0) {
#line 210 "zgbcon.f"
	*info = -2;
#line 211 "zgbcon.f"
    } else if (*kl < 0) {
#line 212 "zgbcon.f"
	*info = -3;
#line 213 "zgbcon.f"
    } else if (*ku < 0) {
#line 214 "zgbcon.f"
	*info = -4;
#line 215 "zgbcon.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 216 "zgbcon.f"
	*info = -6;
#line 217 "zgbcon.f"
    } else if (*anorm < 0.) {
#line 218 "zgbcon.f"
	*info = -8;
#line 219 "zgbcon.f"
    }
#line 220 "zgbcon.f"
    if (*info != 0) {
#line 221 "zgbcon.f"
	i__1 = -(*info);
#line 221 "zgbcon.f"
	xerbla_("ZGBCON", &i__1, (ftnlen)6);
#line 222 "zgbcon.f"
	return 0;
#line 223 "zgbcon.f"
    }

/*     Quick return if possible */

#line 227 "zgbcon.f"
    *rcond = 0.;
#line 228 "zgbcon.f"
    if (*n == 0) {
#line 229 "zgbcon.f"
	*rcond = 1.;
#line 230 "zgbcon.f"
	return 0;
#line 231 "zgbcon.f"
    } else if (*anorm == 0.) {
#line 232 "zgbcon.f"
	return 0;
#line 233 "zgbcon.f"
    }

#line 235 "zgbcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 239 "zgbcon.f"
    ainvnm = 0.;
#line 240 "zgbcon.f"
    *(unsigned char *)normin = 'N';
#line 241 "zgbcon.f"
    if (onenrm) {
#line 242 "zgbcon.f"
	kase1 = 1;
#line 243 "zgbcon.f"
    } else {
#line 244 "zgbcon.f"
	kase1 = 2;
#line 245 "zgbcon.f"
    }
#line 246 "zgbcon.f"
    kd = *kl + *ku + 1;
#line 247 "zgbcon.f"
    lnoti = *kl > 0;
#line 248 "zgbcon.f"
    kase = 0;
#line 249 "zgbcon.f"
L10:
#line 250 "zgbcon.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 251 "zgbcon.f"
    if (kase != 0) {
#line 252 "zgbcon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 256 "zgbcon.f"
	    if (lnoti) {
#line 257 "zgbcon.f"
		i__1 = *n - 1;
#line 257 "zgbcon.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 258 "zgbcon.f"
		    i__2 = *kl, i__3 = *n - j;
#line 258 "zgbcon.f"
		    lm = min(i__2,i__3);
#line 259 "zgbcon.f"
		    jp = ipiv[j];
#line 260 "zgbcon.f"
		    i__2 = jp;
#line 260 "zgbcon.f"
		    t.r = work[i__2].r, t.i = work[i__2].i;
#line 261 "zgbcon.f"
		    if (jp != j) {
#line 262 "zgbcon.f"
			i__2 = jp;
#line 262 "zgbcon.f"
			i__3 = j;
#line 262 "zgbcon.f"
			work[i__2].r = work[i__3].r, work[i__2].i = work[i__3]
				.i;
#line 263 "zgbcon.f"
			i__2 = j;
#line 263 "zgbcon.f"
			work[i__2].r = t.r, work[i__2].i = t.i;
#line 264 "zgbcon.f"
		    }
#line 265 "zgbcon.f"
		    z__1.r = -t.r, z__1.i = -t.i;
#line 265 "zgbcon.f"
		    zaxpy_(&lm, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 266 "zgbcon.f"
/* L20: */
#line 266 "zgbcon.f"
		}
#line 267 "zgbcon.f"
	    }

/*           Multiply by inv(U). */

#line 271 "zgbcon.f"
	    i__1 = *kl + *ku;
#line 271 "zgbcon.f"
	    zlatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
		    ab[ab_offset], ldab, &work[1], &scale, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 273 "zgbcon.f"
	} else {

/*           Multiply by inv(U**H). */

#line 277 "zgbcon.f"
	    i__1 = *kl + *ku;
#line 277 "zgbcon.f"
	    zlatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    i__1, &ab[ab_offset], ldab, &work[1], &scale, &rwork[1], 
		    info, (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**H). */

#line 283 "zgbcon.f"
	    if (lnoti) {
#line 284 "zgbcon.f"
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 285 "zgbcon.f"
		    i__1 = *kl, i__2 = *n - j;
#line 285 "zgbcon.f"
		    lm = min(i__1,i__2);
#line 286 "zgbcon.f"
		    i__1 = j;
#line 286 "zgbcon.f"
		    i__2 = j;
#line 286 "zgbcon.f"
		    zdotc_(&z__2, &lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 286 "zgbcon.f"
		    z__1.r = work[i__2].r - z__2.r, z__1.i = work[i__2].i - 
			    z__2.i;
#line 286 "zgbcon.f"
		    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 288 "zgbcon.f"
		    jp = ipiv[j];
#line 289 "zgbcon.f"
		    if (jp != j) {
#line 290 "zgbcon.f"
			i__1 = jp;
#line 290 "zgbcon.f"
			t.r = work[i__1].r, t.i = work[i__1].i;
#line 291 "zgbcon.f"
			i__1 = jp;
#line 291 "zgbcon.f"
			i__2 = j;
#line 291 "zgbcon.f"
			work[i__1].r = work[i__2].r, work[i__1].i = work[i__2]
				.i;
#line 292 "zgbcon.f"
			i__1 = j;
#line 292 "zgbcon.f"
			work[i__1].r = t.r, work[i__1].i = t.i;
#line 293 "zgbcon.f"
		    }
#line 294 "zgbcon.f"
/* L30: */
#line 294 "zgbcon.f"
		}
#line 295 "zgbcon.f"
	    }
#line 296 "zgbcon.f"
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

#line 300 "zgbcon.f"
	*(unsigned char *)normin = 'Y';
#line 301 "zgbcon.f"
	if (scale != 1.) {
#line 302 "zgbcon.f"
	    ix = izamax_(n, &work[1], &c__1);
#line 303 "zgbcon.f"
	    i__1 = ix;
#line 303 "zgbcon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 303 "zgbcon.f"
		goto L40;
#line 303 "zgbcon.f"
	    }
#line 305 "zgbcon.f"
	    zdrscl_(n, &scale, &work[1], &c__1);
#line 306 "zgbcon.f"
	}
#line 307 "zgbcon.f"
	goto L10;
#line 308 "zgbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 312 "zgbcon.f"
    if (ainvnm != 0.) {
#line 312 "zgbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 312 "zgbcon.f"
    }

#line 315 "zgbcon.f"
L40:
#line 316 "zgbcon.f"
    return 0;

/*     End of ZGBCON */

} /* zgbcon_ */

