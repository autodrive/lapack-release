#line 1 "cgbcon.f"
/* cgbcon.f -- translated by f2c (version 20100827).
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

#line 1 "cgbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CGBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, KL, KU, LDAB, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBCON estimates the reciprocal of the condition number of a complex */
/* > general band matrix A, in either the 1-norm or the infinity-norm, */
/* > using the LU factorization computed by CGBTRF. */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by CGBTRF.  U is stored as an upper triangular band */
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
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
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

/* > \date November 2011 */

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgbcon_(char *norm, integer *n, integer *kl, integer *ku,
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
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical lnoti;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static logical onenrm;
    static char normin[1];
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 205 "cgbcon.f"
    /* Parameter adjustments */
#line 205 "cgbcon.f"
    ab_dim1 = *ldab;
#line 205 "cgbcon.f"
    ab_offset = 1 + ab_dim1;
#line 205 "cgbcon.f"
    ab -= ab_offset;
#line 205 "cgbcon.f"
    --ipiv;
#line 205 "cgbcon.f"
    --work;
#line 205 "cgbcon.f"
    --rwork;
#line 205 "cgbcon.f"

#line 205 "cgbcon.f"
    /* Function Body */
#line 205 "cgbcon.f"
    *info = 0;
#line 206 "cgbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 207 "cgbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 208 "cgbcon.f"
	*info = -1;
#line 209 "cgbcon.f"
    } else if (*n < 0) {
#line 210 "cgbcon.f"
	*info = -2;
#line 211 "cgbcon.f"
    } else if (*kl < 0) {
#line 212 "cgbcon.f"
	*info = -3;
#line 213 "cgbcon.f"
    } else if (*ku < 0) {
#line 214 "cgbcon.f"
	*info = -4;
#line 215 "cgbcon.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 216 "cgbcon.f"
	*info = -6;
#line 217 "cgbcon.f"
    } else if (*anorm < 0.) {
#line 218 "cgbcon.f"
	*info = -8;
#line 219 "cgbcon.f"
    }
#line 220 "cgbcon.f"
    if (*info != 0) {
#line 221 "cgbcon.f"
	i__1 = -(*info);
#line 221 "cgbcon.f"
	xerbla_("CGBCON", &i__1, (ftnlen)6);
#line 222 "cgbcon.f"
	return 0;
#line 223 "cgbcon.f"
    }

/*     Quick return if possible */

#line 227 "cgbcon.f"
    *rcond = 0.;
#line 228 "cgbcon.f"
    if (*n == 0) {
#line 229 "cgbcon.f"
	*rcond = 1.;
#line 230 "cgbcon.f"
	return 0;
#line 231 "cgbcon.f"
    } else if (*anorm == 0.) {
#line 232 "cgbcon.f"
	return 0;
#line 233 "cgbcon.f"
    }

#line 235 "cgbcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 239 "cgbcon.f"
    ainvnm = 0.;
#line 240 "cgbcon.f"
    *(unsigned char *)normin = 'N';
#line 241 "cgbcon.f"
    if (onenrm) {
#line 242 "cgbcon.f"
	kase1 = 1;
#line 243 "cgbcon.f"
    } else {
#line 244 "cgbcon.f"
	kase1 = 2;
#line 245 "cgbcon.f"
    }
#line 246 "cgbcon.f"
    kd = *kl + *ku + 1;
#line 247 "cgbcon.f"
    lnoti = *kl > 0;
#line 248 "cgbcon.f"
    kase = 0;
#line 249 "cgbcon.f"
L10:
#line 250 "cgbcon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 251 "cgbcon.f"
    if (kase != 0) {
#line 252 "cgbcon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 256 "cgbcon.f"
	    if (lnoti) {
#line 257 "cgbcon.f"
		i__1 = *n - 1;
#line 257 "cgbcon.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 258 "cgbcon.f"
		    i__2 = *kl, i__3 = *n - j;
#line 258 "cgbcon.f"
		    lm = min(i__2,i__3);
#line 259 "cgbcon.f"
		    jp = ipiv[j];
#line 260 "cgbcon.f"
		    i__2 = jp;
#line 260 "cgbcon.f"
		    t.r = work[i__2].r, t.i = work[i__2].i;
#line 261 "cgbcon.f"
		    if (jp != j) {
#line 262 "cgbcon.f"
			i__2 = jp;
#line 262 "cgbcon.f"
			i__3 = j;
#line 262 "cgbcon.f"
			work[i__2].r = work[i__3].r, work[i__2].i = work[i__3]
				.i;
#line 263 "cgbcon.f"
			i__2 = j;
#line 263 "cgbcon.f"
			work[i__2].r = t.r, work[i__2].i = t.i;
#line 264 "cgbcon.f"
		    }
#line 265 "cgbcon.f"
		    z__1.r = -t.r, z__1.i = -t.i;
#line 265 "cgbcon.f"
		    caxpy_(&lm, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 266 "cgbcon.f"
/* L20: */
#line 266 "cgbcon.f"
		}
#line 267 "cgbcon.f"
	    }

/*           Multiply by inv(U). */

#line 271 "cgbcon.f"
	    i__1 = *kl + *ku;
#line 271 "cgbcon.f"
	    clatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
		    ab[ab_offset], ldab, &work[1], &scale, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 273 "cgbcon.f"
	} else {

/*           Multiply by inv(U**H). */

#line 277 "cgbcon.f"
	    i__1 = *kl + *ku;
#line 277 "cgbcon.f"
	    clatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    i__1, &ab[ab_offset], ldab, &work[1], &scale, &rwork[1], 
		    info, (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**H). */

#line 283 "cgbcon.f"
	    if (lnoti) {
#line 284 "cgbcon.f"
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 285 "cgbcon.f"
		    i__1 = *kl, i__2 = *n - j;
#line 285 "cgbcon.f"
		    lm = min(i__1,i__2);
#line 286 "cgbcon.f"
		    i__1 = j;
#line 286 "cgbcon.f"
		    i__2 = j;
#line 286 "cgbcon.f"
		    cdotc_(&z__2, &lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
#line 286 "cgbcon.f"
		    z__1.r = work[i__2].r - z__2.r, z__1.i = work[i__2].i - 
			    z__2.i;
#line 286 "cgbcon.f"
		    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 288 "cgbcon.f"
		    jp = ipiv[j];
#line 289 "cgbcon.f"
		    if (jp != j) {
#line 290 "cgbcon.f"
			i__1 = jp;
#line 290 "cgbcon.f"
			t.r = work[i__1].r, t.i = work[i__1].i;
#line 291 "cgbcon.f"
			i__1 = jp;
#line 291 "cgbcon.f"
			i__2 = j;
#line 291 "cgbcon.f"
			work[i__1].r = work[i__2].r, work[i__1].i = work[i__2]
				.i;
#line 292 "cgbcon.f"
			i__1 = j;
#line 292 "cgbcon.f"
			work[i__1].r = t.r, work[i__1].i = t.i;
#line 293 "cgbcon.f"
		    }
#line 294 "cgbcon.f"
/* L30: */
#line 294 "cgbcon.f"
		}
#line 295 "cgbcon.f"
	    }
#line 296 "cgbcon.f"
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

#line 300 "cgbcon.f"
	*(unsigned char *)normin = 'Y';
#line 301 "cgbcon.f"
	if (scale != 1.) {
#line 302 "cgbcon.f"
	    ix = icamax_(n, &work[1], &c__1);
#line 303 "cgbcon.f"
	    i__1 = ix;
#line 303 "cgbcon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 303 "cgbcon.f"
		goto L40;
#line 303 "cgbcon.f"
	    }
#line 305 "cgbcon.f"
	    csrscl_(n, &scale, &work[1], &c__1);
#line 306 "cgbcon.f"
	}
#line 307 "cgbcon.f"
	goto L10;
#line 308 "cgbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 312 "cgbcon.f"
    if (ainvnm != 0.) {
#line 312 "cgbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 312 "cgbcon.f"
    }

#line 315 "cgbcon.f"
L40:
#line 316 "cgbcon.f"
    return 0;

/*     End of CGBCON */

} /* cgbcon_ */

