#line 1 "stbcon.f"
/* stbcon.f -- translated by f2c (version 20100827).
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

#line 1 "stbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b STBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STBCON estimates the reciprocal of the condition number of a */
/* > triangular band matrix A, in either the 1-norm or the infinity-norm. */
/* > */
/* > The norm of A is computed and an estimate is obtained for */
/* > norm(inv(A)), then the reciprocal of the condition number is */
/* > computed as */
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
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals or subdiagonals of the */
/* >          triangular band matrix A.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          The upper or lower triangular band matrix A, stored in the */
/* >          first kd+1 rows of the array. The j-th column of A is stored */
/* >          in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* >          If DIAG = 'U', the diagonal elements of A are not referenced */
/* >          and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int stbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *rcond, 
	doublereal *work, integer *iwork, integer *info, ftnlen norm_len, 
	ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer ix, kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    extern /* Subroutine */ int srscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    static doublereal xnorm;
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern doublereal slantb_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int slatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static logical onenrm;
    static char normin[1];
    static doublereal smlnum;
    static logical nounit;


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 192 "stbcon.f"
    /* Parameter adjustments */
#line 192 "stbcon.f"
    ab_dim1 = *ldab;
#line 192 "stbcon.f"
    ab_offset = 1 + ab_dim1;
#line 192 "stbcon.f"
    ab -= ab_offset;
#line 192 "stbcon.f"
    --work;
#line 192 "stbcon.f"
    --iwork;
#line 192 "stbcon.f"

#line 192 "stbcon.f"
    /* Function Body */
#line 192 "stbcon.f"
    *info = 0;
#line 193 "stbcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 194 "stbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 195 "stbcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 197 "stbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 198 "stbcon.f"
	*info = -1;
#line 199 "stbcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 200 "stbcon.f"
	*info = -2;
#line 201 "stbcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 202 "stbcon.f"
	*info = -3;
#line 203 "stbcon.f"
    } else if (*n < 0) {
#line 204 "stbcon.f"
	*info = -4;
#line 205 "stbcon.f"
    } else if (*kd < 0) {
#line 206 "stbcon.f"
	*info = -5;
#line 207 "stbcon.f"
    } else if (*ldab < *kd + 1) {
#line 208 "stbcon.f"
	*info = -7;
#line 209 "stbcon.f"
    }
#line 210 "stbcon.f"
    if (*info != 0) {
#line 211 "stbcon.f"
	i__1 = -(*info);
#line 211 "stbcon.f"
	xerbla_("STBCON", &i__1, (ftnlen)6);
#line 212 "stbcon.f"
	return 0;
#line 213 "stbcon.f"
    }

/*     Quick return if possible */

#line 217 "stbcon.f"
    if (*n == 0) {
#line 218 "stbcon.f"
	*rcond = 1.;
#line 219 "stbcon.f"
	return 0;
#line 220 "stbcon.f"
    }

#line 222 "stbcon.f"
    *rcond = 0.;
#line 223 "stbcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 227 "stbcon.f"
    anorm = slantb_(norm, uplo, diag, n, kd, &ab[ab_offset], ldab, &work[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 231 "stbcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 235 "stbcon.f"
	ainvnm = 0.;
#line 236 "stbcon.f"
	*(unsigned char *)normin = 'N';
#line 237 "stbcon.f"
	if (onenrm) {
#line 238 "stbcon.f"
	    kase1 = 1;
#line 239 "stbcon.f"
	} else {
#line 240 "stbcon.f"
	    kase1 = 2;
#line 241 "stbcon.f"
	}
#line 242 "stbcon.f"
	kase = 0;
#line 243 "stbcon.f"
L10:
#line 244 "stbcon.f"
	slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 245 "stbcon.f"
	if (kase != 0) {
#line 246 "stbcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 250 "stbcon.f"
		slatbs_(uplo, "No transpose", diag, normin, n, kd, &ab[
			ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 
			1], info, (ftnlen)1, (ftnlen)12, (ftnlen)1, (ftnlen)1)
			;
#line 252 "stbcon.f"
	    } else {

/*              Multiply by inv(A**T). */

#line 256 "stbcon.f"
		slatbs_(uplo, "Transpose", diag, normin, n, kd, &ab[ab_offset]
			, ldab, &work[1], &scale, &work[(*n << 1) + 1], info, 
			(ftnlen)1, (ftnlen)9, (ftnlen)1, (ftnlen)1);
#line 258 "stbcon.f"
	    }
#line 259 "stbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 263 "stbcon.f"
	    if (scale != 1.) {
#line 264 "stbcon.f"
		ix = isamax_(n, &work[1], &c__1);
#line 265 "stbcon.f"
		xnorm = (d__1 = work[ix], abs(d__1));
#line 266 "stbcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 266 "stbcon.f"
		    goto L20;
#line 266 "stbcon.f"
		}
#line 268 "stbcon.f"
		srscl_(n, &scale, &work[1], &c__1);
#line 269 "stbcon.f"
	    }
#line 270 "stbcon.f"
	    goto L10;
#line 271 "stbcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 275 "stbcon.f"
	if (ainvnm != 0.) {
#line 275 "stbcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 275 "stbcon.f"
	}
#line 277 "stbcon.f"
    }

#line 279 "stbcon.f"
L20:
#line 280 "stbcon.f"
    return 0;

/*     End of STBCON */

} /* stbcon_ */

