#line 1 "dtbcon.f"
/* dtbcon.f -- translated by f2c (version 20100827).
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

#line 1 "dtbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DTBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTBCON estimates the reciprocal of the condition number of a */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtbcon_(char *norm, char *uplo, char *diag, integer *n, 
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
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm;
    static logical upper;
    static doublereal xnorm;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern doublereal dlantb_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int dlatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
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

#line 192 "dtbcon.f"
    /* Parameter adjustments */
#line 192 "dtbcon.f"
    ab_dim1 = *ldab;
#line 192 "dtbcon.f"
    ab_offset = 1 + ab_dim1;
#line 192 "dtbcon.f"
    ab -= ab_offset;
#line 192 "dtbcon.f"
    --work;
#line 192 "dtbcon.f"
    --iwork;
#line 192 "dtbcon.f"

#line 192 "dtbcon.f"
    /* Function Body */
#line 192 "dtbcon.f"
    *info = 0;
#line 193 "dtbcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 194 "dtbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 195 "dtbcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 197 "dtbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 198 "dtbcon.f"
	*info = -1;
#line 199 "dtbcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 200 "dtbcon.f"
	*info = -2;
#line 201 "dtbcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 202 "dtbcon.f"
	*info = -3;
#line 203 "dtbcon.f"
    } else if (*n < 0) {
#line 204 "dtbcon.f"
	*info = -4;
#line 205 "dtbcon.f"
    } else if (*kd < 0) {
#line 206 "dtbcon.f"
	*info = -5;
#line 207 "dtbcon.f"
    } else if (*ldab < *kd + 1) {
#line 208 "dtbcon.f"
	*info = -7;
#line 209 "dtbcon.f"
    }
#line 210 "dtbcon.f"
    if (*info != 0) {
#line 211 "dtbcon.f"
	i__1 = -(*info);
#line 211 "dtbcon.f"
	xerbla_("DTBCON", &i__1, (ftnlen)6);
#line 212 "dtbcon.f"
	return 0;
#line 213 "dtbcon.f"
    }

/*     Quick return if possible */

#line 217 "dtbcon.f"
    if (*n == 0) {
#line 218 "dtbcon.f"
	*rcond = 1.;
#line 219 "dtbcon.f"
	return 0;
#line 220 "dtbcon.f"
    }

#line 222 "dtbcon.f"
    *rcond = 0.;
#line 223 "dtbcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 227 "dtbcon.f"
    anorm = dlantb_(norm, uplo, diag, n, kd, &ab[ab_offset], ldab, &work[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 231 "dtbcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 235 "dtbcon.f"
	ainvnm = 0.;
#line 236 "dtbcon.f"
	*(unsigned char *)normin = 'N';
#line 237 "dtbcon.f"
	if (onenrm) {
#line 238 "dtbcon.f"
	    kase1 = 1;
#line 239 "dtbcon.f"
	} else {
#line 240 "dtbcon.f"
	    kase1 = 2;
#line 241 "dtbcon.f"
	}
#line 242 "dtbcon.f"
	kase = 0;
#line 243 "dtbcon.f"
L10:
#line 244 "dtbcon.f"
	dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 245 "dtbcon.f"
	if (kase != 0) {
#line 246 "dtbcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 250 "dtbcon.f"
		dlatbs_(uplo, "No transpose", diag, normin, n, kd, &ab[
			ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 
			1], info, (ftnlen)1, (ftnlen)12, (ftnlen)1, (ftnlen)1)
			;
#line 252 "dtbcon.f"
	    } else {

/*              Multiply by inv(A**T). */

#line 256 "dtbcon.f"
		dlatbs_(uplo, "Transpose", diag, normin, n, kd, &ab[ab_offset]
			, ldab, &work[1], &scale, &work[(*n << 1) + 1], info, 
			(ftnlen)1, (ftnlen)9, (ftnlen)1, (ftnlen)1);
#line 258 "dtbcon.f"
	    }
#line 259 "dtbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 263 "dtbcon.f"
	    if (scale != 1.) {
#line 264 "dtbcon.f"
		ix = idamax_(n, &work[1], &c__1);
#line 265 "dtbcon.f"
		xnorm = (d__1 = work[ix], abs(d__1));
#line 266 "dtbcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 266 "dtbcon.f"
		    goto L20;
#line 266 "dtbcon.f"
		}
#line 268 "dtbcon.f"
		drscl_(n, &scale, &work[1], &c__1);
#line 269 "dtbcon.f"
	    }
#line 270 "dtbcon.f"
	    goto L10;
#line 271 "dtbcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 275 "dtbcon.f"
	if (ainvnm != 0.) {
#line 275 "dtbcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 275 "dtbcon.f"
	}
#line 277 "dtbcon.f"
    }

#line 279 "dtbcon.f"
L20:
#line 280 "dtbcon.f"
    return 0;

/*     End of DTBCON */

} /* dtbcon_ */

