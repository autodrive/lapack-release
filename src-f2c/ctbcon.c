#line 1 "ctbcon.f"
/* ctbcon.f -- translated by f2c (version 20100827).
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

#line 1 "ctbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTBCON estimates the reciprocal of the condition number of a */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	norm_len, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer ix, kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    static doublereal xnorm;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal clantb_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen, ftnlen),
	     slamch_(char *, ftnlen);
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
    static logical nounit;


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

#line 199 "ctbcon.f"
    /* Parameter adjustments */
#line 199 "ctbcon.f"
    ab_dim1 = *ldab;
#line 199 "ctbcon.f"
    ab_offset = 1 + ab_dim1;
#line 199 "ctbcon.f"
    ab -= ab_offset;
#line 199 "ctbcon.f"
    --work;
#line 199 "ctbcon.f"
    --rwork;
#line 199 "ctbcon.f"

#line 199 "ctbcon.f"
    /* Function Body */
#line 199 "ctbcon.f"
    *info = 0;
#line 200 "ctbcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 201 "ctbcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 202 "ctbcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 204 "ctbcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 205 "ctbcon.f"
	*info = -1;
#line 206 "ctbcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 207 "ctbcon.f"
	*info = -2;
#line 208 "ctbcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 209 "ctbcon.f"
	*info = -3;
#line 210 "ctbcon.f"
    } else if (*n < 0) {
#line 211 "ctbcon.f"
	*info = -4;
#line 212 "ctbcon.f"
    } else if (*kd < 0) {
#line 213 "ctbcon.f"
	*info = -5;
#line 214 "ctbcon.f"
    } else if (*ldab < *kd + 1) {
#line 215 "ctbcon.f"
	*info = -7;
#line 216 "ctbcon.f"
    }
#line 217 "ctbcon.f"
    if (*info != 0) {
#line 218 "ctbcon.f"
	i__1 = -(*info);
#line 218 "ctbcon.f"
	xerbla_("CTBCON", &i__1, (ftnlen)6);
#line 219 "ctbcon.f"
	return 0;
#line 220 "ctbcon.f"
    }

/*     Quick return if possible */

#line 224 "ctbcon.f"
    if (*n == 0) {
#line 225 "ctbcon.f"
	*rcond = 1.;
#line 226 "ctbcon.f"
	return 0;
#line 227 "ctbcon.f"
    }

#line 229 "ctbcon.f"
    *rcond = 0.;
#line 230 "ctbcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) * (doublereal) max(*n,1);

/*     Compute the 1-norm of the triangular matrix A or A**H. */

#line 234 "ctbcon.f"
    anorm = clantb_(norm, uplo, diag, n, kd, &ab[ab_offset], ldab, &rwork[1], 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 238 "ctbcon.f"
    if (anorm > 0.) {

/*        Estimate the 1-norm of the inverse of A. */

#line 242 "ctbcon.f"
	ainvnm = 0.;
#line 243 "ctbcon.f"
	*(unsigned char *)normin = 'N';
#line 244 "ctbcon.f"
	if (onenrm) {
#line 245 "ctbcon.f"
	    kase1 = 1;
#line 246 "ctbcon.f"
	} else {
#line 247 "ctbcon.f"
	    kase1 = 2;
#line 248 "ctbcon.f"
	}
#line 249 "ctbcon.f"
	kase = 0;
#line 250 "ctbcon.f"
L10:
#line 251 "ctbcon.f"
	clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 252 "ctbcon.f"
	if (kase != 0) {
#line 253 "ctbcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 257 "ctbcon.f"
		clatbs_(uplo, "No transpose", diag, normin, n, kd, &ab[
			ab_offset], ldab, &work[1], &scale, &rwork[1], info, (
			ftnlen)1, (ftnlen)12, (ftnlen)1, (ftnlen)1);
#line 259 "ctbcon.f"
	    } else {

/*              Multiply by inv(A**H). */

#line 263 "ctbcon.f"
		clatbs_(uplo, "Conjugate transpose", diag, normin, n, kd, &ab[
			ab_offset], ldab, &work[1], &scale, &rwork[1], info, (
			ftnlen)1, (ftnlen)19, (ftnlen)1, (ftnlen)1);
#line 265 "ctbcon.f"
	    }
#line 266 "ctbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 270 "ctbcon.f"
	    if (scale != 1.) {
#line 271 "ctbcon.f"
		ix = icamax_(n, &work[1], &c__1);
#line 272 "ctbcon.f"
		i__1 = ix;
#line 272 "ctbcon.f"
		xnorm = (d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			work[ix]), abs(d__2));
#line 273 "ctbcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 273 "ctbcon.f"
		    goto L20;
#line 273 "ctbcon.f"
		}
#line 275 "ctbcon.f"
		csrscl_(n, &scale, &work[1], &c__1);
#line 276 "ctbcon.f"
	    }
#line 277 "ctbcon.f"
	    goto L10;
#line 278 "ctbcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 282 "ctbcon.f"
	if (ainvnm != 0.) {
#line 282 "ctbcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 282 "ctbcon.f"
	}
#line 284 "ctbcon.f"
    }

#line 286 "ctbcon.f"
L20:
#line 287 "ctbcon.f"
    return 0;

/*     End of CTBCON */

} /* ctbcon_ */

