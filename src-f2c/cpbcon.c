#line 1 "cpbcon.f"
/* cpbcon.f -- translated by f2c (version 20100827).
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

#line 1 "cpbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CPBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       REAL               ANORM, RCOND */
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
/* > CPBCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite band matrix using */
/* > the Cholesky factorization A = U**H*U or A = L*L**H computed by */
/* > CPBTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangular factor stored in AB; */
/* >          = 'L':  Lower triangular factor stored in AB. */
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
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H of the band matrix A, stored in the */
/* >          first KD+1 rows of the array.  The j-th column of U or L is */
/* >          stored in the j-th column of the array AB as follows: */
/* >          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          The 1-norm (or infinity-norm) of the Hermitian band matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* >          estimate of the 1-norm of inv(A) computed in this routine. */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpbcon_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *
	rcond, doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer ix, kase;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    static doublereal scalel;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static doublereal scaleu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
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

#line 189 "cpbcon.f"
    /* Parameter adjustments */
#line 189 "cpbcon.f"
    ab_dim1 = *ldab;
#line 189 "cpbcon.f"
    ab_offset = 1 + ab_dim1;
#line 189 "cpbcon.f"
    ab -= ab_offset;
#line 189 "cpbcon.f"
    --work;
#line 189 "cpbcon.f"
    --rwork;
#line 189 "cpbcon.f"

#line 189 "cpbcon.f"
    /* Function Body */
#line 189 "cpbcon.f"
    *info = 0;
#line 190 "cpbcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 191 "cpbcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 192 "cpbcon.f"
	*info = -1;
#line 193 "cpbcon.f"
    } else if (*n < 0) {
#line 194 "cpbcon.f"
	*info = -2;
#line 195 "cpbcon.f"
    } else if (*kd < 0) {
#line 196 "cpbcon.f"
	*info = -3;
#line 197 "cpbcon.f"
    } else if (*ldab < *kd + 1) {
#line 198 "cpbcon.f"
	*info = -5;
#line 199 "cpbcon.f"
    } else if (*anorm < 0.) {
#line 200 "cpbcon.f"
	*info = -6;
#line 201 "cpbcon.f"
    }
#line 202 "cpbcon.f"
    if (*info != 0) {
#line 203 "cpbcon.f"
	i__1 = -(*info);
#line 203 "cpbcon.f"
	xerbla_("CPBCON", &i__1, (ftnlen)6);
#line 204 "cpbcon.f"
	return 0;
#line 205 "cpbcon.f"
    }

/*     Quick return if possible */

#line 209 "cpbcon.f"
    *rcond = 0.;
#line 210 "cpbcon.f"
    if (*n == 0) {
#line 211 "cpbcon.f"
	*rcond = 1.;
#line 212 "cpbcon.f"
	return 0;
#line 213 "cpbcon.f"
    } else if (*anorm == 0.) {
#line 214 "cpbcon.f"
	return 0;
#line 215 "cpbcon.f"
    }

#line 217 "cpbcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 221 "cpbcon.f"
    kase = 0;
#line 222 "cpbcon.f"
    *(unsigned char *)normin = 'N';
#line 223 "cpbcon.f"
L10:
#line 224 "cpbcon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 225 "cpbcon.f"
    if (kase != 0) {
#line 226 "cpbcon.f"
	if (upper) {

/*           Multiply by inv(U**H). */

#line 230 "cpbcon.f"
	    clatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, kd,
		     &ab[ab_offset], ldab, &work[1], &scalel, &rwork[1], info,
		     (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 233 "cpbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 237 "cpbcon.f"
	    clatbs_("Upper", "No transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scaleu, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 239 "cpbcon.f"
	} else {

/*           Multiply by inv(L). */

#line 243 "cpbcon.f"
	    clatbs_("Lower", "No transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scalel, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 245 "cpbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**H). */

#line 249 "cpbcon.f"
	    clatbs_("Lower", "Conjugate transpose", "Non-unit", normin, n, kd,
		     &ab[ab_offset], ldab, &work[1], &scaleu, &rwork[1], info,
		     (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 252 "cpbcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 256 "cpbcon.f"
	scale = scalel * scaleu;
#line 257 "cpbcon.f"
	if (scale != 1.) {
#line 258 "cpbcon.f"
	    ix = icamax_(n, &work[1], &c__1);
#line 259 "cpbcon.f"
	    i__1 = ix;
#line 259 "cpbcon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 259 "cpbcon.f"
		goto L20;
#line 259 "cpbcon.f"
	    }
#line 261 "cpbcon.f"
	    csrscl_(n, &scale, &work[1], &c__1);
#line 262 "cpbcon.f"
	}
#line 263 "cpbcon.f"
	goto L10;
#line 264 "cpbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 268 "cpbcon.f"
    if (ainvnm != 0.) {
#line 268 "cpbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 268 "cpbcon.f"
    }

#line 271 "cpbcon.f"
L20:

#line 273 "cpbcon.f"
    return 0;

/*     End of CPBCON */

} /* cpbcon_ */

