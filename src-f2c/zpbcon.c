#line 1 "zpbcon.f"
/* zpbcon.f -- translated by f2c (version 20100827).
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

#line 1 "zpbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZPBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPBCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite band matrix using */
/* > the Cholesky factorization A = U**H*U or A = L*L**H computed by */
/* > ZPBTRF. */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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
/* >          ANORM is DOUBLE PRECISION */
/* >          The 1-norm (or infinity-norm) of the Hermitian band matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* >          estimate of the 1-norm of inv(A) computed in this routine. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpbcon_(char *uplo, integer *n, integer *kd, 
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
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalel, scaleu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern integer izamax_(integer *, doublecomplex *, integer *);
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

#line 189 "zpbcon.f"
    /* Parameter adjustments */
#line 189 "zpbcon.f"
    ab_dim1 = *ldab;
#line 189 "zpbcon.f"
    ab_offset = 1 + ab_dim1;
#line 189 "zpbcon.f"
    ab -= ab_offset;
#line 189 "zpbcon.f"
    --work;
#line 189 "zpbcon.f"
    --rwork;
#line 189 "zpbcon.f"

#line 189 "zpbcon.f"
    /* Function Body */
#line 189 "zpbcon.f"
    *info = 0;
#line 190 "zpbcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 191 "zpbcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 192 "zpbcon.f"
	*info = -1;
#line 193 "zpbcon.f"
    } else if (*n < 0) {
#line 194 "zpbcon.f"
	*info = -2;
#line 195 "zpbcon.f"
    } else if (*kd < 0) {
#line 196 "zpbcon.f"
	*info = -3;
#line 197 "zpbcon.f"
    } else if (*ldab < *kd + 1) {
#line 198 "zpbcon.f"
	*info = -5;
#line 199 "zpbcon.f"
    } else if (*anorm < 0.) {
#line 200 "zpbcon.f"
	*info = -6;
#line 201 "zpbcon.f"
    }
#line 202 "zpbcon.f"
    if (*info != 0) {
#line 203 "zpbcon.f"
	i__1 = -(*info);
#line 203 "zpbcon.f"
	xerbla_("ZPBCON", &i__1, (ftnlen)6);
#line 204 "zpbcon.f"
	return 0;
#line 205 "zpbcon.f"
    }

/*     Quick return if possible */

#line 209 "zpbcon.f"
    *rcond = 0.;
#line 210 "zpbcon.f"
    if (*n == 0) {
#line 211 "zpbcon.f"
	*rcond = 1.;
#line 212 "zpbcon.f"
	return 0;
#line 213 "zpbcon.f"
    } else if (*anorm == 0.) {
#line 214 "zpbcon.f"
	return 0;
#line 215 "zpbcon.f"
    }

#line 217 "zpbcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 221 "zpbcon.f"
    kase = 0;
#line 222 "zpbcon.f"
    *(unsigned char *)normin = 'N';
#line 223 "zpbcon.f"
L10:
#line 224 "zpbcon.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 225 "zpbcon.f"
    if (kase != 0) {
#line 226 "zpbcon.f"
	if (upper) {

/*           Multiply by inv(U**H). */

#line 230 "zpbcon.f"
	    zlatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, kd,
		     &ab[ab_offset], ldab, &work[1], &scalel, &rwork[1], info,
		     (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 233 "zpbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 237 "zpbcon.f"
	    zlatbs_("Upper", "No transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scaleu, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 239 "zpbcon.f"
	} else {

/*           Multiply by inv(L). */

#line 243 "zpbcon.f"
	    zlatbs_("Lower", "No transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scalel, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 245 "zpbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**H). */

#line 249 "zpbcon.f"
	    zlatbs_("Lower", "Conjugate transpose", "Non-unit", normin, n, kd,
		     &ab[ab_offset], ldab, &work[1], &scaleu, &rwork[1], info,
		     (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 252 "zpbcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 256 "zpbcon.f"
	scale = scalel * scaleu;
#line 257 "zpbcon.f"
	if (scale != 1.) {
#line 258 "zpbcon.f"
	    ix = izamax_(n, &work[1], &c__1);
#line 259 "zpbcon.f"
	    i__1 = ix;
#line 259 "zpbcon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 259 "zpbcon.f"
		goto L20;
#line 259 "zpbcon.f"
	    }
#line 261 "zpbcon.f"
	    zdrscl_(n, &scale, &work[1], &c__1);
#line 262 "zpbcon.f"
	}
#line 263 "zpbcon.f"
	goto L10;
#line 264 "zpbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 268 "zpbcon.f"
    if (ainvnm != 0.) {
#line 268 "zpbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 268 "zpbcon.f"
    }

#line 271 "zpbcon.f"
L20:

#line 273 "zpbcon.f"
    return 0;

/*     End of ZPBCON */

} /* zpbcon_ */

