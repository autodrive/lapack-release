#line 1 "dpbcon.f"
/* dpbcon.f -- translated by f2c (version 20100827).
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

#line 1 "dpbcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DPBCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
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
/* > DPBCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite band matrix using the */
/* > Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF. */
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
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T of the band matrix A, stored in the */
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
/* >          The 1-norm (or infinity-norm) of the symmetric band matrix A. */
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
/* Subroutine */ int dpbcon_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer ix, kase;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalel;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static doublereal scaleu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 181 "dpbcon.f"
    /* Parameter adjustments */
#line 181 "dpbcon.f"
    ab_dim1 = *ldab;
#line 181 "dpbcon.f"
    ab_offset = 1 + ab_dim1;
#line 181 "dpbcon.f"
    ab -= ab_offset;
#line 181 "dpbcon.f"
    --work;
#line 181 "dpbcon.f"
    --iwork;
#line 181 "dpbcon.f"

#line 181 "dpbcon.f"
    /* Function Body */
#line 181 "dpbcon.f"
    *info = 0;
#line 182 "dpbcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 183 "dpbcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 184 "dpbcon.f"
	*info = -1;
#line 185 "dpbcon.f"
    } else if (*n < 0) {
#line 186 "dpbcon.f"
	*info = -2;
#line 187 "dpbcon.f"
    } else if (*kd < 0) {
#line 188 "dpbcon.f"
	*info = -3;
#line 189 "dpbcon.f"
    } else if (*ldab < *kd + 1) {
#line 190 "dpbcon.f"
	*info = -5;
#line 191 "dpbcon.f"
    } else if (*anorm < 0.) {
#line 192 "dpbcon.f"
	*info = -6;
#line 193 "dpbcon.f"
    }
#line 194 "dpbcon.f"
    if (*info != 0) {
#line 195 "dpbcon.f"
	i__1 = -(*info);
#line 195 "dpbcon.f"
	xerbla_("DPBCON", &i__1, (ftnlen)6);
#line 196 "dpbcon.f"
	return 0;
#line 197 "dpbcon.f"
    }

/*     Quick return if possible */

#line 201 "dpbcon.f"
    *rcond = 0.;
#line 202 "dpbcon.f"
    if (*n == 0) {
#line 203 "dpbcon.f"
	*rcond = 1.;
#line 204 "dpbcon.f"
	return 0;
#line 205 "dpbcon.f"
    } else if (*anorm == 0.) {
#line 206 "dpbcon.f"
	return 0;
#line 207 "dpbcon.f"
    }

#line 209 "dpbcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 213 "dpbcon.f"
    kase = 0;
#line 214 "dpbcon.f"
    *(unsigned char *)normin = 'N';
#line 215 "dpbcon.f"
L10:
#line 216 "dpbcon.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 217 "dpbcon.f"
    if (kase != 0) {
#line 218 "dpbcon.f"
	if (upper) {

/*           Multiply by inv(U**T). */

#line 222 "dpbcon.f"
	    dlatbs_("Upper", "Transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scalel, &work[(*n << 1) + 1],
		     info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 225 "dpbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 229 "dpbcon.f"
	    dlatbs_("Upper", "No transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scaleu, &work[(*n << 1) + 1],
		     info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 232 "dpbcon.f"
	} else {

/*           Multiply by inv(L). */

#line 236 "dpbcon.f"
	    dlatbs_("Lower", "No transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scalel, &work[(*n << 1) + 1],
		     info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 239 "dpbcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**T). */

#line 243 "dpbcon.f"
	    dlatbs_("Lower", "Transpose", "Non-unit", normin, n, kd, &ab[
		    ab_offset], ldab, &work[1], &scaleu, &work[(*n << 1) + 1],
		     info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 246 "dpbcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 250 "dpbcon.f"
	scale = scalel * scaleu;
#line 251 "dpbcon.f"
	if (scale != 1.) {
#line 252 "dpbcon.f"
	    ix = idamax_(n, &work[1], &c__1);
#line 253 "dpbcon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 253 "dpbcon.f"
		goto L20;
#line 253 "dpbcon.f"
	    }
#line 255 "dpbcon.f"
	    drscl_(n, &scale, &work[1], &c__1);
#line 256 "dpbcon.f"
	}
#line 257 "dpbcon.f"
	goto L10;
#line 258 "dpbcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 262 "dpbcon.f"
    if (ainvnm != 0.) {
#line 262 "dpbcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 262 "dpbcon.f"
    }

#line 265 "dpbcon.f"
L20:

#line 267 "dpbcon.f"
    return 0;

/*     End of DPBCON */

} /* dpbcon_ */

