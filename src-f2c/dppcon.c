#line 1 "dppcon.f"
/* dppcon.f -- translated by f2c (version 20100827).
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

#line 1 "dppcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DPPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dppcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dppcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dppcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite packed matrix using */
/* > the Cholesky factorization A = U**T*U or A = L*L**T computed by */
/* > DPPTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, packed columnwise in a linear */
/* >          array.  The j-th column of U or L is stored in the array AP */
/* >          as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          The 1-norm (or infinity-norm) of the symmetric matrix A. */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dppcon_(char *uplo, integer *n, doublereal *ap, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1;
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
    static doublereal scaleu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlatps_(
	    char *, char *, char *, char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
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

#line 167 "dppcon.f"
    /* Parameter adjustments */
#line 167 "dppcon.f"
    --iwork;
#line 167 "dppcon.f"
    --work;
#line 167 "dppcon.f"
    --ap;
#line 167 "dppcon.f"

#line 167 "dppcon.f"
    /* Function Body */
#line 167 "dppcon.f"
    *info = 0;
#line 168 "dppcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 169 "dppcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 170 "dppcon.f"
	*info = -1;
#line 171 "dppcon.f"
    } else if (*n < 0) {
#line 172 "dppcon.f"
	*info = -2;
#line 173 "dppcon.f"
    } else if (*anorm < 0.) {
#line 174 "dppcon.f"
	*info = -4;
#line 175 "dppcon.f"
    }
#line 176 "dppcon.f"
    if (*info != 0) {
#line 177 "dppcon.f"
	i__1 = -(*info);
#line 177 "dppcon.f"
	xerbla_("DPPCON", &i__1, (ftnlen)6);
#line 178 "dppcon.f"
	return 0;
#line 179 "dppcon.f"
    }

/*     Quick return if possible */

#line 183 "dppcon.f"
    *rcond = 0.;
#line 184 "dppcon.f"
    if (*n == 0) {
#line 185 "dppcon.f"
	*rcond = 1.;
#line 186 "dppcon.f"
	return 0;
#line 187 "dppcon.f"
    } else if (*anorm == 0.) {
#line 188 "dppcon.f"
	return 0;
#line 189 "dppcon.f"
    }

#line 191 "dppcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 195 "dppcon.f"
    kase = 0;
#line 196 "dppcon.f"
    *(unsigned char *)normin = 'N';
#line 197 "dppcon.f"
L10:
#line 198 "dppcon.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 199 "dppcon.f"
    if (kase != 0) {
#line 200 "dppcon.f"
	if (upper) {

/*           Multiply by inv(U**T). */

#line 204 "dppcon.f"
	    dlatps_("Upper", "Transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scalel, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 206 "dppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 210 "dppcon.f"
	    dlatps_("Upper", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scaleu, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 212 "dppcon.f"
	} else {

/*           Multiply by inv(L). */

#line 216 "dppcon.f"
	    dlatps_("Lower", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scalel, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 218 "dppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**T). */

#line 222 "dppcon.f"
	    dlatps_("Lower", "Transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scaleu, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 224 "dppcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 228 "dppcon.f"
	scale = scalel * scaleu;
#line 229 "dppcon.f"
	if (scale != 1.) {
#line 230 "dppcon.f"
	    ix = idamax_(n, &work[1], &c__1);
#line 231 "dppcon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 231 "dppcon.f"
		goto L20;
#line 231 "dppcon.f"
	    }
#line 233 "dppcon.f"
	    drscl_(n, &scale, &work[1], &c__1);
#line 234 "dppcon.f"
	}
#line 235 "dppcon.f"
	goto L10;
#line 236 "dppcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 240 "dppcon.f"
    if (ainvnm != 0.) {
#line 240 "dppcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 240 "dppcon.f"
    }

#line 243 "dppcon.f"
L20:
#line 244 "dppcon.f"
    return 0;

/*     End of DPPCON */

} /* dppcon_ */

