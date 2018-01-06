#line 1 "zppcon.f"
/* zppcon.f -- translated by f2c (version 20100827).
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

#line 1 "zppcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZPPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite packed matrix using */
/* > the Cholesky factorization A = U**H*U or A = L*L**H computed by */
/* > ZPPTRF. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H, packed columnwise in a linear */
/* >          array.  The j-th column of U or L is stored in the array AP */
/* >          as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          The 1-norm (or infinity-norm) of the Hermitian matrix A. */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zppcon_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1;
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
    extern /* Subroutine */ int zdrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static char normin[1];
    static doublereal smlnum;
    extern /* Subroutine */ int zlatps_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);


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

#line 174 "zppcon.f"
    /* Parameter adjustments */
#line 174 "zppcon.f"
    --rwork;
#line 174 "zppcon.f"
    --work;
#line 174 "zppcon.f"
    --ap;
#line 174 "zppcon.f"

#line 174 "zppcon.f"
    /* Function Body */
#line 174 "zppcon.f"
    *info = 0;
#line 175 "zppcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 176 "zppcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 177 "zppcon.f"
	*info = -1;
#line 178 "zppcon.f"
    } else if (*n < 0) {
#line 179 "zppcon.f"
	*info = -2;
#line 180 "zppcon.f"
    } else if (*anorm < 0.) {
#line 181 "zppcon.f"
	*info = -4;
#line 182 "zppcon.f"
    }
#line 183 "zppcon.f"
    if (*info != 0) {
#line 184 "zppcon.f"
	i__1 = -(*info);
#line 184 "zppcon.f"
	xerbla_("ZPPCON", &i__1, (ftnlen)6);
#line 185 "zppcon.f"
	return 0;
#line 186 "zppcon.f"
    }

/*     Quick return if possible */

#line 190 "zppcon.f"
    *rcond = 0.;
#line 191 "zppcon.f"
    if (*n == 0) {
#line 192 "zppcon.f"
	*rcond = 1.;
#line 193 "zppcon.f"
	return 0;
#line 194 "zppcon.f"
    } else if (*anorm == 0.) {
#line 195 "zppcon.f"
	return 0;
#line 196 "zppcon.f"
    }

#line 198 "zppcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 202 "zppcon.f"
    kase = 0;
#line 203 "zppcon.f"
    *(unsigned char *)normin = 'N';
#line 204 "zppcon.f"
L10:
#line 205 "zppcon.f"
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 206 "zppcon.f"
    if (kase != 0) {
#line 207 "zppcon.f"
	if (upper) {

/*           Multiply by inv(U**H). */

#line 211 "zppcon.f"
	    zlatps_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    ap[1], &work[1], &scalel, &rwork[1], info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 213 "zppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 217 "zppcon.f"
	    zlatps_("Upper", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scaleu, &rwork[1], info, (ftnlen)5, (ftnlen)12, 
		    (ftnlen)8, (ftnlen)1);
#line 219 "zppcon.f"
	} else {

/*           Multiply by inv(L). */

#line 223 "zppcon.f"
	    zlatps_("Lower", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scalel, &rwork[1], info, (ftnlen)5, (ftnlen)12, 
		    (ftnlen)8, (ftnlen)1);
#line 225 "zppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**H). */

#line 229 "zppcon.f"
	    zlatps_("Lower", "Conjugate transpose", "Non-unit", normin, n, &
		    ap[1], &work[1], &scaleu, &rwork[1], info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 231 "zppcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 235 "zppcon.f"
	scale = scalel * scaleu;
#line 236 "zppcon.f"
	if (scale != 1.) {
#line 237 "zppcon.f"
	    ix = izamax_(n, &work[1], &c__1);
#line 238 "zppcon.f"
	    i__1 = ix;
#line 238 "zppcon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 238 "zppcon.f"
		goto L20;
#line 238 "zppcon.f"
	    }
#line 240 "zppcon.f"
	    zdrscl_(n, &scale, &work[1], &c__1);
#line 241 "zppcon.f"
	}
#line 242 "zppcon.f"
	goto L10;
#line 243 "zppcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 247 "zppcon.f"
    if (ainvnm != 0.) {
#line 247 "zppcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 247 "zppcon.f"
    }

#line 250 "zppcon.f"
L20:
#line 251 "zppcon.f"
    return 0;

/*     End of ZPPCON */

} /* zppcon_ */

