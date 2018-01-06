#line 1 "sppcon.f"
/* sppcon.f -- translated by f2c (version 20100827).
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

#line 1 "sppcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SPPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sppcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sppcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sppcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite packed matrix using */
/* > the Cholesky factorization A = U**T*U or A = L*L**T computed by */
/* > SPPTRF. */
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
/* >          AP is REAL array, dimension (N*(N+1)/2) */
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
/* >          ANORM is REAL */
/* >          The 1-norm (or infinity-norm) of the symmetric matrix A. */
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

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sppcon_(char *uplo, integer *n, doublereal *ap, 
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
    extern /* Subroutine */ int srscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static doublereal scalel;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal scaleu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal ainvnm;
    static char normin[1];
    extern /* Subroutine */ int slatps_(char *, char *, char *, char *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 167 "sppcon.f"
    /* Parameter adjustments */
#line 167 "sppcon.f"
    --iwork;
#line 167 "sppcon.f"
    --work;
#line 167 "sppcon.f"
    --ap;
#line 167 "sppcon.f"

#line 167 "sppcon.f"
    /* Function Body */
#line 167 "sppcon.f"
    *info = 0;
#line 168 "sppcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 169 "sppcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 170 "sppcon.f"
	*info = -1;
#line 171 "sppcon.f"
    } else if (*n < 0) {
#line 172 "sppcon.f"
	*info = -2;
#line 173 "sppcon.f"
    } else if (*anorm < 0.) {
#line 174 "sppcon.f"
	*info = -4;
#line 175 "sppcon.f"
    }
#line 176 "sppcon.f"
    if (*info != 0) {
#line 177 "sppcon.f"
	i__1 = -(*info);
#line 177 "sppcon.f"
	xerbla_("SPPCON", &i__1, (ftnlen)6);
#line 178 "sppcon.f"
	return 0;
#line 179 "sppcon.f"
    }

/*     Quick return if possible */

#line 183 "sppcon.f"
    *rcond = 0.;
#line 184 "sppcon.f"
    if (*n == 0) {
#line 185 "sppcon.f"
	*rcond = 1.;
#line 186 "sppcon.f"
	return 0;
#line 187 "sppcon.f"
    } else if (*anorm == 0.) {
#line 188 "sppcon.f"
	return 0;
#line 189 "sppcon.f"
    }

#line 191 "sppcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 195 "sppcon.f"
    kase = 0;
#line 196 "sppcon.f"
    *(unsigned char *)normin = 'N';
#line 197 "sppcon.f"
L10:
#line 198 "sppcon.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 199 "sppcon.f"
    if (kase != 0) {
#line 200 "sppcon.f"
	if (upper) {

/*           Multiply by inv(U**T). */

#line 204 "sppcon.f"
	    slatps_("Upper", "Transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scalel, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 206 "sppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 210 "sppcon.f"
	    slatps_("Upper", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scaleu, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 212 "sppcon.f"
	} else {

/*           Multiply by inv(L). */

#line 216 "sppcon.f"
	    slatps_("Lower", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scalel, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 218 "sppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**T). */

#line 222 "sppcon.f"
	    slatps_("Lower", "Transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scaleu, &work[(*n << 1) + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 224 "sppcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 228 "sppcon.f"
	scale = scalel * scaleu;
#line 229 "sppcon.f"
	if (scale != 1.) {
#line 230 "sppcon.f"
	    ix = isamax_(n, &work[1], &c__1);
#line 231 "sppcon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 231 "sppcon.f"
		goto L20;
#line 231 "sppcon.f"
	    }
#line 233 "sppcon.f"
	    srscl_(n, &scale, &work[1], &c__1);
#line 234 "sppcon.f"
	}
#line 235 "sppcon.f"
	goto L10;
#line 236 "sppcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 240 "sppcon.f"
    if (ainvnm != 0.) {
#line 240 "sppcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 240 "sppcon.f"
    }

#line 243 "sppcon.f"
L20:
#line 244 "sppcon.f"
    return 0;

/*     End of SPPCON */

} /* sppcon_ */

