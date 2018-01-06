#line 1 "dtpcon.f"
/* dtpcon.f -- translated by f2c (version 20100827).
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

#line 1 "dtpcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DTPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   RCOND */
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
/* > DTPCON estimates the reciprocal of the condition number of a packed */
/* > triangular matrix A, in either the 1-norm or the infinity-norm. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* >          If DIAG = 'U', the diagonal elements of A are not referenced */
/* >          and are assumed to be 1. */
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
/* Subroutine */ int dtpcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *ap, doublereal *rcond, doublereal *work, integer *iwork, 
	integer *info, ftnlen norm_len, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1;
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal dlantp_(char *, char *, char *, integer *, doublereal *,
	     doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int dlatps_(char *, char *, char *, char *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 179 "dtpcon.f"
    /* Parameter adjustments */
#line 179 "dtpcon.f"
    --iwork;
#line 179 "dtpcon.f"
    --work;
#line 179 "dtpcon.f"
    --ap;
#line 179 "dtpcon.f"

#line 179 "dtpcon.f"
    /* Function Body */
#line 179 "dtpcon.f"
    *info = 0;
#line 180 "dtpcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 181 "dtpcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 182 "dtpcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 184 "dtpcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 185 "dtpcon.f"
	*info = -1;
#line 186 "dtpcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 187 "dtpcon.f"
	*info = -2;
#line 188 "dtpcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "dtpcon.f"
	*info = -3;
#line 190 "dtpcon.f"
    } else if (*n < 0) {
#line 191 "dtpcon.f"
	*info = -4;
#line 192 "dtpcon.f"
    }
#line 193 "dtpcon.f"
    if (*info != 0) {
#line 194 "dtpcon.f"
	i__1 = -(*info);
#line 194 "dtpcon.f"
	xerbla_("DTPCON", &i__1, (ftnlen)6);
#line 195 "dtpcon.f"
	return 0;
#line 196 "dtpcon.f"
    }

/*     Quick return if possible */

#line 200 "dtpcon.f"
    if (*n == 0) {
#line 201 "dtpcon.f"
	*rcond = 1.;
#line 202 "dtpcon.f"
	return 0;
#line 203 "dtpcon.f"
    }

#line 205 "dtpcon.f"
    *rcond = 0.;
#line 206 "dtpcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 210 "dtpcon.f"
    anorm = dlantp_(norm, uplo, diag, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)
	    1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 214 "dtpcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 218 "dtpcon.f"
	ainvnm = 0.;
#line 219 "dtpcon.f"
	*(unsigned char *)normin = 'N';
#line 220 "dtpcon.f"
	if (onenrm) {
#line 221 "dtpcon.f"
	    kase1 = 1;
#line 222 "dtpcon.f"
	} else {
#line 223 "dtpcon.f"
	    kase1 = 2;
#line 224 "dtpcon.f"
	}
#line 225 "dtpcon.f"
	kase = 0;
#line 226 "dtpcon.f"
L10:
#line 227 "dtpcon.f"
	dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 228 "dtpcon.f"
	if (kase != 0) {
#line 229 "dtpcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 233 "dtpcon.f"
		dlatps_(uplo, "No transpose", diag, normin, n, &ap[1], &work[
			1], &scale, &work[(*n << 1) + 1], info, (ftnlen)1, (
			ftnlen)12, (ftnlen)1, (ftnlen)1);
#line 235 "dtpcon.f"
	    } else {

/*              Multiply by inv(A**T). */

#line 239 "dtpcon.f"
		dlatps_(uplo, "Transpose", diag, normin, n, &ap[1], &work[1], 
			&scale, &work[(*n << 1) + 1], info, (ftnlen)1, (
			ftnlen)9, (ftnlen)1, (ftnlen)1);
#line 241 "dtpcon.f"
	    }
#line 242 "dtpcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 246 "dtpcon.f"
	    if (scale != 1.) {
#line 247 "dtpcon.f"
		ix = idamax_(n, &work[1], &c__1);
#line 248 "dtpcon.f"
		xnorm = (d__1 = work[ix], abs(d__1));
#line 249 "dtpcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 249 "dtpcon.f"
		    goto L20;
#line 249 "dtpcon.f"
		}
#line 251 "dtpcon.f"
		drscl_(n, &scale, &work[1], &c__1);
#line 252 "dtpcon.f"
	    }
#line 253 "dtpcon.f"
	    goto L10;
#line 254 "dtpcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 258 "dtpcon.f"
	if (ainvnm != 0.) {
#line 258 "dtpcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 258 "dtpcon.f"
	}
#line 260 "dtpcon.f"
    }

#line 262 "dtpcon.f"
L20:
#line 263 "dtpcon.f"
    return 0;

/*     End of DTPCON */

} /* dtpcon_ */

