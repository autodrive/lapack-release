#line 1 "dtrcon.f"
/* dtrcon.f -- translated by f2c (version 20100827).
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

#line 1 "dtrcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DTRCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRCON estimates the reciprocal of the condition number of a */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The triangular matrix A.  If UPLO = 'U', the leading N-by-N */
/* >          upper triangular part of the array A contains the upper */
/* >          triangular matrix, and the strictly lower triangular part of */
/* >          A is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/* >          triangular part of the array A contains the lower triangular */
/* >          matrix, and the strictly upper triangular part of A is not */
/* >          referenced.  If DIAG = 'U', the diagonal elements of A are */
/* >          also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
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
/* Subroutine */ int dtrcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *rcond, doublereal *work, 
	integer *iwork, integer *info, ftnlen norm_len, ftnlen uplo_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
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
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int dlatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 186 "dtrcon.f"
    /* Parameter adjustments */
#line 186 "dtrcon.f"
    a_dim1 = *lda;
#line 186 "dtrcon.f"
    a_offset = 1 + a_dim1;
#line 186 "dtrcon.f"
    a -= a_offset;
#line 186 "dtrcon.f"
    --work;
#line 186 "dtrcon.f"
    --iwork;
#line 186 "dtrcon.f"

#line 186 "dtrcon.f"
    /* Function Body */
#line 186 "dtrcon.f"
    *info = 0;
#line 187 "dtrcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 188 "dtrcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 189 "dtrcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 191 "dtrcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 192 "dtrcon.f"
	*info = -1;
#line 193 "dtrcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 194 "dtrcon.f"
	*info = -2;
#line 195 "dtrcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "dtrcon.f"
	*info = -3;
#line 197 "dtrcon.f"
    } else if (*n < 0) {
#line 198 "dtrcon.f"
	*info = -4;
#line 199 "dtrcon.f"
    } else if (*lda < max(1,*n)) {
#line 200 "dtrcon.f"
	*info = -6;
#line 201 "dtrcon.f"
    }
#line 202 "dtrcon.f"
    if (*info != 0) {
#line 203 "dtrcon.f"
	i__1 = -(*info);
#line 203 "dtrcon.f"
	xerbla_("DTRCON", &i__1, (ftnlen)6);
#line 204 "dtrcon.f"
	return 0;
#line 205 "dtrcon.f"
    }

/*     Quick return if possible */

#line 209 "dtrcon.f"
    if (*n == 0) {
#line 210 "dtrcon.f"
	*rcond = 1.;
#line 211 "dtrcon.f"
	return 0;
#line 212 "dtrcon.f"
    }

#line 214 "dtrcon.f"
    *rcond = 0.;
#line 215 "dtrcon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 219 "dtrcon.f"
    anorm = dlantr_(norm, uplo, diag, n, n, &a[a_offset], lda, &work[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 223 "dtrcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 227 "dtrcon.f"
	ainvnm = 0.;
#line 228 "dtrcon.f"
	*(unsigned char *)normin = 'N';
#line 229 "dtrcon.f"
	if (onenrm) {
#line 230 "dtrcon.f"
	    kase1 = 1;
#line 231 "dtrcon.f"
	} else {
#line 232 "dtrcon.f"
	    kase1 = 2;
#line 233 "dtrcon.f"
	}
#line 234 "dtrcon.f"
	kase = 0;
#line 235 "dtrcon.f"
L10:
#line 236 "dtrcon.f"
	dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 237 "dtrcon.f"
	if (kase != 0) {
#line 238 "dtrcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 242 "dtrcon.f"
		dlatrs_(uplo, "No transpose", diag, normin, n, &a[a_offset], 
			lda, &work[1], &scale, &work[(*n << 1) + 1], info, (
			ftnlen)1, (ftnlen)12, (ftnlen)1, (ftnlen)1);
#line 244 "dtrcon.f"
	    } else {

/*              Multiply by inv(A**T). */

#line 248 "dtrcon.f"
		dlatrs_(uplo, "Transpose", diag, normin, n, &a[a_offset], lda,
			 &work[1], &scale, &work[(*n << 1) + 1], info, (
			ftnlen)1, (ftnlen)9, (ftnlen)1, (ftnlen)1);
#line 250 "dtrcon.f"
	    }
#line 251 "dtrcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 255 "dtrcon.f"
	    if (scale != 1.) {
#line 256 "dtrcon.f"
		ix = idamax_(n, &work[1], &c__1);
#line 257 "dtrcon.f"
		xnorm = (d__1 = work[ix], abs(d__1));
#line 258 "dtrcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 258 "dtrcon.f"
		    goto L20;
#line 258 "dtrcon.f"
		}
#line 260 "dtrcon.f"
		drscl_(n, &scale, &work[1], &c__1);
#line 261 "dtrcon.f"
	    }
#line 262 "dtrcon.f"
	    goto L10;
#line 263 "dtrcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 267 "dtrcon.f"
	if (ainvnm != 0.) {
#line 267 "dtrcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 267 "dtrcon.f"
	}
#line 269 "dtrcon.f"
    }

#line 271 "dtrcon.f"
L20:
#line 272 "dtrcon.f"
    return 0;

/*     End of DTRCON */

} /* dtrcon_ */

