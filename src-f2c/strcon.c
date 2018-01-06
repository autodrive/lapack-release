#line 1 "strcon.f"
/* strcon.f -- translated by f2c (version 20100827).
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

#line 1 "strcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b STRCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRCON estimates the reciprocal of the condition number of a */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          RCOND is REAL */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(norm(A) * norm(inv(A))). */
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

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int strcon_(char *norm, char *uplo, char *diag, integer *n, 
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
    static doublereal anorm;
    extern /* Subroutine */ int srscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    static doublereal xnorm;
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal ainvnm;
    static logical onenrm;
    static char normin[1];
    extern doublereal slantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int slatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 186 "strcon.f"
    /* Parameter adjustments */
#line 186 "strcon.f"
    a_dim1 = *lda;
#line 186 "strcon.f"
    a_offset = 1 + a_dim1;
#line 186 "strcon.f"
    a -= a_offset;
#line 186 "strcon.f"
    --work;
#line 186 "strcon.f"
    --iwork;
#line 186 "strcon.f"

#line 186 "strcon.f"
    /* Function Body */
#line 186 "strcon.f"
    *info = 0;
#line 187 "strcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 188 "strcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 189 "strcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 191 "strcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 192 "strcon.f"
	*info = -1;
#line 193 "strcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 194 "strcon.f"
	*info = -2;
#line 195 "strcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "strcon.f"
	*info = -3;
#line 197 "strcon.f"
    } else if (*n < 0) {
#line 198 "strcon.f"
	*info = -4;
#line 199 "strcon.f"
    } else if (*lda < max(1,*n)) {
#line 200 "strcon.f"
	*info = -6;
#line 201 "strcon.f"
    }
#line 202 "strcon.f"
    if (*info != 0) {
#line 203 "strcon.f"
	i__1 = -(*info);
#line 203 "strcon.f"
	xerbla_("STRCON", &i__1, (ftnlen)6);
#line 204 "strcon.f"
	return 0;
#line 205 "strcon.f"
    }

/*     Quick return if possible */

#line 209 "strcon.f"
    if (*n == 0) {
#line 210 "strcon.f"
	*rcond = 1.;
#line 211 "strcon.f"
	return 0;
#line 212 "strcon.f"
    }

#line 214 "strcon.f"
    *rcond = 0.;
#line 215 "strcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 219 "strcon.f"
    anorm = slantr_(norm, uplo, diag, n, n, &a[a_offset], lda, &work[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 223 "strcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 227 "strcon.f"
	ainvnm = 0.;
#line 228 "strcon.f"
	*(unsigned char *)normin = 'N';
#line 229 "strcon.f"
	if (onenrm) {
#line 230 "strcon.f"
	    kase1 = 1;
#line 231 "strcon.f"
	} else {
#line 232 "strcon.f"
	    kase1 = 2;
#line 233 "strcon.f"
	}
#line 234 "strcon.f"
	kase = 0;
#line 235 "strcon.f"
L10:
#line 236 "strcon.f"
	slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 237 "strcon.f"
	if (kase != 0) {
#line 238 "strcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 242 "strcon.f"
		slatrs_(uplo, "No transpose", diag, normin, n, &a[a_offset], 
			lda, &work[1], &scale, &work[(*n << 1) + 1], info, (
			ftnlen)1, (ftnlen)12, (ftnlen)1, (ftnlen)1);
#line 244 "strcon.f"
	    } else {

/*              Multiply by inv(A**T). */

#line 248 "strcon.f"
		slatrs_(uplo, "Transpose", diag, normin, n, &a[a_offset], lda,
			 &work[1], &scale, &work[(*n << 1) + 1], info, (
			ftnlen)1, (ftnlen)9, (ftnlen)1, (ftnlen)1);
#line 250 "strcon.f"
	    }
#line 251 "strcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 255 "strcon.f"
	    if (scale != 1.) {
#line 256 "strcon.f"
		ix = isamax_(n, &work[1], &c__1);
#line 257 "strcon.f"
		xnorm = (d__1 = work[ix], abs(d__1));
#line 258 "strcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 258 "strcon.f"
		    goto L20;
#line 258 "strcon.f"
		}
#line 260 "strcon.f"
		srscl_(n, &scale, &work[1], &c__1);
#line 261 "strcon.f"
	    }
#line 262 "strcon.f"
	    goto L10;
#line 263 "strcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 267 "strcon.f"
	if (ainvnm != 0.) {
#line 267 "strcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 267 "strcon.f"
	}
#line 269 "strcon.f"
    }

#line 271 "strcon.f"
L20:
#line 272 "strcon.f"
    return 0;

/*     End of STRCON */

} /* strcon_ */

