#line 1 "ctrcon.f"
/* ctrcon.f -- translated by f2c (version 20100827).
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

#line 1 "ctrcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTRCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTRCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRCON estimates the reciprocal of the condition number of a */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* Subroutine */ int ctrcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen norm_len, ftnlen 
	uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
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
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal clantr_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int clatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), csrscl_(integer *, doublereal *, doublecomplex *, 
	    integer *);
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

#line 193 "ctrcon.f"
    /* Parameter adjustments */
#line 193 "ctrcon.f"
    a_dim1 = *lda;
#line 193 "ctrcon.f"
    a_offset = 1 + a_dim1;
#line 193 "ctrcon.f"
    a -= a_offset;
#line 193 "ctrcon.f"
    --work;
#line 193 "ctrcon.f"
    --rwork;
#line 193 "ctrcon.f"

#line 193 "ctrcon.f"
    /* Function Body */
#line 193 "ctrcon.f"
    *info = 0;
#line 194 "ctrcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 195 "ctrcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 196 "ctrcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 198 "ctrcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 199 "ctrcon.f"
	*info = -1;
#line 200 "ctrcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 201 "ctrcon.f"
	*info = -2;
#line 202 "ctrcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 203 "ctrcon.f"
	*info = -3;
#line 204 "ctrcon.f"
    } else if (*n < 0) {
#line 205 "ctrcon.f"
	*info = -4;
#line 206 "ctrcon.f"
    } else if (*lda < max(1,*n)) {
#line 207 "ctrcon.f"
	*info = -6;
#line 208 "ctrcon.f"
    }
#line 209 "ctrcon.f"
    if (*info != 0) {
#line 210 "ctrcon.f"
	i__1 = -(*info);
#line 210 "ctrcon.f"
	xerbla_("CTRCON", &i__1, (ftnlen)6);
#line 211 "ctrcon.f"
	return 0;
#line 212 "ctrcon.f"
    }

/*     Quick return if possible */

#line 216 "ctrcon.f"
    if (*n == 0) {
#line 217 "ctrcon.f"
	*rcond = 1.;
#line 218 "ctrcon.f"
	return 0;
#line 219 "ctrcon.f"
    }

#line 221 "ctrcon.f"
    *rcond = 0.;
#line 222 "ctrcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 226 "ctrcon.f"
    anorm = clantr_(norm, uplo, diag, n, n, &a[a_offset], lda, &rwork[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 230 "ctrcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 234 "ctrcon.f"
	ainvnm = 0.;
#line 235 "ctrcon.f"
	*(unsigned char *)normin = 'N';
#line 236 "ctrcon.f"
	if (onenrm) {
#line 237 "ctrcon.f"
	    kase1 = 1;
#line 238 "ctrcon.f"
	} else {
#line 239 "ctrcon.f"
	    kase1 = 2;
#line 240 "ctrcon.f"
	}
#line 241 "ctrcon.f"
	kase = 0;
#line 242 "ctrcon.f"
L10:
#line 243 "ctrcon.f"
	clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 244 "ctrcon.f"
	if (kase != 0) {
#line 245 "ctrcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 249 "ctrcon.f"
		clatrs_(uplo, "No transpose", diag, normin, n, &a[a_offset], 
			lda, &work[1], &scale, &rwork[1], info, (ftnlen)1, (
			ftnlen)12, (ftnlen)1, (ftnlen)1);
#line 251 "ctrcon.f"
	    } else {

/*              Multiply by inv(A**H). */

#line 255 "ctrcon.f"
		clatrs_(uplo, "Conjugate transpose", diag, normin, n, &a[
			a_offset], lda, &work[1], &scale, &rwork[1], info, (
			ftnlen)1, (ftnlen)19, (ftnlen)1, (ftnlen)1);
#line 257 "ctrcon.f"
	    }
#line 258 "ctrcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 262 "ctrcon.f"
	    if (scale != 1.) {
#line 263 "ctrcon.f"
		ix = icamax_(n, &work[1], &c__1);
#line 264 "ctrcon.f"
		i__1 = ix;
#line 264 "ctrcon.f"
		xnorm = (d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			work[ix]), abs(d__2));
#line 265 "ctrcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 265 "ctrcon.f"
		    goto L20;
#line 265 "ctrcon.f"
		}
#line 267 "ctrcon.f"
		csrscl_(n, &scale, &work[1], &c__1);
#line 268 "ctrcon.f"
	    }
#line 269 "ctrcon.f"
	    goto L10;
#line 270 "ctrcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 274 "ctrcon.f"
	if (ainvnm != 0.) {
#line 274 "ctrcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 274 "ctrcon.f"
	}
#line 276 "ctrcon.f"
    }

#line 278 "ctrcon.f"
L20:
#line 279 "ctrcon.f"
    return 0;

/*     End of CTRCON */

} /* ctrcon_ */

