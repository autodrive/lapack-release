#line 1 "ctpcon.f"
/* ctpcon.f -- translated by f2c (version 20100827).
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

#line 1 "ctpcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            INFO, N */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ) */
/*       COMPLEX            AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPCON estimates the reciprocal of the condition number of a packed */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctpcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info, ftnlen norm_len, ftnlen uplo_len, ftnlen 
	diag_len)
{
    /* System generated locals */
    integer i__1;
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
    extern doublereal clantp_(char *, char *, char *, integer *, 
	    doublecomplex *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int clatps_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 186 "ctpcon.f"
    /* Parameter adjustments */
#line 186 "ctpcon.f"
    --rwork;
#line 186 "ctpcon.f"
    --work;
#line 186 "ctpcon.f"
    --ap;
#line 186 "ctpcon.f"

#line 186 "ctpcon.f"
    /* Function Body */
#line 186 "ctpcon.f"
    *info = 0;
#line 187 "ctpcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 188 "ctpcon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 189 "ctpcon.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 191 "ctpcon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 192 "ctpcon.f"
	*info = -1;
#line 193 "ctpcon.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 194 "ctpcon.f"
	*info = -2;
#line 195 "ctpcon.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "ctpcon.f"
	*info = -3;
#line 197 "ctpcon.f"
    } else if (*n < 0) {
#line 198 "ctpcon.f"
	*info = -4;
#line 199 "ctpcon.f"
    }
#line 200 "ctpcon.f"
    if (*info != 0) {
#line 201 "ctpcon.f"
	i__1 = -(*info);
#line 201 "ctpcon.f"
	xerbla_("CTPCON", &i__1, (ftnlen)6);
#line 202 "ctpcon.f"
	return 0;
#line 203 "ctpcon.f"
    }

/*     Quick return if possible */

#line 207 "ctpcon.f"
    if (*n == 0) {
#line 208 "ctpcon.f"
	*rcond = 1.;
#line 209 "ctpcon.f"
	return 0;
#line 210 "ctpcon.f"
    }

#line 212 "ctpcon.f"
    *rcond = 0.;
#line 213 "ctpcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) * (doublereal) max(1,*n);

/*     Compute the norm of the triangular matrix A. */

#line 217 "ctpcon.f"
    anorm = clantp_(norm, uplo, diag, n, &ap[1], &rwork[1], (ftnlen)1, (
	    ftnlen)1, (ftnlen)1);

/*     Continue only if ANORM > 0. */

#line 221 "ctpcon.f"
    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

#line 225 "ctpcon.f"
	ainvnm = 0.;
#line 226 "ctpcon.f"
	*(unsigned char *)normin = 'N';
#line 227 "ctpcon.f"
	if (onenrm) {
#line 228 "ctpcon.f"
	    kase1 = 1;
#line 229 "ctpcon.f"
	} else {
#line 230 "ctpcon.f"
	    kase1 = 2;
#line 231 "ctpcon.f"
	}
#line 232 "ctpcon.f"
	kase = 0;
#line 233 "ctpcon.f"
L10:
#line 234 "ctpcon.f"
	clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 235 "ctpcon.f"
	if (kase != 0) {
#line 236 "ctpcon.f"
	    if (kase == kase1) {

/*              Multiply by inv(A). */

#line 240 "ctpcon.f"
		clatps_(uplo, "No transpose", diag, normin, n, &ap[1], &work[
			1], &scale, &rwork[1], info, (ftnlen)1, (ftnlen)12, (
			ftnlen)1, (ftnlen)1);
#line 242 "ctpcon.f"
	    } else {

/*              Multiply by inv(A**H). */

#line 246 "ctpcon.f"
		clatps_(uplo, "Conjugate transpose", diag, normin, n, &ap[1], 
			&work[1], &scale, &rwork[1], info, (ftnlen)1, (ftnlen)
			19, (ftnlen)1, (ftnlen)1);
#line 248 "ctpcon.f"
	    }
#line 249 "ctpcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overflow. */

#line 253 "ctpcon.f"
	    if (scale != 1.) {
#line 254 "ctpcon.f"
		ix = icamax_(n, &work[1], &c__1);
#line 255 "ctpcon.f"
		i__1 = ix;
#line 255 "ctpcon.f"
		xnorm = (d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			work[ix]), abs(d__2));
#line 256 "ctpcon.f"
		if (scale < xnorm * smlnum || scale == 0.) {
#line 256 "ctpcon.f"
		    goto L20;
#line 256 "ctpcon.f"
		}
#line 258 "ctpcon.f"
		csrscl_(n, &scale, &work[1], &c__1);
#line 259 "ctpcon.f"
	    }
#line 260 "ctpcon.f"
	    goto L10;
#line 261 "ctpcon.f"
	}

/*        Compute the estimate of the reciprocal condition number. */

#line 265 "ctpcon.f"
	if (ainvnm != 0.) {
#line 265 "ctpcon.f"
	    *rcond = 1. / anorm / ainvnm;
#line 265 "ctpcon.f"
	}
#line 267 "ctpcon.f"
    }

#line 269 "ctpcon.f"
L20:
#line 270 "ctpcon.f"
    return 0;

/*     End of CTPCON */

} /* ctpcon_ */

