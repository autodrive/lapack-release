#line 1 "cpocon.f"
/* cpocon.f -- translated by f2c (version 20100827).
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

#line 1 "cpocon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CPOCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPOCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpocon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpocon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpocon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               ANORM, RCOND */
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
/* > CPOCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite matrix using the */
/* > Cholesky factorization A = U**H*U or A = L*L**H computed by CPOTRF. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H, as computed by CPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          The 1-norm (or infinity-norm) of the Hermitian matrix A. */
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

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpocon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer ix, kase;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    static doublereal scalel;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal scaleu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int clatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), csrscl_(integer *, doublereal *, doublecomplex *, 
	    integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 177 "cpocon.f"
    /* Parameter adjustments */
#line 177 "cpocon.f"
    a_dim1 = *lda;
#line 177 "cpocon.f"
    a_offset = 1 + a_dim1;
#line 177 "cpocon.f"
    a -= a_offset;
#line 177 "cpocon.f"
    --work;
#line 177 "cpocon.f"
    --rwork;
#line 177 "cpocon.f"

#line 177 "cpocon.f"
    /* Function Body */
#line 177 "cpocon.f"
    *info = 0;
#line 178 "cpocon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 179 "cpocon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 180 "cpocon.f"
	*info = -1;
#line 181 "cpocon.f"
    } else if (*n < 0) {
#line 182 "cpocon.f"
	*info = -2;
#line 183 "cpocon.f"
    } else if (*lda < max(1,*n)) {
#line 184 "cpocon.f"
	*info = -4;
#line 185 "cpocon.f"
    } else if (*anorm < 0.) {
#line 186 "cpocon.f"
	*info = -5;
#line 187 "cpocon.f"
    }
#line 188 "cpocon.f"
    if (*info != 0) {
#line 189 "cpocon.f"
	i__1 = -(*info);
#line 189 "cpocon.f"
	xerbla_("CPOCON", &i__1, (ftnlen)6);
#line 190 "cpocon.f"
	return 0;
#line 191 "cpocon.f"
    }

/*     Quick return if possible */

#line 195 "cpocon.f"
    *rcond = 0.;
#line 196 "cpocon.f"
    if (*n == 0) {
#line 197 "cpocon.f"
	*rcond = 1.;
#line 198 "cpocon.f"
	return 0;
#line 199 "cpocon.f"
    } else if (*anorm == 0.) {
#line 200 "cpocon.f"
	return 0;
#line 201 "cpocon.f"
    }

#line 203 "cpocon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of inv(A). */

#line 207 "cpocon.f"
    kase = 0;
#line 208 "cpocon.f"
    *(unsigned char *)normin = 'N';
#line 209 "cpocon.f"
L10:
#line 210 "cpocon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 211 "cpocon.f"
    if (kase != 0) {
#line 212 "cpocon.f"
	if (upper) {

/*           Multiply by inv(U**H). */

#line 216 "cpocon.f"
	    clatrs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scalel, &rwork[1], info, (
		    ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 218 "cpocon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 222 "cpocon.f"
	    clatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scaleu, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 224 "cpocon.f"
	} else {

/*           Multiply by inv(L). */

#line 228 "cpocon.f"
	    clatrs_("Lower", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scalel, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 230 "cpocon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**H). */

#line 234 "cpocon.f"
	    clatrs_("Lower", "Conjugate transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scaleu, &rwork[1], info, (
		    ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 236 "cpocon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 240 "cpocon.f"
	scale = scalel * scaleu;
#line 241 "cpocon.f"
	if (scale != 1.) {
#line 242 "cpocon.f"
	    ix = icamax_(n, &work[1], &c__1);
#line 243 "cpocon.f"
	    i__1 = ix;
#line 243 "cpocon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 243 "cpocon.f"
		goto L20;
#line 243 "cpocon.f"
	    }
#line 245 "cpocon.f"
	    csrscl_(n, &scale, &work[1], &c__1);
#line 246 "cpocon.f"
	}
#line 247 "cpocon.f"
	goto L10;
#line 248 "cpocon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 252 "cpocon.f"
    if (ainvnm != 0.) {
#line 252 "cpocon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 252 "cpocon.f"
    }

#line 255 "cpocon.f"
L20:
#line 256 "cpocon.f"
    return 0;

/*     End of CPOCON */

} /* cpocon_ */

