#line 1 "spocon.f"
/* spocon.f -- translated by f2c (version 20100827).
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

#line 1 "spocon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SPOCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPOCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spocon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spocon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spocon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               ANORM, RCOND */
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
/* > SPOCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite matrix using the */
/* > Cholesky factorization A = U**T*U or A = L*L**T computed by SPOTRF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, as computed by SPOTRF. */
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

/* > \date November 2011 */

/* > \ingroup realPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int spocon_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
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
    extern /* Subroutine */ int slatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 170 "spocon.f"
    /* Parameter adjustments */
#line 170 "spocon.f"
    a_dim1 = *lda;
#line 170 "spocon.f"
    a_offset = 1 + a_dim1;
#line 170 "spocon.f"
    a -= a_offset;
#line 170 "spocon.f"
    --work;
#line 170 "spocon.f"
    --iwork;
#line 170 "spocon.f"

#line 170 "spocon.f"
    /* Function Body */
#line 170 "spocon.f"
    *info = 0;
#line 171 "spocon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 172 "spocon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 173 "spocon.f"
	*info = -1;
#line 174 "spocon.f"
    } else if (*n < 0) {
#line 175 "spocon.f"
	*info = -2;
#line 176 "spocon.f"
    } else if (*lda < max(1,*n)) {
#line 177 "spocon.f"
	*info = -4;
#line 178 "spocon.f"
    } else if (*anorm < 0.) {
#line 179 "spocon.f"
	*info = -5;
#line 180 "spocon.f"
    }
#line 181 "spocon.f"
    if (*info != 0) {
#line 182 "spocon.f"
	i__1 = -(*info);
#line 182 "spocon.f"
	xerbla_("SPOCON", &i__1, (ftnlen)6);
#line 183 "spocon.f"
	return 0;
#line 184 "spocon.f"
    }

/*     Quick return if possible */

#line 188 "spocon.f"
    *rcond = 0.;
#line 189 "spocon.f"
    if (*n == 0) {
#line 190 "spocon.f"
	*rcond = 1.;
#line 191 "spocon.f"
	return 0;
#line 192 "spocon.f"
    } else if (*anorm == 0.) {
#line 193 "spocon.f"
	return 0;
#line 194 "spocon.f"
    }

#line 196 "spocon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of inv(A). */

#line 200 "spocon.f"
    kase = 0;
#line 201 "spocon.f"
    *(unsigned char *)normin = 'N';
#line 202 "spocon.f"
L10:
#line 203 "spocon.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 204 "spocon.f"
    if (kase != 0) {
#line 205 "spocon.f"
	if (upper) {

/*           Multiply by inv(U**T). */

#line 209 "spocon.f"
	    slatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &scalel, &work[(*n << 1) + 1], info, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 211 "spocon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 215 "spocon.f"
	    slatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scaleu, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 217 "spocon.f"
	} else {

/*           Multiply by inv(L). */

#line 221 "spocon.f"
	    slatrs_("Lower", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scalel, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 223 "spocon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**T). */

#line 227 "spocon.f"
	    slatrs_("Lower", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &scaleu, &work[(*n << 1) + 1], info, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
#line 229 "spocon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 233 "spocon.f"
	scale = scalel * scaleu;
#line 234 "spocon.f"
	if (scale != 1.) {
#line 235 "spocon.f"
	    ix = isamax_(n, &work[1], &c__1);
#line 236 "spocon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 236 "spocon.f"
		goto L20;
#line 236 "spocon.f"
	    }
#line 238 "spocon.f"
	    srscl_(n, &scale, &work[1], &c__1);
#line 239 "spocon.f"
	}
#line 240 "spocon.f"
	goto L10;
#line 241 "spocon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 245 "spocon.f"
    if (ainvnm != 0.) {
#line 245 "spocon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 245 "spocon.f"
    }

#line 248 "spocon.f"
L20:
#line 249 "spocon.f"
    return 0;

/*     End of SPOCON */

} /* spocon_ */

