#line 1 "dsycon_rook.f"
/* dsycon_rook.f -- translated by f2c (version 20100827).
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

#line 1 "dsycon_rook.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DSYCON_ROOK </b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYCON_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsycon_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsycon_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsycon_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYCON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND, */
/*                               WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYCON_ROOK estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric matrix A using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by DSYTRF_ROOK. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
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
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSYTRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSYTRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          The 1-norm of the original matrix A. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \date April 2012 */

/* > \ingroup doubleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > \verbatim */
/* > */
/* >   April 2012, Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int dsycon_rook__(char *uplo, integer *n, doublereal *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublereal *work, integer *iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int dsytrs_rook__(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal ainvnm;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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

#line 190 "dsycon_rook.f"
    /* Parameter adjustments */
#line 190 "dsycon_rook.f"
    a_dim1 = *lda;
#line 190 "dsycon_rook.f"
    a_offset = 1 + a_dim1;
#line 190 "dsycon_rook.f"
    a -= a_offset;
#line 190 "dsycon_rook.f"
    --ipiv;
#line 190 "dsycon_rook.f"
    --work;
#line 190 "dsycon_rook.f"
    --iwork;
#line 190 "dsycon_rook.f"

#line 190 "dsycon_rook.f"
    /* Function Body */
#line 190 "dsycon_rook.f"
    *info = 0;
#line 191 "dsycon_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 192 "dsycon_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 193 "dsycon_rook.f"
	*info = -1;
#line 194 "dsycon_rook.f"
    } else if (*n < 0) {
#line 195 "dsycon_rook.f"
	*info = -2;
#line 196 "dsycon_rook.f"
    } else if (*lda < max(1,*n)) {
#line 197 "dsycon_rook.f"
	*info = -4;
#line 198 "dsycon_rook.f"
    } else if (*anorm < 0.) {
#line 199 "dsycon_rook.f"
	*info = -6;
#line 200 "dsycon_rook.f"
    }
#line 201 "dsycon_rook.f"
    if (*info != 0) {
#line 202 "dsycon_rook.f"
	i__1 = -(*info);
#line 202 "dsycon_rook.f"
	xerbla_("DSYCON_ROOK", &i__1, (ftnlen)11);
#line 203 "dsycon_rook.f"
	return 0;
#line 204 "dsycon_rook.f"
    }

/*     Quick return if possible */

#line 208 "dsycon_rook.f"
    *rcond = 0.;
#line 209 "dsycon_rook.f"
    if (*n == 0) {
#line 210 "dsycon_rook.f"
	*rcond = 1.;
#line 211 "dsycon_rook.f"
	return 0;
#line 212 "dsycon_rook.f"
    } else if (*anorm <= 0.) {
#line 213 "dsycon_rook.f"
	return 0;
#line 214 "dsycon_rook.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 218 "dsycon_rook.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 222 "dsycon_rook.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 223 "dsycon_rook.f"
	    if (ipiv[i__] > 0 && a[i__ + i__ * a_dim1] == 0.) {
#line 223 "dsycon_rook.f"
		return 0;
#line 223 "dsycon_rook.f"
	    }
#line 225 "dsycon_rook.f"
/* L10: */
#line 225 "dsycon_rook.f"
	}
#line 226 "dsycon_rook.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 230 "dsycon_rook.f"
	i__1 = *n;
#line 230 "dsycon_rook.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 231 "dsycon_rook.f"
	    if (ipiv[i__] > 0 && a[i__ + i__ * a_dim1] == 0.) {
#line 231 "dsycon_rook.f"
		return 0;
#line 231 "dsycon_rook.f"
	    }
#line 233 "dsycon_rook.f"
/* L20: */
#line 233 "dsycon_rook.f"
	}
#line 234 "dsycon_rook.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 238 "dsycon_rook.f"
    kase = 0;
#line 239 "dsycon_rook.f"
L30:
#line 240 "dsycon_rook.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 241 "dsycon_rook.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**T) or inv(U*D*U**T). */

#line 245 "dsycon_rook.f"
	dsytrs_rook__(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], 
		n, info, (ftnlen)1);
#line 246 "dsycon_rook.f"
	goto L30;
#line 247 "dsycon_rook.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 251 "dsycon_rook.f"
    if (ainvnm != 0.) {
#line 251 "dsycon_rook.f"
	*rcond = 1. / ainvnm / *anorm;
#line 251 "dsycon_rook.f"
    }

#line 254 "dsycon_rook.f"
    return 0;

/*     End of DSYCON_ROOK */

} /* dsycon_rook__ */

