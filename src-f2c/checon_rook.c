#line 1 "checon_rook.f"
/* checon_rook.f -- translated by f2c (version 20100827).
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

#line 1 "checon_rook.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHECON_ROOK estimates the reciprocal of the condition number fort HE matrices using factorizati
on obtained with one of the bounded diagonal pivoting methods (max 2 interchanges) */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHECON_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/checon_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/checon_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/checon_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHECON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, */
/*                               INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHECON_ROOK estimates the reciprocal of the condition number of a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by CHETRF_ROOK. */
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
/* >          = 'U':  Upper triangular, form is A = U*D*U**H; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**H. */
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
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CHETRF_ROOK. */
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
/* >          as determined by CHETRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          The 1-norm of the original matrix A. */
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

/* > \date November 2013 */

/* > \ingroup complexHEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > \verbatim */
/* > */
/* >  November 2013,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int checon_rook__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int chetrs_rook__(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

#line 185 "checon_rook.f"
    /* Parameter adjustments */
#line 185 "checon_rook.f"
    a_dim1 = *lda;
#line 185 "checon_rook.f"
    a_offset = 1 + a_dim1;
#line 185 "checon_rook.f"
    a -= a_offset;
#line 185 "checon_rook.f"
    --ipiv;
#line 185 "checon_rook.f"
    --work;
#line 185 "checon_rook.f"

#line 185 "checon_rook.f"
    /* Function Body */
#line 185 "checon_rook.f"
    *info = 0;
#line 186 "checon_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 187 "checon_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 188 "checon_rook.f"
	*info = -1;
#line 189 "checon_rook.f"
    } else if (*n < 0) {
#line 190 "checon_rook.f"
	*info = -2;
#line 191 "checon_rook.f"
    } else if (*lda < max(1,*n)) {
#line 192 "checon_rook.f"
	*info = -4;
#line 193 "checon_rook.f"
    } else if (*anorm < 0.) {
#line 194 "checon_rook.f"
	*info = -6;
#line 195 "checon_rook.f"
    }
#line 196 "checon_rook.f"
    if (*info != 0) {
#line 197 "checon_rook.f"
	i__1 = -(*info);
#line 197 "checon_rook.f"
	xerbla_("CHECON_ROOK", &i__1, (ftnlen)11);
#line 198 "checon_rook.f"
	return 0;
#line 199 "checon_rook.f"
    }

/*     Quick return if possible */

#line 203 "checon_rook.f"
    *rcond = 0.;
#line 204 "checon_rook.f"
    if (*n == 0) {
#line 205 "checon_rook.f"
	*rcond = 1.;
#line 206 "checon_rook.f"
	return 0;
#line 207 "checon_rook.f"
    } else if (*anorm <= 0.) {
#line 208 "checon_rook.f"
	return 0;
#line 209 "checon_rook.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 213 "checon_rook.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 217 "checon_rook.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 218 "checon_rook.f"
	    i__1 = i__ + i__ * a_dim1;
#line 218 "checon_rook.f"
	    if (ipiv[i__] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 218 "checon_rook.f"
		return 0;
#line 218 "checon_rook.f"
	    }
#line 220 "checon_rook.f"
/* L10: */
#line 220 "checon_rook.f"
	}
#line 221 "checon_rook.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 225 "checon_rook.f"
	i__1 = *n;
#line 225 "checon_rook.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 226 "checon_rook.f"
	    i__2 = i__ + i__ * a_dim1;
#line 226 "checon_rook.f"
	    if (ipiv[i__] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 226 "checon_rook.f"
		return 0;
#line 226 "checon_rook.f"
	    }
#line 228 "checon_rook.f"
/* L20: */
#line 228 "checon_rook.f"
	}
#line 229 "checon_rook.f"
    }

/*     Estimate the 1-norm of the inverse. */

#line 233 "checon_rook.f"
    kase = 0;
#line 234 "checon_rook.f"
L30:
#line 235 "checon_rook.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 236 "checon_rook.f"
    if (kase != 0) {

/*        Multiply by inv(L*D*L**H) or inv(U*D*U**H). */

#line 240 "checon_rook.f"
	chetrs_rook__(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], 
		n, info, (ftnlen)1);
#line 241 "checon_rook.f"
	goto L30;
#line 242 "checon_rook.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 246 "checon_rook.f"
    if (ainvnm != 0.) {
#line 246 "checon_rook.f"
	*rcond = 1. / ainvnm / *anorm;
#line 246 "checon_rook.f"
    }

#line 249 "checon_rook.f"
    return 0;

/*     End of CHECON_ROOK */

} /* checon_rook__ */

