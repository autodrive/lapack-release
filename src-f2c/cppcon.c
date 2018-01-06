#line 1 "cppcon.f"
/* cppcon.f -- translated by f2c (version 20100827).
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

#line 1 "cppcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CPPCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
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
/* > CPPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite packed matrix using */
/* > the Cholesky factorization A = U**H*U or A = L*L**H computed by */
/* > CPPTRF. */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H, packed columnwise in a linear */
/* >          array.  The j-th column of U or L is stored in the array AP */
/* >          as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cppcon_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1;
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), clatps_(
	    char *, char *, char *, char *, integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int csrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static char normin[1];
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 174 "cppcon.f"
    /* Parameter adjustments */
#line 174 "cppcon.f"
    --rwork;
#line 174 "cppcon.f"
    --work;
#line 174 "cppcon.f"
    --ap;
#line 174 "cppcon.f"

#line 174 "cppcon.f"
    /* Function Body */
#line 174 "cppcon.f"
    *info = 0;
#line 175 "cppcon.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 176 "cppcon.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 177 "cppcon.f"
	*info = -1;
#line 178 "cppcon.f"
    } else if (*n < 0) {
#line 179 "cppcon.f"
	*info = -2;
#line 180 "cppcon.f"
    } else if (*anorm < 0.) {
#line 181 "cppcon.f"
	*info = -4;
#line 182 "cppcon.f"
    }
#line 183 "cppcon.f"
    if (*info != 0) {
#line 184 "cppcon.f"
	i__1 = -(*info);
#line 184 "cppcon.f"
	xerbla_("CPPCON", &i__1, (ftnlen)6);
#line 185 "cppcon.f"
	return 0;
#line 186 "cppcon.f"
    }

/*     Quick return if possible */

#line 190 "cppcon.f"
    *rcond = 0.;
#line 191 "cppcon.f"
    if (*n == 0) {
#line 192 "cppcon.f"
	*rcond = 1.;
#line 193 "cppcon.f"
	return 0;
#line 194 "cppcon.f"
    } else if (*anorm == 0.) {
#line 195 "cppcon.f"
	return 0;
#line 196 "cppcon.f"
    }

#line 198 "cppcon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the 1-norm of the inverse. */

#line 202 "cppcon.f"
    kase = 0;
#line 203 "cppcon.f"
    *(unsigned char *)normin = 'N';
#line 204 "cppcon.f"
L10:
#line 205 "cppcon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 206 "cppcon.f"
    if (kase != 0) {
#line 207 "cppcon.f"
	if (upper) {

/*           Multiply by inv(U**H). */

#line 211 "cppcon.f"
	    clatps_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    ap[1], &work[1], &scalel, &rwork[1], info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 213 "cppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

#line 217 "cppcon.f"
	    clatps_("Upper", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scaleu, &rwork[1], info, (ftnlen)5, (ftnlen)12, 
		    (ftnlen)8, (ftnlen)1);
#line 219 "cppcon.f"
	} else {

/*           Multiply by inv(L). */

#line 223 "cppcon.f"
	    clatps_("Lower", "No transpose", "Non-unit", normin, n, &ap[1], &
		    work[1], &scalel, &rwork[1], info, (ftnlen)5, (ftnlen)12, 
		    (ftnlen)8, (ftnlen)1);
#line 225 "cppcon.f"
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L**H). */

#line 229 "cppcon.f"
	    clatps_("Lower", "Conjugate transpose", "Non-unit", normin, n, &
		    ap[1], &work[1], &scaleu, &rwork[1], info, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8, (ftnlen)1);
#line 231 "cppcon.f"
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

#line 235 "cppcon.f"
	scale = scalel * scaleu;
#line 236 "cppcon.f"
	if (scale != 1.) {
#line 237 "cppcon.f"
	    ix = icamax_(n, &work[1], &c__1);
#line 238 "cppcon.f"
	    i__1 = ix;
#line 238 "cppcon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 238 "cppcon.f"
		goto L20;
#line 238 "cppcon.f"
	    }
#line 240 "cppcon.f"
	    csrscl_(n, &scale, &work[1], &c__1);
#line 241 "cppcon.f"
	}
#line 242 "cppcon.f"
	goto L10;
#line 243 "cppcon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 247 "cppcon.f"
    if (ainvnm != 0.) {
#line 247 "cppcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 247 "cppcon.f"
    }

#line 250 "cppcon.f"
L20:
#line 251 "cppcon.f"
    return 0;

/*     End of CPPCON */

} /* cppcon_ */

