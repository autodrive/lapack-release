#line 1 "cgecon.f"
/* cgecon.f -- translated by f2c (version 20100827).
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

#line 1 "cgecon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CGECON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGECON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgecon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgecon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgecon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
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
/* > CGECON estimates the reciprocal of the condition number of a general */
/* > complex matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by CGETRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as */
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
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by CGETRF. */
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
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
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
/* >          RWORK is REAL array, dimension (2*N) */
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

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgecon_(char *norm, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal sl;
    static integer ix;
    static doublereal su;
    static integer kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int clatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), csrscl_(integer *, doublereal *, doublecomplex *, 
	    integer *);
    static logical onenrm;
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

#line 180 "cgecon.f"
    /* Parameter adjustments */
#line 180 "cgecon.f"
    a_dim1 = *lda;
#line 180 "cgecon.f"
    a_offset = 1 + a_dim1;
#line 180 "cgecon.f"
    a -= a_offset;
#line 180 "cgecon.f"
    --work;
#line 180 "cgecon.f"
    --rwork;
#line 180 "cgecon.f"

#line 180 "cgecon.f"
    /* Function Body */
#line 180 "cgecon.f"
    *info = 0;
#line 181 "cgecon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 182 "cgecon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 183 "cgecon.f"
	*info = -1;
#line 184 "cgecon.f"
    } else if (*n < 0) {
#line 185 "cgecon.f"
	*info = -2;
#line 186 "cgecon.f"
    } else if (*lda < max(1,*n)) {
#line 187 "cgecon.f"
	*info = -4;
#line 188 "cgecon.f"
    } else if (*anorm < 0.) {
#line 189 "cgecon.f"
	*info = -5;
#line 190 "cgecon.f"
    }
#line 191 "cgecon.f"
    if (*info != 0) {
#line 192 "cgecon.f"
	i__1 = -(*info);
#line 192 "cgecon.f"
	xerbla_("CGECON", &i__1, (ftnlen)6);
#line 193 "cgecon.f"
	return 0;
#line 194 "cgecon.f"
    }

/*     Quick return if possible */

#line 198 "cgecon.f"
    *rcond = 0.;
#line 199 "cgecon.f"
    if (*n == 0) {
#line 200 "cgecon.f"
	*rcond = 1.;
#line 201 "cgecon.f"
	return 0;
#line 202 "cgecon.f"
    } else if (*anorm == 0.) {
#line 203 "cgecon.f"
	return 0;
#line 204 "cgecon.f"
    }

#line 206 "cgecon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 210 "cgecon.f"
    ainvnm = 0.;
#line 211 "cgecon.f"
    *(unsigned char *)normin = 'N';
#line 212 "cgecon.f"
    if (onenrm) {
#line 213 "cgecon.f"
	kase1 = 1;
#line 214 "cgecon.f"
    } else {
#line 215 "cgecon.f"
	kase1 = 2;
#line 216 "cgecon.f"
    }
#line 217 "cgecon.f"
    kase = 0;
#line 218 "cgecon.f"
L10:
#line 219 "cgecon.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 220 "cgecon.f"
    if (kase != 0) {
#line 221 "cgecon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 225 "cgecon.f"
	    clatrs_("Lower", "No transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &rwork[1], info, (ftnlen)5, (ftnlen)
		    12, (ftnlen)4, (ftnlen)1);

/*           Multiply by inv(U). */

#line 230 "cgecon.f"
	    clatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &su, &rwork[*n + 1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 232 "cgecon.f"
	} else {

/*           Multiply by inv(U**H). */

#line 236 "cgecon.f"
	    clatrs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &su, &rwork[*n + 1], info, (
		    ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**H). */

#line 242 "cgecon.f"
	    clatrs_("Lower", "Conjugate transpose", "Unit", normin, n, &a[
		    a_offset], lda, &work[1], &sl, &rwork[1], info, (ftnlen)5,
		     (ftnlen)19, (ftnlen)4, (ftnlen)1);
#line 244 "cgecon.f"
	}

/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

#line 248 "cgecon.f"
	scale = sl * su;
#line 249 "cgecon.f"
	*(unsigned char *)normin = 'Y';
#line 250 "cgecon.f"
	if (scale != 1.) {
#line 251 "cgecon.f"
	    ix = icamax_(n, &work[1], &c__1);
#line 252 "cgecon.f"
	    i__1 = ix;
#line 252 "cgecon.f"
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
#line 252 "cgecon.f"
		goto L20;
#line 252 "cgecon.f"
	    }
#line 254 "cgecon.f"
	    csrscl_(n, &scale, &work[1], &c__1);
#line 255 "cgecon.f"
	}
#line 256 "cgecon.f"
	goto L10;
#line 257 "cgecon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 261 "cgecon.f"
    if (ainvnm != 0.) {
#line 261 "cgecon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 261 "cgecon.f"
    }

#line 264 "cgecon.f"
L20:
#line 265 "cgecon.f"
    return 0;

/*     End of CGECON */

} /* cgecon_ */

