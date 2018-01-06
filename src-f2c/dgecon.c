#line 1 "dgecon.f"
/* dgecon.f -- translated by f2c (version 20100827).
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

#line 1 "dgecon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGECON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGECON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgecon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgecon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgecon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
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
/* > DGECON estimates the reciprocal of the condition number of a general */
/* > real matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by DGETRF. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by DGETRF. */
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
/* >          ANORM is DOUBLE PRECISION */
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgecon_(char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal sl;
    static integer ix;
    static doublereal su;
    static integer kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dlacn2_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int dlatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 173 "dgecon.f"
    /* Parameter adjustments */
#line 173 "dgecon.f"
    a_dim1 = *lda;
#line 173 "dgecon.f"
    a_offset = 1 + a_dim1;
#line 173 "dgecon.f"
    a -= a_offset;
#line 173 "dgecon.f"
    --work;
#line 173 "dgecon.f"
    --iwork;
#line 173 "dgecon.f"

#line 173 "dgecon.f"
    /* Function Body */
#line 173 "dgecon.f"
    *info = 0;
#line 174 "dgecon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 175 "dgecon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 176 "dgecon.f"
	*info = -1;
#line 177 "dgecon.f"
    } else if (*n < 0) {
#line 178 "dgecon.f"
	*info = -2;
#line 179 "dgecon.f"
    } else if (*lda < max(1,*n)) {
#line 180 "dgecon.f"
	*info = -4;
#line 181 "dgecon.f"
    } else if (*anorm < 0.) {
#line 182 "dgecon.f"
	*info = -5;
#line 183 "dgecon.f"
    }
#line 184 "dgecon.f"
    if (*info != 0) {
#line 185 "dgecon.f"
	i__1 = -(*info);
#line 185 "dgecon.f"
	xerbla_("DGECON", &i__1, (ftnlen)6);
#line 186 "dgecon.f"
	return 0;
#line 187 "dgecon.f"
    }

/*     Quick return if possible */

#line 191 "dgecon.f"
    *rcond = 0.;
#line 192 "dgecon.f"
    if (*n == 0) {
#line 193 "dgecon.f"
	*rcond = 1.;
#line 194 "dgecon.f"
	return 0;
#line 195 "dgecon.f"
    } else if (*anorm == 0.) {
#line 196 "dgecon.f"
	return 0;
#line 197 "dgecon.f"
    }

#line 199 "dgecon.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 203 "dgecon.f"
    ainvnm = 0.;
#line 204 "dgecon.f"
    *(unsigned char *)normin = 'N';
#line 205 "dgecon.f"
    if (onenrm) {
#line 206 "dgecon.f"
	kase1 = 1;
#line 207 "dgecon.f"
    } else {
#line 208 "dgecon.f"
	kase1 = 2;
#line 209 "dgecon.f"
    }
#line 210 "dgecon.f"
    kase = 0;
#line 211 "dgecon.f"
L10:
#line 212 "dgecon.f"
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 213 "dgecon.f"
    if (kase != 0) {
#line 214 "dgecon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 218 "dgecon.f"
	    dlatrs_("Lower", "No transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen)5,
		     (ftnlen)12, (ftnlen)4, (ftnlen)1);

/*           Multiply by inv(U). */

#line 223 "dgecon.f"
	    dlatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &su, &work[*n * 3 + 1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 225 "dgecon.f"
	} else {

/*           Multiply by inv(U**T). */

#line 229 "dgecon.f"
	    dlatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &su, &work[*n * 3 + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**T). */

#line 234 "dgecon.f"
	    dlatrs_("Lower", "Transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen)5,
		     (ftnlen)9, (ftnlen)4, (ftnlen)1);
#line 236 "dgecon.f"
	}

/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

#line 240 "dgecon.f"
	scale = sl * su;
#line 241 "dgecon.f"
	*(unsigned char *)normin = 'Y';
#line 242 "dgecon.f"
	if (scale != 1.) {
#line 243 "dgecon.f"
	    ix = idamax_(n, &work[1], &c__1);
#line 244 "dgecon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 244 "dgecon.f"
		goto L20;
#line 244 "dgecon.f"
	    }
#line 246 "dgecon.f"
	    drscl_(n, &scale, &work[1], &c__1);
#line 247 "dgecon.f"
	}
#line 248 "dgecon.f"
	goto L10;
#line 249 "dgecon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 253 "dgecon.f"
    if (ainvnm != 0.) {
#line 253 "dgecon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 253 "dgecon.f"
    }

#line 256 "dgecon.f"
L20:
#line 257 "dgecon.f"
    return 0;

/*     End of DGECON */

} /* dgecon_ */

