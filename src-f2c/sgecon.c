#line 1 "sgecon.f"
/* sgecon.f -- translated by f2c (version 20100827).
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

#line 1 "sgecon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SGECON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGECON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgecon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgecon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgecon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
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
/* > SGECON estimates the reciprocal of the condition number of a general */
/* > real matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by SGETRF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by SGETRF. */
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
/* >          WORK is REAL array, dimension (4*N) */
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

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgecon_(char *norm, integer *n, doublereal *a, integer *
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
    extern /* Subroutine */ int srscl_(integer *, doublereal *, doublereal *, 
	    integer *), slacn2_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal ainvnm;
    static logical onenrm;
    static char normin[1];
    extern /* Subroutine */ int slatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 173 "sgecon.f"
    /* Parameter adjustments */
#line 173 "sgecon.f"
    a_dim1 = *lda;
#line 173 "sgecon.f"
    a_offset = 1 + a_dim1;
#line 173 "sgecon.f"
    a -= a_offset;
#line 173 "sgecon.f"
    --work;
#line 173 "sgecon.f"
    --iwork;
#line 173 "sgecon.f"

#line 173 "sgecon.f"
    /* Function Body */
#line 173 "sgecon.f"
    *info = 0;
#line 174 "sgecon.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 175 "sgecon.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 176 "sgecon.f"
	*info = -1;
#line 177 "sgecon.f"
    } else if (*n < 0) {
#line 178 "sgecon.f"
	*info = -2;
#line 179 "sgecon.f"
    } else if (*lda < max(1,*n)) {
#line 180 "sgecon.f"
	*info = -4;
#line 181 "sgecon.f"
    } else if (*anorm < 0.) {
#line 182 "sgecon.f"
	*info = -5;
#line 183 "sgecon.f"
    }
#line 184 "sgecon.f"
    if (*info != 0) {
#line 185 "sgecon.f"
	i__1 = -(*info);
#line 185 "sgecon.f"
	xerbla_("SGECON", &i__1, (ftnlen)6);
#line 186 "sgecon.f"
	return 0;
#line 187 "sgecon.f"
    }

/*     Quick return if possible */

#line 191 "sgecon.f"
    *rcond = 0.;
#line 192 "sgecon.f"
    if (*n == 0) {
#line 193 "sgecon.f"
	*rcond = 1.;
#line 194 "sgecon.f"
	return 0;
#line 195 "sgecon.f"
    } else if (*anorm == 0.) {
#line 196 "sgecon.f"
	return 0;
#line 197 "sgecon.f"
    }

#line 199 "sgecon.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

#line 203 "sgecon.f"
    ainvnm = 0.;
#line 204 "sgecon.f"
    *(unsigned char *)normin = 'N';
#line 205 "sgecon.f"
    if (onenrm) {
#line 206 "sgecon.f"
	kase1 = 1;
#line 207 "sgecon.f"
    } else {
#line 208 "sgecon.f"
	kase1 = 2;
#line 209 "sgecon.f"
    }
#line 210 "sgecon.f"
    kase = 0;
#line 211 "sgecon.f"
L10:
#line 212 "sgecon.f"
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
#line 213 "sgecon.f"
    if (kase != 0) {
#line 214 "sgecon.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 218 "sgecon.f"
	    slatrs_("Lower", "No transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen)5,
		     (ftnlen)12, (ftnlen)4, (ftnlen)1);

/*           Multiply by inv(U). */

#line 223 "sgecon.f"
	    slatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &su, &work[*n * 3 + 1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 225 "sgecon.f"
	} else {

/*           Multiply by inv(U**T). */

#line 229 "sgecon.f"
	    slatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &su, &work[*n * 3 + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**T). */

#line 234 "sgecon.f"
	    slatrs_("Lower", "Transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen)5,
		     (ftnlen)9, (ftnlen)4, (ftnlen)1);
#line 236 "sgecon.f"
	}

/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

#line 240 "sgecon.f"
	scale = sl * su;
#line 241 "sgecon.f"
	*(unsigned char *)normin = 'Y';
#line 242 "sgecon.f"
	if (scale != 1.) {
#line 243 "sgecon.f"
	    ix = isamax_(n, &work[1], &c__1);
#line 244 "sgecon.f"
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
#line 244 "sgecon.f"
		goto L20;
#line 244 "sgecon.f"
	    }
#line 246 "sgecon.f"
	    srscl_(n, &scale, &work[1], &c__1);
#line 247 "sgecon.f"
	}
#line 248 "sgecon.f"
	goto L10;
#line 249 "sgecon.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 253 "sgecon.f"
    if (ainvnm != 0.) {
#line 253 "sgecon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 253 "sgecon.f"
    }

#line 256 "sgecon.f"
L20:
#line 257 "sgecon.f"
    return 0;

/*     End of SGECON */

} /* sgecon_ */

