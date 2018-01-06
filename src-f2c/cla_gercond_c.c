#line 1 "cla_gercond_c.f"
/* cla_gercond_c.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gercond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_GERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general mat
rices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GERCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_GERCOND_C( TRANS, N, A, LDA, AF, LDAF, IPIV, C, */
/*                                    CAPPLY, INFO, WORK, RWORK ) */

/*       .. Scalar Aguments .. */
/*       CHARACTER          TRANS */
/*       LOGICAL            CAPPLY */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/*       REAL               C( * ), RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* >    CLA_GERCOND_C computes the infinity norm condition number of */
/* >    op(A) * inv(diag(C)) where C is a REAL vector. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >     Specifies the form of the system of equations: */
/* >       = 'N':  A * X = B     (No transpose) */
/* >       = 'T':  A**T * X = B  (Transpose) */
/* >       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >     The number of linear equations, i.e., the order of the */
/* >     matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >     The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by CGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >     The pivot indices from the factorization A = P*L*U */
/* >     as computed by CGETRF; row i of the matrix was interchanged */
/* >     with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >     The vector C in the formula op(A) * inv(diag(C)). */
/* > \endverbatim */
/* > */
/* > \param[in] CAPPLY */
/* > \verbatim */
/* >          CAPPLY is LOGICAL */
/* >     If .TRUE. then access the vector C in the formula above. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >       = 0:  Successful exit. */
/* >     i > 0:  The ith argument is invalid. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N). */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N). */
/* >     Workspace. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
doublereal cla_gercond_c__(char *trans, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *
	c__, logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen), cgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    static doublereal ainvnm;
    static logical notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Aguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 189 "cla_gercond_c.f"
    /* Parameter adjustments */
#line 189 "cla_gercond_c.f"
    a_dim1 = *lda;
#line 189 "cla_gercond_c.f"
    a_offset = 1 + a_dim1;
#line 189 "cla_gercond_c.f"
    a -= a_offset;
#line 189 "cla_gercond_c.f"
    af_dim1 = *ldaf;
#line 189 "cla_gercond_c.f"
    af_offset = 1 + af_dim1;
#line 189 "cla_gercond_c.f"
    af -= af_offset;
#line 189 "cla_gercond_c.f"
    --ipiv;
#line 189 "cla_gercond_c.f"
    --c__;
#line 189 "cla_gercond_c.f"
    --work;
#line 189 "cla_gercond_c.f"
    --rwork;
#line 189 "cla_gercond_c.f"

#line 189 "cla_gercond_c.f"
    /* Function Body */
#line 189 "cla_gercond_c.f"
    ret_val = 0.;

#line 191 "cla_gercond_c.f"
    *info = 0;
#line 192 "cla_gercond_c.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 193 "cla_gercond_c.f"
    if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 195 "cla_gercond_c.f"
	*info = -1;
#line 196 "cla_gercond_c.f"
    } else if (*n < 0) {
#line 197 "cla_gercond_c.f"
	*info = -2;
#line 198 "cla_gercond_c.f"
    } else if (*lda < max(1,*n)) {
#line 199 "cla_gercond_c.f"
	*info = -4;
#line 200 "cla_gercond_c.f"
    } else if (*ldaf < max(1,*n)) {
#line 201 "cla_gercond_c.f"
	*info = -6;
#line 202 "cla_gercond_c.f"
    }
#line 203 "cla_gercond_c.f"
    if (*info != 0) {
#line 204 "cla_gercond_c.f"
	i__1 = -(*info);
#line 204 "cla_gercond_c.f"
	xerbla_("CLA_GERCOND_C", &i__1, (ftnlen)13);
#line 205 "cla_gercond_c.f"
	return ret_val;
#line 206 "cla_gercond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 210 "cla_gercond_c.f"
    anorm = 0.;
#line 211 "cla_gercond_c.f"
    if (notrans) {
#line 212 "cla_gercond_c.f"
	i__1 = *n;
#line 212 "cla_gercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 213 "cla_gercond_c.f"
	    tmp = 0.;
#line 214 "cla_gercond_c.f"
	    if (*capply) {
#line 215 "cla_gercond_c.f"
		i__2 = *n;
#line 215 "cla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 216 "cla_gercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 216 "cla_gercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 217 "cla_gercond_c.f"
		}
#line 218 "cla_gercond_c.f"
	    } else {
#line 219 "cla_gercond_c.f"
		i__2 = *n;
#line 219 "cla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 220 "cla_gercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 220 "cla_gercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 221 "cla_gercond_c.f"
		}
#line 222 "cla_gercond_c.f"
	    }
#line 223 "cla_gercond_c.f"
	    rwork[i__] = tmp;
#line 224 "cla_gercond_c.f"
	    anorm = max(anorm,tmp);
#line 225 "cla_gercond_c.f"
	}
#line 226 "cla_gercond_c.f"
    } else {
#line 227 "cla_gercond_c.f"
	i__1 = *n;
#line 227 "cla_gercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 228 "cla_gercond_c.f"
	    tmp = 0.;
#line 229 "cla_gercond_c.f"
	    if (*capply) {
#line 230 "cla_gercond_c.f"
		i__2 = *n;
#line 230 "cla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 231 "cla_gercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 231 "cla_gercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 232 "cla_gercond_c.f"
		}
#line 233 "cla_gercond_c.f"
	    } else {
#line 234 "cla_gercond_c.f"
		i__2 = *n;
#line 234 "cla_gercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 235 "cla_gercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 235 "cla_gercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 236 "cla_gercond_c.f"
		}
#line 237 "cla_gercond_c.f"
	    }
#line 238 "cla_gercond_c.f"
	    rwork[i__] = tmp;
#line 239 "cla_gercond_c.f"
	    anorm = max(anorm,tmp);
#line 240 "cla_gercond_c.f"
	}
#line 241 "cla_gercond_c.f"
    }

/*     Quick return if possible. */

#line 245 "cla_gercond_c.f"
    if (*n == 0) {
#line 246 "cla_gercond_c.f"
	ret_val = 1.;
#line 247 "cla_gercond_c.f"
	return ret_val;
#line 248 "cla_gercond_c.f"
    } else if (anorm == 0.) {
#line 249 "cla_gercond_c.f"
	return ret_val;
#line 250 "cla_gercond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 254 "cla_gercond_c.f"
    ainvnm = 0.;

#line 256 "cla_gercond_c.f"
    kase = 0;
#line 257 "cla_gercond_c.f"
L10:
#line 258 "cla_gercond_c.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 259 "cla_gercond_c.f"
    if (kase != 0) {
#line 260 "cla_gercond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 264 "cla_gercond_c.f"
	    i__1 = *n;
#line 264 "cla_gercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 265 "cla_gercond_c.f"
		i__2 = i__;
#line 265 "cla_gercond_c.f"
		i__3 = i__;
#line 265 "cla_gercond_c.f"
		i__4 = i__;
#line 265 "cla_gercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 265 "cla_gercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 266 "cla_gercond_c.f"
	    }

#line 268 "cla_gercond_c.f"
	    if (notrans) {
#line 269 "cla_gercond_c.f"
		cgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 271 "cla_gercond_c.f"
	    } else {
#line 272 "cla_gercond_c.f"
		cgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 274 "cla_gercond_c.f"
	    }

/*           Multiply by inv(C). */

#line 278 "cla_gercond_c.f"
	    if (*capply) {
#line 279 "cla_gercond_c.f"
		i__1 = *n;
#line 279 "cla_gercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "cla_gercond_c.f"
		    i__2 = i__;
#line 280 "cla_gercond_c.f"
		    i__3 = i__;
#line 280 "cla_gercond_c.f"
		    i__4 = i__;
#line 280 "cla_gercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 280 "cla_gercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 281 "cla_gercond_c.f"
		}
#line 282 "cla_gercond_c.f"
	    }
#line 283 "cla_gercond_c.f"
	} else {

/*           Multiply by inv(C**H). */

#line 287 "cla_gercond_c.f"
	    if (*capply) {
#line 288 "cla_gercond_c.f"
		i__1 = *n;
#line 288 "cla_gercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "cla_gercond_c.f"
		    i__2 = i__;
#line 289 "cla_gercond_c.f"
		    i__3 = i__;
#line 289 "cla_gercond_c.f"
		    i__4 = i__;
#line 289 "cla_gercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 289 "cla_gercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 290 "cla_gercond_c.f"
		}
#line 291 "cla_gercond_c.f"
	    }

#line 293 "cla_gercond_c.f"
	    if (notrans) {
#line 294 "cla_gercond_c.f"
		cgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf,
			 &ipiv[1], &work[1], n, info, (ftnlen)19);
#line 296 "cla_gercond_c.f"
	    } else {
#line 297 "cla_gercond_c.f"
		cgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[
			1], &work[1], n, info, (ftnlen)12);
#line 299 "cla_gercond_c.f"
	    }

/*           Multiply by R. */

#line 303 "cla_gercond_c.f"
	    i__1 = *n;
#line 303 "cla_gercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 304 "cla_gercond_c.f"
		i__2 = i__;
#line 304 "cla_gercond_c.f"
		i__3 = i__;
#line 304 "cla_gercond_c.f"
		i__4 = i__;
#line 304 "cla_gercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 304 "cla_gercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 305 "cla_gercond_c.f"
	    }
#line 306 "cla_gercond_c.f"
	}
#line 307 "cla_gercond_c.f"
	goto L10;
#line 308 "cla_gercond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 312 "cla_gercond_c.f"
    if (ainvnm != 0.) {
#line 312 "cla_gercond_c.f"
	ret_val = 1. / ainvnm;
#line 312 "cla_gercond_c.f"
    }

#line 315 "cla_gercond_c.f"
    return ret_val;

} /* cla_gercond_c__ */

