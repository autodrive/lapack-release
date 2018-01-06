#line 1 "cla_hercond_c.f"
/* cla_hercond_c.f -- translated by f2c (version 20100827).
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

#line 1 "cla_hercond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_HERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian i
ndefinite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_HERCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_her
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_her
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_her
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_HERCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C, */
/*                                    CAPPLY, INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       LOGICAL            CAPPLY */
/*       INTEGER            N, LDA, LDAF, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/*       REAL               C ( * ), RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLA_HERCOND_C computes the infinity norm condition number of */
/* >    op(A) * inv(diag(C)) where C is a REAL vector. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >       = 'U':  Upper triangle of A is stored; */
/* >       = 'L':  Lower triangle of A is stored. */
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
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by CHETRF. */
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
/* >     Details of the interchanges and the block structure of D */
/* >     as determined by CHETRF. */
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

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
doublereal cla_hercond_c__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *c__,
	 logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static logical up;
    static doublereal tmp;
    static integer kase;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static doublereal anorm;
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern /* Subroutine */ int chetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
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

#line 186 "cla_hercond_c.f"
    /* Parameter adjustments */
#line 186 "cla_hercond_c.f"
    a_dim1 = *lda;
#line 186 "cla_hercond_c.f"
    a_offset = 1 + a_dim1;
#line 186 "cla_hercond_c.f"
    a -= a_offset;
#line 186 "cla_hercond_c.f"
    af_dim1 = *ldaf;
#line 186 "cla_hercond_c.f"
    af_offset = 1 + af_dim1;
#line 186 "cla_hercond_c.f"
    af -= af_offset;
#line 186 "cla_hercond_c.f"
    --ipiv;
#line 186 "cla_hercond_c.f"
    --c__;
#line 186 "cla_hercond_c.f"
    --work;
#line 186 "cla_hercond_c.f"
    --rwork;
#line 186 "cla_hercond_c.f"

#line 186 "cla_hercond_c.f"
    /* Function Body */
#line 186 "cla_hercond_c.f"
    ret_val = 0.;

#line 188 "cla_hercond_c.f"
    *info = 0;
#line 189 "cla_hercond_c.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 190 "cla_hercond_c.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 191 "cla_hercond_c.f"
	*info = -1;
#line 192 "cla_hercond_c.f"
    } else if (*n < 0) {
#line 193 "cla_hercond_c.f"
	*info = -2;
#line 194 "cla_hercond_c.f"
    } else if (*lda < max(1,*n)) {
#line 195 "cla_hercond_c.f"
	*info = -4;
#line 196 "cla_hercond_c.f"
    } else if (*ldaf < max(1,*n)) {
#line 197 "cla_hercond_c.f"
	*info = -6;
#line 198 "cla_hercond_c.f"
    }
#line 199 "cla_hercond_c.f"
    if (*info != 0) {
#line 200 "cla_hercond_c.f"
	i__1 = -(*info);
#line 200 "cla_hercond_c.f"
	xerbla_("CLA_HERCOND_C", &i__1, (ftnlen)13);
#line 201 "cla_hercond_c.f"
	return ret_val;
#line 202 "cla_hercond_c.f"
    }
#line 203 "cla_hercond_c.f"
    up = FALSE_;
#line 204 "cla_hercond_c.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 204 "cla_hercond_c.f"
	up = TRUE_;
#line 204 "cla_hercond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 208 "cla_hercond_c.f"
    anorm = 0.;
#line 209 "cla_hercond_c.f"
    if (up) {
#line 210 "cla_hercond_c.f"
	i__1 = *n;
#line 210 "cla_hercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "cla_hercond_c.f"
	    tmp = 0.;
#line 212 "cla_hercond_c.f"
	    if (*capply) {
#line 213 "cla_hercond_c.f"
		i__2 = i__;
#line 213 "cla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 214 "cla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 214 "cla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 215 "cla_hercond_c.f"
		}
#line 216 "cla_hercond_c.f"
		i__2 = *n;
#line 216 "cla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 217 "cla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 217 "cla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 218 "cla_hercond_c.f"
		}
#line 219 "cla_hercond_c.f"
	    } else {
#line 220 "cla_hercond_c.f"
		i__2 = i__;
#line 220 "cla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 221 "cla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 221 "cla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 222 "cla_hercond_c.f"
		}
#line 223 "cla_hercond_c.f"
		i__2 = *n;
#line 223 "cla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 224 "cla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 224 "cla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 225 "cla_hercond_c.f"
		}
#line 226 "cla_hercond_c.f"
	    }
#line 227 "cla_hercond_c.f"
	    rwork[i__] = tmp;
#line 228 "cla_hercond_c.f"
	    anorm = max(anorm,tmp);
#line 229 "cla_hercond_c.f"
	}
#line 230 "cla_hercond_c.f"
    } else {
#line 231 "cla_hercond_c.f"
	i__1 = *n;
#line 231 "cla_hercond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "cla_hercond_c.f"
	    tmp = 0.;
#line 233 "cla_hercond_c.f"
	    if (*capply) {
#line 234 "cla_hercond_c.f"
		i__2 = i__;
#line 234 "cla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 235 "cla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 235 "cla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 236 "cla_hercond_c.f"
		}
#line 237 "cla_hercond_c.f"
		i__2 = *n;
#line 237 "cla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 238 "cla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 238 "cla_hercond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 239 "cla_hercond_c.f"
		}
#line 240 "cla_hercond_c.f"
	    } else {
#line 241 "cla_hercond_c.f"
		i__2 = i__;
#line 241 "cla_hercond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 242 "cla_hercond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 242 "cla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 243 "cla_hercond_c.f"
		}
#line 244 "cla_hercond_c.f"
		i__2 = *n;
#line 244 "cla_hercond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 245 "cla_hercond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 245 "cla_hercond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 246 "cla_hercond_c.f"
		}
#line 247 "cla_hercond_c.f"
	    }
#line 248 "cla_hercond_c.f"
	    rwork[i__] = tmp;
#line 249 "cla_hercond_c.f"
	    anorm = max(anorm,tmp);
#line 250 "cla_hercond_c.f"
	}
#line 251 "cla_hercond_c.f"
    }

/*     Quick return if possible. */

#line 255 "cla_hercond_c.f"
    if (*n == 0) {
#line 256 "cla_hercond_c.f"
	ret_val = 1.;
#line 257 "cla_hercond_c.f"
	return ret_val;
#line 258 "cla_hercond_c.f"
    } else if (anorm == 0.) {
#line 259 "cla_hercond_c.f"
	return ret_val;
#line 260 "cla_hercond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 264 "cla_hercond_c.f"
    ainvnm = 0.;

#line 266 "cla_hercond_c.f"
    kase = 0;
#line 267 "cla_hercond_c.f"
L10:
#line 268 "cla_hercond_c.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 269 "cla_hercond_c.f"
    if (kase != 0) {
#line 270 "cla_hercond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 274 "cla_hercond_c.f"
	    i__1 = *n;
#line 274 "cla_hercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "cla_hercond_c.f"
		i__2 = i__;
#line 275 "cla_hercond_c.f"
		i__3 = i__;
#line 275 "cla_hercond_c.f"
		i__4 = i__;
#line 275 "cla_hercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 275 "cla_hercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 276 "cla_hercond_c.f"
	    }

#line 278 "cla_hercond_c.f"
	    if (up) {
#line 279 "cla_hercond_c.f"
		chetrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 281 "cla_hercond_c.f"
	    } else {
#line 282 "cla_hercond_c.f"
		chetrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 284 "cla_hercond_c.f"
	    }

/*           Multiply by inv(C). */

#line 288 "cla_hercond_c.f"
	    if (*capply) {
#line 289 "cla_hercond_c.f"
		i__1 = *n;
#line 289 "cla_hercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "cla_hercond_c.f"
		    i__2 = i__;
#line 290 "cla_hercond_c.f"
		    i__3 = i__;
#line 290 "cla_hercond_c.f"
		    i__4 = i__;
#line 290 "cla_hercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 290 "cla_hercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 291 "cla_hercond_c.f"
		}
#line 292 "cla_hercond_c.f"
	    }
#line 293 "cla_hercond_c.f"
	} else {

/*           Multiply by inv(C**H). */

#line 297 "cla_hercond_c.f"
	    if (*capply) {
#line 298 "cla_hercond_c.f"
		i__1 = *n;
#line 298 "cla_hercond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "cla_hercond_c.f"
		    i__2 = i__;
#line 299 "cla_hercond_c.f"
		    i__3 = i__;
#line 299 "cla_hercond_c.f"
		    i__4 = i__;
#line 299 "cla_hercond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 299 "cla_hercond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 300 "cla_hercond_c.f"
		}
#line 301 "cla_hercond_c.f"
	    }

#line 303 "cla_hercond_c.f"
	    if (up) {
#line 304 "cla_hercond_c.f"
		chetrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 306 "cla_hercond_c.f"
	    } else {
#line 307 "cla_hercond_c.f"
		chetrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 309 "cla_hercond_c.f"
	    }

/*           Multiply by R. */

#line 313 "cla_hercond_c.f"
	    i__1 = *n;
#line 313 "cla_hercond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 314 "cla_hercond_c.f"
		i__2 = i__;
#line 314 "cla_hercond_c.f"
		i__3 = i__;
#line 314 "cla_hercond_c.f"
		i__4 = i__;
#line 314 "cla_hercond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 314 "cla_hercond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 315 "cla_hercond_c.f"
	    }
#line 316 "cla_hercond_c.f"
	}
#line 317 "cla_hercond_c.f"
	goto L10;
#line 318 "cla_hercond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 322 "cla_hercond_c.f"
    if (ainvnm != 0.) {
#line 322 "cla_hercond_c.f"
	ret_val = 1. / ainvnm;
#line 322 "cla_hercond_c.f"
    }

#line 325 "cla_hercond_c.f"
    return ret_val;

} /* cla_hercond_c__ */

