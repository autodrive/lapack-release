#line 1 "cla_syrcond_c.f"
/* cla_syrcond_c.f -- translated by f2c (version 20100827).
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

#line 1 "cla_syrcond_c.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLA_SYRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for symmetric i
ndefinite matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_SYRCOND_C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syr
cond_c.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syr
cond_c.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syr
cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_SYRCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C, */
/*                                    CAPPLY, INFO, WORK, RWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
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
/* >    CLA_SYRCOND_C Computes the infinity norm condition number of */
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
/* >     obtain the factor U or L as computed by CSYTRF. */
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
/* >     as determined by CSYTRF. */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
doublereal cla_syrcond_c__(char *uplo, integer *n, doublecomplex *a, integer *
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
    extern /* Subroutine */ int csytrs_(char *, integer *, integer *, 
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

#line 187 "cla_syrcond_c.f"
    /* Parameter adjustments */
#line 187 "cla_syrcond_c.f"
    a_dim1 = *lda;
#line 187 "cla_syrcond_c.f"
    a_offset = 1 + a_dim1;
#line 187 "cla_syrcond_c.f"
    a -= a_offset;
#line 187 "cla_syrcond_c.f"
    af_dim1 = *ldaf;
#line 187 "cla_syrcond_c.f"
    af_offset = 1 + af_dim1;
#line 187 "cla_syrcond_c.f"
    af -= af_offset;
#line 187 "cla_syrcond_c.f"
    --ipiv;
#line 187 "cla_syrcond_c.f"
    --c__;
#line 187 "cla_syrcond_c.f"
    --work;
#line 187 "cla_syrcond_c.f"
    --rwork;
#line 187 "cla_syrcond_c.f"

#line 187 "cla_syrcond_c.f"
    /* Function Body */
#line 187 "cla_syrcond_c.f"
    ret_val = 0.;

#line 189 "cla_syrcond_c.f"
    *info = 0;
#line 190 "cla_syrcond_c.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 191 "cla_syrcond_c.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 192 "cla_syrcond_c.f"
	*info = -1;
#line 193 "cla_syrcond_c.f"
    } else if (*n < 0) {
#line 194 "cla_syrcond_c.f"
	*info = -2;
#line 195 "cla_syrcond_c.f"
    } else if (*lda < max(1,*n)) {
#line 196 "cla_syrcond_c.f"
	*info = -4;
#line 197 "cla_syrcond_c.f"
    } else if (*ldaf < max(1,*n)) {
#line 198 "cla_syrcond_c.f"
	*info = -6;
#line 199 "cla_syrcond_c.f"
    }
#line 200 "cla_syrcond_c.f"
    if (*info != 0) {
#line 201 "cla_syrcond_c.f"
	i__1 = -(*info);
#line 201 "cla_syrcond_c.f"
	xerbla_("CLA_SYRCOND_C", &i__1, (ftnlen)13);
#line 202 "cla_syrcond_c.f"
	return ret_val;
#line 203 "cla_syrcond_c.f"
    }
#line 204 "cla_syrcond_c.f"
    up = FALSE_;
#line 205 "cla_syrcond_c.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 205 "cla_syrcond_c.f"
	up = TRUE_;
#line 205 "cla_syrcond_c.f"
    }

/*     Compute norm of op(A)*op2(C). */

#line 209 "cla_syrcond_c.f"
    anorm = 0.;
#line 210 "cla_syrcond_c.f"
    if (up) {
#line 211 "cla_syrcond_c.f"
	i__1 = *n;
#line 211 "cla_syrcond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "cla_syrcond_c.f"
	    tmp = 0.;
#line 213 "cla_syrcond_c.f"
	    if (*capply) {
#line 214 "cla_syrcond_c.f"
		i__2 = i__;
#line 214 "cla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 215 "cla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 215 "cla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 216 "cla_syrcond_c.f"
		}
#line 217 "cla_syrcond_c.f"
		i__2 = *n;
#line 217 "cla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 218 "cla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 218 "cla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 219 "cla_syrcond_c.f"
		}
#line 220 "cla_syrcond_c.f"
	    } else {
#line 221 "cla_syrcond_c.f"
		i__2 = i__;
#line 221 "cla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 222 "cla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 222 "cla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 223 "cla_syrcond_c.f"
		}
#line 224 "cla_syrcond_c.f"
		i__2 = *n;
#line 224 "cla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 225 "cla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 225 "cla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 226 "cla_syrcond_c.f"
		}
#line 227 "cla_syrcond_c.f"
	    }
#line 228 "cla_syrcond_c.f"
	    rwork[i__] = tmp;
#line 229 "cla_syrcond_c.f"
	    anorm = max(anorm,tmp);
#line 230 "cla_syrcond_c.f"
	}
#line 231 "cla_syrcond_c.f"
    } else {
#line 232 "cla_syrcond_c.f"
	i__1 = *n;
#line 232 "cla_syrcond_c.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "cla_syrcond_c.f"
	    tmp = 0.;
#line 234 "cla_syrcond_c.f"
	    if (*capply) {
#line 235 "cla_syrcond_c.f"
		i__2 = i__;
#line 235 "cla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 236 "cla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 236 "cla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) / c__[j];
#line 237 "cla_syrcond_c.f"
		}
#line 238 "cla_syrcond_c.f"
		i__2 = *n;
#line 238 "cla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 239 "cla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 239 "cla_syrcond_c.f"
		    tmp += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2))) / c__[j];
#line 240 "cla_syrcond_c.f"
		}
#line 241 "cla_syrcond_c.f"
	    } else {
#line 242 "cla_syrcond_c.f"
		i__2 = i__;
#line 242 "cla_syrcond_c.f"
		for (j = 1; j <= i__2; ++j) {
#line 243 "cla_syrcond_c.f"
		    i__3 = i__ + j * a_dim1;
#line 243 "cla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2));
#line 244 "cla_syrcond_c.f"
		}
#line 245 "cla_syrcond_c.f"
		i__2 = *n;
#line 245 "cla_syrcond_c.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 246 "cla_syrcond_c.f"
		    i__3 = j + i__ * a_dim1;
#line 246 "cla_syrcond_c.f"
		    tmp += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[
			    j + i__ * a_dim1]), abs(d__2));
#line 247 "cla_syrcond_c.f"
		}
#line 248 "cla_syrcond_c.f"
	    }
#line 249 "cla_syrcond_c.f"
	    rwork[i__] = tmp;
#line 250 "cla_syrcond_c.f"
	    anorm = max(anorm,tmp);
#line 251 "cla_syrcond_c.f"
	}
#line 252 "cla_syrcond_c.f"
    }

/*     Quick return if possible. */

#line 256 "cla_syrcond_c.f"
    if (*n == 0) {
#line 257 "cla_syrcond_c.f"
	ret_val = 1.;
#line 258 "cla_syrcond_c.f"
	return ret_val;
#line 259 "cla_syrcond_c.f"
    } else if (anorm == 0.) {
#line 260 "cla_syrcond_c.f"
	return ret_val;
#line 261 "cla_syrcond_c.f"
    }

/*     Estimate the norm of inv(op(A)). */

#line 265 "cla_syrcond_c.f"
    ainvnm = 0.;

#line 267 "cla_syrcond_c.f"
    kase = 0;
#line 268 "cla_syrcond_c.f"
L10:
#line 269 "cla_syrcond_c.f"
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
#line 270 "cla_syrcond_c.f"
    if (kase != 0) {
#line 271 "cla_syrcond_c.f"
	if (kase == 2) {

/*           Multiply by R. */

#line 275 "cla_syrcond_c.f"
	    i__1 = *n;
#line 275 "cla_syrcond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 276 "cla_syrcond_c.f"
		i__2 = i__;
#line 276 "cla_syrcond_c.f"
		i__3 = i__;
#line 276 "cla_syrcond_c.f"
		i__4 = i__;
#line 276 "cla_syrcond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 276 "cla_syrcond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 277 "cla_syrcond_c.f"
	    }

#line 279 "cla_syrcond_c.f"
	    if (up) {
#line 280 "cla_syrcond_c.f"
		csytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 282 "cla_syrcond_c.f"
	    } else {
#line 283 "cla_syrcond_c.f"
		csytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 285 "cla_syrcond_c.f"
	    }

/*           Multiply by inv(C). */

#line 289 "cla_syrcond_c.f"
	    if (*capply) {
#line 290 "cla_syrcond_c.f"
		i__1 = *n;
#line 290 "cla_syrcond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "cla_syrcond_c.f"
		    i__2 = i__;
#line 291 "cla_syrcond_c.f"
		    i__3 = i__;
#line 291 "cla_syrcond_c.f"
		    i__4 = i__;
#line 291 "cla_syrcond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 291 "cla_syrcond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 292 "cla_syrcond_c.f"
		}
#line 293 "cla_syrcond_c.f"
	    }
#line 294 "cla_syrcond_c.f"
	} else {

/*           Multiply by inv(C**T). */

#line 298 "cla_syrcond_c.f"
	    if (*capply) {
#line 299 "cla_syrcond_c.f"
		i__1 = *n;
#line 299 "cla_syrcond_c.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 300 "cla_syrcond_c.f"
		    i__2 = i__;
#line 300 "cla_syrcond_c.f"
		    i__3 = i__;
#line 300 "cla_syrcond_c.f"
		    i__4 = i__;
#line 300 "cla_syrcond_c.f"
		    z__1.r = c__[i__4] * work[i__3].r, z__1.i = c__[i__4] * 
			    work[i__3].i;
#line 300 "cla_syrcond_c.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 301 "cla_syrcond_c.f"
		}
#line 302 "cla_syrcond_c.f"
	    }

#line 304 "cla_syrcond_c.f"
	    if (up) {
#line 305 "cla_syrcond_c.f"
		csytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 307 "cla_syrcond_c.f"
	    } else {
#line 308 "cla_syrcond_c.f"
		csytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 310 "cla_syrcond_c.f"
	    }

/*           Multiply by R. */

#line 314 "cla_syrcond_c.f"
	    i__1 = *n;
#line 314 "cla_syrcond_c.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 315 "cla_syrcond_c.f"
		i__2 = i__;
#line 315 "cla_syrcond_c.f"
		i__3 = i__;
#line 315 "cla_syrcond_c.f"
		i__4 = i__;
#line 315 "cla_syrcond_c.f"
		z__1.r = rwork[i__4] * work[i__3].r, z__1.i = rwork[i__4] * 
			work[i__3].i;
#line 315 "cla_syrcond_c.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 316 "cla_syrcond_c.f"
	    }
#line 317 "cla_syrcond_c.f"
	}
#line 318 "cla_syrcond_c.f"
	goto L10;
#line 319 "cla_syrcond_c.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 323 "cla_syrcond_c.f"
    if (ainvnm != 0.) {
#line 323 "cla_syrcond_c.f"
	ret_val = 1. / ainvnm;
#line 323 "cla_syrcond_c.f"
    }

#line 326 "cla_syrcond_c.f"
    return ret_val;

} /* cla_syrcond_c__ */

