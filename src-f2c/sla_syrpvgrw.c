#line 1 "sla_syrpvgrw.f"
/* sla_syrpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "sla_syrpvgrw.f"
/* > \brief \b SLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefi
nite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_SYRPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syr
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syr
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syr
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, */
/*                                   WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            N, INFO, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > SLA_SYRPVGRW computes the reciprocal pivot growth factor */
/* > norm(A)/norm(U). The "max absolute element" norm is used. If this is */
/* > much less than 1, the stability of the LU factorization of the */
/* > (equilibrated) matrix A could be poor. This also means that the */
/* > solution X, estimated condition numbers, and error bounds could be */
/* > unreliable. */
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
/* > \param[in] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >     The value of INFO returned from SSYTRF, .i.e., the pivot in */
/* >     column INFO is exactly 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A. */
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
/* >          AF is REAL array, dimension (LDAF,N) */
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by SSYTRF. */
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
/* >     as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (2*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
doublereal sla_syrpvgrw__(char *uplo, integer *n, integer *info, doublereal *
	a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, 
	doublereal *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, kp;
    static doublereal tmp, amax, umax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ncols;
    static logical upper;
    static doublereal rpvgrw;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 155 "sla_syrpvgrw.f"
    /* Parameter adjustments */
#line 155 "sla_syrpvgrw.f"
    a_dim1 = *lda;
#line 155 "sla_syrpvgrw.f"
    a_offset = 1 + a_dim1;
#line 155 "sla_syrpvgrw.f"
    a -= a_offset;
#line 155 "sla_syrpvgrw.f"
    af_dim1 = *ldaf;
#line 155 "sla_syrpvgrw.f"
    af_offset = 1 + af_dim1;
#line 155 "sla_syrpvgrw.f"
    af -= af_offset;
#line 155 "sla_syrpvgrw.f"
    --ipiv;
#line 155 "sla_syrpvgrw.f"
    --work;
#line 155 "sla_syrpvgrw.f"

#line 155 "sla_syrpvgrw.f"
    /* Function Body */
#line 155 "sla_syrpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);
#line 156 "sla_syrpvgrw.f"
    if (*info == 0) {
#line 157 "sla_syrpvgrw.f"
	if (upper) {
#line 158 "sla_syrpvgrw.f"
	    ncols = 1;
#line 159 "sla_syrpvgrw.f"
	} else {
#line 160 "sla_syrpvgrw.f"
	    ncols = *n;
#line 161 "sla_syrpvgrw.f"
	}
#line 162 "sla_syrpvgrw.f"
    } else {
#line 163 "sla_syrpvgrw.f"
	ncols = *info;
#line 164 "sla_syrpvgrw.f"
    }
#line 166 "sla_syrpvgrw.f"
    rpvgrw = 1.;
#line 167 "sla_syrpvgrw.f"
    i__1 = *n << 1;
#line 167 "sla_syrpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 168 "sla_syrpvgrw.f"
	work[i__] = 0.;
#line 169 "sla_syrpvgrw.f"
    }

/*     Find the max magnitude entry of each column of A.  Compute the max */
/*     for all N columns so we can apply the pivot permutation while */
/*     looping below.  Assume a full factorization is the common case. */

#line 175 "sla_syrpvgrw.f"
    if (upper) {
#line 176 "sla_syrpvgrw.f"
	i__1 = *n;
#line 176 "sla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 177 "sla_syrpvgrw.f"
	    i__2 = j;
#line 177 "sla_syrpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 178 "sla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + i__];
#line 178 "sla_syrpvgrw.f"
		work[*n + i__] = max(d__2,d__3);
/* Computing MAX */
#line 179 "sla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + j];
#line 179 "sla_syrpvgrw.f"
		work[*n + j] = max(d__2,d__3);
#line 180 "sla_syrpvgrw.f"
	    }
#line 181 "sla_syrpvgrw.f"
	}
#line 182 "sla_syrpvgrw.f"
    } else {
#line 183 "sla_syrpvgrw.f"
	i__1 = *n;
#line 183 "sla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "sla_syrpvgrw.f"
	    i__2 = *n;
#line 184 "sla_syrpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 185 "sla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + i__];
#line 185 "sla_syrpvgrw.f"
		work[*n + i__] = max(d__2,d__3);
/* Computing MAX */
#line 186 "sla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + j];
#line 186 "sla_syrpvgrw.f"
		work[*n + j] = max(d__2,d__3);
#line 187 "sla_syrpvgrw.f"
	    }
#line 188 "sla_syrpvgrw.f"
	}
#line 189 "sla_syrpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of U or L.  Also */
/*     permute the magnitudes of A above so they're in the same order as */
/*     the factor. */

/*     The iteration orders and permutations were copied from ssytrs. */
/*     Calls to SSWAP would be severe overkill. */

#line 198 "sla_syrpvgrw.f"
    if (upper) {
#line 199 "sla_syrpvgrw.f"
	k = *n;
#line 200 "sla_syrpvgrw.f"
	while(k < ncols && k > 0) {
#line 201 "sla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 203 "sla_syrpvgrw.f"
		kp = ipiv[k];
#line 204 "sla_syrpvgrw.f"
		if (kp != k) {
#line 205 "sla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 206 "sla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 207 "sla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 208 "sla_syrpvgrw.f"
		}
#line 209 "sla_syrpvgrw.f"
		i__1 = k;
#line 209 "sla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 210 "sla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 210 "sla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
#line 211 "sla_syrpvgrw.f"
		}
#line 212 "sla_syrpvgrw.f"
		--k;
#line 213 "sla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 215 "sla_syrpvgrw.f"
		kp = -ipiv[k];
#line 216 "sla_syrpvgrw.f"
		tmp = work[*n + k - 1];
#line 217 "sla_syrpvgrw.f"
		work[*n + k - 1] = work[*n + kp];
#line 218 "sla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 219 "sla_syrpvgrw.f"
		i__1 = k - 1;
#line 219 "sla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 220 "sla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 220 "sla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
/* Computing MAX */
#line 221 "sla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + (k - 1) * af_dim1], abs(d__1)), 
			    d__3 = work[k - 1];
#line 221 "sla_syrpvgrw.f"
		    work[k - 1] = max(d__2,d__3);
#line 222 "sla_syrpvgrw.f"
		}
/* Computing MAX */
#line 223 "sla_syrpvgrw.f"
		d__2 = (d__1 = af[k + k * af_dim1], abs(d__1)), d__3 = work[k]
			;
#line 223 "sla_syrpvgrw.f"
		work[k] = max(d__2,d__3);
#line 224 "sla_syrpvgrw.f"
		k += -2;
#line 225 "sla_syrpvgrw.f"
	    }
#line 226 "sla_syrpvgrw.f"
	}
#line 227 "sla_syrpvgrw.f"
	k = ncols;
#line 228 "sla_syrpvgrw.f"
	while(k <= *n) {
#line 229 "sla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 230 "sla_syrpvgrw.f"
		kp = ipiv[k];
#line 231 "sla_syrpvgrw.f"
		if (kp != k) {
#line 232 "sla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 233 "sla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 234 "sla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 235 "sla_syrpvgrw.f"
		}
#line 236 "sla_syrpvgrw.f"
		++k;
#line 237 "sla_syrpvgrw.f"
	    } else {
#line 238 "sla_syrpvgrw.f"
		kp = -ipiv[k];
#line 239 "sla_syrpvgrw.f"
		tmp = work[*n + k];
#line 240 "sla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 241 "sla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 242 "sla_syrpvgrw.f"
		k += 2;
#line 243 "sla_syrpvgrw.f"
	    }
#line 244 "sla_syrpvgrw.f"
	}
#line 245 "sla_syrpvgrw.f"
    } else {
#line 246 "sla_syrpvgrw.f"
	k = 1;
#line 247 "sla_syrpvgrw.f"
	while(k <= ncols) {
#line 248 "sla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 250 "sla_syrpvgrw.f"
		kp = ipiv[k];
#line 251 "sla_syrpvgrw.f"
		if (kp != k) {
#line 252 "sla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 253 "sla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 254 "sla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 255 "sla_syrpvgrw.f"
		}
#line 256 "sla_syrpvgrw.f"
		i__1 = *n;
#line 256 "sla_syrpvgrw.f"
		for (i__ = k; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 257 "sla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 257 "sla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
#line 258 "sla_syrpvgrw.f"
		}
#line 259 "sla_syrpvgrw.f"
		++k;
#line 260 "sla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 262 "sla_syrpvgrw.f"
		kp = -ipiv[k];
#line 263 "sla_syrpvgrw.f"
		tmp = work[*n + k + 1];
#line 264 "sla_syrpvgrw.f"
		work[*n + k + 1] = work[*n + kp];
#line 265 "sla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 266 "sla_syrpvgrw.f"
		i__1 = *n;
#line 266 "sla_syrpvgrw.f"
		for (i__ = k + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 267 "sla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 267 "sla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
/* Computing MAX */
#line 268 "sla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + (k + 1) * af_dim1], abs(d__1)), 
			    d__3 = work[k + 1];
#line 268 "sla_syrpvgrw.f"
		    work[k + 1] = max(d__2,d__3);
#line 269 "sla_syrpvgrw.f"
		}
/* Computing MAX */
#line 270 "sla_syrpvgrw.f"
		d__2 = (d__1 = af[k + k * af_dim1], abs(d__1)), d__3 = work[k]
			;
#line 270 "sla_syrpvgrw.f"
		work[k] = max(d__2,d__3);
#line 271 "sla_syrpvgrw.f"
		k += 2;
#line 272 "sla_syrpvgrw.f"
	    }
#line 273 "sla_syrpvgrw.f"
	}
#line 274 "sla_syrpvgrw.f"
	k = ncols;
#line 275 "sla_syrpvgrw.f"
	while(k >= 1) {
#line 276 "sla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 277 "sla_syrpvgrw.f"
		kp = ipiv[k];
#line 278 "sla_syrpvgrw.f"
		if (kp != k) {
#line 279 "sla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 280 "sla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 281 "sla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 282 "sla_syrpvgrw.f"
		}
#line 283 "sla_syrpvgrw.f"
		--k;
#line 284 "sla_syrpvgrw.f"
	    } else {
#line 285 "sla_syrpvgrw.f"
		kp = -ipiv[k];
#line 286 "sla_syrpvgrw.f"
		tmp = work[*n + k];
#line 287 "sla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 288 "sla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 289 "sla_syrpvgrw.f"
		k += -2;
#line 290 "sla_syrpvgrw.f"
	    }
#line 291 "sla_syrpvgrw.f"
	}
#line 292 "sla_syrpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 301 "sla_syrpvgrw.f"
    if (upper) {
#line 302 "sla_syrpvgrw.f"
	i__1 = *n;
#line 302 "sla_syrpvgrw.f"
	for (i__ = ncols; i__ <= i__1; ++i__) {
#line 303 "sla_syrpvgrw.f"
	    umax = work[i__];
#line 304 "sla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 305 "sla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 306 "sla_syrpvgrw.f"
		d__1 = amax / umax;
#line 306 "sla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 307 "sla_syrpvgrw.f"
	    }
#line 308 "sla_syrpvgrw.f"
	}
#line 309 "sla_syrpvgrw.f"
    } else {
#line 310 "sla_syrpvgrw.f"
	i__1 = ncols;
#line 310 "sla_syrpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 311 "sla_syrpvgrw.f"
	    umax = work[i__];
#line 312 "sla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 313 "sla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 314 "sla_syrpvgrw.f"
		d__1 = amax / umax;
#line 314 "sla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 315 "sla_syrpvgrw.f"
	    }
#line 316 "sla_syrpvgrw.f"
	}
#line 317 "sla_syrpvgrw.f"
    }
#line 319 "sla_syrpvgrw.f"
    ret_val = rpvgrw;
#line 320 "sla_syrpvgrw.f"
    return ret_val;
} /* sla_syrpvgrw__ */

