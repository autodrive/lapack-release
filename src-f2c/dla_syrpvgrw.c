#line 1 "dla_syrpvgrw.f"
/* dla_syrpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "dla_syrpvgrw.f"
/* > \brief \b DLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefi
nite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_SYRPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_syr
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_syr
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_syr
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, */
/*                                               LDAF, IPIV, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            N, INFO, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > DLA_SYRPVGRW computes the reciprocal pivot growth factor */
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
/* >     The value of INFO returned from DSYTRF, .i.e., the pivot in */
/* >     column INFO is exactly 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          AF is DOUBLE PRECISION array, dimension (LDAF,N) */
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by DSYTRF. */
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
/* >     as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
doublereal dla_syrpvgrw__(char *uplo, integer *n, integer *info, doublereal *
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 155 "dla_syrpvgrw.f"
    /* Parameter adjustments */
#line 155 "dla_syrpvgrw.f"
    a_dim1 = *lda;
#line 155 "dla_syrpvgrw.f"
    a_offset = 1 + a_dim1;
#line 155 "dla_syrpvgrw.f"
    a -= a_offset;
#line 155 "dla_syrpvgrw.f"
    af_dim1 = *ldaf;
#line 155 "dla_syrpvgrw.f"
    af_offset = 1 + af_dim1;
#line 155 "dla_syrpvgrw.f"
    af -= af_offset;
#line 155 "dla_syrpvgrw.f"
    --ipiv;
#line 155 "dla_syrpvgrw.f"
    --work;
#line 155 "dla_syrpvgrw.f"

#line 155 "dla_syrpvgrw.f"
    /* Function Body */
#line 155 "dla_syrpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);
#line 156 "dla_syrpvgrw.f"
    if (*info == 0) {
#line 157 "dla_syrpvgrw.f"
	if (upper) {
#line 158 "dla_syrpvgrw.f"
	    ncols = 1;
#line 159 "dla_syrpvgrw.f"
	} else {
#line 160 "dla_syrpvgrw.f"
	    ncols = *n;
#line 161 "dla_syrpvgrw.f"
	}
#line 162 "dla_syrpvgrw.f"
    } else {
#line 163 "dla_syrpvgrw.f"
	ncols = *info;
#line 164 "dla_syrpvgrw.f"
    }
#line 166 "dla_syrpvgrw.f"
    rpvgrw = 1.;
#line 167 "dla_syrpvgrw.f"
    i__1 = *n << 1;
#line 167 "dla_syrpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 168 "dla_syrpvgrw.f"
	work[i__] = 0.;
#line 169 "dla_syrpvgrw.f"
    }

/*     Find the max magnitude entry of each column of A.  Compute the max */
/*     for all N columns so we can apply the pivot permutation while */
/*     looping below.  Assume a full factorization is the common case. */

#line 175 "dla_syrpvgrw.f"
    if (upper) {
#line 176 "dla_syrpvgrw.f"
	i__1 = *n;
#line 176 "dla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 177 "dla_syrpvgrw.f"
	    i__2 = j;
#line 177 "dla_syrpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 178 "dla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + i__];
#line 178 "dla_syrpvgrw.f"
		work[*n + i__] = max(d__2,d__3);
/* Computing MAX */
#line 179 "dla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + j];
#line 179 "dla_syrpvgrw.f"
		work[*n + j] = max(d__2,d__3);
#line 180 "dla_syrpvgrw.f"
	    }
#line 181 "dla_syrpvgrw.f"
	}
#line 182 "dla_syrpvgrw.f"
    } else {
#line 183 "dla_syrpvgrw.f"
	i__1 = *n;
#line 183 "dla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "dla_syrpvgrw.f"
	    i__2 = *n;
#line 184 "dla_syrpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 185 "dla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + i__];
#line 185 "dla_syrpvgrw.f"
		work[*n + i__] = max(d__2,d__3);
/* Computing MAX */
#line 186 "dla_syrpvgrw.f"
		d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)), d__3 = work[*
			n + j];
#line 186 "dla_syrpvgrw.f"
		work[*n + j] = max(d__2,d__3);
#line 187 "dla_syrpvgrw.f"
	    }
#line 188 "dla_syrpvgrw.f"
	}
#line 189 "dla_syrpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of U or L.  Also */
/*     permute the magnitudes of A above so they're in the same order as */
/*     the factor. */

/*     The iteration orders and permutations were copied from dsytrs. */
/*     Calls to SSWAP would be severe overkill. */

#line 198 "dla_syrpvgrw.f"
    if (upper) {
#line 199 "dla_syrpvgrw.f"
	k = *n;
#line 200 "dla_syrpvgrw.f"
	while(k < ncols && k > 0) {
#line 201 "dla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 203 "dla_syrpvgrw.f"
		kp = ipiv[k];
#line 204 "dla_syrpvgrw.f"
		if (kp != k) {
#line 205 "dla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 206 "dla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 207 "dla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 208 "dla_syrpvgrw.f"
		}
#line 209 "dla_syrpvgrw.f"
		i__1 = k;
#line 209 "dla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 210 "dla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 210 "dla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
#line 211 "dla_syrpvgrw.f"
		}
#line 212 "dla_syrpvgrw.f"
		--k;
#line 213 "dla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 215 "dla_syrpvgrw.f"
		kp = -ipiv[k];
#line 216 "dla_syrpvgrw.f"
		tmp = work[*n + k - 1];
#line 217 "dla_syrpvgrw.f"
		work[*n + k - 1] = work[*n + kp];
#line 218 "dla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 219 "dla_syrpvgrw.f"
		i__1 = k - 1;
#line 219 "dla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 220 "dla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 220 "dla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
/* Computing MAX */
#line 221 "dla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + (k - 1) * af_dim1], abs(d__1)), 
			    d__3 = work[k - 1];
#line 221 "dla_syrpvgrw.f"
		    work[k - 1] = max(d__2,d__3);
#line 222 "dla_syrpvgrw.f"
		}
/* Computing MAX */
#line 223 "dla_syrpvgrw.f"
		d__2 = (d__1 = af[k + k * af_dim1], abs(d__1)), d__3 = work[k]
			;
#line 223 "dla_syrpvgrw.f"
		work[k] = max(d__2,d__3);
#line 224 "dla_syrpvgrw.f"
		k += -2;
#line 225 "dla_syrpvgrw.f"
	    }
#line 226 "dla_syrpvgrw.f"
	}
#line 227 "dla_syrpvgrw.f"
	k = ncols;
#line 228 "dla_syrpvgrw.f"
	while(k <= *n) {
#line 229 "dla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 230 "dla_syrpvgrw.f"
		kp = ipiv[k];
#line 231 "dla_syrpvgrw.f"
		if (kp != k) {
#line 232 "dla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 233 "dla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 234 "dla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 235 "dla_syrpvgrw.f"
		}
#line 236 "dla_syrpvgrw.f"
		++k;
#line 237 "dla_syrpvgrw.f"
	    } else {
#line 238 "dla_syrpvgrw.f"
		kp = -ipiv[k];
#line 239 "dla_syrpvgrw.f"
		tmp = work[*n + k];
#line 240 "dla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 241 "dla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 242 "dla_syrpvgrw.f"
		k += 2;
#line 243 "dla_syrpvgrw.f"
	    }
#line 244 "dla_syrpvgrw.f"
	}
#line 245 "dla_syrpvgrw.f"
    } else {
#line 246 "dla_syrpvgrw.f"
	k = 1;
#line 247 "dla_syrpvgrw.f"
	while(k <= ncols) {
#line 248 "dla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 250 "dla_syrpvgrw.f"
		kp = ipiv[k];
#line 251 "dla_syrpvgrw.f"
		if (kp != k) {
#line 252 "dla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 253 "dla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 254 "dla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 255 "dla_syrpvgrw.f"
		}
#line 256 "dla_syrpvgrw.f"
		i__1 = *n;
#line 256 "dla_syrpvgrw.f"
		for (i__ = k; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 257 "dla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 257 "dla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
#line 258 "dla_syrpvgrw.f"
		}
#line 259 "dla_syrpvgrw.f"
		++k;
#line 260 "dla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 262 "dla_syrpvgrw.f"
		kp = -ipiv[k];
#line 263 "dla_syrpvgrw.f"
		tmp = work[*n + k + 1];
#line 264 "dla_syrpvgrw.f"
		work[*n + k + 1] = work[*n + kp];
#line 265 "dla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 266 "dla_syrpvgrw.f"
		i__1 = *n;
#line 266 "dla_syrpvgrw.f"
		for (i__ = k + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 267 "dla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + k * af_dim1], abs(d__1)), d__3 = 
			    work[k];
#line 267 "dla_syrpvgrw.f"
		    work[k] = max(d__2,d__3);
/* Computing MAX */
#line 268 "dla_syrpvgrw.f"
		    d__2 = (d__1 = af[i__ + (k + 1) * af_dim1], abs(d__1)), 
			    d__3 = work[k + 1];
#line 268 "dla_syrpvgrw.f"
		    work[k + 1] = max(d__2,d__3);
#line 269 "dla_syrpvgrw.f"
		}
/* Computing MAX */
#line 270 "dla_syrpvgrw.f"
		d__2 = (d__1 = af[k + k * af_dim1], abs(d__1)), d__3 = work[k]
			;
#line 270 "dla_syrpvgrw.f"
		work[k] = max(d__2,d__3);
#line 271 "dla_syrpvgrw.f"
		k += 2;
#line 272 "dla_syrpvgrw.f"
	    }
#line 273 "dla_syrpvgrw.f"
	}
#line 274 "dla_syrpvgrw.f"
	k = ncols;
#line 275 "dla_syrpvgrw.f"
	while(k >= 1) {
#line 276 "dla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 277 "dla_syrpvgrw.f"
		kp = ipiv[k];
#line 278 "dla_syrpvgrw.f"
		if (kp != k) {
#line 279 "dla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 280 "dla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 281 "dla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 282 "dla_syrpvgrw.f"
		}
#line 283 "dla_syrpvgrw.f"
		--k;
#line 284 "dla_syrpvgrw.f"
	    } else {
#line 285 "dla_syrpvgrw.f"
		kp = -ipiv[k];
#line 286 "dla_syrpvgrw.f"
		tmp = work[*n + k];
#line 287 "dla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 288 "dla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 289 "dla_syrpvgrw.f"
		k += -2;
#line 290 "dla_syrpvgrw.f"
	    }
#line 291 "dla_syrpvgrw.f"
	}
#line 292 "dla_syrpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 301 "dla_syrpvgrw.f"
    if (upper) {
#line 302 "dla_syrpvgrw.f"
	i__1 = *n;
#line 302 "dla_syrpvgrw.f"
	for (i__ = ncols; i__ <= i__1; ++i__) {
#line 303 "dla_syrpvgrw.f"
	    umax = work[i__];
#line 304 "dla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 305 "dla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 306 "dla_syrpvgrw.f"
		d__1 = amax / umax;
#line 306 "dla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 307 "dla_syrpvgrw.f"
	    }
#line 308 "dla_syrpvgrw.f"
	}
#line 309 "dla_syrpvgrw.f"
    } else {
#line 310 "dla_syrpvgrw.f"
	i__1 = ncols;
#line 310 "dla_syrpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 311 "dla_syrpvgrw.f"
	    umax = work[i__];
#line 312 "dla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 313 "dla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 314 "dla_syrpvgrw.f"
		d__1 = amax / umax;
#line 314 "dla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 315 "dla_syrpvgrw.f"
	    }
#line 316 "dla_syrpvgrw.f"
	}
#line 317 "dla_syrpvgrw.f"
    }
#line 319 "dla_syrpvgrw.f"
    ret_val = rpvgrw;
#line 320 "dla_syrpvgrw.f"
    return ret_val;
} /* dla_syrpvgrw__ */

