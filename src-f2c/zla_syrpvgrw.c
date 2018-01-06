#line 1 "zla_syrpvgrw.f"
/* zla_syrpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "zla_syrpvgrw.f"
/* > \brief \b ZLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefi
nite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_SYRPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_syr
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_syr
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_syr
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, */
/*                                               LDAF, IPIV, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            N, INFO, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ) */
/*       DOUBLE PRECISION   WORK( * ) */
/*       INTEGER            IPIV( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > ZLA_SYRPVGRW computes the reciprocal pivot growth factor */
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
/* >     The value of INFO returned from ZSYTRF, .i.e., the pivot in */
/* >     column INFO is exactly 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by ZSYTRF. */
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
/* >     as determined by ZSYTRF. */
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

/* > \date November 2015 */

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
doublereal zla_syrpvgrw__(char *uplo, integer *n, integer *info, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublereal *work, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, k, kp;
    static doublereal tmp, amax, umax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ncols;
    static logical upper;
    static doublereal rpvgrw;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 164 "zla_syrpvgrw.f"
    /* Parameter adjustments */
#line 164 "zla_syrpvgrw.f"
    a_dim1 = *lda;
#line 164 "zla_syrpvgrw.f"
    a_offset = 1 + a_dim1;
#line 164 "zla_syrpvgrw.f"
    a -= a_offset;
#line 164 "zla_syrpvgrw.f"
    af_dim1 = *ldaf;
#line 164 "zla_syrpvgrw.f"
    af_offset = 1 + af_dim1;
#line 164 "zla_syrpvgrw.f"
    af -= af_offset;
#line 164 "zla_syrpvgrw.f"
    --ipiv;
#line 164 "zla_syrpvgrw.f"
    --work;
#line 164 "zla_syrpvgrw.f"

#line 164 "zla_syrpvgrw.f"
    /* Function Body */
#line 164 "zla_syrpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);
#line 165 "zla_syrpvgrw.f"
    if (*info == 0) {
#line 166 "zla_syrpvgrw.f"
	if (upper) {
#line 167 "zla_syrpvgrw.f"
	    ncols = 1;
#line 168 "zla_syrpvgrw.f"
	} else {
#line 169 "zla_syrpvgrw.f"
	    ncols = *n;
#line 170 "zla_syrpvgrw.f"
	}
#line 171 "zla_syrpvgrw.f"
    } else {
#line 172 "zla_syrpvgrw.f"
	ncols = *info;
#line 173 "zla_syrpvgrw.f"
    }
#line 175 "zla_syrpvgrw.f"
    rpvgrw = 1.;
#line 176 "zla_syrpvgrw.f"
    i__1 = *n << 1;
#line 176 "zla_syrpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 177 "zla_syrpvgrw.f"
	work[i__] = 0.;
#line 178 "zla_syrpvgrw.f"
    }

/*     Find the max magnitude entry of each column of A.  Compute the max */
/*     for all N columns so we can apply the pivot permutation while */
/*     looping below.  Assume a full factorization is the common case. */

#line 184 "zla_syrpvgrw.f"
    if (upper) {
#line 185 "zla_syrpvgrw.f"
	i__1 = *n;
#line 185 "zla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 186 "zla_syrpvgrw.f"
	    i__2 = j;
#line 186 "zla_syrpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 187 "zla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 187 "zla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + i__];
#line 187 "zla_syrpvgrw.f"
		work[*n + i__] = max(d__3,d__4);
/* Computing MAX */
#line 188 "zla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 188 "zla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + j];
#line 188 "zla_syrpvgrw.f"
		work[*n + j] = max(d__3,d__4);
#line 189 "zla_syrpvgrw.f"
	    }
#line 190 "zla_syrpvgrw.f"
	}
#line 191 "zla_syrpvgrw.f"
    } else {
#line 192 "zla_syrpvgrw.f"
	i__1 = *n;
#line 192 "zla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 193 "zla_syrpvgrw.f"
	    i__2 = *n;
#line 193 "zla_syrpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 194 "zla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 194 "zla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + i__];
#line 194 "zla_syrpvgrw.f"
		work[*n + i__] = max(d__3,d__4);
/* Computing MAX */
#line 195 "zla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 195 "zla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + j];
#line 195 "zla_syrpvgrw.f"
		work[*n + j] = max(d__3,d__4);
#line 196 "zla_syrpvgrw.f"
	    }
#line 197 "zla_syrpvgrw.f"
	}
#line 198 "zla_syrpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of U or L.  Also */
/*     permute the magnitudes of A above so they're in the same order as */
/*     the factor. */

/*     The iteration orders and permutations were copied from zsytrs. */
/*     Calls to SSWAP would be severe overkill. */

#line 207 "zla_syrpvgrw.f"
    if (upper) {
#line 208 "zla_syrpvgrw.f"
	k = *n;
#line 209 "zla_syrpvgrw.f"
	while(k < ncols && k > 0) {
#line 210 "zla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 212 "zla_syrpvgrw.f"
		kp = ipiv[k];
#line 213 "zla_syrpvgrw.f"
		if (kp != k) {
#line 214 "zla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 215 "zla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 216 "zla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 217 "zla_syrpvgrw.f"
		}
#line 218 "zla_syrpvgrw.f"
		i__1 = k;
#line 218 "zla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 219 "zla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 219 "zla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 219 "zla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
#line 220 "zla_syrpvgrw.f"
		}
#line 221 "zla_syrpvgrw.f"
		--k;
#line 222 "zla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 224 "zla_syrpvgrw.f"
		kp = -ipiv[k];
#line 225 "zla_syrpvgrw.f"
		tmp = work[*n + k - 1];
#line 226 "zla_syrpvgrw.f"
		work[*n + k - 1] = work[*n + kp];
#line 227 "zla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 228 "zla_syrpvgrw.f"
		i__1 = k - 1;
#line 228 "zla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 229 "zla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 229 "zla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 229 "zla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
/* Computing MAX */
#line 230 "zla_syrpvgrw.f"
		    i__2 = i__ + (k - 1) * af_dim1;
#line 230 "zla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + (k - 1) * af_dim1]), abs(d__2)), d__4 = 
			    work[k - 1];
#line 230 "zla_syrpvgrw.f"
		    work[k - 1] = max(d__3,d__4);
#line 232 "zla_syrpvgrw.f"
		}
/* Computing MAX */
#line 233 "zla_syrpvgrw.f"
		i__1 = k + k * af_dim1;
#line 233 "zla_syrpvgrw.f"
		d__3 = (d__1 = af[i__1].r, abs(d__1)) + (d__2 = d_imag(&af[k 
			+ k * af_dim1]), abs(d__2)), d__4 = work[k];
#line 233 "zla_syrpvgrw.f"
		work[k] = max(d__3,d__4);
#line 234 "zla_syrpvgrw.f"
		k += -2;
#line 235 "zla_syrpvgrw.f"
	    }
#line 236 "zla_syrpvgrw.f"
	}
#line 237 "zla_syrpvgrw.f"
	k = ncols;
#line 238 "zla_syrpvgrw.f"
	while(k <= *n) {
#line 239 "zla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 240 "zla_syrpvgrw.f"
		kp = ipiv[k];
#line 241 "zla_syrpvgrw.f"
		if (kp != k) {
#line 242 "zla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 243 "zla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 244 "zla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 245 "zla_syrpvgrw.f"
		}
#line 246 "zla_syrpvgrw.f"
		++k;
#line 247 "zla_syrpvgrw.f"
	    } else {
#line 248 "zla_syrpvgrw.f"
		kp = -ipiv[k];
#line 249 "zla_syrpvgrw.f"
		tmp = work[*n + k];
#line 250 "zla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 251 "zla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 252 "zla_syrpvgrw.f"
		k += 2;
#line 253 "zla_syrpvgrw.f"
	    }
#line 254 "zla_syrpvgrw.f"
	}
#line 255 "zla_syrpvgrw.f"
    } else {
#line 256 "zla_syrpvgrw.f"
	k = 1;
#line 257 "zla_syrpvgrw.f"
	while(k <= ncols) {
#line 258 "zla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 260 "zla_syrpvgrw.f"
		kp = ipiv[k];
#line 261 "zla_syrpvgrw.f"
		if (kp != k) {
#line 262 "zla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 263 "zla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 264 "zla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 265 "zla_syrpvgrw.f"
		}
#line 266 "zla_syrpvgrw.f"
		i__1 = *n;
#line 266 "zla_syrpvgrw.f"
		for (i__ = k; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 267 "zla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 267 "zla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 267 "zla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
#line 268 "zla_syrpvgrw.f"
		}
#line 269 "zla_syrpvgrw.f"
		++k;
#line 270 "zla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 272 "zla_syrpvgrw.f"
		kp = -ipiv[k];
#line 273 "zla_syrpvgrw.f"
		tmp = work[*n + k + 1];
#line 274 "zla_syrpvgrw.f"
		work[*n + k + 1] = work[*n + kp];
#line 275 "zla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 276 "zla_syrpvgrw.f"
		i__1 = *n;
#line 276 "zla_syrpvgrw.f"
		for (i__ = k + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 277 "zla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 277 "zla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 277 "zla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
/* Computing MAX */
#line 278 "zla_syrpvgrw.f"
		    i__2 = i__ + (k + 1) * af_dim1;
#line 278 "zla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + (k + 1) * af_dim1]), abs(d__2)), d__4 = 
			    work[k + 1];
#line 278 "zla_syrpvgrw.f"
		    work[k + 1] = max(d__3,d__4);
#line 280 "zla_syrpvgrw.f"
		}
/* Computing MAX */
#line 281 "zla_syrpvgrw.f"
		i__1 = k + k * af_dim1;
#line 281 "zla_syrpvgrw.f"
		d__3 = (d__1 = af[i__1].r, abs(d__1)) + (d__2 = d_imag(&af[k 
			+ k * af_dim1]), abs(d__2)), d__4 = work[k];
#line 281 "zla_syrpvgrw.f"
		work[k] = max(d__3,d__4);
#line 282 "zla_syrpvgrw.f"
		k += 2;
#line 283 "zla_syrpvgrw.f"
	    }
#line 284 "zla_syrpvgrw.f"
	}
#line 285 "zla_syrpvgrw.f"
	k = ncols;
#line 286 "zla_syrpvgrw.f"
	while(k >= 1) {
#line 287 "zla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 288 "zla_syrpvgrw.f"
		kp = ipiv[k];
#line 289 "zla_syrpvgrw.f"
		if (kp != k) {
#line 290 "zla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 291 "zla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 292 "zla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 293 "zla_syrpvgrw.f"
		}
#line 294 "zla_syrpvgrw.f"
		--k;
#line 295 "zla_syrpvgrw.f"
	    } else {
#line 296 "zla_syrpvgrw.f"
		kp = -ipiv[k];
#line 297 "zla_syrpvgrw.f"
		tmp = work[*n + k];
#line 298 "zla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 299 "zla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 300 "zla_syrpvgrw.f"
		k += -2;
#line 301 "zla_syrpvgrw.f"
	    }
#line 302 "zla_syrpvgrw.f"
	}
#line 303 "zla_syrpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 312 "zla_syrpvgrw.f"
    if (upper) {
#line 313 "zla_syrpvgrw.f"
	i__1 = *n;
#line 313 "zla_syrpvgrw.f"
	for (i__ = ncols; i__ <= i__1; ++i__) {
#line 314 "zla_syrpvgrw.f"
	    umax = work[i__];
#line 315 "zla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 316 "zla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 317 "zla_syrpvgrw.f"
		d__1 = amax / umax;
#line 317 "zla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 318 "zla_syrpvgrw.f"
	    }
#line 319 "zla_syrpvgrw.f"
	}
#line 320 "zla_syrpvgrw.f"
    } else {
#line 321 "zla_syrpvgrw.f"
	i__1 = ncols;
#line 321 "zla_syrpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 322 "zla_syrpvgrw.f"
	    umax = work[i__];
#line 323 "zla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 324 "zla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 325 "zla_syrpvgrw.f"
		d__1 = amax / umax;
#line 325 "zla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 326 "zla_syrpvgrw.f"
	    }
#line 327 "zla_syrpvgrw.f"
	}
#line 328 "zla_syrpvgrw.f"
    }
#line 330 "zla_syrpvgrw.f"
    ret_val = rpvgrw;
#line 331 "zla_syrpvgrw.f"
    return ret_val;
} /* zla_syrpvgrw__ */

