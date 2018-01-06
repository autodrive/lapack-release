#line 1 "cla_herpvgrw.f"
/* cla_herpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "cla_herpvgrw.f"
/* > \brief \b CLA_HERPVGRW */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_HERPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_her
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_her
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_her
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, */
/*                                   WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            N, INFO, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ) */
/*       REAL               WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_HERPVGRW computes the reciprocal pivot growth factor */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* > \param[in] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
doublereal cla_herpvgrw__(char *uplo, integer *n, integer *info, 
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


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 163 "cla_herpvgrw.f"
    /* Parameter adjustments */
#line 163 "cla_herpvgrw.f"
    a_dim1 = *lda;
#line 163 "cla_herpvgrw.f"
    a_offset = 1 + a_dim1;
#line 163 "cla_herpvgrw.f"
    a -= a_offset;
#line 163 "cla_herpvgrw.f"
    af_dim1 = *ldaf;
#line 163 "cla_herpvgrw.f"
    af_offset = 1 + af_dim1;
#line 163 "cla_herpvgrw.f"
    af -= af_offset;
#line 163 "cla_herpvgrw.f"
    --ipiv;
#line 163 "cla_herpvgrw.f"
    --work;
#line 163 "cla_herpvgrw.f"

#line 163 "cla_herpvgrw.f"
    /* Function Body */
#line 163 "cla_herpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);
#line 164 "cla_herpvgrw.f"
    if (*info == 0) {
#line 165 "cla_herpvgrw.f"
	if (upper) {
#line 166 "cla_herpvgrw.f"
	    ncols = 1;
#line 167 "cla_herpvgrw.f"
	} else {
#line 168 "cla_herpvgrw.f"
	    ncols = *n;
#line 169 "cla_herpvgrw.f"
	}
#line 170 "cla_herpvgrw.f"
    } else {
#line 171 "cla_herpvgrw.f"
	ncols = *info;
#line 172 "cla_herpvgrw.f"
    }
#line 174 "cla_herpvgrw.f"
    rpvgrw = 1.;
#line 175 "cla_herpvgrw.f"
    i__1 = *n << 1;
#line 175 "cla_herpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 176 "cla_herpvgrw.f"
	work[i__] = 0.;
#line 177 "cla_herpvgrw.f"
    }

/*     Find the max magnitude entry of each column of A.  Compute the max */
/*     for all N columns so we can apply the pivot permutation while */
/*     looping below.  Assume a full factorization is the common case. */

#line 183 "cla_herpvgrw.f"
    if (upper) {
#line 184 "cla_herpvgrw.f"
	i__1 = *n;
#line 184 "cla_herpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 185 "cla_herpvgrw.f"
	    i__2 = j;
#line 185 "cla_herpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 186 "cla_herpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 186 "cla_herpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + i__];
#line 186 "cla_herpvgrw.f"
		work[*n + i__] = max(d__3,d__4);
/* Computing MAX */
#line 187 "cla_herpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 187 "cla_herpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + j];
#line 187 "cla_herpvgrw.f"
		work[*n + j] = max(d__3,d__4);
#line 188 "cla_herpvgrw.f"
	    }
#line 189 "cla_herpvgrw.f"
	}
#line 190 "cla_herpvgrw.f"
    } else {
#line 191 "cla_herpvgrw.f"
	i__1 = *n;
#line 191 "cla_herpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 192 "cla_herpvgrw.f"
	    i__2 = *n;
#line 192 "cla_herpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 193 "cla_herpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 193 "cla_herpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + i__];
#line 193 "cla_herpvgrw.f"
		work[*n + i__] = max(d__3,d__4);
/* Computing MAX */
#line 194 "cla_herpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 194 "cla_herpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + j];
#line 194 "cla_herpvgrw.f"
		work[*n + j] = max(d__3,d__4);
#line 195 "cla_herpvgrw.f"
	    }
#line 196 "cla_herpvgrw.f"
	}
#line 197 "cla_herpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of U or L.  Also */
/*     permute the magnitudes of A above so they're in the same order as */
/*     the factor. */

/*     The iteration orders and permutations were copied from csytrs. */
/*     Calls to SSWAP would be severe overkill. */

#line 206 "cla_herpvgrw.f"
    if (upper) {
#line 207 "cla_herpvgrw.f"
	k = *n;
#line 208 "cla_herpvgrw.f"
	while(k < ncols && k > 0) {
#line 209 "cla_herpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 211 "cla_herpvgrw.f"
		kp = ipiv[k];
#line 212 "cla_herpvgrw.f"
		if (kp != k) {
#line 213 "cla_herpvgrw.f"
		    tmp = work[*n + k];
#line 214 "cla_herpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 215 "cla_herpvgrw.f"
		    work[*n + kp] = tmp;
#line 216 "cla_herpvgrw.f"
		}
#line 217 "cla_herpvgrw.f"
		i__1 = k;
#line 217 "cla_herpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 218 "cla_herpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 218 "cla_herpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 218 "cla_herpvgrw.f"
		    work[k] = max(d__3,d__4);
#line 219 "cla_herpvgrw.f"
		}
#line 220 "cla_herpvgrw.f"
		--k;
#line 221 "cla_herpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 223 "cla_herpvgrw.f"
		kp = -ipiv[k];
#line 224 "cla_herpvgrw.f"
		tmp = work[*n + k - 1];
#line 225 "cla_herpvgrw.f"
		work[*n + k - 1] = work[*n + kp];
#line 226 "cla_herpvgrw.f"
		work[*n + kp] = tmp;
#line 227 "cla_herpvgrw.f"
		i__1 = k - 1;
#line 227 "cla_herpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 228 "cla_herpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 228 "cla_herpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 228 "cla_herpvgrw.f"
		    work[k] = max(d__3,d__4);
/* Computing MAX */
#line 229 "cla_herpvgrw.f"
		    i__2 = i__ + (k - 1) * af_dim1;
#line 229 "cla_herpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + (k - 1) * af_dim1]), abs(d__2)), d__4 = 
			    work[k - 1];
#line 229 "cla_herpvgrw.f"
		    work[k - 1] = max(d__3,d__4);
#line 231 "cla_herpvgrw.f"
		}
/* Computing MAX */
#line 232 "cla_herpvgrw.f"
		i__1 = k + k * af_dim1;
#line 232 "cla_herpvgrw.f"
		d__3 = (d__1 = af[i__1].r, abs(d__1)) + (d__2 = d_imag(&af[k 
			+ k * af_dim1]), abs(d__2)), d__4 = work[k];
#line 232 "cla_herpvgrw.f"
		work[k] = max(d__3,d__4);
#line 233 "cla_herpvgrw.f"
		k += -2;
#line 234 "cla_herpvgrw.f"
	    }
#line 235 "cla_herpvgrw.f"
	}
#line 236 "cla_herpvgrw.f"
	k = ncols;
#line 237 "cla_herpvgrw.f"
	while(k <= *n) {
#line 238 "cla_herpvgrw.f"
	    if (ipiv[k] > 0) {
#line 239 "cla_herpvgrw.f"
		kp = ipiv[k];
#line 240 "cla_herpvgrw.f"
		if (kp != k) {
#line 241 "cla_herpvgrw.f"
		    tmp = work[*n + k];
#line 242 "cla_herpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 243 "cla_herpvgrw.f"
		    work[*n + kp] = tmp;
#line 244 "cla_herpvgrw.f"
		}
#line 245 "cla_herpvgrw.f"
		++k;
#line 246 "cla_herpvgrw.f"
	    } else {
#line 247 "cla_herpvgrw.f"
		kp = -ipiv[k];
#line 248 "cla_herpvgrw.f"
		tmp = work[*n + k];
#line 249 "cla_herpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 250 "cla_herpvgrw.f"
		work[*n + kp] = tmp;
#line 251 "cla_herpvgrw.f"
		k += 2;
#line 252 "cla_herpvgrw.f"
	    }
#line 253 "cla_herpvgrw.f"
	}
#line 254 "cla_herpvgrw.f"
    } else {
#line 255 "cla_herpvgrw.f"
	k = 1;
#line 256 "cla_herpvgrw.f"
	while(k <= ncols) {
#line 257 "cla_herpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 259 "cla_herpvgrw.f"
		kp = ipiv[k];
#line 260 "cla_herpvgrw.f"
		if (kp != k) {
#line 261 "cla_herpvgrw.f"
		    tmp = work[*n + k];
#line 262 "cla_herpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 263 "cla_herpvgrw.f"
		    work[*n + kp] = tmp;
#line 264 "cla_herpvgrw.f"
		}
#line 265 "cla_herpvgrw.f"
		i__1 = *n;
#line 265 "cla_herpvgrw.f"
		for (i__ = k; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 266 "cla_herpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 266 "cla_herpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 266 "cla_herpvgrw.f"
		    work[k] = max(d__3,d__4);
#line 267 "cla_herpvgrw.f"
		}
#line 268 "cla_herpvgrw.f"
		++k;
#line 269 "cla_herpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 271 "cla_herpvgrw.f"
		kp = -ipiv[k];
#line 272 "cla_herpvgrw.f"
		tmp = work[*n + k + 1];
#line 273 "cla_herpvgrw.f"
		work[*n + k + 1] = work[*n + kp];
#line 274 "cla_herpvgrw.f"
		work[*n + kp] = tmp;
#line 275 "cla_herpvgrw.f"
		i__1 = *n;
#line 275 "cla_herpvgrw.f"
		for (i__ = k + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 276 "cla_herpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 276 "cla_herpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 276 "cla_herpvgrw.f"
		    work[k] = max(d__3,d__4);
/* Computing MAX */
#line 277 "cla_herpvgrw.f"
		    i__2 = i__ + (k + 1) * af_dim1;
#line 277 "cla_herpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + (k + 1) * af_dim1]), abs(d__2)), d__4 = 
			    work[k + 1];
#line 277 "cla_herpvgrw.f"
		    work[k + 1] = max(d__3,d__4);
#line 279 "cla_herpvgrw.f"
		}
/* Computing MAX */
#line 280 "cla_herpvgrw.f"
		i__1 = k + k * af_dim1;
#line 280 "cla_herpvgrw.f"
		d__3 = (d__1 = af[i__1].r, abs(d__1)) + (d__2 = d_imag(&af[k 
			+ k * af_dim1]), abs(d__2)), d__4 = work[k];
#line 280 "cla_herpvgrw.f"
		work[k] = max(d__3,d__4);
#line 281 "cla_herpvgrw.f"
		k += 2;
#line 282 "cla_herpvgrw.f"
	    }
#line 283 "cla_herpvgrw.f"
	}
#line 284 "cla_herpvgrw.f"
	k = ncols;
#line 285 "cla_herpvgrw.f"
	while(k >= 1) {
#line 286 "cla_herpvgrw.f"
	    if (ipiv[k] > 0) {
#line 287 "cla_herpvgrw.f"
		kp = ipiv[k];
#line 288 "cla_herpvgrw.f"
		if (kp != k) {
#line 289 "cla_herpvgrw.f"
		    tmp = work[*n + k];
#line 290 "cla_herpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 291 "cla_herpvgrw.f"
		    work[*n + kp] = tmp;
#line 292 "cla_herpvgrw.f"
		}
#line 293 "cla_herpvgrw.f"
		--k;
#line 294 "cla_herpvgrw.f"
	    } else {
#line 295 "cla_herpvgrw.f"
		kp = -ipiv[k];
#line 296 "cla_herpvgrw.f"
		tmp = work[*n + k];
#line 297 "cla_herpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 298 "cla_herpvgrw.f"
		work[*n + kp] = tmp;
#line 299 "cla_herpvgrw.f"
		k += -2;
#line 300 "cla_herpvgrw.f"
	    }
#line 301 "cla_herpvgrw.f"
	}
#line 302 "cla_herpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 311 "cla_herpvgrw.f"
    if (upper) {
#line 312 "cla_herpvgrw.f"
	i__1 = *n;
#line 312 "cla_herpvgrw.f"
	for (i__ = ncols; i__ <= i__1; ++i__) {
#line 313 "cla_herpvgrw.f"
	    umax = work[i__];
#line 314 "cla_herpvgrw.f"
	    amax = work[*n + i__];
#line 315 "cla_herpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 316 "cla_herpvgrw.f"
		d__1 = amax / umax;
#line 316 "cla_herpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 317 "cla_herpvgrw.f"
	    }
#line 318 "cla_herpvgrw.f"
	}
#line 319 "cla_herpvgrw.f"
    } else {
#line 320 "cla_herpvgrw.f"
	i__1 = ncols;
#line 320 "cla_herpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 321 "cla_herpvgrw.f"
	    umax = work[i__];
#line 322 "cla_herpvgrw.f"
	    amax = work[*n + i__];
#line 323 "cla_herpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 324 "cla_herpvgrw.f"
		d__1 = amax / umax;
#line 324 "cla_herpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 325 "cla_herpvgrw.f"
	    }
#line 326 "cla_herpvgrw.f"
	}
#line 327 "cla_herpvgrw.f"
    }
#line 329 "cla_herpvgrw.f"
    ret_val = rpvgrw;
#line 330 "cla_herpvgrw.f"
    return ret_val;
} /* cla_herpvgrw__ */

