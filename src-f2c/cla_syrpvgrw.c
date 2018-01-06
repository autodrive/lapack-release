#line 1 "cla_syrpvgrw.f"
/* cla_syrpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "cla_syrpvgrw.f"
/* > \brief \b CLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefi
nite matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_SYRPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syr
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syr
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syr
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, */
/*                                   WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*1        UPLO */
/*       INTEGER            N, INFO, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ) */
/*       REAL               WORK( * ) */
/*       INTEGER            IPIV( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_SYRPVGRW computes the reciprocal pivot growth factor */
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
/* >     The value of INFO returned from CSYTRF, .i.e., the pivot in */
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

/* > \date November 2015 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
doublereal cla_syrpvgrw__(char *uplo, integer *n, integer *info, 
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

#line 164 "cla_syrpvgrw.f"
    /* Parameter adjustments */
#line 164 "cla_syrpvgrw.f"
    a_dim1 = *lda;
#line 164 "cla_syrpvgrw.f"
    a_offset = 1 + a_dim1;
#line 164 "cla_syrpvgrw.f"
    a -= a_offset;
#line 164 "cla_syrpvgrw.f"
    af_dim1 = *ldaf;
#line 164 "cla_syrpvgrw.f"
    af_offset = 1 + af_dim1;
#line 164 "cla_syrpvgrw.f"
    af -= af_offset;
#line 164 "cla_syrpvgrw.f"
    --ipiv;
#line 164 "cla_syrpvgrw.f"
    --work;
#line 164 "cla_syrpvgrw.f"

#line 164 "cla_syrpvgrw.f"
    /* Function Body */
#line 164 "cla_syrpvgrw.f"
    upper = lsame_("Upper", uplo, (ftnlen)5, (ftnlen)1);
#line 165 "cla_syrpvgrw.f"
    if (*info == 0) {
#line 166 "cla_syrpvgrw.f"
	if (upper) {
#line 167 "cla_syrpvgrw.f"
	    ncols = 1;
#line 168 "cla_syrpvgrw.f"
	} else {
#line 169 "cla_syrpvgrw.f"
	    ncols = *n;
#line 170 "cla_syrpvgrw.f"
	}
#line 171 "cla_syrpvgrw.f"
    } else {
#line 172 "cla_syrpvgrw.f"
	ncols = *info;
#line 173 "cla_syrpvgrw.f"
    }
#line 175 "cla_syrpvgrw.f"
    rpvgrw = 1.;
#line 176 "cla_syrpvgrw.f"
    i__1 = *n << 1;
#line 176 "cla_syrpvgrw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 177 "cla_syrpvgrw.f"
	work[i__] = 0.;
#line 178 "cla_syrpvgrw.f"
    }

/*     Find the max magnitude entry of each column of A.  Compute the max */
/*     for all N columns so we can apply the pivot permutation while */
/*     looping below.  Assume a full factorization is the common case. */

#line 184 "cla_syrpvgrw.f"
    if (upper) {
#line 185 "cla_syrpvgrw.f"
	i__1 = *n;
#line 185 "cla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 186 "cla_syrpvgrw.f"
	    i__2 = j;
#line 186 "cla_syrpvgrw.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 187 "cla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 187 "cla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + i__];
#line 187 "cla_syrpvgrw.f"
		work[*n + i__] = max(d__3,d__4);
/* Computing MAX */
#line 188 "cla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 188 "cla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + j];
#line 188 "cla_syrpvgrw.f"
		work[*n + j] = max(d__3,d__4);
#line 189 "cla_syrpvgrw.f"
	    }
#line 190 "cla_syrpvgrw.f"
	}
#line 191 "cla_syrpvgrw.f"
    } else {
#line 192 "cla_syrpvgrw.f"
	i__1 = *n;
#line 192 "cla_syrpvgrw.f"
	for (j = 1; j <= i__1; ++j) {
#line 193 "cla_syrpvgrw.f"
	    i__2 = *n;
#line 193 "cla_syrpvgrw.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 194 "cla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 194 "cla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + i__];
#line 194 "cla_syrpvgrw.f"
		work[*n + i__] = max(d__3,d__4);
/* Computing MAX */
#line 195 "cla_syrpvgrw.f"
		i__3 = i__ + j * a_dim1;
#line 195 "cla_syrpvgrw.f"
		d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			+ j * a_dim1]), abs(d__2)), d__4 = work[*n + j];
#line 195 "cla_syrpvgrw.f"
		work[*n + j] = max(d__3,d__4);
#line 196 "cla_syrpvgrw.f"
	    }
#line 197 "cla_syrpvgrw.f"
	}
#line 198 "cla_syrpvgrw.f"
    }

/*     Now find the max magnitude entry of each column of U or L.  Also */
/*     permute the magnitudes of A above so they're in the same order as */
/*     the factor. */

/*     The iteration orders and permutations were copied from csytrs. */
/*     Calls to SSWAP would be severe overkill. */

#line 207 "cla_syrpvgrw.f"
    if (upper) {
#line 208 "cla_syrpvgrw.f"
	k = *n;
#line 209 "cla_syrpvgrw.f"
	while(k < ncols && k > 0) {
#line 210 "cla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 212 "cla_syrpvgrw.f"
		kp = ipiv[k];
#line 213 "cla_syrpvgrw.f"
		if (kp != k) {
#line 214 "cla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 215 "cla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 216 "cla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 217 "cla_syrpvgrw.f"
		}
#line 218 "cla_syrpvgrw.f"
		i__1 = k;
#line 218 "cla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 219 "cla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 219 "cla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 219 "cla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
#line 220 "cla_syrpvgrw.f"
		}
#line 221 "cla_syrpvgrw.f"
		--k;
#line 222 "cla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 224 "cla_syrpvgrw.f"
		kp = -ipiv[k];
#line 225 "cla_syrpvgrw.f"
		tmp = work[*n + k - 1];
#line 226 "cla_syrpvgrw.f"
		work[*n + k - 1] = work[*n + kp];
#line 227 "cla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 228 "cla_syrpvgrw.f"
		i__1 = k - 1;
#line 228 "cla_syrpvgrw.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 229 "cla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 229 "cla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 229 "cla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
/* Computing MAX */
#line 230 "cla_syrpvgrw.f"
		    i__2 = i__ + (k - 1) * af_dim1;
#line 230 "cla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + (k - 1) * af_dim1]), abs(d__2)), d__4 = 
			    work[k - 1];
#line 230 "cla_syrpvgrw.f"
		    work[k - 1] = max(d__3,d__4);
#line 232 "cla_syrpvgrw.f"
		}
/* Computing MAX */
#line 233 "cla_syrpvgrw.f"
		i__1 = k + k * af_dim1;
#line 233 "cla_syrpvgrw.f"
		d__3 = (d__1 = af[i__1].r, abs(d__1)) + (d__2 = d_imag(&af[k 
			+ k * af_dim1]), abs(d__2)), d__4 = work[k];
#line 233 "cla_syrpvgrw.f"
		work[k] = max(d__3,d__4);
#line 234 "cla_syrpvgrw.f"
		k += -2;
#line 235 "cla_syrpvgrw.f"
	    }
#line 236 "cla_syrpvgrw.f"
	}
#line 237 "cla_syrpvgrw.f"
	k = ncols;
#line 238 "cla_syrpvgrw.f"
	while(k <= *n) {
#line 239 "cla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 240 "cla_syrpvgrw.f"
		kp = ipiv[k];
#line 241 "cla_syrpvgrw.f"
		if (kp != k) {
#line 242 "cla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 243 "cla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 244 "cla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 245 "cla_syrpvgrw.f"
		}
#line 246 "cla_syrpvgrw.f"
		++k;
#line 247 "cla_syrpvgrw.f"
	    } else {
#line 248 "cla_syrpvgrw.f"
		kp = -ipiv[k];
#line 249 "cla_syrpvgrw.f"
		tmp = work[*n + k];
#line 250 "cla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 251 "cla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 252 "cla_syrpvgrw.f"
		k += 2;
#line 253 "cla_syrpvgrw.f"
	    }
#line 254 "cla_syrpvgrw.f"
	}
#line 255 "cla_syrpvgrw.f"
    } else {
#line 256 "cla_syrpvgrw.f"
	k = 1;
#line 257 "cla_syrpvgrw.f"
	while(k <= ncols) {
#line 258 "cla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
/*              1x1 pivot */
#line 260 "cla_syrpvgrw.f"
		kp = ipiv[k];
#line 261 "cla_syrpvgrw.f"
		if (kp != k) {
#line 262 "cla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 263 "cla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 264 "cla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 265 "cla_syrpvgrw.f"
		}
#line 266 "cla_syrpvgrw.f"
		i__1 = *n;
#line 266 "cla_syrpvgrw.f"
		for (i__ = k; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 267 "cla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 267 "cla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 267 "cla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
#line 268 "cla_syrpvgrw.f"
		}
#line 269 "cla_syrpvgrw.f"
		++k;
#line 270 "cla_syrpvgrw.f"
	    } else {
/*              2x2 pivot */
#line 272 "cla_syrpvgrw.f"
		kp = -ipiv[k];
#line 273 "cla_syrpvgrw.f"
		tmp = work[*n + k + 1];
#line 274 "cla_syrpvgrw.f"
		work[*n + k + 1] = work[*n + kp];
#line 275 "cla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 276 "cla_syrpvgrw.f"
		i__1 = *n;
#line 276 "cla_syrpvgrw.f"
		for (i__ = k + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 277 "cla_syrpvgrw.f"
		    i__2 = i__ + k * af_dim1;
#line 277 "cla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + k * af_dim1]), abs(d__2)), d__4 = work[k]
			    ;
#line 277 "cla_syrpvgrw.f"
		    work[k] = max(d__3,d__4);
/* Computing MAX */
#line 278 "cla_syrpvgrw.f"
		    i__2 = i__ + (k + 1) * af_dim1;
#line 278 "cla_syrpvgrw.f"
		    d__3 = (d__1 = af[i__2].r, abs(d__1)) + (d__2 = d_imag(&
			    af[i__ + (k + 1) * af_dim1]), abs(d__2)), d__4 = 
			    work[k + 1];
#line 278 "cla_syrpvgrw.f"
		    work[k + 1] = max(d__3,d__4);
#line 280 "cla_syrpvgrw.f"
		}
/* Computing MAX */
#line 281 "cla_syrpvgrw.f"
		i__1 = k + k * af_dim1;
#line 281 "cla_syrpvgrw.f"
		d__3 = (d__1 = af[i__1].r, abs(d__1)) + (d__2 = d_imag(&af[k 
			+ k * af_dim1]), abs(d__2)), d__4 = work[k];
#line 281 "cla_syrpvgrw.f"
		work[k] = max(d__3,d__4);
#line 282 "cla_syrpvgrw.f"
		k += 2;
#line 283 "cla_syrpvgrw.f"
	    }
#line 284 "cla_syrpvgrw.f"
	}
#line 285 "cla_syrpvgrw.f"
	k = ncols;
#line 286 "cla_syrpvgrw.f"
	while(k >= 1) {
#line 287 "cla_syrpvgrw.f"
	    if (ipiv[k] > 0) {
#line 288 "cla_syrpvgrw.f"
		kp = ipiv[k];
#line 289 "cla_syrpvgrw.f"
		if (kp != k) {
#line 290 "cla_syrpvgrw.f"
		    tmp = work[*n + k];
#line 291 "cla_syrpvgrw.f"
		    work[*n + k] = work[*n + kp];
#line 292 "cla_syrpvgrw.f"
		    work[*n + kp] = tmp;
#line 293 "cla_syrpvgrw.f"
		}
#line 294 "cla_syrpvgrw.f"
		--k;
#line 295 "cla_syrpvgrw.f"
	    } else {
#line 296 "cla_syrpvgrw.f"
		kp = -ipiv[k];
#line 297 "cla_syrpvgrw.f"
		tmp = work[*n + k];
#line 298 "cla_syrpvgrw.f"
		work[*n + k] = work[*n + kp];
#line 299 "cla_syrpvgrw.f"
		work[*n + kp] = tmp;
#line 300 "cla_syrpvgrw.f"
		k += -2;
#line 301 "cla_syrpvgrw.f"
	    }
#line 302 "cla_syrpvgrw.f"
	}
#line 303 "cla_syrpvgrw.f"
    }

/*     Compute the *inverse* of the max element growth factor.  Dividing */
/*     by zero would imply the largest entry of the factor's column is */
/*     zero.  Than can happen when either the column of A is zero or */
/*     massive pivots made the factor underflow to zero.  Neither counts */
/*     as growth in itself, so simply ignore terms with zero */
/*     denominators. */

#line 312 "cla_syrpvgrw.f"
    if (upper) {
#line 313 "cla_syrpvgrw.f"
	i__1 = *n;
#line 313 "cla_syrpvgrw.f"
	for (i__ = ncols; i__ <= i__1; ++i__) {
#line 314 "cla_syrpvgrw.f"
	    umax = work[i__];
#line 315 "cla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 316 "cla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 317 "cla_syrpvgrw.f"
		d__1 = amax / umax;
#line 317 "cla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 318 "cla_syrpvgrw.f"
	    }
#line 319 "cla_syrpvgrw.f"
	}
#line 320 "cla_syrpvgrw.f"
    } else {
#line 321 "cla_syrpvgrw.f"
	i__1 = ncols;
#line 321 "cla_syrpvgrw.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 322 "cla_syrpvgrw.f"
	    umax = work[i__];
#line 323 "cla_syrpvgrw.f"
	    amax = work[*n + i__];
#line 324 "cla_syrpvgrw.f"
	    if (umax != 0.) {
/* Computing MIN */
#line 325 "cla_syrpvgrw.f"
		d__1 = amax / umax;
#line 325 "cla_syrpvgrw.f"
		rpvgrw = min(d__1,rpvgrw);
#line 326 "cla_syrpvgrw.f"
	    }
#line 327 "cla_syrpvgrw.f"
	}
#line 328 "cla_syrpvgrw.f"
    }
#line 330 "cla_syrpvgrw.f"
    ret_val = rpvgrw;
#line 331 "cla_syrpvgrw.f"
    return ret_val;
} /* cla_syrpvgrw__ */

