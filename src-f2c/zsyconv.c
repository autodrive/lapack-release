#line 1 "zsyconv.f"
/* zsyconv.f -- translated by f2c (version 20100827).
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

#line 1 "zsyconv.f"
/* > \brief \b ZSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYCONV converts A given by ZHETRF into L and D or vice-versa. */
/* > Get nondiagonal elements of D (returned in workspace) and */
/* > apply or reverse permutation done in TRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N) */
/* >          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1 */
/* >          or 2-by-2 block diagonal matrix D in LDLT. */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsyconv_(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *e, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ip;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


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
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 150 "zsyconv.f"
    /* Parameter adjustments */
#line 150 "zsyconv.f"
    a_dim1 = *lda;
#line 150 "zsyconv.f"
    a_offset = 1 + a_dim1;
#line 150 "zsyconv.f"
    a -= a_offset;
#line 150 "zsyconv.f"
    --ipiv;
#line 150 "zsyconv.f"
    --e;
#line 150 "zsyconv.f"

#line 150 "zsyconv.f"
    /* Function Body */
#line 150 "zsyconv.f"
    *info = 0;
#line 151 "zsyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 152 "zsyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 153 "zsyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 154 "zsyconv.f"
	*info = -1;
#line 155 "zsyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 156 "zsyconv.f"
	*info = -2;
#line 157 "zsyconv.f"
    } else if (*n < 0) {
#line 158 "zsyconv.f"
	*info = -3;
#line 159 "zsyconv.f"
    } else if (*lda < max(1,*n)) {
#line 160 "zsyconv.f"
	*info = -5;
#line 162 "zsyconv.f"
    }
#line 163 "zsyconv.f"
    if (*info != 0) {
#line 164 "zsyconv.f"
	i__1 = -(*info);
#line 164 "zsyconv.f"
	xerbla_("ZSYCONV", &i__1, (ftnlen)7);
#line 165 "zsyconv.f"
	return 0;
#line 166 "zsyconv.f"
    }

/*     Quick return if possible */

#line 170 "zsyconv.f"
    if (*n == 0) {
#line 170 "zsyconv.f"
	return 0;
#line 170 "zsyconv.f"
    }

#line 173 "zsyconv.f"
    if (upper) {

/*        A is UPPER */

#line 177 "zsyconv.f"
	if (convert) {

/*           Convert A (A is upper) */

/*           Convert VALUE */

#line 183 "zsyconv.f"
	    i__ = *n;
#line 184 "zsyconv.f"
	    e[1].r = 0., e[1].i = 0.;
#line 185 "zsyconv.f"
	    while(i__ > 1) {
#line 186 "zsyconv.f"
		if (ipiv[i__] < 0) {
#line 187 "zsyconv.f"
		    i__1 = i__;
#line 187 "zsyconv.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 187 "zsyconv.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 188 "zsyconv.f"
		    i__1 = i__ - 1;
#line 188 "zsyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 189 "zsyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 189 "zsyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 190 "zsyconv.f"
		    --i__;
#line 191 "zsyconv.f"
		} else {
#line 192 "zsyconv.f"
		    i__1 = i__;
#line 192 "zsyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 193 "zsyconv.f"
		}
#line 194 "zsyconv.f"
		--i__;
#line 195 "zsyconv.f"
	    }

/*           Convert PERMUTATIONS */

#line 199 "zsyconv.f"
	    i__ = *n;
#line 200 "zsyconv.f"
	    while(i__ >= 1) {
#line 201 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 202 "zsyconv.f"
		    ip = ipiv[i__];
#line 203 "zsyconv.f"
		    if (i__ < *n) {
#line 204 "zsyconv.f"
			i__1 = *n;
#line 204 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 205 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 205 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 206 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 206 "zsyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 206 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 207 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 207 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 208 "zsyconv.f"
/* L12: */
#line 208 "zsyconv.f"
			}
#line 209 "zsyconv.f"
		    }
#line 210 "zsyconv.f"
		} else {
#line 211 "zsyconv.f"
		    ip = -ipiv[i__];
#line 212 "zsyconv.f"
		    if (i__ < *n) {
#line 213 "zsyconv.f"
			i__1 = *n;
#line 213 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 214 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 214 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 215 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 215 "zsyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 215 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 216 "zsyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 216 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 217 "zsyconv.f"
/* L13: */
#line 217 "zsyconv.f"
			}
#line 218 "zsyconv.f"
		    }
#line 219 "zsyconv.f"
		    --i__;
#line 220 "zsyconv.f"
		}
#line 221 "zsyconv.f"
		--i__;
#line 222 "zsyconv.f"
	    }

#line 224 "zsyconv.f"
	} else {

/*           Revert A (A is upper) */

/*           Revert PERMUTATIONS */

#line 230 "zsyconv.f"
	    i__ = 1;
#line 231 "zsyconv.f"
	    while(i__ <= *n) {
#line 232 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 233 "zsyconv.f"
		    ip = ipiv[i__];
#line 234 "zsyconv.f"
		    if (i__ < *n) {
#line 235 "zsyconv.f"
			i__1 = *n;
#line 235 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 236 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 236 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 237 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 237 "zsyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 237 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 238 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 238 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 239 "zsyconv.f"
			}
#line 240 "zsyconv.f"
		    }
#line 241 "zsyconv.f"
		} else {
#line 242 "zsyconv.f"
		    ip = -ipiv[i__];
#line 243 "zsyconv.f"
		    ++i__;
#line 244 "zsyconv.f"
		    if (i__ < *n) {
#line 245 "zsyconv.f"
			i__1 = *n;
#line 245 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 246 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 246 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 247 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 247 "zsyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 247 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 248 "zsyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 248 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 249 "zsyconv.f"
			}
#line 250 "zsyconv.f"
		    }
#line 251 "zsyconv.f"
		}
#line 252 "zsyconv.f"
		++i__;
#line 253 "zsyconv.f"
	    }

/*           Revert VALUE */

#line 257 "zsyconv.f"
	    i__ = *n;
#line 258 "zsyconv.f"
	    while(i__ > 1) {
#line 259 "zsyconv.f"
		if (ipiv[i__] < 0) {
#line 260 "zsyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 260 "zsyconv.f"
		    i__2 = i__;
#line 260 "zsyconv.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 261 "zsyconv.f"
		    --i__;
#line 262 "zsyconv.f"
		}
#line 263 "zsyconv.f"
		--i__;
#line 264 "zsyconv.f"
	    }
#line 265 "zsyconv.f"
	}

#line 267 "zsyconv.f"
    } else {

/*        A is LOWER */

#line 271 "zsyconv.f"
	if (convert) {

/*           Convert A (A is lower) */

/*           Convert VALUE */

#line 277 "zsyconv.f"
	    i__ = 1;
#line 278 "zsyconv.f"
	    i__1 = *n;
#line 278 "zsyconv.f"
	    e[i__1].r = 0., e[i__1].i = 0.;
#line 279 "zsyconv.f"
	    while(i__ <= *n) {
#line 280 "zsyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 281 "zsyconv.f"
		    i__1 = i__;
#line 281 "zsyconv.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 281 "zsyconv.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 282 "zsyconv.f"
		    i__1 = i__ + 1;
#line 282 "zsyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 283 "zsyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 283 "zsyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 284 "zsyconv.f"
		    ++i__;
#line 285 "zsyconv.f"
		} else {
#line 286 "zsyconv.f"
		    i__1 = i__;
#line 286 "zsyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 287 "zsyconv.f"
		}
#line 288 "zsyconv.f"
		++i__;
#line 289 "zsyconv.f"
	    }

/*           Convert PERMUTATIONS */

#line 293 "zsyconv.f"
	    i__ = 1;
#line 294 "zsyconv.f"
	    while(i__ <= *n) {
#line 295 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 296 "zsyconv.f"
		    ip = ipiv[i__];
#line 297 "zsyconv.f"
		    if (i__ > 1) {
#line 298 "zsyconv.f"
			i__1 = i__ - 1;
#line 298 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 299 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 299 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 300 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 300 "zsyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 300 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 301 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 301 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 302 "zsyconv.f"
/* L22: */
#line 302 "zsyconv.f"
			}
#line 303 "zsyconv.f"
		    }
#line 304 "zsyconv.f"
		} else {
#line 305 "zsyconv.f"
		    ip = -ipiv[i__];
#line 306 "zsyconv.f"
		    if (i__ > 1) {
#line 307 "zsyconv.f"
			i__1 = i__ - 1;
#line 307 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 308 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 308 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 309 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 309 "zsyconv.f"
			    i__3 = i__ + 1 + j * a_dim1;
#line 309 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 310 "zsyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 310 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 311 "zsyconv.f"
/* L23: */
#line 311 "zsyconv.f"
			}
#line 312 "zsyconv.f"
		    }
#line 313 "zsyconv.f"
		    ++i__;
#line 314 "zsyconv.f"
		}
#line 315 "zsyconv.f"
		++i__;
#line 316 "zsyconv.f"
	    }

#line 318 "zsyconv.f"
	} else {

/*           Revert A (A is lower) */

/*           Revert PERMUTATIONS */

#line 324 "zsyconv.f"
	    i__ = *n;
#line 325 "zsyconv.f"
	    while(i__ >= 1) {
#line 326 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 327 "zsyconv.f"
		    ip = ipiv[i__];
#line 328 "zsyconv.f"
		    if (i__ > 1) {
#line 329 "zsyconv.f"
			i__1 = i__ - 1;
#line 329 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 330 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 330 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 331 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 331 "zsyconv.f"
			    i__3 = ip + j * a_dim1;
#line 331 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 332 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 332 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 333 "zsyconv.f"
			}
#line 334 "zsyconv.f"
		    }
#line 335 "zsyconv.f"
		} else {
#line 336 "zsyconv.f"
		    ip = -ipiv[i__];
#line 337 "zsyconv.f"
		    --i__;
#line 338 "zsyconv.f"
		    if (i__ > 1) {
#line 339 "zsyconv.f"
			i__1 = i__ - 1;
#line 339 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 340 "zsyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 340 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 341 "zsyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 341 "zsyconv.f"
			    i__3 = ip + j * a_dim1;
#line 341 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 342 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 342 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 343 "zsyconv.f"
			}
#line 344 "zsyconv.f"
		    }
#line 345 "zsyconv.f"
		}
#line 346 "zsyconv.f"
		--i__;
#line 347 "zsyconv.f"
	    }

/*           Revert VALUE */

#line 351 "zsyconv.f"
	    i__ = 1;
#line 352 "zsyconv.f"
	    while(i__ <= *n - 1) {
#line 353 "zsyconv.f"
		if (ipiv[i__] < 0) {
#line 354 "zsyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 354 "zsyconv.f"
		    i__2 = i__;
#line 354 "zsyconv.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 355 "zsyconv.f"
		    ++i__;
#line 356 "zsyconv.f"
		}
#line 357 "zsyconv.f"
		++i__;
#line 358 "zsyconv.f"
	    }
#line 359 "zsyconv.f"
	}
#line 360 "zsyconv.f"
    }

#line 362 "zsyconv.f"
    return 0;

/*     End of ZSYCONV */

} /* zsyconv_ */

