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

/*       SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
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
/* > \param[in] A */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
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

/* > \date November 2011 */

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsyconv_(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work, 
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


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 148 "zsyconv.f"
    /* Parameter adjustments */
#line 148 "zsyconv.f"
    a_dim1 = *lda;
#line 148 "zsyconv.f"
    a_offset = 1 + a_dim1;
#line 148 "zsyconv.f"
    a -= a_offset;
#line 148 "zsyconv.f"
    --ipiv;
#line 148 "zsyconv.f"
    --work;
#line 148 "zsyconv.f"

#line 148 "zsyconv.f"
    /* Function Body */
#line 148 "zsyconv.f"
    *info = 0;
#line 149 "zsyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 150 "zsyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 151 "zsyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 152 "zsyconv.f"
	*info = -1;
#line 153 "zsyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 154 "zsyconv.f"
	*info = -2;
#line 155 "zsyconv.f"
    } else if (*n < 0) {
#line 156 "zsyconv.f"
	*info = -3;
#line 157 "zsyconv.f"
    } else if (*lda < max(1,*n)) {
#line 158 "zsyconv.f"
	*info = -5;
#line 160 "zsyconv.f"
    }
#line 161 "zsyconv.f"
    if (*info != 0) {
#line 162 "zsyconv.f"
	i__1 = -(*info);
#line 162 "zsyconv.f"
	xerbla_("ZSYCONV", &i__1, (ftnlen)7);
#line 163 "zsyconv.f"
	return 0;
#line 164 "zsyconv.f"
    }

/*     Quick return if possible */

#line 168 "zsyconv.f"
    if (*n == 0) {
#line 168 "zsyconv.f"
	return 0;
#line 168 "zsyconv.f"
    }

#line 171 "zsyconv.f"
    if (upper) {

/*        A is UPPER */

#line 175 "zsyconv.f"
	if (convert) {

/*           Convert A (A is upper) */

/*           Convert VALUE */

#line 181 "zsyconv.f"
	    i__ = *n;
#line 182 "zsyconv.f"
	    work[1].r = 0., work[1].i = 0.;
#line 183 "zsyconv.f"
	    while(i__ > 1) {
#line 184 "zsyconv.f"
		if (ipiv[i__] < 0) {
#line 185 "zsyconv.f"
		    i__1 = i__;
#line 185 "zsyconv.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 185 "zsyconv.f"
		    work[i__1].r = a[i__2].r, work[i__1].i = a[i__2].i;
#line 186 "zsyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 186 "zsyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 187 "zsyconv.f"
		    --i__;
#line 188 "zsyconv.f"
		} else {
#line 189 "zsyconv.f"
		    i__1 = i__;
#line 189 "zsyconv.f"
		    work[i__1].r = 0., work[i__1].i = 0.;
#line 190 "zsyconv.f"
		}
#line 191 "zsyconv.f"
		--i__;
#line 192 "zsyconv.f"
	    }

/*           Convert PERMUTATIONS */

#line 196 "zsyconv.f"
	    i__ = *n;
#line 197 "zsyconv.f"
	    while(i__ >= 1) {
#line 198 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 199 "zsyconv.f"
		    ip = ipiv[i__];
#line 200 "zsyconv.f"
		    if (i__ < *n) {
#line 201 "zsyconv.f"
			i__1 = *n;
#line 201 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 202 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 202 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 203 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 203 "zsyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 203 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 204 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 204 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 205 "zsyconv.f"
/* L12: */
#line 205 "zsyconv.f"
			}
#line 206 "zsyconv.f"
		    }
#line 207 "zsyconv.f"
		} else {
#line 208 "zsyconv.f"
		    ip = -ipiv[i__];
#line 209 "zsyconv.f"
		    if (i__ < *n) {
#line 210 "zsyconv.f"
			i__1 = *n;
#line 210 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 211 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 211 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 212 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 212 "zsyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 212 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 213 "zsyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 213 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 214 "zsyconv.f"
/* L13: */
#line 214 "zsyconv.f"
			}
#line 215 "zsyconv.f"
		    }
#line 216 "zsyconv.f"
		    --i__;
#line 217 "zsyconv.f"
		}
#line 218 "zsyconv.f"
		--i__;
#line 219 "zsyconv.f"
	    }

#line 221 "zsyconv.f"
	} else {

/*           Revert A (A is upper) */

/*           Revert PERMUTATIONS */

#line 227 "zsyconv.f"
	    i__ = 1;
#line 228 "zsyconv.f"
	    while(i__ <= *n) {
#line 229 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 230 "zsyconv.f"
		    ip = ipiv[i__];
#line 231 "zsyconv.f"
		    if (i__ < *n) {
#line 232 "zsyconv.f"
			i__1 = *n;
#line 232 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 233 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 233 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 234 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 234 "zsyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 234 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 235 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 235 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 236 "zsyconv.f"
			}
#line 237 "zsyconv.f"
		    }
#line 238 "zsyconv.f"
		} else {
#line 239 "zsyconv.f"
		    ip = -ipiv[i__];
#line 240 "zsyconv.f"
		    ++i__;
#line 241 "zsyconv.f"
		    if (i__ < *n) {
#line 242 "zsyconv.f"
			i__1 = *n;
#line 242 "zsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 243 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 243 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 244 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 244 "zsyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 244 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 245 "zsyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 245 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 246 "zsyconv.f"
			}
#line 247 "zsyconv.f"
		    }
#line 248 "zsyconv.f"
		}
#line 249 "zsyconv.f"
		++i__;
#line 250 "zsyconv.f"
	    }

/*           Revert VALUE */

#line 254 "zsyconv.f"
	    i__ = *n;
#line 255 "zsyconv.f"
	    while(i__ > 1) {
#line 256 "zsyconv.f"
		if (ipiv[i__] < 0) {
#line 257 "zsyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 257 "zsyconv.f"
		    i__2 = i__;
#line 257 "zsyconv.f"
		    a[i__1].r = work[i__2].r, a[i__1].i = work[i__2].i;
#line 258 "zsyconv.f"
		    --i__;
#line 259 "zsyconv.f"
		}
#line 260 "zsyconv.f"
		--i__;
#line 261 "zsyconv.f"
	    }
#line 262 "zsyconv.f"
	}

#line 264 "zsyconv.f"
    } else {

/*        A is LOWER */

#line 268 "zsyconv.f"
	if (convert) {

/*           Convert A (A is lower) */

/*           Convert VALUE */

#line 274 "zsyconv.f"
	    i__ = 1;
#line 275 "zsyconv.f"
	    i__1 = *n;
#line 275 "zsyconv.f"
	    work[i__1].r = 0., work[i__1].i = 0.;
#line 276 "zsyconv.f"
	    while(i__ <= *n) {
#line 277 "zsyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 278 "zsyconv.f"
		    i__1 = i__;
#line 278 "zsyconv.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 278 "zsyconv.f"
		    work[i__1].r = a[i__2].r, work[i__1].i = a[i__2].i;
#line 279 "zsyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 279 "zsyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 280 "zsyconv.f"
		    ++i__;
#line 281 "zsyconv.f"
		} else {
#line 282 "zsyconv.f"
		    i__1 = i__;
#line 282 "zsyconv.f"
		    work[i__1].r = 0., work[i__1].i = 0.;
#line 283 "zsyconv.f"
		}
#line 284 "zsyconv.f"
		++i__;
#line 285 "zsyconv.f"
	    }

/*           Convert PERMUTATIONS */

#line 289 "zsyconv.f"
	    i__ = 1;
#line 290 "zsyconv.f"
	    while(i__ <= *n) {
#line 291 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 292 "zsyconv.f"
		    ip = ipiv[i__];
#line 293 "zsyconv.f"
		    if (i__ > 1) {
#line 294 "zsyconv.f"
			i__1 = i__ - 1;
#line 294 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 295 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 295 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 296 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 296 "zsyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 296 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 297 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 297 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 298 "zsyconv.f"
/* L22: */
#line 298 "zsyconv.f"
			}
#line 299 "zsyconv.f"
		    }
#line 300 "zsyconv.f"
		} else {
#line 301 "zsyconv.f"
		    ip = -ipiv[i__];
#line 302 "zsyconv.f"
		    if (i__ > 1) {
#line 303 "zsyconv.f"
			i__1 = i__ - 1;
#line 303 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 304 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 304 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 305 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 305 "zsyconv.f"
			    i__3 = i__ + 1 + j * a_dim1;
#line 305 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 306 "zsyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 306 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 307 "zsyconv.f"
/* L23: */
#line 307 "zsyconv.f"
			}
#line 308 "zsyconv.f"
		    }
#line 309 "zsyconv.f"
		    ++i__;
#line 310 "zsyconv.f"
		}
#line 311 "zsyconv.f"
		++i__;
#line 312 "zsyconv.f"
	    }

#line 314 "zsyconv.f"
	} else {

/*           Revert A (A is lower) */

/*           Revert PERMUTATIONS */

#line 320 "zsyconv.f"
	    i__ = *n;
#line 321 "zsyconv.f"
	    while(i__ >= 1) {
#line 322 "zsyconv.f"
		if (ipiv[i__] > 0) {
#line 323 "zsyconv.f"
		    ip = ipiv[i__];
#line 324 "zsyconv.f"
		    if (i__ > 1) {
#line 325 "zsyconv.f"
			i__1 = i__ - 1;
#line 325 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 326 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 326 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 327 "zsyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 327 "zsyconv.f"
			    i__3 = ip + j * a_dim1;
#line 327 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 328 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 328 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 329 "zsyconv.f"
			}
#line 330 "zsyconv.f"
		    }
#line 331 "zsyconv.f"
		} else {
#line 332 "zsyconv.f"
		    ip = -ipiv[i__];
#line 333 "zsyconv.f"
		    --i__;
#line 334 "zsyconv.f"
		    if (i__ > 1) {
#line 335 "zsyconv.f"
			i__1 = i__ - 1;
#line 335 "zsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 336 "zsyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 336 "zsyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 337 "zsyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 337 "zsyconv.f"
			    i__3 = ip + j * a_dim1;
#line 337 "zsyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 338 "zsyconv.f"
			    i__2 = ip + j * a_dim1;
#line 338 "zsyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 339 "zsyconv.f"
			}
#line 340 "zsyconv.f"
		    }
#line 341 "zsyconv.f"
		}
#line 342 "zsyconv.f"
		--i__;
#line 343 "zsyconv.f"
	    }

/*           Revert VALUE */

#line 347 "zsyconv.f"
	    i__ = 1;
#line 348 "zsyconv.f"
	    while(i__ <= *n - 1) {
#line 349 "zsyconv.f"
		if (ipiv[i__] < 0) {
#line 350 "zsyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 350 "zsyconv.f"
		    i__2 = i__;
#line 350 "zsyconv.f"
		    a[i__1].r = work[i__2].r, a[i__1].i = work[i__2].i;
#line 351 "zsyconv.f"
		    ++i__;
#line 352 "zsyconv.f"
		}
#line 353 "zsyconv.f"
		++i__;
#line 354 "zsyconv.f"
	    }
#line 355 "zsyconv.f"
	}
#line 356 "zsyconv.f"
    }

#line 358 "zsyconv.f"
    return 0;

/*     End of ZSYCONV */

} /* zsyconv_ */

