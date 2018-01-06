#line 1 "csyconv.f"
/* csyconv.f -- translated by f2c (version 20100827).
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

#line 1 "csyconv.f"
/* > \brief \b CSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYCONV convert A given by TRF into L and D and vice-versa. */
/* > Get Non-diag elements of D (returned in workspace) and */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSYTRF. */
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
/* >          as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
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

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csyconv_(char *uplo, char *way, integer *n, 
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

#line 148 "csyconv.f"
    /* Parameter adjustments */
#line 148 "csyconv.f"
    a_dim1 = *lda;
#line 148 "csyconv.f"
    a_offset = 1 + a_dim1;
#line 148 "csyconv.f"
    a -= a_offset;
#line 148 "csyconv.f"
    --ipiv;
#line 148 "csyconv.f"
    --work;
#line 148 "csyconv.f"

#line 148 "csyconv.f"
    /* Function Body */
#line 148 "csyconv.f"
    *info = 0;
#line 149 "csyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 150 "csyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 151 "csyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 152 "csyconv.f"
	*info = -1;
#line 153 "csyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 154 "csyconv.f"
	*info = -2;
#line 155 "csyconv.f"
    } else if (*n < 0) {
#line 156 "csyconv.f"
	*info = -3;
#line 157 "csyconv.f"
    } else if (*lda < max(1,*n)) {
#line 158 "csyconv.f"
	*info = -5;
#line 160 "csyconv.f"
    }
#line 161 "csyconv.f"
    if (*info != 0) {
#line 162 "csyconv.f"
	i__1 = -(*info);
#line 162 "csyconv.f"
	xerbla_("CSYCONV", &i__1, (ftnlen)7);
#line 163 "csyconv.f"
	return 0;
#line 164 "csyconv.f"
    }

/*     Quick return if possible */

#line 168 "csyconv.f"
    if (*n == 0) {
#line 168 "csyconv.f"
	return 0;
#line 168 "csyconv.f"
    }

#line 171 "csyconv.f"
    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

#line 179 "csyconv.f"
	if (convert) {
#line 180 "csyconv.f"
	    i__ = *n;
#line 181 "csyconv.f"
	    work[1].r = 0., work[1].i = 0.;
#line 182 "csyconv.f"
	    while(i__ > 1) {
#line 183 "csyconv.f"
		if (ipiv[i__] < 0) {
#line 184 "csyconv.f"
		    i__1 = i__;
#line 184 "csyconv.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 184 "csyconv.f"
		    work[i__1].r = a[i__2].r, work[i__1].i = a[i__2].i;
#line 185 "csyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 185 "csyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 186 "csyconv.f"
		    --i__;
#line 187 "csyconv.f"
		} else {
#line 188 "csyconv.f"
		    i__1 = i__;
#line 188 "csyconv.f"
		    work[i__1].r = 0., work[i__1].i = 0.;
#line 189 "csyconv.f"
		}
#line 190 "csyconv.f"
		--i__;
#line 191 "csyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 195 "csyconv.f"
	    i__ = *n;
#line 196 "csyconv.f"
	    while(i__ >= 1) {
#line 197 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 198 "csyconv.f"
		    ip = ipiv[i__];
#line 199 "csyconv.f"
		    if (i__ < *n) {
#line 200 "csyconv.f"
			i__1 = *n;
#line 200 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 201 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 201 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 202 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 202 "csyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 202 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 203 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 203 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 204 "csyconv.f"
/* L12: */
#line 204 "csyconv.f"
			}
#line 205 "csyconv.f"
		    }
#line 206 "csyconv.f"
		} else {
#line 207 "csyconv.f"
		    ip = -ipiv[i__];
#line 208 "csyconv.f"
		    if (i__ < *n) {
#line 209 "csyconv.f"
			i__1 = *n;
#line 209 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 210 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 210 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 211 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 211 "csyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 211 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 212 "csyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 212 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 213 "csyconv.f"
/* L13: */
#line 213 "csyconv.f"
			}
#line 214 "csyconv.f"
		    }
#line 215 "csyconv.f"
		    --i__;
#line 216 "csyconv.f"
		}
#line 217 "csyconv.f"
		--i__;
#line 218 "csyconv.f"
	    }
#line 220 "csyconv.f"
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

#line 227 "csyconv.f"
	    i__ = 1;
#line 228 "csyconv.f"
	    while(i__ <= *n) {
#line 229 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 230 "csyconv.f"
		    ip = ipiv[i__];
#line 231 "csyconv.f"
		    if (i__ < *n) {
#line 232 "csyconv.f"
			i__1 = *n;
#line 232 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 233 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 233 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 234 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 234 "csyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 234 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 235 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 235 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 236 "csyconv.f"
			}
#line 237 "csyconv.f"
		    }
#line 238 "csyconv.f"
		} else {
#line 239 "csyconv.f"
		    ip = -ipiv[i__];
#line 240 "csyconv.f"
		    ++i__;
#line 241 "csyconv.f"
		    if (i__ < *n) {
#line 242 "csyconv.f"
			i__1 = *n;
#line 242 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 243 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 243 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 244 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 244 "csyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 244 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 245 "csyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 245 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 246 "csyconv.f"
			}
#line 247 "csyconv.f"
		    }
#line 248 "csyconv.f"
		}
#line 249 "csyconv.f"
		++i__;
#line 250 "csyconv.f"
	    }

/*        Revert VALUE */

#line 254 "csyconv.f"
	    i__ = *n;
#line 255 "csyconv.f"
	    while(i__ > 1) {
#line 256 "csyconv.f"
		if (ipiv[i__] < 0) {
#line 257 "csyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 257 "csyconv.f"
		    i__2 = i__;
#line 257 "csyconv.f"
		    a[i__1].r = work[i__2].r, a[i__1].i = work[i__2].i;
#line 258 "csyconv.f"
		    --i__;
#line 259 "csyconv.f"
		}
#line 260 "csyconv.f"
		--i__;
#line 261 "csyconv.f"
	    }
#line 262 "csyconv.f"
	}
#line 263 "csyconv.f"
    } else {

/*      A is LOWER */

#line 267 "csyconv.f"
	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

#line 274 "csyconv.f"
	    i__ = 1;
#line 275 "csyconv.f"
	    i__1 = *n;
#line 275 "csyconv.f"
	    work[i__1].r = 0., work[i__1].i = 0.;
#line 276 "csyconv.f"
	    while(i__ <= *n) {
#line 277 "csyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 278 "csyconv.f"
		    i__1 = i__;
#line 278 "csyconv.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 278 "csyconv.f"
		    work[i__1].r = a[i__2].r, work[i__1].i = a[i__2].i;
#line 279 "csyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 279 "csyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 280 "csyconv.f"
		    ++i__;
#line 281 "csyconv.f"
		} else {
#line 282 "csyconv.f"
		    i__1 = i__;
#line 282 "csyconv.f"
		    work[i__1].r = 0., work[i__1].i = 0.;
#line 283 "csyconv.f"
		}
#line 284 "csyconv.f"
		++i__;
#line 285 "csyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 289 "csyconv.f"
	    i__ = 1;
#line 290 "csyconv.f"
	    while(i__ <= *n) {
#line 291 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 292 "csyconv.f"
		    ip = ipiv[i__];
#line 293 "csyconv.f"
		    if (i__ > 1) {
#line 294 "csyconv.f"
			i__1 = i__ - 1;
#line 294 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 295 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 295 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 296 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 296 "csyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 296 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 297 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 297 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 298 "csyconv.f"
/* L22: */
#line 298 "csyconv.f"
			}
#line 299 "csyconv.f"
		    }
#line 300 "csyconv.f"
		} else {
#line 301 "csyconv.f"
		    ip = -ipiv[i__];
#line 302 "csyconv.f"
		    if (i__ > 1) {
#line 303 "csyconv.f"
			i__1 = i__ - 1;
#line 303 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 304 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 304 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 305 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 305 "csyconv.f"
			    i__3 = i__ + 1 + j * a_dim1;
#line 305 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 306 "csyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 306 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 307 "csyconv.f"
/* L23: */
#line 307 "csyconv.f"
			}
#line 308 "csyconv.f"
		    }
#line 309 "csyconv.f"
		    ++i__;
#line 310 "csyconv.f"
		}
#line 311 "csyconv.f"
		++i__;
#line 312 "csyconv.f"
	    }
#line 313 "csyconv.f"
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

#line 320 "csyconv.f"
	    i__ = *n;
#line 321 "csyconv.f"
	    while(i__ >= 1) {
#line 322 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 323 "csyconv.f"
		    ip = ipiv[i__];
#line 324 "csyconv.f"
		    if (i__ > 1) {
#line 325 "csyconv.f"
			i__1 = i__ - 1;
#line 325 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 326 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 326 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 327 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 327 "csyconv.f"
			    i__3 = ip + j * a_dim1;
#line 327 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 328 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 328 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 329 "csyconv.f"
			}
#line 330 "csyconv.f"
		    }
#line 331 "csyconv.f"
		} else {
#line 332 "csyconv.f"
		    ip = -ipiv[i__];
#line 333 "csyconv.f"
		    --i__;
#line 334 "csyconv.f"
		    if (i__ > 1) {
#line 335 "csyconv.f"
			i__1 = i__ - 1;
#line 335 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 336 "csyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 336 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 337 "csyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 337 "csyconv.f"
			    i__3 = ip + j * a_dim1;
#line 337 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 338 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 338 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 339 "csyconv.f"
			}
#line 340 "csyconv.f"
		    }
#line 341 "csyconv.f"
		}
#line 342 "csyconv.f"
		--i__;
#line 343 "csyconv.f"
	    }

/*        Revert VALUE */

#line 347 "csyconv.f"
	    i__ = 1;
#line 348 "csyconv.f"
	    while(i__ <= *n - 1) {
#line 349 "csyconv.f"
		if (ipiv[i__] < 0) {
#line 350 "csyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 350 "csyconv.f"
		    i__2 = i__;
#line 350 "csyconv.f"
		    a[i__1].r = work[i__2].r, a[i__1].i = work[i__2].i;
#line 351 "csyconv.f"
		    ++i__;
#line 352 "csyconv.f"
		}
#line 353 "csyconv.f"
		++i__;
#line 354 "csyconv.f"
	    }
#line 355 "csyconv.f"
	}
#line 356 "csyconv.f"
    }
#line 358 "csyconv.f"
    return 0;

/*     End of CSYCONV */

} /* csyconv_ */

