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

/*       SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), E( * ) */
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
/* > \param[in,out] A */
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
/* > \param[out] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N) */
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

/* > \date November 2015 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csyconv_(char *uplo, char *way, integer *n, 
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


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 150 "csyconv.f"
    /* Parameter adjustments */
#line 150 "csyconv.f"
    a_dim1 = *lda;
#line 150 "csyconv.f"
    a_offset = 1 + a_dim1;
#line 150 "csyconv.f"
    a -= a_offset;
#line 150 "csyconv.f"
    --ipiv;
#line 150 "csyconv.f"
    --e;
#line 150 "csyconv.f"

#line 150 "csyconv.f"
    /* Function Body */
#line 150 "csyconv.f"
    *info = 0;
#line 151 "csyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 152 "csyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 153 "csyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 154 "csyconv.f"
	*info = -1;
#line 155 "csyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 156 "csyconv.f"
	*info = -2;
#line 157 "csyconv.f"
    } else if (*n < 0) {
#line 158 "csyconv.f"
	*info = -3;
#line 159 "csyconv.f"
    } else if (*lda < max(1,*n)) {
#line 160 "csyconv.f"
	*info = -5;
#line 162 "csyconv.f"
    }
#line 163 "csyconv.f"
    if (*info != 0) {
#line 164 "csyconv.f"
	i__1 = -(*info);
#line 164 "csyconv.f"
	xerbla_("CSYCONV", &i__1, (ftnlen)7);
#line 165 "csyconv.f"
	return 0;
#line 166 "csyconv.f"
    }

/*     Quick return if possible */

#line 170 "csyconv.f"
    if (*n == 0) {
#line 170 "csyconv.f"
	return 0;
#line 170 "csyconv.f"
    }

#line 173 "csyconv.f"
    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

#line 181 "csyconv.f"
	if (convert) {
#line 182 "csyconv.f"
	    i__ = *n;
#line 183 "csyconv.f"
	    e[1].r = 0., e[1].i = 0.;
#line 184 "csyconv.f"
	    while(i__ > 1) {
#line 185 "csyconv.f"
		if (ipiv[i__] < 0) {
#line 186 "csyconv.f"
		    i__1 = i__;
#line 186 "csyconv.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 186 "csyconv.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 187 "csyconv.f"
		    i__1 = i__ - 1;
#line 187 "csyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 188 "csyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 188 "csyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 189 "csyconv.f"
		    --i__;
#line 190 "csyconv.f"
		} else {
#line 191 "csyconv.f"
		    i__1 = i__;
#line 191 "csyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 192 "csyconv.f"
		}
#line 193 "csyconv.f"
		--i__;
#line 194 "csyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 198 "csyconv.f"
	    i__ = *n;
#line 199 "csyconv.f"
	    while(i__ >= 1) {
#line 200 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 201 "csyconv.f"
		    ip = ipiv[i__];
#line 202 "csyconv.f"
		    if (i__ < *n) {
#line 203 "csyconv.f"
			i__1 = *n;
#line 203 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 204 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 204 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 205 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 205 "csyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 205 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 206 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 206 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 207 "csyconv.f"
/* L12: */
#line 207 "csyconv.f"
			}
#line 208 "csyconv.f"
		    }
#line 209 "csyconv.f"
		} else {
#line 210 "csyconv.f"
		    ip = -ipiv[i__];
#line 211 "csyconv.f"
		    if (i__ < *n) {
#line 212 "csyconv.f"
			i__1 = *n;
#line 212 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 213 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 213 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 214 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 214 "csyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 214 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 215 "csyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 215 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 216 "csyconv.f"
/* L13: */
#line 216 "csyconv.f"
			}
#line 217 "csyconv.f"
		    }
#line 218 "csyconv.f"
		    --i__;
#line 219 "csyconv.f"
		}
#line 220 "csyconv.f"
		--i__;
#line 221 "csyconv.f"
	    }
#line 223 "csyconv.f"
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

#line 230 "csyconv.f"
	    i__ = 1;
#line 231 "csyconv.f"
	    while(i__ <= *n) {
#line 232 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 233 "csyconv.f"
		    ip = ipiv[i__];
#line 234 "csyconv.f"
		    if (i__ < *n) {
#line 235 "csyconv.f"
			i__1 = *n;
#line 235 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 236 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 236 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 237 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 237 "csyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 237 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 238 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 238 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 239 "csyconv.f"
			}
#line 240 "csyconv.f"
		    }
#line 241 "csyconv.f"
		} else {
#line 242 "csyconv.f"
		    ip = -ipiv[i__];
#line 243 "csyconv.f"
		    ++i__;
#line 244 "csyconv.f"
		    if (i__ < *n) {
#line 245 "csyconv.f"
			i__1 = *n;
#line 245 "csyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 246 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 246 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 247 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 247 "csyconv.f"
			    i__3 = i__ - 1 + j * a_dim1;
#line 247 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 248 "csyconv.f"
			    i__2 = i__ - 1 + j * a_dim1;
#line 248 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 249 "csyconv.f"
			}
#line 250 "csyconv.f"
		    }
#line 251 "csyconv.f"
		}
#line 252 "csyconv.f"
		++i__;
#line 253 "csyconv.f"
	    }

/*        Revert VALUE */

#line 257 "csyconv.f"
	    i__ = *n;
#line 258 "csyconv.f"
	    while(i__ > 1) {
#line 259 "csyconv.f"
		if (ipiv[i__] < 0) {
#line 260 "csyconv.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 260 "csyconv.f"
		    i__2 = i__;
#line 260 "csyconv.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 261 "csyconv.f"
		    --i__;
#line 262 "csyconv.f"
		}
#line 263 "csyconv.f"
		--i__;
#line 264 "csyconv.f"
	    }
#line 265 "csyconv.f"
	}
#line 266 "csyconv.f"
    } else {

/*      A is LOWER */

#line 270 "csyconv.f"
	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

#line 277 "csyconv.f"
	    i__ = 1;
#line 278 "csyconv.f"
	    i__1 = *n;
#line 278 "csyconv.f"
	    e[i__1].r = 0., e[i__1].i = 0.;
#line 279 "csyconv.f"
	    while(i__ <= *n) {
#line 280 "csyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 281 "csyconv.f"
		    i__1 = i__;
#line 281 "csyconv.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 281 "csyconv.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 282 "csyconv.f"
		    i__1 = i__ + 1;
#line 282 "csyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 283 "csyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 283 "csyconv.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 284 "csyconv.f"
		    ++i__;
#line 285 "csyconv.f"
		} else {
#line 286 "csyconv.f"
		    i__1 = i__;
#line 286 "csyconv.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 287 "csyconv.f"
		}
#line 288 "csyconv.f"
		++i__;
#line 289 "csyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 293 "csyconv.f"
	    i__ = 1;
#line 294 "csyconv.f"
	    while(i__ <= *n) {
#line 295 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 296 "csyconv.f"
		    ip = ipiv[i__];
#line 297 "csyconv.f"
		    if (i__ > 1) {
#line 298 "csyconv.f"
			i__1 = i__ - 1;
#line 298 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 299 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 299 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 300 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 300 "csyconv.f"
			    i__3 = i__ + j * a_dim1;
#line 300 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 301 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 301 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 302 "csyconv.f"
/* L22: */
#line 302 "csyconv.f"
			}
#line 303 "csyconv.f"
		    }
#line 304 "csyconv.f"
		} else {
#line 305 "csyconv.f"
		    ip = -ipiv[i__];
#line 306 "csyconv.f"
		    if (i__ > 1) {
#line 307 "csyconv.f"
			i__1 = i__ - 1;
#line 307 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 308 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 308 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 309 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 309 "csyconv.f"
			    i__3 = i__ + 1 + j * a_dim1;
#line 309 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 310 "csyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 310 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 311 "csyconv.f"
/* L23: */
#line 311 "csyconv.f"
			}
#line 312 "csyconv.f"
		    }
#line 313 "csyconv.f"
		    ++i__;
#line 314 "csyconv.f"
		}
#line 315 "csyconv.f"
		++i__;
#line 316 "csyconv.f"
	    }
#line 317 "csyconv.f"
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

#line 324 "csyconv.f"
	    i__ = *n;
#line 325 "csyconv.f"
	    while(i__ >= 1) {
#line 326 "csyconv.f"
		if (ipiv[i__] > 0) {
#line 327 "csyconv.f"
		    ip = ipiv[i__];
#line 328 "csyconv.f"
		    if (i__ > 1) {
#line 329 "csyconv.f"
			i__1 = i__ - 1;
#line 329 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 330 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 330 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 331 "csyconv.f"
			    i__2 = i__ + j * a_dim1;
#line 331 "csyconv.f"
			    i__3 = ip + j * a_dim1;
#line 331 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 332 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 332 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 333 "csyconv.f"
			}
#line 334 "csyconv.f"
		    }
#line 335 "csyconv.f"
		} else {
#line 336 "csyconv.f"
		    ip = -ipiv[i__];
#line 337 "csyconv.f"
		    --i__;
#line 338 "csyconv.f"
		    if (i__ > 1) {
#line 339 "csyconv.f"
			i__1 = i__ - 1;
#line 339 "csyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 340 "csyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 340 "csyconv.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 341 "csyconv.f"
			    i__2 = i__ + 1 + j * a_dim1;
#line 341 "csyconv.f"
			    i__3 = ip + j * a_dim1;
#line 341 "csyconv.f"
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 342 "csyconv.f"
			    i__2 = ip + j * a_dim1;
#line 342 "csyconv.f"
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 343 "csyconv.f"
			}
#line 344 "csyconv.f"
		    }
#line 345 "csyconv.f"
		}
#line 346 "csyconv.f"
		--i__;
#line 347 "csyconv.f"
	    }

/*        Revert VALUE */

#line 351 "csyconv.f"
	    i__ = 1;
#line 352 "csyconv.f"
	    while(i__ <= *n - 1) {
#line 353 "csyconv.f"
		if (ipiv[i__] < 0) {
#line 354 "csyconv.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 354 "csyconv.f"
		    i__2 = i__;
#line 354 "csyconv.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 355 "csyconv.f"
		    ++i__;
#line 356 "csyconv.f"
		}
#line 357 "csyconv.f"
		++i__;
#line 358 "csyconv.f"
	    }
#line 359 "csyconv.f"
	}
#line 360 "csyconv.f"
    }
#line 362 "csyconv.f"
    return 0;

/*     End of CSYCONV */

} /* csyconv_ */

