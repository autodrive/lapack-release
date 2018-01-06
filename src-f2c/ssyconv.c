#line 1 "ssyconv.f"
/* ssyconv.f -- translated by f2c (version 20100827).
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

#line 1 "ssyconv.f"
/* > \brief \b SSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYCONV convert A given by TRF into L and D and vice-versa. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by SSYTRF. */
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
/* >          as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssyconv_(char *uplo, char *way, integer *n, doublereal *
	a, integer *lda, integer *ipiv, doublereal *e, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, j, ip;
    static doublereal temp;
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

#line 150 "ssyconv.f"
    /* Parameter adjustments */
#line 150 "ssyconv.f"
    a_dim1 = *lda;
#line 150 "ssyconv.f"
    a_offset = 1 + a_dim1;
#line 150 "ssyconv.f"
    a -= a_offset;
#line 150 "ssyconv.f"
    --ipiv;
#line 150 "ssyconv.f"
    --e;
#line 150 "ssyconv.f"

#line 150 "ssyconv.f"
    /* Function Body */
#line 150 "ssyconv.f"
    *info = 0;
#line 151 "ssyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 152 "ssyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 153 "ssyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 154 "ssyconv.f"
	*info = -1;
#line 155 "ssyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 156 "ssyconv.f"
	*info = -2;
#line 157 "ssyconv.f"
    } else if (*n < 0) {
#line 158 "ssyconv.f"
	*info = -3;
#line 159 "ssyconv.f"
    } else if (*lda < max(1,*n)) {
#line 160 "ssyconv.f"
	*info = -5;
#line 162 "ssyconv.f"
    }
#line 163 "ssyconv.f"
    if (*info != 0) {
#line 164 "ssyconv.f"
	i__1 = -(*info);
#line 164 "ssyconv.f"
	xerbla_("SSYCONV", &i__1, (ftnlen)7);
#line 165 "ssyconv.f"
	return 0;
#line 166 "ssyconv.f"
    }

/*     Quick return if possible */

#line 170 "ssyconv.f"
    if (*n == 0) {
#line 170 "ssyconv.f"
	return 0;
#line 170 "ssyconv.f"
    }

#line 173 "ssyconv.f"
    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

#line 181 "ssyconv.f"
	if (convert) {
#line 182 "ssyconv.f"
	    i__ = *n;
#line 183 "ssyconv.f"
	    e[1] = 0.;
#line 184 "ssyconv.f"
	    while(i__ > 1) {
#line 185 "ssyconv.f"
		if (ipiv[i__] < 0) {
#line 186 "ssyconv.f"
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
#line 187 "ssyconv.f"
		    e[i__ - 1] = 0.;
#line 188 "ssyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = 0.;
#line 189 "ssyconv.f"
		    --i__;
#line 190 "ssyconv.f"
		} else {
#line 191 "ssyconv.f"
		    e[i__] = 0.;
#line 192 "ssyconv.f"
		}
#line 193 "ssyconv.f"
		--i__;
#line 194 "ssyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 198 "ssyconv.f"
	    i__ = *n;
#line 199 "ssyconv.f"
	    while(i__ >= 1) {
#line 200 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 201 "ssyconv.f"
		    ip = ipiv[i__];
#line 202 "ssyconv.f"
		    if (i__ < *n) {
#line 203 "ssyconv.f"
			i__1 = *n;
#line 203 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 204 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 205 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 206 "ssyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 207 "ssyconv.f"
/* L12: */
#line 207 "ssyconv.f"
			}
#line 208 "ssyconv.f"
		    }
#line 209 "ssyconv.f"
		} else {
#line 210 "ssyconv.f"
		    ip = -ipiv[i__];
#line 211 "ssyconv.f"
		    if (i__ < *n) {
#line 212 "ssyconv.f"
			i__1 = *n;
#line 212 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 213 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 214 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 215 "ssyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 216 "ssyconv.f"
/* L13: */
#line 216 "ssyconv.f"
			}
#line 217 "ssyconv.f"
		    }
#line 218 "ssyconv.f"
		    --i__;
#line 219 "ssyconv.f"
		}
#line 220 "ssyconv.f"
		--i__;
#line 221 "ssyconv.f"
	    }
#line 223 "ssyconv.f"
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

#line 230 "ssyconv.f"
	    i__ = 1;
#line 231 "ssyconv.f"
	    while(i__ <= *n) {
#line 232 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 233 "ssyconv.f"
		    ip = ipiv[i__];
#line 234 "ssyconv.f"
		    if (i__ < *n) {
#line 235 "ssyconv.f"
			i__1 = *n;
#line 235 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 236 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 237 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 238 "ssyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 239 "ssyconv.f"
			}
#line 240 "ssyconv.f"
		    }
#line 241 "ssyconv.f"
		} else {
#line 242 "ssyconv.f"
		    ip = -ipiv[i__];
#line 243 "ssyconv.f"
		    ++i__;
#line 244 "ssyconv.f"
		    if (i__ < *n) {
#line 245 "ssyconv.f"
			i__1 = *n;
#line 245 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 246 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 247 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 248 "ssyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 249 "ssyconv.f"
			}
#line 250 "ssyconv.f"
		    }
#line 251 "ssyconv.f"
		}
#line 252 "ssyconv.f"
		++i__;
#line 253 "ssyconv.f"
	    }

/*        Revert VALUE */

#line 257 "ssyconv.f"
	    i__ = *n;
#line 258 "ssyconv.f"
	    while(i__ > 1) {
#line 259 "ssyconv.f"
		if (ipiv[i__] < 0) {
#line 260 "ssyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
#line 261 "ssyconv.f"
		    --i__;
#line 262 "ssyconv.f"
		}
#line 263 "ssyconv.f"
		--i__;
#line 264 "ssyconv.f"
	    }
#line 265 "ssyconv.f"
	}
#line 266 "ssyconv.f"
    } else {

/*      A is LOWER */

#line 270 "ssyconv.f"
	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

#line 277 "ssyconv.f"
	    i__ = 1;
#line 278 "ssyconv.f"
	    e[*n] = 0.;
#line 279 "ssyconv.f"
	    while(i__ <= *n) {
#line 280 "ssyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 281 "ssyconv.f"
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 282 "ssyconv.f"
		    e[i__ + 1] = 0.;
#line 283 "ssyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = 0.;
#line 284 "ssyconv.f"
		    ++i__;
#line 285 "ssyconv.f"
		} else {
#line 286 "ssyconv.f"
		    e[i__] = 0.;
#line 287 "ssyconv.f"
		}
#line 288 "ssyconv.f"
		++i__;
#line 289 "ssyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 293 "ssyconv.f"
	    i__ = 1;
#line 294 "ssyconv.f"
	    while(i__ <= *n) {
#line 295 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 296 "ssyconv.f"
		    ip = ipiv[i__];
#line 297 "ssyconv.f"
		    if (i__ > 1) {
#line 298 "ssyconv.f"
			i__1 = i__ - 1;
#line 298 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 299 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 300 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 301 "ssyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 302 "ssyconv.f"
/* L22: */
#line 302 "ssyconv.f"
			}
#line 303 "ssyconv.f"
		    }
#line 304 "ssyconv.f"
		} else {
#line 305 "ssyconv.f"
		    ip = -ipiv[i__];
#line 306 "ssyconv.f"
		    if (i__ > 1) {
#line 307 "ssyconv.f"
			i__1 = i__ - 1;
#line 307 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 308 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 309 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
#line 310 "ssyconv.f"
			    a[i__ + 1 + j * a_dim1] = temp;
#line 311 "ssyconv.f"
/* L23: */
#line 311 "ssyconv.f"
			}
#line 312 "ssyconv.f"
		    }
#line 313 "ssyconv.f"
		    ++i__;
#line 314 "ssyconv.f"
		}
#line 315 "ssyconv.f"
		++i__;
#line 316 "ssyconv.f"
	    }
#line 317 "ssyconv.f"
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

#line 324 "ssyconv.f"
	    i__ = *n;
#line 325 "ssyconv.f"
	    while(i__ >= 1) {
#line 326 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 327 "ssyconv.f"
		    ip = ipiv[i__];
#line 328 "ssyconv.f"
		    if (i__ > 1) {
#line 329 "ssyconv.f"
			i__1 = i__ - 1;
#line 329 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 330 "ssyconv.f"
			    temp = a[i__ + j * a_dim1];
#line 331 "ssyconv.f"
			    a[i__ + j * a_dim1] = a[ip + j * a_dim1];
#line 332 "ssyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 333 "ssyconv.f"
			}
#line 334 "ssyconv.f"
		    }
#line 335 "ssyconv.f"
		} else {
#line 336 "ssyconv.f"
		    ip = -ipiv[i__];
#line 337 "ssyconv.f"
		    --i__;
#line 338 "ssyconv.f"
		    if (i__ > 1) {
#line 339 "ssyconv.f"
			i__1 = i__ - 1;
#line 339 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 340 "ssyconv.f"
			    temp = a[i__ + 1 + j * a_dim1];
#line 341 "ssyconv.f"
			    a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
#line 342 "ssyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 343 "ssyconv.f"
			}
#line 344 "ssyconv.f"
		    }
#line 345 "ssyconv.f"
		}
#line 346 "ssyconv.f"
		--i__;
#line 347 "ssyconv.f"
	    }

/*        Revert VALUE */

#line 351 "ssyconv.f"
	    i__ = 1;
#line 352 "ssyconv.f"
	    while(i__ <= *n - 1) {
#line 353 "ssyconv.f"
		if (ipiv[i__] < 0) {
#line 354 "ssyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 355 "ssyconv.f"
		    ++i__;
#line 356 "ssyconv.f"
		}
#line 357 "ssyconv.f"
		++i__;
#line 358 "ssyconv.f"
	    }
#line 359 "ssyconv.f"
	}
#line 360 "ssyconv.f"
    }
#line 362 "ssyconv.f"
    return 0;

/*     End of SSYCONV */

} /* ssyconv_ */

