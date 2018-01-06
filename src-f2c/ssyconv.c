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

/*       SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), WORK( * ) */
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
/* > \param[in] A */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssyconv_(char *uplo, char *way, integer *n, doublereal *
	a, integer *lda, integer *ipiv, doublereal *work, integer *info, 
	ftnlen uplo_len, ftnlen way_len)
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

#line 148 "ssyconv.f"
    /* Parameter adjustments */
#line 148 "ssyconv.f"
    a_dim1 = *lda;
#line 148 "ssyconv.f"
    a_offset = 1 + a_dim1;
#line 148 "ssyconv.f"
    a -= a_offset;
#line 148 "ssyconv.f"
    --ipiv;
#line 148 "ssyconv.f"
    --work;
#line 148 "ssyconv.f"

#line 148 "ssyconv.f"
    /* Function Body */
#line 148 "ssyconv.f"
    *info = 0;
#line 149 "ssyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 150 "ssyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 151 "ssyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 152 "ssyconv.f"
	*info = -1;
#line 153 "ssyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 154 "ssyconv.f"
	*info = -2;
#line 155 "ssyconv.f"
    } else if (*n < 0) {
#line 156 "ssyconv.f"
	*info = -3;
#line 157 "ssyconv.f"
    } else if (*lda < max(1,*n)) {
#line 158 "ssyconv.f"
	*info = -5;
#line 160 "ssyconv.f"
    }
#line 161 "ssyconv.f"
    if (*info != 0) {
#line 162 "ssyconv.f"
	i__1 = -(*info);
#line 162 "ssyconv.f"
	xerbla_("SSYCONV", &i__1, (ftnlen)7);
#line 163 "ssyconv.f"
	return 0;
#line 164 "ssyconv.f"
    }

/*     Quick return if possible */

#line 168 "ssyconv.f"
    if (*n == 0) {
#line 168 "ssyconv.f"
	return 0;
#line 168 "ssyconv.f"
    }

#line 171 "ssyconv.f"
    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

#line 179 "ssyconv.f"
	if (convert) {
#line 180 "ssyconv.f"
	    i__ = *n;
#line 181 "ssyconv.f"
	    work[1] = 0.;
#line 182 "ssyconv.f"
	    while(i__ > 1) {
#line 183 "ssyconv.f"
		if (ipiv[i__] < 0) {
#line 184 "ssyconv.f"
		    work[i__] = a[i__ - 1 + i__ * a_dim1];
#line 185 "ssyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = 0.;
#line 186 "ssyconv.f"
		    --i__;
#line 187 "ssyconv.f"
		} else {
#line 188 "ssyconv.f"
		    work[i__] = 0.;
#line 189 "ssyconv.f"
		}
#line 190 "ssyconv.f"
		--i__;
#line 191 "ssyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 195 "ssyconv.f"
	    i__ = *n;
#line 196 "ssyconv.f"
	    while(i__ >= 1) {
#line 197 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 198 "ssyconv.f"
		    ip = ipiv[i__];
#line 199 "ssyconv.f"
		    if (i__ < *n) {
#line 200 "ssyconv.f"
			i__1 = *n;
#line 200 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 201 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 202 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 203 "ssyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 204 "ssyconv.f"
/* L12: */
#line 204 "ssyconv.f"
			}
#line 205 "ssyconv.f"
		    }
#line 206 "ssyconv.f"
		} else {
#line 207 "ssyconv.f"
		    ip = -ipiv[i__];
#line 208 "ssyconv.f"
		    if (i__ < *n) {
#line 209 "ssyconv.f"
			i__1 = *n;
#line 209 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 210 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 211 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 212 "ssyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 213 "ssyconv.f"
/* L13: */
#line 213 "ssyconv.f"
			}
#line 214 "ssyconv.f"
		    }
#line 215 "ssyconv.f"
		    --i__;
#line 216 "ssyconv.f"
		}
#line 217 "ssyconv.f"
		--i__;
#line 218 "ssyconv.f"
	    }
#line 220 "ssyconv.f"
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

#line 227 "ssyconv.f"
	    i__ = 1;
#line 228 "ssyconv.f"
	    while(i__ <= *n) {
#line 229 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 230 "ssyconv.f"
		    ip = ipiv[i__];
#line 231 "ssyconv.f"
		    if (i__ < *n) {
#line 232 "ssyconv.f"
			i__1 = *n;
#line 232 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 233 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 234 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 235 "ssyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 236 "ssyconv.f"
			}
#line 237 "ssyconv.f"
		    }
#line 238 "ssyconv.f"
		} else {
#line 239 "ssyconv.f"
		    ip = -ipiv[i__];
#line 240 "ssyconv.f"
		    ++i__;
#line 241 "ssyconv.f"
		    if (i__ < *n) {
#line 242 "ssyconv.f"
			i__1 = *n;
#line 242 "ssyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 243 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 244 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 245 "ssyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 246 "ssyconv.f"
			}
#line 247 "ssyconv.f"
		    }
#line 248 "ssyconv.f"
		}
#line 249 "ssyconv.f"
		++i__;
#line 250 "ssyconv.f"
	    }

/*        Revert VALUE */

#line 254 "ssyconv.f"
	    i__ = *n;
#line 255 "ssyconv.f"
	    while(i__ > 1) {
#line 256 "ssyconv.f"
		if (ipiv[i__] < 0) {
#line 257 "ssyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = work[i__];
#line 258 "ssyconv.f"
		    --i__;
#line 259 "ssyconv.f"
		}
#line 260 "ssyconv.f"
		--i__;
#line 261 "ssyconv.f"
	    }
#line 262 "ssyconv.f"
	}
#line 263 "ssyconv.f"
    } else {

/*      A is LOWER */

#line 267 "ssyconv.f"
	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

#line 274 "ssyconv.f"
	    i__ = 1;
#line 275 "ssyconv.f"
	    work[*n] = 0.;
#line 276 "ssyconv.f"
	    while(i__ <= *n) {
#line 277 "ssyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 278 "ssyconv.f"
		    work[i__] = a[i__ + 1 + i__ * a_dim1];
#line 279 "ssyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = 0.;
#line 280 "ssyconv.f"
		    ++i__;
#line 281 "ssyconv.f"
		} else {
#line 282 "ssyconv.f"
		    work[i__] = 0.;
#line 283 "ssyconv.f"
		}
#line 284 "ssyconv.f"
		++i__;
#line 285 "ssyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 289 "ssyconv.f"
	    i__ = 1;
#line 290 "ssyconv.f"
	    while(i__ <= *n) {
#line 291 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 292 "ssyconv.f"
		    ip = ipiv[i__];
#line 293 "ssyconv.f"
		    if (i__ > 1) {
#line 294 "ssyconv.f"
			i__1 = i__ - 1;
#line 294 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 295 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 296 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 297 "ssyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 298 "ssyconv.f"
/* L22: */
#line 298 "ssyconv.f"
			}
#line 299 "ssyconv.f"
		    }
#line 300 "ssyconv.f"
		} else {
#line 301 "ssyconv.f"
		    ip = -ipiv[i__];
#line 302 "ssyconv.f"
		    if (i__ > 1) {
#line 303 "ssyconv.f"
			i__1 = i__ - 1;
#line 303 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 304 "ssyconv.f"
			    temp = a[ip + j * a_dim1];
#line 305 "ssyconv.f"
			    a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
#line 306 "ssyconv.f"
			    a[i__ + 1 + j * a_dim1] = temp;
#line 307 "ssyconv.f"
/* L23: */
#line 307 "ssyconv.f"
			}
#line 308 "ssyconv.f"
		    }
#line 309 "ssyconv.f"
		    ++i__;
#line 310 "ssyconv.f"
		}
#line 311 "ssyconv.f"
		++i__;
#line 312 "ssyconv.f"
	    }
#line 313 "ssyconv.f"
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

#line 320 "ssyconv.f"
	    i__ = *n;
#line 321 "ssyconv.f"
	    while(i__ >= 1) {
#line 322 "ssyconv.f"
		if (ipiv[i__] > 0) {
#line 323 "ssyconv.f"
		    ip = ipiv[i__];
#line 324 "ssyconv.f"
		    if (i__ > 1) {
#line 325 "ssyconv.f"
			i__1 = i__ - 1;
#line 325 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 326 "ssyconv.f"
			    temp = a[i__ + j * a_dim1];
#line 327 "ssyconv.f"
			    a[i__ + j * a_dim1] = a[ip + j * a_dim1];
#line 328 "ssyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 329 "ssyconv.f"
			}
#line 330 "ssyconv.f"
		    }
#line 331 "ssyconv.f"
		} else {
#line 332 "ssyconv.f"
		    ip = -ipiv[i__];
#line 333 "ssyconv.f"
		    --i__;
#line 334 "ssyconv.f"
		    if (i__ > 1) {
#line 335 "ssyconv.f"
			i__1 = i__ - 1;
#line 335 "ssyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 336 "ssyconv.f"
			    temp = a[i__ + 1 + j * a_dim1];
#line 337 "ssyconv.f"
			    a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
#line 338 "ssyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 339 "ssyconv.f"
			}
#line 340 "ssyconv.f"
		    }
#line 341 "ssyconv.f"
		}
#line 342 "ssyconv.f"
		--i__;
#line 343 "ssyconv.f"
	    }

/*        Revert VALUE */

#line 347 "ssyconv.f"
	    i__ = 1;
#line 348 "ssyconv.f"
	    while(i__ <= *n - 1) {
#line 349 "ssyconv.f"
		if (ipiv[i__] < 0) {
#line 350 "ssyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = work[i__];
#line 351 "ssyconv.f"
		    ++i__;
#line 352 "ssyconv.f"
		}
#line 353 "ssyconv.f"
		++i__;
#line 354 "ssyconv.f"
	    }
#line 355 "ssyconv.f"
	}
#line 356 "ssyconv.f"
    }
#line 358 "ssyconv.f"
    return 0;

/*     End of SSYCONV */

} /* ssyconv_ */

