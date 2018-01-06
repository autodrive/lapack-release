#line 1 "dsyconv.f"
/* dsyconv.f -- translated by f2c (version 20100827).
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

#line 1 "dsyconv.f"
/* > \brief \b DSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYCONV convert A given by TRF into L and D and vice-versa. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSYTRF. */
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
/* >          as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsyconv_(char *uplo, char *way, integer *n, doublereal *
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

#line 148 "dsyconv.f"
    /* Parameter adjustments */
#line 148 "dsyconv.f"
    a_dim1 = *lda;
#line 148 "dsyconv.f"
    a_offset = 1 + a_dim1;
#line 148 "dsyconv.f"
    a -= a_offset;
#line 148 "dsyconv.f"
    --ipiv;
#line 148 "dsyconv.f"
    --work;
#line 148 "dsyconv.f"

#line 148 "dsyconv.f"
    /* Function Body */
#line 148 "dsyconv.f"
    *info = 0;
#line 149 "dsyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 150 "dsyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 151 "dsyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 152 "dsyconv.f"
	*info = -1;
#line 153 "dsyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 154 "dsyconv.f"
	*info = -2;
#line 155 "dsyconv.f"
    } else if (*n < 0) {
#line 156 "dsyconv.f"
	*info = -3;
#line 157 "dsyconv.f"
    } else if (*lda < max(1,*n)) {
#line 158 "dsyconv.f"
	*info = -5;
#line 160 "dsyconv.f"
    }
#line 161 "dsyconv.f"
    if (*info != 0) {
#line 162 "dsyconv.f"
	i__1 = -(*info);
#line 162 "dsyconv.f"
	xerbla_("DSYCONV", &i__1, (ftnlen)7);
#line 163 "dsyconv.f"
	return 0;
#line 164 "dsyconv.f"
    }

/*     Quick return if possible */

#line 168 "dsyconv.f"
    if (*n == 0) {
#line 168 "dsyconv.f"
	return 0;
#line 168 "dsyconv.f"
    }

#line 171 "dsyconv.f"
    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

#line 179 "dsyconv.f"
	if (convert) {
#line 180 "dsyconv.f"
	    i__ = *n;
#line 181 "dsyconv.f"
	    work[1] = 0.;
#line 182 "dsyconv.f"
	    while(i__ > 1) {
#line 183 "dsyconv.f"
		if (ipiv[i__] < 0) {
#line 184 "dsyconv.f"
		    work[i__] = a[i__ - 1 + i__ * a_dim1];
#line 185 "dsyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = 0.;
#line 186 "dsyconv.f"
		    --i__;
#line 187 "dsyconv.f"
		} else {
#line 188 "dsyconv.f"
		    work[i__] = 0.;
#line 189 "dsyconv.f"
		}
#line 190 "dsyconv.f"
		--i__;
#line 191 "dsyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 195 "dsyconv.f"
	    i__ = *n;
#line 196 "dsyconv.f"
	    while(i__ >= 1) {
#line 197 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 198 "dsyconv.f"
		    ip = ipiv[i__];
#line 199 "dsyconv.f"
		    if (i__ < *n) {
#line 200 "dsyconv.f"
			i__1 = *n;
#line 200 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 201 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 202 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 203 "dsyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 204 "dsyconv.f"
/* L12: */
#line 204 "dsyconv.f"
			}
#line 205 "dsyconv.f"
		    }
#line 206 "dsyconv.f"
		} else {
#line 207 "dsyconv.f"
		    ip = -ipiv[i__];
#line 208 "dsyconv.f"
		    if (i__ < *n) {
#line 209 "dsyconv.f"
			i__1 = *n;
#line 209 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 210 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 211 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 212 "dsyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 213 "dsyconv.f"
/* L13: */
#line 213 "dsyconv.f"
			}
#line 214 "dsyconv.f"
		    }
#line 215 "dsyconv.f"
		    --i__;
#line 216 "dsyconv.f"
		}
#line 217 "dsyconv.f"
		--i__;
#line 218 "dsyconv.f"
	    }
#line 220 "dsyconv.f"
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

#line 227 "dsyconv.f"
	    i__ = 1;
#line 228 "dsyconv.f"
	    while(i__ <= *n) {
#line 229 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 230 "dsyconv.f"
		    ip = ipiv[i__];
#line 231 "dsyconv.f"
		    if (i__ < *n) {
#line 232 "dsyconv.f"
			i__1 = *n;
#line 232 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 233 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 234 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 235 "dsyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 236 "dsyconv.f"
			}
#line 237 "dsyconv.f"
		    }
#line 238 "dsyconv.f"
		} else {
#line 239 "dsyconv.f"
		    ip = -ipiv[i__];
#line 240 "dsyconv.f"
		    ++i__;
#line 241 "dsyconv.f"
		    if (i__ < *n) {
#line 242 "dsyconv.f"
			i__1 = *n;
#line 242 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 243 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 244 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 245 "dsyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 246 "dsyconv.f"
			}
#line 247 "dsyconv.f"
		    }
#line 248 "dsyconv.f"
		}
#line 249 "dsyconv.f"
		++i__;
#line 250 "dsyconv.f"
	    }

/*        Revert VALUE */

#line 254 "dsyconv.f"
	    i__ = *n;
#line 255 "dsyconv.f"
	    while(i__ > 1) {
#line 256 "dsyconv.f"
		if (ipiv[i__] < 0) {
#line 257 "dsyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = work[i__];
#line 258 "dsyconv.f"
		    --i__;
#line 259 "dsyconv.f"
		}
#line 260 "dsyconv.f"
		--i__;
#line 261 "dsyconv.f"
	    }
#line 262 "dsyconv.f"
	}
#line 263 "dsyconv.f"
    } else {

/*      A is LOWER */

#line 267 "dsyconv.f"
	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

#line 274 "dsyconv.f"
	    i__ = 1;
#line 275 "dsyconv.f"
	    work[*n] = 0.;
#line 276 "dsyconv.f"
	    while(i__ <= *n) {
#line 277 "dsyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 278 "dsyconv.f"
		    work[i__] = a[i__ + 1 + i__ * a_dim1];
#line 279 "dsyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = 0.;
#line 280 "dsyconv.f"
		    ++i__;
#line 281 "dsyconv.f"
		} else {
#line 282 "dsyconv.f"
		    work[i__] = 0.;
#line 283 "dsyconv.f"
		}
#line 284 "dsyconv.f"
		++i__;
#line 285 "dsyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 289 "dsyconv.f"
	    i__ = 1;
#line 290 "dsyconv.f"
	    while(i__ <= *n) {
#line 291 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 292 "dsyconv.f"
		    ip = ipiv[i__];
#line 293 "dsyconv.f"
		    if (i__ > 1) {
#line 294 "dsyconv.f"
			i__1 = i__ - 1;
#line 294 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 295 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 296 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 297 "dsyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 298 "dsyconv.f"
/* L22: */
#line 298 "dsyconv.f"
			}
#line 299 "dsyconv.f"
		    }
#line 300 "dsyconv.f"
		} else {
#line 301 "dsyconv.f"
		    ip = -ipiv[i__];
#line 302 "dsyconv.f"
		    if (i__ > 1) {
#line 303 "dsyconv.f"
			i__1 = i__ - 1;
#line 303 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 304 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 305 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
#line 306 "dsyconv.f"
			    a[i__ + 1 + j * a_dim1] = temp;
#line 307 "dsyconv.f"
/* L23: */
#line 307 "dsyconv.f"
			}
#line 308 "dsyconv.f"
		    }
#line 309 "dsyconv.f"
		    ++i__;
#line 310 "dsyconv.f"
		}
#line 311 "dsyconv.f"
		++i__;
#line 312 "dsyconv.f"
	    }
#line 313 "dsyconv.f"
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

#line 320 "dsyconv.f"
	    i__ = *n;
#line 321 "dsyconv.f"
	    while(i__ >= 1) {
#line 322 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 323 "dsyconv.f"
		    ip = ipiv[i__];
#line 324 "dsyconv.f"
		    if (i__ > 1) {
#line 325 "dsyconv.f"
			i__1 = i__ - 1;
#line 325 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 326 "dsyconv.f"
			    temp = a[i__ + j * a_dim1];
#line 327 "dsyconv.f"
			    a[i__ + j * a_dim1] = a[ip + j * a_dim1];
#line 328 "dsyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 329 "dsyconv.f"
			}
#line 330 "dsyconv.f"
		    }
#line 331 "dsyconv.f"
		} else {
#line 332 "dsyconv.f"
		    ip = -ipiv[i__];
#line 333 "dsyconv.f"
		    --i__;
#line 334 "dsyconv.f"
		    if (i__ > 1) {
#line 335 "dsyconv.f"
			i__1 = i__ - 1;
#line 335 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 336 "dsyconv.f"
			    temp = a[i__ + 1 + j * a_dim1];
#line 337 "dsyconv.f"
			    a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
#line 338 "dsyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 339 "dsyconv.f"
			}
#line 340 "dsyconv.f"
		    }
#line 341 "dsyconv.f"
		}
#line 342 "dsyconv.f"
		--i__;
#line 343 "dsyconv.f"
	    }

/*        Revert VALUE */

#line 347 "dsyconv.f"
	    i__ = 1;
#line 348 "dsyconv.f"
	    while(i__ <= *n - 1) {
#line 349 "dsyconv.f"
		if ((doublereal) ipiv[i__] < 0.) {
#line 350 "dsyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = work[i__];
#line 351 "dsyconv.f"
		    ++i__;
#line 352 "dsyconv.f"
		}
#line 353 "dsyconv.f"
		++i__;
#line 354 "dsyconv.f"
	    }
#line 355 "dsyconv.f"
	}
#line 356 "dsyconv.f"
    }
#line 358 "dsyconv.f"
    return 0;

/*     End of DSYCONV */

} /* dsyconv_ */

