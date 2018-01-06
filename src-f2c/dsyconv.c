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

/*       SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), E( * ) */
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
/* > \param[in,out] A */
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
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsyconv_(char *uplo, char *way, integer *n, doublereal *
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

#line 150 "dsyconv.f"
    /* Parameter adjustments */
#line 150 "dsyconv.f"
    a_dim1 = *lda;
#line 150 "dsyconv.f"
    a_offset = 1 + a_dim1;
#line 150 "dsyconv.f"
    a -= a_offset;
#line 150 "dsyconv.f"
    --ipiv;
#line 150 "dsyconv.f"
    --e;
#line 150 "dsyconv.f"

#line 150 "dsyconv.f"
    /* Function Body */
#line 150 "dsyconv.f"
    *info = 0;
#line 151 "dsyconv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 152 "dsyconv.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 153 "dsyconv.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 154 "dsyconv.f"
	*info = -1;
#line 155 "dsyconv.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 156 "dsyconv.f"
	*info = -2;
#line 157 "dsyconv.f"
    } else if (*n < 0) {
#line 158 "dsyconv.f"
	*info = -3;
#line 159 "dsyconv.f"
    } else if (*lda < max(1,*n)) {
#line 160 "dsyconv.f"
	*info = -5;
#line 162 "dsyconv.f"
    }
#line 163 "dsyconv.f"
    if (*info != 0) {
#line 164 "dsyconv.f"
	i__1 = -(*info);
#line 164 "dsyconv.f"
	xerbla_("DSYCONV", &i__1, (ftnlen)7);
#line 165 "dsyconv.f"
	return 0;
#line 166 "dsyconv.f"
    }

/*     Quick return if possible */

#line 170 "dsyconv.f"
    if (*n == 0) {
#line 170 "dsyconv.f"
	return 0;
#line 170 "dsyconv.f"
    }

#line 173 "dsyconv.f"
    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

#line 181 "dsyconv.f"
	if (convert) {
#line 182 "dsyconv.f"
	    i__ = *n;
#line 183 "dsyconv.f"
	    e[1] = 0.;
#line 184 "dsyconv.f"
	    while(i__ > 1) {
#line 185 "dsyconv.f"
		if (ipiv[i__] < 0) {
#line 186 "dsyconv.f"
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
#line 187 "dsyconv.f"
		    e[i__ - 1] = 0.;
#line 188 "dsyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = 0.;
#line 189 "dsyconv.f"
		    --i__;
#line 190 "dsyconv.f"
		} else {
#line 191 "dsyconv.f"
		    e[i__] = 0.;
#line 192 "dsyconv.f"
		}
#line 193 "dsyconv.f"
		--i__;
#line 194 "dsyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 198 "dsyconv.f"
	    i__ = *n;
#line 199 "dsyconv.f"
	    while(i__ >= 1) {
#line 200 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 201 "dsyconv.f"
		    ip = ipiv[i__];
#line 202 "dsyconv.f"
		    if (i__ < *n) {
#line 203 "dsyconv.f"
			i__1 = *n;
#line 203 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 204 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 205 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 206 "dsyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 207 "dsyconv.f"
/* L12: */
#line 207 "dsyconv.f"
			}
#line 208 "dsyconv.f"
		    }
#line 209 "dsyconv.f"
		} else {
#line 210 "dsyconv.f"
		    ip = -ipiv[i__];
#line 211 "dsyconv.f"
		    if (i__ < *n) {
#line 212 "dsyconv.f"
			i__1 = *n;
#line 212 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 213 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 214 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 215 "dsyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 216 "dsyconv.f"
/* L13: */
#line 216 "dsyconv.f"
			}
#line 217 "dsyconv.f"
		    }
#line 218 "dsyconv.f"
		    --i__;
#line 219 "dsyconv.f"
		}
#line 220 "dsyconv.f"
		--i__;
#line 221 "dsyconv.f"
	    }
#line 223 "dsyconv.f"
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

#line 230 "dsyconv.f"
	    i__ = 1;
#line 231 "dsyconv.f"
	    while(i__ <= *n) {
#line 232 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 233 "dsyconv.f"
		    ip = ipiv[i__];
#line 234 "dsyconv.f"
		    if (i__ < *n) {
#line 235 "dsyconv.f"
			i__1 = *n;
#line 235 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 236 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 237 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 238 "dsyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 239 "dsyconv.f"
			}
#line 240 "dsyconv.f"
		    }
#line 241 "dsyconv.f"
		} else {
#line 242 "dsyconv.f"
		    ip = -ipiv[i__];
#line 243 "dsyconv.f"
		    ++i__;
#line 244 "dsyconv.f"
		    if (i__ < *n) {
#line 245 "dsyconv.f"
			i__1 = *n;
#line 245 "dsyconv.f"
			for (j = i__ + 1; j <= i__1; ++j) {
#line 246 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 247 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 248 "dsyconv.f"
			    a[i__ - 1 + j * a_dim1] = temp;
#line 249 "dsyconv.f"
			}
#line 250 "dsyconv.f"
		    }
#line 251 "dsyconv.f"
		}
#line 252 "dsyconv.f"
		++i__;
#line 253 "dsyconv.f"
	    }

/*        Revert VALUE */

#line 257 "dsyconv.f"
	    i__ = *n;
#line 258 "dsyconv.f"
	    while(i__ > 1) {
#line 259 "dsyconv.f"
		if (ipiv[i__] < 0) {
#line 260 "dsyconv.f"
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
#line 261 "dsyconv.f"
		    --i__;
#line 262 "dsyconv.f"
		}
#line 263 "dsyconv.f"
		--i__;
#line 264 "dsyconv.f"
	    }
#line 265 "dsyconv.f"
	}
#line 266 "dsyconv.f"
    } else {

/*      A is LOWER */

#line 270 "dsyconv.f"
	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

#line 277 "dsyconv.f"
	    i__ = 1;
#line 278 "dsyconv.f"
	    e[*n] = 0.;
#line 279 "dsyconv.f"
	    while(i__ <= *n) {
#line 280 "dsyconv.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 281 "dsyconv.f"
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 282 "dsyconv.f"
		    e[i__ + 1] = 0.;
#line 283 "dsyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = 0.;
#line 284 "dsyconv.f"
		    ++i__;
#line 285 "dsyconv.f"
		} else {
#line 286 "dsyconv.f"
		    e[i__] = 0.;
#line 287 "dsyconv.f"
		}
#line 288 "dsyconv.f"
		++i__;
#line 289 "dsyconv.f"
	    }

/*        Convert PERMUTATIONS */

#line 293 "dsyconv.f"
	    i__ = 1;
#line 294 "dsyconv.f"
	    while(i__ <= *n) {
#line 295 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 296 "dsyconv.f"
		    ip = ipiv[i__];
#line 297 "dsyconv.f"
		    if (i__ > 1) {
#line 298 "dsyconv.f"
			i__1 = i__ - 1;
#line 298 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 299 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 300 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
#line 301 "dsyconv.f"
			    a[i__ + j * a_dim1] = temp;
#line 302 "dsyconv.f"
/* L22: */
#line 302 "dsyconv.f"
			}
#line 303 "dsyconv.f"
		    }
#line 304 "dsyconv.f"
		} else {
#line 305 "dsyconv.f"
		    ip = -ipiv[i__];
#line 306 "dsyconv.f"
		    if (i__ > 1) {
#line 307 "dsyconv.f"
			i__1 = i__ - 1;
#line 307 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 308 "dsyconv.f"
			    temp = a[ip + j * a_dim1];
#line 309 "dsyconv.f"
			    a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
#line 310 "dsyconv.f"
			    a[i__ + 1 + j * a_dim1] = temp;
#line 311 "dsyconv.f"
/* L23: */
#line 311 "dsyconv.f"
			}
#line 312 "dsyconv.f"
		    }
#line 313 "dsyconv.f"
		    ++i__;
#line 314 "dsyconv.f"
		}
#line 315 "dsyconv.f"
		++i__;
#line 316 "dsyconv.f"
	    }
#line 317 "dsyconv.f"
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

#line 324 "dsyconv.f"
	    i__ = *n;
#line 325 "dsyconv.f"
	    while(i__ >= 1) {
#line 326 "dsyconv.f"
		if (ipiv[i__] > 0) {
#line 327 "dsyconv.f"
		    ip = ipiv[i__];
#line 328 "dsyconv.f"
		    if (i__ > 1) {
#line 329 "dsyconv.f"
			i__1 = i__ - 1;
#line 329 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 330 "dsyconv.f"
			    temp = a[i__ + j * a_dim1];
#line 331 "dsyconv.f"
			    a[i__ + j * a_dim1] = a[ip + j * a_dim1];
#line 332 "dsyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 333 "dsyconv.f"
			}
#line 334 "dsyconv.f"
		    }
#line 335 "dsyconv.f"
		} else {
#line 336 "dsyconv.f"
		    ip = -ipiv[i__];
#line 337 "dsyconv.f"
		    --i__;
#line 338 "dsyconv.f"
		    if (i__ > 1) {
#line 339 "dsyconv.f"
			i__1 = i__ - 1;
#line 339 "dsyconv.f"
			for (j = 1; j <= i__1; ++j) {
#line 340 "dsyconv.f"
			    temp = a[i__ + 1 + j * a_dim1];
#line 341 "dsyconv.f"
			    a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
#line 342 "dsyconv.f"
			    a[ip + j * a_dim1] = temp;
#line 343 "dsyconv.f"
			}
#line 344 "dsyconv.f"
		    }
#line 345 "dsyconv.f"
		}
#line 346 "dsyconv.f"
		--i__;
#line 347 "dsyconv.f"
	    }

/*        Revert VALUE */

#line 351 "dsyconv.f"
	    i__ = 1;
#line 352 "dsyconv.f"
	    while(i__ <= *n - 1) {
#line 353 "dsyconv.f"
		if (ipiv[i__] < 0) {
#line 354 "dsyconv.f"
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 355 "dsyconv.f"
		    ++i__;
#line 356 "dsyconv.f"
		}
#line 357 "dsyconv.f"
		++i__;
#line 358 "dsyconv.f"
	    }
#line 359 "dsyconv.f"
	}
#line 360 "dsyconv.f"
    }
#line 362 "dsyconv.f"
    return 0;

/*     End of DSYCONV */

} /* dsyconv_ */

