#line 1 "slantp.f"
/* slantp.f -- translated by f2c (version 20100827).
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

#line 1 "slantp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANTP( NORM, UPLO, DIAG, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANTP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > triangular matrix A, supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return SLANTP */
/* > \verbatim */
/* > */
/* >    SLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
/* >             ( */
/* >             ( norm1(A),         NORM = '1', 'O' or 'o' */
/* >             ( */
/* >             ( normI(A),         NORM = 'I' or 'i' */
/* >             ( */
/* >             ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/* > normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/* > normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/* > squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies the value to be returned in SLANTP as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower triangular. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A is unit triangular. */
/* >          = 'N':  Non-unit triangular */
/* >          = 'U':  Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANTP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* >          Note that when DIAG = 'U', the elements of the array AP */
/* >          corresponding to the diagonal elements of the matrix A are */
/* >          not referenced, but are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= N when NORM = 'I'; otherwise, WORK is not */
/* >          referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
doublereal slantp_(char *norm, char *uplo, char *diag, integer *n, doublereal 
	*ap, doublereal *work, ftnlen norm_len, ftnlen uplo_len, ftnlen 
	diag_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum, scale;
    static logical udiag;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 163 "slantp.f"
    /* Parameter adjustments */
#line 163 "slantp.f"
    --work;
#line 163 "slantp.f"
    --ap;
#line 163 "slantp.f"

#line 163 "slantp.f"
    /* Function Body */
#line 163 "slantp.f"
    if (*n == 0) {
#line 164 "slantp.f"
	value = 0.;
#line 165 "slantp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 169 "slantp.f"
	k = 1;
#line 170 "slantp.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 171 "slantp.f"
	    value = 1.;
#line 172 "slantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 173 "slantp.f"
		i__1 = *n;
#line 173 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 174 "slantp.f"
		    i__2 = k + j - 2;
#line 174 "slantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 175 "slantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 176 "slantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 176 "slantp.f"
			    value = sum;
#line 176 "slantp.f"
			}
#line 177 "slantp.f"
/* L10: */
#line 177 "slantp.f"
		    }
#line 178 "slantp.f"
		    k += j;
#line 179 "slantp.f"
/* L20: */
#line 179 "slantp.f"
		}
#line 180 "slantp.f"
	    } else {
#line 181 "slantp.f"
		i__1 = *n;
#line 181 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 182 "slantp.f"
		    i__2 = k + *n - j;
#line 182 "slantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 183 "slantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 184 "slantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 184 "slantp.f"
			    value = sum;
#line 184 "slantp.f"
			}
#line 185 "slantp.f"
/* L30: */
#line 185 "slantp.f"
		    }
#line 186 "slantp.f"
		    k = k + *n - j + 1;
#line 187 "slantp.f"
/* L40: */
#line 187 "slantp.f"
		}
#line 188 "slantp.f"
	    }
#line 189 "slantp.f"
	} else {
#line 190 "slantp.f"
	    value = 0.;
#line 191 "slantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 192 "slantp.f"
		i__1 = *n;
#line 192 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 193 "slantp.f"
		    i__2 = k + j - 1;
#line 193 "slantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 194 "slantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 195 "slantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 195 "slantp.f"
			    value = sum;
#line 195 "slantp.f"
			}
#line 196 "slantp.f"
/* L50: */
#line 196 "slantp.f"
		    }
#line 197 "slantp.f"
		    k += j;
#line 198 "slantp.f"
/* L60: */
#line 198 "slantp.f"
		}
#line 199 "slantp.f"
	    } else {
#line 200 "slantp.f"
		i__1 = *n;
#line 200 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 201 "slantp.f"
		    i__2 = k + *n - j;
#line 201 "slantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 202 "slantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 203 "slantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 203 "slantp.f"
			    value = sum;
#line 203 "slantp.f"
			}
#line 204 "slantp.f"
/* L70: */
#line 204 "slantp.f"
		    }
#line 205 "slantp.f"
		    k = k + *n - j + 1;
#line 206 "slantp.f"
/* L80: */
#line 206 "slantp.f"
		}
#line 207 "slantp.f"
	    }
#line 208 "slantp.f"
	}
#line 209 "slantp.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 213 "slantp.f"
	value = 0.;
#line 214 "slantp.f"
	k = 1;
#line 215 "slantp.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 216 "slantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 217 "slantp.f"
	    i__1 = *n;
#line 217 "slantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 218 "slantp.f"
		if (udiag) {
#line 219 "slantp.f"
		    sum = 1.;
#line 220 "slantp.f"
		    i__2 = k + j - 2;
#line 220 "slantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 221 "slantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 222 "slantp.f"
/* L90: */
#line 222 "slantp.f"
		    }
#line 223 "slantp.f"
		} else {
#line 224 "slantp.f"
		    sum = 0.;
#line 225 "slantp.f"
		    i__2 = k + j - 1;
#line 225 "slantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 226 "slantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 227 "slantp.f"
/* L100: */
#line 227 "slantp.f"
		    }
#line 228 "slantp.f"
		}
#line 229 "slantp.f"
		k += j;
#line 230 "slantp.f"
		if (value < sum || sisnan_(&sum)) {
#line 230 "slantp.f"
		    value = sum;
#line 230 "slantp.f"
		}
#line 231 "slantp.f"
/* L110: */
#line 231 "slantp.f"
	    }
#line 232 "slantp.f"
	} else {
#line 233 "slantp.f"
	    i__1 = *n;
#line 233 "slantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 234 "slantp.f"
		if (udiag) {
#line 235 "slantp.f"
		    sum = 1.;
#line 236 "slantp.f"
		    i__2 = k + *n - j;
#line 236 "slantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 237 "slantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 238 "slantp.f"
/* L120: */
#line 238 "slantp.f"
		    }
#line 239 "slantp.f"
		} else {
#line 240 "slantp.f"
		    sum = 0.;
#line 241 "slantp.f"
		    i__2 = k + *n - j;
#line 241 "slantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 242 "slantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 243 "slantp.f"
/* L130: */
#line 243 "slantp.f"
		    }
#line 244 "slantp.f"
		}
#line 245 "slantp.f"
		k = k + *n - j + 1;
#line 246 "slantp.f"
		if (value < sum || sisnan_(&sum)) {
#line 246 "slantp.f"
		    value = sum;
#line 246 "slantp.f"
		}
#line 247 "slantp.f"
/* L140: */
#line 247 "slantp.f"
	    }
#line 248 "slantp.f"
	}
#line 249 "slantp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 253 "slantp.f"
	k = 1;
#line 254 "slantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 255 "slantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "slantp.f"
		i__1 = *n;
#line 256 "slantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 257 "slantp.f"
		    work[i__] = 1.;
#line 258 "slantp.f"
/* L150: */
#line 258 "slantp.f"
		}
#line 259 "slantp.f"
		i__1 = *n;
#line 259 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "slantp.f"
		    i__2 = j - 1;
#line 260 "slantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 261 "slantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 262 "slantp.f"
			++k;
#line 263 "slantp.f"
/* L160: */
#line 263 "slantp.f"
		    }
#line 264 "slantp.f"
		    ++k;
#line 265 "slantp.f"
/* L170: */
#line 265 "slantp.f"
		}
#line 266 "slantp.f"
	    } else {
#line 267 "slantp.f"
		i__1 = *n;
#line 267 "slantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 268 "slantp.f"
		    work[i__] = 0.;
#line 269 "slantp.f"
/* L180: */
#line 269 "slantp.f"
		}
#line 270 "slantp.f"
		i__1 = *n;
#line 270 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 271 "slantp.f"
		    i__2 = j;
#line 271 "slantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 272 "slantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 273 "slantp.f"
			++k;
#line 274 "slantp.f"
/* L190: */
#line 274 "slantp.f"
		    }
#line 275 "slantp.f"
/* L200: */
#line 275 "slantp.f"
		}
#line 276 "slantp.f"
	    }
#line 277 "slantp.f"
	} else {
#line 278 "slantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 279 "slantp.f"
		i__1 = *n;
#line 279 "slantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "slantp.f"
		    work[i__] = 1.;
#line 281 "slantp.f"
/* L210: */
#line 281 "slantp.f"
		}
#line 282 "slantp.f"
		i__1 = *n;
#line 282 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 283 "slantp.f"
		    ++k;
#line 284 "slantp.f"
		    i__2 = *n;
#line 284 "slantp.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 285 "slantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 286 "slantp.f"
			++k;
#line 287 "slantp.f"
/* L220: */
#line 287 "slantp.f"
		    }
#line 288 "slantp.f"
/* L230: */
#line 288 "slantp.f"
		}
#line 289 "slantp.f"
	    } else {
#line 290 "slantp.f"
		i__1 = *n;
#line 290 "slantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "slantp.f"
		    work[i__] = 0.;
#line 292 "slantp.f"
/* L240: */
#line 292 "slantp.f"
		}
#line 293 "slantp.f"
		i__1 = *n;
#line 293 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 294 "slantp.f"
		    i__2 = *n;
#line 294 "slantp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 295 "slantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 296 "slantp.f"
			++k;
#line 297 "slantp.f"
/* L250: */
#line 297 "slantp.f"
		    }
#line 298 "slantp.f"
/* L260: */
#line 298 "slantp.f"
		}
#line 299 "slantp.f"
	    }
#line 300 "slantp.f"
	}
#line 301 "slantp.f"
	value = 0.;
#line 302 "slantp.f"
	i__1 = *n;
#line 302 "slantp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 303 "slantp.f"
	    sum = work[i__];
#line 304 "slantp.f"
	    if (value < sum || sisnan_(&sum)) {
#line 304 "slantp.f"
		value = sum;
#line 304 "slantp.f"
	    }
#line 305 "slantp.f"
/* L270: */
#line 305 "slantp.f"
	}
#line 306 "slantp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 310 "slantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 311 "slantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 312 "slantp.f"
		scale = 1.;
#line 313 "slantp.f"
		sum = (doublereal) (*n);
#line 314 "slantp.f"
		k = 2;
#line 315 "slantp.f"
		i__1 = *n;
#line 315 "slantp.f"
		for (j = 2; j <= i__1; ++j) {
#line 316 "slantp.f"
		    i__2 = j - 1;
#line 316 "slantp.f"
		    slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 317 "slantp.f"
		    k += j;
#line 318 "slantp.f"
/* L280: */
#line 318 "slantp.f"
		}
#line 319 "slantp.f"
	    } else {
#line 320 "slantp.f"
		scale = 0.;
#line 321 "slantp.f"
		sum = 1.;
#line 322 "slantp.f"
		k = 1;
#line 323 "slantp.f"
		i__1 = *n;
#line 323 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 324 "slantp.f"
		    slassq_(&j, &ap[k], &c__1, &scale, &sum);
#line 325 "slantp.f"
		    k += j;
#line 326 "slantp.f"
/* L290: */
#line 326 "slantp.f"
		}
#line 327 "slantp.f"
	    }
#line 328 "slantp.f"
	} else {
#line 329 "slantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 330 "slantp.f"
		scale = 1.;
#line 331 "slantp.f"
		sum = (doublereal) (*n);
#line 332 "slantp.f"
		k = 2;
#line 333 "slantp.f"
		i__1 = *n - 1;
#line 333 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 334 "slantp.f"
		    i__2 = *n - j;
#line 334 "slantp.f"
		    slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 335 "slantp.f"
		    k = k + *n - j + 1;
#line 336 "slantp.f"
/* L300: */
#line 336 "slantp.f"
		}
#line 337 "slantp.f"
	    } else {
#line 338 "slantp.f"
		scale = 0.;
#line 339 "slantp.f"
		sum = 1.;
#line 340 "slantp.f"
		k = 1;
#line 341 "slantp.f"
		i__1 = *n;
#line 341 "slantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 342 "slantp.f"
		    i__2 = *n - j + 1;
#line 342 "slantp.f"
		    slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 343 "slantp.f"
		    k = k + *n - j + 1;
#line 344 "slantp.f"
/* L310: */
#line 344 "slantp.f"
		}
#line 345 "slantp.f"
	    }
#line 346 "slantp.f"
	}
#line 347 "slantp.f"
	value = scale * sqrt(sum);
#line 348 "slantp.f"
    }

#line 350 "slantp.f"
    ret_val = value;
#line 351 "slantp.f"
    return ret_val;

/*     End of SLANTP */

} /* slantp_ */

