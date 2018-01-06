#line 1 "dlantp.f"
/* dlantp.f -- translated by f2c (version 20100827).
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

#line 1 "dlantp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANTP( NORM, UPLO, DIAG, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANTP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > triangular matrix A, supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return DLANTP */
/* > \verbatim */
/* > */
/* >    DLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANTP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANTP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
doublereal dlantp_(char *norm, char *uplo, char *diag, integer *n, doublereal 
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
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
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

#line 163 "dlantp.f"
    /* Parameter adjustments */
#line 163 "dlantp.f"
    --work;
#line 163 "dlantp.f"
    --ap;
#line 163 "dlantp.f"

#line 163 "dlantp.f"
    /* Function Body */
#line 163 "dlantp.f"
    if (*n == 0) {
#line 164 "dlantp.f"
	value = 0.;
#line 165 "dlantp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 169 "dlantp.f"
	k = 1;
#line 170 "dlantp.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 171 "dlantp.f"
	    value = 1.;
#line 172 "dlantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 173 "dlantp.f"
		i__1 = *n;
#line 173 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 174 "dlantp.f"
		    i__2 = k + j - 2;
#line 174 "dlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 175 "dlantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 176 "dlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 176 "dlantp.f"
			    value = sum;
#line 176 "dlantp.f"
			}
#line 177 "dlantp.f"
/* L10: */
#line 177 "dlantp.f"
		    }
#line 178 "dlantp.f"
		    k += j;
#line 179 "dlantp.f"
/* L20: */
#line 179 "dlantp.f"
		}
#line 180 "dlantp.f"
	    } else {
#line 181 "dlantp.f"
		i__1 = *n;
#line 181 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 182 "dlantp.f"
		    i__2 = k + *n - j;
#line 182 "dlantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 183 "dlantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 184 "dlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 184 "dlantp.f"
			    value = sum;
#line 184 "dlantp.f"
			}
#line 185 "dlantp.f"
/* L30: */
#line 185 "dlantp.f"
		    }
#line 186 "dlantp.f"
		    k = k + *n - j + 1;
#line 187 "dlantp.f"
/* L40: */
#line 187 "dlantp.f"
		}
#line 188 "dlantp.f"
	    }
#line 189 "dlantp.f"
	} else {
#line 190 "dlantp.f"
	    value = 0.;
#line 191 "dlantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 192 "dlantp.f"
		i__1 = *n;
#line 192 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 193 "dlantp.f"
		    i__2 = k + j - 1;
#line 193 "dlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 194 "dlantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 195 "dlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 195 "dlantp.f"
			    value = sum;
#line 195 "dlantp.f"
			}
#line 196 "dlantp.f"
/* L50: */
#line 196 "dlantp.f"
		    }
#line 197 "dlantp.f"
		    k += j;
#line 198 "dlantp.f"
/* L60: */
#line 198 "dlantp.f"
		}
#line 199 "dlantp.f"
	    } else {
#line 200 "dlantp.f"
		i__1 = *n;
#line 200 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 201 "dlantp.f"
		    i__2 = k + *n - j;
#line 201 "dlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 202 "dlantp.f"
			sum = (d__1 = ap[i__], abs(d__1));
#line 203 "dlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 203 "dlantp.f"
			    value = sum;
#line 203 "dlantp.f"
			}
#line 204 "dlantp.f"
/* L70: */
#line 204 "dlantp.f"
		    }
#line 205 "dlantp.f"
		    k = k + *n - j + 1;
#line 206 "dlantp.f"
/* L80: */
#line 206 "dlantp.f"
		}
#line 207 "dlantp.f"
	    }
#line 208 "dlantp.f"
	}
#line 209 "dlantp.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 213 "dlantp.f"
	value = 0.;
#line 214 "dlantp.f"
	k = 1;
#line 215 "dlantp.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 216 "dlantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 217 "dlantp.f"
	    i__1 = *n;
#line 217 "dlantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 218 "dlantp.f"
		if (udiag) {
#line 219 "dlantp.f"
		    sum = 1.;
#line 220 "dlantp.f"
		    i__2 = k + j - 2;
#line 220 "dlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 221 "dlantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 222 "dlantp.f"
/* L90: */
#line 222 "dlantp.f"
		    }
#line 223 "dlantp.f"
		} else {
#line 224 "dlantp.f"
		    sum = 0.;
#line 225 "dlantp.f"
		    i__2 = k + j - 1;
#line 225 "dlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 226 "dlantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 227 "dlantp.f"
/* L100: */
#line 227 "dlantp.f"
		    }
#line 228 "dlantp.f"
		}
#line 229 "dlantp.f"
		k += j;
#line 230 "dlantp.f"
		if (value < sum || disnan_(&sum)) {
#line 230 "dlantp.f"
		    value = sum;
#line 230 "dlantp.f"
		}
#line 231 "dlantp.f"
/* L110: */
#line 231 "dlantp.f"
	    }
#line 232 "dlantp.f"
	} else {
#line 233 "dlantp.f"
	    i__1 = *n;
#line 233 "dlantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 234 "dlantp.f"
		if (udiag) {
#line 235 "dlantp.f"
		    sum = 1.;
#line 236 "dlantp.f"
		    i__2 = k + *n - j;
#line 236 "dlantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 237 "dlantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 238 "dlantp.f"
/* L120: */
#line 238 "dlantp.f"
		    }
#line 239 "dlantp.f"
		} else {
#line 240 "dlantp.f"
		    sum = 0.;
#line 241 "dlantp.f"
		    i__2 = k + *n - j;
#line 241 "dlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 242 "dlantp.f"
			sum += (d__1 = ap[i__], abs(d__1));
#line 243 "dlantp.f"
/* L130: */
#line 243 "dlantp.f"
		    }
#line 244 "dlantp.f"
		}
#line 245 "dlantp.f"
		k = k + *n - j + 1;
#line 246 "dlantp.f"
		if (value < sum || disnan_(&sum)) {
#line 246 "dlantp.f"
		    value = sum;
#line 246 "dlantp.f"
		}
#line 247 "dlantp.f"
/* L140: */
#line 247 "dlantp.f"
	    }
#line 248 "dlantp.f"
	}
#line 249 "dlantp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 253 "dlantp.f"
	k = 1;
#line 254 "dlantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 255 "dlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "dlantp.f"
		i__1 = *n;
#line 256 "dlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 257 "dlantp.f"
		    work[i__] = 1.;
#line 258 "dlantp.f"
/* L150: */
#line 258 "dlantp.f"
		}
#line 259 "dlantp.f"
		i__1 = *n;
#line 259 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "dlantp.f"
		    i__2 = j - 1;
#line 260 "dlantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 261 "dlantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 262 "dlantp.f"
			++k;
#line 263 "dlantp.f"
/* L160: */
#line 263 "dlantp.f"
		    }
#line 264 "dlantp.f"
		    ++k;
#line 265 "dlantp.f"
/* L170: */
#line 265 "dlantp.f"
		}
#line 266 "dlantp.f"
	    } else {
#line 267 "dlantp.f"
		i__1 = *n;
#line 267 "dlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 268 "dlantp.f"
		    work[i__] = 0.;
#line 269 "dlantp.f"
/* L180: */
#line 269 "dlantp.f"
		}
#line 270 "dlantp.f"
		i__1 = *n;
#line 270 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 271 "dlantp.f"
		    i__2 = j;
#line 271 "dlantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 272 "dlantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 273 "dlantp.f"
			++k;
#line 274 "dlantp.f"
/* L190: */
#line 274 "dlantp.f"
		    }
#line 275 "dlantp.f"
/* L200: */
#line 275 "dlantp.f"
		}
#line 276 "dlantp.f"
	    }
#line 277 "dlantp.f"
	} else {
#line 278 "dlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 279 "dlantp.f"
		i__1 = *n;
#line 279 "dlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "dlantp.f"
		    work[i__] = 1.;
#line 281 "dlantp.f"
/* L210: */
#line 281 "dlantp.f"
		}
#line 282 "dlantp.f"
		i__1 = *n;
#line 282 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 283 "dlantp.f"
		    ++k;
#line 284 "dlantp.f"
		    i__2 = *n;
#line 284 "dlantp.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 285 "dlantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 286 "dlantp.f"
			++k;
#line 287 "dlantp.f"
/* L220: */
#line 287 "dlantp.f"
		    }
#line 288 "dlantp.f"
/* L230: */
#line 288 "dlantp.f"
		}
#line 289 "dlantp.f"
	    } else {
#line 290 "dlantp.f"
		i__1 = *n;
#line 290 "dlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "dlantp.f"
		    work[i__] = 0.;
#line 292 "dlantp.f"
/* L240: */
#line 292 "dlantp.f"
		}
#line 293 "dlantp.f"
		i__1 = *n;
#line 293 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 294 "dlantp.f"
		    i__2 = *n;
#line 294 "dlantp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 295 "dlantp.f"
			work[i__] += (d__1 = ap[k], abs(d__1));
#line 296 "dlantp.f"
			++k;
#line 297 "dlantp.f"
/* L250: */
#line 297 "dlantp.f"
		    }
#line 298 "dlantp.f"
/* L260: */
#line 298 "dlantp.f"
		}
#line 299 "dlantp.f"
	    }
#line 300 "dlantp.f"
	}
#line 301 "dlantp.f"
	value = 0.;
#line 302 "dlantp.f"
	i__1 = *n;
#line 302 "dlantp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 303 "dlantp.f"
	    sum = work[i__];
#line 304 "dlantp.f"
	    if (value < sum || disnan_(&sum)) {
#line 304 "dlantp.f"
		value = sum;
#line 304 "dlantp.f"
	    }
#line 305 "dlantp.f"
/* L270: */
#line 305 "dlantp.f"
	}
#line 306 "dlantp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 310 "dlantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 311 "dlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 312 "dlantp.f"
		scale = 1.;
#line 313 "dlantp.f"
		sum = (doublereal) (*n);
#line 314 "dlantp.f"
		k = 2;
#line 315 "dlantp.f"
		i__1 = *n;
#line 315 "dlantp.f"
		for (j = 2; j <= i__1; ++j) {
#line 316 "dlantp.f"
		    i__2 = j - 1;
#line 316 "dlantp.f"
		    dlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 317 "dlantp.f"
		    k += j;
#line 318 "dlantp.f"
/* L280: */
#line 318 "dlantp.f"
		}
#line 319 "dlantp.f"
	    } else {
#line 320 "dlantp.f"
		scale = 0.;
#line 321 "dlantp.f"
		sum = 1.;
#line 322 "dlantp.f"
		k = 1;
#line 323 "dlantp.f"
		i__1 = *n;
#line 323 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 324 "dlantp.f"
		    dlassq_(&j, &ap[k], &c__1, &scale, &sum);
#line 325 "dlantp.f"
		    k += j;
#line 326 "dlantp.f"
/* L290: */
#line 326 "dlantp.f"
		}
#line 327 "dlantp.f"
	    }
#line 328 "dlantp.f"
	} else {
#line 329 "dlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 330 "dlantp.f"
		scale = 1.;
#line 331 "dlantp.f"
		sum = (doublereal) (*n);
#line 332 "dlantp.f"
		k = 2;
#line 333 "dlantp.f"
		i__1 = *n - 1;
#line 333 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 334 "dlantp.f"
		    i__2 = *n - j;
#line 334 "dlantp.f"
		    dlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 335 "dlantp.f"
		    k = k + *n - j + 1;
#line 336 "dlantp.f"
/* L300: */
#line 336 "dlantp.f"
		}
#line 337 "dlantp.f"
	    } else {
#line 338 "dlantp.f"
		scale = 0.;
#line 339 "dlantp.f"
		sum = 1.;
#line 340 "dlantp.f"
		k = 1;
#line 341 "dlantp.f"
		i__1 = *n;
#line 341 "dlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 342 "dlantp.f"
		    i__2 = *n - j + 1;
#line 342 "dlantp.f"
		    dlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 343 "dlantp.f"
		    k = k + *n - j + 1;
#line 344 "dlantp.f"
/* L310: */
#line 344 "dlantp.f"
		}
#line 345 "dlantp.f"
	    }
#line 346 "dlantp.f"
	}
#line 347 "dlantp.f"
	value = scale * sqrt(sum);
#line 348 "dlantp.f"
    }

#line 350 "dlantp.f"
    ret_val = value;
#line 351 "dlantp.f"
    return ret_val;

/*     End of DLANTP */

} /* dlantp_ */

