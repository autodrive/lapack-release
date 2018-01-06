#line 1 "zlantp.f"
/* zlantp.f -- translated by f2c (version 20100827).
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

#line 1 "zlantp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlantp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlantp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlantp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANTP( NORM, UPLO, DIAG, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( * ) */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANTP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > triangular matrix A, supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return ZLANTP */
/* > \verbatim */
/* > */
/* >    ZLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANTP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANTP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
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

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlantp_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *work, ftnlen norm_len, ftnlen uplo_len,
	 ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum, scale;
    static logical udiag;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 165 "zlantp.f"
    /* Parameter adjustments */
#line 165 "zlantp.f"
    --work;
#line 165 "zlantp.f"
    --ap;
#line 165 "zlantp.f"

#line 165 "zlantp.f"
    /* Function Body */
#line 165 "zlantp.f"
    if (*n == 0) {
#line 166 "zlantp.f"
	value = 0.;
#line 167 "zlantp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 171 "zlantp.f"
	k = 1;
#line 172 "zlantp.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 173 "zlantp.f"
	    value = 1.;
#line 174 "zlantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 175 "zlantp.f"
		i__1 = *n;
#line 175 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 176 "zlantp.f"
		    i__2 = k + j - 2;
#line 176 "zlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 177 "zlantp.f"
			sum = z_abs(&ap[i__]);
#line 178 "zlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 178 "zlantp.f"
			    value = sum;
#line 178 "zlantp.f"
			}
#line 179 "zlantp.f"
/* L10: */
#line 179 "zlantp.f"
		    }
#line 180 "zlantp.f"
		    k += j;
#line 181 "zlantp.f"
/* L20: */
#line 181 "zlantp.f"
		}
#line 182 "zlantp.f"
	    } else {
#line 183 "zlantp.f"
		i__1 = *n;
#line 183 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 184 "zlantp.f"
		    i__2 = k + *n - j;
#line 184 "zlantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 185 "zlantp.f"
			sum = z_abs(&ap[i__]);
#line 186 "zlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 186 "zlantp.f"
			    value = sum;
#line 186 "zlantp.f"
			}
#line 187 "zlantp.f"
/* L30: */
#line 187 "zlantp.f"
		    }
#line 188 "zlantp.f"
		    k = k + *n - j + 1;
#line 189 "zlantp.f"
/* L40: */
#line 189 "zlantp.f"
		}
#line 190 "zlantp.f"
	    }
#line 191 "zlantp.f"
	} else {
#line 192 "zlantp.f"
	    value = 0.;
#line 193 "zlantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 194 "zlantp.f"
		i__1 = *n;
#line 194 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 195 "zlantp.f"
		    i__2 = k + j - 1;
#line 195 "zlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 196 "zlantp.f"
			sum = z_abs(&ap[i__]);
#line 197 "zlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 197 "zlantp.f"
			    value = sum;
#line 197 "zlantp.f"
			}
#line 198 "zlantp.f"
/* L50: */
#line 198 "zlantp.f"
		    }
#line 199 "zlantp.f"
		    k += j;
#line 200 "zlantp.f"
/* L60: */
#line 200 "zlantp.f"
		}
#line 201 "zlantp.f"
	    } else {
#line 202 "zlantp.f"
		i__1 = *n;
#line 202 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 203 "zlantp.f"
		    i__2 = k + *n - j;
#line 203 "zlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 204 "zlantp.f"
			sum = z_abs(&ap[i__]);
#line 205 "zlantp.f"
			if (value < sum || disnan_(&sum)) {
#line 205 "zlantp.f"
			    value = sum;
#line 205 "zlantp.f"
			}
#line 206 "zlantp.f"
/* L70: */
#line 206 "zlantp.f"
		    }
#line 207 "zlantp.f"
		    k = k + *n - j + 1;
#line 208 "zlantp.f"
/* L80: */
#line 208 "zlantp.f"
		}
#line 209 "zlantp.f"
	    }
#line 210 "zlantp.f"
	}
#line 211 "zlantp.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 215 "zlantp.f"
	value = 0.;
#line 216 "zlantp.f"
	k = 1;
#line 217 "zlantp.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 218 "zlantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 219 "zlantp.f"
	    i__1 = *n;
#line 219 "zlantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 220 "zlantp.f"
		if (udiag) {
#line 221 "zlantp.f"
		    sum = 1.;
#line 222 "zlantp.f"
		    i__2 = k + j - 2;
#line 222 "zlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 223 "zlantp.f"
			sum += z_abs(&ap[i__]);
#line 224 "zlantp.f"
/* L90: */
#line 224 "zlantp.f"
		    }
#line 225 "zlantp.f"
		} else {
#line 226 "zlantp.f"
		    sum = 0.;
#line 227 "zlantp.f"
		    i__2 = k + j - 1;
#line 227 "zlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 228 "zlantp.f"
			sum += z_abs(&ap[i__]);
#line 229 "zlantp.f"
/* L100: */
#line 229 "zlantp.f"
		    }
#line 230 "zlantp.f"
		}
#line 231 "zlantp.f"
		k += j;
#line 232 "zlantp.f"
		if (value < sum || disnan_(&sum)) {
#line 232 "zlantp.f"
		    value = sum;
#line 232 "zlantp.f"
		}
#line 233 "zlantp.f"
/* L110: */
#line 233 "zlantp.f"
	    }
#line 234 "zlantp.f"
	} else {
#line 235 "zlantp.f"
	    i__1 = *n;
#line 235 "zlantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 236 "zlantp.f"
		if (udiag) {
#line 237 "zlantp.f"
		    sum = 1.;
#line 238 "zlantp.f"
		    i__2 = k + *n - j;
#line 238 "zlantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 239 "zlantp.f"
			sum += z_abs(&ap[i__]);
#line 240 "zlantp.f"
/* L120: */
#line 240 "zlantp.f"
		    }
#line 241 "zlantp.f"
		} else {
#line 242 "zlantp.f"
		    sum = 0.;
#line 243 "zlantp.f"
		    i__2 = k + *n - j;
#line 243 "zlantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 244 "zlantp.f"
			sum += z_abs(&ap[i__]);
#line 245 "zlantp.f"
/* L130: */
#line 245 "zlantp.f"
		    }
#line 246 "zlantp.f"
		}
#line 247 "zlantp.f"
		k = k + *n - j + 1;
#line 248 "zlantp.f"
		if (value < sum || disnan_(&sum)) {
#line 248 "zlantp.f"
		    value = sum;
#line 248 "zlantp.f"
		}
#line 249 "zlantp.f"
/* L140: */
#line 249 "zlantp.f"
	    }
#line 250 "zlantp.f"
	}
#line 251 "zlantp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 255 "zlantp.f"
	k = 1;
#line 256 "zlantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 257 "zlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 258 "zlantp.f"
		i__1 = *n;
#line 258 "zlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 259 "zlantp.f"
		    work[i__] = 1.;
#line 260 "zlantp.f"
/* L150: */
#line 260 "zlantp.f"
		}
#line 261 "zlantp.f"
		i__1 = *n;
#line 261 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 262 "zlantp.f"
		    i__2 = j - 1;
#line 262 "zlantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 263 "zlantp.f"
			work[i__] += z_abs(&ap[k]);
#line 264 "zlantp.f"
			++k;
#line 265 "zlantp.f"
/* L160: */
#line 265 "zlantp.f"
		    }
#line 266 "zlantp.f"
		    ++k;
#line 267 "zlantp.f"
/* L170: */
#line 267 "zlantp.f"
		}
#line 268 "zlantp.f"
	    } else {
#line 269 "zlantp.f"
		i__1 = *n;
#line 269 "zlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 270 "zlantp.f"
		    work[i__] = 0.;
#line 271 "zlantp.f"
/* L180: */
#line 271 "zlantp.f"
		}
#line 272 "zlantp.f"
		i__1 = *n;
#line 272 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 273 "zlantp.f"
		    i__2 = j;
#line 273 "zlantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 274 "zlantp.f"
			work[i__] += z_abs(&ap[k]);
#line 275 "zlantp.f"
			++k;
#line 276 "zlantp.f"
/* L190: */
#line 276 "zlantp.f"
		    }
#line 277 "zlantp.f"
/* L200: */
#line 277 "zlantp.f"
		}
#line 278 "zlantp.f"
	    }
#line 279 "zlantp.f"
	} else {
#line 280 "zlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 281 "zlantp.f"
		i__1 = *n;
#line 281 "zlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 282 "zlantp.f"
		    work[i__] = 1.;
#line 283 "zlantp.f"
/* L210: */
#line 283 "zlantp.f"
		}
#line 284 "zlantp.f"
		i__1 = *n;
#line 284 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 285 "zlantp.f"
		    ++k;
#line 286 "zlantp.f"
		    i__2 = *n;
#line 286 "zlantp.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 287 "zlantp.f"
			work[i__] += z_abs(&ap[k]);
#line 288 "zlantp.f"
			++k;
#line 289 "zlantp.f"
/* L220: */
#line 289 "zlantp.f"
		    }
#line 290 "zlantp.f"
/* L230: */
#line 290 "zlantp.f"
		}
#line 291 "zlantp.f"
	    } else {
#line 292 "zlantp.f"
		i__1 = *n;
#line 292 "zlantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 293 "zlantp.f"
		    work[i__] = 0.;
#line 294 "zlantp.f"
/* L240: */
#line 294 "zlantp.f"
		}
#line 295 "zlantp.f"
		i__1 = *n;
#line 295 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "zlantp.f"
		    i__2 = *n;
#line 296 "zlantp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 297 "zlantp.f"
			work[i__] += z_abs(&ap[k]);
#line 298 "zlantp.f"
			++k;
#line 299 "zlantp.f"
/* L250: */
#line 299 "zlantp.f"
		    }
#line 300 "zlantp.f"
/* L260: */
#line 300 "zlantp.f"
		}
#line 301 "zlantp.f"
	    }
#line 302 "zlantp.f"
	}
#line 303 "zlantp.f"
	value = 0.;
#line 304 "zlantp.f"
	i__1 = *n;
#line 304 "zlantp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 305 "zlantp.f"
	    sum = work[i__];
#line 306 "zlantp.f"
	    if (value < sum || disnan_(&sum)) {
#line 306 "zlantp.f"
		value = sum;
#line 306 "zlantp.f"
	    }
#line 307 "zlantp.f"
/* L270: */
#line 307 "zlantp.f"
	}
#line 308 "zlantp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 312 "zlantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 313 "zlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 314 "zlantp.f"
		scale = 1.;
#line 315 "zlantp.f"
		sum = (doublereal) (*n);
#line 316 "zlantp.f"
		k = 2;
#line 317 "zlantp.f"
		i__1 = *n;
#line 317 "zlantp.f"
		for (j = 2; j <= i__1; ++j) {
#line 318 "zlantp.f"
		    i__2 = j - 1;
#line 318 "zlantp.f"
		    zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 319 "zlantp.f"
		    k += j;
#line 320 "zlantp.f"
/* L280: */
#line 320 "zlantp.f"
		}
#line 321 "zlantp.f"
	    } else {
#line 322 "zlantp.f"
		scale = 0.;
#line 323 "zlantp.f"
		sum = 1.;
#line 324 "zlantp.f"
		k = 1;
#line 325 "zlantp.f"
		i__1 = *n;
#line 325 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "zlantp.f"
		    zlassq_(&j, &ap[k], &c__1, &scale, &sum);
#line 327 "zlantp.f"
		    k += j;
#line 328 "zlantp.f"
/* L290: */
#line 328 "zlantp.f"
		}
#line 329 "zlantp.f"
	    }
#line 330 "zlantp.f"
	} else {
#line 331 "zlantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 332 "zlantp.f"
		scale = 1.;
#line 333 "zlantp.f"
		sum = (doublereal) (*n);
#line 334 "zlantp.f"
		k = 2;
#line 335 "zlantp.f"
		i__1 = *n - 1;
#line 335 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 336 "zlantp.f"
		    i__2 = *n - j;
#line 336 "zlantp.f"
		    zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 337 "zlantp.f"
		    k = k + *n - j + 1;
#line 338 "zlantp.f"
/* L300: */
#line 338 "zlantp.f"
		}
#line 339 "zlantp.f"
	    } else {
#line 340 "zlantp.f"
		scale = 0.;
#line 341 "zlantp.f"
		sum = 1.;
#line 342 "zlantp.f"
		k = 1;
#line 343 "zlantp.f"
		i__1 = *n;
#line 343 "zlantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 344 "zlantp.f"
		    i__2 = *n - j + 1;
#line 344 "zlantp.f"
		    zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 345 "zlantp.f"
		    k = k + *n - j + 1;
#line 346 "zlantp.f"
/* L310: */
#line 346 "zlantp.f"
		}
#line 347 "zlantp.f"
	    }
#line 348 "zlantp.f"
	}
#line 349 "zlantp.f"
	value = scale * sqrt(sum);
#line 350 "zlantp.f"
    }

#line 352 "zlantp.f"
    ret_val = value;
#line 353 "zlantp.f"
    return ret_val;

/*     End of ZLANTP */

} /* zlantp_ */

