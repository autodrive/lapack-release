#line 1 "clantp.f"
/* clantp.f -- translated by f2c (version 20100827).
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

#line 1 "clantp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANTP( NORM, UPLO, DIAG, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               WORK( * ) */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANTP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > triangular matrix A, supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return CLANTP */
/* > \verbatim */
/* > */
/* >    CLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANTP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANTP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clantp_(char *norm, char *uplo, char *diag, integer *n, 
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
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    extern logical sisnan_(doublereal *);


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

#line 165 "clantp.f"
    /* Parameter adjustments */
#line 165 "clantp.f"
    --work;
#line 165 "clantp.f"
    --ap;
#line 165 "clantp.f"

#line 165 "clantp.f"
    /* Function Body */
#line 165 "clantp.f"
    if (*n == 0) {
#line 166 "clantp.f"
	value = 0.;
#line 167 "clantp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 171 "clantp.f"
	k = 1;
#line 172 "clantp.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 173 "clantp.f"
	    value = 1.;
#line 174 "clantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 175 "clantp.f"
		i__1 = *n;
#line 175 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 176 "clantp.f"
		    i__2 = k + j - 2;
#line 176 "clantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 177 "clantp.f"
			sum = z_abs(&ap[i__]);
#line 178 "clantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 178 "clantp.f"
			    value = sum;
#line 178 "clantp.f"
			}
#line 179 "clantp.f"
/* L10: */
#line 179 "clantp.f"
		    }
#line 180 "clantp.f"
		    k += j;
#line 181 "clantp.f"
/* L20: */
#line 181 "clantp.f"
		}
#line 182 "clantp.f"
	    } else {
#line 183 "clantp.f"
		i__1 = *n;
#line 183 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 184 "clantp.f"
		    i__2 = k + *n - j;
#line 184 "clantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 185 "clantp.f"
			sum = z_abs(&ap[i__]);
#line 186 "clantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 186 "clantp.f"
			    value = sum;
#line 186 "clantp.f"
			}
#line 187 "clantp.f"
/* L30: */
#line 187 "clantp.f"
		    }
#line 188 "clantp.f"
		    k = k + *n - j + 1;
#line 189 "clantp.f"
/* L40: */
#line 189 "clantp.f"
		}
#line 190 "clantp.f"
	    }
#line 191 "clantp.f"
	} else {
#line 192 "clantp.f"
	    value = 0.;
#line 193 "clantp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 194 "clantp.f"
		i__1 = *n;
#line 194 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 195 "clantp.f"
		    i__2 = k + j - 1;
#line 195 "clantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 196 "clantp.f"
			sum = z_abs(&ap[i__]);
#line 197 "clantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 197 "clantp.f"
			    value = sum;
#line 197 "clantp.f"
			}
#line 198 "clantp.f"
/* L50: */
#line 198 "clantp.f"
		    }
#line 199 "clantp.f"
		    k += j;
#line 200 "clantp.f"
/* L60: */
#line 200 "clantp.f"
		}
#line 201 "clantp.f"
	    } else {
#line 202 "clantp.f"
		i__1 = *n;
#line 202 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 203 "clantp.f"
		    i__2 = k + *n - j;
#line 203 "clantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 204 "clantp.f"
			sum = z_abs(&ap[i__]);
#line 205 "clantp.f"
			if (value < sum || sisnan_(&sum)) {
#line 205 "clantp.f"
			    value = sum;
#line 205 "clantp.f"
			}
#line 206 "clantp.f"
/* L70: */
#line 206 "clantp.f"
		    }
#line 207 "clantp.f"
		    k = k + *n - j + 1;
#line 208 "clantp.f"
/* L80: */
#line 208 "clantp.f"
		}
#line 209 "clantp.f"
	    }
#line 210 "clantp.f"
	}
#line 211 "clantp.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 215 "clantp.f"
	value = 0.;
#line 216 "clantp.f"
	k = 1;
#line 217 "clantp.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 218 "clantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 219 "clantp.f"
	    i__1 = *n;
#line 219 "clantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 220 "clantp.f"
		if (udiag) {
#line 221 "clantp.f"
		    sum = 1.;
#line 222 "clantp.f"
		    i__2 = k + j - 2;
#line 222 "clantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 223 "clantp.f"
			sum += z_abs(&ap[i__]);
#line 224 "clantp.f"
/* L90: */
#line 224 "clantp.f"
		    }
#line 225 "clantp.f"
		} else {
#line 226 "clantp.f"
		    sum = 0.;
#line 227 "clantp.f"
		    i__2 = k + j - 1;
#line 227 "clantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 228 "clantp.f"
			sum += z_abs(&ap[i__]);
#line 229 "clantp.f"
/* L100: */
#line 229 "clantp.f"
		    }
#line 230 "clantp.f"
		}
#line 231 "clantp.f"
		k += j;
#line 232 "clantp.f"
		if (value < sum || sisnan_(&sum)) {
#line 232 "clantp.f"
		    value = sum;
#line 232 "clantp.f"
		}
#line 233 "clantp.f"
/* L110: */
#line 233 "clantp.f"
	    }
#line 234 "clantp.f"
	} else {
#line 235 "clantp.f"
	    i__1 = *n;
#line 235 "clantp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 236 "clantp.f"
		if (udiag) {
#line 237 "clantp.f"
		    sum = 1.;
#line 238 "clantp.f"
		    i__2 = k + *n - j;
#line 238 "clantp.f"
		    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 239 "clantp.f"
			sum += z_abs(&ap[i__]);
#line 240 "clantp.f"
/* L120: */
#line 240 "clantp.f"
		    }
#line 241 "clantp.f"
		} else {
#line 242 "clantp.f"
		    sum = 0.;
#line 243 "clantp.f"
		    i__2 = k + *n - j;
#line 243 "clantp.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 244 "clantp.f"
			sum += z_abs(&ap[i__]);
#line 245 "clantp.f"
/* L130: */
#line 245 "clantp.f"
		    }
#line 246 "clantp.f"
		}
#line 247 "clantp.f"
		k = k + *n - j + 1;
#line 248 "clantp.f"
		if (value < sum || sisnan_(&sum)) {
#line 248 "clantp.f"
		    value = sum;
#line 248 "clantp.f"
		}
#line 249 "clantp.f"
/* L140: */
#line 249 "clantp.f"
	    }
#line 250 "clantp.f"
	}
#line 251 "clantp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 255 "clantp.f"
	k = 1;
#line 256 "clantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 257 "clantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 258 "clantp.f"
		i__1 = *n;
#line 258 "clantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 259 "clantp.f"
		    work[i__] = 1.;
#line 260 "clantp.f"
/* L150: */
#line 260 "clantp.f"
		}
#line 261 "clantp.f"
		i__1 = *n;
#line 261 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 262 "clantp.f"
		    i__2 = j - 1;
#line 262 "clantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 263 "clantp.f"
			work[i__] += z_abs(&ap[k]);
#line 264 "clantp.f"
			++k;
#line 265 "clantp.f"
/* L160: */
#line 265 "clantp.f"
		    }
#line 266 "clantp.f"
		    ++k;
#line 267 "clantp.f"
/* L170: */
#line 267 "clantp.f"
		}
#line 268 "clantp.f"
	    } else {
#line 269 "clantp.f"
		i__1 = *n;
#line 269 "clantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 270 "clantp.f"
		    work[i__] = 0.;
#line 271 "clantp.f"
/* L180: */
#line 271 "clantp.f"
		}
#line 272 "clantp.f"
		i__1 = *n;
#line 272 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 273 "clantp.f"
		    i__2 = j;
#line 273 "clantp.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 274 "clantp.f"
			work[i__] += z_abs(&ap[k]);
#line 275 "clantp.f"
			++k;
#line 276 "clantp.f"
/* L190: */
#line 276 "clantp.f"
		    }
#line 277 "clantp.f"
/* L200: */
#line 277 "clantp.f"
		}
#line 278 "clantp.f"
	    }
#line 279 "clantp.f"
	} else {
#line 280 "clantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 281 "clantp.f"
		i__1 = *n;
#line 281 "clantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 282 "clantp.f"
		    work[i__] = 1.;
#line 283 "clantp.f"
/* L210: */
#line 283 "clantp.f"
		}
#line 284 "clantp.f"
		i__1 = *n;
#line 284 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 285 "clantp.f"
		    ++k;
#line 286 "clantp.f"
		    i__2 = *n;
#line 286 "clantp.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 287 "clantp.f"
			work[i__] += z_abs(&ap[k]);
#line 288 "clantp.f"
			++k;
#line 289 "clantp.f"
/* L220: */
#line 289 "clantp.f"
		    }
#line 290 "clantp.f"
/* L230: */
#line 290 "clantp.f"
		}
#line 291 "clantp.f"
	    } else {
#line 292 "clantp.f"
		i__1 = *n;
#line 292 "clantp.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 293 "clantp.f"
		    work[i__] = 0.;
#line 294 "clantp.f"
/* L240: */
#line 294 "clantp.f"
		}
#line 295 "clantp.f"
		i__1 = *n;
#line 295 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "clantp.f"
		    i__2 = *n;
#line 296 "clantp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 297 "clantp.f"
			work[i__] += z_abs(&ap[k]);
#line 298 "clantp.f"
			++k;
#line 299 "clantp.f"
/* L250: */
#line 299 "clantp.f"
		    }
#line 300 "clantp.f"
/* L260: */
#line 300 "clantp.f"
		}
#line 301 "clantp.f"
	    }
#line 302 "clantp.f"
	}
#line 303 "clantp.f"
	value = 0.;
#line 304 "clantp.f"
	i__1 = *n;
#line 304 "clantp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 305 "clantp.f"
	    sum = work[i__];
#line 306 "clantp.f"
	    if (value < sum || sisnan_(&sum)) {
#line 306 "clantp.f"
		value = sum;
#line 306 "clantp.f"
	    }
#line 307 "clantp.f"
/* L270: */
#line 307 "clantp.f"
	}
#line 308 "clantp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 312 "clantp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 313 "clantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 314 "clantp.f"
		scale = 1.;
#line 315 "clantp.f"
		sum = (doublereal) (*n);
#line 316 "clantp.f"
		k = 2;
#line 317 "clantp.f"
		i__1 = *n;
#line 317 "clantp.f"
		for (j = 2; j <= i__1; ++j) {
#line 318 "clantp.f"
		    i__2 = j - 1;
#line 318 "clantp.f"
		    classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 319 "clantp.f"
		    k += j;
#line 320 "clantp.f"
/* L280: */
#line 320 "clantp.f"
		}
#line 321 "clantp.f"
	    } else {
#line 322 "clantp.f"
		scale = 0.;
#line 323 "clantp.f"
		sum = 1.;
#line 324 "clantp.f"
		k = 1;
#line 325 "clantp.f"
		i__1 = *n;
#line 325 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "clantp.f"
		    classq_(&j, &ap[k], &c__1, &scale, &sum);
#line 327 "clantp.f"
		    k += j;
#line 328 "clantp.f"
/* L290: */
#line 328 "clantp.f"
		}
#line 329 "clantp.f"
	    }
#line 330 "clantp.f"
	} else {
#line 331 "clantp.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 332 "clantp.f"
		scale = 1.;
#line 333 "clantp.f"
		sum = (doublereal) (*n);
#line 334 "clantp.f"
		k = 2;
#line 335 "clantp.f"
		i__1 = *n - 1;
#line 335 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 336 "clantp.f"
		    i__2 = *n - j;
#line 336 "clantp.f"
		    classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 337 "clantp.f"
		    k = k + *n - j + 1;
#line 338 "clantp.f"
/* L300: */
#line 338 "clantp.f"
		}
#line 339 "clantp.f"
	    } else {
#line 340 "clantp.f"
		scale = 0.;
#line 341 "clantp.f"
		sum = 1.;
#line 342 "clantp.f"
		k = 1;
#line 343 "clantp.f"
		i__1 = *n;
#line 343 "clantp.f"
		for (j = 1; j <= i__1; ++j) {
#line 344 "clantp.f"
		    i__2 = *n - j + 1;
#line 344 "clantp.f"
		    classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 345 "clantp.f"
		    k = k + *n - j + 1;
#line 346 "clantp.f"
/* L310: */
#line 346 "clantp.f"
		}
#line 347 "clantp.f"
	    }
#line 348 "clantp.f"
	}
#line 349 "clantp.f"
	value = scale * sqrt(sum);
#line 350 "clantp.f"
    }

#line 352 "clantp.f"
    ret_val = value;
#line 353 "clantp.f"
    return ret_val;

/*     End of CLANTP */

} /* clantp_ */

