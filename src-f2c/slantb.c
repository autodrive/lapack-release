#line 1 "slantb.f"
/* slantb.f -- translated by f2c (version 20100827).
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

#line 1 "slantb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANTB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANTB( NORM, UPLO, DIAG, N, K, AB, */
/*                        LDAB, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            K, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANTB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n triangular band matrix A,  with ( k + 1 ) diagonals. */
/* > \endverbatim */
/* > */
/* > \return SLANTB */
/* > \verbatim */
/* > */
/* >    SLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANTB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANTB is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of super-diagonals of the matrix A if UPLO = 'U', */
/* >          or the number of sub-diagonals of the matrix A if UPLO = 'L'. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          The upper or lower triangular band matrix A, stored in the */
/* >          first k+1 rows of AB.  The j-th column of A is stored */
/* >          in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k). */
/* >          Note that when DIAG = 'U', the elements of the array AB */
/* >          corresponding to the diagonal elements of the matrix A are */
/* >          not referenced, but are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= K+1. */
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
doublereal slantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
	 doublereal *ab, integer *ldab, doublereal *work, ftnlen norm_len, 
	ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l;
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

#line 179 "slantb.f"
    /* Parameter adjustments */
#line 179 "slantb.f"
    ab_dim1 = *ldab;
#line 179 "slantb.f"
    ab_offset = 1 + ab_dim1;
#line 179 "slantb.f"
    ab -= ab_offset;
#line 179 "slantb.f"
    --work;
#line 179 "slantb.f"

#line 179 "slantb.f"
    /* Function Body */
#line 179 "slantb.f"
    if (*n == 0) {
#line 180 "slantb.f"
	value = 0.;
#line 181 "slantb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 185 "slantb.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 186 "slantb.f"
	    value = 1.;
#line 187 "slantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "slantb.f"
		i__1 = *n;
#line 188 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 189 "slantb.f"
		    i__2 = *k + 2 - j;
#line 189 "slantb.f"
		    i__3 = *k;
#line 189 "slantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 190 "slantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 191 "slantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 191 "slantb.f"
			    value = sum;
#line 191 "slantb.f"
			}
#line 192 "slantb.f"
/* L10: */
#line 192 "slantb.f"
		    }
#line 193 "slantb.f"
/* L20: */
#line 193 "slantb.f"
		}
#line 194 "slantb.f"
	    } else {
#line 195 "slantb.f"
		i__1 = *n;
#line 195 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 196 "slantb.f"
		    i__2 = *n + 1 - j, i__4 = *k + 1;
#line 196 "slantb.f"
		    i__3 = min(i__2,i__4);
#line 196 "slantb.f"
		    for (i__ = 2; i__ <= i__3; ++i__) {
#line 197 "slantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 198 "slantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 198 "slantb.f"
			    value = sum;
#line 198 "slantb.f"
			}
#line 199 "slantb.f"
/* L30: */
#line 199 "slantb.f"
		    }
#line 200 "slantb.f"
/* L40: */
#line 200 "slantb.f"
		}
#line 201 "slantb.f"
	    }
#line 202 "slantb.f"
	} else {
#line 203 "slantb.f"
	    value = 0.;
#line 204 "slantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 205 "slantb.f"
		i__1 = *n;
#line 205 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 206 "slantb.f"
		    i__3 = *k + 2 - j;
#line 206 "slantb.f"
		    i__2 = *k + 1;
#line 206 "slantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 207 "slantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 208 "slantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 208 "slantb.f"
			    value = sum;
#line 208 "slantb.f"
			}
#line 209 "slantb.f"
/* L50: */
#line 209 "slantb.f"
		    }
#line 210 "slantb.f"
/* L60: */
#line 210 "slantb.f"
		}
#line 211 "slantb.f"
	    } else {
#line 212 "slantb.f"
		i__1 = *n;
#line 212 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 213 "slantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 213 "slantb.f"
		    i__2 = min(i__3,i__4);
#line 213 "slantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 214 "slantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 215 "slantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 215 "slantb.f"
			    value = sum;
#line 215 "slantb.f"
			}
#line 216 "slantb.f"
/* L70: */
#line 216 "slantb.f"
		    }
#line 217 "slantb.f"
/* L80: */
#line 217 "slantb.f"
		}
#line 218 "slantb.f"
	    }
#line 219 "slantb.f"
	}
#line 220 "slantb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 224 "slantb.f"
	value = 0.;
#line 225 "slantb.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 226 "slantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 227 "slantb.f"
	    i__1 = *n;
#line 227 "slantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 228 "slantb.f"
		if (udiag) {
#line 229 "slantb.f"
		    sum = 1.;
/* Computing MAX */
#line 230 "slantb.f"
		    i__2 = *k + 2 - j;
#line 230 "slantb.f"
		    i__3 = *k;
#line 230 "slantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 231 "slantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 232 "slantb.f"
/* L90: */
#line 232 "slantb.f"
		    }
#line 233 "slantb.f"
		} else {
#line 234 "slantb.f"
		    sum = 0.;
/* Computing MAX */
#line 235 "slantb.f"
		    i__3 = *k + 2 - j;
#line 235 "slantb.f"
		    i__2 = *k + 1;
#line 235 "slantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 236 "slantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 237 "slantb.f"
/* L100: */
#line 237 "slantb.f"
		    }
#line 238 "slantb.f"
		}
#line 239 "slantb.f"
		if (value < sum || sisnan_(&sum)) {
#line 239 "slantb.f"
		    value = sum;
#line 239 "slantb.f"
		}
#line 240 "slantb.f"
/* L110: */
#line 240 "slantb.f"
	    }
#line 241 "slantb.f"
	} else {
#line 242 "slantb.f"
	    i__1 = *n;
#line 242 "slantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 243 "slantb.f"
		if (udiag) {
#line 244 "slantb.f"
		    sum = 1.;
/* Computing MIN */
#line 245 "slantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 245 "slantb.f"
		    i__2 = min(i__3,i__4);
#line 245 "slantb.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 246 "slantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 247 "slantb.f"
/* L120: */
#line 247 "slantb.f"
		    }
#line 248 "slantb.f"
		} else {
#line 249 "slantb.f"
		    sum = 0.;
/* Computing MIN */
#line 250 "slantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 250 "slantb.f"
		    i__2 = min(i__3,i__4);
#line 250 "slantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 251 "slantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 252 "slantb.f"
/* L130: */
#line 252 "slantb.f"
		    }
#line 253 "slantb.f"
		}
#line 254 "slantb.f"
		if (value < sum || sisnan_(&sum)) {
#line 254 "slantb.f"
		    value = sum;
#line 254 "slantb.f"
		}
#line 255 "slantb.f"
/* L140: */
#line 255 "slantb.f"
	    }
#line 256 "slantb.f"
	}
#line 257 "slantb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 261 "slantb.f"
	value = 0.;
#line 262 "slantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 263 "slantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 264 "slantb.f"
		i__1 = *n;
#line 264 "slantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 265 "slantb.f"
		    work[i__] = 1.;
#line 266 "slantb.f"
/* L150: */
#line 266 "slantb.f"
		}
#line 267 "slantb.f"
		i__1 = *n;
#line 267 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "slantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 269 "slantb.f"
		    i__2 = 1, i__3 = j - *k;
#line 269 "slantb.f"
		    i__4 = j - 1;
#line 269 "slantb.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 270 "slantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 271 "slantb.f"
/* L160: */
#line 271 "slantb.f"
		    }
#line 272 "slantb.f"
/* L170: */
#line 272 "slantb.f"
		}
#line 273 "slantb.f"
	    } else {
#line 274 "slantb.f"
		i__1 = *n;
#line 274 "slantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "slantb.f"
		    work[i__] = 0.;
#line 276 "slantb.f"
/* L180: */
#line 276 "slantb.f"
		}
#line 277 "slantb.f"
		i__1 = *n;
#line 277 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 278 "slantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 279 "slantb.f"
		    i__4 = 1, i__2 = j - *k;
#line 279 "slantb.f"
		    i__3 = j;
#line 279 "slantb.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 280 "slantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 281 "slantb.f"
/* L190: */
#line 281 "slantb.f"
		    }
#line 282 "slantb.f"
/* L200: */
#line 282 "slantb.f"
		}
#line 283 "slantb.f"
	    }
#line 284 "slantb.f"
	} else {
#line 285 "slantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 286 "slantb.f"
		i__1 = *n;
#line 286 "slantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "slantb.f"
		    work[i__] = 1.;
#line 288 "slantb.f"
/* L210: */
#line 288 "slantb.f"
		}
#line 289 "slantb.f"
		i__1 = *n;
#line 289 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 290 "slantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 291 "slantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 291 "slantb.f"
		    i__3 = min(i__4,i__2);
#line 291 "slantb.f"
		    for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 292 "slantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 293 "slantb.f"
/* L220: */
#line 293 "slantb.f"
		    }
#line 294 "slantb.f"
/* L230: */
#line 294 "slantb.f"
		}
#line 295 "slantb.f"
	    } else {
#line 296 "slantb.f"
		i__1 = *n;
#line 296 "slantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "slantb.f"
		    work[i__] = 0.;
#line 298 "slantb.f"
/* L240: */
#line 298 "slantb.f"
		}
#line 299 "slantb.f"
		i__1 = *n;
#line 299 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 300 "slantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 301 "slantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 301 "slantb.f"
		    i__3 = min(i__4,i__2);
#line 301 "slantb.f"
		    for (i__ = j; i__ <= i__3; ++i__) {
#line 302 "slantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 303 "slantb.f"
/* L250: */
#line 303 "slantb.f"
		    }
#line 304 "slantb.f"
/* L260: */
#line 304 "slantb.f"
		}
#line 305 "slantb.f"
	    }
#line 306 "slantb.f"
	}
#line 307 "slantb.f"
	i__1 = *n;
#line 307 "slantb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "slantb.f"
	    sum = work[i__];
#line 309 "slantb.f"
	    if (value < sum || sisnan_(&sum)) {
#line 309 "slantb.f"
		value = sum;
#line 309 "slantb.f"
	    }
#line 310 "slantb.f"
/* L270: */
#line 310 "slantb.f"
	}
#line 311 "slantb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 315 "slantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 316 "slantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 317 "slantb.f"
		scale = 1.;
#line 318 "slantb.f"
		sum = (doublereal) (*n);
#line 319 "slantb.f"
		if (*k > 0) {
#line 320 "slantb.f"
		    i__1 = *n;
#line 320 "slantb.f"
		    for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 321 "slantb.f"
			i__4 = j - 1;
#line 321 "slantb.f"
			i__3 = min(i__4,*k);
/* Computing MAX */
#line 321 "slantb.f"
			i__2 = *k + 2 - j;
#line 321 "slantb.f"
			slassq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, 
				&scale, &sum);
#line 324 "slantb.f"
/* L280: */
#line 324 "slantb.f"
		    }
#line 325 "slantb.f"
		}
#line 326 "slantb.f"
	    } else {
#line 327 "slantb.f"
		scale = 0.;
#line 328 "slantb.f"
		sum = 1.;
#line 329 "slantb.f"
		i__1 = *n;
#line 329 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 330 "slantb.f"
		    i__4 = j, i__2 = *k + 1;
#line 330 "slantb.f"
		    i__3 = min(i__4,i__2);
/* Computing MAX */
#line 330 "slantb.f"
		    i__5 = *k + 2 - j;
#line 330 "slantb.f"
		    slassq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 332 "slantb.f"
/* L290: */
#line 332 "slantb.f"
		}
#line 333 "slantb.f"
	    }
#line 334 "slantb.f"
	} else {
#line 335 "slantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 336 "slantb.f"
		scale = 1.;
#line 337 "slantb.f"
		sum = (doublereal) (*n);
#line 338 "slantb.f"
		if (*k > 0) {
#line 339 "slantb.f"
		    i__1 = *n - 1;
#line 339 "slantb.f"
		    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 340 "slantb.f"
			i__4 = *n - j;
#line 340 "slantb.f"
			i__3 = min(i__4,*k);
#line 340 "slantb.f"
			slassq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, &
				sum);
#line 342 "slantb.f"
/* L300: */
#line 342 "slantb.f"
		    }
#line 343 "slantb.f"
		}
#line 344 "slantb.f"
	    } else {
#line 345 "slantb.f"
		scale = 0.;
#line 346 "slantb.f"
		sum = 1.;
#line 347 "slantb.f"
		i__1 = *n;
#line 347 "slantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 348 "slantb.f"
		    i__4 = *n - j + 1, i__2 = *k + 1;
#line 348 "slantb.f"
		    i__3 = min(i__4,i__2);
#line 348 "slantb.f"
		    slassq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
#line 350 "slantb.f"
/* L310: */
#line 350 "slantb.f"
		}
#line 351 "slantb.f"
	    }
#line 352 "slantb.f"
	}
#line 353 "slantb.f"
	value = scale * sqrt(sum);
#line 354 "slantb.f"
    }

#line 356 "slantb.f"
    ret_val = value;
#line 357 "slantb.f"
    return ret_val;

/*     End of SLANTB */

} /* slantb_ */

