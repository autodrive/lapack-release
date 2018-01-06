#line 1 "dlantb.f"
/* dlantb.f -- translated by f2c (version 20100827).
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

#line 1 "dlantb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANTB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANTB( NORM, UPLO, DIAG, N, K, AB, */
/*                        LDAB, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            K, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANTB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n triangular band matrix A,  with ( k + 1 ) diagonals. */
/* > \endverbatim */
/* > */
/* > \return DLANTB */
/* > \verbatim */
/* > */
/* >    DLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANTB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANTB is */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
doublereal dlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
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
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 179 "dlantb.f"
    /* Parameter adjustments */
#line 179 "dlantb.f"
    ab_dim1 = *ldab;
#line 179 "dlantb.f"
    ab_offset = 1 + ab_dim1;
#line 179 "dlantb.f"
    ab -= ab_offset;
#line 179 "dlantb.f"
    --work;
#line 179 "dlantb.f"

#line 179 "dlantb.f"
    /* Function Body */
#line 179 "dlantb.f"
    if (*n == 0) {
#line 180 "dlantb.f"
	value = 0.;
#line 181 "dlantb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 185 "dlantb.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 186 "dlantb.f"
	    value = 1.;
#line 187 "dlantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "dlantb.f"
		i__1 = *n;
#line 188 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 189 "dlantb.f"
		    i__2 = *k + 2 - j;
#line 189 "dlantb.f"
		    i__3 = *k;
#line 189 "dlantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 190 "dlantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 191 "dlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 191 "dlantb.f"
			    value = sum;
#line 191 "dlantb.f"
			}
#line 192 "dlantb.f"
/* L10: */
#line 192 "dlantb.f"
		    }
#line 193 "dlantb.f"
/* L20: */
#line 193 "dlantb.f"
		}
#line 194 "dlantb.f"
	    } else {
#line 195 "dlantb.f"
		i__1 = *n;
#line 195 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 196 "dlantb.f"
		    i__2 = *n + 1 - j, i__4 = *k + 1;
#line 196 "dlantb.f"
		    i__3 = min(i__2,i__4);
#line 196 "dlantb.f"
		    for (i__ = 2; i__ <= i__3; ++i__) {
#line 197 "dlantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 198 "dlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 198 "dlantb.f"
			    value = sum;
#line 198 "dlantb.f"
			}
#line 199 "dlantb.f"
/* L30: */
#line 199 "dlantb.f"
		    }
#line 200 "dlantb.f"
/* L40: */
#line 200 "dlantb.f"
		}
#line 201 "dlantb.f"
	    }
#line 202 "dlantb.f"
	} else {
#line 203 "dlantb.f"
	    value = 0.;
#line 204 "dlantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 205 "dlantb.f"
		i__1 = *n;
#line 205 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 206 "dlantb.f"
		    i__3 = *k + 2 - j;
#line 206 "dlantb.f"
		    i__2 = *k + 1;
#line 206 "dlantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 207 "dlantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 208 "dlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 208 "dlantb.f"
			    value = sum;
#line 208 "dlantb.f"
			}
#line 209 "dlantb.f"
/* L50: */
#line 209 "dlantb.f"
		    }
#line 210 "dlantb.f"
/* L60: */
#line 210 "dlantb.f"
		}
#line 211 "dlantb.f"
	    } else {
#line 212 "dlantb.f"
		i__1 = *n;
#line 212 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 213 "dlantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 213 "dlantb.f"
		    i__2 = min(i__3,i__4);
#line 213 "dlantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 214 "dlantb.f"
			sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 215 "dlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 215 "dlantb.f"
			    value = sum;
#line 215 "dlantb.f"
			}
#line 216 "dlantb.f"
/* L70: */
#line 216 "dlantb.f"
		    }
#line 217 "dlantb.f"
/* L80: */
#line 217 "dlantb.f"
		}
#line 218 "dlantb.f"
	    }
#line 219 "dlantb.f"
	}
#line 220 "dlantb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 224 "dlantb.f"
	value = 0.;
#line 225 "dlantb.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 226 "dlantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 227 "dlantb.f"
	    i__1 = *n;
#line 227 "dlantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 228 "dlantb.f"
		if (udiag) {
#line 229 "dlantb.f"
		    sum = 1.;
/* Computing MAX */
#line 230 "dlantb.f"
		    i__2 = *k + 2 - j;
#line 230 "dlantb.f"
		    i__3 = *k;
#line 230 "dlantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 231 "dlantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 232 "dlantb.f"
/* L90: */
#line 232 "dlantb.f"
		    }
#line 233 "dlantb.f"
		} else {
#line 234 "dlantb.f"
		    sum = 0.;
/* Computing MAX */
#line 235 "dlantb.f"
		    i__3 = *k + 2 - j;
#line 235 "dlantb.f"
		    i__2 = *k + 1;
#line 235 "dlantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 236 "dlantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 237 "dlantb.f"
/* L100: */
#line 237 "dlantb.f"
		    }
#line 238 "dlantb.f"
		}
#line 239 "dlantb.f"
		if (value < sum || disnan_(&sum)) {
#line 239 "dlantb.f"
		    value = sum;
#line 239 "dlantb.f"
		}
#line 240 "dlantb.f"
/* L110: */
#line 240 "dlantb.f"
	    }
#line 241 "dlantb.f"
	} else {
#line 242 "dlantb.f"
	    i__1 = *n;
#line 242 "dlantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 243 "dlantb.f"
		if (udiag) {
#line 244 "dlantb.f"
		    sum = 1.;
/* Computing MIN */
#line 245 "dlantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 245 "dlantb.f"
		    i__2 = min(i__3,i__4);
#line 245 "dlantb.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 246 "dlantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 247 "dlantb.f"
/* L120: */
#line 247 "dlantb.f"
		    }
#line 248 "dlantb.f"
		} else {
#line 249 "dlantb.f"
		    sum = 0.;
/* Computing MIN */
#line 250 "dlantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 250 "dlantb.f"
		    i__2 = min(i__3,i__4);
#line 250 "dlantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 251 "dlantb.f"
			sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 252 "dlantb.f"
/* L130: */
#line 252 "dlantb.f"
		    }
#line 253 "dlantb.f"
		}
#line 254 "dlantb.f"
		if (value < sum || disnan_(&sum)) {
#line 254 "dlantb.f"
		    value = sum;
#line 254 "dlantb.f"
		}
#line 255 "dlantb.f"
/* L140: */
#line 255 "dlantb.f"
	    }
#line 256 "dlantb.f"
	}
#line 257 "dlantb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 261 "dlantb.f"
	value = 0.;
#line 262 "dlantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 263 "dlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 264 "dlantb.f"
		i__1 = *n;
#line 264 "dlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 265 "dlantb.f"
		    work[i__] = 1.;
#line 266 "dlantb.f"
/* L150: */
#line 266 "dlantb.f"
		}
#line 267 "dlantb.f"
		i__1 = *n;
#line 267 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "dlantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 269 "dlantb.f"
		    i__2 = 1, i__3 = j - *k;
#line 269 "dlantb.f"
		    i__4 = j - 1;
#line 269 "dlantb.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 270 "dlantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 271 "dlantb.f"
/* L160: */
#line 271 "dlantb.f"
		    }
#line 272 "dlantb.f"
/* L170: */
#line 272 "dlantb.f"
		}
#line 273 "dlantb.f"
	    } else {
#line 274 "dlantb.f"
		i__1 = *n;
#line 274 "dlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "dlantb.f"
		    work[i__] = 0.;
#line 276 "dlantb.f"
/* L180: */
#line 276 "dlantb.f"
		}
#line 277 "dlantb.f"
		i__1 = *n;
#line 277 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 278 "dlantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 279 "dlantb.f"
		    i__4 = 1, i__2 = j - *k;
#line 279 "dlantb.f"
		    i__3 = j;
#line 279 "dlantb.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 280 "dlantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 281 "dlantb.f"
/* L190: */
#line 281 "dlantb.f"
		    }
#line 282 "dlantb.f"
/* L200: */
#line 282 "dlantb.f"
		}
#line 283 "dlantb.f"
	    }
#line 284 "dlantb.f"
	} else {
#line 285 "dlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 286 "dlantb.f"
		i__1 = *n;
#line 286 "dlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "dlantb.f"
		    work[i__] = 1.;
#line 288 "dlantb.f"
/* L210: */
#line 288 "dlantb.f"
		}
#line 289 "dlantb.f"
		i__1 = *n;
#line 289 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 290 "dlantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 291 "dlantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 291 "dlantb.f"
		    i__3 = min(i__4,i__2);
#line 291 "dlantb.f"
		    for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 292 "dlantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 293 "dlantb.f"
/* L220: */
#line 293 "dlantb.f"
		    }
#line 294 "dlantb.f"
/* L230: */
#line 294 "dlantb.f"
		}
#line 295 "dlantb.f"
	    } else {
#line 296 "dlantb.f"
		i__1 = *n;
#line 296 "dlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "dlantb.f"
		    work[i__] = 0.;
#line 298 "dlantb.f"
/* L240: */
#line 298 "dlantb.f"
		}
#line 299 "dlantb.f"
		i__1 = *n;
#line 299 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 300 "dlantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 301 "dlantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 301 "dlantb.f"
		    i__3 = min(i__4,i__2);
#line 301 "dlantb.f"
		    for (i__ = j; i__ <= i__3; ++i__) {
#line 302 "dlantb.f"
			work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
				d__1));
#line 303 "dlantb.f"
/* L250: */
#line 303 "dlantb.f"
		    }
#line 304 "dlantb.f"
/* L260: */
#line 304 "dlantb.f"
		}
#line 305 "dlantb.f"
	    }
#line 306 "dlantb.f"
	}
#line 307 "dlantb.f"
	i__1 = *n;
#line 307 "dlantb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "dlantb.f"
	    sum = work[i__];
#line 309 "dlantb.f"
	    if (value < sum || disnan_(&sum)) {
#line 309 "dlantb.f"
		value = sum;
#line 309 "dlantb.f"
	    }
#line 310 "dlantb.f"
/* L270: */
#line 310 "dlantb.f"
	}
#line 311 "dlantb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 315 "dlantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 316 "dlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 317 "dlantb.f"
		scale = 1.;
#line 318 "dlantb.f"
		sum = (doublereal) (*n);
#line 319 "dlantb.f"
		if (*k > 0) {
#line 320 "dlantb.f"
		    i__1 = *n;
#line 320 "dlantb.f"
		    for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 321 "dlantb.f"
			i__4 = j - 1;
#line 321 "dlantb.f"
			i__3 = min(i__4,*k);
/* Computing MAX */
#line 321 "dlantb.f"
			i__2 = *k + 2 - j;
#line 321 "dlantb.f"
			dlassq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, 
				&scale, &sum);
#line 324 "dlantb.f"
/* L280: */
#line 324 "dlantb.f"
		    }
#line 325 "dlantb.f"
		}
#line 326 "dlantb.f"
	    } else {
#line 327 "dlantb.f"
		scale = 0.;
#line 328 "dlantb.f"
		sum = 1.;
#line 329 "dlantb.f"
		i__1 = *n;
#line 329 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 330 "dlantb.f"
		    i__4 = j, i__2 = *k + 1;
#line 330 "dlantb.f"
		    i__3 = min(i__4,i__2);
/* Computing MAX */
#line 330 "dlantb.f"
		    i__5 = *k + 2 - j;
#line 330 "dlantb.f"
		    dlassq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 332 "dlantb.f"
/* L290: */
#line 332 "dlantb.f"
		}
#line 333 "dlantb.f"
	    }
#line 334 "dlantb.f"
	} else {
#line 335 "dlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 336 "dlantb.f"
		scale = 1.;
#line 337 "dlantb.f"
		sum = (doublereal) (*n);
#line 338 "dlantb.f"
		if (*k > 0) {
#line 339 "dlantb.f"
		    i__1 = *n - 1;
#line 339 "dlantb.f"
		    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 340 "dlantb.f"
			i__4 = *n - j;
#line 340 "dlantb.f"
			i__3 = min(i__4,*k);
#line 340 "dlantb.f"
			dlassq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, &
				sum);
#line 342 "dlantb.f"
/* L300: */
#line 342 "dlantb.f"
		    }
#line 343 "dlantb.f"
		}
#line 344 "dlantb.f"
	    } else {
#line 345 "dlantb.f"
		scale = 0.;
#line 346 "dlantb.f"
		sum = 1.;
#line 347 "dlantb.f"
		i__1 = *n;
#line 347 "dlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 348 "dlantb.f"
		    i__4 = *n - j + 1, i__2 = *k + 1;
#line 348 "dlantb.f"
		    i__3 = min(i__4,i__2);
#line 348 "dlantb.f"
		    dlassq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
#line 350 "dlantb.f"
/* L310: */
#line 350 "dlantb.f"
		}
#line 351 "dlantb.f"
	    }
#line 352 "dlantb.f"
	}
#line 353 "dlantb.f"
	value = scale * sqrt(sum);
#line 354 "dlantb.f"
    }

#line 356 "dlantb.f"
    ret_val = value;
#line 357 "dlantb.f"
    return ret_val;

/*     End of DLANTB */

} /* dlantb_ */

