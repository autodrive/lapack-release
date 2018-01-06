#line 1 "zlantb.f"
/* zlantb.f -- translated by f2c (version 20100827).
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

#line 1 "zlantb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANTB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlantb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlantb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlantb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANTB( NORM, UPLO, DIAG, N, K, AB, */
/*                        LDAB, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            K, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( * ) */
/*       COMPLEX*16         AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANTB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n triangular band matrix A,  with ( k + 1 ) diagonals. */
/* > \endverbatim */
/* > */
/* > \return ZLANTB */
/* > \verbatim */
/* > */
/* >    ZLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANTB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANTB is */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
	 doublecomplex *ab, integer *ldab, doublereal *work, ftnlen norm_len, 
	ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l;
    static doublereal sum, scale;
    static logical udiag;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 181 "zlantb.f"
    /* Parameter adjustments */
#line 181 "zlantb.f"
    ab_dim1 = *ldab;
#line 181 "zlantb.f"
    ab_offset = 1 + ab_dim1;
#line 181 "zlantb.f"
    ab -= ab_offset;
#line 181 "zlantb.f"
    --work;
#line 181 "zlantb.f"

#line 181 "zlantb.f"
    /* Function Body */
#line 181 "zlantb.f"
    if (*n == 0) {
#line 182 "zlantb.f"
	value = 0.;
#line 183 "zlantb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 187 "zlantb.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "zlantb.f"
	    value = 1.;
#line 189 "zlantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 190 "zlantb.f"
		i__1 = *n;
#line 190 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 191 "zlantb.f"
		    i__2 = *k + 2 - j;
#line 191 "zlantb.f"
		    i__3 = *k;
#line 191 "zlantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 192 "zlantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 193 "zlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 193 "zlantb.f"
			    value = sum;
#line 193 "zlantb.f"
			}
#line 194 "zlantb.f"
/* L10: */
#line 194 "zlantb.f"
		    }
#line 195 "zlantb.f"
/* L20: */
#line 195 "zlantb.f"
		}
#line 196 "zlantb.f"
	    } else {
#line 197 "zlantb.f"
		i__1 = *n;
#line 197 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 198 "zlantb.f"
		    i__2 = *n + 1 - j, i__4 = *k + 1;
#line 198 "zlantb.f"
		    i__3 = min(i__2,i__4);
#line 198 "zlantb.f"
		    for (i__ = 2; i__ <= i__3; ++i__) {
#line 199 "zlantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 200 "zlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 200 "zlantb.f"
			    value = sum;
#line 200 "zlantb.f"
			}
#line 201 "zlantb.f"
/* L30: */
#line 201 "zlantb.f"
		    }
#line 202 "zlantb.f"
/* L40: */
#line 202 "zlantb.f"
		}
#line 203 "zlantb.f"
	    }
#line 204 "zlantb.f"
	} else {
#line 205 "zlantb.f"
	    value = 0.;
#line 206 "zlantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 207 "zlantb.f"
		i__1 = *n;
#line 207 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 208 "zlantb.f"
		    i__3 = *k + 2 - j;
#line 208 "zlantb.f"
		    i__2 = *k + 1;
#line 208 "zlantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 209 "zlantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 210 "zlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 210 "zlantb.f"
			    value = sum;
#line 210 "zlantb.f"
			}
#line 211 "zlantb.f"
/* L50: */
#line 211 "zlantb.f"
		    }
#line 212 "zlantb.f"
/* L60: */
#line 212 "zlantb.f"
		}
#line 213 "zlantb.f"
	    } else {
#line 214 "zlantb.f"
		i__1 = *n;
#line 214 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 215 "zlantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 215 "zlantb.f"
		    i__2 = min(i__3,i__4);
#line 215 "zlantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 216 "zlantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 217 "zlantb.f"
			if (value < sum || disnan_(&sum)) {
#line 217 "zlantb.f"
			    value = sum;
#line 217 "zlantb.f"
			}
#line 218 "zlantb.f"
/* L70: */
#line 218 "zlantb.f"
		    }
#line 219 "zlantb.f"
/* L80: */
#line 219 "zlantb.f"
		}
#line 220 "zlantb.f"
	    }
#line 221 "zlantb.f"
	}
#line 222 "zlantb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 226 "zlantb.f"
	value = 0.;
#line 227 "zlantb.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 228 "zlantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 229 "zlantb.f"
	    i__1 = *n;
#line 229 "zlantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 230 "zlantb.f"
		if (udiag) {
#line 231 "zlantb.f"
		    sum = 1.;
/* Computing MAX */
#line 232 "zlantb.f"
		    i__2 = *k + 2 - j;
#line 232 "zlantb.f"
		    i__3 = *k;
#line 232 "zlantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 233 "zlantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 234 "zlantb.f"
/* L90: */
#line 234 "zlantb.f"
		    }
#line 235 "zlantb.f"
		} else {
#line 236 "zlantb.f"
		    sum = 0.;
/* Computing MAX */
#line 237 "zlantb.f"
		    i__3 = *k + 2 - j;
#line 237 "zlantb.f"
		    i__2 = *k + 1;
#line 237 "zlantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 238 "zlantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 239 "zlantb.f"
/* L100: */
#line 239 "zlantb.f"
		    }
#line 240 "zlantb.f"
		}
#line 241 "zlantb.f"
		if (value < sum || disnan_(&sum)) {
#line 241 "zlantb.f"
		    value = sum;
#line 241 "zlantb.f"
		}
#line 242 "zlantb.f"
/* L110: */
#line 242 "zlantb.f"
	    }
#line 243 "zlantb.f"
	} else {
#line 244 "zlantb.f"
	    i__1 = *n;
#line 244 "zlantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "zlantb.f"
		if (udiag) {
#line 246 "zlantb.f"
		    sum = 1.;
/* Computing MIN */
#line 247 "zlantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 247 "zlantb.f"
		    i__2 = min(i__3,i__4);
#line 247 "zlantb.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 248 "zlantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 249 "zlantb.f"
/* L120: */
#line 249 "zlantb.f"
		    }
#line 250 "zlantb.f"
		} else {
#line 251 "zlantb.f"
		    sum = 0.;
/* Computing MIN */
#line 252 "zlantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 252 "zlantb.f"
		    i__2 = min(i__3,i__4);
#line 252 "zlantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 253 "zlantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 254 "zlantb.f"
/* L130: */
#line 254 "zlantb.f"
		    }
#line 255 "zlantb.f"
		}
#line 256 "zlantb.f"
		if (value < sum || disnan_(&sum)) {
#line 256 "zlantb.f"
		    value = sum;
#line 256 "zlantb.f"
		}
#line 257 "zlantb.f"
/* L140: */
#line 257 "zlantb.f"
	    }
#line 258 "zlantb.f"
	}
#line 259 "zlantb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 263 "zlantb.f"
	value = 0.;
#line 264 "zlantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 265 "zlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 266 "zlantb.f"
		i__1 = *n;
#line 266 "zlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "zlantb.f"
		    work[i__] = 1.;
#line 268 "zlantb.f"
/* L150: */
#line 268 "zlantb.f"
		}
#line 269 "zlantb.f"
		i__1 = *n;
#line 269 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 270 "zlantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 271 "zlantb.f"
		    i__2 = 1, i__3 = j - *k;
#line 271 "zlantb.f"
		    i__4 = j - 1;
#line 271 "zlantb.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 272 "zlantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 273 "zlantb.f"
/* L160: */
#line 273 "zlantb.f"
		    }
#line 274 "zlantb.f"
/* L170: */
#line 274 "zlantb.f"
		}
#line 275 "zlantb.f"
	    } else {
#line 276 "zlantb.f"
		i__1 = *n;
#line 276 "zlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 277 "zlantb.f"
		    work[i__] = 0.;
#line 278 "zlantb.f"
/* L180: */
#line 278 "zlantb.f"
		}
#line 279 "zlantb.f"
		i__1 = *n;
#line 279 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 280 "zlantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 281 "zlantb.f"
		    i__4 = 1, i__2 = j - *k;
#line 281 "zlantb.f"
		    i__3 = j;
#line 281 "zlantb.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 282 "zlantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 283 "zlantb.f"
/* L190: */
#line 283 "zlantb.f"
		    }
#line 284 "zlantb.f"
/* L200: */
#line 284 "zlantb.f"
		}
#line 285 "zlantb.f"
	    }
#line 286 "zlantb.f"
	} else {
#line 287 "zlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 288 "zlantb.f"
		i__1 = *n;
#line 288 "zlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "zlantb.f"
		    work[i__] = 1.;
#line 290 "zlantb.f"
/* L210: */
#line 290 "zlantb.f"
		}
#line 291 "zlantb.f"
		i__1 = *n;
#line 291 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "zlantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 293 "zlantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 293 "zlantb.f"
		    i__3 = min(i__4,i__2);
#line 293 "zlantb.f"
		    for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 294 "zlantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 295 "zlantb.f"
/* L220: */
#line 295 "zlantb.f"
		    }
#line 296 "zlantb.f"
/* L230: */
#line 296 "zlantb.f"
		}
#line 297 "zlantb.f"
	    } else {
#line 298 "zlantb.f"
		i__1 = *n;
#line 298 "zlantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "zlantb.f"
		    work[i__] = 0.;
#line 300 "zlantb.f"
/* L240: */
#line 300 "zlantb.f"
		}
#line 301 "zlantb.f"
		i__1 = *n;
#line 301 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 302 "zlantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 303 "zlantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 303 "zlantb.f"
		    i__3 = min(i__4,i__2);
#line 303 "zlantb.f"
		    for (i__ = j; i__ <= i__3; ++i__) {
#line 304 "zlantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 305 "zlantb.f"
/* L250: */
#line 305 "zlantb.f"
		    }
#line 306 "zlantb.f"
/* L260: */
#line 306 "zlantb.f"
		}
#line 307 "zlantb.f"
	    }
#line 308 "zlantb.f"
	}
#line 309 "zlantb.f"
	i__1 = *n;
#line 309 "zlantb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 310 "zlantb.f"
	    sum = work[i__];
#line 311 "zlantb.f"
	    if (value < sum || disnan_(&sum)) {
#line 311 "zlantb.f"
		value = sum;
#line 311 "zlantb.f"
	    }
#line 312 "zlantb.f"
/* L270: */
#line 312 "zlantb.f"
	}
#line 313 "zlantb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 317 "zlantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 318 "zlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 319 "zlantb.f"
		scale = 1.;
#line 320 "zlantb.f"
		sum = (doublereal) (*n);
#line 321 "zlantb.f"
		if (*k > 0) {
#line 322 "zlantb.f"
		    i__1 = *n;
#line 322 "zlantb.f"
		    for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 323 "zlantb.f"
			i__4 = j - 1;
#line 323 "zlantb.f"
			i__3 = min(i__4,*k);
/* Computing MAX */
#line 323 "zlantb.f"
			i__2 = *k + 2 - j;
#line 323 "zlantb.f"
			zlassq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, 
				&scale, &sum);
#line 326 "zlantb.f"
/* L280: */
#line 326 "zlantb.f"
		    }
#line 327 "zlantb.f"
		}
#line 328 "zlantb.f"
	    } else {
#line 329 "zlantb.f"
		scale = 0.;
#line 330 "zlantb.f"
		sum = 1.;
#line 331 "zlantb.f"
		i__1 = *n;
#line 331 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 332 "zlantb.f"
		    i__4 = j, i__2 = *k + 1;
#line 332 "zlantb.f"
		    i__3 = min(i__4,i__2);
/* Computing MAX */
#line 332 "zlantb.f"
		    i__5 = *k + 2 - j;
#line 332 "zlantb.f"
		    zlassq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 334 "zlantb.f"
/* L290: */
#line 334 "zlantb.f"
		}
#line 335 "zlantb.f"
	    }
#line 336 "zlantb.f"
	} else {
#line 337 "zlantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 338 "zlantb.f"
		scale = 1.;
#line 339 "zlantb.f"
		sum = (doublereal) (*n);
#line 340 "zlantb.f"
		if (*k > 0) {
#line 341 "zlantb.f"
		    i__1 = *n - 1;
#line 341 "zlantb.f"
		    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 342 "zlantb.f"
			i__4 = *n - j;
#line 342 "zlantb.f"
			i__3 = min(i__4,*k);
#line 342 "zlantb.f"
			zlassq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, &
				sum);
#line 344 "zlantb.f"
/* L300: */
#line 344 "zlantb.f"
		    }
#line 345 "zlantb.f"
		}
#line 346 "zlantb.f"
	    } else {
#line 347 "zlantb.f"
		scale = 0.;
#line 348 "zlantb.f"
		sum = 1.;
#line 349 "zlantb.f"
		i__1 = *n;
#line 349 "zlantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 350 "zlantb.f"
		    i__4 = *n - j + 1, i__2 = *k + 1;
#line 350 "zlantb.f"
		    i__3 = min(i__4,i__2);
#line 350 "zlantb.f"
		    zlassq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
#line 352 "zlantb.f"
/* L310: */
#line 352 "zlantb.f"
		}
#line 353 "zlantb.f"
	    }
#line 354 "zlantb.f"
	}
#line 355 "zlantb.f"
	value = scale * sqrt(sum);
#line 356 "zlantb.f"
    }

#line 358 "zlantb.f"
    ret_val = value;
#line 359 "zlantb.f"
    return ret_val;

/*     End of ZLANTB */

} /* zlantb_ */

