#line 1 "clantb.f"
/* clantb.f -- translated by f2c (version 20100827).
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

#line 1 "clantb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a triangular band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANTB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANTB( NORM, UPLO, DIAG, N, K, AB, */
/*                        LDAB, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            K, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               WORK( * ) */
/*       COMPLEX            AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANTB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n triangular band matrix A,  with ( k + 1 ) diagonals. */
/* > \endverbatim */
/* > */
/* > \return CLANTB */
/* > \verbatim */
/* > */
/* >    CLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANTB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANTB is */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
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

#line 181 "clantb.f"
    /* Parameter adjustments */
#line 181 "clantb.f"
    ab_dim1 = *ldab;
#line 181 "clantb.f"
    ab_offset = 1 + ab_dim1;
#line 181 "clantb.f"
    ab -= ab_offset;
#line 181 "clantb.f"
    --work;
#line 181 "clantb.f"

#line 181 "clantb.f"
    /* Function Body */
#line 181 "clantb.f"
    if (*n == 0) {
#line 182 "clantb.f"
	value = 0.;
#line 183 "clantb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 187 "clantb.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "clantb.f"
	    value = 1.;
#line 189 "clantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 190 "clantb.f"
		i__1 = *n;
#line 190 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 191 "clantb.f"
		    i__2 = *k + 2 - j;
#line 191 "clantb.f"
		    i__3 = *k;
#line 191 "clantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 192 "clantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 193 "clantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 193 "clantb.f"
			    value = sum;
#line 193 "clantb.f"
			}
#line 194 "clantb.f"
/* L10: */
#line 194 "clantb.f"
		    }
#line 195 "clantb.f"
/* L20: */
#line 195 "clantb.f"
		}
#line 196 "clantb.f"
	    } else {
#line 197 "clantb.f"
		i__1 = *n;
#line 197 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 198 "clantb.f"
		    i__2 = *n + 1 - j, i__4 = *k + 1;
#line 198 "clantb.f"
		    i__3 = min(i__2,i__4);
#line 198 "clantb.f"
		    for (i__ = 2; i__ <= i__3; ++i__) {
#line 199 "clantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 200 "clantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 200 "clantb.f"
			    value = sum;
#line 200 "clantb.f"
			}
#line 201 "clantb.f"
/* L30: */
#line 201 "clantb.f"
		    }
#line 202 "clantb.f"
/* L40: */
#line 202 "clantb.f"
		}
#line 203 "clantb.f"
	    }
#line 204 "clantb.f"
	} else {
#line 205 "clantb.f"
	    value = 0.;
#line 206 "clantb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 207 "clantb.f"
		i__1 = *n;
#line 207 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 208 "clantb.f"
		    i__3 = *k + 2 - j;
#line 208 "clantb.f"
		    i__2 = *k + 1;
#line 208 "clantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 209 "clantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 210 "clantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 210 "clantb.f"
			    value = sum;
#line 210 "clantb.f"
			}
#line 211 "clantb.f"
/* L50: */
#line 211 "clantb.f"
		    }
#line 212 "clantb.f"
/* L60: */
#line 212 "clantb.f"
		}
#line 213 "clantb.f"
	    } else {
#line 214 "clantb.f"
		i__1 = *n;
#line 214 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 215 "clantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 215 "clantb.f"
		    i__2 = min(i__3,i__4);
#line 215 "clantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 216 "clantb.f"
			sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 217 "clantb.f"
			if (value < sum || sisnan_(&sum)) {
#line 217 "clantb.f"
			    value = sum;
#line 217 "clantb.f"
			}
#line 218 "clantb.f"
/* L70: */
#line 218 "clantb.f"
		    }
#line 219 "clantb.f"
/* L80: */
#line 219 "clantb.f"
		}
#line 220 "clantb.f"
	    }
#line 221 "clantb.f"
	}
#line 222 "clantb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 226 "clantb.f"
	value = 0.;
#line 227 "clantb.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 228 "clantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 229 "clantb.f"
	    i__1 = *n;
#line 229 "clantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 230 "clantb.f"
		if (udiag) {
#line 231 "clantb.f"
		    sum = 1.;
/* Computing MAX */
#line 232 "clantb.f"
		    i__2 = *k + 2 - j;
#line 232 "clantb.f"
		    i__3 = *k;
#line 232 "clantb.f"
		    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 233 "clantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 234 "clantb.f"
/* L90: */
#line 234 "clantb.f"
		    }
#line 235 "clantb.f"
		} else {
#line 236 "clantb.f"
		    sum = 0.;
/* Computing MAX */
#line 237 "clantb.f"
		    i__3 = *k + 2 - j;
#line 237 "clantb.f"
		    i__2 = *k + 1;
#line 237 "clantb.f"
		    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 238 "clantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 239 "clantb.f"
/* L100: */
#line 239 "clantb.f"
		    }
#line 240 "clantb.f"
		}
#line 241 "clantb.f"
		if (value < sum || sisnan_(&sum)) {
#line 241 "clantb.f"
		    value = sum;
#line 241 "clantb.f"
		}
#line 242 "clantb.f"
/* L110: */
#line 242 "clantb.f"
	    }
#line 243 "clantb.f"
	} else {
#line 244 "clantb.f"
	    i__1 = *n;
#line 244 "clantb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "clantb.f"
		if (udiag) {
#line 246 "clantb.f"
		    sum = 1.;
/* Computing MIN */
#line 247 "clantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 247 "clantb.f"
		    i__2 = min(i__3,i__4);
#line 247 "clantb.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 248 "clantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 249 "clantb.f"
/* L120: */
#line 249 "clantb.f"
		    }
#line 250 "clantb.f"
		} else {
#line 251 "clantb.f"
		    sum = 0.;
/* Computing MIN */
#line 252 "clantb.f"
		    i__3 = *n + 1 - j, i__4 = *k + 1;
#line 252 "clantb.f"
		    i__2 = min(i__3,i__4);
#line 252 "clantb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 253 "clantb.f"
			sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 254 "clantb.f"
/* L130: */
#line 254 "clantb.f"
		    }
#line 255 "clantb.f"
		}
#line 256 "clantb.f"
		if (value < sum || sisnan_(&sum)) {
#line 256 "clantb.f"
		    value = sum;
#line 256 "clantb.f"
		}
#line 257 "clantb.f"
/* L140: */
#line 257 "clantb.f"
	    }
#line 258 "clantb.f"
	}
#line 259 "clantb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 263 "clantb.f"
	value = 0.;
#line 264 "clantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 265 "clantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 266 "clantb.f"
		i__1 = *n;
#line 266 "clantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "clantb.f"
		    work[i__] = 1.;
#line 268 "clantb.f"
/* L150: */
#line 268 "clantb.f"
		}
#line 269 "clantb.f"
		i__1 = *n;
#line 269 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 270 "clantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 271 "clantb.f"
		    i__2 = 1, i__3 = j - *k;
#line 271 "clantb.f"
		    i__4 = j - 1;
#line 271 "clantb.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 272 "clantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 273 "clantb.f"
/* L160: */
#line 273 "clantb.f"
		    }
#line 274 "clantb.f"
/* L170: */
#line 274 "clantb.f"
		}
#line 275 "clantb.f"
	    } else {
#line 276 "clantb.f"
		i__1 = *n;
#line 276 "clantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 277 "clantb.f"
		    work[i__] = 0.;
#line 278 "clantb.f"
/* L180: */
#line 278 "clantb.f"
		}
#line 279 "clantb.f"
		i__1 = *n;
#line 279 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 280 "clantb.f"
		    l = *k + 1 - j;
/* Computing MAX */
#line 281 "clantb.f"
		    i__4 = 1, i__2 = j - *k;
#line 281 "clantb.f"
		    i__3 = j;
#line 281 "clantb.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 282 "clantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 283 "clantb.f"
/* L190: */
#line 283 "clantb.f"
		    }
#line 284 "clantb.f"
/* L200: */
#line 284 "clantb.f"
		}
#line 285 "clantb.f"
	    }
#line 286 "clantb.f"
	} else {
#line 287 "clantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 288 "clantb.f"
		i__1 = *n;
#line 288 "clantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "clantb.f"
		    work[i__] = 1.;
#line 290 "clantb.f"
/* L210: */
#line 290 "clantb.f"
		}
#line 291 "clantb.f"
		i__1 = *n;
#line 291 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "clantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 293 "clantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 293 "clantb.f"
		    i__3 = min(i__4,i__2);
#line 293 "clantb.f"
		    for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 294 "clantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 295 "clantb.f"
/* L220: */
#line 295 "clantb.f"
		    }
#line 296 "clantb.f"
/* L230: */
#line 296 "clantb.f"
		}
#line 297 "clantb.f"
	    } else {
#line 298 "clantb.f"
		i__1 = *n;
#line 298 "clantb.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "clantb.f"
		    work[i__] = 0.;
#line 300 "clantb.f"
/* L240: */
#line 300 "clantb.f"
		}
#line 301 "clantb.f"
		i__1 = *n;
#line 301 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
#line 302 "clantb.f"
		    l = 1 - j;
/* Computing MIN */
#line 303 "clantb.f"
		    i__4 = *n, i__2 = j + *k;
#line 303 "clantb.f"
		    i__3 = min(i__4,i__2);
#line 303 "clantb.f"
		    for (i__ = j; i__ <= i__3; ++i__) {
#line 304 "clantb.f"
			work[i__] += z_abs(&ab[l + i__ + j * ab_dim1]);
#line 305 "clantb.f"
/* L250: */
#line 305 "clantb.f"
		    }
#line 306 "clantb.f"
/* L260: */
#line 306 "clantb.f"
		}
#line 307 "clantb.f"
	    }
#line 308 "clantb.f"
	}
#line 309 "clantb.f"
	i__1 = *n;
#line 309 "clantb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 310 "clantb.f"
	    sum = work[i__];
#line 311 "clantb.f"
	    if (value < sum || sisnan_(&sum)) {
#line 311 "clantb.f"
		value = sum;
#line 311 "clantb.f"
	    }
#line 312 "clantb.f"
/* L270: */
#line 312 "clantb.f"
	}
#line 313 "clantb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 317 "clantb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 318 "clantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 319 "clantb.f"
		scale = 1.;
#line 320 "clantb.f"
		sum = (doublereal) (*n);
#line 321 "clantb.f"
		if (*k > 0) {
#line 322 "clantb.f"
		    i__1 = *n;
#line 322 "clantb.f"
		    for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 323 "clantb.f"
			i__4 = j - 1;
#line 323 "clantb.f"
			i__3 = min(i__4,*k);
/* Computing MAX */
#line 323 "clantb.f"
			i__2 = *k + 2 - j;
#line 323 "clantb.f"
			classq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, 
				&scale, &sum);
#line 326 "clantb.f"
/* L280: */
#line 326 "clantb.f"
		    }
#line 327 "clantb.f"
		}
#line 328 "clantb.f"
	    } else {
#line 329 "clantb.f"
		scale = 0.;
#line 330 "clantb.f"
		sum = 1.;
#line 331 "clantb.f"
		i__1 = *n;
#line 331 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 332 "clantb.f"
		    i__4 = j, i__2 = *k + 1;
#line 332 "clantb.f"
		    i__3 = min(i__4,i__2);
/* Computing MAX */
#line 332 "clantb.f"
		    i__5 = *k + 2 - j;
#line 332 "clantb.f"
		    classq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 334 "clantb.f"
/* L290: */
#line 334 "clantb.f"
		}
#line 335 "clantb.f"
	    }
#line 336 "clantb.f"
	} else {
#line 337 "clantb.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 338 "clantb.f"
		scale = 1.;
#line 339 "clantb.f"
		sum = (doublereal) (*n);
#line 340 "clantb.f"
		if (*k > 0) {
#line 341 "clantb.f"
		    i__1 = *n - 1;
#line 341 "clantb.f"
		    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 342 "clantb.f"
			i__4 = *n - j;
#line 342 "clantb.f"
			i__3 = min(i__4,*k);
#line 342 "clantb.f"
			classq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, &
				sum);
#line 344 "clantb.f"
/* L300: */
#line 344 "clantb.f"
		    }
#line 345 "clantb.f"
		}
#line 346 "clantb.f"
	    } else {
#line 347 "clantb.f"
		scale = 0.;
#line 348 "clantb.f"
		sum = 1.;
#line 349 "clantb.f"
		i__1 = *n;
#line 349 "clantb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 350 "clantb.f"
		    i__4 = *n - j + 1, i__2 = *k + 1;
#line 350 "clantb.f"
		    i__3 = min(i__4,i__2);
#line 350 "clantb.f"
		    classq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
#line 352 "clantb.f"
/* L310: */
#line 352 "clantb.f"
		}
#line 353 "clantb.f"
	    }
#line 354 "clantb.f"
	}
#line 355 "clantb.f"
	value = scale * sqrt(sum);
#line 356 "clantb.f"
    }

#line 358 "clantb.f"
    ret_val = value;
#line 359 "clantb.f"
    return ret_val;

/*     End of CLANTB */

} /* clantb_ */

