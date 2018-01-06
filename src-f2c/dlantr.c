#line 1 "dlantr.f"
/* dlantr.f -- translated by f2c (version 20100827).
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

#line 1 "dlantr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a trapezoidal or triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANTR  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > trapezoidal or triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANTR */
/* > \verbatim */
/* > */
/* >    DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANTR as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower trapezoidal. */
/* >          = 'U':  Upper trapezoidal */
/* >          = 'L':  Lower trapezoidal */
/* >          Note that A is triangular instead of trapezoidal if M = N. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A has unit diagonal. */
/* >          = 'N':  Non-unit diagonal */
/* >          = 'U':  Unit diagonal */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0, and if */
/* >          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0, and if */
/* >          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The trapezoidal matrix A (A is triangular if M = N). */
/* >          If UPLO = 'U', the leading m by n upper trapezoidal part of */
/* >          the array A contains the upper trapezoidal matrix, and the */
/* >          strictly lower triangular part of A is not referenced. */
/* >          If UPLO = 'L', the leading m by n lower trapezoidal part of */
/* >          the array A contains the lower trapezoidal matrix, and the */
/* >          strictly upper triangular part of A is not referenced.  Note */
/* >          that when DIAG = 'U', the diagonal elements of A are not */
/* >          referenced and are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
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
doublereal dlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
	 doublereal *a, integer *lda, doublereal *work, ftnlen norm_len, 
	ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
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

#line 180 "dlantr.f"
    /* Parameter adjustments */
#line 180 "dlantr.f"
    a_dim1 = *lda;
#line 180 "dlantr.f"
    a_offset = 1 + a_dim1;
#line 180 "dlantr.f"
    a -= a_offset;
#line 180 "dlantr.f"
    --work;
#line 180 "dlantr.f"

#line 180 "dlantr.f"
    /* Function Body */
#line 180 "dlantr.f"
    if (min(*m,*n) == 0) {
#line 181 "dlantr.f"
	value = 0.;
#line 182 "dlantr.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 186 "dlantr.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 187 "dlantr.f"
	    value = 1.;
#line 188 "dlantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "dlantr.f"
		i__1 = *n;
#line 189 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 190 "dlantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 190 "dlantr.f"
		    i__2 = min(i__3,i__4);
#line 190 "dlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 191 "dlantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 192 "dlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 192 "dlantr.f"
			    value = sum;
#line 192 "dlantr.f"
			}
#line 193 "dlantr.f"
/* L10: */
#line 193 "dlantr.f"
		    }
#line 194 "dlantr.f"
/* L20: */
#line 194 "dlantr.f"
		}
#line 195 "dlantr.f"
	    } else {
#line 196 "dlantr.f"
		i__1 = *n;
#line 196 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 197 "dlantr.f"
		    i__2 = *m;
#line 197 "dlantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 198 "dlantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 199 "dlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 199 "dlantr.f"
			    value = sum;
#line 199 "dlantr.f"
			}
#line 200 "dlantr.f"
/* L30: */
#line 200 "dlantr.f"
		    }
#line 201 "dlantr.f"
/* L40: */
#line 201 "dlantr.f"
		}
#line 202 "dlantr.f"
	    }
#line 203 "dlantr.f"
	} else {
#line 204 "dlantr.f"
	    value = 0.;
#line 205 "dlantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 206 "dlantr.f"
		i__1 = *n;
#line 206 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 207 "dlantr.f"
		    i__2 = min(*m,j);
#line 207 "dlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 208 "dlantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 209 "dlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 209 "dlantr.f"
			    value = sum;
#line 209 "dlantr.f"
			}
#line 210 "dlantr.f"
/* L50: */
#line 210 "dlantr.f"
		    }
#line 211 "dlantr.f"
/* L60: */
#line 211 "dlantr.f"
		}
#line 212 "dlantr.f"
	    } else {
#line 213 "dlantr.f"
		i__1 = *n;
#line 213 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 214 "dlantr.f"
		    i__2 = *m;
#line 214 "dlantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 215 "dlantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 216 "dlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 216 "dlantr.f"
			    value = sum;
#line 216 "dlantr.f"
			}
#line 217 "dlantr.f"
/* L70: */
#line 217 "dlantr.f"
		    }
#line 218 "dlantr.f"
/* L80: */
#line 218 "dlantr.f"
		}
#line 219 "dlantr.f"
	    }
#line 220 "dlantr.f"
	}
#line 221 "dlantr.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 225 "dlantr.f"
	value = 0.;
#line 226 "dlantr.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 227 "dlantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 228 "dlantr.f"
	    i__1 = *n;
#line 228 "dlantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 229 "dlantr.f"
		if (udiag && j <= *m) {
#line 230 "dlantr.f"
		    sum = 1.;
#line 231 "dlantr.f"
		    i__2 = j - 1;
#line 231 "dlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 232 "dlantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 233 "dlantr.f"
/* L90: */
#line 233 "dlantr.f"
		    }
#line 234 "dlantr.f"
		} else {
#line 235 "dlantr.f"
		    sum = 0.;
#line 236 "dlantr.f"
		    i__2 = min(*m,j);
#line 236 "dlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 237 "dlantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 238 "dlantr.f"
/* L100: */
#line 238 "dlantr.f"
		    }
#line 239 "dlantr.f"
		}
#line 240 "dlantr.f"
		if (value < sum || disnan_(&sum)) {
#line 240 "dlantr.f"
		    value = sum;
#line 240 "dlantr.f"
		}
#line 241 "dlantr.f"
/* L110: */
#line 241 "dlantr.f"
	    }
#line 242 "dlantr.f"
	} else {
#line 243 "dlantr.f"
	    i__1 = *n;
#line 243 "dlantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 244 "dlantr.f"
		if (udiag) {
#line 245 "dlantr.f"
		    sum = 1.;
#line 246 "dlantr.f"
		    i__2 = *m;
#line 246 "dlantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 247 "dlantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 248 "dlantr.f"
/* L120: */
#line 248 "dlantr.f"
		    }
#line 249 "dlantr.f"
		} else {
#line 250 "dlantr.f"
		    sum = 0.;
#line 251 "dlantr.f"
		    i__2 = *m;
#line 251 "dlantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 252 "dlantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 253 "dlantr.f"
/* L130: */
#line 253 "dlantr.f"
		    }
#line 254 "dlantr.f"
		}
#line 255 "dlantr.f"
		if (value < sum || disnan_(&sum)) {
#line 255 "dlantr.f"
		    value = sum;
#line 255 "dlantr.f"
		}
#line 256 "dlantr.f"
/* L140: */
#line 256 "dlantr.f"
	    }
#line 257 "dlantr.f"
	}
#line 258 "dlantr.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 262 "dlantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 263 "dlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 264 "dlantr.f"
		i__1 = *m;
#line 264 "dlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 265 "dlantr.f"
		    work[i__] = 1.;
#line 266 "dlantr.f"
/* L150: */
#line 266 "dlantr.f"
		}
#line 267 "dlantr.f"
		i__1 = *n;
#line 267 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 268 "dlantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 268 "dlantr.f"
		    i__2 = min(i__3,i__4);
#line 268 "dlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 269 "dlantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 270 "dlantr.f"
/* L160: */
#line 270 "dlantr.f"
		    }
#line 271 "dlantr.f"
/* L170: */
#line 271 "dlantr.f"
		}
#line 272 "dlantr.f"
	    } else {
#line 273 "dlantr.f"
		i__1 = *m;
#line 273 "dlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "dlantr.f"
		    work[i__] = 0.;
#line 275 "dlantr.f"
/* L180: */
#line 275 "dlantr.f"
		}
#line 276 "dlantr.f"
		i__1 = *n;
#line 276 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 277 "dlantr.f"
		    i__2 = min(*m,j);
#line 277 "dlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "dlantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 279 "dlantr.f"
/* L190: */
#line 279 "dlantr.f"
		    }
#line 280 "dlantr.f"
/* L200: */
#line 280 "dlantr.f"
		}
#line 281 "dlantr.f"
	    }
#line 282 "dlantr.f"
	} else {
#line 283 "dlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 284 "dlantr.f"
		i__1 = *n;
#line 284 "dlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "dlantr.f"
		    work[i__] = 1.;
#line 286 "dlantr.f"
/* L210: */
#line 286 "dlantr.f"
		}
#line 287 "dlantr.f"
		i__1 = *m;
#line 287 "dlantr.f"
		for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 288 "dlantr.f"
		    work[i__] = 0.;
#line 289 "dlantr.f"
/* L220: */
#line 289 "dlantr.f"
		}
#line 290 "dlantr.f"
		i__1 = *n;
#line 290 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 291 "dlantr.f"
		    i__2 = *m;
#line 291 "dlantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 292 "dlantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 293 "dlantr.f"
/* L230: */
#line 293 "dlantr.f"
		    }
#line 294 "dlantr.f"
/* L240: */
#line 294 "dlantr.f"
		}
#line 295 "dlantr.f"
	    } else {
#line 296 "dlantr.f"
		i__1 = *m;
#line 296 "dlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "dlantr.f"
		    work[i__] = 0.;
#line 298 "dlantr.f"
/* L250: */
#line 298 "dlantr.f"
		}
#line 299 "dlantr.f"
		i__1 = *n;
#line 299 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 300 "dlantr.f"
		    i__2 = *m;
#line 300 "dlantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 301 "dlantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 302 "dlantr.f"
/* L260: */
#line 302 "dlantr.f"
		    }
#line 303 "dlantr.f"
/* L270: */
#line 303 "dlantr.f"
		}
#line 304 "dlantr.f"
	    }
#line 305 "dlantr.f"
	}
#line 306 "dlantr.f"
	value = 0.;
#line 307 "dlantr.f"
	i__1 = *m;
#line 307 "dlantr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "dlantr.f"
	    sum = work[i__];
#line 309 "dlantr.f"
	    if (value < sum || disnan_(&sum)) {
#line 309 "dlantr.f"
		value = sum;
#line 309 "dlantr.f"
	    }
#line 310 "dlantr.f"
/* L280: */
#line 310 "dlantr.f"
	}
#line 311 "dlantr.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 315 "dlantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 316 "dlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 317 "dlantr.f"
		scale = 1.;
#line 318 "dlantr.f"
		sum = (doublereal) min(*m,*n);
#line 319 "dlantr.f"
		i__1 = *n;
#line 319 "dlantr.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 320 "dlantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 320 "dlantr.f"
		    i__2 = min(i__3,i__4);
#line 320 "dlantr.f"
		    dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 321 "dlantr.f"
/* L290: */
#line 321 "dlantr.f"
		}
#line 322 "dlantr.f"
	    } else {
#line 323 "dlantr.f"
		scale = 0.;
#line 324 "dlantr.f"
		sum = 1.;
#line 325 "dlantr.f"
		i__1 = *n;
#line 325 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "dlantr.f"
		    i__2 = min(*m,j);
#line 326 "dlantr.f"
		    dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 327 "dlantr.f"
/* L300: */
#line 327 "dlantr.f"
		}
#line 328 "dlantr.f"
	    }
#line 329 "dlantr.f"
	} else {
#line 330 "dlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 331 "dlantr.f"
		scale = 1.;
#line 332 "dlantr.f"
		sum = (doublereal) min(*m,*n);
#line 333 "dlantr.f"
		i__1 = *n;
#line 333 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 334 "dlantr.f"
		    i__2 = *m - j;
/* Computing MIN */
#line 334 "dlantr.f"
		    i__3 = *m, i__4 = j + 1;
#line 334 "dlantr.f"
		    dlassq_(&i__2, &a[min(i__3,i__4) + j * a_dim1], &c__1, &
			    scale, &sum);
#line 336 "dlantr.f"
/* L310: */
#line 336 "dlantr.f"
		}
#line 337 "dlantr.f"
	    } else {
#line 338 "dlantr.f"
		scale = 0.;
#line 339 "dlantr.f"
		sum = 1.;
#line 340 "dlantr.f"
		i__1 = *n;
#line 340 "dlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 341 "dlantr.f"
		    i__2 = *m - j + 1;
#line 341 "dlantr.f"
		    dlassq_(&i__2, &a[j + j * a_dim1], &c__1, &scale, &sum);
#line 342 "dlantr.f"
/* L320: */
#line 342 "dlantr.f"
		}
#line 343 "dlantr.f"
	    }
#line 344 "dlantr.f"
	}
#line 345 "dlantr.f"
	value = scale * sqrt(sum);
#line 346 "dlantr.f"
    }

#line 348 "dlantr.f"
    ret_val = value;
#line 349 "dlantr.f"
    return ret_val;

/*     End of DLANTR */

} /* dlantr_ */

