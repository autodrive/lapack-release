#line 1 "slantr.f"
/* slantr.f -- translated by f2c (version 20100827).
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

#line 1 "slantr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a trapezoidal or triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANTR( NORM, UPLO, DIAG, M, N, A, LDA, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANTR  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > trapezoidal or triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANTR */
/* > \verbatim */
/* > */
/* >    SLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANTR as described */
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
/* >          UPLO = 'U', M <= N.  When M = 0, SLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0, and if */
/* >          UPLO = 'L', N <= M.  When N = 0, SLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
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
doublereal slantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
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

#line 180 "slantr.f"
    /* Parameter adjustments */
#line 180 "slantr.f"
    a_dim1 = *lda;
#line 180 "slantr.f"
    a_offset = 1 + a_dim1;
#line 180 "slantr.f"
    a -= a_offset;
#line 180 "slantr.f"
    --work;
#line 180 "slantr.f"

#line 180 "slantr.f"
    /* Function Body */
#line 180 "slantr.f"
    if (min(*m,*n) == 0) {
#line 181 "slantr.f"
	value = 0.;
#line 182 "slantr.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 186 "slantr.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 187 "slantr.f"
	    value = 1.;
#line 188 "slantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "slantr.f"
		i__1 = *n;
#line 189 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 190 "slantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 190 "slantr.f"
		    i__2 = min(i__3,i__4);
#line 190 "slantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 191 "slantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 192 "slantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 192 "slantr.f"
			    value = sum;
#line 192 "slantr.f"
			}
#line 193 "slantr.f"
/* L10: */
#line 193 "slantr.f"
		    }
#line 194 "slantr.f"
/* L20: */
#line 194 "slantr.f"
		}
#line 195 "slantr.f"
	    } else {
#line 196 "slantr.f"
		i__1 = *n;
#line 196 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 197 "slantr.f"
		    i__2 = *m;
#line 197 "slantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 198 "slantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 199 "slantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 199 "slantr.f"
			    value = sum;
#line 199 "slantr.f"
			}
#line 200 "slantr.f"
/* L30: */
#line 200 "slantr.f"
		    }
#line 201 "slantr.f"
/* L40: */
#line 201 "slantr.f"
		}
#line 202 "slantr.f"
	    }
#line 203 "slantr.f"
	} else {
#line 204 "slantr.f"
	    value = 0.;
#line 205 "slantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 206 "slantr.f"
		i__1 = *n;
#line 206 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 207 "slantr.f"
		    i__2 = min(*m,j);
#line 207 "slantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 208 "slantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 209 "slantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 209 "slantr.f"
			    value = sum;
#line 209 "slantr.f"
			}
#line 210 "slantr.f"
/* L50: */
#line 210 "slantr.f"
		    }
#line 211 "slantr.f"
/* L60: */
#line 211 "slantr.f"
		}
#line 212 "slantr.f"
	    } else {
#line 213 "slantr.f"
		i__1 = *n;
#line 213 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 214 "slantr.f"
		    i__2 = *m;
#line 214 "slantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 215 "slantr.f"
			sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 216 "slantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 216 "slantr.f"
			    value = sum;
#line 216 "slantr.f"
			}
#line 217 "slantr.f"
/* L70: */
#line 217 "slantr.f"
		    }
#line 218 "slantr.f"
/* L80: */
#line 218 "slantr.f"
		}
#line 219 "slantr.f"
	    }
#line 220 "slantr.f"
	}
#line 221 "slantr.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 225 "slantr.f"
	value = 0.;
#line 226 "slantr.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 227 "slantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 228 "slantr.f"
	    i__1 = *n;
#line 228 "slantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 229 "slantr.f"
		if (udiag && j <= *m) {
#line 230 "slantr.f"
		    sum = 1.;
#line 231 "slantr.f"
		    i__2 = j - 1;
#line 231 "slantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 232 "slantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 233 "slantr.f"
/* L90: */
#line 233 "slantr.f"
		    }
#line 234 "slantr.f"
		} else {
#line 235 "slantr.f"
		    sum = 0.;
#line 236 "slantr.f"
		    i__2 = min(*m,j);
#line 236 "slantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 237 "slantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 238 "slantr.f"
/* L100: */
#line 238 "slantr.f"
		    }
#line 239 "slantr.f"
		}
#line 240 "slantr.f"
		if (value < sum || sisnan_(&sum)) {
#line 240 "slantr.f"
		    value = sum;
#line 240 "slantr.f"
		}
#line 241 "slantr.f"
/* L110: */
#line 241 "slantr.f"
	    }
#line 242 "slantr.f"
	} else {
#line 243 "slantr.f"
	    i__1 = *n;
#line 243 "slantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 244 "slantr.f"
		if (udiag) {
#line 245 "slantr.f"
		    sum = 1.;
#line 246 "slantr.f"
		    i__2 = *m;
#line 246 "slantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 247 "slantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 248 "slantr.f"
/* L120: */
#line 248 "slantr.f"
		    }
#line 249 "slantr.f"
		} else {
#line 250 "slantr.f"
		    sum = 0.;
#line 251 "slantr.f"
		    i__2 = *m;
#line 251 "slantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 252 "slantr.f"
			sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 253 "slantr.f"
/* L130: */
#line 253 "slantr.f"
		    }
#line 254 "slantr.f"
		}
#line 255 "slantr.f"
		if (value < sum || sisnan_(&sum)) {
#line 255 "slantr.f"
		    value = sum;
#line 255 "slantr.f"
		}
#line 256 "slantr.f"
/* L140: */
#line 256 "slantr.f"
	    }
#line 257 "slantr.f"
	}
#line 258 "slantr.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 262 "slantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 263 "slantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 264 "slantr.f"
		i__1 = *m;
#line 264 "slantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 265 "slantr.f"
		    work[i__] = 1.;
#line 266 "slantr.f"
/* L150: */
#line 266 "slantr.f"
		}
#line 267 "slantr.f"
		i__1 = *n;
#line 267 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 268 "slantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 268 "slantr.f"
		    i__2 = min(i__3,i__4);
#line 268 "slantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 269 "slantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 270 "slantr.f"
/* L160: */
#line 270 "slantr.f"
		    }
#line 271 "slantr.f"
/* L170: */
#line 271 "slantr.f"
		}
#line 272 "slantr.f"
	    } else {
#line 273 "slantr.f"
		i__1 = *m;
#line 273 "slantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "slantr.f"
		    work[i__] = 0.;
#line 275 "slantr.f"
/* L180: */
#line 275 "slantr.f"
		}
#line 276 "slantr.f"
		i__1 = *n;
#line 276 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 277 "slantr.f"
		    i__2 = min(*m,j);
#line 277 "slantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "slantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 279 "slantr.f"
/* L190: */
#line 279 "slantr.f"
		    }
#line 280 "slantr.f"
/* L200: */
#line 280 "slantr.f"
		}
#line 281 "slantr.f"
	    }
#line 282 "slantr.f"
	} else {
#line 283 "slantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 284 "slantr.f"
		i__1 = *n;
#line 284 "slantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "slantr.f"
		    work[i__] = 1.;
#line 286 "slantr.f"
/* L210: */
#line 286 "slantr.f"
		}
#line 287 "slantr.f"
		i__1 = *m;
#line 287 "slantr.f"
		for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 288 "slantr.f"
		    work[i__] = 0.;
#line 289 "slantr.f"
/* L220: */
#line 289 "slantr.f"
		}
#line 290 "slantr.f"
		i__1 = *n;
#line 290 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 291 "slantr.f"
		    i__2 = *m;
#line 291 "slantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 292 "slantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 293 "slantr.f"
/* L230: */
#line 293 "slantr.f"
		    }
#line 294 "slantr.f"
/* L240: */
#line 294 "slantr.f"
		}
#line 295 "slantr.f"
	    } else {
#line 296 "slantr.f"
		i__1 = *m;
#line 296 "slantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "slantr.f"
		    work[i__] = 0.;
#line 298 "slantr.f"
/* L250: */
#line 298 "slantr.f"
		}
#line 299 "slantr.f"
		i__1 = *n;
#line 299 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 300 "slantr.f"
		    i__2 = *m;
#line 300 "slantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 301 "slantr.f"
			work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 302 "slantr.f"
/* L260: */
#line 302 "slantr.f"
		    }
#line 303 "slantr.f"
/* L270: */
#line 303 "slantr.f"
		}
#line 304 "slantr.f"
	    }
#line 305 "slantr.f"
	}
#line 306 "slantr.f"
	value = 0.;
#line 307 "slantr.f"
	i__1 = *m;
#line 307 "slantr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "slantr.f"
	    sum = work[i__];
#line 309 "slantr.f"
	    if (value < sum || sisnan_(&sum)) {
#line 309 "slantr.f"
		value = sum;
#line 309 "slantr.f"
	    }
#line 310 "slantr.f"
/* L280: */
#line 310 "slantr.f"
	}
#line 311 "slantr.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 315 "slantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 316 "slantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 317 "slantr.f"
		scale = 1.;
#line 318 "slantr.f"
		sum = (doublereal) min(*m,*n);
#line 319 "slantr.f"
		i__1 = *n;
#line 319 "slantr.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 320 "slantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 320 "slantr.f"
		    i__2 = min(i__3,i__4);
#line 320 "slantr.f"
		    slassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 321 "slantr.f"
/* L290: */
#line 321 "slantr.f"
		}
#line 322 "slantr.f"
	    } else {
#line 323 "slantr.f"
		scale = 0.;
#line 324 "slantr.f"
		sum = 1.;
#line 325 "slantr.f"
		i__1 = *n;
#line 325 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "slantr.f"
		    i__2 = min(*m,j);
#line 326 "slantr.f"
		    slassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 327 "slantr.f"
/* L300: */
#line 327 "slantr.f"
		}
#line 328 "slantr.f"
	    }
#line 329 "slantr.f"
	} else {
#line 330 "slantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 331 "slantr.f"
		scale = 1.;
#line 332 "slantr.f"
		sum = (doublereal) min(*m,*n);
#line 333 "slantr.f"
		i__1 = *n;
#line 333 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 334 "slantr.f"
		    i__2 = *m - j;
/* Computing MIN */
#line 334 "slantr.f"
		    i__3 = *m, i__4 = j + 1;
#line 334 "slantr.f"
		    slassq_(&i__2, &a[min(i__3,i__4) + j * a_dim1], &c__1, &
			    scale, &sum);
#line 336 "slantr.f"
/* L310: */
#line 336 "slantr.f"
		}
#line 337 "slantr.f"
	    } else {
#line 338 "slantr.f"
		scale = 0.;
#line 339 "slantr.f"
		sum = 1.;
#line 340 "slantr.f"
		i__1 = *n;
#line 340 "slantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 341 "slantr.f"
		    i__2 = *m - j + 1;
#line 341 "slantr.f"
		    slassq_(&i__2, &a[j + j * a_dim1], &c__1, &scale, &sum);
#line 342 "slantr.f"
/* L320: */
#line 342 "slantr.f"
		}
#line 343 "slantr.f"
	    }
#line 344 "slantr.f"
	}
#line 345 "slantr.f"
	value = scale * sqrt(sum);
#line 346 "slantr.f"
    }

#line 348 "slantr.f"
    ret_val = value;
#line 349 "slantr.f"
    return ret_val;

/*     End of SLANTR */

} /* slantr_ */

