#line 1 "zlantr.f"
/* zlantr.f -- translated by f2c (version 20100827).
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

#line 1 "zlantr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a trapezoidal or triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlantr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlantr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlantr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANTR( NORM, UPLO, DIAG, M, N, A, LDA, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANTR  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > trapezoidal or triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \return ZLANTR */
/* > \verbatim */
/* > */
/* >    ZLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANTR as described */
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
/* >          UPLO = 'U', M <= N.  When M = 0, ZLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0, and if */
/* >          UPLO = 'L', N <= M.  When N = 0, ZLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
	 doublecomplex *a, integer *lda, doublereal *work, ftnlen norm_len, 
	ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
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

#line 182 "zlantr.f"
    /* Parameter adjustments */
#line 182 "zlantr.f"
    a_dim1 = *lda;
#line 182 "zlantr.f"
    a_offset = 1 + a_dim1;
#line 182 "zlantr.f"
    a -= a_offset;
#line 182 "zlantr.f"
    --work;
#line 182 "zlantr.f"

#line 182 "zlantr.f"
    /* Function Body */
#line 182 "zlantr.f"
    if (min(*m,*n) == 0) {
#line 183 "zlantr.f"
	value = 0.;
#line 184 "zlantr.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 188 "zlantr.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "zlantr.f"
	    value = 1.;
#line 190 "zlantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 191 "zlantr.f"
		i__1 = *n;
#line 191 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 192 "zlantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 192 "zlantr.f"
		    i__2 = min(i__3,i__4);
#line 192 "zlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 193 "zlantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 194 "zlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 194 "zlantr.f"
			    value = sum;
#line 194 "zlantr.f"
			}
#line 195 "zlantr.f"
/* L10: */
#line 195 "zlantr.f"
		    }
#line 196 "zlantr.f"
/* L20: */
#line 196 "zlantr.f"
		}
#line 197 "zlantr.f"
	    } else {
#line 198 "zlantr.f"
		i__1 = *n;
#line 198 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 199 "zlantr.f"
		    i__2 = *m;
#line 199 "zlantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 200 "zlantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 201 "zlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 201 "zlantr.f"
			    value = sum;
#line 201 "zlantr.f"
			}
#line 202 "zlantr.f"
/* L30: */
#line 202 "zlantr.f"
		    }
#line 203 "zlantr.f"
/* L40: */
#line 203 "zlantr.f"
		}
#line 204 "zlantr.f"
	    }
#line 205 "zlantr.f"
	} else {
#line 206 "zlantr.f"
	    value = 0.;
#line 207 "zlantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 208 "zlantr.f"
		i__1 = *n;
#line 208 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 209 "zlantr.f"
		    i__2 = min(*m,j);
#line 209 "zlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 210 "zlantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 211 "zlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 211 "zlantr.f"
			    value = sum;
#line 211 "zlantr.f"
			}
#line 212 "zlantr.f"
/* L50: */
#line 212 "zlantr.f"
		    }
#line 213 "zlantr.f"
/* L60: */
#line 213 "zlantr.f"
		}
#line 214 "zlantr.f"
	    } else {
#line 215 "zlantr.f"
		i__1 = *n;
#line 215 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 216 "zlantr.f"
		    i__2 = *m;
#line 216 "zlantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 217 "zlantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 218 "zlantr.f"
			if (value < sum || disnan_(&sum)) {
#line 218 "zlantr.f"
			    value = sum;
#line 218 "zlantr.f"
			}
#line 219 "zlantr.f"
/* L70: */
#line 219 "zlantr.f"
		    }
#line 220 "zlantr.f"
/* L80: */
#line 220 "zlantr.f"
		}
#line 221 "zlantr.f"
	    }
#line 222 "zlantr.f"
	}
#line 223 "zlantr.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 227 "zlantr.f"
	value = 0.;
#line 228 "zlantr.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 229 "zlantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "zlantr.f"
	    i__1 = *n;
#line 230 "zlantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 231 "zlantr.f"
		if (udiag && j <= *m) {
#line 232 "zlantr.f"
		    sum = 1.;
#line 233 "zlantr.f"
		    i__2 = j - 1;
#line 233 "zlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 234 "zlantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 235 "zlantr.f"
/* L90: */
#line 235 "zlantr.f"
		    }
#line 236 "zlantr.f"
		} else {
#line 237 "zlantr.f"
		    sum = 0.;
#line 238 "zlantr.f"
		    i__2 = min(*m,j);
#line 238 "zlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 239 "zlantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 240 "zlantr.f"
/* L100: */
#line 240 "zlantr.f"
		    }
#line 241 "zlantr.f"
		}
#line 242 "zlantr.f"
		if (value < sum || disnan_(&sum)) {
#line 242 "zlantr.f"
		    value = sum;
#line 242 "zlantr.f"
		}
#line 243 "zlantr.f"
/* L110: */
#line 243 "zlantr.f"
	    }
#line 244 "zlantr.f"
	} else {
#line 245 "zlantr.f"
	    i__1 = *n;
#line 245 "zlantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 246 "zlantr.f"
		if (udiag) {
#line 247 "zlantr.f"
		    sum = 1.;
#line 248 "zlantr.f"
		    i__2 = *m;
#line 248 "zlantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 249 "zlantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 250 "zlantr.f"
/* L120: */
#line 250 "zlantr.f"
		    }
#line 251 "zlantr.f"
		} else {
#line 252 "zlantr.f"
		    sum = 0.;
#line 253 "zlantr.f"
		    i__2 = *m;
#line 253 "zlantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 254 "zlantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 255 "zlantr.f"
/* L130: */
#line 255 "zlantr.f"
		    }
#line 256 "zlantr.f"
		}
#line 257 "zlantr.f"
		if (value < sum || disnan_(&sum)) {
#line 257 "zlantr.f"
		    value = sum;
#line 257 "zlantr.f"
		}
#line 258 "zlantr.f"
/* L140: */
#line 258 "zlantr.f"
	    }
#line 259 "zlantr.f"
	}
#line 260 "zlantr.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 264 "zlantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 265 "zlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 266 "zlantr.f"
		i__1 = *m;
#line 266 "zlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "zlantr.f"
		    work[i__] = 1.;
#line 268 "zlantr.f"
/* L150: */
#line 268 "zlantr.f"
		}
#line 269 "zlantr.f"
		i__1 = *n;
#line 269 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 270 "zlantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 270 "zlantr.f"
		    i__2 = min(i__3,i__4);
#line 270 "zlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 271 "zlantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 272 "zlantr.f"
/* L160: */
#line 272 "zlantr.f"
		    }
#line 273 "zlantr.f"
/* L170: */
#line 273 "zlantr.f"
		}
#line 274 "zlantr.f"
	    } else {
#line 275 "zlantr.f"
		i__1 = *m;
#line 275 "zlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 276 "zlantr.f"
		    work[i__] = 0.;
#line 277 "zlantr.f"
/* L180: */
#line 277 "zlantr.f"
		}
#line 278 "zlantr.f"
		i__1 = *n;
#line 278 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 279 "zlantr.f"
		    i__2 = min(*m,j);
#line 279 "zlantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 280 "zlantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 281 "zlantr.f"
/* L190: */
#line 281 "zlantr.f"
		    }
#line 282 "zlantr.f"
/* L200: */
#line 282 "zlantr.f"
		}
#line 283 "zlantr.f"
	    }
#line 284 "zlantr.f"
	} else {
#line 285 "zlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 286 "zlantr.f"
		i__1 = *n;
#line 286 "zlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "zlantr.f"
		    work[i__] = 1.;
#line 288 "zlantr.f"
/* L210: */
#line 288 "zlantr.f"
		}
#line 289 "zlantr.f"
		i__1 = *m;
#line 289 "zlantr.f"
		for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 290 "zlantr.f"
		    work[i__] = 0.;
#line 291 "zlantr.f"
/* L220: */
#line 291 "zlantr.f"
		}
#line 292 "zlantr.f"
		i__1 = *n;
#line 292 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 293 "zlantr.f"
		    i__2 = *m;
#line 293 "zlantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 294 "zlantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 295 "zlantr.f"
/* L230: */
#line 295 "zlantr.f"
		    }
#line 296 "zlantr.f"
/* L240: */
#line 296 "zlantr.f"
		}
#line 297 "zlantr.f"
	    } else {
#line 298 "zlantr.f"
		i__1 = *m;
#line 298 "zlantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "zlantr.f"
		    work[i__] = 0.;
#line 300 "zlantr.f"
/* L250: */
#line 300 "zlantr.f"
		}
#line 301 "zlantr.f"
		i__1 = *n;
#line 301 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 302 "zlantr.f"
		    i__2 = *m;
#line 302 "zlantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 303 "zlantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 304 "zlantr.f"
/* L260: */
#line 304 "zlantr.f"
		    }
#line 305 "zlantr.f"
/* L270: */
#line 305 "zlantr.f"
		}
#line 306 "zlantr.f"
	    }
#line 307 "zlantr.f"
	}
#line 308 "zlantr.f"
	value = 0.;
#line 309 "zlantr.f"
	i__1 = *m;
#line 309 "zlantr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 310 "zlantr.f"
	    sum = work[i__];
#line 311 "zlantr.f"
	    if (value < sum || disnan_(&sum)) {
#line 311 "zlantr.f"
		value = sum;
#line 311 "zlantr.f"
	    }
#line 312 "zlantr.f"
/* L280: */
#line 312 "zlantr.f"
	}
#line 313 "zlantr.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 317 "zlantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 318 "zlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 319 "zlantr.f"
		scale = 1.;
#line 320 "zlantr.f"
		sum = (doublereal) min(*m,*n);
#line 321 "zlantr.f"
		i__1 = *n;
#line 321 "zlantr.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 322 "zlantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 322 "zlantr.f"
		    i__2 = min(i__3,i__4);
#line 322 "zlantr.f"
		    zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 323 "zlantr.f"
/* L290: */
#line 323 "zlantr.f"
		}
#line 324 "zlantr.f"
	    } else {
#line 325 "zlantr.f"
		scale = 0.;
#line 326 "zlantr.f"
		sum = 1.;
#line 327 "zlantr.f"
		i__1 = *n;
#line 327 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 328 "zlantr.f"
		    i__2 = min(*m,j);
#line 328 "zlantr.f"
		    zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 329 "zlantr.f"
/* L300: */
#line 329 "zlantr.f"
		}
#line 330 "zlantr.f"
	    }
#line 331 "zlantr.f"
	} else {
#line 332 "zlantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 333 "zlantr.f"
		scale = 1.;
#line 334 "zlantr.f"
		sum = (doublereal) min(*m,*n);
#line 335 "zlantr.f"
		i__1 = *n;
#line 335 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 336 "zlantr.f"
		    i__2 = *m - j;
/* Computing MIN */
#line 336 "zlantr.f"
		    i__3 = *m, i__4 = j + 1;
#line 336 "zlantr.f"
		    zlassq_(&i__2, &a[min(i__3,i__4) + j * a_dim1], &c__1, &
			    scale, &sum);
#line 338 "zlantr.f"
/* L310: */
#line 338 "zlantr.f"
		}
#line 339 "zlantr.f"
	    } else {
#line 340 "zlantr.f"
		scale = 0.;
#line 341 "zlantr.f"
		sum = 1.;
#line 342 "zlantr.f"
		i__1 = *n;
#line 342 "zlantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 343 "zlantr.f"
		    i__2 = *m - j + 1;
#line 343 "zlantr.f"
		    zlassq_(&i__2, &a[j + j * a_dim1], &c__1, &scale, &sum);
#line 344 "zlantr.f"
/* L320: */
#line 344 "zlantr.f"
		}
#line 345 "zlantr.f"
	    }
#line 346 "zlantr.f"
	}
#line 347 "zlantr.f"
	value = scale * sqrt(sum);
#line 348 "zlantr.f"
    }

#line 350 "zlantr.f"
    ret_val = value;
#line 351 "zlantr.f"
    return ret_val;

/*     End of ZLANTR */

} /* zlantr_ */

