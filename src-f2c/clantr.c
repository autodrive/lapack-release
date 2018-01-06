#line 1 "clantr.f"
/* clantr.f -- translated by f2c (version 20100827).
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

#line 1 "clantr.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a trapezoidal or triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANTR( NORM, UPLO, DIAG, M, N, A, LDA, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               WORK( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANTR  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > trapezoidal or triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANTR */
/* > \verbatim */
/* > */
/* >    CLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANTR as described */
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
/* >          UPLO = 'U', M <= N.  When M = 0, CLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0, and if */
/* >          UPLO = 'L', N <= M.  When N = 0, CLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
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

#line 182 "clantr.f"
    /* Parameter adjustments */
#line 182 "clantr.f"
    a_dim1 = *lda;
#line 182 "clantr.f"
    a_offset = 1 + a_dim1;
#line 182 "clantr.f"
    a -= a_offset;
#line 182 "clantr.f"
    --work;
#line 182 "clantr.f"

#line 182 "clantr.f"
    /* Function Body */
#line 182 "clantr.f"
    if (min(*m,*n) == 0) {
#line 183 "clantr.f"
	value = 0.;
#line 184 "clantr.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 188 "clantr.f"
	if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "clantr.f"
	    value = 1.;
#line 190 "clantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 191 "clantr.f"
		i__1 = *n;
#line 191 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 192 "clantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 192 "clantr.f"
		    i__2 = min(i__3,i__4);
#line 192 "clantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 193 "clantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 194 "clantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 194 "clantr.f"
			    value = sum;
#line 194 "clantr.f"
			}
#line 195 "clantr.f"
/* L10: */
#line 195 "clantr.f"
		    }
#line 196 "clantr.f"
/* L20: */
#line 196 "clantr.f"
		}
#line 197 "clantr.f"
	    } else {
#line 198 "clantr.f"
		i__1 = *n;
#line 198 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 199 "clantr.f"
		    i__2 = *m;
#line 199 "clantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 200 "clantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 201 "clantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 201 "clantr.f"
			    value = sum;
#line 201 "clantr.f"
			}
#line 202 "clantr.f"
/* L30: */
#line 202 "clantr.f"
		    }
#line 203 "clantr.f"
/* L40: */
#line 203 "clantr.f"
		}
#line 204 "clantr.f"
	    }
#line 205 "clantr.f"
	} else {
#line 206 "clantr.f"
	    value = 0.;
#line 207 "clantr.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 208 "clantr.f"
		i__1 = *n;
#line 208 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 209 "clantr.f"
		    i__2 = min(*m,j);
#line 209 "clantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 210 "clantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 211 "clantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 211 "clantr.f"
			    value = sum;
#line 211 "clantr.f"
			}
#line 212 "clantr.f"
/* L50: */
#line 212 "clantr.f"
		    }
#line 213 "clantr.f"
/* L60: */
#line 213 "clantr.f"
		}
#line 214 "clantr.f"
	    } else {
#line 215 "clantr.f"
		i__1 = *n;
#line 215 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 216 "clantr.f"
		    i__2 = *m;
#line 216 "clantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 217 "clantr.f"
			sum = z_abs(&a[i__ + j * a_dim1]);
#line 218 "clantr.f"
			if (value < sum || sisnan_(&sum)) {
#line 218 "clantr.f"
			    value = sum;
#line 218 "clantr.f"
			}
#line 219 "clantr.f"
/* L70: */
#line 219 "clantr.f"
		    }
#line 220 "clantr.f"
/* L80: */
#line 220 "clantr.f"
		}
#line 221 "clantr.f"
	    }
#line 222 "clantr.f"
	}
#line 223 "clantr.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 227 "clantr.f"
	value = 0.;
#line 228 "clantr.f"
	udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
#line 229 "clantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "clantr.f"
	    i__1 = *n;
#line 230 "clantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 231 "clantr.f"
		if (udiag && j <= *m) {
#line 232 "clantr.f"
		    sum = 1.;
#line 233 "clantr.f"
		    i__2 = j - 1;
#line 233 "clantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 234 "clantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 235 "clantr.f"
/* L90: */
#line 235 "clantr.f"
		    }
#line 236 "clantr.f"
		} else {
#line 237 "clantr.f"
		    sum = 0.;
#line 238 "clantr.f"
		    i__2 = min(*m,j);
#line 238 "clantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 239 "clantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 240 "clantr.f"
/* L100: */
#line 240 "clantr.f"
		    }
#line 241 "clantr.f"
		}
#line 242 "clantr.f"
		if (value < sum || sisnan_(&sum)) {
#line 242 "clantr.f"
		    value = sum;
#line 242 "clantr.f"
		}
#line 243 "clantr.f"
/* L110: */
#line 243 "clantr.f"
	    }
#line 244 "clantr.f"
	} else {
#line 245 "clantr.f"
	    i__1 = *n;
#line 245 "clantr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 246 "clantr.f"
		if (udiag) {
#line 247 "clantr.f"
		    sum = 1.;
#line 248 "clantr.f"
		    i__2 = *m;
#line 248 "clantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 249 "clantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 250 "clantr.f"
/* L120: */
#line 250 "clantr.f"
		    }
#line 251 "clantr.f"
		} else {
#line 252 "clantr.f"
		    sum = 0.;
#line 253 "clantr.f"
		    i__2 = *m;
#line 253 "clantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 254 "clantr.f"
			sum += z_abs(&a[i__ + j * a_dim1]);
#line 255 "clantr.f"
/* L130: */
#line 255 "clantr.f"
		    }
#line 256 "clantr.f"
		}
#line 257 "clantr.f"
		if (value < sum || sisnan_(&sum)) {
#line 257 "clantr.f"
		    value = sum;
#line 257 "clantr.f"
		}
#line 258 "clantr.f"
/* L140: */
#line 258 "clantr.f"
	    }
#line 259 "clantr.f"
	}
#line 260 "clantr.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 264 "clantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 265 "clantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 266 "clantr.f"
		i__1 = *m;
#line 266 "clantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "clantr.f"
		    work[i__] = 1.;
#line 268 "clantr.f"
/* L150: */
#line 268 "clantr.f"
		}
#line 269 "clantr.f"
		i__1 = *n;
#line 269 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 270 "clantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 270 "clantr.f"
		    i__2 = min(i__3,i__4);
#line 270 "clantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 271 "clantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 272 "clantr.f"
/* L160: */
#line 272 "clantr.f"
		    }
#line 273 "clantr.f"
/* L170: */
#line 273 "clantr.f"
		}
#line 274 "clantr.f"
	    } else {
#line 275 "clantr.f"
		i__1 = *m;
#line 275 "clantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 276 "clantr.f"
		    work[i__] = 0.;
#line 277 "clantr.f"
/* L180: */
#line 277 "clantr.f"
		}
#line 278 "clantr.f"
		i__1 = *n;
#line 278 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 279 "clantr.f"
		    i__2 = min(*m,j);
#line 279 "clantr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 280 "clantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 281 "clantr.f"
/* L190: */
#line 281 "clantr.f"
		    }
#line 282 "clantr.f"
/* L200: */
#line 282 "clantr.f"
		}
#line 283 "clantr.f"
	    }
#line 284 "clantr.f"
	} else {
#line 285 "clantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 286 "clantr.f"
		i__1 = *n;
#line 286 "clantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "clantr.f"
		    work[i__] = 1.;
#line 288 "clantr.f"
/* L210: */
#line 288 "clantr.f"
		}
#line 289 "clantr.f"
		i__1 = *m;
#line 289 "clantr.f"
		for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 290 "clantr.f"
		    work[i__] = 0.;
#line 291 "clantr.f"
/* L220: */
#line 291 "clantr.f"
		}
#line 292 "clantr.f"
		i__1 = *n;
#line 292 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 293 "clantr.f"
		    i__2 = *m;
#line 293 "clantr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 294 "clantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 295 "clantr.f"
/* L230: */
#line 295 "clantr.f"
		    }
#line 296 "clantr.f"
/* L240: */
#line 296 "clantr.f"
		}
#line 297 "clantr.f"
	    } else {
#line 298 "clantr.f"
		i__1 = *m;
#line 298 "clantr.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 299 "clantr.f"
		    work[i__] = 0.;
#line 300 "clantr.f"
/* L250: */
#line 300 "clantr.f"
		}
#line 301 "clantr.f"
		i__1 = *n;
#line 301 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 302 "clantr.f"
		    i__2 = *m;
#line 302 "clantr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 303 "clantr.f"
			work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 304 "clantr.f"
/* L260: */
#line 304 "clantr.f"
		    }
#line 305 "clantr.f"
/* L270: */
#line 305 "clantr.f"
		}
#line 306 "clantr.f"
	    }
#line 307 "clantr.f"
	}
#line 308 "clantr.f"
	value = 0.;
#line 309 "clantr.f"
	i__1 = *m;
#line 309 "clantr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 310 "clantr.f"
	    sum = work[i__];
#line 311 "clantr.f"
	    if (value < sum || sisnan_(&sum)) {
#line 311 "clantr.f"
		value = sum;
#line 311 "clantr.f"
	    }
#line 312 "clantr.f"
/* L280: */
#line 312 "clantr.f"
	}
#line 313 "clantr.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 317 "clantr.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 318 "clantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 319 "clantr.f"
		scale = 1.;
#line 320 "clantr.f"
		sum = (doublereal) min(*m,*n);
#line 321 "clantr.f"
		i__1 = *n;
#line 321 "clantr.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 322 "clantr.f"
		    i__3 = *m, i__4 = j - 1;
#line 322 "clantr.f"
		    i__2 = min(i__3,i__4);
#line 322 "clantr.f"
		    classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 323 "clantr.f"
/* L290: */
#line 323 "clantr.f"
		}
#line 324 "clantr.f"
	    } else {
#line 325 "clantr.f"
		scale = 0.;
#line 326 "clantr.f"
		sum = 1.;
#line 327 "clantr.f"
		i__1 = *n;
#line 327 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 328 "clantr.f"
		    i__2 = min(*m,j);
#line 328 "clantr.f"
		    classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 329 "clantr.f"
/* L300: */
#line 329 "clantr.f"
		}
#line 330 "clantr.f"
	    }
#line 331 "clantr.f"
	} else {
#line 332 "clantr.f"
	    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 333 "clantr.f"
		scale = 1.;
#line 334 "clantr.f"
		sum = (doublereal) min(*m,*n);
#line 335 "clantr.f"
		i__1 = *n;
#line 335 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 336 "clantr.f"
		    i__2 = *m - j;
/* Computing MIN */
#line 336 "clantr.f"
		    i__3 = *m, i__4 = j + 1;
#line 336 "clantr.f"
		    classq_(&i__2, &a[min(i__3,i__4) + j * a_dim1], &c__1, &
			    scale, &sum);
#line 338 "clantr.f"
/* L310: */
#line 338 "clantr.f"
		}
#line 339 "clantr.f"
	    } else {
#line 340 "clantr.f"
		scale = 0.;
#line 341 "clantr.f"
		sum = 1.;
#line 342 "clantr.f"
		i__1 = *n;
#line 342 "clantr.f"
		for (j = 1; j <= i__1; ++j) {
#line 343 "clantr.f"
		    i__2 = *m - j + 1;
#line 343 "clantr.f"
		    classq_(&i__2, &a[j + j * a_dim1], &c__1, &scale, &sum);
#line 344 "clantr.f"
/* L320: */
#line 344 "clantr.f"
		}
#line 345 "clantr.f"
	    }
#line 346 "clantr.f"
	}
#line 347 "clantr.f"
	value = scale * sqrt(sum);
#line 348 "clantr.f"
    }

#line 350 "clantr.f"
    ret_val = value;
#line 351 "clantr.f"
    return ret_val;

/*     End of CLANTR */

} /* clantr_ */

