#line 1 "zlanhe.f"
/* zlanhe.f -- translated by f2c (version 20100827).
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

#line 1 "zlanhe.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex Hermitian matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANHE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhe.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhe.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhe.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
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
/* > ZLANHE  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex hermitian matrix A. */
/* > \endverbatim */
/* > */
/* > \return ZLANHE */
/* > \verbatim */
/* > */
/* >    ZLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANHE as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          hermitian matrix A is to be referenced. */
/* >          = 'U':  Upper triangular part of A is referenced */
/* >          = 'L':  Lower triangular part of A is referenced */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANHE is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The hermitian matrix A.  If UPLO = 'U', the leading n by n */
/* >          upper triangular part of A contains the upper triangular part */
/* >          of the matrix A, and the strictly lower triangular part of A */
/* >          is not referenced.  If UPLO = 'L', the leading n by n lower */
/* >          triangular part of A contains the lower triangular part of */
/* >          the matrix A, and the strictly upper triangular part of A is */
/* >          not referenced. Note that the imaginary parts of the diagonal */
/* >          elements need not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise, */
/* >          WORK is not referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16HEauxiliary */

/*  ===================================================================== */
doublereal zlanhe_(char *norm, char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal sum, absa, scale;
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

#line 163 "zlanhe.f"
    /* Parameter adjustments */
#line 163 "zlanhe.f"
    a_dim1 = *lda;
#line 163 "zlanhe.f"
    a_offset = 1 + a_dim1;
#line 163 "zlanhe.f"
    a -= a_offset;
#line 163 "zlanhe.f"
    --work;
#line 163 "zlanhe.f"

#line 163 "zlanhe.f"
    /* Function Body */
#line 163 "zlanhe.f"
    if (*n == 0) {
#line 164 "zlanhe.f"
	value = 0.;
#line 165 "zlanhe.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 169 "zlanhe.f"
	value = 0.;
#line 170 "zlanhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 171 "zlanhe.f"
	    i__1 = *n;
#line 171 "zlanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 172 "zlanhe.f"
		i__2 = j - 1;
#line 172 "zlanhe.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 173 "zlanhe.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 174 "zlanhe.f"
		    if (value < sum || disnan_(&sum)) {
#line 174 "zlanhe.f"
			value = sum;
#line 174 "zlanhe.f"
		    }
#line 175 "zlanhe.f"
/* L10: */
#line 175 "zlanhe.f"
		}
#line 176 "zlanhe.f"
		i__2 = j + j * a_dim1;
#line 176 "zlanhe.f"
		sum = (d__1 = a[i__2].r, abs(d__1));
#line 177 "zlanhe.f"
		if (value < sum || disnan_(&sum)) {
#line 177 "zlanhe.f"
		    value = sum;
#line 177 "zlanhe.f"
		}
#line 178 "zlanhe.f"
/* L20: */
#line 178 "zlanhe.f"
	    }
#line 179 "zlanhe.f"
	} else {
#line 180 "zlanhe.f"
	    i__1 = *n;
#line 180 "zlanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 181 "zlanhe.f"
		i__2 = j + j * a_dim1;
#line 181 "zlanhe.f"
		sum = (d__1 = a[i__2].r, abs(d__1));
#line 182 "zlanhe.f"
		if (value < sum || disnan_(&sum)) {
#line 182 "zlanhe.f"
		    value = sum;
#line 182 "zlanhe.f"
		}
#line 183 "zlanhe.f"
		i__2 = *n;
#line 183 "zlanhe.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 184 "zlanhe.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 185 "zlanhe.f"
		    if (value < sum || disnan_(&sum)) {
#line 185 "zlanhe.f"
			value = sum;
#line 185 "zlanhe.f"
		    }
#line 186 "zlanhe.f"
/* L30: */
#line 186 "zlanhe.f"
		}
#line 187 "zlanhe.f"
/* L40: */
#line 187 "zlanhe.f"
	    }
#line 188 "zlanhe.f"
	}
#line 189 "zlanhe.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is hermitian). */

#line 194 "zlanhe.f"
	value = 0.;
#line 195 "zlanhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "zlanhe.f"
	    i__1 = *n;
#line 196 "zlanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 197 "zlanhe.f"
		sum = 0.;
#line 198 "zlanhe.f"
		i__2 = j - 1;
#line 198 "zlanhe.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 199 "zlanhe.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 200 "zlanhe.f"
		    sum += absa;
#line 201 "zlanhe.f"
		    work[i__] += absa;
#line 202 "zlanhe.f"
/* L50: */
#line 202 "zlanhe.f"
		}
#line 203 "zlanhe.f"
		i__2 = j + j * a_dim1;
#line 203 "zlanhe.f"
		work[j] = sum + (d__1 = a[i__2].r, abs(d__1));
#line 204 "zlanhe.f"
/* L60: */
#line 204 "zlanhe.f"
	    }
#line 205 "zlanhe.f"
	    i__1 = *n;
#line 205 "zlanhe.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "zlanhe.f"
		sum = work[i__];
#line 207 "zlanhe.f"
		if (value < sum || disnan_(&sum)) {
#line 207 "zlanhe.f"
		    value = sum;
#line 207 "zlanhe.f"
		}
#line 208 "zlanhe.f"
/* L70: */
#line 208 "zlanhe.f"
	    }
#line 209 "zlanhe.f"
	} else {
#line 210 "zlanhe.f"
	    i__1 = *n;
#line 210 "zlanhe.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "zlanhe.f"
		work[i__] = 0.;
#line 212 "zlanhe.f"
/* L80: */
#line 212 "zlanhe.f"
	    }
#line 213 "zlanhe.f"
	    i__1 = *n;
#line 213 "zlanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 214 "zlanhe.f"
		i__2 = j + j * a_dim1;
#line 214 "zlanhe.f"
		sum = work[j] + (d__1 = a[i__2].r, abs(d__1));
#line 215 "zlanhe.f"
		i__2 = *n;
#line 215 "zlanhe.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 216 "zlanhe.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 217 "zlanhe.f"
		    sum += absa;
#line 218 "zlanhe.f"
		    work[i__] += absa;
#line 219 "zlanhe.f"
/* L90: */
#line 219 "zlanhe.f"
		}
#line 220 "zlanhe.f"
		if (value < sum || disnan_(&sum)) {
#line 220 "zlanhe.f"
		    value = sum;
#line 220 "zlanhe.f"
		}
#line 221 "zlanhe.f"
/* L100: */
#line 221 "zlanhe.f"
	    }
#line 222 "zlanhe.f"
	}
#line 223 "zlanhe.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 227 "zlanhe.f"
	scale = 0.;
#line 228 "zlanhe.f"
	sum = 1.;
#line 229 "zlanhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "zlanhe.f"
	    i__1 = *n;
#line 230 "zlanhe.f"
	    for (j = 2; j <= i__1; ++j) {
#line 231 "zlanhe.f"
		i__2 = j - 1;
#line 231 "zlanhe.f"
		zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 232 "zlanhe.f"
/* L110: */
#line 232 "zlanhe.f"
	    }
#line 233 "zlanhe.f"
	} else {
#line 234 "zlanhe.f"
	    i__1 = *n - 1;
#line 234 "zlanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 235 "zlanhe.f"
		i__2 = *n - j;
#line 235 "zlanhe.f"
		zlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
#line 236 "zlanhe.f"
/* L120: */
#line 236 "zlanhe.f"
	    }
#line 237 "zlanhe.f"
	}
#line 238 "zlanhe.f"
	sum *= 2;
#line 239 "zlanhe.f"
	i__1 = *n;
#line 239 "zlanhe.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "zlanhe.f"
	    i__2 = i__ + i__ * a_dim1;
#line 240 "zlanhe.f"
	    if (a[i__2].r != 0.) {
#line 241 "zlanhe.f"
		i__2 = i__ + i__ * a_dim1;
#line 241 "zlanhe.f"
		absa = (d__1 = a[i__2].r, abs(d__1));
#line 242 "zlanhe.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 243 "zlanhe.f"
		    d__1 = scale / absa;
#line 243 "zlanhe.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 244 "zlanhe.f"
		    scale = absa;
#line 245 "zlanhe.f"
		} else {
/* Computing 2nd power */
#line 246 "zlanhe.f"
		    d__1 = absa / scale;
#line 246 "zlanhe.f"
		    sum += d__1 * d__1;
#line 247 "zlanhe.f"
		}
#line 248 "zlanhe.f"
	    }
#line 249 "zlanhe.f"
/* L130: */
#line 249 "zlanhe.f"
	}
#line 250 "zlanhe.f"
	value = scale * sqrt(sum);
#line 251 "zlanhe.f"
    }

#line 253 "zlanhe.f"
    ret_val = value;
#line 254 "zlanhe.f"
    return ret_val;

/*     End of ZLANHE */

} /* zlanhe_ */

