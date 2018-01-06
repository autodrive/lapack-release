#line 1 "clanhe.f"
/* clanhe.f -- translated by f2c (version 20100827).
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

#line 1 "clanhe.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex Hermitian matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANHE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhe.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhe.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhe.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
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
/* > CLANHE  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex hermitian matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANHE */
/* > \verbatim */
/* > */
/* >    CLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANHE as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANHE is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
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

/* > \ingroup complexHEauxiliary */

/*  ===================================================================== */
doublereal clanhe_(char *norm, char *uplo, integer *n, doublecomplex *a, 
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

#line 163 "clanhe.f"
    /* Parameter adjustments */
#line 163 "clanhe.f"
    a_dim1 = *lda;
#line 163 "clanhe.f"
    a_offset = 1 + a_dim1;
#line 163 "clanhe.f"
    a -= a_offset;
#line 163 "clanhe.f"
    --work;
#line 163 "clanhe.f"

#line 163 "clanhe.f"
    /* Function Body */
#line 163 "clanhe.f"
    if (*n == 0) {
#line 164 "clanhe.f"
	value = 0.;
#line 165 "clanhe.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 169 "clanhe.f"
	value = 0.;
#line 170 "clanhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 171 "clanhe.f"
	    i__1 = *n;
#line 171 "clanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 172 "clanhe.f"
		i__2 = j - 1;
#line 172 "clanhe.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 173 "clanhe.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 174 "clanhe.f"
		    if (value < sum || sisnan_(&sum)) {
#line 174 "clanhe.f"
			value = sum;
#line 174 "clanhe.f"
		    }
#line 175 "clanhe.f"
/* L10: */
#line 175 "clanhe.f"
		}
#line 176 "clanhe.f"
		i__2 = j + j * a_dim1;
#line 176 "clanhe.f"
		sum = (d__1 = a[i__2].r, abs(d__1));
#line 177 "clanhe.f"
		if (value < sum || sisnan_(&sum)) {
#line 177 "clanhe.f"
		    value = sum;
#line 177 "clanhe.f"
		}
#line 178 "clanhe.f"
/* L20: */
#line 178 "clanhe.f"
	    }
#line 179 "clanhe.f"
	} else {
#line 180 "clanhe.f"
	    i__1 = *n;
#line 180 "clanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 181 "clanhe.f"
		i__2 = j + j * a_dim1;
#line 181 "clanhe.f"
		sum = (d__1 = a[i__2].r, abs(d__1));
#line 182 "clanhe.f"
		if (value < sum || sisnan_(&sum)) {
#line 182 "clanhe.f"
		    value = sum;
#line 182 "clanhe.f"
		}
#line 183 "clanhe.f"
		i__2 = *n;
#line 183 "clanhe.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 184 "clanhe.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 185 "clanhe.f"
		    if (value < sum || sisnan_(&sum)) {
#line 185 "clanhe.f"
			value = sum;
#line 185 "clanhe.f"
		    }
#line 186 "clanhe.f"
/* L30: */
#line 186 "clanhe.f"
		}
#line 187 "clanhe.f"
/* L40: */
#line 187 "clanhe.f"
	    }
#line 188 "clanhe.f"
	}
#line 189 "clanhe.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is hermitian). */

#line 194 "clanhe.f"
	value = 0.;
#line 195 "clanhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "clanhe.f"
	    i__1 = *n;
#line 196 "clanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 197 "clanhe.f"
		sum = 0.;
#line 198 "clanhe.f"
		i__2 = j - 1;
#line 198 "clanhe.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 199 "clanhe.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 200 "clanhe.f"
		    sum += absa;
#line 201 "clanhe.f"
		    work[i__] += absa;
#line 202 "clanhe.f"
/* L50: */
#line 202 "clanhe.f"
		}
#line 203 "clanhe.f"
		i__2 = j + j * a_dim1;
#line 203 "clanhe.f"
		work[j] = sum + (d__1 = a[i__2].r, abs(d__1));
#line 204 "clanhe.f"
/* L60: */
#line 204 "clanhe.f"
	    }
#line 205 "clanhe.f"
	    i__1 = *n;
#line 205 "clanhe.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "clanhe.f"
		sum = work[i__];
#line 207 "clanhe.f"
		if (value < sum || sisnan_(&sum)) {
#line 207 "clanhe.f"
		    value = sum;
#line 207 "clanhe.f"
		}
#line 208 "clanhe.f"
/* L70: */
#line 208 "clanhe.f"
	    }
#line 209 "clanhe.f"
	} else {
#line 210 "clanhe.f"
	    i__1 = *n;
#line 210 "clanhe.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "clanhe.f"
		work[i__] = 0.;
#line 212 "clanhe.f"
/* L80: */
#line 212 "clanhe.f"
	    }
#line 213 "clanhe.f"
	    i__1 = *n;
#line 213 "clanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 214 "clanhe.f"
		i__2 = j + j * a_dim1;
#line 214 "clanhe.f"
		sum = work[j] + (d__1 = a[i__2].r, abs(d__1));
#line 215 "clanhe.f"
		i__2 = *n;
#line 215 "clanhe.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 216 "clanhe.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 217 "clanhe.f"
		    sum += absa;
#line 218 "clanhe.f"
		    work[i__] += absa;
#line 219 "clanhe.f"
/* L90: */
#line 219 "clanhe.f"
		}
#line 220 "clanhe.f"
		if (value < sum || sisnan_(&sum)) {
#line 220 "clanhe.f"
		    value = sum;
#line 220 "clanhe.f"
		}
#line 221 "clanhe.f"
/* L100: */
#line 221 "clanhe.f"
	    }
#line 222 "clanhe.f"
	}
#line 223 "clanhe.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 227 "clanhe.f"
	scale = 0.;
#line 228 "clanhe.f"
	sum = 1.;
#line 229 "clanhe.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "clanhe.f"
	    i__1 = *n;
#line 230 "clanhe.f"
	    for (j = 2; j <= i__1; ++j) {
#line 231 "clanhe.f"
		i__2 = j - 1;
#line 231 "clanhe.f"
		classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 232 "clanhe.f"
/* L110: */
#line 232 "clanhe.f"
	    }
#line 233 "clanhe.f"
	} else {
#line 234 "clanhe.f"
	    i__1 = *n - 1;
#line 234 "clanhe.f"
	    for (j = 1; j <= i__1; ++j) {
#line 235 "clanhe.f"
		i__2 = *n - j;
#line 235 "clanhe.f"
		classq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
#line 236 "clanhe.f"
/* L120: */
#line 236 "clanhe.f"
	    }
#line 237 "clanhe.f"
	}
#line 238 "clanhe.f"
	sum *= 2;
#line 239 "clanhe.f"
	i__1 = *n;
#line 239 "clanhe.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "clanhe.f"
	    i__2 = i__ + i__ * a_dim1;
#line 240 "clanhe.f"
	    if (a[i__2].r != 0.) {
#line 241 "clanhe.f"
		i__2 = i__ + i__ * a_dim1;
#line 241 "clanhe.f"
		absa = (d__1 = a[i__2].r, abs(d__1));
#line 242 "clanhe.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 243 "clanhe.f"
		    d__1 = scale / absa;
#line 243 "clanhe.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 244 "clanhe.f"
		    scale = absa;
#line 245 "clanhe.f"
		} else {
/* Computing 2nd power */
#line 246 "clanhe.f"
		    d__1 = absa / scale;
#line 246 "clanhe.f"
		    sum += d__1 * d__1;
#line 247 "clanhe.f"
		}
#line 248 "clanhe.f"
	    }
#line 249 "clanhe.f"
/* L130: */
#line 249 "clanhe.f"
	}
#line 250 "clanhe.f"
	value = scale * sqrt(sum);
#line 251 "clanhe.f"
    }

#line 253 "clanhe.f"
    ret_val = value;
#line 254 "clanhe.f"
    return ret_val;

/*     End of CLANHE */

} /* clanhe_ */

