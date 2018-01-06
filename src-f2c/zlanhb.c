#line 1 "zlanhb.f"
/* zlanhb.f -- translated by f2c (version 20100827).
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

#line 1 "zlanhb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a Hermitian band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANHB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANHB( NORM, UPLO, N, K, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
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
/* > ZLANHB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n hermitian band matrix A,  with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return ZLANHB */
/* > \verbatim */
/* > */
/* >    ZLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANHB as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          band matrix A is supplied. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANHB is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of super-diagonals or sub-diagonals of the */
/* >          band matrix A.  K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          The upper or lower triangle of the hermitian band matrix A, */
/* >          stored in the first K+1 rows of AB.  The j-th column of A is */
/* >          stored in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k). */
/* >          Note that the imaginary parts of the diagonal elements need */
/* >          not be set and are assumed to be zero. */
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

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlanhb_(char *norm, char *uplo, integer *n, integer *k, 
	doublecomplex *ab, integer *ldab, doublereal *work, ftnlen norm_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l;
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

#line 171 "zlanhb.f"
    /* Parameter adjustments */
#line 171 "zlanhb.f"
    ab_dim1 = *ldab;
#line 171 "zlanhb.f"
    ab_offset = 1 + ab_dim1;
#line 171 "zlanhb.f"
    ab -= ab_offset;
#line 171 "zlanhb.f"
    --work;
#line 171 "zlanhb.f"

#line 171 "zlanhb.f"
    /* Function Body */
#line 171 "zlanhb.f"
    if (*n == 0) {
#line 172 "zlanhb.f"
	value = 0.;
#line 173 "zlanhb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 177 "zlanhb.f"
	value = 0.;
#line 178 "zlanhb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 179 "zlanhb.f"
	    i__1 = *n;
#line 179 "zlanhb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 180 "zlanhb.f"
		i__2 = *k + 2 - j;
#line 180 "zlanhb.f"
		i__3 = *k;
#line 180 "zlanhb.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 181 "zlanhb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 182 "zlanhb.f"
		    if (value < sum || disnan_(&sum)) {
#line 182 "zlanhb.f"
			value = sum;
#line 182 "zlanhb.f"
		    }
#line 183 "zlanhb.f"
/* L10: */
#line 183 "zlanhb.f"
		}
#line 184 "zlanhb.f"
		i__3 = *k + 1 + j * ab_dim1;
#line 184 "zlanhb.f"
		sum = (d__1 = ab[i__3].r, abs(d__1));
#line 185 "zlanhb.f"
		if (value < sum || disnan_(&sum)) {
#line 185 "zlanhb.f"
		    value = sum;
#line 185 "zlanhb.f"
		}
#line 186 "zlanhb.f"
/* L20: */
#line 186 "zlanhb.f"
	    }
#line 187 "zlanhb.f"
	} else {
#line 188 "zlanhb.f"
	    i__1 = *n;
#line 188 "zlanhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 189 "zlanhb.f"
		i__3 = j * ab_dim1 + 1;
#line 189 "zlanhb.f"
		sum = (d__1 = ab[i__3].r, abs(d__1));
#line 190 "zlanhb.f"
		if (value < sum || disnan_(&sum)) {
#line 190 "zlanhb.f"
		    value = sum;
#line 190 "zlanhb.f"
		}
/* Computing MIN */
#line 191 "zlanhb.f"
		i__2 = *n + 1 - j, i__4 = *k + 1;
#line 191 "zlanhb.f"
		i__3 = min(i__2,i__4);
#line 191 "zlanhb.f"
		for (i__ = 2; i__ <= i__3; ++i__) {
#line 192 "zlanhb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 193 "zlanhb.f"
		    if (value < sum || disnan_(&sum)) {
#line 193 "zlanhb.f"
			value = sum;
#line 193 "zlanhb.f"
		    }
#line 194 "zlanhb.f"
/* L30: */
#line 194 "zlanhb.f"
		}
#line 195 "zlanhb.f"
/* L40: */
#line 195 "zlanhb.f"
	    }
#line 196 "zlanhb.f"
	}
#line 197 "zlanhb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is hermitian). */

#line 202 "zlanhb.f"
	value = 0.;
#line 203 "zlanhb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 204 "zlanhb.f"
	    i__1 = *n;
#line 204 "zlanhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "zlanhb.f"
		sum = 0.;
#line 206 "zlanhb.f"
		l = *k + 1 - j;
/* Computing MAX */
#line 207 "zlanhb.f"
		i__3 = 1, i__2 = j - *k;
#line 207 "zlanhb.f"
		i__4 = j - 1;
#line 207 "zlanhb.f"
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
#line 208 "zlanhb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 209 "zlanhb.f"
		    sum += absa;
#line 210 "zlanhb.f"
		    work[i__] += absa;
#line 211 "zlanhb.f"
/* L50: */
#line 211 "zlanhb.f"
		}
#line 212 "zlanhb.f"
		i__4 = *k + 1 + j * ab_dim1;
#line 212 "zlanhb.f"
		work[j] = sum + (d__1 = ab[i__4].r, abs(d__1));
#line 213 "zlanhb.f"
/* L60: */
#line 213 "zlanhb.f"
	    }
#line 214 "zlanhb.f"
	    i__1 = *n;
#line 214 "zlanhb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 215 "zlanhb.f"
		sum = work[i__];
#line 216 "zlanhb.f"
		if (value < sum || disnan_(&sum)) {
#line 216 "zlanhb.f"
		    value = sum;
#line 216 "zlanhb.f"
		}
#line 217 "zlanhb.f"
/* L70: */
#line 217 "zlanhb.f"
	    }
#line 218 "zlanhb.f"
	} else {
#line 219 "zlanhb.f"
	    i__1 = *n;
#line 219 "zlanhb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "zlanhb.f"
		work[i__] = 0.;
#line 221 "zlanhb.f"
/* L80: */
#line 221 "zlanhb.f"
	    }
#line 222 "zlanhb.f"
	    i__1 = *n;
#line 222 "zlanhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 223 "zlanhb.f"
		i__4 = j * ab_dim1 + 1;
#line 223 "zlanhb.f"
		sum = work[j] + (d__1 = ab[i__4].r, abs(d__1));
#line 224 "zlanhb.f"
		l = 1 - j;
/* Computing MIN */
#line 225 "zlanhb.f"
		i__3 = *n, i__2 = j + *k;
#line 225 "zlanhb.f"
		i__4 = min(i__3,i__2);
#line 225 "zlanhb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 226 "zlanhb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 227 "zlanhb.f"
		    sum += absa;
#line 228 "zlanhb.f"
		    work[i__] += absa;
#line 229 "zlanhb.f"
/* L90: */
#line 229 "zlanhb.f"
		}
#line 230 "zlanhb.f"
		if (value < sum || disnan_(&sum)) {
#line 230 "zlanhb.f"
		    value = sum;
#line 230 "zlanhb.f"
		}
#line 231 "zlanhb.f"
/* L100: */
#line 231 "zlanhb.f"
	    }
#line 232 "zlanhb.f"
	}
#line 233 "zlanhb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 237 "zlanhb.f"
	scale = 0.;
#line 238 "zlanhb.f"
	sum = 1.;
#line 239 "zlanhb.f"
	if (*k > 0) {
#line 240 "zlanhb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 241 "zlanhb.f"
		i__1 = *n;
#line 241 "zlanhb.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 242 "zlanhb.f"
		    i__3 = j - 1;
#line 242 "zlanhb.f"
		    i__4 = min(i__3,*k);
/* Computing MAX */
#line 242 "zlanhb.f"
		    i__2 = *k + 2 - j;
#line 242 "zlanhb.f"
		    zlassq_(&i__4, &ab[max(i__2,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 244 "zlanhb.f"
/* L110: */
#line 244 "zlanhb.f"
		}
#line 245 "zlanhb.f"
		l = *k + 1;
#line 246 "zlanhb.f"
	    } else {
#line 247 "zlanhb.f"
		i__1 = *n - 1;
#line 247 "zlanhb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 248 "zlanhb.f"
		    i__3 = *n - j;
#line 248 "zlanhb.f"
		    i__4 = min(i__3,*k);
#line 248 "zlanhb.f"
		    zlassq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
#line 250 "zlanhb.f"
/* L120: */
#line 250 "zlanhb.f"
		}
#line 251 "zlanhb.f"
		l = 1;
#line 252 "zlanhb.f"
	    }
#line 253 "zlanhb.f"
	    sum *= 2;
#line 254 "zlanhb.f"
	} else {
#line 255 "zlanhb.f"
	    l = 1;
#line 256 "zlanhb.f"
	}
#line 257 "zlanhb.f"
	i__1 = *n;
#line 257 "zlanhb.f"
	for (j = 1; j <= i__1; ++j) {
#line 258 "zlanhb.f"
	    i__4 = l + j * ab_dim1;
#line 258 "zlanhb.f"
	    if (ab[i__4].r != 0.) {
#line 259 "zlanhb.f"
		i__4 = l + j * ab_dim1;
#line 259 "zlanhb.f"
		absa = (d__1 = ab[i__4].r, abs(d__1));
#line 260 "zlanhb.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 261 "zlanhb.f"
		    d__1 = scale / absa;
#line 261 "zlanhb.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 262 "zlanhb.f"
		    scale = absa;
#line 263 "zlanhb.f"
		} else {
/* Computing 2nd power */
#line 264 "zlanhb.f"
		    d__1 = absa / scale;
#line 264 "zlanhb.f"
		    sum += d__1 * d__1;
#line 265 "zlanhb.f"
		}
#line 266 "zlanhb.f"
	    }
#line 267 "zlanhb.f"
/* L130: */
#line 267 "zlanhb.f"
	}
#line 268 "zlanhb.f"
	value = scale * sqrt(sum);
#line 269 "zlanhb.f"
    }

#line 271 "zlanhb.f"
    ret_val = value;
#line 272 "zlanhb.f"
    return ret_val;

/*     End of ZLANHB */

} /* zlanhb_ */

