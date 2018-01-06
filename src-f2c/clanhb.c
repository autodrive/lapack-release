#line 1 "clanhb.f"
/* clanhb.f -- translated by f2c (version 20100827).
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

#line 1 "clanhb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a Hermitian band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANHB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
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
/* > CLANHB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n hermitian band matrix A,  with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return CLANHB */
/* > \verbatim */
/* > */
/* >    CLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANHB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANHB is */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clanhb_(char *norm, char *uplo, integer *n, integer *k, 
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

#line 171 "clanhb.f"
    /* Parameter adjustments */
#line 171 "clanhb.f"
    ab_dim1 = *ldab;
#line 171 "clanhb.f"
    ab_offset = 1 + ab_dim1;
#line 171 "clanhb.f"
    ab -= ab_offset;
#line 171 "clanhb.f"
    --work;
#line 171 "clanhb.f"

#line 171 "clanhb.f"
    /* Function Body */
#line 171 "clanhb.f"
    if (*n == 0) {
#line 172 "clanhb.f"
	value = 0.;
#line 173 "clanhb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 177 "clanhb.f"
	value = 0.;
#line 178 "clanhb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 179 "clanhb.f"
	    i__1 = *n;
#line 179 "clanhb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 180 "clanhb.f"
		i__2 = *k + 2 - j;
#line 180 "clanhb.f"
		i__3 = *k;
#line 180 "clanhb.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 181 "clanhb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 182 "clanhb.f"
		    if (value < sum || sisnan_(&sum)) {
#line 182 "clanhb.f"
			value = sum;
#line 182 "clanhb.f"
		    }
#line 183 "clanhb.f"
/* L10: */
#line 183 "clanhb.f"
		}
#line 184 "clanhb.f"
		i__3 = *k + 1 + j * ab_dim1;
#line 184 "clanhb.f"
		sum = (d__1 = ab[i__3].r, abs(d__1));
#line 185 "clanhb.f"
		if (value < sum || sisnan_(&sum)) {
#line 185 "clanhb.f"
		    value = sum;
#line 185 "clanhb.f"
		}
#line 186 "clanhb.f"
/* L20: */
#line 186 "clanhb.f"
	    }
#line 187 "clanhb.f"
	} else {
#line 188 "clanhb.f"
	    i__1 = *n;
#line 188 "clanhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 189 "clanhb.f"
		i__3 = j * ab_dim1 + 1;
#line 189 "clanhb.f"
		sum = (d__1 = ab[i__3].r, abs(d__1));
#line 190 "clanhb.f"
		if (value < sum || sisnan_(&sum)) {
#line 190 "clanhb.f"
		    value = sum;
#line 190 "clanhb.f"
		}
/* Computing MIN */
#line 191 "clanhb.f"
		i__2 = *n + 1 - j, i__4 = *k + 1;
#line 191 "clanhb.f"
		i__3 = min(i__2,i__4);
#line 191 "clanhb.f"
		for (i__ = 2; i__ <= i__3; ++i__) {
#line 192 "clanhb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 193 "clanhb.f"
		    if (value < sum || sisnan_(&sum)) {
#line 193 "clanhb.f"
			value = sum;
#line 193 "clanhb.f"
		    }
#line 194 "clanhb.f"
/* L30: */
#line 194 "clanhb.f"
		}
#line 195 "clanhb.f"
/* L40: */
#line 195 "clanhb.f"
	    }
#line 196 "clanhb.f"
	}
#line 197 "clanhb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is hermitian). */

#line 202 "clanhb.f"
	value = 0.;
#line 203 "clanhb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 204 "clanhb.f"
	    i__1 = *n;
#line 204 "clanhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "clanhb.f"
		sum = 0.;
#line 206 "clanhb.f"
		l = *k + 1 - j;
/* Computing MAX */
#line 207 "clanhb.f"
		i__3 = 1, i__2 = j - *k;
#line 207 "clanhb.f"
		i__4 = j - 1;
#line 207 "clanhb.f"
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
#line 208 "clanhb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 209 "clanhb.f"
		    sum += absa;
#line 210 "clanhb.f"
		    work[i__] += absa;
#line 211 "clanhb.f"
/* L50: */
#line 211 "clanhb.f"
		}
#line 212 "clanhb.f"
		i__4 = *k + 1 + j * ab_dim1;
#line 212 "clanhb.f"
		work[j] = sum + (d__1 = ab[i__4].r, abs(d__1));
#line 213 "clanhb.f"
/* L60: */
#line 213 "clanhb.f"
	    }
#line 214 "clanhb.f"
	    i__1 = *n;
#line 214 "clanhb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 215 "clanhb.f"
		sum = work[i__];
#line 216 "clanhb.f"
		if (value < sum || sisnan_(&sum)) {
#line 216 "clanhb.f"
		    value = sum;
#line 216 "clanhb.f"
		}
#line 217 "clanhb.f"
/* L70: */
#line 217 "clanhb.f"
	    }
#line 218 "clanhb.f"
	} else {
#line 219 "clanhb.f"
	    i__1 = *n;
#line 219 "clanhb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "clanhb.f"
		work[i__] = 0.;
#line 221 "clanhb.f"
/* L80: */
#line 221 "clanhb.f"
	    }
#line 222 "clanhb.f"
	    i__1 = *n;
#line 222 "clanhb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 223 "clanhb.f"
		i__4 = j * ab_dim1 + 1;
#line 223 "clanhb.f"
		sum = work[j] + (d__1 = ab[i__4].r, abs(d__1));
#line 224 "clanhb.f"
		l = 1 - j;
/* Computing MIN */
#line 225 "clanhb.f"
		i__3 = *n, i__2 = j + *k;
#line 225 "clanhb.f"
		i__4 = min(i__3,i__2);
#line 225 "clanhb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 226 "clanhb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 227 "clanhb.f"
		    sum += absa;
#line 228 "clanhb.f"
		    work[i__] += absa;
#line 229 "clanhb.f"
/* L90: */
#line 229 "clanhb.f"
		}
#line 230 "clanhb.f"
		if (value < sum || sisnan_(&sum)) {
#line 230 "clanhb.f"
		    value = sum;
#line 230 "clanhb.f"
		}
#line 231 "clanhb.f"
/* L100: */
#line 231 "clanhb.f"
	    }
#line 232 "clanhb.f"
	}
#line 233 "clanhb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 237 "clanhb.f"
	scale = 0.;
#line 238 "clanhb.f"
	sum = 1.;
#line 239 "clanhb.f"
	if (*k > 0) {
#line 240 "clanhb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 241 "clanhb.f"
		i__1 = *n;
#line 241 "clanhb.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 242 "clanhb.f"
		    i__3 = j - 1;
#line 242 "clanhb.f"
		    i__4 = min(i__3,*k);
/* Computing MAX */
#line 242 "clanhb.f"
		    i__2 = *k + 2 - j;
#line 242 "clanhb.f"
		    classq_(&i__4, &ab[max(i__2,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 244 "clanhb.f"
/* L110: */
#line 244 "clanhb.f"
		}
#line 245 "clanhb.f"
		l = *k + 1;
#line 246 "clanhb.f"
	    } else {
#line 247 "clanhb.f"
		i__1 = *n - 1;
#line 247 "clanhb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 248 "clanhb.f"
		    i__3 = *n - j;
#line 248 "clanhb.f"
		    i__4 = min(i__3,*k);
#line 248 "clanhb.f"
		    classq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
#line 250 "clanhb.f"
/* L120: */
#line 250 "clanhb.f"
		}
#line 251 "clanhb.f"
		l = 1;
#line 252 "clanhb.f"
	    }
#line 253 "clanhb.f"
	    sum *= 2;
#line 254 "clanhb.f"
	} else {
#line 255 "clanhb.f"
	    l = 1;
#line 256 "clanhb.f"
	}
#line 257 "clanhb.f"
	i__1 = *n;
#line 257 "clanhb.f"
	for (j = 1; j <= i__1; ++j) {
#line 258 "clanhb.f"
	    i__4 = l + j * ab_dim1;
#line 258 "clanhb.f"
	    if (ab[i__4].r != 0.) {
#line 259 "clanhb.f"
		i__4 = l + j * ab_dim1;
#line 259 "clanhb.f"
		absa = (d__1 = ab[i__4].r, abs(d__1));
#line 260 "clanhb.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 261 "clanhb.f"
		    d__1 = scale / absa;
#line 261 "clanhb.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 262 "clanhb.f"
		    scale = absa;
#line 263 "clanhb.f"
		} else {
/* Computing 2nd power */
#line 264 "clanhb.f"
		    d__1 = absa / scale;
#line 264 "clanhb.f"
		    sum += d__1 * d__1;
#line 265 "clanhb.f"
		}
#line 266 "clanhb.f"
	    }
#line 267 "clanhb.f"
/* L130: */
#line 267 "clanhb.f"
	}
#line 268 "clanhb.f"
	value = scale * sqrt(sum);
#line 269 "clanhb.f"
    }

#line 271 "clanhb.f"
    ret_val = value;
#line 272 "clanhb.f"
    return ret_val;

/*     End of CLANHB */

} /* clanhb_ */

