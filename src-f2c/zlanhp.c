#line 1 "zlanhp.f"
/* zlanhp.f -- translated by f2c (version 20100827).
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

#line 1 "zlanhp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANHP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex Hermitian matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANHP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( * ) */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANHP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex hermitian matrix A,  supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return ZLANHP */
/* > \verbatim */
/* > */
/* >    ZLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANHP as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          hermitian matrix A is supplied. */
/* >          = 'U':  Upper triangular part of A is supplied */
/* >          = 'L':  Lower triangular part of A is supplied */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANHP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the hermitian matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* >          Note that the  imaginary parts of the diagonal elements need */
/* >          not be set and are assumed to be zero. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlanhp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
	doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum, absa, scale;
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

#line 156 "zlanhp.f"
    /* Parameter adjustments */
#line 156 "zlanhp.f"
    --work;
#line 156 "zlanhp.f"
    --ap;
#line 156 "zlanhp.f"

#line 156 "zlanhp.f"
    /* Function Body */
#line 156 "zlanhp.f"
    if (*n == 0) {
#line 157 "zlanhp.f"
	value = 0.;
#line 158 "zlanhp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 162 "zlanhp.f"
	value = 0.;
#line 163 "zlanhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 164 "zlanhp.f"
	    k = 0;
#line 165 "zlanhp.f"
	    i__1 = *n;
#line 165 "zlanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 166 "zlanhp.f"
		i__2 = k + j - 1;
#line 166 "zlanhp.f"
		for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 167 "zlanhp.f"
		    sum = z_abs(&ap[i__]);
#line 168 "zlanhp.f"
		    if (value < sum || disnan_(&sum)) {
#line 168 "zlanhp.f"
			value = sum;
#line 168 "zlanhp.f"
		    }
#line 169 "zlanhp.f"
/* L10: */
#line 169 "zlanhp.f"
		}
#line 170 "zlanhp.f"
		k += j;
#line 171 "zlanhp.f"
		i__2 = k;
#line 171 "zlanhp.f"
		sum = (d__1 = ap[i__2].r, abs(d__1));
#line 172 "zlanhp.f"
		if (value < sum || disnan_(&sum)) {
#line 172 "zlanhp.f"
		    value = sum;
#line 172 "zlanhp.f"
		}
#line 173 "zlanhp.f"
/* L20: */
#line 173 "zlanhp.f"
	    }
#line 174 "zlanhp.f"
	} else {
#line 175 "zlanhp.f"
	    k = 1;
#line 176 "zlanhp.f"
	    i__1 = *n;
#line 176 "zlanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 177 "zlanhp.f"
		i__2 = k;
#line 177 "zlanhp.f"
		sum = (d__1 = ap[i__2].r, abs(d__1));
#line 178 "zlanhp.f"
		if (value < sum || disnan_(&sum)) {
#line 178 "zlanhp.f"
		    value = sum;
#line 178 "zlanhp.f"
		}
#line 179 "zlanhp.f"
		i__2 = k + *n - j;
#line 179 "zlanhp.f"
		for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 180 "zlanhp.f"
		    sum = z_abs(&ap[i__]);
#line 181 "zlanhp.f"
		    if (value < sum || disnan_(&sum)) {
#line 181 "zlanhp.f"
			value = sum;
#line 181 "zlanhp.f"
		    }
#line 182 "zlanhp.f"
/* L30: */
#line 182 "zlanhp.f"
		}
#line 183 "zlanhp.f"
		k = k + *n - j + 1;
#line 184 "zlanhp.f"
/* L40: */
#line 184 "zlanhp.f"
	    }
#line 185 "zlanhp.f"
	}
#line 186 "zlanhp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is hermitian). */

#line 191 "zlanhp.f"
	value = 0.;
#line 192 "zlanhp.f"
	k = 1;
#line 193 "zlanhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 194 "zlanhp.f"
	    i__1 = *n;
#line 194 "zlanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 195 "zlanhp.f"
		sum = 0.;
#line 196 "zlanhp.f"
		i__2 = j - 1;
#line 196 "zlanhp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 197 "zlanhp.f"
		    absa = z_abs(&ap[k]);
#line 198 "zlanhp.f"
		    sum += absa;
#line 199 "zlanhp.f"
		    work[i__] += absa;
#line 200 "zlanhp.f"
		    ++k;
#line 201 "zlanhp.f"
/* L50: */
#line 201 "zlanhp.f"
		}
#line 202 "zlanhp.f"
		i__2 = k;
#line 202 "zlanhp.f"
		work[j] = sum + (d__1 = ap[i__2].r, abs(d__1));
#line 203 "zlanhp.f"
		++k;
#line 204 "zlanhp.f"
/* L60: */
#line 204 "zlanhp.f"
	    }
#line 205 "zlanhp.f"
	    i__1 = *n;
#line 205 "zlanhp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "zlanhp.f"
		sum = work[i__];
#line 207 "zlanhp.f"
		if (value < sum || disnan_(&sum)) {
#line 207 "zlanhp.f"
		    value = sum;
#line 207 "zlanhp.f"
		}
#line 208 "zlanhp.f"
/* L70: */
#line 208 "zlanhp.f"
	    }
#line 209 "zlanhp.f"
	} else {
#line 210 "zlanhp.f"
	    i__1 = *n;
#line 210 "zlanhp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "zlanhp.f"
		work[i__] = 0.;
#line 212 "zlanhp.f"
/* L80: */
#line 212 "zlanhp.f"
	    }
#line 213 "zlanhp.f"
	    i__1 = *n;
#line 213 "zlanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 214 "zlanhp.f"
		i__2 = k;
#line 214 "zlanhp.f"
		sum = work[j] + (d__1 = ap[i__2].r, abs(d__1));
#line 215 "zlanhp.f"
		++k;
#line 216 "zlanhp.f"
		i__2 = *n;
#line 216 "zlanhp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 217 "zlanhp.f"
		    absa = z_abs(&ap[k]);
#line 218 "zlanhp.f"
		    sum += absa;
#line 219 "zlanhp.f"
		    work[i__] += absa;
#line 220 "zlanhp.f"
		    ++k;
#line 221 "zlanhp.f"
/* L90: */
#line 221 "zlanhp.f"
		}
#line 222 "zlanhp.f"
		if (value < sum || disnan_(&sum)) {
#line 222 "zlanhp.f"
		    value = sum;
#line 222 "zlanhp.f"
		}
#line 223 "zlanhp.f"
/* L100: */
#line 223 "zlanhp.f"
	    }
#line 224 "zlanhp.f"
	}
#line 225 "zlanhp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 229 "zlanhp.f"
	scale = 0.;
#line 230 "zlanhp.f"
	sum = 1.;
#line 231 "zlanhp.f"
	k = 2;
#line 232 "zlanhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 233 "zlanhp.f"
	    i__1 = *n;
#line 233 "zlanhp.f"
	    for (j = 2; j <= i__1; ++j) {
#line 234 "zlanhp.f"
		i__2 = j - 1;
#line 234 "zlanhp.f"
		zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 235 "zlanhp.f"
		k += j;
#line 236 "zlanhp.f"
/* L110: */
#line 236 "zlanhp.f"
	    }
#line 237 "zlanhp.f"
	} else {
#line 238 "zlanhp.f"
	    i__1 = *n - 1;
#line 238 "zlanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "zlanhp.f"
		i__2 = *n - j;
#line 239 "zlanhp.f"
		zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 240 "zlanhp.f"
		k = k + *n - j + 1;
#line 241 "zlanhp.f"
/* L120: */
#line 241 "zlanhp.f"
	    }
#line 242 "zlanhp.f"
	}
#line 243 "zlanhp.f"
	sum *= 2;
#line 244 "zlanhp.f"
	k = 1;
#line 245 "zlanhp.f"
	i__1 = *n;
#line 245 "zlanhp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "zlanhp.f"
	    i__2 = k;
#line 246 "zlanhp.f"
	    if (ap[i__2].r != 0.) {
#line 247 "zlanhp.f"
		i__2 = k;
#line 247 "zlanhp.f"
		absa = (d__1 = ap[i__2].r, abs(d__1));
#line 248 "zlanhp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 249 "zlanhp.f"
		    d__1 = scale / absa;
#line 249 "zlanhp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 250 "zlanhp.f"
		    scale = absa;
#line 251 "zlanhp.f"
		} else {
/* Computing 2nd power */
#line 252 "zlanhp.f"
		    d__1 = absa / scale;
#line 252 "zlanhp.f"
		    sum += d__1 * d__1;
#line 253 "zlanhp.f"
		}
#line 254 "zlanhp.f"
	    }
#line 255 "zlanhp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "zlanhp.f"
		k = k + i__ + 1;
#line 257 "zlanhp.f"
	    } else {
#line 258 "zlanhp.f"
		k = k + *n - i__ + 1;
#line 259 "zlanhp.f"
	    }
#line 260 "zlanhp.f"
/* L130: */
#line 260 "zlanhp.f"
	}
#line 261 "zlanhp.f"
	value = scale * sqrt(sum);
#line 262 "zlanhp.f"
    }

#line 264 "zlanhp.f"
    ret_val = value;
#line 265 "zlanhp.f"
    return ret_val;

/*     End of ZLANHP */

} /* zlanhp_ */

