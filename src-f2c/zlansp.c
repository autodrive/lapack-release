#line 1 "zlansp.f"
/* zlansp.f -- translated by f2c (version 20100827).
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

#line 1 "zlansp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlansp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlansp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlansp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANSP( NORM, UPLO, N, AP, WORK ) */

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
/* > ZLANSP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex symmetric matrix A,  supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return ZLANSP */
/* > \verbatim */
/* > */
/* >    ZLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANSP as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is supplied. */
/* >          = 'U':  Upper triangular part of A is supplied */
/* >          = 'L':  Lower triangular part of A is supplied */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANSP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
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
doublereal zlansp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
	doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
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

#line 154 "zlansp.f"
    /* Parameter adjustments */
#line 154 "zlansp.f"
    --work;
#line 154 "zlansp.f"
    --ap;
#line 154 "zlansp.f"

#line 154 "zlansp.f"
    /* Function Body */
#line 154 "zlansp.f"
    if (*n == 0) {
#line 155 "zlansp.f"
	value = 0.;
#line 156 "zlansp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 160 "zlansp.f"
	value = 0.;
#line 161 "zlansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 162 "zlansp.f"
	    k = 1;
#line 163 "zlansp.f"
	    i__1 = *n;
#line 163 "zlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 164 "zlansp.f"
		i__2 = k + j - 1;
#line 164 "zlansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 165 "zlansp.f"
		    sum = z_abs(&ap[i__]);
#line 166 "zlansp.f"
		    if (value < sum || disnan_(&sum)) {
#line 166 "zlansp.f"
			value = sum;
#line 166 "zlansp.f"
		    }
#line 167 "zlansp.f"
/* L10: */
#line 167 "zlansp.f"
		}
#line 168 "zlansp.f"
		k += j;
#line 169 "zlansp.f"
/* L20: */
#line 169 "zlansp.f"
	    }
#line 170 "zlansp.f"
	} else {
#line 171 "zlansp.f"
	    k = 1;
#line 172 "zlansp.f"
	    i__1 = *n;
#line 172 "zlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 173 "zlansp.f"
		i__2 = k + *n - j;
#line 173 "zlansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 174 "zlansp.f"
		    sum = z_abs(&ap[i__]);
#line 175 "zlansp.f"
		    if (value < sum || disnan_(&sum)) {
#line 175 "zlansp.f"
			value = sum;
#line 175 "zlansp.f"
		    }
#line 176 "zlansp.f"
/* L30: */
#line 176 "zlansp.f"
		}
#line 177 "zlansp.f"
		k = k + *n - j + 1;
#line 178 "zlansp.f"
/* L40: */
#line 178 "zlansp.f"
	    }
#line 179 "zlansp.f"
	}
#line 180 "zlansp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 185 "zlansp.f"
	value = 0.;
#line 186 "zlansp.f"
	k = 1;
#line 187 "zlansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "zlansp.f"
	    i__1 = *n;
#line 188 "zlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 189 "zlansp.f"
		sum = 0.;
#line 190 "zlansp.f"
		i__2 = j - 1;
#line 190 "zlansp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 191 "zlansp.f"
		    absa = z_abs(&ap[k]);
#line 192 "zlansp.f"
		    sum += absa;
#line 193 "zlansp.f"
		    work[i__] += absa;
#line 194 "zlansp.f"
		    ++k;
#line 195 "zlansp.f"
/* L50: */
#line 195 "zlansp.f"
		}
#line 196 "zlansp.f"
		work[j] = sum + z_abs(&ap[k]);
#line 197 "zlansp.f"
		++k;
#line 198 "zlansp.f"
/* L60: */
#line 198 "zlansp.f"
	    }
#line 199 "zlansp.f"
	    i__1 = *n;
#line 199 "zlansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 200 "zlansp.f"
		sum = work[i__];
#line 201 "zlansp.f"
		if (value < sum || disnan_(&sum)) {
#line 201 "zlansp.f"
		    value = sum;
#line 201 "zlansp.f"
		}
#line 202 "zlansp.f"
/* L70: */
#line 202 "zlansp.f"
	    }
#line 203 "zlansp.f"
	} else {
#line 204 "zlansp.f"
	    i__1 = *n;
#line 204 "zlansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 205 "zlansp.f"
		work[i__] = 0.;
#line 206 "zlansp.f"
/* L80: */
#line 206 "zlansp.f"
	    }
#line 207 "zlansp.f"
	    i__1 = *n;
#line 207 "zlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 208 "zlansp.f"
		sum = work[j] + z_abs(&ap[k]);
#line 209 "zlansp.f"
		++k;
#line 210 "zlansp.f"
		i__2 = *n;
#line 210 "zlansp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 211 "zlansp.f"
		    absa = z_abs(&ap[k]);
#line 212 "zlansp.f"
		    sum += absa;
#line 213 "zlansp.f"
		    work[i__] += absa;
#line 214 "zlansp.f"
		    ++k;
#line 215 "zlansp.f"
/* L90: */
#line 215 "zlansp.f"
		}
#line 216 "zlansp.f"
		if (value < sum || disnan_(&sum)) {
#line 216 "zlansp.f"
		    value = sum;
#line 216 "zlansp.f"
		}
#line 217 "zlansp.f"
/* L100: */
#line 217 "zlansp.f"
	    }
#line 218 "zlansp.f"
	}
#line 219 "zlansp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 223 "zlansp.f"
	scale = 0.;
#line 224 "zlansp.f"
	sum = 1.;
#line 225 "zlansp.f"
	k = 2;
#line 226 "zlansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 227 "zlansp.f"
	    i__1 = *n;
#line 227 "zlansp.f"
	    for (j = 2; j <= i__1; ++j) {
#line 228 "zlansp.f"
		i__2 = j - 1;
#line 228 "zlansp.f"
		zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 229 "zlansp.f"
		k += j;
#line 230 "zlansp.f"
/* L110: */
#line 230 "zlansp.f"
	    }
#line 231 "zlansp.f"
	} else {
#line 232 "zlansp.f"
	    i__1 = *n - 1;
#line 232 "zlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 233 "zlansp.f"
		i__2 = *n - j;
#line 233 "zlansp.f"
		zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 234 "zlansp.f"
		k = k + *n - j + 1;
#line 235 "zlansp.f"
/* L120: */
#line 235 "zlansp.f"
	    }
#line 236 "zlansp.f"
	}
#line 237 "zlansp.f"
	sum *= 2;
#line 238 "zlansp.f"
	k = 1;
#line 239 "zlansp.f"
	i__1 = *n;
#line 239 "zlansp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "zlansp.f"
	    i__2 = k;
#line 240 "zlansp.f"
	    if (ap[i__2].r != 0.) {
#line 241 "zlansp.f"
		i__2 = k;
#line 241 "zlansp.f"
		absa = (d__1 = ap[i__2].r, abs(d__1));
#line 242 "zlansp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 243 "zlansp.f"
		    d__1 = scale / absa;
#line 243 "zlansp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 244 "zlansp.f"
		    scale = absa;
#line 245 "zlansp.f"
		} else {
/* Computing 2nd power */
#line 246 "zlansp.f"
		    d__1 = absa / scale;
#line 246 "zlansp.f"
		    sum += d__1 * d__1;
#line 247 "zlansp.f"
		}
#line 248 "zlansp.f"
	    }
#line 249 "zlansp.f"
	    if (d_imag(&ap[k]) != 0.) {
#line 250 "zlansp.f"
		absa = (d__1 = d_imag(&ap[k]), abs(d__1));
#line 251 "zlansp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 252 "zlansp.f"
		    d__1 = scale / absa;
#line 252 "zlansp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 253 "zlansp.f"
		    scale = absa;
#line 254 "zlansp.f"
		} else {
/* Computing 2nd power */
#line 255 "zlansp.f"
		    d__1 = absa / scale;
#line 255 "zlansp.f"
		    sum += d__1 * d__1;
#line 256 "zlansp.f"
		}
#line 257 "zlansp.f"
	    }
#line 258 "zlansp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 259 "zlansp.f"
		k = k + i__ + 1;
#line 260 "zlansp.f"
	    } else {
#line 261 "zlansp.f"
		k = k + *n - i__ + 1;
#line 262 "zlansp.f"
	    }
#line 263 "zlansp.f"
/* L130: */
#line 263 "zlansp.f"
	}
#line 264 "zlansp.f"
	value = scale * sqrt(sum);
#line 265 "zlansp.f"
    }

#line 267 "zlansp.f"
    ret_val = value;
#line 268 "zlansp.f"
    return ret_val;

/*     End of ZLANSP */

} /* zlansp_ */

