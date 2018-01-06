#line 1 "slansp.f"
/* slansp.f -- translated by f2c (version 20100827).
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

#line 1 "slansp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANSP( NORM, UPLO, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANSP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real symmetric matrix A,  supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return SLANSP */
/* > \verbatim */
/* > */
/* >    SLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANSP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANSP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
doublereal slansp_(char *norm, char *uplo, integer *n, doublereal *ap, 
	doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum, absa, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
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

#line 152 "slansp.f"
    /* Parameter adjustments */
#line 152 "slansp.f"
    --work;
#line 152 "slansp.f"
    --ap;
#line 152 "slansp.f"

#line 152 "slansp.f"
    /* Function Body */
#line 152 "slansp.f"
    if (*n == 0) {
#line 153 "slansp.f"
	value = 0.;
#line 154 "slansp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 158 "slansp.f"
	value = 0.;
#line 159 "slansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 160 "slansp.f"
	    k = 1;
#line 161 "slansp.f"
	    i__1 = *n;
#line 161 "slansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 162 "slansp.f"
		i__2 = k + j - 1;
#line 162 "slansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 163 "slansp.f"
		    sum = (d__1 = ap[i__], abs(d__1));
#line 164 "slansp.f"
		    if (value < sum || sisnan_(&sum)) {
#line 164 "slansp.f"
			value = sum;
#line 164 "slansp.f"
		    }
#line 165 "slansp.f"
/* L10: */
#line 165 "slansp.f"
		}
#line 166 "slansp.f"
		k += j;
#line 167 "slansp.f"
/* L20: */
#line 167 "slansp.f"
	    }
#line 168 "slansp.f"
	} else {
#line 169 "slansp.f"
	    k = 1;
#line 170 "slansp.f"
	    i__1 = *n;
#line 170 "slansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 171 "slansp.f"
		i__2 = k + *n - j;
#line 171 "slansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 172 "slansp.f"
		    sum = (d__1 = ap[i__], abs(d__1));
#line 173 "slansp.f"
		    if (value < sum || sisnan_(&sum)) {
#line 173 "slansp.f"
			value = sum;
#line 173 "slansp.f"
		    }
#line 174 "slansp.f"
/* L30: */
#line 174 "slansp.f"
		}
#line 175 "slansp.f"
		k = k + *n - j + 1;
#line 176 "slansp.f"
/* L40: */
#line 176 "slansp.f"
	    }
#line 177 "slansp.f"
	}
#line 178 "slansp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 183 "slansp.f"
	value = 0.;
#line 184 "slansp.f"
	k = 1;
#line 185 "slansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 186 "slansp.f"
	    i__1 = *n;
#line 186 "slansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 187 "slansp.f"
		sum = 0.;
#line 188 "slansp.f"
		i__2 = j - 1;
#line 188 "slansp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 189 "slansp.f"
		    absa = (d__1 = ap[k], abs(d__1));
#line 190 "slansp.f"
		    sum += absa;
#line 191 "slansp.f"
		    work[i__] += absa;
#line 192 "slansp.f"
		    ++k;
#line 193 "slansp.f"
/* L50: */
#line 193 "slansp.f"
		}
#line 194 "slansp.f"
		work[j] = sum + (d__1 = ap[k], abs(d__1));
#line 195 "slansp.f"
		++k;
#line 196 "slansp.f"
/* L60: */
#line 196 "slansp.f"
	    }
#line 197 "slansp.f"
	    i__1 = *n;
#line 197 "slansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 198 "slansp.f"
		sum = work[i__];
#line 199 "slansp.f"
		if (value < sum || sisnan_(&sum)) {
#line 199 "slansp.f"
		    value = sum;
#line 199 "slansp.f"
		}
#line 200 "slansp.f"
/* L70: */
#line 200 "slansp.f"
	    }
#line 201 "slansp.f"
	} else {
#line 202 "slansp.f"
	    i__1 = *n;
#line 202 "slansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "slansp.f"
		work[i__] = 0.;
#line 204 "slansp.f"
/* L80: */
#line 204 "slansp.f"
	    }
#line 205 "slansp.f"
	    i__1 = *n;
#line 205 "slansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 206 "slansp.f"
		sum = work[j] + (d__1 = ap[k], abs(d__1));
#line 207 "slansp.f"
		++k;
#line 208 "slansp.f"
		i__2 = *n;
#line 208 "slansp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 209 "slansp.f"
		    absa = (d__1 = ap[k], abs(d__1));
#line 210 "slansp.f"
		    sum += absa;
#line 211 "slansp.f"
		    work[i__] += absa;
#line 212 "slansp.f"
		    ++k;
#line 213 "slansp.f"
/* L90: */
#line 213 "slansp.f"
		}
#line 214 "slansp.f"
		if (value < sum || sisnan_(&sum)) {
#line 214 "slansp.f"
		    value = sum;
#line 214 "slansp.f"
		}
#line 215 "slansp.f"
/* L100: */
#line 215 "slansp.f"
	    }
#line 216 "slansp.f"
	}
#line 217 "slansp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 221 "slansp.f"
	scale = 0.;
#line 222 "slansp.f"
	sum = 1.;
#line 223 "slansp.f"
	k = 2;
#line 224 "slansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 225 "slansp.f"
	    i__1 = *n;
#line 225 "slansp.f"
	    for (j = 2; j <= i__1; ++j) {
#line 226 "slansp.f"
		i__2 = j - 1;
#line 226 "slansp.f"
		slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 227 "slansp.f"
		k += j;
#line 228 "slansp.f"
/* L110: */
#line 228 "slansp.f"
	    }
#line 229 "slansp.f"
	} else {
#line 230 "slansp.f"
	    i__1 = *n - 1;
#line 230 "slansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 231 "slansp.f"
		i__2 = *n - j;
#line 231 "slansp.f"
		slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 232 "slansp.f"
		k = k + *n - j + 1;
#line 233 "slansp.f"
/* L120: */
#line 233 "slansp.f"
	    }
#line 234 "slansp.f"
	}
#line 235 "slansp.f"
	sum *= 2;
#line 236 "slansp.f"
	k = 1;
#line 237 "slansp.f"
	i__1 = *n;
#line 237 "slansp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 238 "slansp.f"
	    if (ap[k] != 0.) {
#line 239 "slansp.f"
		absa = (d__1 = ap[k], abs(d__1));
#line 240 "slansp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 241 "slansp.f"
		    d__1 = scale / absa;
#line 241 "slansp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 242 "slansp.f"
		    scale = absa;
#line 243 "slansp.f"
		} else {
/* Computing 2nd power */
#line 244 "slansp.f"
		    d__1 = absa / scale;
#line 244 "slansp.f"
		    sum += d__1 * d__1;
#line 245 "slansp.f"
		}
#line 246 "slansp.f"
	    }
#line 247 "slansp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 248 "slansp.f"
		k = k + i__ + 1;
#line 249 "slansp.f"
	    } else {
#line 250 "slansp.f"
		k = k + *n - i__ + 1;
#line 251 "slansp.f"
	    }
#line 252 "slansp.f"
/* L130: */
#line 252 "slansp.f"
	}
#line 253 "slansp.f"
	value = scale * sqrt(sum);
#line 254 "slansp.f"
    }

#line 256 "slansp.f"
    ret_val = value;
#line 257 "slansp.f"
    return ret_val;

/*     End of SLANSP */

} /* slansp_ */

