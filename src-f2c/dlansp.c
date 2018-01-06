#line 1 "dlansp.f"
/* dlansp.f -- translated by f2c (version 20100827).
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

#line 1 "dlansp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANSP( NORM, UPLO, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANSP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real symmetric matrix A,  supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return DLANSP */
/* > \verbatim */
/* > */
/* >    DLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANSP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANSP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
doublereal dlansp_(char *norm, char *uplo, integer *n, doublereal *ap, 
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
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
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

#line 152 "dlansp.f"
    /* Parameter adjustments */
#line 152 "dlansp.f"
    --work;
#line 152 "dlansp.f"
    --ap;
#line 152 "dlansp.f"

#line 152 "dlansp.f"
    /* Function Body */
#line 152 "dlansp.f"
    if (*n == 0) {
#line 153 "dlansp.f"
	value = 0.;
#line 154 "dlansp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 158 "dlansp.f"
	value = 0.;
#line 159 "dlansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 160 "dlansp.f"
	    k = 1;
#line 161 "dlansp.f"
	    i__1 = *n;
#line 161 "dlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 162 "dlansp.f"
		i__2 = k + j - 1;
#line 162 "dlansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 163 "dlansp.f"
		    sum = (d__1 = ap[i__], abs(d__1));
#line 164 "dlansp.f"
		    if (value < sum || disnan_(&sum)) {
#line 164 "dlansp.f"
			value = sum;
#line 164 "dlansp.f"
		    }
#line 165 "dlansp.f"
/* L10: */
#line 165 "dlansp.f"
		}
#line 166 "dlansp.f"
		k += j;
#line 167 "dlansp.f"
/* L20: */
#line 167 "dlansp.f"
	    }
#line 168 "dlansp.f"
	} else {
#line 169 "dlansp.f"
	    k = 1;
#line 170 "dlansp.f"
	    i__1 = *n;
#line 170 "dlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 171 "dlansp.f"
		i__2 = k + *n - j;
#line 171 "dlansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 172 "dlansp.f"
		    sum = (d__1 = ap[i__], abs(d__1));
#line 173 "dlansp.f"
		    if (value < sum || disnan_(&sum)) {
#line 173 "dlansp.f"
			value = sum;
#line 173 "dlansp.f"
		    }
#line 174 "dlansp.f"
/* L30: */
#line 174 "dlansp.f"
		}
#line 175 "dlansp.f"
		k = k + *n - j + 1;
#line 176 "dlansp.f"
/* L40: */
#line 176 "dlansp.f"
	    }
#line 177 "dlansp.f"
	}
#line 178 "dlansp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 183 "dlansp.f"
	value = 0.;
#line 184 "dlansp.f"
	k = 1;
#line 185 "dlansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 186 "dlansp.f"
	    i__1 = *n;
#line 186 "dlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 187 "dlansp.f"
		sum = 0.;
#line 188 "dlansp.f"
		i__2 = j - 1;
#line 188 "dlansp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 189 "dlansp.f"
		    absa = (d__1 = ap[k], abs(d__1));
#line 190 "dlansp.f"
		    sum += absa;
#line 191 "dlansp.f"
		    work[i__] += absa;
#line 192 "dlansp.f"
		    ++k;
#line 193 "dlansp.f"
/* L50: */
#line 193 "dlansp.f"
		}
#line 194 "dlansp.f"
		work[j] = sum + (d__1 = ap[k], abs(d__1));
#line 195 "dlansp.f"
		++k;
#line 196 "dlansp.f"
/* L60: */
#line 196 "dlansp.f"
	    }
#line 197 "dlansp.f"
	    i__1 = *n;
#line 197 "dlansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 198 "dlansp.f"
		sum = work[i__];
#line 199 "dlansp.f"
		if (value < sum || disnan_(&sum)) {
#line 199 "dlansp.f"
		    value = sum;
#line 199 "dlansp.f"
		}
#line 200 "dlansp.f"
/* L70: */
#line 200 "dlansp.f"
	    }
#line 201 "dlansp.f"
	} else {
#line 202 "dlansp.f"
	    i__1 = *n;
#line 202 "dlansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "dlansp.f"
		work[i__] = 0.;
#line 204 "dlansp.f"
/* L80: */
#line 204 "dlansp.f"
	    }
#line 205 "dlansp.f"
	    i__1 = *n;
#line 205 "dlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 206 "dlansp.f"
		sum = work[j] + (d__1 = ap[k], abs(d__1));
#line 207 "dlansp.f"
		++k;
#line 208 "dlansp.f"
		i__2 = *n;
#line 208 "dlansp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 209 "dlansp.f"
		    absa = (d__1 = ap[k], abs(d__1));
#line 210 "dlansp.f"
		    sum += absa;
#line 211 "dlansp.f"
		    work[i__] += absa;
#line 212 "dlansp.f"
		    ++k;
#line 213 "dlansp.f"
/* L90: */
#line 213 "dlansp.f"
		}
#line 214 "dlansp.f"
		if (value < sum || disnan_(&sum)) {
#line 214 "dlansp.f"
		    value = sum;
#line 214 "dlansp.f"
		}
#line 215 "dlansp.f"
/* L100: */
#line 215 "dlansp.f"
	    }
#line 216 "dlansp.f"
	}
#line 217 "dlansp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 221 "dlansp.f"
	scale = 0.;
#line 222 "dlansp.f"
	sum = 1.;
#line 223 "dlansp.f"
	k = 2;
#line 224 "dlansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 225 "dlansp.f"
	    i__1 = *n;
#line 225 "dlansp.f"
	    for (j = 2; j <= i__1; ++j) {
#line 226 "dlansp.f"
		i__2 = j - 1;
#line 226 "dlansp.f"
		dlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 227 "dlansp.f"
		k += j;
#line 228 "dlansp.f"
/* L110: */
#line 228 "dlansp.f"
	    }
#line 229 "dlansp.f"
	} else {
#line 230 "dlansp.f"
	    i__1 = *n - 1;
#line 230 "dlansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 231 "dlansp.f"
		i__2 = *n - j;
#line 231 "dlansp.f"
		dlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 232 "dlansp.f"
		k = k + *n - j + 1;
#line 233 "dlansp.f"
/* L120: */
#line 233 "dlansp.f"
	    }
#line 234 "dlansp.f"
	}
#line 235 "dlansp.f"
	sum *= 2;
#line 236 "dlansp.f"
	k = 1;
#line 237 "dlansp.f"
	i__1 = *n;
#line 237 "dlansp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 238 "dlansp.f"
	    if (ap[k] != 0.) {
#line 239 "dlansp.f"
		absa = (d__1 = ap[k], abs(d__1));
#line 240 "dlansp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 241 "dlansp.f"
		    d__1 = scale / absa;
#line 241 "dlansp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 242 "dlansp.f"
		    scale = absa;
#line 243 "dlansp.f"
		} else {
/* Computing 2nd power */
#line 244 "dlansp.f"
		    d__1 = absa / scale;
#line 244 "dlansp.f"
		    sum += d__1 * d__1;
#line 245 "dlansp.f"
		}
#line 246 "dlansp.f"
	    }
#line 247 "dlansp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 248 "dlansp.f"
		k = k + i__ + 1;
#line 249 "dlansp.f"
	    } else {
#line 250 "dlansp.f"
		k = k + *n - i__ + 1;
#line 251 "dlansp.f"
	    }
#line 252 "dlansp.f"
/* L130: */
#line 252 "dlansp.f"
	}
#line 253 "dlansp.f"
	value = scale * sqrt(sum);
#line 254 "dlansp.f"
    }

#line 256 "dlansp.f"
    ret_val = value;
#line 257 "dlansp.f"
    return ret_val;

/*     End of DLANSP */

} /* dlansp_ */

