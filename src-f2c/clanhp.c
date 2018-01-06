#line 1 "clanhp.f"
/* clanhp.f -- translated by f2c (version 20100827).
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

#line 1 "clanhp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANHP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex Hermitian matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANHP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANHP( NORM, UPLO, N, AP, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               WORK( * ) */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANHP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex hermitian matrix A,  supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return CLANHP */
/* > \verbatim */
/* > */
/* >    CLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANHP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANHP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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
doublereal clanhp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
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

#line 156 "clanhp.f"
    /* Parameter adjustments */
#line 156 "clanhp.f"
    --work;
#line 156 "clanhp.f"
    --ap;
#line 156 "clanhp.f"

#line 156 "clanhp.f"
    /* Function Body */
#line 156 "clanhp.f"
    if (*n == 0) {
#line 157 "clanhp.f"
	value = 0.;
#line 158 "clanhp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 162 "clanhp.f"
	value = 0.;
#line 163 "clanhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 164 "clanhp.f"
	    k = 0;
#line 165 "clanhp.f"
	    i__1 = *n;
#line 165 "clanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 166 "clanhp.f"
		i__2 = k + j - 1;
#line 166 "clanhp.f"
		for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 167 "clanhp.f"
		    sum = z_abs(&ap[i__]);
#line 168 "clanhp.f"
		    if (value < sum || sisnan_(&sum)) {
#line 168 "clanhp.f"
			value = sum;
#line 168 "clanhp.f"
		    }
#line 169 "clanhp.f"
/* L10: */
#line 169 "clanhp.f"
		}
#line 170 "clanhp.f"
		k += j;
#line 171 "clanhp.f"
		i__2 = k;
#line 171 "clanhp.f"
		sum = (d__1 = ap[i__2].r, abs(d__1));
#line 172 "clanhp.f"
		if (value < sum || sisnan_(&sum)) {
#line 172 "clanhp.f"
		    value = sum;
#line 172 "clanhp.f"
		}
#line 173 "clanhp.f"
/* L20: */
#line 173 "clanhp.f"
	    }
#line 174 "clanhp.f"
	} else {
#line 175 "clanhp.f"
	    k = 1;
#line 176 "clanhp.f"
	    i__1 = *n;
#line 176 "clanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 177 "clanhp.f"
		i__2 = k;
#line 177 "clanhp.f"
		sum = (d__1 = ap[i__2].r, abs(d__1));
#line 178 "clanhp.f"
		if (value < sum || sisnan_(&sum)) {
#line 178 "clanhp.f"
		    value = sum;
#line 178 "clanhp.f"
		}
#line 179 "clanhp.f"
		i__2 = k + *n - j;
#line 179 "clanhp.f"
		for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 180 "clanhp.f"
		    sum = z_abs(&ap[i__]);
#line 181 "clanhp.f"
		    if (value < sum || sisnan_(&sum)) {
#line 181 "clanhp.f"
			value = sum;
#line 181 "clanhp.f"
		    }
#line 182 "clanhp.f"
/* L30: */
#line 182 "clanhp.f"
		}
#line 183 "clanhp.f"
		k = k + *n - j + 1;
#line 184 "clanhp.f"
/* L40: */
#line 184 "clanhp.f"
	    }
#line 185 "clanhp.f"
	}
#line 186 "clanhp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is hermitian). */

#line 191 "clanhp.f"
	value = 0.;
#line 192 "clanhp.f"
	k = 1;
#line 193 "clanhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 194 "clanhp.f"
	    i__1 = *n;
#line 194 "clanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 195 "clanhp.f"
		sum = 0.;
#line 196 "clanhp.f"
		i__2 = j - 1;
#line 196 "clanhp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 197 "clanhp.f"
		    absa = z_abs(&ap[k]);
#line 198 "clanhp.f"
		    sum += absa;
#line 199 "clanhp.f"
		    work[i__] += absa;
#line 200 "clanhp.f"
		    ++k;
#line 201 "clanhp.f"
/* L50: */
#line 201 "clanhp.f"
		}
#line 202 "clanhp.f"
		i__2 = k;
#line 202 "clanhp.f"
		work[j] = sum + (d__1 = ap[i__2].r, abs(d__1));
#line 203 "clanhp.f"
		++k;
#line 204 "clanhp.f"
/* L60: */
#line 204 "clanhp.f"
	    }
#line 205 "clanhp.f"
	    i__1 = *n;
#line 205 "clanhp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "clanhp.f"
		sum = work[i__];
#line 207 "clanhp.f"
		if (value < sum || sisnan_(&sum)) {
#line 207 "clanhp.f"
		    value = sum;
#line 207 "clanhp.f"
		}
#line 208 "clanhp.f"
/* L70: */
#line 208 "clanhp.f"
	    }
#line 209 "clanhp.f"
	} else {
#line 210 "clanhp.f"
	    i__1 = *n;
#line 210 "clanhp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "clanhp.f"
		work[i__] = 0.;
#line 212 "clanhp.f"
/* L80: */
#line 212 "clanhp.f"
	    }
#line 213 "clanhp.f"
	    i__1 = *n;
#line 213 "clanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 214 "clanhp.f"
		i__2 = k;
#line 214 "clanhp.f"
		sum = work[j] + (d__1 = ap[i__2].r, abs(d__1));
#line 215 "clanhp.f"
		++k;
#line 216 "clanhp.f"
		i__2 = *n;
#line 216 "clanhp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 217 "clanhp.f"
		    absa = z_abs(&ap[k]);
#line 218 "clanhp.f"
		    sum += absa;
#line 219 "clanhp.f"
		    work[i__] += absa;
#line 220 "clanhp.f"
		    ++k;
#line 221 "clanhp.f"
/* L90: */
#line 221 "clanhp.f"
		}
#line 222 "clanhp.f"
		if (value < sum || sisnan_(&sum)) {
#line 222 "clanhp.f"
		    value = sum;
#line 222 "clanhp.f"
		}
#line 223 "clanhp.f"
/* L100: */
#line 223 "clanhp.f"
	    }
#line 224 "clanhp.f"
	}
#line 225 "clanhp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 229 "clanhp.f"
	scale = 0.;
#line 230 "clanhp.f"
	sum = 1.;
#line 231 "clanhp.f"
	k = 2;
#line 232 "clanhp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 233 "clanhp.f"
	    i__1 = *n;
#line 233 "clanhp.f"
	    for (j = 2; j <= i__1; ++j) {
#line 234 "clanhp.f"
		i__2 = j - 1;
#line 234 "clanhp.f"
		classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 235 "clanhp.f"
		k += j;
#line 236 "clanhp.f"
/* L110: */
#line 236 "clanhp.f"
	    }
#line 237 "clanhp.f"
	} else {
#line 238 "clanhp.f"
	    i__1 = *n - 1;
#line 238 "clanhp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "clanhp.f"
		i__2 = *n - j;
#line 239 "clanhp.f"
		classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 240 "clanhp.f"
		k = k + *n - j + 1;
#line 241 "clanhp.f"
/* L120: */
#line 241 "clanhp.f"
	    }
#line 242 "clanhp.f"
	}
#line 243 "clanhp.f"
	sum *= 2;
#line 244 "clanhp.f"
	k = 1;
#line 245 "clanhp.f"
	i__1 = *n;
#line 245 "clanhp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "clanhp.f"
	    i__2 = k;
#line 246 "clanhp.f"
	    if (ap[i__2].r != 0.) {
#line 247 "clanhp.f"
		i__2 = k;
#line 247 "clanhp.f"
		absa = (d__1 = ap[i__2].r, abs(d__1));
#line 248 "clanhp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 249 "clanhp.f"
		    d__1 = scale / absa;
#line 249 "clanhp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 250 "clanhp.f"
		    scale = absa;
#line 251 "clanhp.f"
		} else {
/* Computing 2nd power */
#line 252 "clanhp.f"
		    d__1 = absa / scale;
#line 252 "clanhp.f"
		    sum += d__1 * d__1;
#line 253 "clanhp.f"
		}
#line 254 "clanhp.f"
	    }
#line 255 "clanhp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "clanhp.f"
		k = k + i__ + 1;
#line 257 "clanhp.f"
	    } else {
#line 258 "clanhp.f"
		k = k + *n - i__ + 1;
#line 259 "clanhp.f"
	    }
#line 260 "clanhp.f"
/* L130: */
#line 260 "clanhp.f"
	}
#line 261 "clanhp.f"
	value = scale * sqrt(sum);
#line 262 "clanhp.f"
    }

#line 264 "clanhp.f"
    ret_val = value;
#line 265 "clanhp.f"
    return ret_val;

/*     End of CLANHP */

} /* clanhp_ */

