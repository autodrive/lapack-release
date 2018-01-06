#line 1 "clansp.f"
/* clansp.f -- translated by f2c (version 20100827).
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

#line 1 "clansp.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric matrix supplied in packed form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANSP( NORM, UPLO, N, AP, WORK ) */

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
/* > CLANSP  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex symmetric matrix A,  supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return CLANSP */
/* > \verbatim */
/* > */
/* >    CLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANSP as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANSP is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clansp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
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

#line 154 "clansp.f"
    /* Parameter adjustments */
#line 154 "clansp.f"
    --work;
#line 154 "clansp.f"
    --ap;
#line 154 "clansp.f"

#line 154 "clansp.f"
    /* Function Body */
#line 154 "clansp.f"
    if (*n == 0) {
#line 155 "clansp.f"
	value = 0.;
#line 156 "clansp.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 160 "clansp.f"
	value = 0.;
#line 161 "clansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 162 "clansp.f"
	    k = 1;
#line 163 "clansp.f"
	    i__1 = *n;
#line 163 "clansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 164 "clansp.f"
		i__2 = k + j - 1;
#line 164 "clansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 165 "clansp.f"
		    sum = z_abs(&ap[i__]);
#line 166 "clansp.f"
		    if (value < sum || sisnan_(&sum)) {
#line 166 "clansp.f"
			value = sum;
#line 166 "clansp.f"
		    }
#line 167 "clansp.f"
/* L10: */
#line 167 "clansp.f"
		}
#line 168 "clansp.f"
		k += j;
#line 169 "clansp.f"
/* L20: */
#line 169 "clansp.f"
	    }
#line 170 "clansp.f"
	} else {
#line 171 "clansp.f"
	    k = 1;
#line 172 "clansp.f"
	    i__1 = *n;
#line 172 "clansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 173 "clansp.f"
		i__2 = k + *n - j;
#line 173 "clansp.f"
		for (i__ = k; i__ <= i__2; ++i__) {
#line 174 "clansp.f"
		    sum = z_abs(&ap[i__]);
#line 175 "clansp.f"
		    if (value < sum || sisnan_(&sum)) {
#line 175 "clansp.f"
			value = sum;
#line 175 "clansp.f"
		    }
#line 176 "clansp.f"
/* L30: */
#line 176 "clansp.f"
		}
#line 177 "clansp.f"
		k = k + *n - j + 1;
#line 178 "clansp.f"
/* L40: */
#line 178 "clansp.f"
	    }
#line 179 "clansp.f"
	}
#line 180 "clansp.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 185 "clansp.f"
	value = 0.;
#line 186 "clansp.f"
	k = 1;
#line 187 "clansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "clansp.f"
	    i__1 = *n;
#line 188 "clansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 189 "clansp.f"
		sum = 0.;
#line 190 "clansp.f"
		i__2 = j - 1;
#line 190 "clansp.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 191 "clansp.f"
		    absa = z_abs(&ap[k]);
#line 192 "clansp.f"
		    sum += absa;
#line 193 "clansp.f"
		    work[i__] += absa;
#line 194 "clansp.f"
		    ++k;
#line 195 "clansp.f"
/* L50: */
#line 195 "clansp.f"
		}
#line 196 "clansp.f"
		work[j] = sum + z_abs(&ap[k]);
#line 197 "clansp.f"
		++k;
#line 198 "clansp.f"
/* L60: */
#line 198 "clansp.f"
	    }
#line 199 "clansp.f"
	    i__1 = *n;
#line 199 "clansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 200 "clansp.f"
		sum = work[i__];
#line 201 "clansp.f"
		if (value < sum || sisnan_(&sum)) {
#line 201 "clansp.f"
		    value = sum;
#line 201 "clansp.f"
		}
#line 202 "clansp.f"
/* L70: */
#line 202 "clansp.f"
	    }
#line 203 "clansp.f"
	} else {
#line 204 "clansp.f"
	    i__1 = *n;
#line 204 "clansp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 205 "clansp.f"
		work[i__] = 0.;
#line 206 "clansp.f"
/* L80: */
#line 206 "clansp.f"
	    }
#line 207 "clansp.f"
	    i__1 = *n;
#line 207 "clansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 208 "clansp.f"
		sum = work[j] + z_abs(&ap[k]);
#line 209 "clansp.f"
		++k;
#line 210 "clansp.f"
		i__2 = *n;
#line 210 "clansp.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 211 "clansp.f"
		    absa = z_abs(&ap[k]);
#line 212 "clansp.f"
		    sum += absa;
#line 213 "clansp.f"
		    work[i__] += absa;
#line 214 "clansp.f"
		    ++k;
#line 215 "clansp.f"
/* L90: */
#line 215 "clansp.f"
		}
#line 216 "clansp.f"
		if (value < sum || sisnan_(&sum)) {
#line 216 "clansp.f"
		    value = sum;
#line 216 "clansp.f"
		}
#line 217 "clansp.f"
/* L100: */
#line 217 "clansp.f"
	    }
#line 218 "clansp.f"
	}
#line 219 "clansp.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 223 "clansp.f"
	scale = 0.;
#line 224 "clansp.f"
	sum = 1.;
#line 225 "clansp.f"
	k = 2;
#line 226 "clansp.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 227 "clansp.f"
	    i__1 = *n;
#line 227 "clansp.f"
	    for (j = 2; j <= i__1; ++j) {
#line 228 "clansp.f"
		i__2 = j - 1;
#line 228 "clansp.f"
		classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 229 "clansp.f"
		k += j;
#line 230 "clansp.f"
/* L110: */
#line 230 "clansp.f"
	    }
#line 231 "clansp.f"
	} else {
#line 232 "clansp.f"
	    i__1 = *n - 1;
#line 232 "clansp.f"
	    for (j = 1; j <= i__1; ++j) {
#line 233 "clansp.f"
		i__2 = *n - j;
#line 233 "clansp.f"
		classq_(&i__2, &ap[k], &c__1, &scale, &sum);
#line 234 "clansp.f"
		k = k + *n - j + 1;
#line 235 "clansp.f"
/* L120: */
#line 235 "clansp.f"
	    }
#line 236 "clansp.f"
	}
#line 237 "clansp.f"
	sum *= 2;
#line 238 "clansp.f"
	k = 1;
#line 239 "clansp.f"
	i__1 = *n;
#line 239 "clansp.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "clansp.f"
	    i__2 = k;
#line 240 "clansp.f"
	    if (ap[i__2].r != 0.) {
#line 241 "clansp.f"
		i__2 = k;
#line 241 "clansp.f"
		absa = (d__1 = ap[i__2].r, abs(d__1));
#line 242 "clansp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 243 "clansp.f"
		    d__1 = scale / absa;
#line 243 "clansp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 244 "clansp.f"
		    scale = absa;
#line 245 "clansp.f"
		} else {
/* Computing 2nd power */
#line 246 "clansp.f"
		    d__1 = absa / scale;
#line 246 "clansp.f"
		    sum += d__1 * d__1;
#line 247 "clansp.f"
		}
#line 248 "clansp.f"
	    }
#line 249 "clansp.f"
	    if (d_imag(&ap[k]) != 0.) {
#line 250 "clansp.f"
		absa = (d__1 = d_imag(&ap[k]), abs(d__1));
#line 251 "clansp.f"
		if (scale < absa) {
/* Computing 2nd power */
#line 252 "clansp.f"
		    d__1 = scale / absa;
#line 252 "clansp.f"
		    sum = sum * (d__1 * d__1) + 1.;
#line 253 "clansp.f"
		    scale = absa;
#line 254 "clansp.f"
		} else {
/* Computing 2nd power */
#line 255 "clansp.f"
		    d__1 = absa / scale;
#line 255 "clansp.f"
		    sum += d__1 * d__1;
#line 256 "clansp.f"
		}
#line 257 "clansp.f"
	    }
#line 258 "clansp.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 259 "clansp.f"
		k = k + i__ + 1;
#line 260 "clansp.f"
	    } else {
#line 261 "clansp.f"
		k = k + *n - i__ + 1;
#line 262 "clansp.f"
	    }
#line 263 "clansp.f"
/* L130: */
#line 263 "clansp.f"
	}
#line 264 "clansp.f"
	value = scale * sqrt(sum);
#line 265 "clansp.f"
    }

#line 267 "clansp.f"
    ret_val = value;
#line 268 "clansp.f"
    return ret_val;

/*     End of CLANSP */

} /* clansp_ */

