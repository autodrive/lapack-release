#line 1 "dlansb.f"
/* dlansb.f -- translated by f2c (version 20100827).
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

#line 1 "dlansb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANSB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANSB( NORM, UPLO, N, K, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            K, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANSB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n symmetric band matrix A,  with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return DLANSB */
/* > \verbatim */
/* > */
/* >    DLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANSB as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          band matrix A is supplied. */
/* >          = 'U':  Upper triangular part is supplied */
/* >          = 'L':  Lower triangular part is supplied */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANSB is */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          The upper or lower triangle of the symmetric band matrix A, */
/* >          stored in the first K+1 rows of AB.  The j-th column of A is */
/* >          stored in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k). */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
doublereal dlansb_(char *norm, char *uplo, integer *n, integer *k, doublereal 
	*ab, integer *ldab, doublereal *work, ftnlen norm_len, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l;
    static doublereal sum, absa, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
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

#line 167 "dlansb.f"
    /* Parameter adjustments */
#line 167 "dlansb.f"
    ab_dim1 = *ldab;
#line 167 "dlansb.f"
    ab_offset = 1 + ab_dim1;
#line 167 "dlansb.f"
    ab -= ab_offset;
#line 167 "dlansb.f"
    --work;
#line 167 "dlansb.f"

#line 167 "dlansb.f"
    /* Function Body */
#line 167 "dlansb.f"
    if (*n == 0) {
#line 168 "dlansb.f"
	value = 0.;
#line 169 "dlansb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 173 "dlansb.f"
	value = 0.;
#line 174 "dlansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 175 "dlansb.f"
	    i__1 = *n;
#line 175 "dlansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 176 "dlansb.f"
		i__2 = *k + 2 - j;
#line 176 "dlansb.f"
		i__3 = *k + 1;
#line 176 "dlansb.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 177 "dlansb.f"
		    sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 178 "dlansb.f"
		    if (value < sum || disnan_(&sum)) {
#line 178 "dlansb.f"
			value = sum;
#line 178 "dlansb.f"
		    }
#line 179 "dlansb.f"
/* L10: */
#line 179 "dlansb.f"
		}
#line 180 "dlansb.f"
/* L20: */
#line 180 "dlansb.f"
	    }
#line 181 "dlansb.f"
	} else {
#line 182 "dlansb.f"
	    i__1 = *n;
#line 182 "dlansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 183 "dlansb.f"
		i__2 = *n + 1 - j, i__4 = *k + 1;
#line 183 "dlansb.f"
		i__3 = min(i__2,i__4);
#line 183 "dlansb.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 184 "dlansb.f"
		    sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 185 "dlansb.f"
		    if (value < sum || disnan_(&sum)) {
#line 185 "dlansb.f"
			value = sum;
#line 185 "dlansb.f"
		    }
#line 186 "dlansb.f"
/* L30: */
#line 186 "dlansb.f"
		}
#line 187 "dlansb.f"
/* L40: */
#line 187 "dlansb.f"
	    }
#line 188 "dlansb.f"
	}
#line 189 "dlansb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 194 "dlansb.f"
	value = 0.;
#line 195 "dlansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "dlansb.f"
	    i__1 = *n;
#line 196 "dlansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 197 "dlansb.f"
		sum = 0.;
#line 198 "dlansb.f"
		l = *k + 1 - j;
/* Computing MAX */
#line 199 "dlansb.f"
		i__3 = 1, i__2 = j - *k;
#line 199 "dlansb.f"
		i__4 = j - 1;
#line 199 "dlansb.f"
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
#line 200 "dlansb.f"
		    absa = (d__1 = ab[l + i__ + j * ab_dim1], abs(d__1));
#line 201 "dlansb.f"
		    sum += absa;
#line 202 "dlansb.f"
		    work[i__] += absa;
#line 203 "dlansb.f"
/* L50: */
#line 203 "dlansb.f"
		}
#line 204 "dlansb.f"
		work[j] = sum + (d__1 = ab[*k + 1 + j * ab_dim1], abs(d__1));
#line 205 "dlansb.f"
/* L60: */
#line 205 "dlansb.f"
	    }
#line 206 "dlansb.f"
	    i__1 = *n;
#line 206 "dlansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 207 "dlansb.f"
		sum = work[i__];
#line 208 "dlansb.f"
		if (value < sum || disnan_(&sum)) {
#line 208 "dlansb.f"
		    value = sum;
#line 208 "dlansb.f"
		}
#line 209 "dlansb.f"
/* L70: */
#line 209 "dlansb.f"
	    }
#line 210 "dlansb.f"
	} else {
#line 211 "dlansb.f"
	    i__1 = *n;
#line 211 "dlansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "dlansb.f"
		work[i__] = 0.;
#line 213 "dlansb.f"
/* L80: */
#line 213 "dlansb.f"
	    }
#line 214 "dlansb.f"
	    i__1 = *n;
#line 214 "dlansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 215 "dlansb.f"
		sum = work[j] + (d__1 = ab[j * ab_dim1 + 1], abs(d__1));
#line 216 "dlansb.f"
		l = 1 - j;
/* Computing MIN */
#line 217 "dlansb.f"
		i__3 = *n, i__2 = j + *k;
#line 217 "dlansb.f"
		i__4 = min(i__3,i__2);
#line 217 "dlansb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 218 "dlansb.f"
		    absa = (d__1 = ab[l + i__ + j * ab_dim1], abs(d__1));
#line 219 "dlansb.f"
		    sum += absa;
#line 220 "dlansb.f"
		    work[i__] += absa;
#line 221 "dlansb.f"
/* L90: */
#line 221 "dlansb.f"
		}
#line 222 "dlansb.f"
		if (value < sum || disnan_(&sum)) {
#line 222 "dlansb.f"
		    value = sum;
#line 222 "dlansb.f"
		}
#line 223 "dlansb.f"
/* L100: */
#line 223 "dlansb.f"
	    }
#line 224 "dlansb.f"
	}
#line 225 "dlansb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 229 "dlansb.f"
	scale = 0.;
#line 230 "dlansb.f"
	sum = 1.;
#line 231 "dlansb.f"
	if (*k > 0) {
#line 232 "dlansb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 233 "dlansb.f"
		i__1 = *n;
#line 233 "dlansb.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 234 "dlansb.f"
		    i__3 = j - 1;
#line 234 "dlansb.f"
		    i__4 = min(i__3,*k);
/* Computing MAX */
#line 234 "dlansb.f"
		    i__2 = *k + 2 - j;
#line 234 "dlansb.f"
		    dlassq_(&i__4, &ab[max(i__2,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 236 "dlansb.f"
/* L110: */
#line 236 "dlansb.f"
		}
#line 237 "dlansb.f"
		l = *k + 1;
#line 238 "dlansb.f"
	    } else {
#line 239 "dlansb.f"
		i__1 = *n - 1;
#line 239 "dlansb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 240 "dlansb.f"
		    i__3 = *n - j;
#line 240 "dlansb.f"
		    i__4 = min(i__3,*k);
#line 240 "dlansb.f"
		    dlassq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
#line 242 "dlansb.f"
/* L120: */
#line 242 "dlansb.f"
		}
#line 243 "dlansb.f"
		l = 1;
#line 244 "dlansb.f"
	    }
#line 245 "dlansb.f"
	    sum *= 2;
#line 246 "dlansb.f"
	} else {
#line 247 "dlansb.f"
	    l = 1;
#line 248 "dlansb.f"
	}
#line 249 "dlansb.f"
	dlassq_(n, &ab[l + ab_dim1], ldab, &scale, &sum);
#line 250 "dlansb.f"
	value = scale * sqrt(sum);
#line 251 "dlansb.f"
    }

#line 253 "dlansb.f"
    ret_val = value;
#line 254 "dlansb.f"
    return ret_val;

/*     End of DLANSB */

} /* dlansb_ */

