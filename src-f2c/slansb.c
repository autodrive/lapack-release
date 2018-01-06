#line 1 "slansb.f"
/* slansb.f -- translated by f2c (version 20100827).
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

#line 1 "slansb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANSB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANSB( NORM, UPLO, N, K, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            K, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANSB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n symmetric band matrix A,  with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return SLANSB */
/* > \verbatim */
/* > */
/* >    SLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANSB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANSB is */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
doublereal slansb_(char *norm, char *uplo, integer *n, integer *k, doublereal 
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
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
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

#line 167 "slansb.f"
    /* Parameter adjustments */
#line 167 "slansb.f"
    ab_dim1 = *ldab;
#line 167 "slansb.f"
    ab_offset = 1 + ab_dim1;
#line 167 "slansb.f"
    ab -= ab_offset;
#line 167 "slansb.f"
    --work;
#line 167 "slansb.f"

#line 167 "slansb.f"
    /* Function Body */
#line 167 "slansb.f"
    if (*n == 0) {
#line 168 "slansb.f"
	value = 0.;
#line 169 "slansb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 173 "slansb.f"
	value = 0.;
#line 174 "slansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 175 "slansb.f"
	    i__1 = *n;
#line 175 "slansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 176 "slansb.f"
		i__2 = *k + 2 - j;
#line 176 "slansb.f"
		i__3 = *k + 1;
#line 176 "slansb.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 177 "slansb.f"
		    sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 178 "slansb.f"
		    if (value < sum || sisnan_(&sum)) {
#line 178 "slansb.f"
			value = sum;
#line 178 "slansb.f"
		    }
#line 179 "slansb.f"
/* L10: */
#line 179 "slansb.f"
		}
#line 180 "slansb.f"
/* L20: */
#line 180 "slansb.f"
	    }
#line 181 "slansb.f"
	} else {
#line 182 "slansb.f"
	    i__1 = *n;
#line 182 "slansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 183 "slansb.f"
		i__2 = *n + 1 - j, i__4 = *k + 1;
#line 183 "slansb.f"
		i__3 = min(i__2,i__4);
#line 183 "slansb.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 184 "slansb.f"
		    sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 185 "slansb.f"
		    if (value < sum || sisnan_(&sum)) {
#line 185 "slansb.f"
			value = sum;
#line 185 "slansb.f"
		    }
#line 186 "slansb.f"
/* L30: */
#line 186 "slansb.f"
		}
#line 187 "slansb.f"
/* L40: */
#line 187 "slansb.f"
	    }
#line 188 "slansb.f"
	}
#line 189 "slansb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 194 "slansb.f"
	value = 0.;
#line 195 "slansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 196 "slansb.f"
	    i__1 = *n;
#line 196 "slansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 197 "slansb.f"
		sum = 0.;
#line 198 "slansb.f"
		l = *k + 1 - j;
/* Computing MAX */
#line 199 "slansb.f"
		i__3 = 1, i__2 = j - *k;
#line 199 "slansb.f"
		i__4 = j - 1;
#line 199 "slansb.f"
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
#line 200 "slansb.f"
		    absa = (d__1 = ab[l + i__ + j * ab_dim1], abs(d__1));
#line 201 "slansb.f"
		    sum += absa;
#line 202 "slansb.f"
		    work[i__] += absa;
#line 203 "slansb.f"
/* L50: */
#line 203 "slansb.f"
		}
#line 204 "slansb.f"
		work[j] = sum + (d__1 = ab[*k + 1 + j * ab_dim1], abs(d__1));
#line 205 "slansb.f"
/* L60: */
#line 205 "slansb.f"
	    }
#line 206 "slansb.f"
	    i__1 = *n;
#line 206 "slansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 207 "slansb.f"
		sum = work[i__];
#line 208 "slansb.f"
		if (value < sum || sisnan_(&sum)) {
#line 208 "slansb.f"
		    value = sum;
#line 208 "slansb.f"
		}
#line 209 "slansb.f"
/* L70: */
#line 209 "slansb.f"
	    }
#line 210 "slansb.f"
	} else {
#line 211 "slansb.f"
	    i__1 = *n;
#line 211 "slansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "slansb.f"
		work[i__] = 0.;
#line 213 "slansb.f"
/* L80: */
#line 213 "slansb.f"
	    }
#line 214 "slansb.f"
	    i__1 = *n;
#line 214 "slansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 215 "slansb.f"
		sum = work[j] + (d__1 = ab[j * ab_dim1 + 1], abs(d__1));
#line 216 "slansb.f"
		l = 1 - j;
/* Computing MIN */
#line 217 "slansb.f"
		i__3 = *n, i__2 = j + *k;
#line 217 "slansb.f"
		i__4 = min(i__3,i__2);
#line 217 "slansb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 218 "slansb.f"
		    absa = (d__1 = ab[l + i__ + j * ab_dim1], abs(d__1));
#line 219 "slansb.f"
		    sum += absa;
#line 220 "slansb.f"
		    work[i__] += absa;
#line 221 "slansb.f"
/* L90: */
#line 221 "slansb.f"
		}
#line 222 "slansb.f"
		if (value < sum || sisnan_(&sum)) {
#line 222 "slansb.f"
		    value = sum;
#line 222 "slansb.f"
		}
#line 223 "slansb.f"
/* L100: */
#line 223 "slansb.f"
	    }
#line 224 "slansb.f"
	}
#line 225 "slansb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 229 "slansb.f"
	scale = 0.;
#line 230 "slansb.f"
	sum = 1.;
#line 231 "slansb.f"
	if (*k > 0) {
#line 232 "slansb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 233 "slansb.f"
		i__1 = *n;
#line 233 "slansb.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 234 "slansb.f"
		    i__3 = j - 1;
#line 234 "slansb.f"
		    i__4 = min(i__3,*k);
/* Computing MAX */
#line 234 "slansb.f"
		    i__2 = *k + 2 - j;
#line 234 "slansb.f"
		    slassq_(&i__4, &ab[max(i__2,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 236 "slansb.f"
/* L110: */
#line 236 "slansb.f"
		}
#line 237 "slansb.f"
		l = *k + 1;
#line 238 "slansb.f"
	    } else {
#line 239 "slansb.f"
		i__1 = *n - 1;
#line 239 "slansb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 240 "slansb.f"
		    i__3 = *n - j;
#line 240 "slansb.f"
		    i__4 = min(i__3,*k);
#line 240 "slansb.f"
		    slassq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
#line 242 "slansb.f"
/* L120: */
#line 242 "slansb.f"
		}
#line 243 "slansb.f"
		l = 1;
#line 244 "slansb.f"
	    }
#line 245 "slansb.f"
	    sum *= 2;
#line 246 "slansb.f"
	} else {
#line 247 "slansb.f"
	    l = 1;
#line 248 "slansb.f"
	}
#line 249 "slansb.f"
	slassq_(n, &ab[l + ab_dim1], ldab, &scale, &sum);
#line 250 "slansb.f"
	value = scale * sqrt(sum);
#line 251 "slansb.f"
    }

#line 253 "slansb.f"
    ret_val = value;
#line 254 "slansb.f"
    return ret_val;

/*     End of SLANSB */

} /* slansb_ */

