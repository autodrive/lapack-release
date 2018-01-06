#line 1 "zlansb.f"
/* zlansb.f -- translated by f2c (version 20100827).
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

#line 1 "zlansb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANSB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlansb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlansb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlansb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANSB( NORM, UPLO, N, K, AB, LDAB, */
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
/* > ZLANSB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n symmetric band matrix A,  with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return ZLANSB */
/* > \verbatim */
/* > */
/* >    ZLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANSB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANSB is */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlansb_(char *norm, char *uplo, integer *n, integer *k, 
	doublecomplex *ab, integer *ldab, doublereal *work, ftnlen norm_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val;

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

#line 169 "zlansb.f"
    /* Parameter adjustments */
#line 169 "zlansb.f"
    ab_dim1 = *ldab;
#line 169 "zlansb.f"
    ab_offset = 1 + ab_dim1;
#line 169 "zlansb.f"
    ab -= ab_offset;
#line 169 "zlansb.f"
    --work;
#line 169 "zlansb.f"

#line 169 "zlansb.f"
    /* Function Body */
#line 169 "zlansb.f"
    if (*n == 0) {
#line 170 "zlansb.f"
	value = 0.;
#line 171 "zlansb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 175 "zlansb.f"
	value = 0.;
#line 176 "zlansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 177 "zlansb.f"
	    i__1 = *n;
#line 177 "zlansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 178 "zlansb.f"
		i__2 = *k + 2 - j;
#line 178 "zlansb.f"
		i__3 = *k + 1;
#line 178 "zlansb.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 179 "zlansb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 180 "zlansb.f"
		    if (value < sum || disnan_(&sum)) {
#line 180 "zlansb.f"
			value = sum;
#line 180 "zlansb.f"
		    }
#line 181 "zlansb.f"
/* L10: */
#line 181 "zlansb.f"
		}
#line 182 "zlansb.f"
/* L20: */
#line 182 "zlansb.f"
	    }
#line 183 "zlansb.f"
	} else {
#line 184 "zlansb.f"
	    i__1 = *n;
#line 184 "zlansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 185 "zlansb.f"
		i__2 = *n + 1 - j, i__4 = *k + 1;
#line 185 "zlansb.f"
		i__3 = min(i__2,i__4);
#line 185 "zlansb.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 186 "zlansb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 187 "zlansb.f"
		    if (value < sum || disnan_(&sum)) {
#line 187 "zlansb.f"
			value = sum;
#line 187 "zlansb.f"
		    }
#line 188 "zlansb.f"
/* L30: */
#line 188 "zlansb.f"
		}
#line 189 "zlansb.f"
/* L40: */
#line 189 "zlansb.f"
	    }
#line 190 "zlansb.f"
	}
#line 191 "zlansb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 196 "zlansb.f"
	value = 0.;
#line 197 "zlansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 198 "zlansb.f"
	    i__1 = *n;
#line 198 "zlansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 199 "zlansb.f"
		sum = 0.;
#line 200 "zlansb.f"
		l = *k + 1 - j;
/* Computing MAX */
#line 201 "zlansb.f"
		i__3 = 1, i__2 = j - *k;
#line 201 "zlansb.f"
		i__4 = j - 1;
#line 201 "zlansb.f"
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
#line 202 "zlansb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 203 "zlansb.f"
		    sum += absa;
#line 204 "zlansb.f"
		    work[i__] += absa;
#line 205 "zlansb.f"
/* L50: */
#line 205 "zlansb.f"
		}
#line 206 "zlansb.f"
		work[j] = sum + z_abs(&ab[*k + 1 + j * ab_dim1]);
#line 207 "zlansb.f"
/* L60: */
#line 207 "zlansb.f"
	    }
#line 208 "zlansb.f"
	    i__1 = *n;
#line 208 "zlansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "zlansb.f"
		sum = work[i__];
#line 210 "zlansb.f"
		if (value < sum || disnan_(&sum)) {
#line 210 "zlansb.f"
		    value = sum;
#line 210 "zlansb.f"
		}
#line 211 "zlansb.f"
/* L70: */
#line 211 "zlansb.f"
	    }
#line 212 "zlansb.f"
	} else {
#line 213 "zlansb.f"
	    i__1 = *n;
#line 213 "zlansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "zlansb.f"
		work[i__] = 0.;
#line 215 "zlansb.f"
/* L80: */
#line 215 "zlansb.f"
	    }
#line 216 "zlansb.f"
	    i__1 = *n;
#line 216 "zlansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 217 "zlansb.f"
		sum = work[j] + z_abs(&ab[j * ab_dim1 + 1]);
#line 218 "zlansb.f"
		l = 1 - j;
/* Computing MIN */
#line 219 "zlansb.f"
		i__3 = *n, i__2 = j + *k;
#line 219 "zlansb.f"
		i__4 = min(i__3,i__2);
#line 219 "zlansb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 220 "zlansb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 221 "zlansb.f"
		    sum += absa;
#line 222 "zlansb.f"
		    work[i__] += absa;
#line 223 "zlansb.f"
/* L90: */
#line 223 "zlansb.f"
		}
#line 224 "zlansb.f"
		if (value < sum || disnan_(&sum)) {
#line 224 "zlansb.f"
		    value = sum;
#line 224 "zlansb.f"
		}
#line 225 "zlansb.f"
/* L100: */
#line 225 "zlansb.f"
	    }
#line 226 "zlansb.f"
	}
#line 227 "zlansb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 231 "zlansb.f"
	scale = 0.;
#line 232 "zlansb.f"
	sum = 1.;
#line 233 "zlansb.f"
	if (*k > 0) {
#line 234 "zlansb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 235 "zlansb.f"
		i__1 = *n;
#line 235 "zlansb.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 236 "zlansb.f"
		    i__3 = j - 1;
#line 236 "zlansb.f"
		    i__4 = min(i__3,*k);
/* Computing MAX */
#line 236 "zlansb.f"
		    i__2 = *k + 2 - j;
#line 236 "zlansb.f"
		    zlassq_(&i__4, &ab[max(i__2,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 238 "zlansb.f"
/* L110: */
#line 238 "zlansb.f"
		}
#line 239 "zlansb.f"
		l = *k + 1;
#line 240 "zlansb.f"
	    } else {
#line 241 "zlansb.f"
		i__1 = *n - 1;
#line 241 "zlansb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 242 "zlansb.f"
		    i__3 = *n - j;
#line 242 "zlansb.f"
		    i__4 = min(i__3,*k);
#line 242 "zlansb.f"
		    zlassq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
#line 244 "zlansb.f"
/* L120: */
#line 244 "zlansb.f"
		}
#line 245 "zlansb.f"
		l = 1;
#line 246 "zlansb.f"
	    }
#line 247 "zlansb.f"
	    sum *= 2;
#line 248 "zlansb.f"
	} else {
#line 249 "zlansb.f"
	    l = 1;
#line 250 "zlansb.f"
	}
#line 251 "zlansb.f"
	zlassq_(n, &ab[l + ab_dim1], ldab, &scale, &sum);
#line 252 "zlansb.f"
	value = scale * sqrt(sum);
#line 253 "zlansb.f"
    }

#line 255 "zlansb.f"
    ret_val = value;
#line 256 "zlansb.f"
    return ret_val;

/*     End of ZLANSB */

} /* zlansb_ */

