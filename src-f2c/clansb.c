#line 1 "clansb.f"
/* clansb.f -- translated by f2c (version 20100827).
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

#line 1 "clansb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANSB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANSB( NORM, UPLO, N, K, AB, LDAB, */
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
/* > CLANSB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n symmetric band matrix A,  with k super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return CLANSB */
/* > \verbatim */
/* > */
/* >    CLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANSB as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANSB is */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clansb_(char *norm, char *uplo, integer *n, integer *k, 
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
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    extern logical sisnan_(doublereal *);


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

#line 169 "clansb.f"
    /* Parameter adjustments */
#line 169 "clansb.f"
    ab_dim1 = *ldab;
#line 169 "clansb.f"
    ab_offset = 1 + ab_dim1;
#line 169 "clansb.f"
    ab -= ab_offset;
#line 169 "clansb.f"
    --work;
#line 169 "clansb.f"

#line 169 "clansb.f"
    /* Function Body */
#line 169 "clansb.f"
    if (*n == 0) {
#line 170 "clansb.f"
	value = 0.;
#line 171 "clansb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 175 "clansb.f"
	value = 0.;
#line 176 "clansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 177 "clansb.f"
	    i__1 = *n;
#line 177 "clansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 178 "clansb.f"
		i__2 = *k + 2 - j;
#line 178 "clansb.f"
		i__3 = *k + 1;
#line 178 "clansb.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 179 "clansb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 180 "clansb.f"
		    if (value < sum || sisnan_(&sum)) {
#line 180 "clansb.f"
			value = sum;
#line 180 "clansb.f"
		    }
#line 181 "clansb.f"
/* L10: */
#line 181 "clansb.f"
		}
#line 182 "clansb.f"
/* L20: */
#line 182 "clansb.f"
	    }
#line 183 "clansb.f"
	} else {
#line 184 "clansb.f"
	    i__1 = *n;
#line 184 "clansb.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 185 "clansb.f"
		i__2 = *n + 1 - j, i__4 = *k + 1;
#line 185 "clansb.f"
		i__3 = min(i__2,i__4);
#line 185 "clansb.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 186 "clansb.f"
		    sum = z_abs(&ab[i__ + j * ab_dim1]);
#line 187 "clansb.f"
		    if (value < sum || sisnan_(&sum)) {
#line 187 "clansb.f"
			value = sum;
#line 187 "clansb.f"
		    }
#line 188 "clansb.f"
/* L30: */
#line 188 "clansb.f"
		}
#line 189 "clansb.f"
/* L40: */
#line 189 "clansb.f"
	    }
#line 190 "clansb.f"
	}
#line 191 "clansb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 196 "clansb.f"
	value = 0.;
#line 197 "clansb.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 198 "clansb.f"
	    i__1 = *n;
#line 198 "clansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 199 "clansb.f"
		sum = 0.;
#line 200 "clansb.f"
		l = *k + 1 - j;
/* Computing MAX */
#line 201 "clansb.f"
		i__3 = 1, i__2 = j - *k;
#line 201 "clansb.f"
		i__4 = j - 1;
#line 201 "clansb.f"
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
#line 202 "clansb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 203 "clansb.f"
		    sum += absa;
#line 204 "clansb.f"
		    work[i__] += absa;
#line 205 "clansb.f"
/* L50: */
#line 205 "clansb.f"
		}
#line 206 "clansb.f"
		work[j] = sum + z_abs(&ab[*k + 1 + j * ab_dim1]);
#line 207 "clansb.f"
/* L60: */
#line 207 "clansb.f"
	    }
#line 208 "clansb.f"
	    i__1 = *n;
#line 208 "clansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "clansb.f"
		sum = work[i__];
#line 210 "clansb.f"
		if (value < sum || sisnan_(&sum)) {
#line 210 "clansb.f"
		    value = sum;
#line 210 "clansb.f"
		}
#line 211 "clansb.f"
/* L70: */
#line 211 "clansb.f"
	    }
#line 212 "clansb.f"
	} else {
#line 213 "clansb.f"
	    i__1 = *n;
#line 213 "clansb.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "clansb.f"
		work[i__] = 0.;
#line 215 "clansb.f"
/* L80: */
#line 215 "clansb.f"
	    }
#line 216 "clansb.f"
	    i__1 = *n;
#line 216 "clansb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 217 "clansb.f"
		sum = work[j] + z_abs(&ab[j * ab_dim1 + 1]);
#line 218 "clansb.f"
		l = 1 - j;
/* Computing MIN */
#line 219 "clansb.f"
		i__3 = *n, i__2 = j + *k;
#line 219 "clansb.f"
		i__4 = min(i__3,i__2);
#line 219 "clansb.f"
		for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 220 "clansb.f"
		    absa = z_abs(&ab[l + i__ + j * ab_dim1]);
#line 221 "clansb.f"
		    sum += absa;
#line 222 "clansb.f"
		    work[i__] += absa;
#line 223 "clansb.f"
/* L90: */
#line 223 "clansb.f"
		}
#line 224 "clansb.f"
		if (value < sum || sisnan_(&sum)) {
#line 224 "clansb.f"
		    value = sum;
#line 224 "clansb.f"
		}
#line 225 "clansb.f"
/* L100: */
#line 225 "clansb.f"
	    }
#line 226 "clansb.f"
	}
#line 227 "clansb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 231 "clansb.f"
	scale = 0.;
#line 232 "clansb.f"
	sum = 1.;
#line 233 "clansb.f"
	if (*k > 0) {
#line 234 "clansb.f"
	    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 235 "clansb.f"
		i__1 = *n;
#line 235 "clansb.f"
		for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 236 "clansb.f"
		    i__3 = j - 1;
#line 236 "clansb.f"
		    i__4 = min(i__3,*k);
/* Computing MAX */
#line 236 "clansb.f"
		    i__2 = *k + 2 - j;
#line 236 "clansb.f"
		    classq_(&i__4, &ab[max(i__2,1) + j * ab_dim1], &c__1, &
			    scale, &sum);
#line 238 "clansb.f"
/* L110: */
#line 238 "clansb.f"
		}
#line 239 "clansb.f"
		l = *k + 1;
#line 240 "clansb.f"
	    } else {
#line 241 "clansb.f"
		i__1 = *n - 1;
#line 241 "clansb.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 242 "clansb.f"
		    i__3 = *n - j;
#line 242 "clansb.f"
		    i__4 = min(i__3,*k);
#line 242 "clansb.f"
		    classq_(&i__4, &ab[j * ab_dim1 + 2], &c__1, &scale, &sum);
#line 244 "clansb.f"
/* L120: */
#line 244 "clansb.f"
		}
#line 245 "clansb.f"
		l = 1;
#line 246 "clansb.f"
	    }
#line 247 "clansb.f"
	    sum *= 2;
#line 248 "clansb.f"
	} else {
#line 249 "clansb.f"
	    l = 1;
#line 250 "clansb.f"
	}
#line 251 "clansb.f"
	classq_(n, &ab[l + ab_dim1], ldab, &scale, &sum);
#line 252 "clansb.f"
	value = scale * sqrt(sum);
#line 253 "clansb.f"
    }

#line 255 "clansb.f"
    ret_val = value;
#line 256 "clansb.f"
    return ret_val;

/*     End of CLANSB */

} /* clansb_ */

