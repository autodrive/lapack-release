#line 1 "zlansy.f"
/* zlansy.f -- translated by f2c (version 20100827).
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

#line 1 "zlansy.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlansy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlansy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlansy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANSY( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANSY  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex symmetric matrix A. */
/* > \endverbatim */
/* > */
/* > \return ZLANSY */
/* > \verbatim */
/* > */
/* >    ZLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANSY as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is to be referenced. */
/* >          = 'U':  Upper triangular part of A is referenced */
/* >          = 'L':  Lower triangular part of A is referenced */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANSY is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The symmetric matrix A.  If UPLO = 'U', the leading n by n */
/* >          upper triangular part of A contains the upper triangular part */
/* >          of the matrix A, and the strictly lower triangular part of A */
/* >          is not referenced.  If UPLO = 'L', the leading n by n lower */
/* >          triangular part of A contains the lower triangular part of */
/* >          the matrix A, and the strictly upper triangular part of A is */
/* >          not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(N,1). */
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

/* > \ingroup complex16SYauxiliary */

/*  ===================================================================== */
doublereal zlansy_(char *norm, char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
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

#line 162 "zlansy.f"
    /* Parameter adjustments */
#line 162 "zlansy.f"
    a_dim1 = *lda;
#line 162 "zlansy.f"
    a_offset = 1 + a_dim1;
#line 162 "zlansy.f"
    a -= a_offset;
#line 162 "zlansy.f"
    --work;
#line 162 "zlansy.f"

#line 162 "zlansy.f"
    /* Function Body */
#line 162 "zlansy.f"
    if (*n == 0) {
#line 163 "zlansy.f"
	value = 0.;
#line 164 "zlansy.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 168 "zlansy.f"
	value = 0.;
#line 169 "zlansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 170 "zlansy.f"
	    i__1 = *n;
#line 170 "zlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 171 "zlansy.f"
		i__2 = j;
#line 171 "zlansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 172 "zlansy.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 173 "zlansy.f"
		    if (value < sum || disnan_(&sum)) {
#line 173 "zlansy.f"
			value = sum;
#line 173 "zlansy.f"
		    }
#line 174 "zlansy.f"
/* L10: */
#line 174 "zlansy.f"
		}
#line 175 "zlansy.f"
/* L20: */
#line 175 "zlansy.f"
	    }
#line 176 "zlansy.f"
	} else {
#line 177 "zlansy.f"
	    i__1 = *n;
#line 177 "zlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 178 "zlansy.f"
		i__2 = *n;
#line 178 "zlansy.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 179 "zlansy.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 180 "zlansy.f"
		    if (value < sum || disnan_(&sum)) {
#line 180 "zlansy.f"
			value = sum;
#line 180 "zlansy.f"
		    }
#line 181 "zlansy.f"
/* L30: */
#line 181 "zlansy.f"
		}
#line 182 "zlansy.f"
/* L40: */
#line 182 "zlansy.f"
	    }
#line 183 "zlansy.f"
	}
#line 184 "zlansy.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 189 "zlansy.f"
	value = 0.;
#line 190 "zlansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 191 "zlansy.f"
	    i__1 = *n;
#line 191 "zlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 192 "zlansy.f"
		sum = 0.;
#line 193 "zlansy.f"
		i__2 = j - 1;
#line 193 "zlansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 194 "zlansy.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 195 "zlansy.f"
		    sum += absa;
#line 196 "zlansy.f"
		    work[i__] += absa;
#line 197 "zlansy.f"
/* L50: */
#line 197 "zlansy.f"
		}
#line 198 "zlansy.f"
		work[j] = sum + z_abs(&a[j + j * a_dim1]);
#line 199 "zlansy.f"
/* L60: */
#line 199 "zlansy.f"
	    }
#line 200 "zlansy.f"
	    i__1 = *n;
#line 200 "zlansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 201 "zlansy.f"
		sum = work[i__];
#line 202 "zlansy.f"
		if (value < sum || disnan_(&sum)) {
#line 202 "zlansy.f"
		    value = sum;
#line 202 "zlansy.f"
		}
#line 203 "zlansy.f"
/* L70: */
#line 203 "zlansy.f"
	    }
#line 204 "zlansy.f"
	} else {
#line 205 "zlansy.f"
	    i__1 = *n;
#line 205 "zlansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "zlansy.f"
		work[i__] = 0.;
#line 207 "zlansy.f"
/* L80: */
#line 207 "zlansy.f"
	    }
#line 208 "zlansy.f"
	    i__1 = *n;
#line 208 "zlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 209 "zlansy.f"
		sum = work[j] + z_abs(&a[j + j * a_dim1]);
#line 210 "zlansy.f"
		i__2 = *n;
#line 210 "zlansy.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 211 "zlansy.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 212 "zlansy.f"
		    sum += absa;
#line 213 "zlansy.f"
		    work[i__] += absa;
#line 214 "zlansy.f"
/* L90: */
#line 214 "zlansy.f"
		}
#line 215 "zlansy.f"
		if (value < sum || disnan_(&sum)) {
#line 215 "zlansy.f"
		    value = sum;
#line 215 "zlansy.f"
		}
#line 216 "zlansy.f"
/* L100: */
#line 216 "zlansy.f"
	    }
#line 217 "zlansy.f"
	}
#line 218 "zlansy.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 222 "zlansy.f"
	scale = 0.;
#line 223 "zlansy.f"
	sum = 1.;
#line 224 "zlansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 225 "zlansy.f"
	    i__1 = *n;
#line 225 "zlansy.f"
	    for (j = 2; j <= i__1; ++j) {
#line 226 "zlansy.f"
		i__2 = j - 1;
#line 226 "zlansy.f"
		zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 227 "zlansy.f"
/* L110: */
#line 227 "zlansy.f"
	    }
#line 228 "zlansy.f"
	} else {
#line 229 "zlansy.f"
	    i__1 = *n - 1;
#line 229 "zlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 230 "zlansy.f"
		i__2 = *n - j;
#line 230 "zlansy.f"
		zlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
#line 231 "zlansy.f"
/* L120: */
#line 231 "zlansy.f"
	    }
#line 232 "zlansy.f"
	}
#line 233 "zlansy.f"
	sum *= 2;
#line 234 "zlansy.f"
	i__1 = *lda + 1;
#line 234 "zlansy.f"
	zlassq_(n, &a[a_offset], &i__1, &scale, &sum);
#line 235 "zlansy.f"
	value = scale * sqrt(sum);
#line 236 "zlansy.f"
    }

#line 238 "zlansy.f"
    ret_val = value;
#line 239 "zlansy.f"
    return ret_val;

/*     End of ZLANSY */

} /* zlansy_ */

