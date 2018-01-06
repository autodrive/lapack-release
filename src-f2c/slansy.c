#line 1 "slansy.f"
/* slansy.f -- translated by f2c (version 20100827).
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

#line 1 "slansy.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a real symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANSY( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANSY  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real symmetric matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANSY */
/* > \verbatim */
/* > */
/* >    SLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANSY as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANSY is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \date November 2015 */

/* > \ingroup realSYauxiliary */

/*  ===================================================================== */
doublereal slansy_(char *norm, char *uplo, integer *n, doublereal *a, integer 
	*lda, doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal sum, absa, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 160 "slansy.f"
    /* Parameter adjustments */
#line 160 "slansy.f"
    a_dim1 = *lda;
#line 160 "slansy.f"
    a_offset = 1 + a_dim1;
#line 160 "slansy.f"
    a -= a_offset;
#line 160 "slansy.f"
    --work;
#line 160 "slansy.f"

#line 160 "slansy.f"
    /* Function Body */
#line 160 "slansy.f"
    if (*n == 0) {
#line 161 "slansy.f"
	value = 0.;
#line 162 "slansy.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 166 "slansy.f"
	value = 0.;
#line 167 "slansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 168 "slansy.f"
	    i__1 = *n;
#line 168 "slansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 169 "slansy.f"
		i__2 = j;
#line 169 "slansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 170 "slansy.f"
		    sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 171 "slansy.f"
		    if (value < sum || sisnan_(&sum)) {
#line 171 "slansy.f"
			value = sum;
#line 171 "slansy.f"
		    }
#line 172 "slansy.f"
/* L10: */
#line 172 "slansy.f"
		}
#line 173 "slansy.f"
/* L20: */
#line 173 "slansy.f"
	    }
#line 174 "slansy.f"
	} else {
#line 175 "slansy.f"
	    i__1 = *n;
#line 175 "slansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 176 "slansy.f"
		i__2 = *n;
#line 176 "slansy.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 177 "slansy.f"
		    sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 178 "slansy.f"
		    if (value < sum || sisnan_(&sum)) {
#line 178 "slansy.f"
			value = sum;
#line 178 "slansy.f"
		    }
#line 179 "slansy.f"
/* L30: */
#line 179 "slansy.f"
		}
#line 180 "slansy.f"
/* L40: */
#line 180 "slansy.f"
	    }
#line 181 "slansy.f"
	}
#line 182 "slansy.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 187 "slansy.f"
	value = 0.;
#line 188 "slansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "slansy.f"
	    i__1 = *n;
#line 189 "slansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 190 "slansy.f"
		sum = 0.;
#line 191 "slansy.f"
		i__2 = j - 1;
#line 191 "slansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 192 "slansy.f"
		    absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 193 "slansy.f"
		    sum += absa;
#line 194 "slansy.f"
		    work[i__] += absa;
#line 195 "slansy.f"
/* L50: */
#line 195 "slansy.f"
		}
#line 196 "slansy.f"
		work[j] = sum + (d__1 = a[j + j * a_dim1], abs(d__1));
#line 197 "slansy.f"
/* L60: */
#line 197 "slansy.f"
	    }
#line 198 "slansy.f"
	    i__1 = *n;
#line 198 "slansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 199 "slansy.f"
		sum = work[i__];
#line 200 "slansy.f"
		if (value < sum || sisnan_(&sum)) {
#line 200 "slansy.f"
		    value = sum;
#line 200 "slansy.f"
		}
#line 201 "slansy.f"
/* L70: */
#line 201 "slansy.f"
	    }
#line 202 "slansy.f"
	} else {
#line 203 "slansy.f"
	    i__1 = *n;
#line 203 "slansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "slansy.f"
		work[i__] = 0.;
#line 205 "slansy.f"
/* L80: */
#line 205 "slansy.f"
	    }
#line 206 "slansy.f"
	    i__1 = *n;
#line 206 "slansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 207 "slansy.f"
		sum = work[j] + (d__1 = a[j + j * a_dim1], abs(d__1));
#line 208 "slansy.f"
		i__2 = *n;
#line 208 "slansy.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 209 "slansy.f"
		    absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 210 "slansy.f"
		    sum += absa;
#line 211 "slansy.f"
		    work[i__] += absa;
#line 212 "slansy.f"
/* L90: */
#line 212 "slansy.f"
		}
#line 213 "slansy.f"
		if (value < sum || sisnan_(&sum)) {
#line 213 "slansy.f"
		    value = sum;
#line 213 "slansy.f"
		}
#line 214 "slansy.f"
/* L100: */
#line 214 "slansy.f"
	    }
#line 215 "slansy.f"
	}
#line 216 "slansy.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 220 "slansy.f"
	scale = 0.;
#line 221 "slansy.f"
	sum = 1.;
#line 222 "slansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 223 "slansy.f"
	    i__1 = *n;
#line 223 "slansy.f"
	    for (j = 2; j <= i__1; ++j) {
#line 224 "slansy.f"
		i__2 = j - 1;
#line 224 "slansy.f"
		slassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 225 "slansy.f"
/* L110: */
#line 225 "slansy.f"
	    }
#line 226 "slansy.f"
	} else {
#line 227 "slansy.f"
	    i__1 = *n - 1;
#line 227 "slansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 228 "slansy.f"
		i__2 = *n - j;
#line 228 "slansy.f"
		slassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
#line 229 "slansy.f"
/* L120: */
#line 229 "slansy.f"
	    }
#line 230 "slansy.f"
	}
#line 231 "slansy.f"
	sum *= 2;
#line 232 "slansy.f"
	i__1 = *lda + 1;
#line 232 "slansy.f"
	slassq_(n, &a[a_offset], &i__1, &scale, &sum);
#line 233 "slansy.f"
	value = scale * sqrt(sum);
#line 234 "slansy.f"
    }

#line 236 "slansy.f"
    ret_val = value;
#line 237 "slansy.f"
    return ret_val;

/*     End of SLANSY */

} /* slansy_ */

