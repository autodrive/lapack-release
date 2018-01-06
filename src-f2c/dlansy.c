#line 1 "dlansy.f"
/* dlansy.f -- translated by f2c (version 20100827).
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

#line 1 "dlansy.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a real symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANSY  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real symmetric matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANSY */
/* > \verbatim */
/* > */
/* >    DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANSY as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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

/* > \date September 2012 */

/* > \ingroup doubleSYauxiliary */

/*  ===================================================================== */
doublereal dlansy_(char *norm, char *uplo, integer *n, doublereal *a, integer 
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

#line 160 "dlansy.f"
    /* Parameter adjustments */
#line 160 "dlansy.f"
    a_dim1 = *lda;
#line 160 "dlansy.f"
    a_offset = 1 + a_dim1;
#line 160 "dlansy.f"
    a -= a_offset;
#line 160 "dlansy.f"
    --work;
#line 160 "dlansy.f"

#line 160 "dlansy.f"
    /* Function Body */
#line 160 "dlansy.f"
    if (*n == 0) {
#line 161 "dlansy.f"
	value = 0.;
#line 162 "dlansy.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 166 "dlansy.f"
	value = 0.;
#line 167 "dlansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 168 "dlansy.f"
	    i__1 = *n;
#line 168 "dlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 169 "dlansy.f"
		i__2 = j;
#line 169 "dlansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 170 "dlansy.f"
		    sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 171 "dlansy.f"
		    if (value < sum || disnan_(&sum)) {
#line 171 "dlansy.f"
			value = sum;
#line 171 "dlansy.f"
		    }
#line 172 "dlansy.f"
/* L10: */
#line 172 "dlansy.f"
		}
#line 173 "dlansy.f"
/* L20: */
#line 173 "dlansy.f"
	    }
#line 174 "dlansy.f"
	} else {
#line 175 "dlansy.f"
	    i__1 = *n;
#line 175 "dlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 176 "dlansy.f"
		i__2 = *n;
#line 176 "dlansy.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 177 "dlansy.f"
		    sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 178 "dlansy.f"
		    if (value < sum || disnan_(&sum)) {
#line 178 "dlansy.f"
			value = sum;
#line 178 "dlansy.f"
		    }
#line 179 "dlansy.f"
/* L30: */
#line 179 "dlansy.f"
		}
#line 180 "dlansy.f"
/* L40: */
#line 180 "dlansy.f"
	    }
#line 181 "dlansy.f"
	}
#line 182 "dlansy.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 187 "dlansy.f"
	value = 0.;
#line 188 "dlansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 189 "dlansy.f"
	    i__1 = *n;
#line 189 "dlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 190 "dlansy.f"
		sum = 0.;
#line 191 "dlansy.f"
		i__2 = j - 1;
#line 191 "dlansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 192 "dlansy.f"
		    absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 193 "dlansy.f"
		    sum += absa;
#line 194 "dlansy.f"
		    work[i__] += absa;
#line 195 "dlansy.f"
/* L50: */
#line 195 "dlansy.f"
		}
#line 196 "dlansy.f"
		work[j] = sum + (d__1 = a[j + j * a_dim1], abs(d__1));
#line 197 "dlansy.f"
/* L60: */
#line 197 "dlansy.f"
	    }
#line 198 "dlansy.f"
	    i__1 = *n;
#line 198 "dlansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 199 "dlansy.f"
		sum = work[i__];
#line 200 "dlansy.f"
		if (value < sum || disnan_(&sum)) {
#line 200 "dlansy.f"
		    value = sum;
#line 200 "dlansy.f"
		}
#line 201 "dlansy.f"
/* L70: */
#line 201 "dlansy.f"
	    }
#line 202 "dlansy.f"
	} else {
#line 203 "dlansy.f"
	    i__1 = *n;
#line 203 "dlansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "dlansy.f"
		work[i__] = 0.;
#line 205 "dlansy.f"
/* L80: */
#line 205 "dlansy.f"
	    }
#line 206 "dlansy.f"
	    i__1 = *n;
#line 206 "dlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 207 "dlansy.f"
		sum = work[j] + (d__1 = a[j + j * a_dim1], abs(d__1));
#line 208 "dlansy.f"
		i__2 = *n;
#line 208 "dlansy.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 209 "dlansy.f"
		    absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 210 "dlansy.f"
		    sum += absa;
#line 211 "dlansy.f"
		    work[i__] += absa;
#line 212 "dlansy.f"
/* L90: */
#line 212 "dlansy.f"
		}
#line 213 "dlansy.f"
		if (value < sum || disnan_(&sum)) {
#line 213 "dlansy.f"
		    value = sum;
#line 213 "dlansy.f"
		}
#line 214 "dlansy.f"
/* L100: */
#line 214 "dlansy.f"
	    }
#line 215 "dlansy.f"
	}
#line 216 "dlansy.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 220 "dlansy.f"
	scale = 0.;
#line 221 "dlansy.f"
	sum = 1.;
#line 222 "dlansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 223 "dlansy.f"
	    i__1 = *n;
#line 223 "dlansy.f"
	    for (j = 2; j <= i__1; ++j) {
#line 224 "dlansy.f"
		i__2 = j - 1;
#line 224 "dlansy.f"
		dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 225 "dlansy.f"
/* L110: */
#line 225 "dlansy.f"
	    }
#line 226 "dlansy.f"
	} else {
#line 227 "dlansy.f"
	    i__1 = *n - 1;
#line 227 "dlansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 228 "dlansy.f"
		i__2 = *n - j;
#line 228 "dlansy.f"
		dlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
#line 229 "dlansy.f"
/* L120: */
#line 229 "dlansy.f"
	    }
#line 230 "dlansy.f"
	}
#line 231 "dlansy.f"
	sum *= 2;
#line 232 "dlansy.f"
	i__1 = *lda + 1;
#line 232 "dlansy.f"
	dlassq_(n, &a[a_offset], &i__1, &scale, &sum);
#line 233 "dlansy.f"
	value = scale * sqrt(sum);
#line 234 "dlansy.f"
    }

#line 236 "dlansy.f"
    ret_val = value;
#line 237 "dlansy.f"
    return ret_val;

/*     End of DLANSY */

} /* dlansy_ */

