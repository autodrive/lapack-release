#line 1 "clansy.f"
/* clansy.f -- translated by f2c (version 20100827).
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

#line 1 "clansy.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANSY( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               WORK( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANSY  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex symmetric matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANSY */
/* > \verbatim */
/* > */
/* >    CLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANSY as described */
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
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANSY is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup complexSYauxiliary */

/*  ===================================================================== */
doublereal clansy_(char *norm, char *uplo, integer *n, doublecomplex *a, 
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

#line 162 "clansy.f"
    /* Parameter adjustments */
#line 162 "clansy.f"
    a_dim1 = *lda;
#line 162 "clansy.f"
    a_offset = 1 + a_dim1;
#line 162 "clansy.f"
    a -= a_offset;
#line 162 "clansy.f"
    --work;
#line 162 "clansy.f"

#line 162 "clansy.f"
    /* Function Body */
#line 162 "clansy.f"
    if (*n == 0) {
#line 163 "clansy.f"
	value = 0.;
#line 164 "clansy.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 168 "clansy.f"
	value = 0.;
#line 169 "clansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 170 "clansy.f"
	    i__1 = *n;
#line 170 "clansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 171 "clansy.f"
		i__2 = j;
#line 171 "clansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 172 "clansy.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 173 "clansy.f"
		    if (value < sum || sisnan_(&sum)) {
#line 173 "clansy.f"
			value = sum;
#line 173 "clansy.f"
		    }
#line 174 "clansy.f"
/* L10: */
#line 174 "clansy.f"
		}
#line 175 "clansy.f"
/* L20: */
#line 175 "clansy.f"
	    }
#line 176 "clansy.f"
	} else {
#line 177 "clansy.f"
	    i__1 = *n;
#line 177 "clansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 178 "clansy.f"
		i__2 = *n;
#line 178 "clansy.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 179 "clansy.f"
		    sum = z_abs(&a[i__ + j * a_dim1]);
#line 180 "clansy.f"
		    if (value < sum || sisnan_(&sum)) {
#line 180 "clansy.f"
			value = sum;
#line 180 "clansy.f"
		    }
#line 181 "clansy.f"
/* L30: */
#line 181 "clansy.f"
		}
#line 182 "clansy.f"
/* L40: */
#line 182 "clansy.f"
	    }
#line 183 "clansy.f"
	}
#line 184 "clansy.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 189 "clansy.f"
	value = 0.;
#line 190 "clansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 191 "clansy.f"
	    i__1 = *n;
#line 191 "clansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 192 "clansy.f"
		sum = 0.;
#line 193 "clansy.f"
		i__2 = j - 1;
#line 193 "clansy.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 194 "clansy.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 195 "clansy.f"
		    sum += absa;
#line 196 "clansy.f"
		    work[i__] += absa;
#line 197 "clansy.f"
/* L50: */
#line 197 "clansy.f"
		}
#line 198 "clansy.f"
		work[j] = sum + z_abs(&a[j + j * a_dim1]);
#line 199 "clansy.f"
/* L60: */
#line 199 "clansy.f"
	    }
#line 200 "clansy.f"
	    i__1 = *n;
#line 200 "clansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 201 "clansy.f"
		sum = work[i__];
#line 202 "clansy.f"
		if (value < sum || sisnan_(&sum)) {
#line 202 "clansy.f"
		    value = sum;
#line 202 "clansy.f"
		}
#line 203 "clansy.f"
/* L70: */
#line 203 "clansy.f"
	    }
#line 204 "clansy.f"
	} else {
#line 205 "clansy.f"
	    i__1 = *n;
#line 205 "clansy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 206 "clansy.f"
		work[i__] = 0.;
#line 207 "clansy.f"
/* L80: */
#line 207 "clansy.f"
	    }
#line 208 "clansy.f"
	    i__1 = *n;
#line 208 "clansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 209 "clansy.f"
		sum = work[j] + z_abs(&a[j + j * a_dim1]);
#line 210 "clansy.f"
		i__2 = *n;
#line 210 "clansy.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 211 "clansy.f"
		    absa = z_abs(&a[i__ + j * a_dim1]);
#line 212 "clansy.f"
		    sum += absa;
#line 213 "clansy.f"
		    work[i__] += absa;
#line 214 "clansy.f"
/* L90: */
#line 214 "clansy.f"
		}
#line 215 "clansy.f"
		if (value < sum || sisnan_(&sum)) {
#line 215 "clansy.f"
		    value = sum;
#line 215 "clansy.f"
		}
#line 216 "clansy.f"
/* L100: */
#line 216 "clansy.f"
	    }
#line 217 "clansy.f"
	}
#line 218 "clansy.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 222 "clansy.f"
	scale = 0.;
#line 223 "clansy.f"
	sum = 1.;
#line 224 "clansy.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 225 "clansy.f"
	    i__1 = *n;
#line 225 "clansy.f"
	    for (j = 2; j <= i__1; ++j) {
#line 226 "clansy.f"
		i__2 = j - 1;
#line 226 "clansy.f"
		classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 227 "clansy.f"
/* L110: */
#line 227 "clansy.f"
	    }
#line 228 "clansy.f"
	} else {
#line 229 "clansy.f"
	    i__1 = *n - 1;
#line 229 "clansy.f"
	    for (j = 1; j <= i__1; ++j) {
#line 230 "clansy.f"
		i__2 = *n - j;
#line 230 "clansy.f"
		classq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
#line 231 "clansy.f"
/* L120: */
#line 231 "clansy.f"
	    }
#line 232 "clansy.f"
	}
#line 233 "clansy.f"
	sum *= 2;
#line 234 "clansy.f"
	i__1 = *lda + 1;
#line 234 "clansy.f"
	classq_(n, &a[a_offset], &i__1, &scale, &sum);
#line 235 "clansy.f"
	value = scale * sqrt(sum);
#line 236 "clansy.f"
    }

#line 238 "clansy.f"
    ret_val = value;
#line 239 "clansy.f"
    return ret_val;

/*     End of CLANSY */

} /* clansy_ */

