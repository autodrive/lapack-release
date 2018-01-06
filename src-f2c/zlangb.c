#line 1 "zlangb.f"
/* zlangb.f -- translated by f2c (version 20100827).
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

#line 1 "zlangb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of general band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlangb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlangb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlangb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANGB( NORM, N, KL, KU, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            KL, KU, LDAB, N */
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
/* > ZLANGB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return ZLANGB */
/* > \verbatim */
/* > */
/* >    ZLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANGB as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANGB is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of sub-diagonals of the matrix A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of super-diagonals of the matrix A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th */
/* >          column of A is stored in the j-th column of the array AB as */
/* >          follows: */
/* >          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= N when NORM = 'I'; otherwise, WORK is not */
/* >          referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GBauxiliary */

/*  ===================================================================== */
doublereal zlangb_(char *norm, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublereal *work, ftnlen norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal sum, temp, scale;
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

#line 164 "zlangb.f"
    /* Parameter adjustments */
#line 164 "zlangb.f"
    ab_dim1 = *ldab;
#line 164 "zlangb.f"
    ab_offset = 1 + ab_dim1;
#line 164 "zlangb.f"
    ab -= ab_offset;
#line 164 "zlangb.f"
    --work;
#line 164 "zlangb.f"

#line 164 "zlangb.f"
    /* Function Body */
#line 164 "zlangb.f"
    if (*n == 0) {
#line 165 "zlangb.f"
	value = 0.;
#line 166 "zlangb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 170 "zlangb.f"
	value = 0.;
#line 171 "zlangb.f"
	i__1 = *n;
#line 171 "zlangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 172 "zlangb.f"
	    i__2 = *ku + 2 - j;
/* Computing MIN */
#line 172 "zlangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 172 "zlangb.f"
	    i__3 = min(i__4,i__5);
#line 172 "zlangb.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 173 "zlangb.f"
		temp = z_abs(&ab[i__ + j * ab_dim1]);
#line 174 "zlangb.f"
		if (value < temp || disnan_(&temp)) {
#line 174 "zlangb.f"
		    value = temp;
#line 174 "zlangb.f"
		}
#line 175 "zlangb.f"
/* L10: */
#line 175 "zlangb.f"
	    }
#line 176 "zlangb.f"
/* L20: */
#line 176 "zlangb.f"
	}
#line 177 "zlangb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 181 "zlangb.f"
	value = 0.;
#line 182 "zlangb.f"
	i__1 = *n;
#line 182 "zlangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 183 "zlangb.f"
	    sum = 0.;
/* Computing MAX */
#line 184 "zlangb.f"
	    i__3 = *ku + 2 - j;
/* Computing MIN */
#line 184 "zlangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 184 "zlangb.f"
	    i__2 = min(i__4,i__5);
#line 184 "zlangb.f"
	    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 185 "zlangb.f"
		sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 186 "zlangb.f"
/* L30: */
#line 186 "zlangb.f"
	    }
#line 187 "zlangb.f"
	    if (value < sum || disnan_(&sum)) {
#line 187 "zlangb.f"
		value = sum;
#line 187 "zlangb.f"
	    }
#line 188 "zlangb.f"
/* L40: */
#line 188 "zlangb.f"
	}
#line 189 "zlangb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 193 "zlangb.f"
	i__1 = *n;
#line 193 "zlangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "zlangb.f"
	    work[i__] = 0.;
#line 195 "zlangb.f"
/* L50: */
#line 195 "zlangb.f"
	}
#line 196 "zlangb.f"
	i__1 = *n;
#line 196 "zlangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 197 "zlangb.f"
	    k = *ku + 1 - j;
/* Computing MAX */
#line 198 "zlangb.f"
	    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 198 "zlangb.f"
	    i__5 = *n, i__6 = j + *kl;
#line 198 "zlangb.f"
	    i__4 = min(i__5,i__6);
#line 198 "zlangb.f"
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 199 "zlangb.f"
		work[i__] += z_abs(&ab[k + i__ + j * ab_dim1]);
#line 200 "zlangb.f"
/* L60: */
#line 200 "zlangb.f"
	    }
#line 201 "zlangb.f"
/* L70: */
#line 201 "zlangb.f"
	}
#line 202 "zlangb.f"
	value = 0.;
#line 203 "zlangb.f"
	i__1 = *n;
#line 203 "zlangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "zlangb.f"
	    temp = work[i__];
#line 205 "zlangb.f"
	    if (value < temp || disnan_(&temp)) {
#line 205 "zlangb.f"
		value = temp;
#line 205 "zlangb.f"
	    }
#line 206 "zlangb.f"
/* L80: */
#line 206 "zlangb.f"
	}
#line 207 "zlangb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 211 "zlangb.f"
	scale = 0.;
#line 212 "zlangb.f"
	sum = 1.;
#line 213 "zlangb.f"
	i__1 = *n;
#line 213 "zlangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 214 "zlangb.f"
	    i__4 = 1, i__2 = j - *ku;
#line 214 "zlangb.f"
	    l = max(i__4,i__2);
#line 215 "zlangb.f"
	    k = *ku + 1 - j + l;
/* Computing MIN */
#line 216 "zlangb.f"
	    i__2 = *n, i__3 = j + *kl;
#line 216 "zlangb.f"
	    i__4 = min(i__2,i__3) - l + 1;
#line 216 "zlangb.f"
	    zlassq_(&i__4, &ab[k + j * ab_dim1], &c__1, &scale, &sum);
#line 217 "zlangb.f"
/* L90: */
#line 217 "zlangb.f"
	}
#line 218 "zlangb.f"
	value = scale * sqrt(sum);
#line 219 "zlangb.f"
    }

#line 221 "zlangb.f"
    ret_val = value;
#line 222 "zlangb.f"
    return ret_val;

/*     End of ZLANGB */

} /* zlangb_ */

