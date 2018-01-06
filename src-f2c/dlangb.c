#line 1 "dlangb.f"
/* dlangb.f -- translated by f2c (version 20100827).
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

#line 1 "dlangb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of general band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlangb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlangb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlangb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANGB( NORM, N, KL, KU, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            KL, KU, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANGB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return DLANGB */
/* > \verbatim */
/* > */
/* >    DLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANGB as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANGB is */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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

/* > \ingroup doubleGBauxiliary */

/*  ===================================================================== */
doublereal dlangb_(char *norm, integer *n, integer *kl, integer *ku, 
	doublereal *ab, integer *ldab, doublereal *work, ftnlen norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal sum, temp, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
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

#line 163 "dlangb.f"
    /* Parameter adjustments */
#line 163 "dlangb.f"
    ab_dim1 = *ldab;
#line 163 "dlangb.f"
    ab_offset = 1 + ab_dim1;
#line 163 "dlangb.f"
    ab -= ab_offset;
#line 163 "dlangb.f"
    --work;
#line 163 "dlangb.f"

#line 163 "dlangb.f"
    /* Function Body */
#line 163 "dlangb.f"
    if (*n == 0) {
#line 164 "dlangb.f"
	value = 0.;
#line 165 "dlangb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 169 "dlangb.f"
	value = 0.;
#line 170 "dlangb.f"
	i__1 = *n;
#line 170 "dlangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 171 "dlangb.f"
	    i__2 = *ku + 2 - j;
/* Computing MIN */
#line 171 "dlangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 171 "dlangb.f"
	    i__3 = min(i__4,i__5);
#line 171 "dlangb.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 172 "dlangb.f"
		temp = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 173 "dlangb.f"
		if (value < temp || disnan_(&temp)) {
#line 173 "dlangb.f"
		    value = temp;
#line 173 "dlangb.f"
		}
#line 174 "dlangb.f"
/* L10: */
#line 174 "dlangb.f"
	    }
#line 175 "dlangb.f"
/* L20: */
#line 175 "dlangb.f"
	}
#line 176 "dlangb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 180 "dlangb.f"
	value = 0.;
#line 181 "dlangb.f"
	i__1 = *n;
#line 181 "dlangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 182 "dlangb.f"
	    sum = 0.;
/* Computing MAX */
#line 183 "dlangb.f"
	    i__3 = *ku + 2 - j;
/* Computing MIN */
#line 183 "dlangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 183 "dlangb.f"
	    i__2 = min(i__4,i__5);
#line 183 "dlangb.f"
	    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 184 "dlangb.f"
		sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 185 "dlangb.f"
/* L30: */
#line 185 "dlangb.f"
	    }
#line 186 "dlangb.f"
	    if (value < sum || disnan_(&sum)) {
#line 186 "dlangb.f"
		value = sum;
#line 186 "dlangb.f"
	    }
#line 187 "dlangb.f"
/* L40: */
#line 187 "dlangb.f"
	}
#line 188 "dlangb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 192 "dlangb.f"
	i__1 = *n;
#line 192 "dlangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 193 "dlangb.f"
	    work[i__] = 0.;
#line 194 "dlangb.f"
/* L50: */
#line 194 "dlangb.f"
	}
#line 195 "dlangb.f"
	i__1 = *n;
#line 195 "dlangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 196 "dlangb.f"
	    k = *ku + 1 - j;
/* Computing MAX */
#line 197 "dlangb.f"
	    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 197 "dlangb.f"
	    i__5 = *n, i__6 = j + *kl;
#line 197 "dlangb.f"
	    i__4 = min(i__5,i__6);
#line 197 "dlangb.f"
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 198 "dlangb.f"
		work[i__] += (d__1 = ab[k + i__ + j * ab_dim1], abs(d__1));
#line 199 "dlangb.f"
/* L60: */
#line 199 "dlangb.f"
	    }
#line 200 "dlangb.f"
/* L70: */
#line 200 "dlangb.f"
	}
#line 201 "dlangb.f"
	value = 0.;
#line 202 "dlangb.f"
	i__1 = *n;
#line 202 "dlangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "dlangb.f"
	    temp = work[i__];
#line 204 "dlangb.f"
	    if (value < temp || disnan_(&temp)) {
#line 204 "dlangb.f"
		value = temp;
#line 204 "dlangb.f"
	    }
#line 205 "dlangb.f"
/* L80: */
#line 205 "dlangb.f"
	}
#line 206 "dlangb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 210 "dlangb.f"
	scale = 0.;
#line 211 "dlangb.f"
	sum = 1.;
#line 212 "dlangb.f"
	i__1 = *n;
#line 212 "dlangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 213 "dlangb.f"
	    i__4 = 1, i__2 = j - *ku;
#line 213 "dlangb.f"
	    l = max(i__4,i__2);
#line 214 "dlangb.f"
	    k = *ku + 1 - j + l;
/* Computing MIN */
#line 215 "dlangb.f"
	    i__2 = *n, i__3 = j + *kl;
#line 215 "dlangb.f"
	    i__4 = min(i__2,i__3) - l + 1;
#line 215 "dlangb.f"
	    dlassq_(&i__4, &ab[k + j * ab_dim1], &c__1, &scale, &sum);
#line 216 "dlangb.f"
/* L90: */
#line 216 "dlangb.f"
	}
#line 217 "dlangb.f"
	value = scale * sqrt(sum);
#line 218 "dlangb.f"
    }

#line 220 "dlangb.f"
    ret_val = value;
#line 221 "dlangb.f"
    return ret_val;

/*     End of DLANGB */

} /* dlangb_ */

